module module_sf_noahmp_init

!==========================================================================
! this module contains major data structure for noahmp initialized driver,
! and the subroutine noahmp_init
! modified from WRF   cheyz 2019/10/26 
!==========================================================================

 use grist_domain_types,  only: global_domain
 !use grist_static_vars_module 
 !use grist_data_types,    only: scalar_field, scalar_2d_field, scalar_3d_field
 !use grist_static_vars_module
 use noahmp_nml_control,  only: nsoil ,dzs1  ! cheyz: change dzs to dzs1
 !---cheyz--
 use grist_lsm_noahmp_vars
 use grist_constants, only: r8
 use grist_mpi
 
 implicit none
 public :: noahmp_init
   

   contains
   
   subroutine noahmp_init ( mminlu, snow , snowh , canwat , isltyp ,   ivgtyp,   & !MMINLU !, isltyp ,   ivgtyp=就是土壤类型和植被类型  lu_index and SOILTOP in grsist
      tslb , smois , sh2o , dzs , fndsoilw , fndsnowh ,                          &
      tsk, isnowxy , tvxy     ,tgxy     ,canicexy ,         tmn,     xice,       &
      canliqxy ,eahxy    ,tahxy    ,cmxy     ,chxy     ,                         &
      fwetxy   ,sneqvoxy ,alboldxy ,qsnowxy  ,wslakexy ,zwtxy    ,waxy     ,     &
      wtxy     ,tsnoxy   ,zsnsoxy  ,snicexy  ,snliqxy  ,lfmassxy ,rtmassxy ,     &
      stmassxy ,woodxy   ,stblcpxy ,fastcpxy ,xsaixy   ,lai      ,               &
      grainxy  ,gddxy    ,                                                       & !jref:start
      t2mvxy   ,t2mbxy   ,chstarxy,                                              & !jref:end       
      nsoil, restart,                                                            &
      iopt_run,                                                                  &
      ids,ide, jds,jde, kds,kde,                                                 &
      ims,ime, jms,jme, kms,kme,                                                 &
      its,ite, jts,jte, kts,kte,                                                 &
      smoiseq  ,smcwtdxy ,rechxy   ,deeprechxy, areaxy, dx, dy, msftx, msfty,    &     ! optional groundwater
      wtddt    ,stepwtd  ,dt       ,qrfsxy     ,qspringsxy  , qslatxy    ,       &      ! optional groundwater
      fdepthxy ,ht     ,riverbedxy ,eqzwt     ,rivercondxy ,pexpxy            )    ! optional groundwater

      use noahmp_tables

      implicit none

      ! initializing canopy air temperature to 287 k seems dangerous to me [kwm].

      integer, intent(in   )    ::   ids,ide, jds,jde, kds,kde,  &
         &                           ims,ime, jms,jme, kms,kme,  &
         &                           its,ite, jts,jte, kts,kte
      integer, intent(in)       ::   nsoil, iopt_run

      logical, intent(in)       ::   restart   !,                    &
         !&                           !allowed_to_read

      real(r8)  ,    dimension( nsoil), intent(in)    ::     dzs  ! thickness of the soil layers [m]
      real(r8)  ,    intent(in) , optional ::     dx, dy
      real(r8)  ,    dimension( ims:ime, jms:jme ) ,  intent(in) , optional :: msftx,msfty

      real(r8)  ,    dimension( ims:ime, nsoil, jms:jme ) ,    &
         &   intent(inout)    ::     smois,                      &  !TSLB SMOIS分别是4层土壤温度和土壤湿度 
         &                           sh2o,                       &  !一般初始值跟smois相等 是土壤液态水
         &                           tslb

      real(r8)  ,    dimension( ims:ime, jms:jme ) ,                     &
         &   intent(inout)    ::     snow,                       &
         &                           snowh,                      &
         &                           canwat  !0.0 ???

      integer, dimension( ims:ime, jms:jme ),                      &
         &   intent(in)         ::     isltyp,  &
                                       ivgtyp

      logical, intent(in)       ::     fndsoilw,                   &
         &                             fndsnowh

      real(r8)  , dimension(ims:ime,jms:jme), intent(in) :: tsk            !skin temperature (k)
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: tmn         !deep soil temperature (k)
      real(r8)  , dimension(ims:ime,jms:jme), intent(in) :: xice           !sea ice fraction
      integer   , dimension(ims:ime,jms:jme), intent(inout) :: isnowxy        !actual no. of snow layers
      real(r8)  , dimension(ims:ime,-2:nsoil,jms:jme), intent(inout) :: zsnsoxy  !snow layer depth [m]
      real(r8)  , dimension(ims:ime,-2:              0,jms:jme), intent(inout) :: tsnoxy   !snow temperature [k]
      real(r8)  , dimension(ims:ime,-2:              0,jms:jme), intent(inout) :: snicexy  !snow layer ice [mm]
      real(r8)  , dimension(ims:ime,-2:              0,jms:jme), intent(inout) :: snliqxy  !snow layer liquid water [mm]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: tvxy        !vegetation canopy temperature
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: tgxy        !ground surface temperature
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: canicexy    !canopy-intercepted ice (mm)
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: canliqxy    !canopy-intercepted liquid water (mm)
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: eahxy       !canopy air vapor pressure (pa)
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: tahxy       !canopy air temperature (k)
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: cmxy        !momentum drag coefficient
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: chxy        !sensible heat exchange coefficient
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: fwetxy      !wetted or snowed fraction of the canopy (-)
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: sneqvoxy    !snow mass at last time step(mm h2o)
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: alboldxy    !snow albedo at last time step (-)
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: qsnowxy     !snowfall on the ground [mm/s]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: wslakexy    !lake water storage [mm]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: zwtxy       !water table depth [m]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: waxy        !water in the "aquifer" [mm]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: wtxy        !groundwater storage [mm]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: lfmassxy    !leaf mass [g/m2]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: rtmassxy    !mass of fine roots [g/m2]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: stmassxy    !stem mass [g/m2]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: woodxy      !mass of wood (incl. woody roots) [g/m2]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: grainxy     !mass of grain [g/m2] !xing
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: gddxy       !growing degree days !xing
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: stblcpxy    !stable carbon in deep soil [g/m2]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: fastcpxy    !short-lived carbon, shallow soil [g/m2]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: xsaixy      !stem area index
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: lai         !leaf area index

      ! iopt_run = 5 option

      real(r8)  , dimension(ims:ime,1:nsoil,jms:jme), intent(inout) , optional :: smoiseq !equilibrium soil moisture content [m3m-3]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) , optional :: smcwtdxy    !deep soil moisture content [m3m-3]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) , optional :: deeprechxy  !deep recharge [m]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) , optional :: rechxy      !accumulated recharge [mm]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) , optional :: qrfsxy      !accumulated flux from groundwater to rivers [mm]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) , optional :: qspringsxy  !accumulated seeping water [mm]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) , optional :: qslatxy     !accumulated lateral flow [mm]
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) , optional :: areaxy      !grid cell area [m2]
      real(r8)  , dimension(ims:ime,jms:jme), intent(in) , optional :: fdepthxy    !efolding depth for transmissivity (m)
      real(r8)  , dimension(ims:ime,jms:jme), intent(in) , optional :: ht          !terrain height (m)
      real(r8)  , dimension(ims:ime,jms:jme), intent(in) , optional :: riverbedxy  !riverbed depth (m)
      real(r8)  , dimension(ims:ime,jms:jme), intent(in) , optional :: eqzwt       !equilibrium water table depth (m)
      real(r8)  , dimension(ims:ime,jms:jme), intent(in) , optional :: rivercondxy !river conductance
      real(r8)  , dimension(ims:ime,jms:jme), intent(in) , optional :: pexpxy      !factor for river conductance

      integer,  intent(out) , optional :: stepwtd
      real(r8)  , intent(in) , optional :: dt, wtddt

      !jref:start
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: t2mvxy        !2m temperature vegetation part (k)
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: t2mbxy        !2m temperature bare ground part (k)
      real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: chstarxy        !dummy
      !jref:end


      real(r8)  , dimension(1:nsoil)  :: zsoil      ! depth of the soil layer bottom (m) from 
      !                                                   the surface (negative)

      real(r8)            :: bexp, smcmax, psisat
      real(r8)            :: fk, masslai,masssai

      real(r8)  , parameter           :: blim  = 5.5
      real(r8)  , parameter           :: hlice = 3.335e5
      real(r8)  , parameter           :: grav = 9.81
      real(r8)  , parameter           :: t0 = 273.15

      integer                   :: errflag, i,j,itf,jtf,ns

      character(len=240) :: err_message
      character(len=4)  :: mminsl
      character(len=*), intent(in) :: mminlu
      !    !character(len=*), intent(in) :: mminlu  !for simple test
   !    character(len=4)        :: mminlu

   !    mminlu='igbp'
   !    mminsl='stas'

      mminsl='stas'

      call read_mp_veg_parameters(trim(mminlu))
      call read_mp_soil_parameters()
      call read_mp_rad_parameters()
      call read_mp_global_parameters()
      call read_mp_crop_parameters()

      if( .not. restart ) then

         itf=min0(ite,ide-1)
         jtf=min0(jte,jde-1)

         !
         ! initialize physical snow height snowh
         !
         if(.not.fndsnowh)then
            ! if no snowh do the following
            call wrf_message( 'snow height not found - value defined in lsminit' )
            do j = jts,jtf
               do i = its,itf
                  snowh(i,j)=snow(i,j)*0.005               ! snow in mm and snowh in m
               enddo
            enddo
         endif


         ! check if snow/snowh are consistent and cap swe at 2000mm;
         !  the noah-mp code does it internally but if we don't do it here, problems ensue
         do j = jts,jtf
            do i = its,itf
               if ( snow(i,j) > 0. .and. snowh(i,j) == 0. .or. snowh(i,j) > 0. .and. snow(i,j) == 0.) then
               write(err_message,*)"problem with initial snow fields: snow/snowh>0 while snowh/snow=0 at i,j" &
                                       ,i,j,snow(i,j),snowh(i,j)
               call wrf_message(err_message)
               endif
               if ( snow( i,j ) > 2000. ) then
               snowh(i,j) = snowh(i,j) * 2000. / snow(i,j)      ! snow in mm and snowh in m
               snow (i,j) = 2000.                               ! cap snow at 2000, maintain density
               endif
            enddo
         enddo

         errflag = 0
         do j = jts,jtf
            do i = its,itf
               if ( isltyp( i,j ) .lt. 1 ) then
                  errflag = 1
                  write(err_message,*)"module_sf_noahlsm.f: lsminit: out of range isltyp ",i,j,isltyp( i,j )
                  call wrf_message(err_message)
               endif
            enddo
         enddo
         if ( errflag .eq. 1 ) then
            call wrf_error_fatal( "module_sf_noahlsm.f: lsminit: out of range value "// &
               "of isltyp. is this field in the input?" )
         endif
      ! gac-->lateralflow
      ! 20130219 - no longer need this - see module_data_gocart_dust
      !#if ( wrf_chem == 1 )
      !       !
      !       ! need this parameter for dust parameterization in wrf/chem
      !       !
      !       do i=1,nsltype
      !          porosity(i)=maxsmc(i)
      !       enddo
      !#endif
      ! <--gac

      ! initialize soil liquid water content sh2o

         do j = jts , jtf
            do i = its , itf
         if(ivgtyp(i,j)==isice_table .and. xice(i,j) <= 0.0) then
               do ns=1, nsoil
            smois(i,ns,j) = 1.0                     ! glacier starts all frozen
            sh2o(i,ns,j) = 0.0
            tslb(i,ns,j) = min(tslb(i,ns,j),263.15) ! set glacier temp to at most -10c
               end do
            !tmn(i,j) = min(tmn(i,j),263.15)         ! set deep temp to at most -10c
      snow(i,j) = max(snow(i,j), 10.0)        ! set swe to at least 10mm
                  snowh(i,j)=snow(i,j)*0.01               ! snow in mm and snowh in m
         else
         
               bexp   =   bexp_table(isltyp(i,j))
               smcmax = smcmax_table(isltyp(i,j))
               psisat = psisat_table(isltyp(i,j))

               do ns=1, nsoil
            if ( smois(i,ns,j) > smcmax )  smois(i,ns,j) = smcmax
               end do
               if ( ( bexp > 0.0 ) .and. ( smcmax > 0.0 ) .and. ( psisat > 0.0 ) ) then
                  do ns=1, nsoil
                     if ( tslb(i,ns,j) < 273.149 ) then    ! use explicit as initial soil ice
                        fk=(( (hlice/(grav*(-psisat))) *                              &
                           ((tslb(i,ns,j)-t0)/tslb(i,ns,j)) )**(-1/bexp) )*smcmax
                        fk = max(fk, 0.02)
                        sh2o(i,ns,j) = min( fk, smois(i,ns,j) )
                     else
                        sh2o(i,ns,j)=smois(i,ns,j)
                     endif
                  end do
               else
                  do ns=1, nsoil
                     sh2o(i,ns,j)=smois(i,ns,j)
                  end do
               endif
            endif
            enddo
         enddo
      !  endif


         do j = jts,jtf
            do i = its,itf
               tvxy       (i,j) = tsk(i,j)
            if(snow(i,j) > 0.0 .and. tsk(i,j) > 273.15) tvxy(i,j) = 273.15
               tgxy       (i,j) = tsk(i,j)
            if(snow(i,j) > 0.0 .and. tsk(i,j) > 273.15) tgxy(i,j) = 273.15
               canwat     (i,j) = 0.0
               canliqxy   (i,j) = canwat(i,j)
               canicexy   (i,j) = 0.
               eahxy      (i,j) = 2000. 
               tahxy      (i,j) = tsk(i,j)
            if(snow(i,j) > 0.0 .and. tsk(i,j) > 273.15) tahxy(i,j) = 273.15
      !             tahxy      (i,j) = 287.
      !jref:start
               t2mvxy     (i,j) = tsk(i,j)
            if(snow(i,j) > 0.0 .and. tsk(i,j) > 273.15) t2mvxy(i,j) = 273.15
               t2mbxy     (i,j) = tsk(i,j)
            if(snow(i,j) > 0.0 .and. tsk(i,j) > 273.15) t2mbxy(i,j) = 273.15
               chstarxy     (i,j) = 0.1
      !jref:end

               cmxy       (i,j) = 0.0
               chxy       (i,j) = 0.0
               fwetxy     (i,j) = 0.0
               sneqvoxy   (i,j) = 0.0
               alboldxy   (i,j) = 0.65
               qsnowxy    (i,j) = 0.0
               wslakexy   (i,j) = 0.0

               if(iopt_run.ne.5) then 
                     waxy       (i,j) = 4900.                                       !???
                     wtxy       (i,j) = waxy(i,j)                                   !???
                     zwtxy      (i,j) = (25. + 2.0) - waxy(i,j)/1000/0.2            !???
               else
                     waxy       (i,j) = 0.
                     wtxy       (i,j) = 0.
                     areaxy     (i,j) = (dx * dy) / ( msftx(i,j) * msfty(i,j) )
               endif

            if(ivgtyp(i,j) == isbarren_table .or. ivgtyp(i,j) == isice_table .or. &
         ivgtyp(i,j) == isurban_table  .or. ivgtyp(i,j) == iswater_table ) then
         
         lai        (i,j) = 0.0
               xsaixy     (i,j) = 0.0
               lfmassxy   (i,j) = 0.0
               stmassxy   (i,j) = 0.0
               rtmassxy   (i,j) = 0.0
               woodxy     (i,j) = 0.0
               stblcpxy   (i,j) = 0.0
               fastcpxy   (i,j) = 0.0

      else
         
         lai        (i,j) = max(lai(i,j),0.05)             ! at least start with 0.05 for arbitrary initialization (v3.7)
               xsaixy     (i,j) = max(0.1*lai(i,j),0.05)         ! mb: arbitrarily initialize sai using input lai (v3.7)
               masslai = 1000. / max(sla_table(ivgtyp(i,j)),1.0) ! conversion from lai to mass  (v3.7)
               lfmassxy   (i,j) = lai(i,j)*masslai               ! use lai to initialize (v3.7)
               masssai = 1000. / 3.0                             ! conversion from lai to mass (v3.7)
               stmassxy   (i,j) = xsaixy(i,j)*masssai            ! use sai to initialize (v3.7)
               rtmassxy   (i,j) = 500.0                          ! these are all arbitrary and probably should be
               woodxy     (i,j) = 500.0                          ! in the table or read from initialization
               stblcpxy   (i,j) = 1000.0                         !
               fastcpxy   (i,j) = 1000.0                         !
               grainxy    (i,j) = 1e-10         ! add by xing
               gddxy      (i,j) = 0          ! add by xing

      end if

            enddo
         enddo

         ! given the soil layer thicknesses (in dzs), initialize the soil layer
         ! depths from the surface.
         zsoil(1)         = -dzs(1)          ! negative
         do ns=2, nsoil
            zsoil(ns)       = zsoil(ns-1) - dzs(ns)
         end do

         ! initialize snow/soil layer arrays zsnsoxy, tsnoxy, snicexy, snliqxy, 
         ! and isnowxy
         call snow_init ( ims , ime , jms , jme , its , itf , jts , jtf , 3 , &
            &           nsoil , zsoil , snow , tgxy , snowh ,     &
            &           zsnsoxy , tsnoxy , snicexy , snliqxy , isnowxy )

         !initialize arrays for groundwater dynamics iopt_run=5

         if(iopt_run.eq.5) then
            if ( present(smoiseq)     .and. &
            present(smcwtdxy)    .and. &
            present(rechxy)      .and. &
            present(deeprechxy)  .and. &
            present(areaxy)      .and. &
            present(dx)          .and. &
            present(dy)          .and. &
            present(msftx)       .and. &
            present(msfty)       .and. &
            present(wtddt)       .and. &
            present(stepwtd)     .and. &
            present(dt)          .and. &
            present(qrfsxy)      .and. &
            present(qspringsxy)  .and. &
            present(qslatxy)     .and. &
            present(fdepthxy)    .and. &
            present(ht)          .and. &
            present(riverbedxy)  .and. &
            present(eqzwt)       .and. &
            present(rivercondxy) .and. &
            present(pexpxy)            ) then

               stepwtd = nint(wtddt*60./dt)
               stepwtd = max(stepwtd,1)

               call groundwater_init ( & 
      &       nsoil, zsoil , dzs  ,isltyp, ivgtyp,wtddt , &
      &       fdepthxy, ht, riverbedxy, eqzwt, rivercondxy, pexpxy , areaxy, zwtxy,   &
      &       smois,sh2o, smoiseq, smcwtdxy, deeprechxy, rechxy, qslatxy, qrfsxy, qspringsxy, &
      &       ids,ide, jds,jde, kds,kde,                    &
      &       ims,ime, jms,jme, kms,kme,                    &
      &       its,ite, jts,jte, kts,kte                     )

            else
               call wrf_error_fatal ('not enough fields to use groundwater option in noah-mp')
            end if
         endif

      endif
   end subroutine noahmp_init

   !------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

  subroutine snow_init ( ims , ime , jms , jme , its , itf , jts , jtf ,                  &
   &                 nsnow , nsoil , zsoil , swe , tgxy , snodep ,                    &
   &                 zsnsoxy , tsnoxy , snicexy ,snliqxy , isnowxy )
!------------------------------------------------------------------------------------------
!   initialize snow arrays for noah-mp lsm, based in input snowdep, nsnow
!   isnowxy is an index array, indicating the index of the top snow layer.  valid indices
!           for snow layers range from 0 (no snow) and -1 (shallow snow) to (-nsnow)+1 (deep snow).
!   tsnoxy holds the temperature of the snow layer.  snow layers are initialized with 
!          temperature = ground temperature [?].  snow-free levels in the array have value 0.0
!   snicexy is the frozen content of a snow layer.  initial estimate based on snodep and swe
!   snliqxy is the liquid content of a snow layer.  initialized to 0.0
!   znsnoxy is the layer depth from the surface.  
!------------------------------------------------------------------------------------------
implicit none
!------------------------------------------------------------------------------------------
integer, intent(in)                              :: ims, ime, jms, jme
integer, intent(in)                              :: its, itf, jts, jtf
integer, intent(in)                              :: nsnow
integer, intent(in)                              :: nsoil
real(r8)  ,    intent(in), dimension(ims:ime, jms:jme) :: swe 
real(r8)  ,    intent(in), dimension(ims:ime, jms:jme) :: snodep
real(r8)  ,    intent(in), dimension(ims:ime, jms:jme) :: tgxy
real(r8)  ,    intent(in), dimension(1:nsoil)          :: zsoil

integer, intent(out), dimension(ims:ime, jms:jme)                :: isnowxy ! top snow layer index
real(r8)  ,    intent(out), dimension(ims:ime, -nsnow+1:nsoil,jms:jme) :: zsnsoxy ! snow/soil layer depth from surface [m]
real(r8)  ,    intent(out), dimension(ims:ime, -nsnow+1:    0,jms:jme) :: tsnoxy  ! snow layer temperature [k]
real(r8)  ,    intent(out), dimension(ims:ime, -nsnow+1:    0,jms:jme) :: snicexy ! snow layer ice content [mm]
real(r8)  ,    intent(out), dimension(ims:ime, -nsnow+1:    0,jms:jme) :: snliqxy ! snow layer liquid content [mm]

! local variables:
!   dzsno   holds the thicknesses of the various snow layers.
!   dzsnoso holds the thicknesses of the various soil/snow layers.
integer                           :: i,j,iz
real(r8)  ,   dimension(-nsnow+1:    0) :: dzsno
real(r8)  ,   dimension(-nsnow+1:nsoil) :: dzsnso

!------------------------------------------------------------------------------------------

do j = jts , jtf
   do i = its , itf
      if ( snodep(i,j) < 0.025 ) then
         isnowxy(i,j) = 0
         dzsno(-nsnow+1:0) = 0.
      else
         if ( ( snodep(i,j) >= 0.025 ) .and. ( snodep(i,j) <= 0.05 ) ) then
            isnowxy(i,j)    = -1
            dzsno(0)  = snodep(i,j)
         else if ( ( snodep(i,j) > 0.05 ) .and. ( snodep(i,j) <= 0.10 ) ) then
            isnowxy(i,j)    = -2
            dzsno(-1) = snodep(i,j)/2.
            dzsno( 0) = snodep(i,j)/2.
         else if ( (snodep(i,j) > 0.10 ) .and. ( snodep(i,j) <= 0.25 ) ) then
            isnowxy(i,j)    = -2
            dzsno(-1) = 0.05
            dzsno( 0) = snodep(i,j) - dzsno(-1)
         else if ( ( snodep(i,j) > 0.25 ) .and. ( snodep(i,j) <= 0.45 ) ) then
            isnowxy(i,j)    = -3
            dzsno(-2) = 0.05
            dzsno(-1) = 0.5*(snodep(i,j)-dzsno(-2))
            dzsno( 0) = 0.5*(snodep(i,j)-dzsno(-2))
         else if ( snodep(i,j) > 0.45 ) then
            isnowxy(i,j)     = -3
            dzsno(-2) = 0.05
            dzsno(-1) = 0.20
            dzsno( 0) = snodep(i,j) - dzsno(-1) - dzsno(-2)
         else
            call wrf_error_fatal("problem with the logic assigning snow layers.")
         end if
      end if

      tsnoxy (i,-nsnow+1:0,j) = 0.
      snicexy(i,-nsnow+1:0,j) = 0.
      snliqxy(i,-nsnow+1:0,j) = 0.
      do iz = isnowxy(i,j)+1 , 0
         tsnoxy(i,iz,j)  = tgxy(i,j)  ! [k]
         snliqxy(i,iz,j) = 0.00
         snicexy(i,iz,j) = 1.00 * dzsno(iz) * (swe(i,j)/snodep(i,j))  ! [kg/m3]
      end do

      ! assign local variable dzsnso, the soil/snow layer thicknesses, for snow layers
      do iz = isnowxy(i,j)+1 , 0
         dzsnso(iz) = -dzsno(iz)
      end do

      ! assign local variable dzsnso, the soil/snow layer thicknesses, for soil layers
      dzsnso(1) = zsoil(1)
      do iz = 2 , nsoil
         dzsnso(iz) = (zsoil(iz) - zsoil(iz-1))
      end do

      ! assign zsnsoxy, the layer depths, for soil and snow layers
      zsnsoxy(i,isnowxy(i,j)+1,j) = dzsnso(isnowxy(i,j)+1)
      do iz = isnowxy(i,j)+2 , nsoil
         zsnsoxy(i,iz,j) = zsnsoxy(i,iz-1,j) + dzsnso(iz)
      enddo

   end do
end do

end subroutine snow_init
! ==================================================================================================
! ----------------------------------------------------------------------
subroutine groundwater_init (   &
        &            nsoil , zsoil , dzs, isltyp, ivgtyp, wtddt , &
        &            fdepth, topo, riverbed, eqwtd, rivercond, pexp , area ,wtd ,  &
        &            smois,sh2o, smoiseq, smcwtdxy, deeprechxy, rechxy ,  &
        &            qslatxy, qrfsxy, qspringsxy,                  &
        &            ids,ide, jds,jde, kds,kde,                    &
        &            ims,ime, jms,jme, kms,kme,                    &
        &            its,ite, jts,jte, kts,kte                     )


use noahmp_tables, only : bexp_table,smcmax_table,psisat_table,smcwlt_table,dwsat_table,dksat_table, &
                            isurban_table, isice_table ,iswater_table
use module_sf_noahmp_groundwater, only : lateralflow

! ----------------------------------------------------------------------
implicit none
! ----------------------------------------------------------------------

integer,  intent(in   )   ::     ids,ide, jds,jde, kds,kde,  &
     &                           ims,ime, jms,jme, kms,kme,  &
     &                           its,ite, jts,jte, kts,kte
integer, intent(in)                              :: nsoil
real(r8)  ,   intent(in)                               ::     wtddt
real(r8)  ,    intent(in), dimension(1:nsoil)          :: zsoil,dzs
integer, intent(in), dimension(ims:ime, jms:jme) :: isltyp, ivgtyp
real(r8)  ,    intent(in), dimension(ims:ime, jms:jme) :: fdepth, topo, riverbed, eqwtd, rivercond, pexp , area
real(r8)  ,    intent(inout), dimension(ims:ime, jms:jme) :: wtd
real(r8)  ,     dimension( ims:ime , 1:nsoil, jms:jme ), &
     &    intent(inout)   ::                          smois, &
     &                                                 sh2o, &
     &                                                 smoiseq
real(r8)  ,    intent(inout), dimension(ims:ime, jms:jme) ::  &
                                                       smcwtdxy, &
                                                       deeprechxy, &
                                                       rechxy, &
                                                       qslatxy, &
                                                       qrfsxy, &
                                                       qspringsxy  
! local
integer  :: i,j,k,iter,itf,jtf
real(r8)            :: bexp,smcmax,psisat,smcwlt,dwsat,dksat
real(r8)            :: frliq,smceqdeep
real(r8)            :: deltat,rcond
real(r8)            :: aa,bbb,cc,dd,dx,func,dfunc,ddz,expon,smc,flux
real(r8)  , dimension(1:nsoil) :: smceq
real(r8)  ,      dimension( ims:ime, jms:jme )    :: qlat, qrf
integer,   dimension( ims:ime, jms:jme )    :: landmask !-1 for water (ice or no ice) and glacial areas, 1 for land where the lsm does its soil moisture calculations


   itf=min0(ite,ide-1)
   jtf=min0(jte,jde-1)

!first compute lateral flow and flow to rivers to initialize deep soil moisture

deltat = wtddt * 60. !timestep in seconds for this calculation

where(ivgtyp.ne.iswater_table.and.ivgtyp.ne.isice_table)
     landmask=1
elsewhere
     landmask=-1
endwhere

!calculate lateral flow

qlat = 0.
call lateralflow(isltyp,wtd,qlat,fdepth,topo,landmask,deltat,area       &
                    ,ids,ide,jds,jde,kds,kde                      & 
                    ,ims,ime,jms,jme,kms,kme                      &
                    ,its,ite,jts,jte,kts,kte                      )
                    

!compute flux from grounwater to rivers in the cell

do j=jts,jtf
   do i=its,itf
      if(landmask(i,j).gt.0)then
         if(wtd(i,j) .gt. riverbed(i,j) .and.  eqwtd(i,j) .gt. riverbed(i,j)) then
           rcond = rivercond(i,j) * exp(pexp(i,j)*(wtd(i,j)-eqwtd(i,j)))
         else    
           rcond = rivercond(i,j)
         endif
         qrf(i,j) = rcond * (wtd(i,j)-riverbed(i,j)) * deltat/area(i,j)
!for now, dont allow it to go from river to groundwater
         qrf(i,j) = max(qrf(i,j),0.) 
      else
         qrf(i,j) = 0.
      endif
   enddo
enddo

!now compute eq. soil moisture, change soil moisture to be compatible with the water table and compute deep soil moisture

   do j = jts,jtf
      do i = its,itf
         bexp   =   bexp_table(isltyp(i,j))
         smcmax = smcmax_table(isltyp(i,j))
         smcwlt = smcwlt_table(isltyp(i,j))
         if(ivgtyp(i,j)==isurban_table)then
             smcmax = 0.45         
             smcwlt = 0.40         
         endif 
         dwsat  =   dwsat_table(isltyp(i,j))
         dksat  =   dksat_table(isltyp(i,j))
         psisat = -psisat_table(isltyp(i,j))
       if ( ( bexp > 0.0 ) .and. ( smcmax > 0.0 ) .and. ( -psisat > 0.0 ) ) then
         !initialize equilibrium soil moisture for water table diagnostic
                call eqsmoisture(nsoil ,  zsoil , smcmax , smcwlt ,dwsat, dksat  ,bexp  , & !in
                                 smceq                          )  !out

         smoiseq (i,1:nsoil,j) = smceq (1:nsoil)


          !make sure that below the water table the layers are saturated and initialize the deep soil moisture
         if(wtd(i,j) < zsoil(nsoil)-dzs(nsoil)) then

!initialize deep soil moisture so that the flux compensates qlat+qrf
!use newton-raphson method to find soil moisture

                     expon = 2. * bexp + 3.
                     ddz = zsoil(nsoil) - wtd(i,j)
                     cc = psisat/ddz
                     flux = (qlat(i,j)-qrf(i,j))/deltat

                     smc = 0.5 * smcmax

                     do iter = 1, 100
                       dd = (smc+smcmax)/(2.*smcmax)
                       aa = -dksat * dd  ** expon
                       bbb = cc * ( (smcmax/smc)**bexp - 1. ) + 1. 
                       func =  aa * bbb - flux
                       dfunc = -dksat * (expon/(2.*smcmax)) * dd ** (expon - 1.) * bbb &
                               + aa * cc * (-bexp) * smcmax ** bexp * smc ** (-bexp-1.)

                       dx = func/dfunc
                       smc = smc - dx
                       if ( abs (dx) < 1.e-6)exit
                     enddo

              smcwtdxy(i,j) = max(smc,1.e-4)

         elseif(wtd(i,j) < zsoil(nsoil))then
              smceqdeep = smcmax * ( psisat / ( psisat - dzs(nsoil) ) ) ** (1./bexp)
!                  smceqdeep = max(smceqdeep,smcwlt)
              smceqdeep = max(smceqdeep,1.e-4)
              smcwtdxy(i,j) = smcmax * ( wtd(i,j) -  (zsoil(nsoil)-dzs(nsoil))) + &
                              smceqdeep * (zsoil(nsoil) - wtd(i,j))

         else !water table within the resolved layers
              smcwtdxy(i,j) = smcmax
              do k=nsoil,2,-1
                 if(wtd(i,j) .ge. zsoil(k-1))then
                      frliq = sh2o(i,k,j) / smois(i,k,j)
                      smois(i,k,j) = smcmax
                      sh2o(i,k,j) = smcmax * frliq
                 else
                      if(smois(i,k,j).lt.smceq(k))then
                          wtd(i,j) = zsoil(k)
                      else
                          wtd(i,j) = ( smois(i,k,j)*dzs(k) - smceq(k)*zsoil(k-1) + smcmax*zsoil(k) ) / &
                                     (smcmax - smceq(k))   
                      endif
                      exit
                 endif
              enddo
         endif
        else
          smoiseq (i,1:nsoil,j) = smcmax
          smcwtdxy(i,j) = smcmax
          wtd(i,j) = 0.
        endif

!zero out some arrays

         deeprechxy(i,j) = 0.
         rechxy(i,j) = 0.
         qslatxy(i,j) = 0.
         qrfsxy(i,j) = 0.
         qspringsxy(i,j) = 0.

      enddo
   enddo




end  subroutine groundwater_init
! ==================================================================================================
! ----------------------------------------------------------------------
subroutine eqsmoisture(nsoil  ,  zsoil , smcmax , smcwlt, dwsat , dksat ,bexp , & !in
                     smceq                          )  !out
! ----------------------------------------------------------------------
implicit none
! ----------------------------------------------------------------------
! input
integer,                         intent(in) :: nsoil !no. of soil layers
real(r8)  , dimension(       1:nsoil), intent(in) :: zsoil !depth of soil layer-bottom [m]
real(r8)  ,                            intent(in) :: smcmax , smcwlt, bexp , dwsat, dksat
!output
real(r8)  ,  dimension(      1:nsoil), intent(out) :: smceq  !equilibrium soil water  content [m3/m3]
!local
integer                                     :: k , iter
real(r8)            :: ddz , smc, func, dfunc , aa, bb , expon, dx

!gmmcompute equilibrium soil moisture content for the layer when wtd=zsoil(k)


do k=1,nsoil

        if ( k == 1 )then
            ddz = -zsoil(k+1) * 0.5
        elseif ( k < nsoil ) then
            ddz = ( zsoil(k-1) - zsoil(k+1) ) * 0.5
        else
            ddz = zsoil(k-1) - zsoil(k)
        endif

!use newton-raphson method to find eq soil moisture

        expon = bexp +1.
        aa = dwsat/ddz
        bb = dksat / smcmax ** expon

        smc = 0.5 * smcmax

     do iter = 1, 100
        func = (smc - smcmax) * aa +  bb * smc ** expon
        dfunc = aa + bb * expon * smc ** bexp 

        dx = func/dfunc
        smc = smc - dx
        if ( abs (dx) < 1.e-6)exit
     enddo

!             smceq(k) = min(max(smc,smcwlt),smcmax*0.99)
         smceq(k) = min(max(smc,1.e-4),smcmax*0.99)
enddo

end  subroutine eqsmoisture

   ! !------------------------------------------------------------------------------------------

   subroutine wrf_message( str )
   implicit none

   character*(*) str
   if(mpi_rank()==0)then
   print *,trim(str)
   end if

   end subroutine wrf_message
   subroutine wrf_error_fatal( str )
   implicit none
   character*(*) str
   call wrf_message( '-------------- fatal called ---------------' )
   call wrf_message( str )
   call wrf_message( '-------------------------------------------' )

   call wrf_abort
   end subroutine wrf_error_fatal

   subroutine wrf_abort
         stop 'wrf_abort'
   end subroutine wrf_abort

end module module_sf_noahmp_init
