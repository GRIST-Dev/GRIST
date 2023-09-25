
module noahmp_glacier_globals
   use grist_constants, only: r8, zero

  implicit none

  real(r8),   parameter :: grav   = 9.80616   
  real(r8),   parameter :: sb     = 5.67e-08  
  real(r8),   parameter :: vkc    = 0.40      
  real(r8),   parameter :: tfrz   = 273.16    
  real(r8),   parameter :: hsub   = 2.8440e06 
  real(r8),   parameter :: hvap   = 2.5104e06 
  real(r8),   parameter :: hfus   = 0.3336e06 
  real(r8),   parameter :: cwat   = 4.188e06  
  real(r8),   parameter :: cice   = 2.094e06  
  real(r8),   parameter :: cpair  = 1004.64   
  real(r8),   parameter :: tkwat  = 0.6       
  real(r8),   parameter :: tkice  = 2.2       
  real(r8),   parameter :: tkair  = 0.023     
  real(r8),   parameter :: rair   = 287.04    
  real(r8),   parameter :: rw     = 461.269   
  real(r8),   parameter :: denh2o = 1000.     
  real(r8),   parameter :: denice = 917.      
  integer :: opt_alb 
  integer :: opt_snf 
  integer :: opt_tbot 
  integer :: opt_stc 
  integer :: opt_gla
  real(r8),   parameter :: z0sno  = 0.002  
  real(r8),   parameter :: ssi    = 0.03   
  real(r8),   parameter :: swemx  = 1.00   
                                     

end module noahmp_glacier_globals


module noahmp_glacier_routines
  use noahmp_glacier_globals
  use grist_constants, only: r8
  implicit none

  public  :: noahmp_options_glacier
  public  :: noahmp_glacier

  private :: atm_glacier
  private :: energy_glacier
  private ::       thermoprop_glacier
  private ::               csnow_glacier
  private ::       radiation_glacier
  private ::               snow_age_glacier
  private ::               snowalb_bats_glacier  
  private ::               snowalb_class_glacier
  private ::       glacier_flux
  private ::               sfcdif1_glacier                  
  private ::       tsnosoi_glacier
  private ::               hrt_glacier
  private ::               hstep_glacier   
  private ::                         rosr12_glacier
  private ::       phasechange_glacier

  private :: water_glacier
  private ::       snowwater_glacier
  private ::               snowfall_glacier
  private ::               combine_glacier
  private ::               divide_glacier
  private ::                         combo_glacier
  private ::               compact_glacier
  private ::               snowh2o_glacier

  private :: error_glacier

contains



  subroutine noahmp_glacier (&
                   iloc    ,jloc    ,cosz    ,nsnow   ,nsoil   ,dt      , & 
                   sfctmp  ,sfcprs  ,uu      ,vv      ,q2      ,soldn   , & 
                   prcp    ,lwdn    ,tbot    ,zlvl    ,ficeold ,zsoil   , & 
                   qsnow   ,sneqvo  ,albold  ,cm      ,ch      ,isnow   , & 
                   sneqv   ,smc     ,zsnso   ,snowh   ,snice   ,snliq   , & 
                   tg      ,stc     ,sh2o    ,tauss   ,qsfc    ,          & 
                   fsa     ,fsr     ,fira    ,fsh     ,fgev    ,ssoil   , & 
                   trad    ,edir    ,runsrf  ,runsub  ,sag     ,albedo  , & 
                   qsnbot  ,ponding ,ponding1,ponding2,t2m     ,q2e     , & 
                   emissi,  fpice,    ch2b,                                & 
                   taux, tauy,fire, albd_out,albi_out)  !--cheyz





  implicit none


  integer                        , intent(in)    :: iloc   
  integer                        , intent(in)    :: jloc   
  real(r8)                   , intent(in)    :: cosz   
  integer                        , intent(in)    :: nsnow  
  integer                        , intent(in)    :: nsoil  
  real(r8)                   , intent(in)    :: dt     
  real(r8)                   , intent(in)    :: sfctmp 
  real(r8)                   , intent(in)    :: sfcprs 
  real(r8)                   , intent(in)    :: uu     
  real(r8)                   , intent(in)    :: vv     
  real(r8)                   , intent(in)    :: q2     
  real(r8)                   , intent(in)    :: soldn  
  real(r8)                   , intent(in)    :: prcp   
  real(r8)                   , intent(in)    :: lwdn   
  real(r8)                   , intent(in)    :: tbot   
  real(r8)                   , intent(in)    :: zlvl   
  real(r8),   dimension(-nsnow+1:    0), intent(in)    :: ficeold
  real(r8),   dimension(       1:nsoil), intent(in)    :: zsoil  



  real(r8)                   , intent(inout) :: qsnow  
  real(r8)                   , intent(inout) :: sneqvo 
  real(r8)                   , intent(inout) :: albold 
  real(r8)                   , intent(inout) :: cm     
  real(r8)                   , intent(inout) :: ch     


  integer                        , intent(inout) :: isnow  
  real(r8)                   , intent(inout) :: sneqv  
  real(r8),   dimension(       1:nsoil), intent(inout) :: smc    
  real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: zsnso  
  real(r8)                   , intent(inout) :: snowh  
  real(r8),   dimension(-nsnow+1:    0), intent(inout) :: snice  
  real(r8),   dimension(-nsnow+1:    0), intent(inout) :: snliq  
  real(r8)                   , intent(inout) :: tg     
  real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: stc    
  real(r8),   dimension(       1:nsoil), intent(inout) :: sh2o   
  real(r8)                   , intent(inout) :: tauss  
  real(r8)                   , intent(inout) :: qsfc   


  real(r8)                   , intent(out)   :: fsa    
  real(r8)                   , intent(out)   :: fsr    
  real(r8)                   , intent(out)   :: fira   
  real(r8)                   , intent(out)   :: fsh    
  real(r8)                   , intent(out)   :: fgev   
  real(r8)                   , intent(out)   :: ssoil  
  real(r8)                   , intent(out)   :: trad   
  real(r8)                   , intent(out)   :: edir   
  real(r8)                   , intent(out)   :: runsrf 
  real(r8)                   , intent(out)   :: runsub 
  real(r8)                   , intent(out)   :: sag    
  real(r8)                   , intent(out)   :: albedo 
  real(r8)                   , intent(out)   :: qsnbot 
  real(r8)                   , intent(out)   :: ponding
  real(r8)                   , intent(out)   :: ponding1
  real(r8)                   , intent(out)   :: ponding2
  real(r8)                   , intent(out)   :: t2m     
  real(r8)                   , intent(out)   :: q2e
  real(r8)                   , intent(out)   :: emissi
  real(r8)                   , intent(out)   :: fpice
  real(r8)                   , intent(out)   :: ch2b
  real(r8)                   , intent(out)   :: taux, tauy, fire !--cheyz 20200427
  real(r8),dimension(1:2)     , intent(out)     :: albd_out,albi_out !--cheyz

  integer                                        :: iz, ib     
  integer, dimension(-nsnow+1:nsoil)             :: imelt  
  real(r8)  :: rhoair 
  real(r8),   dimension(-nsnow+1:nsoil)                :: dzsnso 
  real(r8)  :: thair  
  real(r8)  :: qair   
  real(r8)  :: eair   
  real(r8),   dimension(       1:    2)                :: solad  
  real(r8),   dimension(       1:    2)                :: solai  
  real(r8),   dimension(       1:nsoil)                :: sice   
  real(r8),   dimension(-nsnow+1:    0)                :: snicev 
  real(r8),   dimension(-nsnow+1:    0)                :: snliqv 
  real(r8),   dimension(-nsnow+1:    0)                :: epore  
  real(r8)  :: qdew   
  real(r8)  :: qvap   
  real(r8)  :: lathea 
  real(r8)  :: qmelt  
  real(r8)  :: swdown 
  real(r8)  :: beg_wb 
  real(r8)  :: zbot = -8.0 

  character*256 message




   call atm_glacier (sfcprs ,sfctmp ,q2     ,soldn  ,cosz   ,thair  , & 
                     qair   ,eair   ,rhoair ,solad  ,solai  ,swdown )

   beg_wb = sneqv



     do iz = isnow+1, nsoil
         if(iz == isnow+1) then
           dzsnso(iz) = - zsnso(iz)
         else
           dzsnso(iz) = zsnso(iz-1) - zsnso(iz)
         end if
     end do



    call energy_glacier (nsnow  ,nsoil  ,isnow  ,dt     ,qsnow  ,rhoair , & 
                         eair   ,sfcprs ,qair   ,sfctmp ,lwdn   ,uu     , & 
                         vv     ,solad  ,solai  ,cosz   ,zlvl   ,         & 
                         tbot   ,zbot   ,zsnso  ,dzsnso ,                 & 
                         tg     ,stc    ,snowh  ,sneqv  ,sneqvo ,sh2o   , & 
                         smc    ,snice  ,snliq  ,albold ,cm     ,ch     , & 
                         tauss  ,qsfc   ,                                 & 
                         imelt  ,snicev ,snliqv ,epore  ,qmelt  ,ponding, & 
		         sag    ,fsa    ,fsr    ,fira   ,fsh    ,fgev   , & 
               trad   ,t2m    ,ssoil  ,lathea ,q2e    ,emissi, ch2b,  &
               taux, tauy, fire ,albd_out,albi_out  )    !--cheyz 20200427

    sice = max(0.0, smc - sh2o)   
    sneqvo  = sneqv

    qvap = max( fgev/lathea, 0.)       
    qdew = abs( min(fgev/lathea, 0.))  
    edir = qvap - qdew



     call water_glacier (nsnow  ,nsoil  ,imelt  ,dt     ,prcp   ,sfctmp , & 
                         qvap   ,qdew   ,ficeold,zsoil  ,                 & 
                         isnow  ,snowh  ,sneqv  ,snice  ,snliq  ,stc    , & 
                         dzsnso ,sh2o   ,sice   ,ponding,zsnso  ,fsh    , & 
                         runsrf ,runsub ,qsnow  ,ponding1       ,ponding2,qsnbot,fpice &  
                        )

     if(opt_gla == 2) then
       edir = qvap - qdew
       fgev = edir * lathea
     end if

     if(maxval(sice) < 0.0001) then
       write(message,*) "glacier has melted at:",iloc,jloc," are you sure this should be a glacier point?"
      ! call wrf_debug(10,trim(message))
     end if



     call error_glacier (iloc   ,jloc   ,swdown ,fsa    ,fsr    ,fira   , &
                         fsh    ,fgev   ,ssoil  ,sag    ,prcp   ,edir   , &
		         runsrf ,runsub ,sneqv  ,dt     ,beg_wb )

    if(snowh <= 1.e-6 .or. sneqv <= 1.e-3) then
     snowh = 0.0
     sneqv = 0.0
    end if

    if(swdown.ne.0.) then
      albedo = fsr / swdown
    else
      albedo = -999.9
    end if
    
   !--cheyz  2020/04/27 should be tested!
! yizhang 2021/01/26 
#ifndef REGLSM
   do ib=1, 2 ! nband
      !if(swdown .ne.0. )then ! this is not the swdown from radiation, lead to
      !inconsitency, use cosz, and this value should better come from the same
      !one used by other (sea, seaice, radiation)
      ! comment this if or not does not affect bit solutions; no zero for
      ! diagnose
      if(cosz.gt.zero)then
         albd_out(ib)= albd_out(ib) !/swdown
         albi_out(ib)= albi_out(ib) !/swdown
      else
         albd_out(ib)= 0._r8 !-999.9
         albi_out(ib)= 0._r8 !-999.9
      end if
   end do
#else
   do ib=1, 2 ! nband
      if(swdown .ne.0. )then
         albd_out(ib)= albd_out(ib) /swdown
         albi_out(ib)= albi_out(ib) /swdown
      else
         albd_out(ib)= 0._r8 !-999.9
         albi_out(ib)= 0._r8 !-999.9
      end if
   end do
#endif

  end subroutine noahmp_glacier

  subroutine atm_glacier (sfcprs ,sfctmp ,q2     ,soldn  ,cosz   ,thair  , &
                          qair   ,eair   ,rhoair ,solad  ,solai  , &
                          swdown )     



  implicit none



  real(r8)                   , intent(in)  :: sfcprs 
  real(r8)                   , intent(in)  :: sfctmp 
  real(r8)                   , intent(in)  :: q2     
  real(r8)                   , intent(in)  :: soldn  
  real(r8)                   , intent(in)  :: cosz   



  real(r8)                   , intent(out) :: thair  
  real(r8)                   , intent(out) :: qair   
  real(r8)                   , intent(out) :: eair   
  real(r8),   dimension(       1:   2), intent(out) :: solad  
  real(r8),   dimension(       1:   2), intent(out) :: solai  
  real(r8)                   , intent(out) :: rhoair 
  real(r8)                   , intent(out) :: swdown 



  real(r8)                              :: pair   


       pair   = sfcprs                   
       thair  = sfctmp * (sfcprs/pair)**(rair/cpair) 

       qair   = q2                       

       eair   = qair*sfcprs / (0.622+0.378*qair)
       rhoair = (sfcprs-0.378*eair) / (rair*sfctmp)

       if(cosz <= 0.) then 
          swdown = 0.
       else
          swdown = soldn
       end if 

       solad(1) = swdown*0.7*0.5     
       solad(2) = swdown*0.7*0.5     
       solai(1) = swdown*0.3*0.5     
       solai(2) = swdown*0.3*0.5     

  end subroutine atm_glacier


  subroutine energy_glacier (nsnow  ,nsoil  ,isnow  ,dt     ,qsnow  ,rhoair , & 
                             eair   ,sfcprs ,qair   ,sfctmp ,lwdn   ,uu     , & 
                             vv     ,solad  ,solai  ,cosz   ,zref   ,         & 
                             tbot   ,zbot   ,zsnso  ,dzsnso ,                 & 
                             tg     ,stc    ,snowh  ,sneqv  ,sneqvo ,sh2o   , & 
                             smc    ,snice  ,snliq  ,albold ,cm     ,ch     , & 
                             tauss  ,qsfc   ,                                 & 
                             imelt  ,snicev ,snliqv ,epore  ,qmelt  ,ponding, & 
                             sag    ,fsa    ,fsr    ,fira   ,fsh    ,fgev   , & 
                             trad   ,t2m    ,ssoil  ,lathea ,q2e    ,emissi, ch2b, &
                             taux, tauy,fire , albd_out, albi_out)  !--cheyz   

  implicit none


  integer                           , intent(in)    :: nsnow  
  integer                           , intent(in)    :: nsoil  
  integer                           , intent(in)    :: isnow  
  real(r8)                      , intent(in)    :: dt     
  real(r8)                      , intent(in)    :: qsnow  
  real(r8)                      , intent(in)    :: rhoair 
  real(r8)                      , intent(in)    :: eair   
  real(r8)                      , intent(in)    :: sfcprs 
  real(r8)                      , intent(in)    :: qair   
  real(r8)                      , intent(in)    :: sfctmp 
  real(r8)                      , intent(in)    :: lwdn   
  real(r8)                      , intent(in)    :: uu     
  real(r8)                      , intent(in)    :: vv     
  real(r8)   , dimension(       1:    2), intent(in)    :: solad  
  real(r8)   , dimension(       1:    2), intent(in)    :: solai  
  real(r8)                      , intent(in)    :: cosz   
  real(r8)                      , intent(in)    :: zref   
  real(r8)                      , intent(in)    :: tbot   
  real(r8)                      , intent(in)    :: zbot   
  real(r8)   , dimension(-nsnow+1:nsoil), intent(in)    :: zsnso  
  real(r8)   , dimension(-nsnow+1:nsoil), intent(in)    :: dzsnso 


  real(r8)                      , intent(inout) :: tg     
  real(r8)   , dimension(-nsnow+1:nsoil), intent(inout) :: stc    
  real(r8)                      , intent(inout) :: snowh  
  real(r8)                      , intent(inout) :: sneqv  
  real(r8)                      , intent(inout) :: sneqvo 
  real(r8)   , dimension(       1:nsoil), intent(inout) :: sh2o   
  real(r8)   , dimension(       1:nsoil), intent(inout) :: smc    
  real(r8)   , dimension(-nsnow+1:    0), intent(inout) :: snice  
  real(r8)                      , dimension(-nsnow+1:    0), intent(inout) :: snliq  
  real(r8)                      , intent(inout) :: albold 
  real(r8)                      , intent(inout) :: cm     
  real(r8)                      , intent(inout) :: ch     
  real(r8)                      , intent(inout) :: tauss  
  real(r8)                      , intent(inout) :: qsfc   
  integer, dimension(-nsnow+1:nsoil), intent(out)   :: imelt  
  real(r8)                      , dimension(-nsnow+1:    0), intent(out)   :: snicev 
  real(r8)                      , dimension(-nsnow+1:    0), intent(out)   :: snliqv 
  real(r8)                      , dimension(-nsnow+1:    0), intent(out)   :: epore  
  real(r8)                      , intent(out)   :: qmelt  
  real(r8)                      , intent(out)   :: ponding
  real(r8)                      , intent(out)   :: sag    
  real(r8)                      , intent(out)   :: fsa    
  real(r8)                      , intent(out)   :: fsr    
  real(r8)                      , intent(out)   :: fira   
  real(r8)                      , intent(out)   :: fsh    
  real(r8)                      , intent(out)   :: fgev   
  real(r8)                      , intent(out)   :: trad   
  real(r8)                      , intent(out)   :: t2m    
  real(r8)                      , intent(out)   :: ssoil  
  real(r8)                      , intent(out)   :: lathea 
  real(r8)                      , intent(out)   :: q2e
  real(r8)                      , intent(out)   :: emissi
  real(r8)                      , intent(out)   :: ch2b
  real(r8)                      , intent(out)   :: taux, tauy, fire   !---cheyz 20200427   
  real(r8)                        , dimension(       1:    2), intent(out)   :: albi_out,albd_out 
  real(r8)                      :: ur     
  real(r8)                      :: zlvl   
  real(r8)                      :: rsurf  
  real(r8)                      :: zpd    
  real(r8)                      :: z0mg   
  real(r8)                      :: emg    
  !real(r8)                      :: fire      !
  real(r8),   dimension(-nsnow+1:nsoil)                   :: fact   
  real(r8),   dimension(-nsnow+1:nsoil)                   :: df     
  real(r8),   dimension(-nsnow+1:nsoil)                   :: hcpct  
  real(r8)                      :: gamma  
  real(r8)                      :: rhsur  



    ur = max( sqrt(uu**2.+vv**2.), 1. )

     z0mg = z0sno
     zpd  = snowh

     zlvl = zpd + zref



  call thermoprop_glacier (nsoil   ,nsnow   ,isnow   ,dzsnso  ,          & 
                           dt      ,snowh   ,snice   ,snliq   ,          & 
                           df      ,hcpct   ,snicev  ,snliqv  ,epore   , & 
                           fact    )                                       


  call  radiation_glacier (dt      ,tg      ,sneqvo  ,sneqv   ,cosz    , & 
                           qsnow   ,solad   ,solai   ,                   & 
                           albold  ,tauss   ,                            & 
                           sag     ,fsr     ,fsa,  albd_out, albi_out   )  !--cheyz                              



     emg = 0.98



     rhsur = 1.0
     rsurf = 1.0



     lathea = hsub
     gamma = cpair*sfcprs/(0.622*lathea)



    call glacier_flux (nsoil   ,nsnow   ,emg     ,isnow   ,df      ,dzsnso  ,z0mg    , & 
                       zlvl    ,zpd     ,qair    ,sfctmp  ,rhoair  ,sfcprs  , & 
		       ur      ,gamma   ,rsurf   ,lwdn    ,rhsur   ,smc     , & 
		       eair    ,stc     ,sag     ,snowh   ,lathea  ,sh2o    , & 
		       cm      ,ch      ,tg      ,qsfc    ,          & 
		       fira    ,fsh     ,fgev    ,ssoil   ,          & 
		       t2m     ,q2e     ,ch2b)                         



    fire = lwdn + fira

    if(fire <=0.) print*, (&
"stop in noah-mp: emitted longwave <0")

    
    emissi = emg

    
    trad = ( ( fire - (1-emissi)*lwdn ) / (emissi*sb) ) ** 0.25



    call tsnosoi_glacier (nsoil   ,nsnow   ,isnow   ,dt      ,tbot    , & 
                          ssoil   ,snowh   ,zbot    ,zsnso   ,df      , & 
		          hcpct   ,                                     & 
                          stc     )                                       


     if(opt_stc == 2) then
      if (snowh > 0.05 .and. tg > tfrz) tg = tfrz
     end if



 call phasechange_glacier (nsnow   ,nsoil   ,isnow   ,dt      ,fact    , & 
                           dzsnso  ,                                     & 
                           stc     ,snice   ,snliq   ,sneqv   ,snowh   , & 
                           smc     ,sh2o    ,                            & 
                           qmelt   ,imelt   ,ponding )                     

   taux = -rhoair*cm*ur*uu  !--cheyz 20200427
   tauy = -rhoair*cm*ur*vv  !--cheyz 20200427
  end subroutine energy_glacier

  subroutine thermoprop_glacier (nsoil   ,nsnow   ,isnow   ,dzsnso  , & 
                                 dt      ,snowh   ,snice   ,snliq   , & 
                                 df      ,hcpct   ,snicev  ,snliqv  ,epore   , & 
                                 fact    )                                       


  implicit none


  integer                        , intent(in)  :: nsoil   
  integer                        , intent(in)  :: nsnow   
  integer                        , intent(in)  :: isnow   
  real(r8)                   , intent(in)  :: dt      
  real(r8),   dimension(-nsnow+1:    0), intent(in)  :: snice   
  real(r8),   dimension(-nsnow+1:    0), intent(in)  :: snliq   
  real(r8),   dimension(-nsnow+1:nsoil), intent(in)  :: dzsnso  
  real(r8)                   , intent(in)  :: snowh   


  real(r8),   dimension(-nsnow+1:nsoil), intent(out) :: df      
  real(r8),   dimension(-nsnow+1:nsoil), intent(out) :: hcpct   
  real(r8),   dimension(-nsnow+1:    0), intent(out) :: snicev  
  real(r8),   dimension(-nsnow+1:    0), intent(out) :: snliqv  
  real(r8),   dimension(-nsnow+1:    0), intent(out) :: epore   
  real(r8),   dimension(-nsnow+1:nsoil), intent(out) :: fact    



  integer :: iz, iz2
  real(r8),   dimension(-nsnow+1:    0)              :: cvsno   
  real(r8),   dimension(-nsnow+1:    0)              :: tksno   
  real(r8)                      :: zmid    




    call csnow_glacier (isnow   ,nsnow   ,nsoil   ,snice   ,snliq   ,dzsnso  , & 
                        tksno   ,cvsno   ,snicev  ,snliqv  ,epore   )   

    do iz = isnow+1, 0
      df   (iz) = tksno(iz)
      hcpct(iz) = cvsno(iz)
    end do



    do  iz = 1, nsoil
       zmid      = 0.5 * (dzsnso(iz))
       do iz2 = 1, iz-1
         zmid = zmid + dzsnso(iz2)
       end do
       hcpct(iz) = 1.e6 * ( 0.8194 + 0.1309*zmid )
       df(iz)    = 0.32333 + ( 0.10073 * zmid )
    end do
       


    do iz = isnow+1,nsoil
     fact(iz) = dt/(hcpct(iz)*dzsnso(iz))
    end do



    if(isnow == 0) then
       df(1) = (df(1)*dzsnso(1)+0.35*snowh)      / (snowh    +dzsnso(1)) 
    else
       df(1) = (df(1)*dzsnso(1)+df(0)*dzsnso(0)) / (dzsnso(0)+dzsnso(1))
    end if


  end subroutine thermoprop_glacier


  subroutine csnow_glacier (isnow   ,nsnow   ,nsoil   ,snice   ,snliq   ,dzsnso  , & 
                            tksno   ,cvsno   ,snicev  ,snliqv  ,epore   )   



  implicit none



  integer,                          intent(in) :: isnow  
  integer                        ,  intent(in) :: nsnow  
  integer                        ,  intent(in) :: nsoil  
  real(r8),   dimension(-nsnow+1:    0),  intent(in) :: snice  
  real(r8),   dimension(-nsnow+1:    0),  intent(in) :: snliq  
  real(r8),   dimension(-nsnow+1:nsoil),  intent(in) :: dzsnso 



  real(r8),   dimension(-nsnow+1:    0), intent(out) :: cvsno  
  real(r8),   dimension(-nsnow+1:    0), intent(out) :: tksno  
  real(r8),   dimension(-nsnow+1:    0), intent(out) :: snicev 
  real(r8),   dimension(-nsnow+1:    0), intent(out) :: snliqv 
  real(r8),   dimension(-nsnow+1:    0), intent(out) :: epore  



  integer :: iz
  real(r8),   dimension(-nsnow+1:    0) :: bdsnoi  




  do iz = isnow+1, 0
      snicev(iz)   = min(1., snice(iz)/(dzsnso(iz)*denice) )
      epore(iz)    = 1. - snicev(iz)
      snliqv(iz)   = min(epore(iz),snliq(iz)/(dzsnso(iz)*denh2o))
  enddo

  do iz = isnow+1, 0
      bdsnoi(iz) = (snice(iz)+snliq(iz))/dzsnso(iz)
      cvsno(iz) = cice*snicev(iz)+cwat*snliqv(iz)

  enddo



  do iz = isnow+1, 0
     tksno(iz) = 3.2217e-6*bdsnoi(iz)**2.           




  enddo

  end subroutine csnow_glacier

  subroutine radiation_glacier (dt      ,tg      ,sneqvo  ,sneqv   ,cosz    , & 
                                qsnow   ,solad   ,solai   ,                   & 
                                albold  ,tauss   ,                            & 
                                sag     ,fsr     ,fsa,  albd_out, albi_out   )  !--cheyz                              

  implicit none


  real(r8),   intent(in)                     :: dt     
  real(r8),   intent(in)                     :: tg     
  real(r8),   intent(in)                     :: sneqvo 
  real(r8),   intent(in)                     :: sneqv  
  real(r8),   intent(in)                     :: cosz   
  real(r8),   intent(in)                     :: qsnow  
  real(r8),   dimension(1:2)    , intent(in) :: solad  
  real(r8),   dimension(1:2)    , intent(in) :: solai  


  real(r8),                    intent(inout) :: albold 
  real(r8),                    intent(inout) :: tauss  

  real(r8)  , dimension(1:2), intent(out)                 :: albd_out    !--cheyz
  real(r8)  , dimension(1:2), intent(out)                 :: albi_out    !--cheyz

  real(r8),   intent(out)                    :: sag    
  real(r8),   intent(out)                    :: fsr    
  real(r8),   intent(out)                    :: fsa    


  integer                              :: ib     
  integer                              :: nband  
  real(r8)                      :: fage   
  real(r8),   dimension(1:2)                 :: albsnd 
  real(r8),   dimension(1:2)                 :: albsni 
  real(r8)                      :: alb    
  real(r8)                      :: abs    
  real(r8)                      :: ref    
  real(r8)                      :: fsno   
  real(r8),   dimension(1:2)                 :: albice 

  real(r8),  parameter :: mpe = 1.e-6



  nband = 2
  albsnd = 0.0
  albsni = 0.0
  albice(1) = 0.80    
  albice(2) = 0.55


  call snow_age_glacier (dt,tg,sneqvo,sneqv,tauss,fage)



  if(opt_alb == 1) &
     call snowalb_bats_glacier (nband,cosz,fage,albsnd,albsni)
  if(opt_alb == 2) then
     call snowalb_class_glacier(nband,qsnow,dt,alb,albold,albsnd,albsni)
     albold = alb
  end if

#ifdef AMIPC_PHYSICS
  !--cheyz  2021/6/25-- increase snow alb based on ERA5 reanalysis!
  albsnd=albsnd+0.19
  albsni=albsni+0.19
  albsnd=min(albsnd, 1.0)
  albsni=min(albsni, 1.0)
#endif


   sag = 0.
   fsa = 0.
   fsr = 0.
   
   fsno = 0.0
   if(sneqv > 0.0) fsno = 1.0

  do ib = 1, nband

    albsnd(ib) = albice(ib)*(1.-fsno) + albsnd(ib)*fsno
    albsni(ib) = albice(ib)*(1.-fsno) + albsni(ib)*fsno



    abs = solad(ib)*(1.-albsnd(ib)) + solai(ib)*(1.-albsni(ib))
    sag = sag + abs
    fsa = fsa + abs
    
    ref = solad(ib)*albsnd(ib) + solai(ib)*albsni(ib)
    fsr = fsr + ref
    
  end do
   
  !--cheyz  2020/04/27 should be tested!
#ifndef REGLSM
  do ib=1, nband
   !albd_out(ib)=solad(ib)*albsnd(ib)
   !albi_out(ib)=solai(ib)*albsni(ib)
   albd_out(ib)=albsnd(ib) !yizhang
   albi_out(ib)=albsni(ib) !yizhang 2021/01/26
  end do
#else
  do ib=1, nband
   albd_out(ib)=solad(ib)*albsnd(ib)
   albi_out(ib)=solai(ib)*albsni(ib)
  end do
#endif

  end subroutine radiation_glacier

  subroutine snow_age_glacier (dt,tg,sneqvo,sneqv,tauss,fage)

  implicit none




   real(r8),   intent(in) :: dt        
   real(r8),   intent(in) :: tg        
   real(r8),   intent(in) :: sneqvo    
   real(r8),   intent(in) :: sneqv     


  real(r8),    intent(inout) :: tauss  


   real(r8),   intent(out) :: fage     


   real(r8)                      :: tage       
   real(r8)                      :: age1       
   real(r8)                      :: age2       
   real(r8)                      :: age3       
   real(r8)                      :: dela       
   real(r8)                      :: sge        
   real(r8)                      :: dels       
   real(r8)                      :: dela0      
   real(r8)                      :: arg        



   if(sneqv.le.0.0) then
          tauss = 0.
   else if (sneqv.gt.800.) then
          tauss = 0.
   else

          dela0 = 1.e-6*dt
          arg   = 5.e3*(1./tfrz-1./tg)
          age1  = exp(arg)
          age2  = exp(amin1(0.,10.*arg))
          age3  = 0.3
          tage  = age1+age2+age3
          dela  = dela0*tage
          dels  = amax1(0.0,sneqv-sneqvo) / swemx
          sge   = (tauss+dela)*(1.0-dels)
          tauss = amax1(0.,sge)
   endif

   fage= tauss/(tauss+1.)

  end subroutine snow_age_glacier


  subroutine snowalb_bats_glacier (nband,cosz,fage,albsnd,albsni)

  implicit none



  integer,intent(in) :: nband  

  real(r8),  intent(in) :: cosz    
  real(r8),  intent(in) :: fage    



  real(r8),   dimension(1:2),intent(out) :: albsnd 
  real(r8),   dimension(1:2),intent(out) :: albsni 


  real(r8)                      :: fzen                 
  real(r8)                      :: cf1                  
  real(r8)                      :: sl2                  
  real(r8)                      :: sl1                  
  real(r8)                      :: sl                   
  real(r8),   parameter :: c1 = 0.2  
  real(r8),   parameter :: c2 = 0.5  





        albsnd(1: nband) = 0.
        albsni(1: nband) = 0.



        sl=2.0
        sl1=1./sl
        sl2=2.*sl
        cf1=((1.+sl1)/(1.+sl2*cosz)-sl1)
        fzen=amax1(cf1,0.)

        albsni(1)=0.95*(1.-c1*fage)         
        albsni(2)=0.65*(1.-c2*fage)        
        albsnd(1)=albsni(1)+0.4*fzen*(1.-albsni(1))    
        albsnd(2)=albsni(2)+0.4*fzen*(1.-albsni(2))    

  end subroutine snowalb_bats_glacier


  subroutine snowalb_class_glacier (nband,qsnow,dt,alb,albold,albsnd,albsni)

  implicit none



  integer,intent(in) :: nband  

  real(r8),  intent(in) :: qsnow     
  real(r8),  intent(in) :: dt        
  real(r8),  intent(in) :: albold    



  real(r8),                  intent(inout) :: alb        


  real(r8),   dimension(1:2),intent(out) :: albsnd 
  real(r8),   dimension(1:2),intent(out) :: albsni 





        albsnd(1: nband) = 0.
        albsni(1: nband) = 0.



         alb = 0.55 + (albold-0.55) * exp(-0.01*dt/3600.)



         if (qsnow > 0.) then
           alb = alb + min(qsnow*dt,swemx) * (0.84-alb)/(swemx)
         endif

         albsni(1)= alb
         albsni(2)= alb
         albsnd(1)= alb
         albsnd(2)= alb

  end subroutine snowalb_class_glacier

  subroutine glacier_flux (nsoil   ,nsnow   ,emg     ,isnow   ,df      ,dzsnso  ,z0m , &
                           zlvl    ,zpd     ,qair    ,sfctmp  ,rhoair  ,sfcprs  , &
           ur      ,gamma   ,rsurf   ,lwdn    ,rhsur   ,smc     , &
           eair    ,stc     ,sag     ,snowh   ,lathea  ,sh2o    , &
                           cm      ,ch      ,tgb     ,qsfc    ,       &
                           irb     ,shb     ,evb     ,ghb     ,       &
                           t2mb    ,q2b     ,ehb2)                         










  implicit none


  integer, intent(in)                         :: nsnow  
  integer, intent(in)                         :: nsoil  
  real(r8),                              intent(in) :: emg    
  integer,                         intent(in) :: isnow  
  real(r8),   dimension(-nsnow+1:nsoil), intent(in) :: df     
  real(r8),   dimension(-nsnow+1:nsoil), intent(in) :: dzsnso 
  real(r8),                              intent(in) :: z0m    
  real(r8),                              intent(in) :: zlvl   
  real(r8),                              intent(in) :: zpd    
  real(r8),                              intent(in) :: qair   
  real(r8),                              intent(in) :: sfctmp 
  real(r8),                              intent(in) :: rhoair 
  real(r8),                              intent(in) :: sfcprs 
  real(r8),                              intent(in) :: ur     
  real(r8),                              intent(in) :: gamma  
  real(r8),                              intent(in) :: rsurf  
  real(r8),                              intent(in) :: lwdn   
  real(r8),                              intent(in) :: rhsur  
  real(r8),                              intent(in) :: eair   
  real(r8),   dimension(-nsnow+1:nsoil), intent(in) :: stc    
  real(r8),   dimension(       1:nsoil), intent(in) :: smc    
  real(r8),   dimension(       1:nsoil), intent(in) :: sh2o   
  real(r8),                              intent(in) :: sag    
  real(r8),                              intent(in) :: snowh  
  real(r8),                              intent(in) :: lathea 


  real(r8),                           intent(inout) :: cm     
  real(r8),                           intent(inout) :: ch     
  real(r8),                           intent(inout) :: tgb    
  real(r8),                           intent(inout) :: qsfc   



  real(r8),                             intent(out) :: irb    
  real(r8),                             intent(out) :: shb    
  real(r8),                             intent(out) :: evb    
  real(r8),                             intent(out) :: ghb    
  real(r8),                             intent(out) :: t2mb   
  real(r8),                             intent(out) :: q2b    
  real(r8),                             intent(out) :: ehb2   



  integer :: niterb  
  real(r8)                      :: mpe     
  real(r8)                      :: dtg        
  integer :: mozsgn  
  real(r8)                      :: mozold     
  real(r8)                      :: fm2          
  real(r8)                      :: fh2          
  real(r8)                      :: ch2          
  real(r8)                      :: h          
  real(r8)                      :: fv         
  real(r8)                      :: cir        
  real(r8)                      :: cgh        
  real(r8)                      :: csh        
  real(r8)                      :: cev        
  real(r8)                      :: cq2b       
  integer :: iter    
  real(r8)                      :: z0h        
  real(r8)                      :: moz        
  real(r8)                      :: fm         
  real(r8)                      :: fh         
  real(r8)                      :: ramb       
  real(r8)                      :: rahb       
  real(r8)                      :: rawb       
  real(r8)                      :: estg       
  real(r8)                      :: destg      
  real(r8)                      :: esatw      
  real(r8)                      :: esati      
  real(r8)                      :: dsatw      
  real(r8)                      :: dsati      
  real(r8)                      :: a          
  real(r8)                      :: b          
  real(r8)                      :: t, tdc     
  real(r8),   dimension(       1:nsoil) :: sice   

  tdc(t)   = min( 50., max(-50.,(t-tfrz)) )




        niterb = 5
        mpe    = 1e-6
        dtg    = 0.
        mozsgn = 0
        mozold = 0.
        h      = 0.
        fv     = 0.1

        cir = emg*sb
        cgh = 2.*df(isnow+1)/dzsnso(isnow+1)


      loop3: do iter = 1, niterb  

        z0h = z0m 



        call sfcdif1_glacier(iter   ,zlvl   ,zpd    ,z0h    ,z0m    , & 
                     qair   ,sfctmp ,h      ,rhoair ,mpe    ,ur     , & 
       &             moz    ,mozsgn ,fm     ,fh     ,fm2    ,fh2    , & 
       &             fv     ,cm     ,ch     ,ch2)                       

        ramb = max(1.,1./(cm*ur))
        rahb = max(1.,1./(ch*ur))
        rawb = rahb



        t = tdc(tgb)
        call esat(t, esatw, esati, dsatw, dsati)
        if (t .gt. 0.) then
            estg  = esatw
            destg = dsatw
        else
            estg  = esati
            destg = dsati
        end if

        csh = rhoair*cpair/rahb
	if(snowh > 0.0 .or. opt_gla == 1) then
          cev = rhoair*cpair/gamma/(rsurf+rawb)
	else
	  cev = 0.0   
	end if



        irb   = cir * tgb**4 - emg*lwdn
        shb   = csh * (tgb        - sfctmp      )
        evb   = cev * (estg*rhsur - eair        )
        ghb   = cgh * (tgb        - stc(isnow+1))

        b     = sag-irb-shb-evb-ghb
        a     = 4.*cir*tgb**3 + csh + cev*destg + cgh
        dtg   = b/a

        irb = irb + 4.*cir*tgb**3*dtg
        shb = shb + csh*dtg
        evb = evb + cev*destg*dtg
        ghb = ghb + cgh*dtg


        tgb = tgb + dtg


        h = csh * (tgb - sfctmp)

        t = tdc(tgb)
        call esat(t, esatw, esati, dsatw, dsati)
        if (t .gt. 0.) then
            estg  = esatw
        else
            estg  = esati
        end if
#ifndef AMIPW_PHYSICS
        qsfc = 0.622*(estg*rhsur)/(sfcprs-0.378*(estg*rhsur))
#endif

     end do loop3 




     sice = smc - sh2o
     if(opt_stc == 1 .or. opt_stc ==3) then
     if ((maxval(sice) > 0.0 .or. snowh > 0.0) .and. tgb > tfrz .and. opt_gla == 1) then
          tgb = tfrz
          t = tdc(tgb)                              
          call esat(t, esatw, esati, dsatw, dsati)
          estg  = esati
#ifndef AMIPW_PHYSICS
          qsfc = 0.622*(estg*rhsur)/(sfcprs-0.378*(estg*rhsur))
#endif
          irb = cir * tgb**4 - emg*lwdn
          shb = csh * (tgb        - sfctmp)
          evb = cev * (estg*rhsur - eair )          
          ghb = sag - (irb+shb+evb)
     end if
     end if


     ehb2  = fv*vkc/(log((2.+z0h)/z0h)-fh2)
     cq2b  = ehb2
     if (ehb2.lt.1.e-5 ) then
       t2mb  = tgb
       q2b   = qsfc
     else
       t2mb  = tgb - shb/(rhoair*cpair) * 1./ehb2
       q2b   = qsfc - evb/(lathea*rhoair)*(1./cq2b + rsurf)
     endif


     ch = 1./rahb

  end subroutine glacier_flux

  subroutine esat(t, esw, esi, desw, desi)



  implicit none



  real(r8),   intent(in)  :: t              



  real(r8),   intent(out) :: esw            
  real(r8),   intent(out) :: esi            
  real(r8),   intent(out) :: desw           
  real(r8),   intent(out) :: desi           



  real(r8)                      :: a0,a1,a2,a3,a4,a5,a6  
  real(r8)                      :: b0,b1,b2,b3,b4,b5,b6  
  real(r8)                      :: c0,c1,c2,c3,c4,c5,c6  
  real(r8)                      :: d0,d1,d2,d3,d4,d5,d6  

  parameter (a0=6.107799961    , a1=4.436518521e-01,  &
             a2=1.428945805e-02, a3=2.650648471e-04,  &
             a4=3.031240396e-06, a5=2.034080948e-08,  &
             a6=6.136820929e-11)

  parameter (b0=6.109177956    , b1=5.034698970e-01,  &
             b2=1.886013408e-02, b3=4.176223716e-04,  &
             b4=5.824720280e-06, b5=4.838803174e-08,  &
             b6=1.838826904e-10)

  parameter (c0= 4.438099984e-01, c1=2.857002636e-02,  &
             c2= 7.938054040e-04, c3=1.215215065e-05,  &
             c4= 1.036561403e-07, c5=3.532421810e-10,  &
             c6=-7.090244804e-13)

  parameter (d0=5.030305237e-01, d1=3.773255020e-02,  &
             d2=1.267995369e-03, d3=2.477563108e-05,  &
             d4=3.005693132e-07, d5=2.158542548e-09,  &
             d6=7.131097725e-12)

  esw  = 100.*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
  esi  = 100.*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
  desw = 100.*(c0+t*(c1+t*(c2+t*(c3+t*(c4+t*(c5+t*c6))))))
  desi = 100.*(d0+t*(d1+t*(d2+t*(d3+t*(d4+t*(d5+t*d6))))))

  end subroutine esat


  subroutine sfcdif1_glacier(iter   ,zlvl   ,zpd    ,z0h    ,z0m    , & 
                     qair   ,sfctmp ,h      ,rhoair ,mpe    ,ur     , & 
       &             moz    ,mozsgn ,fm     ,fh     ,fm2    ,fh2    , & 
       &             fv     ,cm     ,ch     ,ch2     )                  



    implicit none


    integer,              intent(in) :: iter   
    real(r8),                   intent(in) :: zlvl   
    real(r8),                   intent(in) :: zpd    
    real(r8),                   intent(in) :: z0h    
    real(r8),                   intent(in) :: z0m    
    real(r8),                   intent(in) :: qair   
    real(r8),                   intent(in) :: sfctmp 
    real(r8),                   intent(in) :: h      
    real(r8),                   intent(in) :: rhoair 
    real(r8),                   intent(in) :: mpe    
    real(r8),                   intent(in) :: ur     


    real(r8),                intent(inout) :: moz    
    integer,           intent(inout) :: mozsgn 
    real(r8),                intent(inout) :: fm     
    real(r8),                intent(inout) :: fh     
    real(r8),                intent(inout) :: fm2    
    real(r8),                intent(inout) :: fh2    


    real(r8),                  intent(out) :: fv     
    real(r8),                  intent(out) :: cm     
    real(r8),                  intent(out) :: ch     
    real(r8),                  intent(out) :: ch2    


    real(r8)                      :: mozold                   
    real(r8)                      :: tmpcm                    
    real(r8)                      :: tmpch                    
    real(r8)                      :: mol                      
    real(r8)                      :: tvir                     
    real(r8)                      :: tmp1,tmp2,tmp3           
    real(r8)                      :: fmnew                    
    real(r8)                      :: fhnew                    
    real(r8)                      :: moz2                     
    real(r8)                      :: tmpcm2                   
    real(r8)                      :: tmpch2                   
    real(r8)                      :: fm2new                   
    real(r8)                      :: fh2new                   
    real(r8)                      :: tmp12,tmp22,tmp32        

    real(r8)                      :: cmfm, chfh, cm2fm2, ch2fh2





    mozold = moz
   
 
    if(zlvl <= zpd) then
       write(*,*) 'critical glacier problem: zlvl <= zpd; model stops', zlvl, zpd
       print*, (&
"stop in noah-mp glacier")
    endif

    tmpcm = log((zlvl-zpd) / z0m)
    tmpch = log((zlvl-zpd) / z0h)
    tmpcm2 = log((2.0 + z0m) / z0m)
    tmpch2 = log((2.0 + z0h) / z0h)

    if(iter == 1) then
       fv   = 0.0
       moz  = 0.0
       mol  = 0.0
       moz2 = 0.0
    else
       tvir = (1. + 0.61*qair) * sfctmp
       tmp1 = vkc * (grav/tvir) * h/(rhoair*cpair)
       if (abs(tmp1) .le. mpe) tmp1 = mpe
       mol  = -1. * fv**3 / tmp1
       moz  = min( (zlvl-zpd)/mol, 1.)
       moz2  = min( (2.0 + z0h)/mol, 1.)
    endif



    if (mozold*moz .lt. 0.) mozsgn = mozsgn+1
    if (mozsgn .ge. 2) then
       moz = 0.
       fm = 0.
       fh = 0.
       moz2 = 0.
       fm2 = 0.
       fh2 = 0.
    endif


    if (moz .lt. 0.) then
       tmp1 = (1. - 16.*moz)**0.25
       tmp2 = log((1.+tmp1*tmp1)/2.)
       tmp3 = log((1.+tmp1)/2.)
       fmnew = 2.*tmp3 + tmp2 - 2.*atan(tmp1) + 1.5707963
       fhnew = 2*tmp2


       tmp12 = (1. - 16.*moz2)**0.25
       tmp22 = log((1.+tmp12*tmp12)/2.)
       tmp32 = log((1.+tmp12)/2.)
       fm2new = 2.*tmp32 + tmp22 - 2.*atan(tmp12) + 1.5707963
       fh2new = 2*tmp22
    else
       fmnew = -5.*moz
       fhnew = fmnew
       fm2new = -5.*moz2
       fh2new = fm2new
    endif




    if (iter == 1) then
       fm = fmnew
       fh = fhnew
       fm2 = fm2new
       fh2 = fh2new
    else
       fm = 0.5 * (fm+fmnew)
       fh = 0.5 * (fh+fhnew)
       fm2 = 0.5 * (fm2+fm2new)
       fh2 = 0.5 * (fh2+fh2new)
    endif



    fh = min(fh,0.9*tmpch)
    fm = min(fm,0.9*tmpcm)
    fh2 = min(fh2,0.9*tmpch2)
    fm2 = min(fm2,0.9*tmpcm2)

    cmfm = tmpcm-fm
    chfh = tmpch-fh
    cm2fm2 = tmpcm2-fm2
    ch2fh2 = tmpch2-fh2
    if(abs(cmfm) <= mpe) cmfm = mpe
    if(abs(chfh) <= mpe) chfh = mpe
    if(abs(cm2fm2) <= mpe) cm2fm2 = mpe
    if(abs(ch2fh2) <= mpe) ch2fh2 = mpe
    cm  = vkc*vkc/(cmfm*cmfm)
    ch  = vkc*vkc/(cmfm*chfh)
    ch2  = vkc*vkc/(cm2fm2*ch2fh2)
        


    fv = ur * sqrt(cm)
    ch2  = vkc*fv/ch2fh2

  end subroutine sfcdif1_glacier

  subroutine tsnosoi_glacier (nsoil   ,nsnow   ,isnow   ,dt      ,tbot    , & 
                              ssoil   ,snowh   ,zbot    ,zsnso   ,df      , & 
			      hcpct   ,                                     & 
                              stc     )                                       





  implicit none



    integer,                         intent(in)  :: nsoil  
    integer,                         intent(in)  :: nsnow  
    integer,                         intent(in)  :: isnow  

    real(r8),                              intent(in)  :: dt     
    real(r8),                              intent(in)  :: tbot   
    real(r8),                              intent(in)  :: ssoil  
    real(r8),                              intent(in)  :: snowh  
    real(r8),                              intent(in)  :: zbot   
    real(r8),   dimension(-nsnow+1:nsoil), intent(in)  :: zsnso  
    real(r8),   dimension(-nsnow+1:nsoil), intent(in)  :: df     
    real(r8),   dimension(-nsnow+1:nsoil), intent(in)  :: hcpct  



    real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: stc



    integer                                      :: iz
    real(r8)                      :: zbotsno   
    real(r8),   dimension(-nsnow+1:nsoil)              :: ai, bi, ci, rhsts
    real(r8)                      :: eflxb 
    real(r8),   dimension(-nsnow+1:nsoil)              :: phi   





    phi(isnow+1:nsoil) = 0.



    zbotsno = zbot - snowh    



      call hrt_glacier   (nsnow     ,nsoil     ,isnow     ,zsnso     , &
                          stc       ,tbot      ,zbotsno   ,df        , &
                          hcpct     ,ssoil     ,phi       ,            &
                          ai        ,bi        ,ci        ,rhsts     , &
                          eflxb     )

      call hstep_glacier (nsnow     ,nsoil     ,isnow     ,dt        , &
                          ai        ,bi        ,ci        ,rhsts     , &
                          stc       ) 

  end subroutine tsnosoi_glacier


  subroutine hrt_glacier (nsnow     ,nsoil     ,isnow     ,zsnso     , & 
                          stc       ,tbot      ,zbot      ,df        , & 
                          hcpct     ,ssoil     ,phi       ,            & 
                          ai        ,bi        ,ci        ,rhsts     , & 
                          botflx    )                                    






    implicit none



    integer,                         intent(in)  :: nsoil  
    integer,                         intent(in)  :: nsnow  
    integer,                         intent(in)  :: isnow  
    real(r8),                              intent(in)  :: tbot   
    real(r8),                              intent(in)  :: zbot   
                                                           
    real(r8),                              intent(in)  :: ssoil  
    real(r8),   dimension(-nsnow+1:nsoil), intent(in)  :: zsnso  
    real(r8),   dimension(-nsnow+1:nsoil), intent(in)  :: stc    
    real(r8),   dimension(-nsnow+1:nsoil), intent(in)  :: df     
    real(r8),   dimension(-nsnow+1:nsoil), intent(in)  :: hcpct  
    real(r8),   dimension(-nsnow+1:nsoil), intent(in)  :: phi    



    real(r8),   dimension(-nsnow+1:nsoil), intent(out) :: rhsts  
    real(r8),   dimension(-nsnow+1:nsoil), intent(out) :: ai     
    real(r8),   dimension(-nsnow+1:nsoil), intent(out) :: bi     
    real(r8),   dimension(-nsnow+1:nsoil), intent(out) :: ci     
    real(r8),                              intent(out) :: botflx 



    integer                                      :: k
    real(r8),   dimension(-nsnow+1:nsoil)              :: ddz
    real(r8),   dimension(-nsnow+1:nsoil)              :: denom
    real(r8),   dimension(-nsnow+1:nsoil)              :: dtsdz
    real(r8),   dimension(-nsnow+1:nsoil)              :: eflux
    real(r8)                      :: temp1


    do k = isnow+1, nsoil
        if (k == isnow+1) then
           denom(k)  = - zsnso(k) * hcpct(k)
           temp1     = - zsnso(k+1)
           ddz(k)    = 2.0 / temp1
           dtsdz(k)  = 2.0 * (stc(k) - stc(k+1)) / temp1
           eflux(k)  = df(k) * dtsdz(k) - ssoil - phi(k)
        else if (k < nsoil) then
           denom(k)  = (zsnso(k-1) - zsnso(k)) * hcpct(k)
           temp1     = zsnso(k-1) - zsnso(k+1)
           ddz(k)    = 2.0 / temp1
           dtsdz(k)  = 2.0 * (stc(k) - stc(k+1)) / temp1
           eflux(k)  = (df(k)*dtsdz(k) - df(k-1)*dtsdz(k-1)) - phi(k)
        else if (k == nsoil) then
           denom(k)  = (zsnso(k-1) - zsnso(k)) * hcpct(k)
           temp1     =  zsnso(k-1) - zsnso(k)
           if(opt_tbot == 1) then
               botflx     = 0. 
           end if
           if(opt_tbot == 2) then
               dtsdz(k)  = (stc(k) - tbot) / ( 0.5*(zsnso(k-1)+zsnso(k)) - zbot)
               botflx    = -df(k) * dtsdz(k)
           end if
           eflux(k)  = (-botflx - df(k-1)*dtsdz(k-1) ) - phi(k)
        end if
    end do

    do k = isnow+1, nsoil
        if (k == isnow+1) then
           ai(k)    =   0.0
           ci(k)    = - df(k)   * ddz(k) / denom(k)
           if (opt_stc == 1 .or. opt_stc == 3) then
              bi(k) = - ci(k)
           end if                                        
           if (opt_stc == 2) then
              bi(k) = - ci(k) + df(k)/(0.5*zsnso(k)*zsnso(k)*hcpct(k))
           end if
        else if (k < nsoil) then
           ai(k)    = - df(k-1) * ddz(k-1) / denom(k) 
           ci(k)    = - df(k  ) * ddz(k  ) / denom(k) 
           bi(k)    = - (ai(k) + ci (k))
        else if (k == nsoil) then
           ai(k)    = - df(k-1) * ddz(k-1) / denom(k) 
           ci(k)    = 0.0
           bi(k)    = - (ai(k) + ci(k))
        end if
           rhsts(k)  = eflux(k)/ (-denom(k))
    end do

  end subroutine hrt_glacier


  subroutine hstep_glacier (nsnow     ,nsoil     ,isnow     ,dt        ,  & 
                            ai        ,bi        ,ci        ,rhsts     ,  & 
                            stc       )                                     



    implicit none



    integer,                         intent(in)    :: nsoil
    integer,                         intent(in)    :: nsnow
    integer,                         intent(in)    :: isnow
    real(r8),                              intent(in)    :: dt


    real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: ai
    real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: bi
    real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: ci
    real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: stc
    real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: rhsts


    integer                                        :: k
    real(r8),   dimension(-nsnow+1:nsoil)                :: rhstsin
    real(r8),   dimension(-nsnow+1:nsoil)                :: ciin


    do k = isnow+1,nsoil
       rhsts(k) =   rhsts(k) * dt
       ai(k)    =      ai(k) * dt
       bi(k)    = 1. + bi(k) * dt
       ci(k)    =      ci(k) * dt
    end do



    do k = isnow+1,nsoil
       rhstsin(k) = rhsts(k)
       ciin(k)    = ci(k)
    end do



    call rosr12_glacier (ci,ai,bi,ciin,rhstsin,rhsts,isnow+1,nsoil,nsnow)



    do k = isnow+1,nsoil
       stc (k) = stc (k) + ci (k)
    end do

  end subroutine hstep_glacier

  subroutine rosr12_glacier (p,a,b,c,d,delta,ntop,nsoil,nsnow)


















    implicit none

    integer, intent(in)   :: ntop           
    integer, intent(in)   :: nsoil,nsnow
    integer               :: k, kk

    real(r8),   dimension(-nsnow+1:nsoil),intent(in):: a, b, d
    real(r8),   dimension(-nsnow+1:nsoil),intent(inout):: c,p,delta




    c (nsoil) = 0.0
    p (ntop) = - c (ntop) / b (ntop)



    delta (ntop) = d (ntop) / b (ntop)



    do k = ntop+1,nsoil
       p (k) = - c (k) * ( 1.0 / (b (k) + a (k) * p (k -1)) )
       delta (k) = (d (k) - a (k)* delta (k -1))* (1.0/ (b (k) + a (k)&
            * p (k -1)))
    end do



    p (nsoil) = delta (nsoil)



    do k = ntop+1,nsoil
       kk = nsoil - k + (ntop-1) + 1
       p (kk) = p (kk) * p (kk +1) + delta (kk)
    end do

  end subroutine rosr12_glacier


  subroutine phasechange_glacier (nsnow   ,nsoil   ,isnow   ,dt      ,fact    , & 
                                  dzsnso  ,                                     & 
                                  stc     ,snice   ,snliq   ,sneqv   ,snowh   , & 
                                  smc     ,sh2o    ,                            & 
                                  qmelt   ,imelt   ,ponding )                     



  implicit none



  integer, intent(in)                             :: nsnow  
  integer, intent(in)                             :: nsoil  
  integer, intent(in)                             :: isnow  
  real(r8),   intent(in)                                :: dt     
  real(r8),   dimension(-nsnow+1:nsoil), intent(in)     :: fact   
  real(r8),   dimension(-nsnow+1:nsoil), intent(in)     :: dzsnso 



  real(r8),   dimension(-nsnow+1:nsoil), intent(inout)  :: stc    
  real(r8),   dimension(-nsnow+1:0)    , intent(inout)  :: snice  
  real(r8),   dimension(-nsnow+1:0)    , intent(inout)  :: snliq  
  real(r8),   intent(inout)                             :: sneqv
  real(r8),   intent(inout)                             :: snowh
  real(r8),   dimension(       1:nsoil), intent(inout)  :: sh2o   
  real(r8),   dimension(       1:nsoil), intent(inout)  :: smc    


  real(r8),                                 intent(out) :: qmelt  
  integer, dimension(-nsnow+1:nsoil), intent(out) :: imelt  
  real(r8),                                 intent(out) :: ponding



  integer                         :: j,k         
  real(r8),   dimension(-nsnow+1:nsoil) :: hm        
  real(r8),   dimension(-nsnow+1:nsoil) :: xm        
  real(r8),   dimension(-nsnow+1:nsoil) :: wmass0
  real(r8),   dimension(-nsnow+1:nsoil) :: wice0 
  real(r8),   dimension(-nsnow+1:nsoil) :: wliq0 
  real(r8),   dimension(-nsnow+1:nsoil) :: mice      
  real(r8),   dimension(-nsnow+1:nsoil) :: mliq      
  real(r8),   dimension(-nsnow+1:nsoil) :: heatr     
  real                            :: temp1     
  real                            :: propor
  real                            :: xmf       




    qmelt   = 0.
    ponding = 0.
    xmf     = 0.

    do j = isnow+1,0           
         mice(j) = snice(j)
         mliq(j) = snliq(j)
    end do

    do j = isnow+1,0           
         imelt(j)    = 0
         hm(j)       = 0.
         xm(j)       = 0.
         wice0(j)    = mice(j)
         wliq0(j)    = mliq(j)
         wmass0(j)   = mice(j) + mliq(j)
    enddo
    
    do j = isnow+1,0
         if (mice(j) > 0. .and. stc(j) >= tfrz) then  
             imelt(j) = 1
         endif
         if (mliq(j) > 0. .and. stc(j)  < tfrz) then  
             imelt(j) = 2
         endif

    enddo



    do j = isnow+1,0
         if (imelt(j) > 0) then
             hm(j) = (stc(j)-tfrz)/fact(j)
             stc(j) = tfrz
         endif

         if (imelt(j) == 1 .and. hm(j) < 0.) then
            hm(j) = 0.
            imelt(j) = 0
         endif
         if (imelt(j) == 2 .and. hm(j) > 0.) then
            hm(j) = 0.
            imelt(j) = 0
         endif
         xm(j) = hm(j)*dt/hfus                           
    enddo



if (opt_gla == 2) then 

    if (isnow == 0 .and. sneqv > 0. .and. stc(1) >= tfrz) then  
        hm(1)    = (stc(1)-tfrz)/fact(1)             
        stc(1)   = tfrz                              
        xm(1)    = hm(1)*dt/hfus                     

        temp1  = sneqv
        sneqv  = max(0.,temp1-xm(1))                 
        propor = sneqv/temp1                         
        snowh  = max(0.,propor * snowh)              
        heatr(1)  = hm(1) - hfus*(temp1-sneqv)/dt    
        if (heatr(1) > 0.) then
              xm(1)  = heatr(1)*dt/hfus             
              stc(1) = stc(1) + fact(1)*heatr(1)     
        else
              xm(1) = 0.                             
              hm(1) = 0.
        endif
        qmelt   = max(0.,(temp1-sneqv))/dt           
        xmf     = hfus*qmelt                         
        ponding = temp1-sneqv                        
    endif

end if  



    do j = isnow+1,0
      if (imelt(j) > 0 .and. abs(hm(j)) > 0.) then

         heatr(j) = 0.
         if (xm(j) > 0.) then                            
            mice(j) = max(0., wice0(j)-xm(j))
            heatr(j) = hm(j) - hfus*(wice0(j)-mice(j))/dt
         else if (xm(j) < 0.) then                      
            mice(j) = min(wmass0(j), wice0(j)-xm(j))  
            heatr(j) = hm(j) - hfus*(wice0(j)-mice(j))/dt
         endif

         mliq(j) = max(0.,wmass0(j)-mice(j))

         if (abs(heatr(j)) > 0.) then
            stc(j) = stc(j) + fact(j)*heatr(j)
            if (mliq(j)*mice(j)>0.) stc(j) = tfrz
         endif

         qmelt = qmelt + max(0.,(wice0(j)-mice(j)))/dt

      endif
    enddo

if (opt_gla == 1) then     

    do j = 1, nsoil            
         mliq(j) =  sh2o(j)            * dzsnso(j) * 1000.
         mice(j) = (smc(j) - sh2o(j))  * dzsnso(j) * 1000.
    end do

    do j = 1,nsoil       
         imelt(j)    = 0
         hm(j)       = 0.
         xm(j)       = 0.
         wice0(j)    = mice(j)
         wliq0(j)    = mliq(j)
         wmass0(j)   = mice(j) + mliq(j)
    enddo
    
    do j = 1,nsoil
         if (mice(j) > 0. .and. stc(j) >= tfrz) then  
             imelt(j) = 1
         endif
         if (mliq(j) > 0. .and. stc(j)  < tfrz) then  
             imelt(j) = 2
         endif

         
         if (isnow == 0 .and. sneqv > 0. .and. j == 1) then
             if (stc(j) >= tfrz) then
                imelt(j) = 1
             endif
         endif
    enddo



    do j = 1,nsoil
         if (imelt(j) > 0) then
             hm(j) = (stc(j)-tfrz)/fact(j)
             stc(j) = tfrz
         endif

         if (imelt(j) == 1 .and. hm(j) < 0.) then
            hm(j) = 0.
            imelt(j) = 0
         endif
         if (imelt(j) == 2 .and. hm(j) > 0.) then
            hm(j) = 0.
            imelt(j) = 0
         endif
         xm(j) = hm(j)*dt/hfus                           
    enddo



    if (isnow == 0 .and. sneqv > 0. .and. xm(1) > 0.) then  
        temp1  = sneqv
        sneqv  = max(0.,temp1-xm(1))  
        propor = sneqv/temp1
        snowh  = max(0.,propor * snowh)
        heatr(1)  = hm(1) - hfus*(temp1-sneqv)/dt  
        if (heatr(1) > 0.) then
              xm(1) = heatr(1)*dt/hfus             
              hm(1) = heatr(1) 
	      imelt(1) = 1                   
        else
              xm(1) = 0.
              hm(1) = 0.
	      imelt(1) = 0                   
        endif
        qmelt   = max(0.,(temp1-sneqv))/dt
        xmf     = hfus*qmelt
        ponding = temp1-sneqv
    endif



    do j = 1,nsoil
      if (imelt(j) > 0 .and. abs(hm(j)) > 0.) then

         heatr(j) = 0.
         if (xm(j) > 0.) then                            
            mice(j) = max(0., wice0(j)-xm(j))
            heatr(j) = hm(j) - hfus*(wice0(j)-mice(j))/dt
         else if (xm(j) < 0.) then                      
            mice(j) = min(wmass0(j), wice0(j)-xm(j))  
            heatr(j) = hm(j) - hfus*(wice0(j)-mice(j))/dt
         endif

         mliq(j) = max(0.,wmass0(j)-mice(j))

         if (abs(heatr(j)) > 0.) then
            stc(j) = stc(j) + fact(j)*heatr(j)
            if (j <= 0) then                             
               if (mliq(j)*mice(j)>0.) stc(j) = tfrz
            end if
         endif

         if (j > 0) xmf = xmf + hfus * (wice0(j)-mice(j))/dt

         if (j < 1) then
            qmelt = qmelt + max(0.,(wice0(j)-mice(j)))/dt
         endif
      endif
    enddo
    heatr = 0.0
    xm = 0.0





    if (any(stc(1:4) > tfrz) .and. any(stc(1:4) < tfrz)) then
      do j = 1,nsoil
        if ( stc(j) > tfrz ) then                                       
	  heatr(j) = (stc(j)-tfrz)/fact(j)
          do k = 1,nsoil
	    if (j .ne. k .and. stc(k) < tfrz .and. heatr(j) > 0.1) then
	      heatr(k) = (stc(k)-tfrz)/fact(k)
	      if (abs(heatr(k)) > heatr(j)) then  
	        heatr(k) = heatr(k) + heatr(j)
		stc(k) = tfrz + heatr(k)*fact(k)
		heatr(j) = 0.0
              else
	        heatr(j) = heatr(j) + heatr(k)
		heatr(k) = 0.0
		stc(k) = tfrz
              end if
	    end if
	  end do
          stc(j) = tfrz + heatr(j)*fact(j)
        end if
      end do
    end if



    if (any(stc(1:4) > tfrz) .and. any(stc(1:4) < tfrz)) then
      do j = 1,nsoil
        if ( stc(j) < tfrz ) then                                       
	  heatr(j) = (stc(j)-tfrz)/fact(j)
          do k = 1,nsoil
	    if (j .ne. k .and. stc(k) > tfrz .and. heatr(j) < -0.1) then
	      heatr(k) = (stc(k)-tfrz)/fact(k)
	      if (heatr(k) > abs(heatr(j))) then  
	        heatr(k) = heatr(k) + heatr(j)
		stc(k) = tfrz + heatr(k)*fact(k)
		heatr(j) = 0.0
              else
	        heatr(j) = heatr(j) + heatr(k)
		heatr(k) = 0.0
		stc(k) = tfrz
              end if
	    end if
	  end do
          stc(j) = tfrz + heatr(j)*fact(j)
        end if
      end do
    end if



    if (any(stc(1:4) > tfrz) .and. any(mice(1:4) > 0.)) then
      do j = 1,nsoil
        if ( stc(j) > tfrz ) then                                       
	  heatr(j) = (stc(j)-tfrz)/fact(j)
          xm(j) = heatr(j)*dt/hfus                           
          do k = 1,nsoil
	    if (j .ne. k .and. mice(k) > 0. .and. xm(j) > 0.1) then
	      if (mice(k) > xm(j)) then  
	        mice(k) = mice(k) - xm(j)
		xmf = xmf + hfus * xm(j)/dt
		stc(k) = tfrz
		xm(j) = 0.0
              else
	        xm(j) = xm(j) - mice(k)
		xmf = xmf + hfus * mice(k)/dt
		mice(k) = 0.0
		stc(k) = tfrz
              end if
              mliq(k) = max(0.,wmass0(k)-mice(k))
	    end if
	  end do
	  heatr(j) = xm(j)*hfus/dt
          stc(j) = tfrz + heatr(j)*fact(j)
        end if
      end do
    end if



    if (any(stc(1:4) < tfrz) .and. any(mliq(1:4) > 0.)) then
      do j = 1,nsoil
        if ( stc(j) < tfrz ) then                                       
	  heatr(j) = (stc(j)-tfrz)/fact(j)
          xm(j) = heatr(j)*dt/hfus                           
          do k = 1,nsoil
	    if (j .ne. k .and. mliq(k) > 0. .and. xm(j) < -0.1) then
	      if (mliq(k) > abs(xm(j))) then  
	        mice(k) = mice(k) - xm(j)
		xmf = xmf + hfus * xm(j)/dt
		stc(k) = tfrz
		xm(j) = 0.0
              else
	        xm(j) = xm(j) + mliq(k)
		xmf = xmf - hfus * mliq(k)/dt
		mice(k) = wmass0(k)
		stc(k) = tfrz
              end if
              mliq(k) = max(0.,wmass0(k)-mice(k))
	    end if
	  end do
	  heatr(j) = xm(j)*hfus/dt
          stc(j) = tfrz + heatr(j)*fact(j)
        end if
      end do
    end if
    
end if   

    do j = isnow+1,0             
       snliq(j) = mliq(j)
       snice(j) = mice(j)
    end do

    do j = 1, nsoil              
      if(opt_gla == 1) then 
       sh2o(j) =  mliq(j)            / (1000. * dzsnso(j))
       sh2o(j) =  max(0.0,min(1.0,sh2o(j)))

      elseif(opt_gla == 2) then 
       sh2o(j) = 0.0             
      end if
      smc(j)  = 1.0 
    end do
   
  end subroutine phasechange_glacier

  subroutine water_glacier (nsnow  ,nsoil  ,imelt  ,dt     ,prcp   ,sfctmp , & 
                            qvap   ,qdew   ,ficeold,zsoil  ,                 & 
                            isnow  ,snowh  ,sneqv  ,snice  ,snliq  ,stc    , & 
                            dzsnso ,sh2o   ,sice   ,ponding,zsnso  ,fsh    , & 
                            runsrf ,runsub ,qsnow  ,ponding1 ,ponding2,qsnbot,fpice     &   
                            )  




  implicit none


  integer,                         intent(in)    :: nsnow   
  integer,                         intent(in)    :: nsoil   
  integer, dimension(-nsnow+1:0) , intent(in)    :: imelt   
  real(r8),                              intent(in)    :: dt      
  real(r8),                              intent(in)    :: prcp    
  real(r8),                              intent(in)    :: sfctmp  
  real(r8),                              intent(inout)    :: qvap    
  real(r8),                              intent(inout)    :: qdew    
  real(r8),   dimension(-nsnow+1:    0), intent(in)    :: ficeold 
  real(r8),   dimension(       1:nsoil), intent(in)    :: zsoil  


  integer,                         intent(inout) :: isnow   
  real(r8),                              intent(inout) :: snowh   
  real(r8),                              intent(inout) :: sneqv   
  real(r8),   dimension(-nsnow+1:    0), intent(inout) :: snice   
  real(r8),   dimension(-nsnow+1:    0), intent(inout) :: snliq   
  real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: stc     
  real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso  
  real(r8),   dimension(       1:nsoil), intent(inout) :: sh2o    
  real(r8),   dimension(       1:nsoil), intent(inout) :: sice    
  real(r8)                   , intent(inout) :: ponding 
  real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: zsnso   
  real(r8)                   , intent(inout) :: fsh     


  real(r8),                              intent(out)   :: runsrf  
  real(r8),                              intent(out)   :: runsub  
  real(r8),                              intent(out)   :: qsnow   
  real(r8),                              intent(out)   :: ponding1
  real(r8),                              intent(out)   :: ponding2
  real(r8),                              intent(out)   :: qsnbot  
  real(r8),                              intent(out)   :: fpice   


  real(r8)  :: qrain   
  real(r8)  :: qseva   
  real(r8)  :: qsdew   
  real(r8)  :: qsnfro  
  real(r8)  :: qsnsub  
  real(r8)  :: snowhin 
  real(r8)  :: snoflow 
  real(r8)  :: bdfall  
  real(r8)  :: replace 
  real(r8),   dimension(       1:nsoil)                :: sice_save  
  real(r8),   dimension(       1:nsoil)                :: sh2o_save  
  integer :: ilev





   snoflow         = 0.
   runsub          = 0.
   runsrf          = 0.
   sice_save       = sice
   sh2o_save       = sh2o






     if(opt_snf == 1 .or. opt_snf == 4) then
       if(sfctmp > tfrz+2.5)then
           fpice = 0.
       else
         if(sfctmp <= tfrz+0.5)then
           fpice = 1.0
         else if(sfctmp <= tfrz+2.)then
           fpice = 1.-(-54.632 + 0.2*sfctmp)
         else
           fpice = 0.6
         endif
       endif
     endif

     if(opt_snf == 2) then
       if(sfctmp >= tfrz+2.2) then
           fpice = 0.
       else
           fpice = 1.0
       endif
     endif

     if(opt_snf == 3) then
       if(sfctmp >= tfrz) then
           fpice = 0.
       else
           fpice = 1.0
       endif
     endif





     bdfall = min(120.,67.92+51.25*exp((sfctmp-tfrz)/2.59)) 

     qrain   = prcp * (1.-fpice)
     qsnow   = prcp * fpice
     snowhin = qsnow/bdfall




     qsnsub = qvap  
     qsnfro = qdew

     call snowwater_glacier (nsnow  ,nsoil  ,imelt  ,dt     ,sfctmp , & 
                             snowhin,qsnow  ,qsnfro ,qsnsub ,qrain  , & 
                             ficeold,zsoil  ,                         & 
                             isnow  ,snowh  ,sneqv  ,snice  ,snliq  , & 
                             sh2o   ,sice   ,stc    ,dzsnso ,zsnso  , & 
                             fsh    ,                                 & 
                             qsnbot ,snoflow,ponding1       ,ponding2)  

    
    
    runsrf = (ponding+ponding1+ponding2)/dt

    if(isnow == 0) then
      runsrf = runsrf + qsnbot + qrain
    else
      runsrf = runsrf + qsnbot
    endif

    
    if(opt_gla == 1) then
      replace = 0.0
      do ilev = 1,nsoil
       replace = replace + dzsnso(ilev)*(sice(ilev) - sice_save(ilev) + sh2o(ilev) - sh2o_save(ilev))
      end do
      replace = replace * 1000.0 / dt     
    
      sice = min(1.0,sice_save)
    elseif(opt_gla == 2) then
      sice = 1.0
    end if
    sh2o = 1.0 - sice
    
    
    

    if(opt_gla == 1) then
      runsub       = snoflow + replace
    elseif(opt_gla == 2) then
      runsub       = snoflow
      qvap = qsnsub
      qdew = qsnfro
    end if

  end subroutine water_glacier


  subroutine snowwater_glacier (nsnow  ,nsoil  ,imelt  ,dt     ,sfctmp , & 
                                snowhin,qsnow  ,qsnfro ,qsnsub ,qrain  , & 
                                ficeold,zsoil  ,                         & 
                                isnow  ,snowh  ,sneqv  ,snice  ,snliq  , & 
                                sh2o   ,sice   ,stc    ,dzsnso ,zsnso  , & 
				fsh    ,                                 & 
                                qsnbot ,snoflow,ponding1       ,ponding2)  

  implicit none


  integer,                         intent(in)    :: nsnow  
  integer,                         intent(in)    :: nsoil  
  integer, dimension(-nsnow+1:0) , intent(in)    :: imelt  
  real(r8),                              intent(in)    :: dt     
  real(r8),                              intent(in)    :: sfctmp 
  real(r8),                              intent(in)    :: snowhin
  real(r8),                              intent(in)    :: qsnow  
  real(r8),                              intent(inout)    :: qsnfro 
  real(r8),                              intent(inout)    :: qsnsub 
  real(r8),                              intent(in)    :: qrain  
  real(r8),   dimension(-nsnow+1:0)    , intent(in)    :: ficeold
  real(r8),   dimension(       1:nsoil), intent(in)    :: zsoil  


  integer,                         intent(inout) :: isnow  
  real(r8),                              intent(inout) :: snowh  
  real(r8),                              intent(inout) :: sneqv  
  real(r8),   dimension(-nsnow+1:    0), intent(inout) :: snice  
  real(r8),   dimension(-nsnow+1:    0), intent(inout) :: snliq  
  real(r8),   dimension(       1:nsoil), intent(inout) :: sh2o   
  real(r8),   dimension(       1:nsoil), intent(inout) :: sice   
  real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: stc    
  real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso 
  real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: zsnso  
  real(r8)                   , intent(inout) :: fsh     


  real(r8),                                intent(out) :: qsnbot 
  real(r8),                                intent(out) :: snoflow
  real(r8),                                intent(out) :: ponding1
  real(r8),                                intent(out) :: ponding2


  integer :: iz
  real(r8)                      :: bdsnow  

   snoflow = 0.0
   ponding1 = 0.0
   ponding2 = 0.0

   call snowfall_glacier (nsoil  ,nsnow  ,dt     ,qsnow  ,snowhin, & 
                          sfctmp ,                                 & 
                          isnow  ,snowh  ,dzsnso ,stc    ,snice  , & 
                          snliq  ,sneqv  )                           

   if(isnow < 0) then        
     call  compact_glacier (nsnow  ,nsoil  ,dt     ,stc    ,snice  , & 
                            snliq  ,imelt  ,ficeold,                 & 
                            isnow  ,dzsnso )                           

     call  combine_glacier (nsnow  ,nsoil  ,                         & 
                            isnow  ,sh2o   ,stc    ,snice  ,snliq  , & 
                            dzsnso ,sice   ,snowh  ,sneqv  ,         & 
                            ponding1       ,ponding2)                  

     call   divide_glacier (nsnow  ,nsoil  ,                         & 
                            isnow  ,stc    ,snice  ,snliq  ,dzsnso )   
   end if



   do iz = -nsnow+1, isnow
        snice(iz) = 0.
        snliq(iz) = 0.
        stc(iz)   = 0.
        dzsnso(iz)= 0.
        zsnso(iz) = 0.
   enddo

   call  snowh2o_glacier (nsnow  ,nsoil  ,dt     ,qsnfro ,qsnsub , & 
                          qrain  ,                                 & 
                          isnow  ,dzsnso ,snowh  ,sneqv  ,snice  , & 
                          snliq  ,sh2o   ,sice   ,stc    ,         & 
			  ponding1       ,ponding2       ,fsh    , & 
                          qsnbot )                                   


       
   if(sneqv > 2000.) then   
      bdsnow      = snice(0) / dzsnso(0)
      snoflow     = (sneqv - 2000.)
      snice(0)    = snice(0)  - snoflow 
      dzsnso(0)   = dzsnso(0) - snoflow/bdsnow
      snoflow     = snoflow / dt
   end if



   if(isnow /= 0) then
       sneqv = 0.
       do iz = isnow+1,0
             sneqv = sneqv + snice(iz) + snliq(iz)
       enddo
   end if



   do iz = isnow+1, 0
        dzsnso(iz) = -dzsnso(iz)
   end do

   dzsnso(1) = zsoil(1)
   do iz = 2,nsoil
        dzsnso(iz) = (zsoil(iz) - zsoil(iz-1))
   end do

   zsnso(isnow+1) = dzsnso(isnow+1)
   do iz = isnow+2 ,nsoil
       zsnso(iz) = zsnso(iz-1) + dzsnso(iz)
   enddo

   do iz = isnow+1 ,nsoil
       dzsnso(iz) = -dzsnso(iz)
   end do

  end subroutine snowwater_glacier

  subroutine snowfall_glacier (nsoil  ,nsnow  ,dt     ,qsnow  ,snowhin , & 
                               sfctmp ,                                  & 
                               isnow  ,snowh  ,dzsnso ,stc    ,snice   , & 
                               snliq  ,sneqv  )                            




    implicit none



  integer,                            intent(in) :: nsoil  
  integer,                            intent(in) :: nsnow  
  real(r8),                                 intent(in) :: dt     
  real(r8),                                 intent(in) :: qsnow  
  real(r8),                                 intent(in) :: snowhin
  real(r8),                                 intent(in) :: sfctmp 



  integer,                         intent(inout) :: isnow  
  real(r8),                              intent(inout) :: snowh  
  real(r8),                              intent(inout) :: sneqv  
  real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso 
  real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: stc    
  real(r8),   dimension(-nsnow+1:    0), intent(inout) :: snice  
  real(r8),   dimension(-nsnow+1:    0), intent(inout) :: snliq  



  integer :: newnode            

    newnode  = 0



    if(isnow == 0 .and. qsnow > 0.)  then
      snowh = snowh + snowhin * dt
      sneqv = sneqv + qsnow * dt
    end if


 
    if(isnow == 0  .and. qsnow>0. .and. snowh >= 0.05) then
      isnow    = -1
      newnode  =  1
      dzsnso(0)= snowh
      snowh    = 0.
      stc(0)   = min(273.16, sfctmp)   
      snice(0) = sneqv
      snliq(0) = 0.
    end if



    if(isnow <  0 .and. newnode == 0 .and. qsnow > 0.) then
         snice(isnow+1)  = snice(isnow+1)   + qsnow   * dt
         dzsnso(isnow+1) = dzsnso(isnow+1)  + snowhin * dt
    endif


  end subroutine snowfall_glacier


  subroutine compact_glacier (nsnow  ,nsoil  ,dt     ,stc    ,snice , & 
                              snliq  ,imelt  ,ficeold,                & 
                              isnow  ,dzsnso )                          


  implicit none


   integer,                         intent(in)    :: nsoil  
   integer,                         intent(in)    :: nsnow  
   integer, dimension(-nsnow+1:0) , intent(in)    :: imelt  
   real(r8),                              intent(in)    :: dt     
   real(r8),   dimension(-nsnow+1:nsoil), intent(in)    :: stc    
   real(r8),   dimension(-nsnow+1:    0), intent(in)    :: snice  
   real(r8),   dimension(-nsnow+1:    0), intent(in)    :: snliq  
   real(r8),   dimension(-nsnow+1:    0), intent(in)    :: ficeold


   integer,                         intent(inout) :: isnow  
   real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso 


   real(r8),   parameter     :: c2 = 21.e-3   
   real(r8),   parameter     :: c3 = 2.5e-6   
   real(r8),   parameter     :: c4 = 0.04     
   real(r8),   parameter     :: c5 = 2.0      
   real(r8),   parameter     :: dm = 100.0    
   real(r8),   parameter     :: eta0 = 0.8e+6 
                                        
   real(r8)                      :: burden 
   real(r8)                      :: ddz1   
   real(r8)                      :: ddz2   
   real(r8)                      :: ddz3   
   real(r8)                      :: dexpf  
   real(r8)                      :: td     
   real(r8)                      :: pdzdtc 
   real(r8)                      :: void   
   real(r8)                      :: wx     
   real(r8)                      :: bi     
   real(r8),   dimension(-nsnow+1:0) :: fice   

   integer  :: j


    burden = 0.0

    do j = isnow+1, 0

        wx      = snice(j) + snliq(j)
        fice(j) = snice(j) / wx
        void    = 1. - (snice(j)/denice + snliq(j)/denh2o) / dzsnso(j)

        
        if (void > 0.001 .and. snice(j) > 0.1) then
           bi = snice(j) / dzsnso(j)
           td = max(0.,tfrz-stc(j))
           dexpf = exp(-c4*td)

           

           ddz1 = -c3*dexpf

           if (bi > dm) ddz1 = ddz1*exp(-46.0e-3*(bi-dm))

           

           if (snliq(j) > 0.01*dzsnso(j)) ddz1=ddz1*c5

           

           ddz2 = -(burden+0.5*wx)*exp(-0.08*td-c2*bi)/eta0 

           

           if (imelt(j) == 1) then
              ddz3 = max(0.,(ficeold(j) - fice(j))/max(1.e-6,ficeold(j)))
              ddz3 = - ddz3/dt           
           else
              ddz3 = 0.
           end if

           

           pdzdtc = (ddz1 + ddz2 + ddz3)*dt
           pdzdtc = max(-0.5,pdzdtc)

           

           dzsnso(j) = dzsnso(j)*(1.+pdzdtc)
        end if

        

        burden = burden + wx

    end do

  end subroutine compact_glacier

  subroutine combine_glacier (nsnow  ,nsoil  ,                         & 
                              isnow  ,sh2o   ,stc    ,snice  ,snliq  , & 
                              dzsnso ,sice   ,snowh  ,sneqv  ,         & 
                              ponding1       ,ponding2)                  

    implicit none



    integer, intent(in)     :: nsnow                        
    integer, intent(in)     :: nsoil                        



    integer,                         intent(inout) :: isnow 
    real(r8),   dimension(       1:nsoil), intent(inout) :: sh2o  
    real(r8),   dimension(       1:nsoil), intent(inout) :: sice  
    real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: stc   
    real(r8),   dimension(-nsnow+1:    0), intent(inout) :: snice 
    real(r8),   dimension(-nsnow+1:    0), intent(inout) :: snliq 
    real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso
    real(r8),                              intent(inout) :: sneqv 
    real(r8),                              intent(inout) :: snowh 
    real(r8),                              intent(inout) :: ponding1
    real(r8),                              intent(inout) :: ponding2



    integer :: i,j,k,l               
    integer :: isnow_old             
    integer :: mssi                  
    integer :: neibor                
    real(r8)                      :: zwice                 
    real(r8)                      :: zwliq                 
    real(r8)                      :: dzmin(3)              
    data dzmin /0.045, 0.05, 0.2/



       isnow_old = isnow

       do j = isnow_old+1,0
          if (snice(j) <= .1) then
             if(j /= 0) then
                snliq(j+1) = snliq(j+1) + snliq(j)
                snice(j+1) = snice(j+1) + snice(j)
             else
               if (isnow_old < -1) then
                snliq(j-1) = snliq(j-1) + snliq(j)
                snice(j-1) = snice(j-1) + snice(j)
               else
                ponding1 = ponding1 + snliq(j)       
                sneqv = snice(j)                     
                snowh = dzsnso(j)                    
                snliq(j) = 0.0                       
                snice(j) = 0.0                       
                dzsnso(j) = 0.0
               endif


             endif

             
             if (j > isnow+1 .and. isnow < -1) then
                do i = j, isnow+2, -1
                   stc(i)   = stc(i-1)
                   snliq(i) = snliq(i-1)
                   snice(i) = snice(i-1)
                   dzsnso(i)= dzsnso(i-1)
                end do
             end if
             isnow = isnow + 1
          end if
       end do



       if(sice(1) < 0.) then
          sh2o(1) = sh2o(1) + sice(1)
          sice(1) = 0.
       end if

       if(isnow ==0) return   

       sneqv  = 0.
       snowh  = 0.
       zwice  = 0.
       zwliq  = 0.

       do j = isnow+1,0
             sneqv = sneqv + snice(j) + snliq(j)
             snowh = snowh + dzsnso(j)
             zwice = zwice + snice(j)
             zwliq = zwliq + snliq(j)
       end do





       if (snowh < 0.05 .and. isnow < 0 ) then
          isnow  = 0
          sneqv = zwice
          ponding2 = ponding2 + zwliq           
          if(sneqv <= 0.) snowh = 0.            
       end if










       if (isnow < -1) then

          isnow_old = isnow
          mssi     = 1

          do i = isnow_old+1,0
             if (dzsnso(i) < dzmin(mssi)) then

                if (i == isnow+1) then
                   neibor = i + 1
                else if (i == 0) then
                   neibor = i - 1
                else
                   neibor = i + 1
                   if ((dzsnso(i-1)+dzsnso(i)) < (dzsnso(i+1)+dzsnso(i))) neibor = i-1
                end if

                
                if (neibor > i) then
                   j = neibor
                   l = i
                else
                   j = i
                   l = neibor
                end if

                call combo_glacier (dzsnso(j), snliq(j), snice(j), &
                   stc(j), dzsnso(l), snliq(l), snice(l), stc(l) )

                
                if (j-1 > isnow+1) then
                   do k = j-1, isnow+2, -1
                      stc(k)   = stc(k-1)
                      snice(k) = snice(k-1)
                      snliq(k) = snliq(k-1)
                      dzsnso(k) = dzsnso(k-1)
                   end do
                end if

                
                isnow = isnow + 1
                if (isnow >= -1) exit
             else

                
                mssi = mssi + 1

             end if
          end do

       end if

  end subroutine combine_glacier



  subroutine combo_glacier(dz,  wliq,  wice, t, dz2, wliq2, wice2, t2)

    implicit none





    real(r8),   intent(in)    :: dz2   
    real(r8),   intent(in)    :: wliq2 
    real(r8),   intent(in)    :: wice2 
    real(r8),   intent(in)    :: t2    
    real(r8),   intent(inout) :: dz    
    real(r8),   intent(inout) :: wliq  
    real(r8),   intent(inout) :: wice  
    real(r8),   intent(inout) :: t     



    real                :: dzc   
    real                :: wliqc 
    real                :: wicec 
    real                :: tc    
    real                :: h     
    real                :: h2    
    real                :: hc    



    dzc = dz+dz2
    wicec = (wice+wice2)
    wliqc = (wliq+wliq2)
    h = (cice*wice+cwat*wliq) * (t-tfrz)+hfus*wliq
    h2= (cice*wice2+cwat*wliq2) * (t2-tfrz)+hfus*wliq2

    hc = h + h2
    if(hc < 0.)then
       tc = tfrz + hc/(cice*wicec + cwat*wliqc)
    else if (hc.le.hfus*wliqc) then
       tc = tfrz
    else
       tc = tfrz + (hc - hfus*wliqc) / (cice*wicec + cwat*wliqc)
    end if

    dz = dzc
    wice = wicec
    wliq = wliqc
    t = tc

  end subroutine combo_glacier

  subroutine divide_glacier (nsnow  ,nsoil  ,                         & 
                             isnow  ,stc    ,snice  ,snliq  ,dzsnso  )  

    implicit none



    integer, intent(in)                            :: nsnow 
    integer, intent(in)                            :: nsoil 



    integer                        , intent(inout) :: isnow 
    real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: stc   
    real(r8),   dimension(-nsnow+1:    0), intent(inout) :: snice 
    real(r8),   dimension(-nsnow+1:    0), intent(inout) :: snliq 
    real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso



    integer                                        :: j     
    integer                                        :: msno  
    real(r8)  :: drr   
    real(r8),   dimension(       1:nsnow)                :: dz    
    real(r8),   dimension(       1:nsnow)                :: swice 
    real(r8),   dimension(       1:nsnow)                :: swliq 
    real(r8),   dimension(       1:nsnow)                :: tsno  
    real(r8)  :: zwice 
    real(r8)  :: zwliq 
    real(r8)  :: propor
    real(r8)  :: dtdz  


    do j = 1,nsnow
          if (j <= abs(isnow)) then
             dz(j)    = dzsnso(j+isnow)
             swice(j) = snice(j+isnow)
             swliq(j) = snliq(j+isnow)
             tsno(j)  = stc(j+isnow)
          end if
    end do

       msno = abs(isnow)

       if (msno == 1) then
          
          if (dz(1) > 0.05) then
             msno = 2
             dz(1)    = dz(1)/2.
             swice(1) = swice(1)/2.
             swliq(1) = swliq(1)/2.
             dz(2)    = dz(1)
             swice(2) = swice(1)
             swliq(2) = swliq(1)
             tsno(2)  = tsno(1)
          end if
       end if

       if (msno > 1) then
          if (dz(1) > 0.05) then
             drr      = dz(1) - 0.05
             propor   = drr/dz(1)
             zwice    = propor*swice(1)
             zwliq    = propor*swliq(1)
             propor   = 0.05/dz(1)
             swice(1) = propor*swice(1)
             swliq(1) = propor*swliq(1)
             dz(1)    = 0.05

             call combo_glacier (dz(2), swliq(2), swice(2), tsno(2), drr, &
                  zwliq, zwice, tsno(1))

             

             if (msno <= 2 .and. dz(2) > 0.10) then
                msno = 3
                dtdz = (tsno(1) - tsno(2))/((dz(1)+dz(2))/2.)
                dz(2)    = dz(2)/2.
                swice(2) = swice(2)/2.
                swliq(2) = swliq(2)/2.
                dz(3)    = dz(2)
                swice(3) = swice(2)
                swliq(3) = swliq(2)
                tsno(3) = tsno(2) - dtdz*dz(2)/2.
                if (tsno(3) >= tfrz) then
                   tsno(3)  = tsno(2)
                else
                   tsno(2) = tsno(2) + dtdz*dz(2)/2.
                endif

             end if
          end if
       end if

       if (msno > 2) then
          if (dz(2) > 0.2) then
             drr = dz(2) - 0.2
             propor   = drr/dz(2)
             zwice    = propor*swice(2)
             zwliq    = propor*swliq(2)
             propor   = 0.2/dz(2)
             swice(2) = propor*swice(2)
             swliq(2) = propor*swliq(2)
             dz(2)    = 0.2
             call combo_glacier (dz(3), swliq(3), swice(3), tsno(3), drr, &
                  zwliq, zwice, tsno(2))
          end if
       end if

       isnow = -msno

    do j = isnow+1,0
             dzsnso(j) = dz(j-isnow)
             snice(j) = swice(j-isnow)
             snliq(j) = swliq(j-isnow)
             stc(j)   = tsno(j-isnow)
    end do






  end subroutine divide_glacier

  subroutine snowh2o_glacier (nsnow  ,nsoil  ,dt     ,qsnfro ,qsnsub , & 
                              qrain  ,                                 & 
                              isnow  ,dzsnso ,snowh  ,sneqv  ,snice  , & 
                              snliq  ,sh2o   ,sice   ,stc    ,         & 
                              ponding1       ,ponding2       ,fsh    , & 
                              qsnbot )                                   




   implicit none



   integer,                         intent(in)    :: nsnow  
   integer,                         intent(in)    :: nsoil  
   real(r8),                              intent(in)    :: dt     
   real(r8),                              intent(inout)    :: qsnfro 
   real(r8),                              intent(inout)    :: qsnsub 
   real(r8),                              intent(in)    :: qrain  



   real(r8),                              intent(out)   :: qsnbot 



   integer,                         intent(inout) :: isnow  
   real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso 
   real(r8),                              intent(inout) :: snowh  
   real(r8),                              intent(inout) :: sneqv  
   real(r8),   dimension(-nsnow+1:0),     intent(inout) :: snice  
   real(r8),   dimension(-nsnow+1:0),     intent(inout) :: snliq  
   real(r8),   dimension(       1:nsoil), intent(inout) :: sh2o   
   real(r8),   dimension(       1:nsoil), intent(inout) :: sice   
   real(r8),   dimension(-nsnow+1:nsoil), intent(inout) :: stc    
   real(r8),                              intent(inout) :: ponding1
   real(r8),                              intent(inout) :: ponding2
   real(r8),                              intent(inout) :: fsh     



   integer                     :: j         
   real                        :: qin       
   real                        :: qout      
   real                        :: wgdif     
   real(r8),   dimension(-nsnow+1:0) :: vol_liq   
   real(r8),   dimension(-nsnow+1:0) :: vol_ice   
   real(r8),   dimension(-nsnow+1:0) :: epore     
   real(r8)                      :: propor, temp




   if(sneqv == 0.) then
     if(opt_gla == 1) then
       sice(1) =  sice(1) + (qsnfro-qsnsub)*dt/(dzsnso(1)*1000.)
     elseif(opt_gla == 2) then
       fsh = fsh - (qsnfro-qsnsub)*hsub
       qsnfro = 0.0
       qsnsub = 0.0
     end if
   end if






   if(isnow == 0 .and. sneqv > 0.) then
      if(opt_gla == 1) then
        temp   = sneqv
        sneqv  = sneqv - qsnsub*dt + qsnfro*dt
        propor = sneqv/temp
        snowh  = max(0.,propor * snowh)
      elseif(opt_gla == 2) then
        fsh = fsh - (qsnfro-qsnsub)*hsub
        qsnfro = 0.0
        qsnsub = 0.0
      end if

      if(sneqv < 0.) then
         sice(1) = sice(1) + sneqv/(dzsnso(1)*1000.)
         sneqv   = 0.
         snowh   = 0.
      end if
      if(sice(1) < 0.) then
         sh2o(1) = sh2o(1) + sice(1)
         sice(1) = 0.
      end if
   end if

   if(snowh <= 1.e-8 .or. sneqv <= 1.e-6) then
     snowh = 0.0
     sneqv = 0.0
   end if



   if ( isnow < 0 ) then 

      wgdif = snice(isnow+1) - qsnsub*dt + qsnfro*dt
      snice(isnow+1) = wgdif
      if (wgdif < 1.e-6 .and. isnow <0) then
         call  combine_glacier (nsnow  ,nsoil  ,                         & 
                                isnow  ,sh2o   ,stc    ,snice  ,snliq  , & 
                                dzsnso ,sice   ,snowh  ,sneqv  ,         & 
                               ponding1, ponding2 )                        
      endif
      
      if ( isnow < 0 ) then 
         snliq(isnow+1) = snliq(isnow+1) + qrain * dt
         snliq(isnow+1) = max(0., snliq(isnow+1))
      endif
      
   endif 



   

   do j = -nsnow+1, 0
      if (j >= isnow+1) then
         vol_ice(j)      = min(1., snice(j)/(dzsnso(j)*denice))
         epore(j)        = 1. - vol_ice(j)
         vol_liq(j)      = min(epore(j),snliq(j)/(dzsnso(j)*denh2o))
      end if
   end do

   qin = 0.
   qout = 0.

   

   do j = -nsnow+1, 0
      if (j >= isnow+1) then
         snliq(j) = snliq(j) + qin
         if (j <= -1) then
            if (epore(j) < 0.05 .or. epore(j+1) < 0.05) then
               qout = 0.
            else
               qout = max(0.,(vol_liq(j)-ssi*epore(j))*dzsnso(j))
               qout = min(qout,(1.-vol_ice(j+1)-vol_liq(j+1))*dzsnso(j+1))
            end if
         else
            qout = max(0.,(vol_liq(j) - ssi*epore(j))*dzsnso(j))
         end if
         qout = qout*1000.
         snliq(j) = snliq(j) - qout
         qin = qout
      end if
   end do



   qsnbot = qout / dt           

  end subroutine snowh2o_glacier


  subroutine error_glacier (iloc   ,jloc   ,swdown ,fsa    ,fsr    ,fira   , &
                            fsh    ,fgev   ,ssoil  ,sag    ,prcp   ,edir   , &
		            runsrf ,runsub ,sneqv  ,dt     ,beg_wb )



  implicit none


  integer                        , intent(in) :: iloc   
  integer                        , intent(in) :: jloc   
  real(r8)                   , intent(in) :: swdown 
  real(r8)                   , intent(in) :: fsa    
  real(r8)                   , intent(in) :: fsr    
  real(r8)                   , intent(in) :: fira   
  real(r8)                   , intent(in) :: fsh    
  real(r8)                   , intent(in) :: fgev   
  real(r8)                   , intent(in) :: ssoil  
  real(r8)                   , intent(in) :: sag

  real(r8)                   , intent(in) :: prcp   
  real(r8)                   , intent(in) :: edir   
  real(r8)                   , intent(in) :: runsrf 
  real(r8)                   , intent(in) :: runsub 
  real(r8)                   , intent(in) :: sneqv  
  real(r8)                   , intent(in) :: dt     
  real(r8)                   , intent(in) :: beg_wb 

  real(r8)                              :: end_wb 
  real(r8)                              :: errwat 
  real(r8)                              :: erreng 
  real(r8)                              :: errsw  
  character(len=256)                          :: message

   errsw   = swdown - (fsa + fsr)
   if (errsw > 0.01) then            
     write(*,*) "sag    =",sag
     write(*,*) "fsa    =",fsa
     write(*,*) "fsr    =",fsr
     write(message,*) 'errsw =',errsw
     print*, (trim(message))
     print*, (&
"radiation budget problem in noahmp glacier")
   end if

   erreng = sag-(fira+fsh+fgev+ssoil)
   if(erreng > 0.01) then
      write(message,*) 'erreng =',erreng
      print*, (trim(message))
      write(message,'(i6,1x,i6,1x,5f10.4)')iloc,jloc,sag,fira,fsh,fgev,ssoil
      print*, (trim(message))
      print*, (&
"energy budget problem in noahmp glacier")
   end if

   end_wb = sneqv
   errwat = end_wb-beg_wb-(prcp-edir-runsrf-runsub)*dt

   if(abs(errwat) > 0.1) then
      if (errwat > 0) then
         print*,  ('the model is gaining water (errwat is positive)')
      else
         print*, ('the model is losing water (errwat is negative)')
      endif
      write(message, *) 'errwat =',errwat, "kg m{-2} timestep{-1}"
      print*, (trim(message))
      write(message,'("    i      j     end_wb     beg_wb       prcp       edir      runsrf     runsub")')
           print*, (trim(message))
           write(message,'(i6,1x,i6,1x,2f15.3,4f11.5)')iloc,jloc,end_wb,beg_wb,prcp*dt,&
                edir*dt,runsrf*dt,runsub*dt
           print*, (trim(message))
           print*, (&
"water budget problem in noahmp glacier")
        end if

 end subroutine error_glacier


  subroutine noahmp_options_glacier(iopt_alb  ,iopt_snf  ,iopt_tbot, iopt_stc, iopt_gla )

  implicit none

  integer,  intent(in) :: iopt_alb  
  integer,  intent(in) :: iopt_snf  
  integer,  intent(in) :: iopt_tbot 
  integer,  intent(in) :: iopt_stc  
                                    
  integer,  intent(in) :: iopt_gla  



  opt_alb  = iopt_alb  
  opt_snf  = iopt_snf  
  opt_tbot = iopt_tbot 
  opt_stc  = iopt_stc
  opt_gla  = iopt_gla
  
  end subroutine noahmp_options_glacier
 
end module noahmp_glacier_routines


module module_sf_noahmp_glacier

  use noahmp_glacier_routines
  use noahmp_glacier_globals

end module module_sf_noahmp_glacier
