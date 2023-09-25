
!----------------------------------------------------------------------------
! Created on 2020
! Author: Yi Zhang
! Version 1.0
! Description: This module contains vars needed for restarting LSM,
!              might be over-inserted, may need tune
! Revision history: 
!----------------------------------------------------------------------------

 module grist_lsm_noahmp_resVars

   use grist_constants,       only: r8, i4, zero
   use grist_data_types,      only: scalar_1d_field, scalar_2d_field
   use grist_lsm_noahmp_vars
 
   implicit none
   private

   public :: pstate_lsm
   public :: grist_lsm_resVars_construct, &
             grist_lsm_resVars_destruct , &
             fill_lsm_resVars_for_write , &
             fill_lsm_resVars_after_read

    type physics_lsm_resvars
! 59 vars
         type(scalar_2d_field) :: tslb
         type(scalar_2d_field) :: smois
         type(scalar_2d_field) :: smoiseq
         type(scalar_2d_field) :: sh2o
         type(scalar_2d_field) :: tsnoxy
         type(scalar_2d_field) :: zsnsoxy
         type(scalar_2d_field) :: snicexy
         type(scalar_2d_field) :: snliqxy

         type(scalar_1d_field) :: tsk     
         type(scalar_1d_field) :: hfx
         type(scalar_1d_field) :: qfx
         type(scalar_1d_field) :: lh
         type(scalar_1d_field) :: grdflx
         type(scalar_1d_field) :: smstav
         type(scalar_1d_field) :: smstot
         type(scalar_1d_field) :: sfcrunoff 
         type(scalar_1d_field) :: udrunoff  
         type(scalar_1d_field) :: albedo    
         type(scalar_1d_field) :: snowc
         type(scalar_1d_field) :: snow
         type(scalar_1d_field) :: snowh
         type(scalar_1d_field) :: canwat
         type(scalar_1d_field) :: acsnom
         type(scalar_1d_field) :: acsnow
         type(scalar_1d_field) :: emiss
         type(scalar_1d_field) :: qsfc
         type(scalar_1d_field) :: z0
         type(scalar_1d_field) :: znt
         type(scalar_1d_field) :: isnowxy   ! must be integer, needs converting
         type(scalar_1d_field) :: tvxy
         type(scalar_1d_field) :: tgxy
         type(scalar_1d_field) :: canicexy
         type(scalar_1d_field) :: canliqxy
         type(scalar_1d_field) :: eahxy
         type(scalar_1d_field) :: tahxy
         type(scalar_1d_field) :: cmxy
         type(scalar_1d_field) :: chxy
         type(scalar_1d_field) :: fwetxy
         type(scalar_1d_field) :: sneqvoxy
         type(scalar_1d_field) :: alboldxy
         type(scalar_1d_field) :: qsnowxy
         type(scalar_1d_field) :: wslakexy
         type(scalar_1d_field) :: zwtxy
         type(scalar_1d_field) :: waxy
         type(scalar_1d_field) :: wtxy
         type(scalar_1d_field) :: lfmassxy  
         type(scalar_1d_field) :: rtmassxy  
         type(scalar_1d_field) :: stmassxy  
         type(scalar_1d_field) :: woodxy    
         type(scalar_1d_field) :: grainxy   
         type(scalar_1d_field) :: gddxy     
         type(scalar_1d_field) :: stblcpxy  
         type(scalar_1d_field) :: fastcpxy  
         type(scalar_1d_field) :: xlaixy    
         type(scalar_1d_field) :: xsaixy    
         type(scalar_1d_field) :: taussxy
         type(scalar_1d_field) :: smcwtdxy  
         type(scalar_1d_field) :: deeprechxy 
         type(scalar_1d_field) :: rechxy
    
    end type physics_lsm_resvars

    type(physics_lsm_resvars)    :: pstate_lsm

  contains

    subroutine grist_lsm_resVars_construct(ncell)
      integer(i4),  intent(in)  :: ncell
! 2d var
        call wrap_allocate_data2d(ncell,nsoil,pstate_lsm%tslb)  
        call wrap_allocate_data2d(ncell,nsoil,pstate_lsm%smois)
        call wrap_allocate_data2d(ncell,nsoil,pstate_lsm%smoiseq)
        call wrap_allocate_data2d(ncell,nsoil,pstate_lsm%sh2o)
        call wrap_allocate_data2d(ncell,3    ,pstate_lsm%tsnoxy)
        call wrap_allocate_data2d(ncell,7    ,pstate_lsm%zsnsoxy)
        call wrap_allocate_data2d(ncell,3    ,pstate_lsm%snicexy)
        call wrap_allocate_data2d(ncell,3    ,pstate_lsm%snliqxy)
! 1d var
        call wrap_allocate_data1d(ncell, pstate_lsm%tsk)   
        call wrap_allocate_data1d(ncell, pstate_lsm%hfx)
        call wrap_allocate_data1d(ncell, pstate_lsm%qfx)
        call wrap_allocate_data1d(ncell, pstate_lsm%lh)
        call wrap_allocate_data1d(ncell, pstate_lsm%grdflx)
        call wrap_allocate_data1d(ncell, pstate_lsm%smstav)
        call wrap_allocate_data1d(ncell, pstate_lsm%smstot)
        call wrap_allocate_data1d(ncell, pstate_lsm%sfcrunoff) 
        call wrap_allocate_data1d(ncell, pstate_lsm%udrunoff)  
        call wrap_allocate_data1d(ncell, pstate_lsm%albedo)   
        call wrap_allocate_data1d(ncell, pstate_lsm%snowc)
        call wrap_allocate_data1d(ncell, pstate_lsm%snow)
        call wrap_allocate_data1d(ncell, pstate_lsm%snowh)
        call wrap_allocate_data1d(ncell, pstate_lsm%canwat)
        call wrap_allocate_data1d(ncell, pstate_lsm%acsnom)
        call wrap_allocate_data1d(ncell, pstate_lsm%acsnow)
        call wrap_allocate_data1d(ncell, pstate_lsm%emiss)
        call wrap_allocate_data1d(ncell, pstate_lsm%qsfc)
        call wrap_allocate_data1d(ncell, pstate_lsm%z0)
        call wrap_allocate_data1d(ncell, pstate_lsm%znt)
        call wrap_allocate_data1d(ncell, pstate_lsm%isnowxy)   ! must be integer, needs converting
        call wrap_allocate_data1d(ncell, pstate_lsm%tvxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%tgxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%canicexy)
        call wrap_allocate_data1d(ncell, pstate_lsm%canliqxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%eahxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%tahxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%cmxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%chxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%fwetxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%sneqvoxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%alboldxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%qsnowxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%wslakexy)
        call wrap_allocate_data1d(ncell, pstate_lsm%zwtxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%waxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%wtxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%lfmassxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%rtmassxy) 
        call wrap_allocate_data1d(ncell, pstate_lsm%stmassxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%woodxy) 
        call wrap_allocate_data1d(ncell, pstate_lsm%grainxy)   
        call wrap_allocate_data1d(ncell, pstate_lsm%gddxy)  
        call wrap_allocate_data1d(ncell, pstate_lsm%stblcpxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%fastcpxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%xlaixy)
        call wrap_allocate_data1d(ncell, pstate_lsm%xsaixy)
        call wrap_allocate_data1d(ncell, pstate_lsm%taussxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%smcwtdxy)  
        call wrap_allocate_data1d(ncell, pstate_lsm%deeprechxy)
        call wrap_allocate_data1d(ncell, pstate_lsm%rechxy)

      return
    end subroutine grist_lsm_resVars_construct

    subroutine grist_lsm_resVars_destruct()
! 2d var
        call wrap_deallocate_data2d(pstate_lsm%tslb)  
        call wrap_deallocate_data2d(pstate_lsm%smois)
        call wrap_deallocate_data2d(pstate_lsm%smoiseq)
        call wrap_deallocate_data2d(pstate_lsm%sh2o)
        call wrap_deallocate_data2d(pstate_lsm%tsnoxy)
        call wrap_deallocate_data2d(pstate_lsm%zsnsoxy)
        call wrap_deallocate_data2d(pstate_lsm%snicexy)
        call wrap_deallocate_data2d(pstate_lsm%snliqxy)
! 1d var
        call wrap_deallocate_data1d(pstate_lsm%tsk)   
        call wrap_deallocate_data1d(pstate_lsm%hfx)
        call wrap_deallocate_data1d(pstate_lsm%qfx)
        call wrap_deallocate_data1d(pstate_lsm%lh)
        call wrap_deallocate_data1d(pstate_lsm%grdflx)
        call wrap_deallocate_data1d(pstate_lsm%smstav)
        call wrap_deallocate_data1d(pstate_lsm%smstot)
        call wrap_deallocate_data1d(pstate_lsm%sfcrunoff) 
        call wrap_deallocate_data1d(pstate_lsm%udrunoff)  
        call wrap_deallocate_data1d(pstate_lsm%albedo)   
        call wrap_deallocate_data1d(pstate_lsm%snowc)
        call wrap_deallocate_data1d(pstate_lsm%snow)
        call wrap_deallocate_data1d(pstate_lsm%snowh)
        call wrap_deallocate_data1d(pstate_lsm%canwat)
        call wrap_deallocate_data1d(pstate_lsm%acsnom)
        call wrap_deallocate_data1d(pstate_lsm%acsnow)
        call wrap_deallocate_data1d(pstate_lsm%emiss)
        call wrap_deallocate_data1d(pstate_lsm%qsfc)
        call wrap_deallocate_data1d(pstate_lsm%z0)
        call wrap_deallocate_data1d(pstate_lsm%znt)
        call wrap_deallocate_data1d(pstate_lsm%isnowxy)   ! must be integer, needs converting
        call wrap_deallocate_data1d(pstate_lsm%tvxy)
        call wrap_deallocate_data1d(pstate_lsm%tgxy)
        call wrap_deallocate_data1d(pstate_lsm%canicexy)
        call wrap_deallocate_data1d(pstate_lsm%canliqxy)
        call wrap_deallocate_data1d(pstate_lsm%eahxy)
        call wrap_deallocate_data1d(pstate_lsm%tahxy)
        call wrap_deallocate_data1d(pstate_lsm%cmxy)
        call wrap_deallocate_data1d(pstate_lsm%chxy)
        call wrap_deallocate_data1d(pstate_lsm%fwetxy)
        call wrap_deallocate_data1d(pstate_lsm%sneqvoxy)
        call wrap_deallocate_data1d(pstate_lsm%alboldxy)
        call wrap_deallocate_data1d(pstate_lsm%qsnowxy)
        call wrap_deallocate_data1d(pstate_lsm%wslakexy)
        call wrap_deallocate_data1d(pstate_lsm%zwtxy)
        call wrap_deallocate_data1d(pstate_lsm%waxy)
        call wrap_deallocate_data1d(pstate_lsm%wtxy)
        call wrap_deallocate_data1d(pstate_lsm%lfmassxy)
        call wrap_deallocate_data1d(pstate_lsm%rtmassxy) 
        call wrap_deallocate_data1d(pstate_lsm%stmassxy)
        call wrap_deallocate_data1d(pstate_lsm%woodxy) 
        call wrap_deallocate_data1d(pstate_lsm%grainxy)   
        call wrap_deallocate_data1d(pstate_lsm%gddxy)  
        call wrap_deallocate_data1d(pstate_lsm%stblcpxy)
        call wrap_deallocate_data1d(pstate_lsm%fastcpxy)
        call wrap_deallocate_data1d(pstate_lsm%xlaixy)
        call wrap_deallocate_data1d(pstate_lsm%xsaixy)
        call wrap_deallocate_data1d(pstate_lsm%taussxy)
        call wrap_deallocate_data1d(pstate_lsm%smcwtdxy)  
        call wrap_deallocate_data1d(pstate_lsm%deeprechxy)
        call wrap_deallocate_data1d(pstate_lsm%rechxy)

      return
    end subroutine grist_lsm_resVars_destruct

    subroutine fill_lsm_resVars_for_write(ncol)
       integer(i4) ,  intent(in) :: ncol       ! not full level but as in lsm-len(i.e., halo(1))

! 2d var
        pstate_lsm%tslb%f   (1:nsoil,1:ncol) = tslb   (1:nsoil, 1:ncol)
        pstate_lsm%smois%f  (1:nsoil,1:ncol) = smois  (1:nsoil, 1:ncol)
        pstate_lsm%smoiseq%f(1:nsoil,1:ncol) = smoiseq(1:nsoil, 1:ncol)
        pstate_lsm%sh2o%f   (1:nsoil,1:ncol) = sh2o   (1:nsoil, 1:ncol)
        pstate_lsm%tsnoxy%f (1:3,1:ncol)     = tsnoxy (-2:0    ,1:ncol)
        pstate_lsm%zsnsoxy%f(1:7,1:ncol)     = zsnsoxy(-2:nsoil,1:ncol)
        pstate_lsm%snicexy%f(1:3,1:ncol)     = snicexy(-2:0    ,1:ncol)
        pstate_lsm%snliqxy%f(1:3,1:ncol)     = snliqxy(-2:0    ,1:ncol)
! 1d var
        pstate_lsm%tsk%f(1:ncol)             = tsk(1:ncol)
        pstate_lsm%hfx%f(1:ncol)             = hfx(1:ncol)
        pstate_lsm%qfx%f(1:ncol)             = qfx(1:ncol)
        pstate_lsm%lh%f(1:ncol)              = lh(1:ncol)
        pstate_lsm%grdflx%f(1:ncol)          = grdflx(1:ncol)
        pstate_lsm%smstav%f(1:ncol)          = smstav(1:ncol)
        pstate_lsm%smstot%f(1:ncol)          = smstot(1:ncol)
        pstate_lsm%sfcrunoff%f(1:ncol)       = sfcrunoff(1:ncol)
        pstate_lsm%udrunoff%f(1:ncol)        = udrunoff(1:ncol)
        pstate_lsm%albedo%f(1:ncol)          = albedo(1:ncol)
        pstate_lsm%snowc%f(1:ncol)           = snowc(1:ncol)
        pstate_lsm%snow%f(1:ncol)            = snow(1:ncol)
        pstate_lsm%snowh%f(1:ncol)           = snowh(1:ncol)
        pstate_lsm%canwat%f(1:ncol)          = canwat(1:ncol)
        pstate_lsm%acsnom%f(1:ncol)          = acsnom(1:ncol)
        pstate_lsm%acsnow%f(1:ncol)          = acsnow(1:ncol)
        pstate_lsm%emiss%f(1:ncol)           = emiss(1:ncol)
        pstate_lsm%qsfc%f(1:ncol)            = qsfc(1:ncol)
        pstate_lsm%z0%f(1:ncol)              = z0(1:ncol)
        pstate_lsm%znt%f(1:ncol)             = znt(1:ncol)
        pstate_lsm%isnowxy%f(1:ncol)         = isnowxy(1:ncol)
        pstate_lsm%tvxy%f(1:ncol)            = tvxy(1:ncol)
        pstate_lsm%tgxy%f(1:ncol)            = tgxy(1:ncol)
        pstate_lsm%canicexy%f(1:ncol)        = canicexy(1:ncol)
        pstate_lsm%canliqxy%f(1:ncol)        = canliqxy(1:ncol)
        pstate_lsm%eahxy%f(1:ncol)           = eahxy(1:ncol)
        pstate_lsm%tahxy%f(1:ncol)           = tahxy(1:ncol)
        pstate_lsm%cmxy%f(1:ncol)            = cmxy(1:ncol)
        pstate_lsm%chxy%f(1:ncol)            = chxy(1:ncol)
        pstate_lsm%fwetxy%f(1:ncol)          = fwetxy(1:ncol)
        pstate_lsm%sneqvoxy%f(1:ncol)        = sneqvoxy(1:ncol)
        pstate_lsm%alboldxy%f(1:ncol)        = alboldxy(1:ncol)
        pstate_lsm%qsnowxy%f(1:ncol)         = qsnowxy(1:ncol)
        pstate_lsm%wslakexy%f(1:ncol)        = wslakexy(1:ncol)
        pstate_lsm%zwtxy%f(1:ncol)           = zwtxy(1:ncol)
        pstate_lsm%waxy%f(1:ncol)            = waxy(1:ncol)
        pstate_lsm%wtxy%f(1:ncol)            = wtxy(1:ncol)
        pstate_lsm%lfmassxy%f(1:ncol)        = lfmassxy(1:ncol)
        pstate_lsm%rtmassxy%f(1:ncol)        = rtmassxy(1:ncol)
        pstate_lsm%stmassxy%f(1:ncol)        = stmassxy(1:ncol)
        pstate_lsm%woodxy%f(1:ncol)          = woodxy(1:ncol)
        pstate_lsm%grainxy%f(1:ncol)         = grainxy(1:ncol)
        pstate_lsm%gddxy%f(1:ncol)           = gddxy(1:ncol)
        pstate_lsm%stblcpxy%f(1:ncol)        = stblcpxy(1:ncol)
        pstate_lsm%fastcpxy%f(1:ncol)        = fastcpxy(1:ncol)
        pstate_lsm%xlaixy%f(1:ncol)          = xlaixy(1:ncol)
        pstate_lsm%xsaixy%f(1:ncol)          = xsaixy(1:ncol)
        pstate_lsm%taussxy%f(1:ncol)         = taussxy(1:ncol)
        pstate_lsm%smcwtdxy%f(1:ncol)        = smcwtdxy(1:ncol)
        pstate_lsm%deeprechxy%f(1:ncol)      = deeprechxy(1:ncol)
        pstate_lsm%rechxy%f(1:ncol)          = rechxy(1:ncol)

        return
    end subroutine fill_lsm_resVars_for_write

    subroutine fill_lsm_resVars_after_read(ncol)
      integer(i4) ,  intent(in) :: ncol

! 2d var
         tslb   (1:nsoil, 1:ncol)  =   pstate_lsm%tslb%f   (1:nsoil,1:ncol)
         smois  (1:nsoil, 1:ncol)  =   pstate_lsm%smois%f  (1:nsoil,1:ncol)
         smoiseq(1:nsoil, 1:ncol)  =   pstate_lsm%smoiseq%f(1:nsoil,1:ncol)
         sh2o   (1:nsoil, 1:ncol)  =   pstate_lsm%sh2o%f   (1:nsoil,1:ncol)
         tsnoxy (-2:0    ,1:ncol)  =   pstate_lsm%tsnoxy%f (1:3,1:ncol)    
         zsnsoxy(-2:nsoil,1:ncol)  =   pstate_lsm%zsnsoxy%f(1:7,1:ncol)    
         snicexy(-2:0    ,1:ncol)  =   pstate_lsm%snicexy%f(1:3,1:ncol)    
         snliqxy(-2:0    ,1:ncol)  =   pstate_lsm%snliqxy%f(1:3,1:ncol)    
! 1d var                               
         tsk(1:ncol)               =   pstate_lsm%tsk%f(1:ncol)
         hfx(1:ncol)               =   pstate_lsm%hfx%f(1:ncol)
         qfx(1:ncol)               =   pstate_lsm%qfx%f(1:ncol)
         lh(1:ncol)                =   pstate_lsm%lh%f(1:ncol)
         grdflx(1:ncol)            =   pstate_lsm%grdflx%f(1:ncol)
         smstav(1:ncol)            =   pstate_lsm%smstav%f(1:ncol)
         smstot(1:ncol)            =   pstate_lsm%smstot%f(1:ncol)
         sfcrunoff(1:ncol)         =   pstate_lsm%sfcrunoff%f(1:ncol)
         udrunoff(1:ncol)          =   pstate_lsm%udrunoff%f(1:ncol)
         albedo(1:ncol)            =   pstate_lsm%albedo%f(1:ncol)
         snowc(1:ncol)             =   pstate_lsm%snowc%f(1:ncol)
         snow(1:ncol)              =   pstate_lsm%snow%f(1:ncol)
         snowh(1:ncol)             =   pstate_lsm%snowh%f(1:ncol)
         canwat(1:ncol)            =   pstate_lsm%canwat%f(1:ncol)
         acsnom(1:ncol)            =   pstate_lsm%acsnom%f(1:ncol)
         acsnow(1:ncol)            =   pstate_lsm%acsnow%f(1:ncol)
         emiss(1:ncol)             =   pstate_lsm%emiss%f(1:ncol)
         qsfc(1:ncol)              =   pstate_lsm%qsfc%f(1:ncol)
         z0(1:ncol)                =   pstate_lsm%z0%f(1:ncol)
         znt(1:ncol)               =   pstate_lsm%znt%f(1:ncol)
         isnowxy(1:ncol)           =   int(pstate_lsm%isnowxy%f(1:ncol))
         tvxy(1:ncol)              =   pstate_lsm%tvxy%f(1:ncol)
         tgxy(1:ncol)              =   pstate_lsm%tgxy%f(1:ncol)
         canicexy(1:ncol)          =   pstate_lsm%canicexy%f(1:ncol)
         canliqxy(1:ncol)          =   pstate_lsm%canliqxy%f(1:ncol)
         eahxy(1:ncol)             =   pstate_lsm%eahxy%f(1:ncol)
         tahxy(1:ncol)             =   pstate_lsm%tahxy%f(1:ncol)
         cmxy(1:ncol)              =   pstate_lsm%cmxy%f(1:ncol)
         chxy(1:ncol)              =   pstate_lsm%chxy%f(1:ncol)
         fwetxy(1:ncol)            =   pstate_lsm%fwetxy%f(1:ncol)
         sneqvoxy(1:ncol)          =   pstate_lsm%sneqvoxy%f(1:ncol)
         alboldxy(1:ncol)          =   pstate_lsm%alboldxy%f(1:ncol)
         qsnowxy(1:ncol)           =   pstate_lsm%qsnowxy%f(1:ncol)
         wslakexy(1:ncol)          =   pstate_lsm%wslakexy%f(1:ncol)
         zwtxy(1:ncol)             =   pstate_lsm%zwtxy%f(1:ncol)
         waxy(1:ncol)              =   pstate_lsm%waxy%f(1:ncol)
         wtxy(1:ncol)              =   pstate_lsm%wtxy%f(1:ncol)
         lfmassxy(1:ncol)          =   pstate_lsm%lfmassxy%f(1:ncol)
         rtmassxy(1:ncol)          =   pstate_lsm%rtmassxy%f(1:ncol)
         stmassxy(1:ncol)          =   pstate_lsm%stmassxy%f(1:ncol)
         woodxy(1:ncol)            =   pstate_lsm%woodxy%f(1:ncol)
         grainxy(1:ncol)           =   pstate_lsm%grainxy%f(1:ncol)
         gddxy(1:ncol)             =   pstate_lsm%gddxy%f(1:ncol)
         stblcpxy(1:ncol)          =   pstate_lsm%stblcpxy%f(1:ncol)
         fastcpxy(1:ncol)          =   pstate_lsm%fastcpxy%f(1:ncol)
         xlaixy(1:ncol)            =   pstate_lsm%xlaixy%f(1:ncol)
         xsaixy(1:ncol)            =   pstate_lsm%xsaixy%f(1:ncol)
         taussxy(1:ncol)           =   pstate_lsm%taussxy%f(1:ncol)
         smcwtdxy(1:ncol)          =   pstate_lsm%smcwtdxy%f(1:ncol)
         deeprechxy(1:ncol)        =   pstate_lsm%deeprechxy%f(1:ncol)
         rechxy(1:ncol)            =   pstate_lsm%rechxy%f(1:ncol)

      return
    end subroutine fill_lsm_resVars_after_read

    subroutine wrap_allocate_data1d(ncell,var)
       integer(i4),           intent(in)    :: ncell
       type(scalar_1d_field), intent(inout) :: var
       if(.not.allocated(var%f)) allocate(var%f(ncell))
       var%f    = zero
       var%pos  = 0
       return
    end subroutine wrap_allocate_data1d

    subroutine wrap_deallocate_data1d(var)
       type(scalar_1d_field), intent(inout) :: var
       if(allocated(var%f)) deallocate(var%f)
       return
    end subroutine wrap_deallocate_data1d

    subroutine wrap_allocate_data2d(ncell,nLevel,var)
       integer(i4),           intent(in)    :: ncell
       integer(i4),           intent(in)    :: nLevel
       type(scalar_2d_field), intent(inout) :: var
       if(.not.allocated(var%f)) allocate(var%f(nLevel,ncell))
       var%f    = zero
       var%pos  = 0
       return
    end subroutine wrap_allocate_data2d

    subroutine wrap_deallocate_data2d(var)
       type(scalar_2d_field), intent(inout) :: var
       if(allocated(var%f)) deallocate(var%f)
       return
    end subroutine wrap_deallocate_data2d

 end module grist_lsm_noahmp_resVars
