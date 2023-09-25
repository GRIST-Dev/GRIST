
!----------------------------------------------------------------------------
! Created on 2020
! Author: Yi Zhang
! Version 1.0
! Description:  This slab ocean module is taken from WRF-V381; we use this module
!               's scheme to modify the staticData_sst in static_data module, 
!               which will be used by PDC coupling
!               pollard, rhines and thompson (prt)
! Revision history:
!----------------------------------------------------------------------------

 module grist_slab_ocean_module

   use grist_constants,     only: i4, r8, one, zero, g=>gravity, stbolt=>stebol, omega
   use grist_handle_error,  only: showinfo
   use grist_mpi

   implicit none

     public :: grist_prt_oml_run   , &
               grist_prt_oml_init  , &
               grist_prt_oml_final

    real(r8), allocatable  :: tml(:)    ! ocean mixed layer temperature (k)
    real(r8), allocatable  :: t0ml(:)   ! ocean mixed layer temperature (k) at initial time
    real(r8), allocatable  :: hml(:)    ! ocean mixed layer depth (m)
    real(r8), allocatable  :: h0ml(:)   ! ocean mixed layer depth (m) at initial time
    real(r8), allocatable  :: huml(:)   ! ocean mixed layer u component of wind
    real(r8), allocatable  :: hvml(:)   ! ocean mixed layer v component of wind
    real(r8), allocatable  :: tmoml(:)  ! top 200 m ocean mean temperature (k) at initial time
    real(r8)               :: oml_gamma = 0.14_r8
    real(r8)               :: oml_hml0  = 50_r8   !

contains

!----------------------------------------------------------------

   subroutine grist_prt_oml_run(nv_full,ncell,nlev,dtime,lats, &
                          sstData,ust,uair,vair,xland,hfx,lh,gsw,glw,emiss)
!----------------------------------------------------------------
!----------------------------------------------------------------
!
!  subroutine oceanml calculates the sea surface temperature (tsk)
!  from a simple ocean mixed layer model based on
!  (pollard, rhines and thompson (1973).
!
!-- tml         ocean mixed layer temperature (k)
!-- t0ml        ocean mixed layer temperature (k) at initial time
!-- tmoml       top 200 m ocean mean temperature (k) at initial time
!-- hml         ocean mixed layer depth (m)
!-- h0ml        ocean mixed layer depth (m) at initial time
!-- huml        ocean mixed layer u component of wind
!-- hvml        ocean mixed layer v component of wind
!===============================================
!-- oml_gamma   deep water lapse rate (k m-1)
!-- uair,vair   lowest model level wind component
!-- ust         frictional velocity
!-- hfx         upward heat flux at the surface (w/m^2), sensible heat
!-- lh          latent heat flux at the surface (w/m^2)
!-- tsk         surface temperature (k)
!-- gsw         downward short wave flux at ground surface (w/m^2)
!-- glw         downward long wave flux at ground surface (w/m^2)
!-- emiss       emissivity of the surface
!-- xland       land mask (1 for land, 2 for water)
!-- stbolt      stefan-boltzmann constant (w/m^2/k^4)
!-- f           coriolis parameter
!-- dt          time step (second)
!-- g           acceleration due to gravity

    integer(i4)        , intent(in)    :: nv_full, ncell ! size of array and size of compute
    integer(i4)        , intent(in)    :: nlev
    real(r8)           , intent(in)    :: dtime
    real(r8)           , intent(in)    :: lats(nv_full)
    real(r8)           , intent(inout) :: sstData(nv_full)  ! return back to model
    real(r8)           , intent(in)    :: ust(nv_full),   uair(nv_full), vair(nv_full), &
                                          xland(nv_full), hfx(nv_full),  lh(nv_full), &
                                          gsw(nv_full),   glw(nv_full),  emiss(nv_full)
! local
    real(r8)     :: fff(ncell)
    integer(i4)  :: icell

    do icell = 1, ncell
       if (xland(icell).ge.1.5)then
            fff(icell)   = 2._r8*omega*sin(lats(icell))
            call oml1d(sstData(icell),t0ml(icell),hml(icell),h0ml(icell),huml(icell), &
                       !hvml(icell),sstData(icell),hfx(icell),      &
                       hvml(icell),hfx(icell),      &
                       lh(icell),gsw(icell),glw(icell),tmoml(icell),       &
                       uair(icell),vair(icell),ust(icell),fff(icell),emiss(icell),stbolt,g,dtime,oml_gamma, &
                       -10000._r8) ! no relaxation
       end if
    end do

    return
   end subroutine grist_prt_oml_run

!================================================================

   SUBROUTINE grist_prt_oml_init(oml_hml0, tsk, ncell)
                      !tml,t0ml,hml,h0ml,huml,hvml,tmoml,       &
                      !allowed_to_read, start_of_simulation,    &
                      !ids,ide, jds,jde, kds,kde,               &
                      !ims,ime, jms,jme, kms,kme,               &
                      !its,ite, jts,jte, kts,kte                )
!----------------------------------------------------------------
   IMPLICIT NONE
!----------------------------------------------------------------
   !LOGICAL , INTENT(IN)      ::      allowed_to_read
   !LOGICAL , INTENT(IN)      ::      start_of_simulation
   !INTEGER, INTENT(IN   )    ::      ids,ide, jds,jde, kds,kde, &
   !                                  ims,ime, jms,jme, kms,kme, &
   !                                  its,ite, jts,jte, kts,kte

   REAL(r8), DIMENSION(ncell), INTENT(IN) :: TSK

   !REAL(r8),    DIMENSION(ncell)                     , &
   !         INTENT(INOUT)    ::     TML, T0ML, HML, H0ML, HUML, HVML, TMOML
   REAL(r8)   , INTENT(IN)    ::     oml_hml0
   integer(i4), intent(in)    ::     ncell
!  LOCAR VAR

   INTEGER                   ::      L,J,I,itf,jtf
   CHARACTER*1024 message

    if(.not.allocated(tml))    allocate(tml(ncell))    ! ocean mixed layer temperature (k)
    if(.not.allocated(t0ml))   allocate(t0ml(ncell))   ! ocean mixed layer temperature (k) at initial time
    if(.not.allocated(hml))    allocate(hml(ncell))    ! ocean mixed layer depth (m)
    if(.not.allocated(h0ml))   allocate(h0ml(ncell))   ! ocean mixed layer depth (m) at initial time
    if(.not.allocated(huml))   allocate(huml(ncell))   ! ocean mixed layer u component of wind
    if(.not.allocated(hvml))   allocate(hvml(ncell))   ! ocean mixed layer v component of wind
    if(.not.allocated(tmoml))  allocate(tmoml(ncell))  ! top 200 m ocean mean temperature (k) at initial time

!----------------------------------------------------------------
 
   !itf=min0(ite,ide-1)
   !jtf=min0(jte,jde-1)

   !IF(start_of_simulation) THEN
     !DO J=jts,jtf
     DO I=1, ncell
       TML(I) =TSK(I)
       T0ML(I)=TSK(I)
     !ENDDO
     ENDDO
     IF (oml_hml0 .gt. 0.) THEN
        if(mpi_rank().eq.0) print*, 'Initializing OML with HML0 = ', oml_hml0
!#if (defined AMIPW_PHYSICS) || (defined AMIPC_PHYSICS)
        if(mpi_rank().eq.0) CALL showinfo(TRIM(message),0)
!#else
!        CALL wrf_debug (0, TRIM(message))
!#endif
        !DO J=jts,jtf
        DO I=1,ncell
          HML(I)=oml_hml0
          H0ML(I)=HML(I)
          HUML(I)=0.
          HVML(I)=0.
          TMOML(I)=TSK(I)-5._r8
        ENDDO
        !ENDDO
     ELSE
        if(mpi_rank().eq.0) print*, 'We should not go here for GRIST'
!#if (defined AMIPW_PHYSICS) || (defined AMIPC_PHYSICS)
        if(mpi_rank().eq.0) call showinfo(trim(message),0)
!#else
!        CALL wrf_debug (0, TRIM(message))
!#endif
        !DO J=jts,jtf
        !DO I=1,ncell
        !  HML(I)=H0ML(I)
! fill in near coast area with SST: 200 K was set as missing value in ocean pre-processing code
        !  IF(TMOML(I).GT.200. .and. TMOML(I).LE.201.) TMOML(I)=TSK(I)
        !ENDDO
        !ENDDO
     ENDIF
   !ENDIF

   END SUBROUTINE grist_prt_oml_init

   subroutine grist_prt_oml_final
      if(allocated(tml  ))   deallocate(tml  )
      if(allocated(t0ml ))   deallocate(t0ml )
      if(allocated(hml  ))   deallocate(hml  )
      if(allocated(h0ml ))   deallocate(h0ml )
      if(allocated(huml ))   deallocate(huml )
      if(allocated(hvml ))   deallocate(hvml )
      if(allocated(tmoml))   deallocate(tmoml)
      return
   end subroutine grist_prt_oml_final

!----------------------------------------------------------------

!----------------------------------------------------------------
!   SUBROUTINE OML1D(I,J,TML,T0ML,H,H0,HUML,                              &
!                    HVML,TSK,HFX,                                        &
!                    LH,GSW,GLW,TMOML,                                    &
!                    UAIR,VAIR,UST,F,EMISS,STBOLT,G,DT,OML_GAMMA,         &
!                    OML_RELAXATION_TIME,                                 &
!                    ids,ide, jds,jde, kds,kde,                           &
!                    ims,ime, jms,jme, kms,kme,                           &
!                    its,ite, jts,jte, kts,kte                            )

   SUBROUTINE OML1D(TML,T0ML,H,H0,HUML,                              &
                    !HVML,TSK,HFX,                                        &
                    HVML,HFX,                                        &
                    LH,GSW,GLW,TMOML,                                    &
                    UAIR,VAIR,UST,F,EMISS,STBOLT,G,DT,OML_GAMMA,         &
                    OML_RELAXATION_TIME  )

!----------------------------------------------------------------
   IMPLICIT NONE
!----------------------------------------------------------------
!
!  SUBROUTINE OCEANML CALCULATES THE SEA SURFACE TEMPERATURE (TSK) 
!  FROM A SIMPLE OCEAN MIXED LAYER MODEL BASED ON 
!  (Pollard, Rhines and Thompson (1973).
!
!-- TML         ocean mixed layer temperature (K)
!-- T0ML        ocean mixed layer temperature (K) at initial time
!-- TMOML       top 200 m ocean mean temperature (K) at initial time
!-- H           ocean mixed layer depth (m)
!-- H0          ocean mixed layer depth (m) at initial time
!-- HUML        ocean mixed layer u component of wind
!-- HVML        ocean mixed layer v component of wind
!-- OML_GAMMA   deep water lapse rate (K m-1)
!-- SF_OCEAN_PHYSICS     whether to call oml model
!-- UAIR,VAIR   lowest model level wind component
!-- UST         frictional velocity
!-- HFX         upward heat flux at the surface (W/m^2)
!-- LH          latent heat flux at the surface (W/m^2)
!-- TSK         surface temperature (K)
!-- GSW         downward short wave flux at ground surface (W/m^2)
!-- GLW         downward long wave flux at ground surface (W/m^2)
!-- EMISS       emissivity of the surface
!-- STBOLT      Stefan-Boltzmann constant (W/m^2/K^4)
!-- F           Coriolis parameter
!-- DT          time step (second)
!-- G           acceleration due to gravity
!-- OML_RELAXATION_TIME  time scale (s) to relax TML to T0ML, H to H0,
!                        HUML and HVML to 0; value <=0 means no relaxation
!
!----------------------------------------------------------------
   !INTEGER(i4), INTENT(IN   )    ::      I, J
   !INTEGER(i4), INTENT(IN   )    ::      ids,ide, jds,jde, kds,kde, &
   !                                  ims,ime, jms,jme, kms,kme, &
   !                                  its,ite, jts,jte, kts,kte

   REAL(r8),    INTENT(INOUT)    :: TML, H, HUML, HVML
   REAL(r8)    :: TSK

   REAL(r8),    INTENT(IN   )    :: T0ML, H0, HFX, LH, GSW, GLW,        &
                                UAIR, VAIR, UST, F, EMISS, TMOML

   REAL(r8),    INTENT(IN) :: STBOLT, G, DT, OML_GAMMA, OML_RELAXATION_TIME

! Local
   REAL(r8) :: rhoair, rhowater, Gam, alp, BV2, A1, A2, B2, u, v, wspd, &
           hu1, hv1, hu2, hv2, taux, tauy, tauxair, tauyair, q, hold, &
           hsqrd, thp, cwater, ust2
   CHARACTER(LEN=120) :: time_series

      hu1=huml
      hv1=hvml
      rhoair=1.
      rhowater=1000.
      cwater=4200.
! Deep ocean lapse rate (K/m) - from Rich
      Gam=oml_gamma
!     if(i.eq.1 .and. j.eq.1 .or. i.eq.105.and.j.eq.105) print *, 'gamma = ', gam
!     Gam=0.14
!     Gam=5.6/40.
!     Gam=5./100.
! Thermal expansion coeff (/K)
!     alp=.0002
!     temp dependence (/K)
      alp=max((tml-273.15)*1.e-5, 1.e-6)
      BV2=alp*g*Gam
      thp=t0ml-Gam*(h-h0)
      A1=(tml-thp)*h - 0.5*Gam*h*h
      if(h.ne.0.)then
        u=hu1/h
        v=hv1/h
      else
        u=0.
        v=0.
      endif

!  time step

        q=(-hfx-lh+gsw+glw*emiss-stbolt*emiss*tml*tml*tml*tml)/(rhowater*cwater)
!       wspd=max(sqrt(uair*uair+vair*vair),0.1)
        wspd=sqrt(uair*uair+vair*vair)
        if (wspd .lt. 1.e-10 ) then
!          print *, 'i,j,wspd are ', i,j,wspd
           wspd = 1.e-10
        endif
! limit ust to 1.6 to give a value of ust for water of 0.05
!       ust2=min(ust, 1.6)
! new limit for ust: reduce atmospheric ust by half for ocean
        ust2=0.5*ust
        tauxair=ust2*ust2*uair/wspd
        taux=rhoair/rhowater*tauxair
        tauyair=ust2*ust2*vair/wspd
        tauy=rhoair/rhowater*tauyair
! note: forward-backward coriolis force for effective time-centering
        hu2=hu1+dt*( f*hv1 + taux)
        hv2=hv1+dt*(-f*hu2 + tauy)
! consider the flux effect
        A2=A1+q*dt

        huml=hu2
        hvml=hv2

        hold=h
        B2=hu2*hu2+hv2*hv2
        hsqrd=-A2/Gam + sqrt(A2*A2/(Gam*Gam) + 2.*B2/BV2)
        h=sqrt(max(hsqrd,0.0))
! limit to positive h change
        if(h.lt.hold)h=hold

! no change unless tml is warmer than layer mean temp tmol or tsk-5 (see omlinit)
        if(tml.ge.tmoml .and. h.ne.0.)then

! no change unless tml is warmer than layer mean temp tmoml or tsk-5 (see omlinit)
          if(tml.ge.tmoml)then
            tml=max(t0ml - Gam*(h-h0) + 0.5*Gam*h + A2/h, tmoml)
          else 
            tml=tmoml
          endif
          u=hu2/h
          v=hv2/h
        else
          tml=t0ml
          u=0.
          v=0.
        endif

! relax TML T0ML and H to H0, HUML and HVML to 0

        if (oml_relaxation_time .gt. 0.) then
          tml = tml - (tml-t0ml)*dt/oml_relaxation_time
          h = h - (h-h0)*dt/oml_relaxation_time
          huml = huml - huml*dt/oml_relaxation_time
          hvml = hvml - hvml*dt/oml_relaxation_time
        end if

        tsk=tml
!        if(h.gt.100.)print *,i,j,h,tml,' h,tml'

! ww: output point data
!     if( (i.eq.190 .and. j.eq.115) .or. (i.eq.170 .and. j.eq.125) ) then
!        write(jtime,fmt='("TS ",f10.0)') float(itimestep)
!        CALL wrf_message ( TRIM(jtime) )
!        write(time_series,fmt='("OML",2I4,2F9.5,2F8.2,2E15.5,F8.3)') &
!              i,j,u,v,tml,h,taux,tauy,a2
!        CALL wrf_message ( TRIM(time_series) )
!     end if

   END SUBROUTINE OML1D

!================================================================

end module grist_slab_ocean_module
