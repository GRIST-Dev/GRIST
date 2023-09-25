
module module_sf_noahmp_groundwater
   use grist_constants, only: r8

contains

  subroutine wtable_mmf_noahmp (nsoil     ,xland    ,xice    ,xice_threshold  ,isice ,& 
                                isltyp    ,smoiseq  ,dzs     ,wtddt                  ,& 
                                fdepth    ,area     ,topo    ,isurban ,ivgtyp        ,& 
                                rivercond ,riverbed ,eqwtd   ,pexp                   ,& 
                                smois     ,sh2oxy   ,smcwtd  ,wtd  ,qrf              ,& 
                                deeprech  ,qspring  ,qslat   ,qrfs ,qsprings  ,rech  ,& 
                                ids,ide, jds,jde, kds,kde,                    &
                                ims,ime, jms,jme, kms,kme,                    &
                                its,ite, jts,jte, kts,kte                     )


  use noahmp_tables, only: bexp_table, dksat_table, smcmax_table,psisat_table, smcwlt_table

  implicit none


  integer,  intent(in   )     ::   ids,ide, jds,jde, kds,kde,  &
       &                           ims,ime, jms,jme, kms,kme,  &
       &                           its,ite, jts,jte, kts,kte
    real(r8)   ,   intent(in)        ::     wtddt
    real(r8)   ,   intent(in)        ::     xice_threshold
    integer,  intent(in   )   ::     isice
    real(r8)   ,    dimension( ims:ime, jms:jme )                     , &
         &   intent(in   )    ::                          xland, &
                                                           xice
    integer, dimension( ims:ime, jms:jme )                     , &
             intent(in   )    ::                         isltyp, &
                                                         ivgtyp
    integer, intent(in)       ::     nsoil
    integer, intent(in)       ::     isurban
    real(r8)   ,     dimension( ims:ime , 1:nsoil, jms:jme ), &
         &    intent(in)      ::                        smoiseq
    real(r8)   ,     dimension(1:nsoil), intent(in)     ::         dzs
    real(r8)   ,    dimension( ims:ime, jms:jme )                     , &
         &   intent(in)       ::                         fdepth, &
                                                           area, &
                                                           topo, &
                                                          eqwtd, &
                                                           pexp, &
                                                       riverbed, &
                                                      rivercond



    real(r8)   ,     dimension( ims:ime , 1:nsoil, jms:jme ), &
         &    intent(inout)   ::                          smois, &
         &                                                sh2oxy 


    real(r8)   ,    dimension( ims:ime, jms:jme )                     , &
         &   intent(inout)    ::                            wtd, &
                                                         smcwtd, &
                                                       deeprech, &
                                                          qslat, &
                                                           qrfs, &
                                                       qsprings, &
                                                           rech



    real(r8)   ,    dimension( ims:ime, jms:jme )                     , &
         &   intent(out)      ::                            qrf, &  
                                                        qspring     


  
  integer                          :: i,j,k  
  real(r8)   , dimension(       0:nsoil)  :: zsoil 
  real(r8)   ,  dimension(      1:nsoil)  :: smceq  
  real(r8)   ,  dimension(      1:nsoil)  :: smc,sh2o
  real(r8)                              :: deltat,rcond,totwater,psi &
                                                ,wfluxdeep,wcnddeep,ddz,smcwtdmid &
                                                ,wplus,wminus
  real(r8)   ,      dimension( ims:ime, jms:jme )    :: qlat
  integer,   dimension( ims:ime, jms:jme )    :: landmask 
  
  real(r8)                              :: bexp,dksat,psisat,smcmax,smcwlt

    deltat = wtddt * 60. 

    zsoil(0) = 0.
    zsoil(1) = -dzs(1)
    do k = 2, nsoil
       zsoil(k)         = -dzs(k) + zsoil(k-1)
    end do

    where(xland-1.5.lt.0..and.xice.lt. xice_threshold.and.ivgtyp.ne.isice)
         landmask=1
    elsewhere
         landmask=-1
    endwhere


    qlat = 0.
call lateralflow(isltyp,wtd,qlat,fdepth,topo,landmask,deltat,area       &
                        ,ids,ide,jds,jde,kds,kde                      &
                        ,ims,ime,jms,jme,kms,kme                      &
                        ,its,ite,jts,jte,kts,kte                      )


    do j=jts,jte
       do i=its,ite
          if(landmask(i,j).gt.0)then
             if(wtd(i,j) .gt. riverbed(i,j) .and.  eqwtd(i,j) .gt. riverbed(i,j)) then
               rcond = rivercond(i,j) * exp(pexp(i,j)*(wtd(i,j)-eqwtd(i,j)))
             else    
               rcond = rivercond(i,j)       
             endif
             qrf(i,j) = rcond * (wtd(i,j)-riverbed(i,j)) * deltat/area(i,j)

             qrf(i,j) = max(qrf(i,j),0.)
          else
             qrf(i,j) = 0.
          endif
       enddo
    enddo


    do j=jts,jte
       do i=its,ite
          if(landmask(i,j).gt.0)then

            bexp   = bexp_table   (isltyp(i,j))
            dksat  = dksat_table  (isltyp(i,j))
            psisat = psisat_table (isltyp(i,j))
            smcmax = smcmax_table (isltyp(i,j))
            smcwlt = smcwlt_table (isltyp(i,j))

             if(ivgtyp(i,j)==isurban)then
                 smcmax = 0.45
                 smcwlt = 0.40
             endif


             if(wtd(i,j) < zsoil(nsoil)-dzs(nsoil))then

                ddz = zsoil(nsoil)-wtd(i,j)
                smcwtdmid = 0.5 * (smcwtd(i,j) + smcmax )
                psi = psisat * ( smcmax / smcwtd(i,j) ) ** bexp
                wcnddeep = dksat * ( smcwtdmid / smcmax ) ** (2.0*bexp + 3.0)
                wfluxdeep =  - deltat * wcnddeep * ( (psisat-psi) / ddz - 1.)

                smcwtd(i,j) = smcwtd(i,j)  + (deeprech(i,j) -  wfluxdeep)  / ddz
                wplus       = max((smcwtd(i,j)-smcmax), 0.0) * ddz
                wminus       = max((1.e-4-smcwtd(i,j)), 0.0) * ddz
                smcwtd(i,j) = max( min(smcwtd(i,j),smcmax) , 1.e-4)
                wfluxdeep = wfluxdeep + wplus - wminus
                deeprech(i,j) = wfluxdeep
              endif

             totwater = qlat(i,j) - qrf(i,j) + deeprech(i,j)

             smc(1:nsoil) = smois(i,1:nsoil,j)
             sh2o(1:nsoil) = sh2oxy(i,1:nsoil,j)
             smceq(1:nsoil) = smoiseq(i,1:nsoil,j)


             call updatewtd ( nsoil, dzs , zsoil, smceq, smcmax, smcwlt, psisat, bexp ,i , j , &
                              totwater, wtd(i,j), smc, sh2o, smcwtd(i,j)      , &
                              qspring(i,j) ) 


             smois(i,1:nsoil,j) = smc(1:nsoil)
             sh2oxy(i,1:nsoil,j) = sh2o(1:nsoil)

           endif
       enddo
    enddo


    do j=jts,jte
       do i=its,ite
           qslat(i,j) = qslat(i,j) + qlat(i,j)*1.e3
           qrfs(i,j) = qrfs(i,j) + qrf(i,j)*1.e3
           qsprings(i,j) = qsprings(i,j) + qspring(i,j)*1.e3
           rech(i,j) = rech(i,j) + deeprech(i,j)*1.e3

           deeprech(i,j) =0.
       enddo
    enddo


end  subroutine wtable_mmf_noahmp


subroutine lateralflow  (isltyp,wtd,qlat,fdepth,topo,landmask,deltat,area &
                           ,ids,ide,jds,jde,kds,kde                      &
                           ,ims,ime,jms,jme,kms,kme                      &
                           ,its,ite,jts,jte,kts,kte                      )

  use noahmp_tables, only : dksat_table

  implicit none

  integer,  intent(in   )   ::     ids,ide, jds,jde, kds,kde,  &
       &                           ims,ime, jms,jme, kms,kme,  &
       &                           its,ite, jts,jte, kts,kte
  real(r8)                                  , intent(in) :: deltat                                 
  integer, dimension( ims:ime, jms:jme ), intent(in) :: isltyp, landmask
  real(r8)   ,    dimension( ims:ime, jms:jme ), intent(in) :: fdepth,wtd,topo,area
  real(r8)   , dimension( ims:ime , jms:jme ), intent(out) :: qlat
  integer                              :: i, j, itsh,iteh,jtsh,jteh
  real(r8)                              :: q,klat
  real(r8)   , dimension( ims:ime , jms:jme ) :: kcell, head
  real(r8)   , dimension(19)      :: klatfactor
  data klatfactor /2.,3.,4.,10.,10.,12.,14.,20.,24.,28.,40.,48.,2.,0.,10.,0.,20.,2.,2./

  real(r8)   ,    parameter :: pi = 3.14159265 
  real(r8)   ,    parameter :: fangle = 0.45508986056   

   itsh=max(its-1,ids)
   iteh=min(ite+1,ide-1)
   jtsh=max(jts-1,jds)
   jteh=min(jte+1,jde-1)

    do j=jtsh,jteh
       do i=itsh,iteh
           if(fdepth(i,j).gt.0.)then
                 klat = dksat_table(isltyp(i,j)) * klatfactor(isltyp(i,j))
                 if(wtd(i,j) < -1.5)then
                     kcell(i,j) = fdepth(i,j) * klat * exp( (wtd(i,j) + 1.5) / fdepth(i,j) )
                 else
                     kcell(i,j) = klat * ( wtd(i,j) + 1.5 + fdepth(i,j) )  
                 endif
           else
                 kcell(i,j) = 0.
           endif

           head(i,j) = topo(i,j) + wtd(i,j)
       enddo
    enddo

   itsh=max(its,ids+1)
   iteh=min(ite,ide-2)
   jtsh=max(jts,jds+1)
   jteh=min(jte,jde-2)

    do j=jtsh,jteh
       do i=itsh,iteh
          if(landmask(i,j).gt.0)then
                 q=0.
                             
                 q  = q + (kcell(i-1,j+1)+kcell(i,j)) &
                        * (head(i-1,j+1)-head(i,j))/sqrt(2.)         
                 q  = q +  (kcell(i-1,j)+kcell(i,j)) &
                        *  (head(i-1,j)-head(i,j))
                 q  = q +  (kcell(i-1,j-1)+kcell(i,j)) &
                        * (head(i-1,j-1)-head(i,j))/sqrt(2.)
                 q  = q +  (kcell(i,j+1)+kcell(i,j)) &
                        * (head(i,j+1)-head(i,j))
                 q  = q +  (kcell(i,j-1)+kcell(i,j)) &
                        * (head(i,j-1)-head(i,j))
                 q  = q +  (kcell(i+1,j+1)+kcell(i,j)) &
                        * (head(i+1,j+1)-head(i,j))/sqrt(2.)
                 q  = q +  (kcell(i+1,j)+kcell(i,j)) &
                        * (head(i+1,j)-head(i,j))
                 q  = q +  (kcell(i+1,j-1)+kcell(i,j)) &
                        * (head(i+1,j-1)-head(i,j))/sqrt(2.)

                 qlat(i,j) = fangle* q * deltat / area(i,j)
          endif
       enddo
    enddo


end  subroutine lateralflow


  subroutine updatewtd  (nsoil,  dzs,  zsoil ,smceq                ,& 
                         smcmax, smcwlt, psisat, bexp ,iloc ,jloc  ,& 
                         totwater, wtd ,smc, sh2o ,smcwtd          ,& 
                         qspring                                 )  

  implicit none


  integer,                         intent(in) :: nsoil 
  integer,                         intent(in) :: iloc, jloc
  real(r8)   ,                         intent(in)    :: smcmax
  real(r8)   ,                         intent(in)    :: smcwlt
  real(r8)   ,                         intent(in)    :: psisat
  real(r8)   ,                         intent(in)    :: bexp
  real(r8)   ,  dimension(       0:nsoil), intent(in) :: zsoil 
  real(r8)   ,  dimension(       1:nsoil), intent(in) :: smceq  
  real(r8)   ,  dimension(       1:nsoil), intent(in) :: dzs 

  real(r8)                            , intent(inout) :: totwater
  real(r8)                            , intent(inout) :: wtd
  real(r8)                            , intent(inout) :: smcwtd
  real(r8)   , dimension(       1:nsoil), intent(inout) :: smc
  real(r8)   , dimension(       1:nsoil), intent(inout) :: sh2o

  real(r8)                            , intent(out) :: qspring

  integer                                     :: k
  integer                                     :: k1
  integer                                     :: iwtd
  integer                                     :: kwtd
  real(r8)                              :: maxwatup, maxwatdw ,wtdold
  real(r8)                              :: wgpmid
  real(r8)                              :: syielddw
  real(r8)                              :: dzup
  real(r8)                              :: smceqdeep
  real(r8)   , dimension(       1:nsoil)             :: sice


qspring=0.

sice = smc - sh2o

iwtd=1


if(totwater.gt.0.)then
         if(wtd.ge.zsoil(nsoil))then

            do k=nsoil-1,1,-1
              if(wtd.lt.zsoil(k))exit
            enddo
            iwtd=k
            kwtd=iwtd+1
            maxwatup=dzs(kwtd)*(smcmax-smc(kwtd))

            if(totwater.le.maxwatup)then
               smc(kwtd) = smc(kwtd) + totwater / dzs(kwtd)
               smc(kwtd) = min(smc(kwtd),smcmax)
               if(smc(kwtd).gt.smceq(kwtd))wtd = min ( ( smc(kwtd)*dzs(kwtd) &
                 - smceq(kwtd)*zsoil(iwtd) + smcmax*zsoil(kwtd) ) / &
                     ( smcmax-smceq(kwtd) ) , zsoil(iwtd) )
               totwater=0.
            else   
              smc(kwtd) = smcmax
              totwater=totwater-maxwatup
              k1=iwtd
              do k=k1,0,-1
                 wtd = zsoil(k)
                 iwtd=k-1
                 if(k.eq.0)exit
                 maxwatup=dzs(k)*(smcmax-smc(k))
                 if(totwater.le.maxwatup)then
                   smc(k) = smc(k) + totwater / dzs(k)
                   smc(k) = min(smc(k),smcmax)
                   if(smc(k).gt.smceq(k))wtd = min ( ( smc(k)*dzs(k) &
                     - smceq(k)*zsoil(iwtd) + smcmax*zsoil(k) ) / &
                     ( smcmax-smceq(k) ) , zsoil(iwtd) )
                   totwater=0.
                   exit
                 else
                    smc(k) = smcmax
                    totwater=totwater-maxwatup
                 endif

              enddo

            endif

         elseif(wtd.ge.zsoil(nsoil)-dzs(nsoil))then 

            
               smceqdeep = smcmax * ( psisat / &
                           (psisat - dzs(nsoil)) ) ** (1./bexp)

               smceqdeep = max(smceqdeep,1.e-4)

            maxwatup=(smcmax-smcwtd)*dzs(nsoil)

            if(totwater.le.maxwatup)then
                smcwtd = smcwtd + totwater / dzs(nsoil)
                smcwtd = min(smcwtd,smcmax)
                if(smcwtd.gt.smceqdeep)wtd = min( ( smcwtd*dzs(nsoil) &
                 - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                     ( smcmax-smceqdeep ) , zsoil(nsoil) )
                totwater=0.
            else
                smcwtd=smcmax
                totwater=totwater-maxwatup
                do k=nsoil,0,-1
                    wtd=zsoil(k)
                    iwtd=k-1
                    if(k.eq.0)exit
                    maxwatup=dzs(k)*(smcmax-smc(k))
                    if(totwater.le.maxwatup)then
                     smc(k) = min(smc(k) + totwater / dzs(k),smcmax)
                     if(smc(k).gt.smceq(k))wtd = min ( ( smc(k)*dzs(k) &
                        - smceq(k)*zsoil(iwtd) + smcmax*zsoil(k) ) / &
                           ( smcmax-smceq(k) ) , zsoil(iwtd) )
                     totwater=0.
                     exit
                    else
                     smc(k) = smcmax
                     totwater=totwater-maxwatup
                    endif
                enddo
             endif


       else

            maxwatup=(smcmax-smcwtd)*(zsoil(nsoil)-dzs(nsoil)-wtd)
            if(totwater.le.maxwatup)then
               wtd = wtd + totwater/(smcmax-smcwtd)
               totwater=0.
            else
               totwater=totwater-maxwatup
               wtd=zsoil(nsoil)-dzs(nsoil)
               maxwatup=(smcmax-smcwtd)*dzs(nsoil)
              if(totwater.le.maxwatup)then

            
               smceqdeep = smcmax * ( psisat / &
                           (psisat - dzs(nsoil)) ) ** (1./bexp)

               smceqdeep = max(smceqdeep,1.e-4)

                smcwtd = smcwtd + totwater / dzs(nsoil)
                smcwtd = min(smcwtd,smcmax)
                wtd = ( smcwtd*dzs(nsoil) &
                 - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                     ( smcmax-smceqdeep )
                totwater=0.
              else
                smcwtd=smcmax
                totwater=totwater-maxwatup
                do k=nsoil,0,-1
                    wtd=zsoil(k)
                    iwtd=k-1
                    if(k.eq.0)exit
                    maxwatup=dzs(k)*(smcmax-smc(k))

                    if(totwater.le.maxwatup)then
                     smc(k) = smc(k) + totwater / dzs(k)
                     smc(k) = min(smc(k),smcmax)
                     if(smc(k).gt.smceq(k))wtd = ( smc(k)*dzs(k) &
                        - smceq(k)*zsoil(iwtd) + smcmax*zsoil(k) ) / &
                           ( smcmax-smceq(k) )
                     totwater=0.
                     exit
                    else
                     smc(k) = smcmax
                     totwater=totwater-maxwatup
                    endif
                   enddo
               endif
             endif
         endif

        qspring=totwater

elseif(totwater.lt.0.)then


         if(wtd.ge.zsoil(nsoil))then 

            do k=nsoil-1,1,-1
               if(wtd.lt.zsoil(k))exit
            enddo
            iwtd=k

               k1=iwtd+1
               do kwtd=k1,nsoil


                  maxwatdw=dzs(kwtd)*(smc(kwtd)-max(smceq(kwtd),sice(kwtd)))

                  if(-totwater.le.maxwatdw)then
                        smc(kwtd) = smc(kwtd) + totwater / dzs(kwtd)
                        if(smc(kwtd).gt.smceq(kwtd))then
                              wtd = ( smc(kwtd)*dzs(kwtd) &
                                 - smceq(kwtd)*zsoil(iwtd) + smcmax*zsoil(kwtd) ) / &
                                 ( smcmax-smceq(kwtd) )
                         else
                              wtd=zsoil(kwtd)
                              iwtd=iwtd+1
                         endif
                         totwater=0.
                         exit
                   else
                         wtd = zsoil(kwtd)
                         iwtd=iwtd+1
                         if(maxwatdw.ge.0.)then
                            smc(kwtd) = smc(kwtd) + maxwatdw / dzs(kwtd)
                            totwater = totwater + maxwatdw
                         endif
                   endif

                enddo

               if(iwtd.eq.nsoil.and.totwater.lt.0.)then
            
               smceqdeep = smcmax * ( psisat / &
                           (psisat - dzs(nsoil)) ) ** (1./bexp)

               smceqdeep = max(smceqdeep,1.e-4)

                  maxwatdw=dzs(nsoil)*(smcwtd-smceqdeep)

                  if(-totwater.le.maxwatdw)then

                       smcwtd = smcwtd + totwater / dzs(nsoil)
                       wtd = max( ( smcwtd*dzs(nsoil) &
                           - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                            ( smcmax-smceqdeep ) , zsoil(nsoil)-dzs(nsoil) )

                  else

                       wtd=zsoil(nsoil)-dzs(nsoil)
                       smcwtd = smcwtd + totwater / dzs(nsoil)

                       dzup=(smceqdeep-smcwtd)*dzs(nsoil)/(smcmax-smceqdeep)
                       wtd=wtd-dzup
                       smcwtd=smceqdeep

                  endif

                endif



        elseif(wtd.ge.zsoil(nsoil)-dzs(nsoil))then


            
               smceqdeep = smcmax * ( psisat / &
                           (psisat - dzs(nsoil)) ) ** (1./bexp)

               smceqdeep = max(smceqdeep,1.e-4)

            maxwatdw=dzs(nsoil)*(smcwtd-smceqdeep)

            if(-totwater.le.maxwatdw)then

               smcwtd = smcwtd + totwater / dzs(nsoil)
               wtd = max( ( smcwtd*dzs(nsoil) &
                    - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                    ( smcmax-smceqdeep ) , zsoil(nsoil)-dzs(nsoil) )

            else

               wtd=zsoil(nsoil)-dzs(nsoil)
               smcwtd = smcwtd + totwater / dzs(nsoil)

               dzup=(smceqdeep-smcwtd)*dzs(nsoil)/(smcmax-smceqdeep)
               wtd=wtd-dzup
               smcwtd=smceqdeep

             endif

         else

               wgpmid = smcmax * ( psisat / &
                    (psisat - (zsoil(nsoil)-wtd)) ) ** (1./bexp)

               wgpmid=max(wgpmid,1.e-4)
               syielddw=smcmax-wgpmid
               wtdold=wtd
               wtd = wtdold + totwater/syielddw

               smcwtd = (smcwtd*(zsoil(nsoil)-wtdold)+wgpmid*(wtdold-wtd) ) / (zsoil(nsoil)-wtd)

          endif

          qspring=0.

   endif

         sh2o = smc - sice

end  subroutine updatewtd

end module module_sf_noahmp_groundwater
