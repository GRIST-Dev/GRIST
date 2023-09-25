module noahmp_tables
  use grist_constants, only: r8
  use grist_mpi
    implicit none

    integer, private, parameter :: mvt   = 27
    integer, private, parameter :: mband = 2
    integer, private, parameter :: msc   = 8
    integer, private, parameter :: max_soiltyp = 30
    integer, private, parameter :: ncrop = 5
    integer, private, parameter :: nstage = 8



    integer :: isurban_table
    integer :: iswater_table
    integer :: isbarren_table
    integer :: isice_table
    integer :: eblforest_table
    integer :: low_density_residential_table
    integer :: high_density_residential_table
    integer :: high_intensity_industrial_table

    real(r8)   :: ch2op_table(mvt)       
    real(r8)   :: dleaf_table(mvt)       
    real(r8)   :: z0mvt_table(mvt)       
    real(r8)   :: hvt_table(mvt)         
    real(r8)   :: hvb_table(mvt)         
    real(r8)   :: den_table(mvt)         
    real(r8)   :: rc_table(mvt)          
    real(r8)   :: mfsno_table(mvt)       
    real(r8)   :: saim_table(mvt,12)     
    real(r8)   :: laim_table(mvt,12)     
    real(r8)   :: sla_table(mvt)         
    real(r8)   :: dilefc_table(mvt)      
    real(r8)   :: dilefw_table(mvt)      
    real(r8)   :: fragr_table(mvt)       
    real(r8)   :: ltovrc_table(mvt)      

    real(r8)   :: c3psn_table(mvt)       
    real(r8)   :: kc25_table(mvt)        
    real(r8)   :: akc_table(mvt)         
    real(r8)   :: ko25_table(mvt)        
    real(r8)   :: ako_table(mvt)         
    real(r8)   :: vcmx25_table(mvt)      
    real(r8)   :: avcmx_table(mvt)       
    real(r8)   :: bp_table(mvt)          
    real(r8)   :: mp_table(mvt)          
    real(r8)   :: qe25_table(mvt)        
    real(r8)   :: aqe_table(mvt)         
    real(r8)   :: rmf25_table(mvt)       
    real(r8)   :: rms25_table(mvt)       
    real(r8)   :: rmr25_table(mvt)       
    real(r8)   :: arm_table(mvt)         
    real(r8)   :: folnmx_table(mvt)      
    real(r8)   :: tmin_table(mvt)        

    real(r8)   :: xl_table(mvt)          
    real(r8)   :: rhol_table(mvt,mband)  
    real(r8)   :: rhos_table(mvt,mband)  
    real(r8)   :: taul_table(mvt,mband)  
    real(r8)   :: taus_table(mvt,mband)  

    real(r8)   :: mrp_table(mvt)         
    real(r8)   :: cwpvt_table(mvt)       

    real(r8)   :: wrrat_table(mvt)       
    real(r8)   :: wdpool_table(mvt)      
    real(r8)   :: tdlef_table(mvt)       

    real(r8)   :: nroot_table(mvt)       
    real(r8)   :: rgl_table(mvt)         
    real(r8)   :: rs_table(mvt)          
    real(r8)   :: hs_table(mvt)          
    real(r8)   :: topt_table(mvt)        
    real(r8)   :: rsmax_table(mvt)       



    integer            :: slcats

    real(r8)   :: bexp_table(max_soiltyp)          
    real(r8)   :: smcdry_table(max_soiltyp)      
    real(r8)   :: f1_table(max_soiltyp)         
    real(r8)   :: smcmax_table(max_soiltyp)      
    real(r8)   :: smcref_table(max_soiltyp)      
    real(r8)   :: psisat_table(max_soiltyp)      
    real(r8)   :: dksat_table(max_soiltyp)       
    real(r8)   :: dwsat_table(max_soiltyp)       
    real(r8)   :: smcwlt_table(max_soiltyp)      
    real(r8)   :: quartz_table(max_soiltyp)         



    real(r8)   :: slope_table(9)    
    
    real(r8)   :: csoil_table       
    real(r8)   :: refdk_table       
    real(r8)   :: refkdt_table      
    real(r8)   :: frzk_table        
    real(r8)   :: zbot_table        
    real(r8)   :: czil_table        



    real(r8)   :: albsat_table(msc,mband)   
    real(r8)   :: albdry_table(msc,mband)   
    real(r8)   :: albice_table(mband)       
    real(r8)   :: alblak_table(mband)       
    real(r8)   :: omegas_table(mband)       
    real(r8)   :: betads_table              
    real(r8)   :: betais_table              
    real(r8)   :: eg_table(2)               



    real(r8)   :: co2_table      
    real(r8)   :: o2_table       
    real(r8)   :: timean_table   
    real(r8)   :: fsatmx_table   
    real(r8)   :: z0sno_table    
    real(r8)   :: ssi_table      
    real(r8)   :: swemx_table    
    real(r8)   :: rsurf_snow_table    



 integer :: pltday_table(ncrop)         
 integer :: hsday_table(ncrop)          
    real(r8)   :: plantpop_table(ncrop)       
    real(r8)   :: irri_table(ncrop)           

    real(r8)   :: gddtbase_table(ncrop)       
    real(r8)   :: gddtcut_table(ncrop)        
    real(r8)   :: gdds1_table(ncrop)          
    real(r8)   :: gdds2_table(ncrop)          
    real(r8)   :: gdds3_table(ncrop)          
    real(r8)   :: gdds4_table(ncrop)          
    real(r8)   :: gdds5_table(ncrop)          

 integer :: c3c4_table(ncrop)           
    real(r8)   :: aref_table(ncrop)           
    real(r8)   :: psnrf_table(ncrop)          
    real(r8)   :: i2par_table(ncrop)          
    real(r8)   :: tassim0_table(ncrop)        
    real(r8)   :: tassim1_table(ncrop)        
    real(r8)   :: tassim2_table(ncrop)        
    real(r8)   :: k_table(ncrop)              
    real(r8)   :: epsi_table(ncrop)           

    real(r8)   :: q10mr_table(ncrop)          
    real(r8)   :: foln_mx_table(ncrop)        
    real(r8)   :: lefreez_table(ncrop)        

    real(r8)   :: dile_fc_table(ncrop,nstage) 
    real(r8)   :: dile_fw_table(ncrop,nstage) 
    real(r8)   :: fra_gr_table(ncrop)         

    real(r8)   :: lf_ovrc_table(ncrop,nstage) 
    real(r8)   :: st_ovrc_table(ncrop,nstage) 
    real(r8)   :: rt_ovrc_table(ncrop,nstage) 
    real(r8)   :: lfmr25_table(ncrop)         
    real(r8)   :: stmr25_table(ncrop)         
    real(r8)   :: rtmr25_table(ncrop)         
    real(r8)   :: grainmr25_table(ncrop)      

    real(r8)   :: lfpt_table(ncrop,nstage)    
    real(r8)   :: stpt_table(ncrop,nstage)    
    real(r8)   :: rtpt_table(ncrop,nstage)    
    real(r8)   :: grainpt_table(ncrop,nstage) 
    real(r8)   :: bio2lai_table(ncrop)        

    integer, parameter                              :: max_cats = 100 , max_seas = 12, luseas=2, lucats=33
    real(r8), dimension(max_cats, max_seas), public :: ALBD, SLMO, SFEM, SFZ0, THERIN, SFHC
    real(r8), dimension(max_cats), public           :: SCFX, li

contains
  
  subroutine read_landuse_tbl
    integer :: ierr, iskip, ls, lc

    open(299, file='LANDUSE.TBL',form='formatted',status='old',iostat=ierr)          
        do iskip=1, 103 ! 215: modis ! 103: mminlu='modified_igbp_modis_noah'
            read (299,*)
        end do

        ls=1    !--summer
        do lc=1,lucats
            read (299,*)li(lc),albd(lc,ls),slmo(lc,ls),sfem(lc,ls),        &
                        sfz0(lc,ls),therin(lc,ls),scfx(lc),sfhc(lc,ls)
        end do
        read (299,*)

        ls=2   !--winter
        do lc=1,lucats
            read (299,*)li(lc),albd(lc,ls),slmo(lc,ls),sfem(lc,ls),        &
                        sfz0(lc,ls),therin(lc,ls),scfx(lc),sfhc(lc,ls)
        end do
    close(299)
  end subroutine read_landuse_tbl

  subroutine read_mp_veg_parameters(dataset_identifier)
    implicit none
    character(len=*), intent(in) :: dataset_identifier
    integer :: ierr
    integer :: ik,im

    integer :: nveg
    character(len=256) :: veg_dataset_description

    integer :: isurban
    integer :: iswater
    integer :: isbarren
    integer :: isice
    integer :: eblforest
    integer :: low_density_residential
    integer :: high_density_residential
    integer :: high_intensity_industrial

    real(r8)  , dimension(mvt) :: sai_jan,sai_feb,sai_mar,sai_apr,sai_may,sai_jun, &
                                     sai_jul,sai_aug,sai_sep,sai_oct,sai_nov,sai_dec
    real(r8)  , dimension(mvt) :: lai_jan,lai_feb,lai_mar,lai_apr,lai_may,lai_jun, &
                                     lai_jul,lai_aug,lai_sep,lai_oct,lai_nov,lai_dec
    real(r8)  , dimension(mvt) :: rhol_vis, rhol_nir, rhos_vis, rhos_nir, &
                                     taul_vis, taul_nir, taus_vis, taus_nir
    real(r8)  , dimension(mvt) :: ch2op, dleaf, z0mvt, hvt, hvb, den, rc, mfsno, xl, cwpvt, c3psn, kc25, akc, ko25, ako, &
                     avcmx, aqe, ltovrc,  dilefc,  dilefw,  rmf25 ,  sla   ,  fragr ,  tmin  ,  vcmx25,  tdlef ,  &
                     bp, mp, qe25, rms25, rmr25, arm, folnmx, wdpool, wrrat, mrp, nroot, rgl, rs, hs, topt, rsmax, &
		     slarea, eps1, eps2, eps3, eps4, eps5
			     
    namelist / noahmp_usgs_veg_categories / veg_dataset_description, nveg
    namelist / noahmp_usgs_parameters / isurban, iswater, isbarren, isice, eblforest, &
         low_density_residential, high_density_residential, high_intensity_industrial, &
         ch2op, dleaf, z0mvt, hvt, hvb, den, rc, mfsno, xl, cwpvt, c3psn, kc25, akc, ko25, ako, avcmx, aqe, &
         ltovrc,  dilefc,  dilefw,  rmf25 ,  sla   ,  fragr ,  tmin  ,  vcmx25,  tdlef ,  bp, mp, qe25, rms25, rmr25, arm, &
         folnmx, wdpool, wrrat, mrp, nroot, rgl, rs, hs, topt, rsmax, &
         sai_jan, sai_feb, sai_mar, sai_apr, sai_may, sai_jun,sai_jul,sai_aug,sai_sep,sai_oct,sai_nov,sai_dec, &
         lai_jan, lai_feb, lai_mar, lai_apr, lai_may, lai_jun,lai_jul,lai_aug,lai_sep,lai_oct,lai_nov,lai_dec, &
         rhol_vis, rhol_nir, rhos_vis, rhos_nir, taul_vis, taul_nir, taus_vis, taus_nir, slarea, eps1, eps2, eps3, eps4, eps5
	 
    namelist / noahmp_modis_veg_categories / veg_dataset_description, nveg
    namelist / noahmp_modis_parameters / isurban, iswater, isbarren, isice, eblforest, &
         low_density_residential, high_density_residential, high_intensity_industrial, &
         ch2op, dleaf, z0mvt, hvt, hvb, den, rc, mfsno, xl, cwpvt, c3psn, kc25, akc, ko25, ako, avcmx, aqe, &
         ltovrc,  dilefc,  dilefw,  rmf25 ,  sla   ,  fragr ,  tmin  ,  vcmx25,  tdlef ,  bp, mp, qe25, rms25, rmr25, arm, &
         folnmx, wdpool, wrrat, mrp, nroot, rgl, rs, hs, topt, rsmax, &
         sai_jan, sai_feb, sai_mar, sai_apr, sai_may, sai_jun,sai_jul,sai_aug,sai_sep,sai_oct,sai_nov,sai_dec, &
         lai_jan, lai_feb, lai_mar, lai_apr, lai_may, lai_jun,lai_jul,lai_aug,lai_sep,lai_oct,lai_nov,lai_dec, &
         rhol_vis, rhol_nir, rhos_vis, rhos_nir, taul_vis, taul_nir, taus_vis, taus_nir, slarea, eps1, eps2, eps3, eps4, eps5

    
    ch2op_table  = -1.e36
    dleaf_table  = -1.e36
    z0mvt_table  = -1.e36
    hvt_table    = -1.e36
    hvb_table    = -1.e36
    den_table    = -1.e36
    rc_table     = -1.e36
    mfsno_table  = -1.e36
    rhol_table   = -1.e36
    rhos_table   = -1.e36
    taul_table   = -1.e36
    taus_table   = -1.e36
    xl_table     = -1.e36
    cwpvt_table  = -1.e36
    c3psn_table  = -1.e36
    kc25_table   = -1.e36
    akc_table    = -1.e36
    ko25_table   = -1.e36
    ako_table    = -1.e36
    avcmx_table  = -1.e36
    aqe_table    = -1.e36
    ltovrc_table = -1.e36
    dilefc_table = -1.e36
    dilefw_table = -1.e36
    rmf25_table  = -1.e36
    sla_table    = -1.e36
    fragr_table  = -1.e36
    tmin_table   = -1.e36
    vcmx25_table = -1.e36
    tdlef_table  = -1.e36
    bp_table     = -1.e36
    mp_table     = -1.e36
    qe25_table   = -1.e36
    rms25_table  = -1.e36
    rmr25_table  = -1.e36
    arm_table    = -1.e36
    folnmx_table = -1.e36
    wdpool_table = -1.e36
    wrrat_table  = -1.e36
    mrp_table    = -1.e36
    saim_table   = -1.e36
    laim_table   = -1.e36
    nroot_table  = -1.e36
    rgl_table    = -1.e36
    rs_table     = -1.e36
    hs_table     = -1.e36
    topt_table   = -1.e36
    rsmax_table  = -1.e36
    isurban_table      = -99999
    iswater_table      = -99999
    isbarren_table     = -99999
    isice_table        = -99999
    eblforest_table    = -99999
    low_density_residential_table   = -99999
    high_density_residential_table  = -99999
    high_intensity_industrial_table = -99999

    
    !open(15, file="/home/cheyz/WORK/CAMS/grist-2019-08-19/local_grist3_20180819/par_grist/src/atmosphere/physics/surface_phys/MPTABLE.TBL", status='old', form='formatted', action='read', iostat=ierr)
    open(15, file="./MPTABLE.TBL", status='old', form='formatted', action='read', iostat=ierr)

    if (ierr /= 0) then
       write(*,'("****** error ******************************************************")')
       write(*,'("cannot find file MPTABLE.TBL")')
       write(*,'("stop")')
       write(*,'("*******************************************************************")')
       print*, (&
"stop in noah-mp read_mp_veg_parameters")
    endif
    
    if ( trim(dataset_identifier) == "usgs" ) then
       read(15,noahmp_usgs_veg_categories)
       read(15,noahmp_usgs_parameters)
    else if ( trim(dataset_identifier) == "modified_igbp_modis_noah" ) then
       read(15,noahmp_modis_veg_categories)
       read(15,noahmp_modis_parameters)
    else
       write(*,'("unrecognized dataset_identifier in subroutine read_mp_veg_parameters")')
       write(*,'("dataset_identifier = ''", a, "''")') trim(dataset_identifier)
       print*, (&
"stop in noah-mp read_mp_veg_parameters")
    endif
    close(15)

                      isurban_table   = isurban
                      iswater_table   = iswater
                     isbarren_table   = isbarren
                        isice_table   = isice
                    eblforest_table   = eblforest
      low_density_residential_table   = low_density_residential
     high_density_residential_table   = high_density_residential
    high_intensity_industrial_table   = high_intensity_industrial

     ch2op_table(1:nveg)  = ch2op(1:nveg)
     dleaf_table(1:nveg)  = dleaf(1:nveg)
     z0mvt_table(1:nveg)  = z0mvt(1:nveg)
       hvt_table(1:nveg)  = hvt(1:nveg)
       hvb_table(1:nveg)  = hvb(1:nveg)
       den_table(1:nveg)  = den(1:nveg)
        rc_table(1:nveg)  = rc(1:nveg)
     mfsno_table(1:nveg)  = mfsno(1:nveg)
        xl_table(1:nveg)  = xl(1:nveg)
     cwpvt_table(1:nveg)  = cwpvt(1:nveg)
     c3psn_table(1:nveg)  = c3psn(1:nveg)
      kc25_table(1:nveg)  = kc25(1:nveg)
       akc_table(1:nveg)  = akc(1:nveg)
      ko25_table(1:nveg)  = ko25(1:nveg)
       ako_table(1:nveg)  = ako(1:nveg)
     avcmx_table(1:nveg)  = avcmx(1:nveg)
       aqe_table(1:nveg)  = aqe(1:nveg)
    ltovrc_table(1:nveg)  = ltovrc(1:nveg)
    dilefc_table(1:nveg)  = dilefc(1:nveg)
    dilefw_table(1:nveg)  = dilefw(1:nveg)
     rmf25_table(1:nveg)  = rmf25(1:nveg)
       sla_table(1:nveg)  = sla(1:nveg)
     fragr_table(1:nveg)  = fragr(1:nveg)
      tmin_table(1:nveg)  = tmin(1:nveg)
    vcmx25_table(1:nveg)  = vcmx25(1:nveg)
     tdlef_table(1:nveg)  = tdlef(1:nveg)
        bp_table(1:nveg)  = bp(1:nveg)
        mp_table(1:nveg)  = mp(1:nveg)
      qe25_table(1:nveg)  = qe25(1:nveg)
     rms25_table(1:nveg)  = rms25(1:nveg)
     rmr25_table(1:nveg)  = rmr25(1:nveg)
       arm_table(1:nveg)  = arm(1:nveg)
    folnmx_table(1:nveg)  = folnmx(1:nveg)
    wdpool_table(1:nveg)  = wdpool(1:nveg)
     wrrat_table(1:nveg)  = wrrat(1:nveg)
       mrp_table(1:nveg)  = mrp(1:nveg)
     nroot_table(1:nveg)  = nroot(1:nveg)
       rgl_table(1:nveg)  = rgl(1:nveg)
        rs_table(1:nveg)  = rs(1:nveg)
        hs_table(1:nveg)  = hs(1:nveg)
      topt_table(1:nveg)  = topt(1:nveg)
     rsmax_table(1:nveg)  = rsmax(1:nveg)
    
    

    saim_table(1:nveg, 1) = sai_jan(1:nveg)
    saim_table(1:nveg, 2) = sai_feb(1:nveg)
    saim_table(1:nveg, 3) = sai_mar(1:nveg)
    saim_table(1:nveg, 4) = sai_apr(1:nveg)
    saim_table(1:nveg, 5) = sai_may(1:nveg)
    saim_table(1:nveg, 6) = sai_jun(1:nveg)
    saim_table(1:nveg, 7) = sai_jul(1:nveg)
    saim_table(1:nveg, 8) = sai_aug(1:nveg)
    saim_table(1:nveg, 9) = sai_sep(1:nveg)
    saim_table(1:nveg,10) = sai_oct(1:nveg)
    saim_table(1:nveg,11) = sai_nov(1:nveg)
    saim_table(1:nveg,12) = sai_dec(1:nveg)

    laim_table(1:nveg, 1) = lai_jan(1:nveg)
    laim_table(1:nveg, 2) = lai_feb(1:nveg)
    laim_table(1:nveg, 3) = lai_mar(1:nveg)
    laim_table(1:nveg, 4) = lai_apr(1:nveg)
    laim_table(1:nveg, 5) = lai_may(1:nveg)
    laim_table(1:nveg, 6) = lai_jun(1:nveg)
    laim_table(1:nveg, 7) = lai_jul(1:nveg)
    laim_table(1:nveg, 8) = lai_aug(1:nveg)
    laim_table(1:nveg, 9) = lai_sep(1:nveg)
    laim_table(1:nveg,10) = lai_oct(1:nveg)
    laim_table(1:nveg,11) = lai_nov(1:nveg)
    laim_table(1:nveg,12) = lai_dec(1:nveg)

    rhol_table(1:nveg,1)  = rhol_vis(1:nveg) 
    rhol_table(1:nveg,2)  = rhol_nir(1:nveg) 
    rhos_table(1:nveg,1)  = rhos_vis(1:nveg) 
    rhos_table(1:nveg,2)  = rhos_nir(1:nveg) 
    taul_table(1:nveg,1)  = taul_vis(1:nveg) 
    taul_table(1:nveg,2)  = taul_nir(1:nveg) 
    taus_table(1:nveg,1)  = taus_vis(1:nveg) 
    taus_table(1:nveg,2)  = taus_nir(1:nveg) 

  end subroutine read_mp_veg_parameters

  subroutine read_mp_soil_parameters()
    implicit none
    integer :: ierr
    character*4         :: sltype
    integer             :: itmp, num_slope, lc
    character(len=256)  :: message
    

    
       bexp_table = -1.e36
     smcdry_table = -1.e36
         f1_table = -1.e36
     smcmax_table = -1.e36
     smcref_table = -1.e36
     psisat_table = -1.e36
      dksat_table = -1.e36
      dwsat_table = -1.e36
     smcwlt_table = -1.e36
     quartz_table = -1.e36
      slope_table = -1.e36
      csoil_table = -1.e36
      refdk_table = -1.e36
     refkdt_table = -1.e36
       frzk_table = -1.e36
       zbot_table = -1.e36
       czil_table = -1.e36




    !open(19, file='/home/cheyz/WORK/CAMS/grist-2019-08-19/local_grist3_20180819/par_grist/src/atmosphere/physics/surface_phys/SOILPARM.TBL',form='formatted',status='old',iostat=ierr)
    open(19, file='./SOILPARM.TBL',form='formatted',status='old',iostat=ierr)

    if(ierr .ne. 0 ) then
    if(mpi_rank()==0)then
      write(message,fmt='(a)') 'module_sf_noahmpdrv.f: read_mp_soil_parameters: failure opening SOILPARM.TBL'
    !  print*, (&
    !           message )
    end if
    end if

    read (19,*)
    read (19,*) sltype
    read (19,*) slcats
    if(mpi_rank()==0)then
    write( message , * ) 'soil texture classification = ', trim ( sltype ) , ' found', &
               slcats,' categories'
    !print*,  ( message )
    end if

    do lc=1,slcats
      read (19,*) itmp,bexp_table(lc),smcdry_table(lc),f1_table(lc),smcmax_table(lc),   &
                  smcref_table(lc),psisat_table(lc),dksat_table(lc), dwsat_table(lc),   &
                  smcwlt_table(lc), quartz_table(lc)
    enddo

    close (19)




    !open(19, file='/home/cheyz/WORK/CAMS/grist-2019-08-19/local_grist3_20180819/par_grist/src/atmosphere/physics/surface_phys/GENPARM.TBL',form='formatted',status='old',iostat=ierr)
    open(19, file='./GENPARM.TBL',form='formatted',status='old',iostat=ierr)
    if(ierr .ne. 0 ) then
    if(mpi_rank()==0)then
      write(message,fmt='(a)') 'module_sf_noahlsm.f: read_mp_soil_parameters: failure opening genparm.tbl'
!      print*, (&
!message )
    end if
    end if

    read (19,*)
    read (19,*)
    read (19,*) num_slope

    do lc=1,num_slope
        read (19,*) slope_table(lc)
    enddo

    read (19,*)
    read (19,*)
    read (19,*)
    read (19,*)
    read (19,*)
    read (19,*) csoil_table
    read (19,*)
    read (19,*)
    read (19,*)
    read (19,*) refdk_table
    read (19,*)
    read (19,*) refkdt_table
    read (19,*)
    read (19,*) frzk_table
    read (19,*)
    read (19,*) zbot_table
    read (19,*)
    read (19,*) czil_table
    read (19,*)
    read (19,*)
    read (19,*)
    read (19,*)

    close (19)

  end subroutine read_mp_soil_parameters

  subroutine read_mp_rad_parameters()
    implicit none
    integer :: ierr

    real(r8)   :: albice(mband),alblak(mband),omegas(mband),betads,betais,eg(2)
    real(r8)   :: albsat_vis(msc)
    real(r8)   :: albsat_nir(msc)
    real(r8)   :: albdry_vis(msc)
    real(r8)   :: albdry_nir(msc)

    namelist / noahmp_rad_parameters / albsat_vis,albsat_nir,albdry_vis,albdry_nir,albice,alblak,omegas,betads,betais,eg


    
    albsat_table     = -1.e36
    albdry_table     = -1.e36
    albice_table     = -1.e36
    alblak_table     = -1.e36
    omegas_table     = -1.e36
    betads_table     = -1.e36
    betais_table     = -1.e36
    eg_table         = -1.e36

    !open(15, file="/home/cheyz/WORK/CAMS/grist-2019-08-19/local_grist3_20180819/par_grist/src/atmosphere/physics/surface_phys/MPTABLE.TBL", status='old', form='formatted', action='read', iostat=ierr)
    open(15, file="./MPTABLE.TBL", status='old', form='formatted', action='read', iostat=ierr)

    if (ierr /= 0) then
       write(*,'("****** error ******************************************************")')
       write(*,'("cannot find file mptable.tbl")')
       write(*,'("stop")')
       write(*,'("*******************************************************************")')
       print*, (&
"stop in noah-mp read_mp_rad_parameters")
    endif

    read(15,noahmp_rad_parameters)
    close(15)

    albsat_table(:,1) = albsat_vis 
    albsat_table(:,2) = albsat_nir 
    albdry_table(:,1) = albdry_vis 
    albdry_table(:,2) = albdry_nir 
    albice_table      = albice
    alblak_table      = alblak
    omegas_table      = omegas
    betads_table      = betads
    betais_table      = betais
    eg_table          = eg

  end subroutine read_mp_rad_parameters

  subroutine read_mp_global_parameters()
    implicit none
    integer :: ierr

    real(r8)   :: co2,o2,timean,fsatmx,z0sno,ssi,swemx,rsurf_snow

    namelist / noahmp_global_parameters / co2,o2,timean,fsatmx,z0sno,ssi,swemx,rsurf_snow


    
       co2_table     = -1.e36
        o2_table     = -1.e36
    timean_table     = -1.e36
    fsatmx_table     = -1.e36
     z0sno_table     = -1.e36
       ssi_table     = -1.e36
     swemx_table     = -1.e36
rsurf_snow_table     = -1.e36

    open(15, file="./MPTABLE.TBL", status='old', form='formatted', action='read', iostat=ierr)
    if (ierr /= 0) then
       write(*,'("****** error ******************************************************")')
       write(*,'("cannot find file mptable.tbl")')
       write(*,'("stop")')
       write(*,'("*******************************************************************")')
       print*, (&
               "stop in noah-mp read_mp_global_parameters")
    endif

    read(15,noahmp_global_parameters)
    close(15)

       co2_table     = co2
        o2_table     = o2
    timean_table     = timean
    fsatmx_table     = fsatmx
     z0sno_table     = z0sno
       ssi_table     = ssi
     swemx_table     = swemx
rsurf_snow_table     = rsurf_snow

  end subroutine read_mp_global_parameters

  subroutine read_mp_crop_parameters()
    implicit none
    integer :: ierr

 integer, dimension(ncrop) :: pltday
 integer, dimension(ncrop) :: hsday
    real(r8)  , dimension(ncrop) :: plantpop
    real(r8)  , dimension(ncrop) :: irri
    real(r8)  , dimension(ncrop) :: gddtbase
    real(r8)  , dimension(ncrop) :: gddtcut
    real(r8)  , dimension(ncrop) :: gdds1
    real(r8)  , dimension(ncrop) :: gdds2
    real(r8)  , dimension(ncrop) :: gdds3
    real(r8)  , dimension(ncrop) :: gdds4
    real(r8)  , dimension(ncrop) :: gdds5
 integer, dimension(ncrop) :: c3c4
    real(r8)  , dimension(ncrop) :: aref
    real(r8)  , dimension(ncrop) :: psnrf
    real(r8)  , dimension(ncrop) :: i2par
    real(r8)  , dimension(ncrop) :: tassim0
    real(r8)  , dimension(ncrop) :: tassim1
    real(r8)  , dimension(ncrop) :: tassim2
    real(r8)  , dimension(ncrop) :: k
    real(r8)  , dimension(ncrop) :: epsi
    real(r8)  , dimension(ncrop) :: q10mr
    real(r8)  , dimension(ncrop) :: foln_mx
    real(r8)  , dimension(ncrop) :: lefreez
    real(r8)  , dimension(ncrop) :: dile_fc_s1,dile_fc_s2,dile_fc_s3,dile_fc_s4,dile_fc_s5,dile_fc_s6,dile_fc_s7,dile_fc_s8
    real(r8)  , dimension(ncrop) :: dile_fw_s1,dile_fw_s2,dile_fw_s3,dile_fw_s4,dile_fw_s5,dile_fw_s6,dile_fw_s7,dile_fw_s8
    real(r8)  , dimension(ncrop) :: fra_gr
    real(r8)  , dimension(ncrop) :: lf_ovrc_s1,lf_ovrc_s2,lf_ovrc_s3,lf_ovrc_s4,lf_ovrc_s5,lf_ovrc_s6,lf_ovrc_s7,lf_ovrc_s8
    real(r8)  , dimension(ncrop) :: st_ovrc_s1,st_ovrc_s2,st_ovrc_s3,st_ovrc_s4,st_ovrc_s5,st_ovrc_s6,st_ovrc_s7,st_ovrc_s8
    real(r8)  , dimension(ncrop) :: rt_ovrc_s1,rt_ovrc_s2,rt_ovrc_s3,rt_ovrc_s4,rt_ovrc_s5,rt_ovrc_s6,rt_ovrc_s7,rt_ovrc_s8
    real(r8)  , dimension(ncrop) :: lfmr25
    real(r8)  , dimension(ncrop) :: stmr25
    real(r8)  , dimension(ncrop) :: rtmr25
    real(r8)  , dimension(ncrop) :: grainmr25
    real(r8)  , dimension(ncrop) :: lfpt_s1,lfpt_s2,lfpt_s3,lfpt_s4,lfpt_s5,lfpt_s6,lfpt_s7,lfpt_s8
    real(r8)  , dimension(ncrop) :: stpt_s1,stpt_s2,stpt_s3,stpt_s4,stpt_s5,stpt_s6,stpt_s7,stpt_s8
    real(r8)  , dimension(ncrop) :: rtpt_s1,rtpt_s2,rtpt_s3,rtpt_s4,rtpt_s5,rtpt_s6,rtpt_s7,rtpt_s8
    real(r8)  , dimension(ncrop) :: grainpt_s1,grainpt_s2,grainpt_s3,grainpt_s4,grainpt_s5,grainpt_s6,grainpt_s7,grainpt_s8
    real(r8)  , dimension(ncrop) :: bio2lai


    namelist / noahmp_crop_parameters /     pltday,     hsday,  plantpop,      irri,  gddtbase,   gddtcut,     gdds1,     gdds2, &
                                             gdds3,     gdds4,     gdds5,      c3c4,      aref,     psnrf,     i2par,   tassim0, &
                                           tassim1,   tassim2,         k,      epsi,     q10mr,   foln_mx,   lefreez,            &
                                        dile_fc_s1,dile_fc_s2,dile_fc_s3,dile_fc_s4,dile_fc_s5,dile_fc_s6,dile_fc_s7,dile_fc_s8, &
                                        dile_fw_s1,dile_fw_s2,dile_fw_s3,dile_fw_s4,dile_fw_s5,dile_fw_s6,dile_fw_s7,dile_fw_s8, &
                                            fra_gr,                                                                              &
                                        lf_ovrc_s1,lf_ovrc_s2,lf_ovrc_s3,lf_ovrc_s4,lf_ovrc_s5,lf_ovrc_s6,lf_ovrc_s7,lf_ovrc_s8, &
                                        st_ovrc_s1,st_ovrc_s2,st_ovrc_s3,st_ovrc_s4,st_ovrc_s5,st_ovrc_s6,st_ovrc_s7,st_ovrc_s8, &
                                        rt_ovrc_s1,rt_ovrc_s2,rt_ovrc_s3,rt_ovrc_s4,rt_ovrc_s5,rt_ovrc_s6,rt_ovrc_s7,rt_ovrc_s8, &
                                            lfmr25,    stmr25,    rtmr25, grainmr25,                                             &
                                           lfpt_s1,   lfpt_s2,   lfpt_s3,   lfpt_s4,   lfpt_s5,   lfpt_s6,   lfpt_s7,   lfpt_s8, &
                                           stpt_s1,   stpt_s2,   stpt_s3,   stpt_s4,   stpt_s5,   stpt_s6,   stpt_s7,   stpt_s8, &
                                           rtpt_s1,   rtpt_s2,   rtpt_s3,   rtpt_s4,   rtpt_s5,   rtpt_s6,   rtpt_s7,   rtpt_s8, &
                                        grainpt_s1,grainpt_s2,grainpt_s3,grainpt_s4,grainpt_s5,grainpt_s6,grainpt_s7,grainpt_s8, &
                                           bio2lai


    
       pltday_table     = -99999
        hsday_table     = -99999
     plantpop_table     = -1.e36
         irri_table     = -1.e36
     gddtbase_table     = -1.e36
      gddtcut_table     = -1.e36
        gdds1_table     = -1.e36
        gdds2_table     = -1.e36
        gdds3_table     = -1.e36
        gdds4_table     = -1.e36
        gdds5_table     = -1.e36
         c3c4_table     = -99999
         aref_table     = -1.e36
        psnrf_table     = -1.e36
        i2par_table     = -1.e36
      tassim0_table     = -1.e36
      tassim1_table     = -1.e36
      tassim2_table     = -1.e36
            k_table     = -1.e36
         epsi_table     = -1.e36
        q10mr_table     = -1.e36
      foln_mx_table     = -1.e36
      lefreez_table     = -1.e36
      dile_fc_table     = -1.e36
      dile_fw_table     = -1.e36
       fra_gr_table     = -1.e36
      lf_ovrc_table     = -1.e36
      st_ovrc_table     = -1.e36
      rt_ovrc_table     = -1.e36
       lfmr25_table     = -1.e36
       stmr25_table     = -1.e36
       rtmr25_table     = -1.e36
    grainmr25_table     = -1.e36
         lfpt_table     = -1.e36
         stpt_table     = -1.e36
         rtpt_table     = -1.e36
      grainpt_table     = -1.e36
      bio2lai_table     = -1.e36


    open(15, file="./MPTABLE.TBL", status='old', form='formatted', action='read', iostat=ierr)
    if (ierr /= 0) then
       write(*,'("****** error ******************************************************")')
       write(*,'("cannot find file mptable.tbl")')
       write(*,'("stop")')
       write(*,'("*******************************************************************")')
       print*, (&
"stop in noah-mp read_mp_crop_parameters")
    endif

    read(15,noahmp_crop_parameters)
    close(15)

       pltday_table      = pltday
        hsday_table      = hsday
     plantpop_table      = plantpop
         irri_table      = irri
     gddtbase_table      = gddtbase
      gddtcut_table      = gddtcut
        gdds1_table      = gdds1
        gdds2_table      = gdds2
        gdds3_table      = gdds3
        gdds4_table      = gdds4
        gdds5_table      = gdds5
         c3c4_table      = c3c4
         aref_table      = aref
        psnrf_table      = psnrf
        i2par_table      = i2par
      tassim0_table      = tassim0
      tassim1_table      = tassim1
      tassim2_table      = tassim2
            k_table      = k
         epsi_table      = epsi
        q10mr_table      = q10mr
      foln_mx_table      = foln_mx
      lefreez_table      = lefreez
      dile_fc_table(:,1) = dile_fc_s1
      dile_fc_table(:,2) = dile_fc_s2
      dile_fc_table(:,3) = dile_fc_s3
      dile_fc_table(:,4) = dile_fc_s4
      dile_fc_table(:,5) = dile_fc_s5
      dile_fc_table(:,6) = dile_fc_s6
      dile_fc_table(:,7) = dile_fc_s7
      dile_fc_table(:,8) = dile_fc_s8
      dile_fw_table(:,1) = dile_fw_s1
      dile_fw_table(:,2) = dile_fw_s2
      dile_fw_table(:,3) = dile_fw_s3
      dile_fw_table(:,4) = dile_fw_s4
      dile_fw_table(:,5) = dile_fw_s5
      dile_fw_table(:,6) = dile_fw_s6
      dile_fw_table(:,7) = dile_fw_s7
      dile_fw_table(:,8) = dile_fw_s8
       fra_gr_table      = fra_gr
      lf_ovrc_table(:,1) = lf_ovrc_s1
      lf_ovrc_table(:,2) = lf_ovrc_s2
      lf_ovrc_table(:,3) = lf_ovrc_s3
      lf_ovrc_table(:,4) = lf_ovrc_s4
      lf_ovrc_table(:,5) = lf_ovrc_s5
      lf_ovrc_table(:,6) = lf_ovrc_s6
      lf_ovrc_table(:,7) = lf_ovrc_s7
      lf_ovrc_table(:,8) = lf_ovrc_s8
      st_ovrc_table(:,1) = st_ovrc_s1
      st_ovrc_table(:,2) = st_ovrc_s2
      st_ovrc_table(:,3) = st_ovrc_s3
      st_ovrc_table(:,4) = st_ovrc_s4
      st_ovrc_table(:,5) = st_ovrc_s5
      st_ovrc_table(:,6) = st_ovrc_s6
      st_ovrc_table(:,7) = st_ovrc_s7
      st_ovrc_table(:,8) = st_ovrc_s8
      rt_ovrc_table(:,1) = rt_ovrc_s1
      rt_ovrc_table(:,2) = rt_ovrc_s2
      rt_ovrc_table(:,3) = rt_ovrc_s3
      rt_ovrc_table(:,4) = rt_ovrc_s4
      rt_ovrc_table(:,5) = rt_ovrc_s5
      rt_ovrc_table(:,6) = rt_ovrc_s6
      rt_ovrc_table(:,7) = rt_ovrc_s7
      rt_ovrc_table(:,8) = rt_ovrc_s8
       lfmr25_table      = lfmr25
       stmr25_table      = stmr25
       rtmr25_table      = rtmr25
    grainmr25_table      = grainmr25
         lfpt_table(:,1) = lfpt_s1
         lfpt_table(:,2) = lfpt_s2
         lfpt_table(:,3) = lfpt_s3
         lfpt_table(:,4) = lfpt_s4
         lfpt_table(:,5) = lfpt_s5
         lfpt_table(:,6) = lfpt_s6
         lfpt_table(:,7) = lfpt_s7
         lfpt_table(:,8) = lfpt_s8
         stpt_table(:,1) = stpt_s1
         stpt_table(:,2) = stpt_s2
         stpt_table(:,3) = stpt_s3
         stpt_table(:,4) = stpt_s4
         stpt_table(:,5) = stpt_s5
         stpt_table(:,6) = stpt_s6
         stpt_table(:,7) = stpt_s7
         stpt_table(:,8) = stpt_s8
         rtpt_table(:,1) = rtpt_s1
         rtpt_table(:,2) = rtpt_s2
         rtpt_table(:,3) = rtpt_s3
         rtpt_table(:,4) = rtpt_s4
         rtpt_table(:,5) = rtpt_s5
         rtpt_table(:,6) = rtpt_s6
         rtpt_table(:,7) = rtpt_s7
         rtpt_table(:,8) = rtpt_s8
      grainpt_table(:,1) = grainpt_s1
      grainpt_table(:,2) = grainpt_s2
      grainpt_table(:,3) = grainpt_s3
      grainpt_table(:,4) = grainpt_s4
      grainpt_table(:,5) = grainpt_s5
      grainpt_table(:,6) = grainpt_s6
      grainpt_table(:,7) = grainpt_s7
      grainpt_table(:,8) = grainpt_s8
      bio2lai_table      = bio2lai

  end subroutine read_mp_crop_parameters

end module noahmp_tables
