!=============================================================
!  Created by LiXiaohan on 19/5/13.
!  Compute cosine of solar zenith angle
!  and earths orbital parameters using Dave Threshers formula
!=============================================================

 module grist_zenith
    use grist_constants,        only: i4, r8, pi, one, zero
    use grist_handle_error,     only: endrun
    use grist_mpi


    implicit none
    private

    public :: zenith,       &
              orb_decl,     &  
              orb_params  

    real(r8) :: eccen       ! Earth's eccentricity factor (unitless) (typically 0 to 0.1)
    real(r8) :: mvelp       ! moving vernal equinox long
    real(r8) :: mvelpp      ! Earth's moving vernal equinox longitude
                            ! of perihelion plus pi (radians)
    real(r8) :: obliqr      ! Earth's obliquity in radians
    real(r8) :: obliq       ! obliquity in degrees
    real(r8) :: lambm0      ! Mean longitude of perihelion at the vernal equinox (radians)

    real   (r8),public,parameter :: SHR_ORB_UNDEF_REAL = 1.e36_r8          ! undefined real 
    integer(i4),public,parameter :: SHR_ORB_UNDEF_INT  = 2000000000        ! undefined int

    ! PRIVATE: by default everything else is private to this module
    real   (r8),parameter :: SHR_ORB_ECCEN_MIN  =   0.0_r8 ! min value for eccen
    real   (r8),parameter :: SHR_ORB_ECCEN_MAX  =   0.1_r8 ! max value for eccen
    real   (r8),parameter :: SHR_ORB_OBLIQ_MIN  = -90.0_r8 ! min value for obliq
    real   (r8),parameter :: SHR_ORB_OBLIQ_MAX  = +90.0_r8 ! max value for obliq
    real   (r8),parameter :: SHR_ORB_MVELP_MIN  =   0.0_r8 ! min value for mvelp
    real   (r8),parameter :: SHR_ORB_MVELP_MAX  = 360.0_r8 ! max value for mvelp


 contains
   subroutine zenith(calday, coszrs, ncol, lon, lat,doAquaPlanet)

   ! io
   integer, intent(in)   :: ncol                ! number of positions
   real(r8), intent(in)  :: calday              ! Calendar day, including fraction
   real(r8), intent(in)  :: lon(ncol)           ! longitude (radians)
   real(r8), intent(in)  :: lat(ncol)           ! latitude (radians)
   real(r8), intent(out) :: coszrs(ncol)        ! Cosine solar zenith angle
   logical,intent(in), optional :: doAquaPlanet
   ! local
   integer  :: i        ! Position loop index
   real(r8) :: delta    ! Solar declination angle in radians
   real(r8) :: eccf     ! Earth orbit eccentricity factor

   if(present(doAquaPlanet))then
      call orb_decl(calday  ,delta   ,eccf, doAquaPlanet )
   else
      call orb_decl(calday  ,delta   ,eccf)
   end if

   do i=1,ncol
      coszrs(i) = shr_orb_cosz( calday, lat(i), lon(i), delta)
   end do

   end subroutine zenith


   ! FUNCTION to return the cosine of the solar zenith angle.
   real(r8) function shr_orb_cosz(jday,lat,lon,declin)

   real   (r8),intent(in) :: jday   ! Julian cal day (1.xx to 365.xx)
   real   (r8),intent(in) :: lat    ! Centered latitude (radians)
   real   (r8),intent(in) :: lon    ! Centered longitude (radians)
   real   (r8),intent(in) :: declin ! Solar declination (radians)

   shr_orb_cosz = sin(lat)*sin(declin) - &
                  cos(lat)*cos(declin)*cos(jday*2.0_r8*pi + lon)

   end function shr_orb_cosz


   subroutine orb_params( iyear, log_print )

   integer(i4),intent(in)    :: iyear     ! Year to calculate orbit for
   logical    ,intent(in)    :: log_print ! Flags print of status/error

   !------------------------------ Parameters ----------------------------------
   integer(i4),parameter :: poblen =47 ! # of elements in series wrt obliquity
   integer(i4),parameter :: pecclen=19 ! # of elements in series wrt eccentricity
   integer(i4),parameter :: pmvelen=78 ! # of elements in series wrt vernal equinox
   real   (r8),parameter :: psecdeg = 1.0_r8/3600.0_r8 ! arc sec to deg conversion

   real   (r8) :: degrad = pi/180._r8   ! degree to radian conversion factor
   real   (r8) :: yb4_1950AD         ! number of years before 1950 AD

   character(len=*),parameter :: subname = '(orb_params)'
 
   ! Cosine series data for computation of obliquity: amplitude (arc seconds),
   ! rate (arc seconds/year), phase (degrees).
 
   real   (r8), parameter :: obamp(poblen) =  & ! amplitudes for obliquity cos series
          (/   -2462.2214466_r8, -857.3232075_r8, -629.3231835_r8,   &
                -414.2804924_r8, -311.7632587_r8,  308.9408604_r8,   &
                -162.5533601_r8, -116.1077911_r8,  101.1189923_r8,   &
                 -67.6856209_r8,   24.9079067_r8,   22.5811241_r8,   &
                 -21.1648355_r8,  -15.6549876_r8,   15.3936813_r8,   &
                  14.6660938_r8,  -11.7273029_r8,   10.2742696_r8,   &
                   6.4914588_r8,    5.8539148_r8,   -5.4872205_r8,   &
                  -5.4290191_r8,    5.1609570_r8,    5.0786314_r8,   &
                  -4.0735782_r8,    3.7227167_r8,    3.3971932_r8,   &
                  -2.8347004_r8,   -2.6550721_r8,   -2.5717867_r8,   &
                  -2.4712188_r8,    2.4625410_r8,    2.2464112_r8,   &
                  -2.0755511_r8,   -1.9713669_r8,   -1.8813061_r8,   &
                  -1.8468785_r8,    1.8186742_r8,    1.7601888_r8,   &
                  -1.5428851_r8,    1.4738838_r8,   -1.4593669_r8,   &
                   1.4192259_r8,   -1.1818980_r8,    1.1756474_r8,   &
                  -1.1316126_r8,    1.0896928_r8/)
 
   real   (r8), parameter :: obrate(poblen) = & ! rates for obliquity cosine series
            (/  31.609974_r8, 32.620504_r8, 24.172203_r8,   &
                31.983787_r8, 44.828336_r8, 30.973257_r8,   &
                43.668246_r8, 32.246691_r8, 30.599444_r8,   &
                42.681324_r8, 43.836462_r8, 47.439436_r8,   &
                63.219948_r8, 64.230478_r8,  1.010530_r8,   &
                 7.437771_r8, 55.782177_r8,  0.373813_r8,   &
                13.218362_r8, 62.583231_r8, 63.593761_r8,   &
                76.438310_r8, 45.815258_r8,  8.448301_r8,   &
                56.792707_r8, 49.747842_r8, 12.058272_r8,   &
                75.278220_r8, 65.241008_r8, 64.604291_r8,   &
                 1.647247_r8,  7.811584_r8, 12.207832_r8,   &
                63.856665_r8, 56.155990_r8, 77.448840_r8,   &
                 6.801054_r8, 62.209418_r8, 20.656133_r8,   &
                48.344406_r8, 55.145460_r8, 69.000539_r8,   &
                11.071350_r8, 74.291298_r8, 11.047742_r8,   &
                 0.636717_r8, 12.844549_r8/)
 
   real   (r8), parameter :: obphas(poblen) = & ! phases for obliquity cosine series
          (/    251.9025_r8, 280.8325_r8, 128.3057_r8,   &
                292.7252_r8,  15.3747_r8, 263.7951_r8,   &
                308.4258_r8, 240.0099_r8, 222.9725_r8,   &
                268.7809_r8, 316.7998_r8, 319.6024_r8,   &
                143.8050_r8, 172.7351_r8,  28.9300_r8,   &
                123.5968_r8,  20.2082_r8,  40.8226_r8,   &
                123.4722_r8, 155.6977_r8, 184.6277_r8,   &
                267.2772_r8,  55.0196_r8, 152.5268_r8,   &
                 49.1382_r8, 204.6609_r8,  56.5233_r8,   &
                200.3284_r8, 201.6651_r8, 213.5577_r8,   &
                 17.0374_r8, 164.4194_r8,  94.5422_r8,   &
                131.9124_r8,  61.0309_r8, 296.2073_r8,   &
                135.4894_r8, 114.8750_r8, 247.0691_r8,   &
                256.6114_r8,  32.1008_r8, 143.6804_r8,   &
                 16.8784_r8, 160.6835_r8,  27.5932_r8,   &
                348.1074_r8,  82.6496_r8/)
 
   ! Cosine/sine series data for computation of eccentricity and fixed vernal 
   ! equinox longitude of perihelion (fvelp): amplitude, 
   ! rate (arc seconds/year), phase (degrees).
 
   real   (r8), parameter :: ecamp (pecclen) = & ! ampl for eccen/fvelp cos/sin series
          (/   0.01860798_r8,  0.01627522_r8, -0.01300660_r8,   &
               0.00988829_r8, -0.00336700_r8,  0.00333077_r8,   &
              -0.00235400_r8,  0.00140015_r8,  0.00100700_r8,   &
               0.00085700_r8,  0.00064990_r8,  0.00059900_r8,   &
               0.00037800_r8, -0.00033700_r8,  0.00027600_r8,   &
               0.00018200_r8, -0.00017400_r8, -0.00012400_r8,   &
               0.00001250_r8/)
 
   real   (r8), parameter :: ecrate(pecclen) = & ! rates for eccen/fvelp cos/sin series
          (/    4.2072050_r8,  7.3460910_r8, 17.8572630_r8,  &
               17.2205460_r8, 16.8467330_r8,  5.1990790_r8,  &
               18.2310760_r8, 26.2167580_r8,  6.3591690_r8,  &
               16.2100160_r8,  3.0651810_r8, 16.5838290_r8,  &
               18.4939800_r8,  6.1909530_r8, 18.8677930_r8,  &
               17.4255670_r8,  6.1860010_r8, 18.4174410_r8,  &
                0.6678630_r8/)
 
   real   (r8), parameter :: ecphas(pecclen) = & ! phases for eccen/fvelp cos/sin series
          (/    28.620089_r8, 193.788772_r8, 308.307024_r8,  &
               320.199637_r8, 279.376984_r8,  87.195000_r8,  &
               349.129677_r8, 128.443387_r8, 154.143880_r8,  &
               291.269597_r8, 114.860583_r8, 332.092251_r8,  &
               296.414411_r8, 145.769910_r8, 337.237063_r8,  &
               152.092288_r8, 126.839891_r8, 210.667199_r8,  &
                72.108838_r8/)
 
   ! Sine series data for computation of moving vernal equinox longitude of 
   ! perihelion: amplitude (arc seconds), rate (arc sec/year), phase (degrees).      
 
   real   (r8), parameter :: mvamp (pmvelen) = & ! amplitudes for mvelp sine series 
          (/   7391.0225890_r8, 2555.1526947_r8, 2022.7629188_r8,  &
              -1973.6517951_r8, 1240.2321818_r8,  953.8679112_r8,  &
               -931.7537108_r8,  872.3795383_r8,  606.3544732_r8,  &
               -496.0274038_r8,  456.9608039_r8,  346.9462320_r8,  &
               -305.8412902_r8,  249.6173246_r8, -199.1027200_r8,  &
                191.0560889_r8, -175.2936572_r8,  165.9068833_r8,  &
                161.1285917_r8,  139.7878093_r8, -133.5228399_r8,  &
                117.0673811_r8,  104.6907281_r8,   95.3227476_r8,  &
                 86.7824524_r8,   86.0857729_r8,   70.5893698_r8,  &
                -69.9719343_r8,  -62.5817473_r8,   61.5450059_r8,  &
                -57.9364011_r8,   57.1899832_r8,  -57.0236109_r8,  &
                -54.2119253_r8,   53.2834147_r8,   52.1223575_r8,  &
                -49.0059908_r8,  -48.3118757_r8,  -45.4191685_r8,  &
                -42.2357920_r8,  -34.7971099_r8,   34.4623613_r8,  &
                -33.8356643_r8,   33.6689362_r8,  -31.2521586_r8,  &
                -30.8798701_r8,   28.4640769_r8,  -27.1960802_r8,  &
                 27.0860736_r8,  -26.3437456_r8,   24.7253740_r8,  &
                 24.6732126_r8,   24.4272733_r8,   24.0127327_r8,  &
                 21.7150294_r8,  -21.5375347_r8,   18.1148363_r8,  &
                -16.9603104_r8,  -16.1765215_r8,   15.5567653_r8,  &
                 15.4846529_r8,   15.2150632_r8,   14.5047426_r8,  &
                -14.3873316_r8,   13.1351419_r8,   12.8776311_r8,  &
                 11.9867234_r8,   11.9385578_r8,   11.7030822_r8,  &
                 11.6018181_r8,  -11.2617293_r8,  -10.4664199_r8,  &
                 10.4333970_r8,  -10.2377466_r8,   10.1934446_r8,  &
                -10.1280191_r8,   10.0289441_r8,  -10.0034259_r8/)
 
   real   (r8), parameter :: mvrate(pmvelen) = & ! rates for mvelp sine series 
          (/    31.609974_r8, 32.620504_r8, 24.172203_r8,   &
                 0.636717_r8, 31.983787_r8,  3.138886_r8,   &
                30.973257_r8, 44.828336_r8,  0.991874_r8,   &
                 0.373813_r8, 43.668246_r8, 32.246691_r8,   &
                30.599444_r8,  2.147012_r8, 10.511172_r8,   &
                42.681324_r8, 13.650058_r8,  0.986922_r8,   &
                 9.874455_r8, 13.013341_r8,  0.262904_r8,   &
                 0.004952_r8,  1.142024_r8, 63.219948_r8,   &
                 0.205021_r8,  2.151964_r8, 64.230478_r8,   &
                43.836462_r8, 47.439436_r8,  1.384343_r8,   &
                 7.437771_r8, 18.829299_r8,  9.500642_r8,   &
                 0.431696_r8,  1.160090_r8, 55.782177_r8,   &
                12.639528_r8,  1.155138_r8,  0.168216_r8,   &
                 1.647247_r8, 10.884985_r8,  5.610937_r8,   &
                12.658184_r8,  1.010530_r8,  1.983748_r8,   &
                14.023871_r8,  0.560178_r8,  1.273434_r8,   &
                12.021467_r8, 62.583231_r8, 63.593761_r8,   &
                76.438310_r8,  4.280910_r8, 13.218362_r8,   &
                17.818769_r8,  8.359495_r8, 56.792707_r8,   &
                8.448301_r8,  1.978796_r8,  8.863925_r8,   &
                 0.186365_r8,  8.996212_r8,  6.771027_r8,   &
                45.815258_r8, 12.002811_r8, 75.278220_r8,   &
                65.241008_r8, 18.870667_r8, 22.009553_r8,   &
                64.604291_r8, 11.498094_r8,  0.578834_r8,   &
                 9.237738_r8, 49.747842_r8,  2.147012_r8,   &
                 1.196895_r8,  2.133898_r8,  0.173168_r8/)

   real   (r8), parameter :: mvphas(pmvelen) = & ! phases for mvelp sine series
          (/    251.9025_r8, 280.8325_r8, 128.3057_r8,   &
                348.1074_r8, 292.7252_r8, 165.1686_r8,   &
                263.7951_r8,  15.3747_r8,  58.5749_r8,   &
                 40.8226_r8, 308.4258_r8, 240.0099_r8,   &
                222.9725_r8, 106.5937_r8, 114.5182_r8,   &
                268.7809_r8, 279.6869_r8,  39.6448_r8,   &
                126.4108_r8, 291.5795_r8, 307.2848_r8,   &
                 18.9300_r8, 273.7596_r8, 143.8050_r8,   &
                191.8927_r8, 125.5237_r8, 172.7351_r8,   &
                316.7998_r8, 319.6024_r8,  69.7526_r8,   &
                123.5968_r8, 217.6432_r8,  85.5882_r8,   &
                156.2147_r8,  66.9489_r8,  20.2082_r8,   &
                250.7568_r8,  48.0188_r8,   8.3739_r8,   &
                 17.0374_r8, 155.3409_r8,  94.1709_r8,   &
                221.1120_r8,  28.9300_r8, 117.1498_r8,   &
                320.5095_r8, 262.3602_r8, 336.2148_r8,   &
                233.0046_r8, 155.6977_r8, 184.6277_r8,   &
                267.2772_r8,  78.9281_r8, 123.4722_r8,   &
                188.7132_r8, 180.1364_r8,  49.1382_r8,   &
                152.5268_r8,  98.2198_r8,  97.4808_r8,   &
                221.5376_r8, 168.2438_r8, 161.1199_r8,   &
                 55.0196_r8, 262.6495_r8, 200.3284_r8,   &
                201.6651_r8, 294.6547_r8,  99.8233_r8,   &
                213.5577_r8, 154.1631_r8, 232.7153_r8,   &
                138.3034_r8, 204.6609_r8, 106.5938_r8,   &
                250.4676_r8, 332.3345_r8,  27.3039_r8/)
 
   integer(i4) :: i       ! Index for series summations
   real   (r8) :: obsum   ! Obliquity series summation
   real   (r8) :: cossum  ! Cos series summation for eccentricity/fvelp
   real   (r8) :: sinsum  ! Sin series summation for eccentricity/fvelp
   real   (r8) :: fvelp   ! Fixed vernal equinox long of perihelion
   real   (r8) :: mvsum   ! mvelp series summation
   real   (r8) :: beta    ! Intermediate argument for lambm0
   real   (r8) :: years   ! Years to time of interest ( pos <=> future)
   real   (r8) :: eccen2  ! eccentricity squared
   real   (r8) :: eccen3  ! eccentricity cubed

   !-------------------------- Formats -----------------------------------------
   character(*),parameter :: svnID  = "SVN " // &
   "$Id: shr_orb_mod.F90 25434 2010-11-04 22:46:24Z tcraig $"
   character(*),parameter :: svnURL = "SVN <unknown URL>" 
!  character(*),parameter :: svnURL = "SVN " // &
!  "$URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/release_tags/cesm1_2_x_n02_share3_130715/shr/shr_orb_mod.F90 $"
   character(len=*),parameter :: F00 = "('(orb_params) ',4a)"
   character(len=*),parameter :: F01 = "('(orb_params) ',a,i9)"
   character(len=*),parameter :: F02 = "('(orb_params) ',a,f6.3)"
   character(len=*),parameter :: F03 = "('(orb_params) ',a,es14.6)"

   !----------------------------------------------------------------------------
   ! radinp and algorithms below will need a degree to radian conversion factor
 
   if ( log_print .and. mpi_rank()==0 ) then
     write(*,F00) 'Calculate characteristics of the orbit:'
!     write(*,F00) svnID
!    write(*,F00) svnURL
   end if
 
   IF ( iyear == SHR_ORB_UNDEF_INT ) THEN

      ! Check input obliq, eccen, and mvelp to ensure reasonable
 
      if( obliq == SHR_ORB_UNDEF_REAL )then
         write(*,F00) trim(subname)//' Have to specify orbital parameters:'
         write(*,F00) 'Either set: iyear, OR [obliq, eccen, and mvelp]:'
         write(*,F00) 'iyear is the year to simulate orbit for (ie. 1950): '
         write(*,F00) 'obliq, eccen, mvelp specify the orbit directly:'
         write(*,F00) 'The AMIP II settings (for a 1995 orbit) are: '
         write(*,F00) ' obliq =  23.4441'
         write(*,F00) ' eccen =   0.016715'
         write(*,F00) ' mvelp = 102.7'
         call endrun(subname//' ERROR: unreasonable obliq')
      else if ( log_print .and. mpi_rank()==0 ) then
         write(*,F00) 'Use input orbital parameters: '
      end if
      if( (obliq < SHR_ORB_OBLIQ_MIN).or.(obliq > SHR_ORB_OBLIQ_MAX) ) then
         write(*,F03) 'Input obliquity unreasonable: ', obliq
         call endrun(subname//' ERROR: unreasonable obliq')
      end if
      if( (eccen < SHR_ORB_ECCEN_MIN).or.(eccen > SHR_ORB_ECCEN_MAX) ) then
         write(*,F03) 'Input eccentricity unreasonable: ', eccen
         call endrun(subname//' ERROR: unreasonable eccen')
      end if
      if( (mvelp < SHR_ORB_MVELP_MIN).or.(mvelp > SHR_ORB_MVELP_MAX) ) then
         write(*,F03) 'Input mvelp unreasonable: ' , mvelp
         call endrun(subname//' ERROR: unreasonable mvelp')
      end if
      eccen2 = eccen*eccen
      eccen3 = eccen2*eccen

   ELSE  ! Otherwise calculate based on years before present
 
      if ( log_print .and. mpi_rank()==0 ) then
         write(*,F01) 'Calculate orbit for year: ' , iyear
      end if
      yb4_1950AD = 1950.0_r8 - real(iyear,r8)
      if ( abs(yb4_1950AD) .gt. 1000000.0_r8 )then
         write(*,F00) 'orbit only valid for years+-1000000'
         write(*,F00) 'Relative to 1950 AD'
         write(*,F03) '# of years before 1950: ',yb4_1950AD
         write(*,F01) 'Year to simulate was  : ',iyear
         call endrun(subname//' ERROR: unreasonable year')
      end if
 
      ! The following calculates the earths obliquity, orbital eccentricity
      ! (and various powers of it) and vernal equinox mean longitude of
      ! perihelion for years in the past (future = negative of years past),
      ! using constants (see parameter section) given in the program of:
      !
      ! Berger, Andre.  1978  A Simple Algorithm to Compute Long-Term Variations
      ! of Daily Insolation.  Contribution 18, Institute of Astronomy and
      ! Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium.
      !
      ! and formulas given in the paper (where less precise constants are also
      ! given):
      !
      ! Berger, Andre.  1978.  Long-Term Variations of Daily Insolation and
      ! Quaternary Climatic Changes.  J. of the Atmo. Sci. 35:2362-2367
      !
      ! The algorithm is valid only to 1,000,000 years past or hence.
      ! For a solution valid to 5-10 million years past see the above author.
      ! Algorithm below is better for years closer to present than is the
      ! 5-10 million year solution.
      !
      ! Years to time of interest must be negative of years before present
      ! (1950) in formulas that follow. 
 
      years = - yb4_1950AD
 
      ! In the summations below, cosine or sine arguments, which end up in
      ! degrees, must be converted to radians via multiplication by degrad.
      !
      ! Summation of cosine series for obliquity (epsilon in Berger 1978) in
      ! degrees. Convert the amplitudes and rates, which are in arc secs, into
      ! degrees via multiplication by psecdeg (arc seconds to degrees conversion
      ! factor).  For obliq, first term is Berger 1978 epsilon star; second
      ! term is series summation in degrees.
 
      obsum = 0.0_r8
      do i = 1, poblen
         obsum = obsum + obamp(i)*psecdeg*cos((obrate(i)*psecdeg*years + obphas(i))*degrad)
      end do
      obliq = 23.320556_r8 + obsum
 
      ! Summation of cosine and sine series for computation of eccentricity 
      ! (eccen; e in Berger 1978) and fixed vernal equinox longitude of 
      ! perihelion (fvelp; pi in Berger 1978), which is used for computation 
      ! of moving vernal equinox longitude of perihelion.  Convert the rates, 
      ! which are in arc seconds, into degrees via multiplication by psecdeg.
 
      cossum = 0.0_r8
      do i = 1, pecclen
        cossum = cossum+ecamp(i)*cos((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do
 
      sinsum = 0.0_r8
      do i = 1, pecclen
        sinsum = sinsum+ecamp(i)*sin((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do
 
      ! Use summations to calculate eccentricity
 
      eccen2 = cossum*cossum + sinsum*sinsum
      eccen  = sqrt(eccen2)
      eccen3 = eccen2*eccen
 
      ! A series of cases for fvelp, which is in radians.
         
      if (abs(cossum) .le. 1.0E-8_r8) then
        if (sinsum .eq. 0.0_r8) then
          fvelp = 0.0_r8
        else if (sinsum .lt. 0.0_r8) then
          fvelp = 1.5_r8*pi
        else if (sinsum .gt. 0.0_r8) then
          fvelp = .5_r8*pi
        endif
      else if (cossum .lt. 0.0_r8) then
        fvelp = atan(sinsum/cossum) + pi
      else if (cossum .gt. 0.0_r8) then
        if (sinsum .lt. 0.0_r8) then
          fvelp = atan(sinsum/cossum) + 2.0_r8*pi
        else
          fvelp = atan(sinsum/cossum)
        endif
      endif
 
      ! Summation of sin series for computation of moving vernal equinox long
      ! of perihelion (mvelp; omega bar in Berger 1978) in degrees.  For mvelp,
      ! first term is fvelp in degrees; second term is Berger 1978 psi bar 
      ! times years and in degrees; third term is Berger 1978 zeta; fourth 
      ! term is series summation in degrees.  Convert the amplitudes and rates,
      ! which are in arc seconds, into degrees via multiplication by psecdeg.  
      ! Series summation plus second and third terms constitute Berger 1978
      ! psi, which is the general precession.
 
      mvsum = 0.0_r8
      do i = 1, pmvelen
        mvsum = mvsum + mvamp(i)*psecdeg*sin((mvrate(i)*psecdeg*years + mvphas(i))*degrad)
      end do
      mvelp = fvelp/degrad + 50.439273_r8*psecdeg*years + 3.392506_r8 + mvsum
 
      ! Cases to make sure mvelp is between 0 and 360.
 
      do while (mvelp .lt. 0.0_r8)
        mvelp = mvelp + 360.0_r8
      end do
      do while (mvelp .ge. 360.0_r8)
        mvelp = mvelp - 360.0_r8
      end do

   END IF  ! end of test on whether to calculate or use input orbital params
 
   ! Orbit needs the obliquity in radians
 
   obliqr = obliq*degrad
 
   ! 180 degrees must be added to mvelp since observations are made from the
   ! earth and the sun is considered (wrongly for the algorithm) to go around
   ! the earth. For a more graphic explanation see Appendix B in:
   !
   ! A. Berger, M. Loutre and C. Tricot. 1993.  Insolation and Earth Orbital
   ! Periods.  J. of Geophysical Research 98:10,341-10,362.
   !
   ! Additionally, orbit will need this value in radians. So mvelp becomes
   ! mvelpp (mvelp plus pi)
    mvelpp = (mvelp + 180._r8)*degrad
 
   ! Set up an argument used several times in lambm0 calculation ahead.
 
   beta = sqrt(1._r8 - eccen2)
 
   ! The mean longitude at the vernal equinox (lambda m nought in Berger
   ! 1978; in radians) is calculated from the following formula given in 
   ! Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
   ! 1978) is 0.

   lambm0 = 2._r8*((.5_r8*eccen + .125_r8*eccen3)*(1._r8 + beta)*sin(mvelpp)  &
         - .250_r8*eccen2*(.5_r8    + beta)*sin(2._r8*mvelpp)            &
         + .125_r8*eccen3*(1._r8/3._r8 + beta)*sin(3._r8*mvelpp))
 
   if ( log_print .and. mpi_rank()==0) then
     write(*,F03) '------ Computed Orbital Parameters ------'
     write(*,F03) 'Eccentricity      = ',eccen
     write(*,F03) 'Obliquity (deg)   = ',obliq
     write(*,F03) 'Obliquity (rad)   = ',obliqr
     write(*,F03) 'Long of perh(deg) = ',mvelp
     write(*,F03) 'Long of perh(rad) = ',mvelpp
     write(*,F03) 'Long at v.e.(rad) = ',lambm0
     write(*,F03) '-----------------------------------------'
   end if
 
   end subroutine orb_params


! Compute earth/orbit parameters using formula suggested by Duane Thresher.
   subroutine orb_decl(calday , delta , eccf, doAquaPlanet)

   real   (r8),intent(in)  :: calday ! Calendar day, including fraction
   real   (r8),intent(out) :: delta  ! Solar declination angle in rad
   real   (r8),intent(out) :: eccf   ! Earth-sun distance factor (ie. (1/r)**2)
   logical, optional, intent(in)  :: doAquaPlanet  ! yiz added for doAquaPlanet
 
   !---------------------------Local variables-----------------------------
   real   (r8),parameter :: dayspy = 365.0_r8  ! days per year
   real   (r8),parameter :: ve     = 80.5_r8   ! Calday of vernal equinox
                                               ! assumes Jan 1 = calday 1
 
   real   (r8) ::   lambm  ! Lambda m, mean long of perihelion (rad)
   real   (r8) ::   lmm    ! Intermediate argument involving lambm
   real   (r8) ::   lamb   ! Lambda, the earths long of perihelion
   real   (r8) ::   invrho ! Inverse normalized sun/earth distance
   real   (r8) ::   sinl   ! Sine of lmm
 
   ! Compute eccentricity factor and solar declination using
   ! day value where a round day (such as 213.0) refers to 0z at
   ! Greenwich longitude.
   !
   ! Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
   ! Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
   ! 35:2362-2367.
   !
   ! To get the earths true longitude (position in orbit; lambda in Berger 
   ! 1978) which is necessary to find the eccentricity factor and declination,
   ! must first calculate the mean longitude (lambda m in Berger 1978) at
   ! the present day.  This is done by adding to lambm0 (the mean longitude
   ! at the vernal equinox, set as March 21 at noon, when lambda=0; in radians)
   ! an increment (delta lambda m in Berger 1978) that is the number of
   ! days past or before (a negative increment) the vernal equinox divided by
   ! the days in a model year times the 2*pi radians in a complete orbit.
 
   lambm = lambm0 + (calday - ve)*2._r8*pi/dayspy
   lmm   = lambm  - mvelpp
 
   ! The earths true longitude, in radians, is then found from
   ! the formula in Berger 1978:
 
   sinl  = sin(lmm)
   lamb  = lambm  + eccen*(2._r8*sinl + eccen*(1.25_r8*sin(2._r8*lmm)  &
         + eccen*((13.0_r8/12.0_r8)*sin(3._r8*lmm) - 0.25_r8*sinl)))
 
   ! Using the obliquity, eccentricity, moving vernal equinox longitude of
   ! perihelion (plus), and earths true longitude, the declination (delta)
   ! and the normalized earth/sun distance (rho in Berger 1978; actually inverse
   ! rho will be used), and thus the eccentricity factor (eccf), can be 
   ! calculated from formulas given in Berger 1978.
 
   invrho = (1._r8 + eccen*cos(lamb - mvelpp)) / (1._r8 - eccen*eccen)
 
   ! Set solar declination and eccentricity factor
 
   delta  = asin(sin(obliqr)*sin(lamb))
   eccf   = invrho*invrho
! perprtual equinox
   if(present(doAquaPlanet).and.doAquaPlanet)then
      delta = zero
      eccf  = one ! eccen is zero
   end if

   return
 
   end subroutine orb_decl

 end module grist_zenith