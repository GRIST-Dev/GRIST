!wrf:model_layer:constants
!

 module module_wrfmodel_constants

   !  2. following are constants for use in defining real number bounds.

   !  a really small number.

   real    , parameter :: epsilon         = 1.e-15

   !  4. following is information related to the physical constants.

   !  these are the physical constants used within the model.

! jm note -- can we name this grav instead?
   real    , parameter :: g = 9.81  ! acceleration due to gravity (m {s}^-2)

#if ( nmm_core == 1 )
   real    , parameter :: r_d          = 287.04
   real    , parameter :: cp           = 1004.6
#else
   real    , parameter :: r_d          = 287.
   real    , parameter :: cp           = 7.*r_d/2.
   real    , parameter :: cpd          = 7.*r_d/2.
#endif
   real    , parameter :: r_v          = 461.6
   real    , parameter :: cpv          = 4.*r_v
   real    , parameter :: cv           = cp-r_d
   real    , parameter :: cvd          = cp-r_d
   real    , parameter :: cvv          = cpv-r_v
   real    , parameter :: cvpm         = -cv/cp
   real    , parameter :: cliq         = 4190.
   real    , parameter :: cice         = 2106.
   real    , parameter :: psat         = 610.78
   real    , parameter :: rcv          = r_d/cv
   real    , parameter :: rcp          = r_d/cp
   real    , parameter :: rovg         = r_d/g
   real    , parameter :: c2           = cp * rcv

   real    , parameter :: p1000mb      = 100000.
   real    , parameter :: t0           = 300.
   real    , parameter :: p0           = p1000mb
   real    , parameter :: cpovcv       = cp/(cp-r_d)
   real    , parameter :: cvovcp       = 1./cpovcv
   real    , parameter :: rvovrd       = r_v/r_d

   real    , parameter :: reradius     = 1./6370.e03 

   real    , parameter :: asselin      = .025
!   real    , parameter :: asselin      = .0
   real    , parameter :: cb           = 25.

   real    , parameter :: xlv0         = 3.15e6
   real    , parameter :: xlv1         = 2370.
   real    , parameter :: xls0         = 2.905e6
   real    , parameter :: xls1         = 259.532

   real    , parameter :: xls          = 2.85e6
   real    , parameter :: xlv          = 2.5e6
   real    , parameter :: xlf          = 3.50e5

   real    , parameter :: rhowater     = 1000.
   real    , parameter :: rhosnow      = 100.
   real    , parameter :: rhoair0      = 1.28

   real    , parameter :: degrad       = 3.1415926/180.
   real    , parameter :: dpd          = 360./365.

   real    , parameter ::  svp1=0.6112
   real    , parameter ::  svp2=17.67
   real    , parameter ::  svp3=29.65
   real    , parameter ::  svpt0=273.15
   real    , parameter ::  ep_1=r_v/r_d-1.
   real    , parameter ::  ep_2=r_d/r_v
   real    , parameter ::  karman=0.4
   real    , parameter ::  eomeg=7.2921e-5
   real    , parameter ::  stbolt=5.67051e-8

                                      ! proportionality constants for eddy viscosity coefficient calc
   real    , parameter ::  c_s = .25  ! turbulence parameterization constant, for smagorinsky
   real    , parameter ::  c_k = .25  ! turbulence parameterization constant, for tke
   real    , parameter ::  prandtl = 1./3.0
                                         ! constants for w-damping option
   real    , parameter ::  w_alpha = 0.3 ! strength m/s/s
   real    , parameter ::  w_beta  = 1.0 ! activation cfl number

       real , parameter ::  pq0=379.90516
       real , parameter ::  epsq2=0.2
       real , parameter ::  a2=17.2693882
       real , parameter ::  a3=273.16
       real , parameter ::  a4=35.86
       real , parameter ::  epsq=1.e-12
       real , parameter ::  p608=rvovrd-1.
#if ( nmm_core == 1 )
       real , parameter ::  climit=1.e-20
       real , parameter ::  cm1=2937.4
       real , parameter ::  cm2=4.9283
       real , parameter ::  cm3=23.5518
       real , parameter ::  defc=8.0
       real , parameter ::  defm=32.0
       real , parameter ::  epsfc=1./1.05
       real , parameter ::  epswet=0.0
       real , parameter ::  fcdif=1./3.
       real , parameter ::  fcm=0.003
       real , parameter ::  gma=-r_d*(1.-rcp)*0.5
       real , parameter ::  p400=40000.0
       real , parameter ::  phitp=15000.0
       real , parameter ::  pi2=2.*3.1415926
       real , parameter ::  plbtm=105000.0
       real , parameter ::  plomd=64200.0
       real , parameter ::  pmdhi=35000.0
       real , parameter ::  q2ini=0.50
       real , parameter ::  rfcp=0.25/cp
       real , parameter ::  rhcrit_land=0.75
       real , parameter ::  rhcrit_sea=0.80
       real , parameter ::  rlag=14.8125
       real , parameter ::  rlx=0.90
       real , parameter ::  scq2=50.0
       real , parameter ::  slopht=0.001
       real , parameter ::  tlc=2.*0.703972477
       real , parameter ::  wa=0.15
       real , parameter ::  wght=0.35
       real , parameter ::  wpc=0.075
       real , parameter ::  z0land=0.10
       real , parameter ::  z0max=0.01
       real , parameter ::  z0sea=0.001
#endif


 contains
   subroutine init_module_model_constants
   end subroutine init_module_model_constants
 end module module_wrfmodel_constants
