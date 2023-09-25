 module cldwat2m_macro_lin

   use grist_constants,                only: i4, r8,                    &
                                             cp, latvap, latice, rdry
   use cloud_fraction,                 only: top_lev => trop_cloud_top_lev 
   use grist_wv_saturation,            only: qsat_water, svp_water, svp_ice
   use grist_handle_error,             only: endrun
   use grist_mpi

   !---------------------Lin Macrop-------------------------
   use cldfrac_conden_mod, only: subgrid_var,cldfrac_conden_single_only
   !---------------------Lin Macrop-------------------------

   implicit none
   private
   save

   public :: ini_macro_lin, mmacro_pcond_lin   

! Parameters for Ice Stratus
   real(r8), public,  parameter :: rhmini       = 0.80_r8       ! Minimum rh for ice cloud fraction > 0.
   real(r8), public,  parameter :: rhmaxi       = 1.1_r8        ! rhi at which ice cloud fraction = 1.
   real(r8), parameter          :: qist_min     = 1.e-6_r8      ! default: 1.e-7_r8, LiXH, Minimum in-stratus ice IWC constraint [ kg/kg ]
   real(r8), parameter          :: qist_max     = 5.e-3_r8      ! Maximum in-stratus ice IWC constraint [ kg/kg ]

                                                                ! use Slingo (triangular PDF-based) liquid stratus fraction
   logical,  parameter          :: freeze_dry   = .false.       ! If .true., use 'freeze dry' in liquid stratus fraction formula
   real(r8), parameter          :: qlst_min     = 2.e-5_r8      ! Minimum in-stratus LWC constraint [ kg/kg ]
   real(r8), parameter          :: qlst_max     = 3.e-3_r8      ! Maximum in-stratus LWC constraint [ kg/kg ]
   real(r8), parameter          :: cc           = 0.1_r8        ! For newly formed/dissipated in-stratus CWC ( 0 <= cc <= 1 )
   real(r8), parameter          :: ramda        = 0.5_r8        ! Explicit : ramda = 0, Implicit : ramda = 1 ( 0<= ramda <= 1 )
   real(r8), private            :: rhminl                       ! Critical RH for low-level  liquid stratus clouds
   real(r8), private            :: rhminl_adj_land              ! rhminl adjustment for snowfree land
   real(r8), private            :: rhminh                       ! Critical RH for high-level liquid stratus clouds
   real(r8), private            :: premit                       ! Top    height for mid-level liquid stratus fraction
   real(r8), private            :: premib                       ! Bottom height for mid-level liquid stratus fraction
   integer,  private            :: iceopt                       ! option for ice cloud closure 
                                                                ! 1=wang & sassen 2=schiller (iciwc)  
                                                                ! 3=wood & field, 4=Wilson (based on smith)
                                                                ! 5=modified slingo (ssat & empyt cloud)        
   real(r8), private            :: icecrit                      ! Critical RH for ice clouds in Wilson & Ballard closure
                                                                ! ( smaller = more ice clouds )

   integer :: pver, pverp, pcnst


   contains

! Purpose: Initialize constants for the liquid stratiform macrophysics 
!          called from stratiform.F90.
   subroutine ini_macro_lin

   use cloud_fraction,   only: cldfrc_getparams
   use grist_nml_module, only: nlev, nlevp, ntracer

   pver  = nlev
   pverp = nlevp
   pcnst = ntracer

   call cldfrc_getparams(rhminl_out = rhminl, rhminl_adj_land_out = rhminl_adj_land,    &
                         rhminh_out = rhminh, premit_out = premit,   premib_out = premib, &
                         iceopt_out = iceopt, icecrit_out = icecrit)

!   if(mpi_rank()==0)then
!       write(*,*) 'ini_macro: tuning parameters : rhminl = ', rhminl, &
!                                                 'rhminl_adj_land = ', rhminl_adj_land, & 
!                                                 'rhminh = ', rhminh, & 
!                                                 'premit = ', premit, & 
!                                                 'premib = ',  premib,  &
!                                                 'iceopt = ',  iceopt,  &
!                                                 'icecrit = ', icecrit
!   end if

   return
   end subroutine ini_macro_lin


! Purpose: Stratiform Liquid Macrophysics (Gauss-PDF)
   subroutine mmacro_pcond_lin( ncol       , dt         , p          , dp           ,                           &
                            T0         , qv0        , ql0        , qi0          , nl0        , ni0        , &
                            a_cu0      , landfrac   , snowh      ,              & 
                            s_tendout  , qv_tendout , ql_tendout , qi_tendout   , nl_tendout , ni_tendout , &
                            qme        ,  &
                            cld        , al_st_star , ai_st_star , ql_st_star   , qi_st_star , do_cldice  , &
                            rpdel, zm, zi,  &
                            lengi, shi, wstarPBL, &
                            qtu_shal, umf_shal, cnt_shal, cnb_shal, thlu_shal, &
                            alst_gp, cond_gp, N1, &
                            sgm, &
                            delta_q, deltaq_sat, deltaq_uns, Q1_sat, Q1_uns,&
                            adjust_factor)

   use grist_physics_data_structure,   only: phy_tracer_info
   use grist_wv_saturation,            only: findsp_vc

   implicit none

! io
   integer,  intent(in)    :: ncol                         ! Number of active columns

   real(r8), intent(in) :: rpdel(pver, ncol)
   real(r8), intent(in) :: zm(pver, ncol)
   real(r8), intent(in) :: zi(pverp, ncol)
 
   real(r8), intent(in) :: lengi(pverp, ncol)
   real(r8), intent(in) :: shi(pverp, ncol)
   real(r8), intent(in) :: wstarPBL(ncol)

   ! shallow
   real(r8), intent(in) :: qtu_shal(pverp, ncol)
   real(r8), intent(in) :: umf_shal(pverp, ncol)
   real(r8), intent(in) :: cnt_shal(ncol)
   real(r8), intent(in) :: cnb_shal(ncol)
   real(r8), intent(in) :: thlu_shal(pverp, ncol)

   ! Input-Output variables

   real(r8), intent(inout) :: T0(pver, ncol)               ! Temperature [K]
   real(r8), intent(inout) :: qv0(pver, ncol)              ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(inout) :: ql0(pver, ncol)              ! Grid-mean liquid water content [kg/kg]
   real(r8), intent(inout) :: qi0(pver, ncol)              ! Grid-mean ice water content [kg/kg]
   real(r8), intent(inout) :: nl0(pver, ncol)              ! Grid-mean number concentration of cloud liquid droplet [#/kg]
   real(r8), intent(inout) :: ni0(pver, ncol)              ! Grid-mean number concentration of cloud ice    droplet [#/kg]

   ! Input variables

   real(r8), intent(in)    :: dt                           ! Model integration time step [s]
   real(r8), intent(in)    :: p(pver, ncol)                ! Pressure at the layer mid-point [Pa]
   real(r8), intent(in)    :: dp(pver, ncol)               ! Pressure thickness [Pa] > 0

   !real(r8), intent(in)    :: a_cud(pver, ncol)            ! Old cumulus fraction before update
   real(r8), intent(in)    :: a_cu0(pver, ncol)            ! New cumulus fraction after update
   
   real(r8), intent(in)    :: landfrac(ncol)              ! Land fraction
   real(r8), intent(in)    :: snowh(ncol)                 ! Snow depth (liquid water equivalent)
   logical, intent(in)     :: do_cldice                    ! Whether or not cldice should be prognosed

   ! Output variables

   real(r8), intent(out)   :: s_tendout(pver, ncol)        ! Net tendency of grid-mean s  from 'Micro+Macro' processes [J/kg/s]
   real(r8), intent(out)   :: qv_tendout(pver, ncol)       ! Net tendency of grid-mean qv from 'Micro+Macro' processes [kg/kg/s]
   real(r8), intent(out)   :: ql_tendout(pver, ncol)       ! Net tendency of grid-mean ql from 'Micro+Macro' processes [kg/kg/s]
   real(r8), intent(out)   :: qi_tendout(pver, ncol)       ! Net tendency of grid-mean qi from 'Micro+Macro' processes [kg/kg/s]
   real(r8), intent(out)   :: nl_tendout(pver, ncol)       ! Net tendency of grid-mean nl from 'Micro+Macro' processes [#/kg/s]
   real(r8), intent(out)   :: ni_tendout(pver, ncol)       ! Net tendency of grid-mean ni from 'Micro+Macro' processes [#/kg/s]
   real(r8), intent(out)   :: qme  (pver, ncol)            ! Net condensation rate [kg/kg/s]
   real(r8), intent(out)   :: cld(pver, ncol)              ! Net cloud fraction ( 0 <= cld <= 1 )
   real(r8), intent(out)   :: al_st_star(pver, ncol)       ! Physical liquid stratus fraction
   real(r8), intent(out)   :: ai_st_star(pver, ncol)       ! Physical ice stratus fraction
   real(r8), intent(out)   :: ql_st_star(pver, ncol)       ! In-stratus LWC [kg/kg] 
   real(r8), intent(out)   :: qi_st_star(pver, ncol)       ! In-stratus IWC [kg/kg] 

   real(r8), intent(out) :: alst_gp(pver, ncol)            ! cloud fraction from Gauss-PDF
   real(r8), intent(out) :: cond_gp(pver, ncol)            ! cloud condensate from Gauss-PDF
   real(r8), intent(out) :: N1(pver, ncol)                 ! use to count the occurrence of over-saturated situation
 
   real(r8), intent(out) :: sgm(pver, ncol)
 
   real(r8), intent(out) :: delta_q(pver, ncol)
   real(r8), intent(out) :: deltaq_sat(pver, ncol)
   real(r8), intent(out) :: deltaq_uns(pver, ncol)
   real(r8), intent(out) :: Q1_sat(pver, ncol)
   real(r8), intent(out) :: Q1_uns(pver, ncol)
 
   real(r8), intent(out) :: adjust_factor(pver, ncol)

! local
   integer :: ixcldliq, ixcldice
   integer :: i, j, k, iter, ii, jj, m                     ! Loop indexes

   ! Thermodynamic state variables
   real(r8) T(pver, ncol)                                  ! Temperature of equilibrium reference state
                                                           ! from which 'Micro & Macro' are computed [K]
   real(r8) T_prime0(pver, ncol)                           ! Temperature after 'Macrophysics (QQ)' on T_05star

   real(r8) qv(pver, ncol)                                 ! Grid-mean qv of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) qv_prime0(pver, ncol)                          ! Grid-mean qv after 'Macrophysics (QQ)' on qv_05star

   real(r8) ql(pver, ncol)                                 ! Grid-mean ql of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) ql_prime0(pver, ncol)                          ! Grid-mean ql after 'Macrophysics (QQ)' on ql_05star

   real(r8) qi(pver, ncol)                                 ! Grid-mean qi of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) qi_prime0(pver, ncol)                          ! Grid-mean qi after 'Macrophysics (QQ)' on qi_05star

   real(r8) nl(pver, ncol)                                 ! Grid-mean nl of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) nl_prime0(pver, ncol)                          ! Grid-mean nl after 'Macrophysics (QQ)' on nl_05star

   real(r8) ni(pver, ncol)                                 ! Grid-mean ni of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) ni_prime0(pver, ncol)                          ! Grid-mean ni after 'Macrophysics (QQ)' on ni_05star

   real(r8) a_st_star(pver, ncol)                          ! Stratus fraction at '_star' state

   real(r8) al_st_nc(pver, ncol)                           ! Non-physical liquid stratus fraction in the non-cumulus pixels

   real(r8) ai_st_nc(pver, ncol)                           ! Non-physical ice stratus fraction in the non-cumulus pixels

 ! Adjustment tendency in association with 'positive_moisture'

   real(r8) Tten_pwi1(pver, ncol)                          ! Pre-process T  tendency of input equilibrium state [K/s] 
   real(r8) qvten_pwi1(pver, ncol)                         ! Pre-process qv tendency of input equilibrium state [kg/kg/s]
   real(r8) qlten_pwi1(pver, ncol)                         ! Pre-process ql tendency of input equilibrium state [kg/kg/s]
   real(r8) qiten_pwi1(pver, ncol)                         ! Pre-process qi tendency of input equilibrium state [kg/kg/s]
   real(r8) nlten_pwi1(pver, ncol)                         ! Pre-process nl tendency of input equilibrium state [#/kg/s]
   real(r8) niten_pwi1(pver, ncol)                         ! Pre-process ni tendency of input equilibrium state [#/kg/s] 

   real(r8) Tten_pwi2(pver, ncol)                          ! Post-process T  tendency of provisional equilibrium state [K/s] 
   real(r8) qvten_pwi2(pver, ncol)                         ! Post-process qv tendency of provisional equilibrium state [kg/kg/s]
   real(r8) qlten_pwi2(pver, ncol)                         ! Post-process ql tendency of provisional equilibrium state [kg/kg/s]
   real(r8) qiten_pwi2(pver, ncol)                         ! Post-process qi tendency of provisional equilibrium state [kg/kg/s]
   real(r8) nlten_pwi2(pver, ncol)                         ! Post-process nl tendency of provisoonal equilibrium state [#/kg/s]
   real(r8) niten_pwi2(pver, ncol)                         ! Post-process ni tendency of provisional equilibrium state [#/kg/s] 


 ! Macrophysical process tendency variables

   real(r8) QQ(pver, ncol)             ! Net condensation rate into water+ice           [kg/kg/s] 
   real(r8) QQw(pver, ncol)            ! Net condensation rate into water               [kg/kg/s] 
   real(r8) QQi(pver, ncol)            ! Net condensation rate into ice                 [kg/kg/s]
   real(r8) QQnl(pver, ncol)           ! Tendency of nl associated with QQw both for condensation and evaporation [#/kg/s]
   real(r8) QQni(pver, ncol)           ! Tendency of ni associated with QQi both for condensation and evaporation [#/kg/s]

 ! Coefficient for computing QQ and related processes

   real(r8) U_nc(pver, ncol)                               ! Mean RH of non-cumulus pixels
   real(r8) G_nc(pver, ncol)                               ! d(U_nc)/d(a_st_nc)

   real(r8) zeros(pver, ncol)

   real(r8) qmin1(pver, ncol)
   real(r8) qmin2(pver, ncol)
   real(r8) qmin3(pver, ncol)

   real(r8) esat_a(ncol)                                  ! Saturation water vapor pressure [Pa]
   real(r8) qsat_a(ncol)                                  ! Saturation water vapor specific humidity [kg/kg]
   real(r8) Twb_aw(pver, ncol)                             ! Wet-bulb temperature [K]
   real(r8) qvwb_aw(pver, ncol)                            ! Wet-bulb water vapor specific humidity [kg/kg]

   real(r8) QQmax,QQmin,QQwmin,QQimin                      ! For limiting QQ
   real(r8) cone                                           ! Number close to but smaller than 1
   real(r8) qsmall                                         ! Smallest mixing ratio considered in the macrophysics

   real(r8) conden_st_nc(pver, ncol)

   conden_st_nc(:,:ncol) = 0._r8
 
   deltaq_sat(:,:ncol) = 0._r8
   deltaq_uns(:,:ncol) = 0._r8
   Q1_sat(:,:ncol) = 0._r8
   Q1_uns(:,:ncol) = 0._r8
   N1(:,:ncol) = 0._r8
   adjust_factor(:,:ncol) = 0._r8
 
   cone            = 0.999_r8
   qsmall          = 1.e-18_r8
   zeros(:,:ncol)  = 0._r8

   ! ------------------------------------ !
   ! Global initialization of main output !
   ! ------------------------------------ !

     s_tendout(:,:ncol)     = 0._r8
     qv_tendout(:,:ncol)    = 0._r8
     ql_tendout(:,:ncol)    = 0._r8
     qi_tendout(:,:ncol)    = 0._r8
     nl_tendout(:,:ncol)    = 0._r8
     ni_tendout(:,:ncol)    = 0._r8

     qme(:,:ncol)           = 0._r8

     cld(:,:ncol)           = 0._r8
     al_st_star(:,:ncol)    = 0._r8
     ai_st_star(:,:ncol)    = 0._r8
     ql_st_star(:,:ncol)    = 0._r8
     qi_st_star(:,:ncol)    = 0._r8

   !---------------------Lin Macrop-------------------------
     alst_gp(:,:ncol)       = 0._r8
     cond_gp(:,:ncol)       = 0._r8
     delta_q(:,:ncol)       = 0._r8
   !---------------------Lin Macrop-------------------------
   
   ! --------------------------------------- !
   ! Initialization of internal 2D variables !
   ! --------------------------------------- !

     T(:,:ncol)             = 0._r8
     T_prime0(:,:ncol)      = 0._r8

     qv(:,:ncol)            = 0._r8
     qv_prime0(:,:ncol)     = 0._r8

     ql(:,:ncol)            = 0._r8
     ql_prime0(:,:ncol)     = 0._r8

     qi(:,:ncol)            = 0._r8
     qi_prime0(:,:ncol)     = 0._r8

     nl(:,:ncol)            = 0._r8
     nl_prime0(:,:ncol)     = 0._r8

     ni(:,:ncol)            = 0._r8
     ni_prime0(:,:ncol)     = 0._r8

     a_st_star(:,:ncol)     = 0._r8
     al_st_nc(:,:ncol)      = 0._r8
     ai_st_nc(:,:ncol)      = 0._r8

 ! Adjustment tendency in association with 'positive_moisture'

     Tten_pwi1(:,:ncol)     = 0._r8
     qvten_pwi1(:,:ncol)    = 0._r8
     qlten_pwi1(:,:ncol)    = 0._r8
     qiten_pwi1(:,:ncol)    = 0._r8
     nlten_pwi1(:,:ncol)    = 0._r8
     niten_pwi1(:,:ncol)    = 0._r8

     Tten_pwi2(:,:ncol)     = 0._r8
     qvten_pwi2(:,:ncol)    = 0._r8
     qlten_pwi2(:,:ncol)    = 0._r8
     qiten_pwi2(:,:ncol)    = 0._r8
     nlten_pwi2(:,:ncol)    = 0._r8
     niten_pwi2(:,:ncol)    = 0._r8

 ! Macrophysical process tendency variables

     QQ(:,:ncol)            = 0._r8
     QQw(:,:ncol)           = 0._r8
     QQi(:,:ncol)           = 0._r8
     QQnl(:,:ncol)          = 0._r8
     QQni(:,:ncol)          = 0._r8

 ! Coefficient for computing QQ and related processes

     U_nc(:,:ncol)          = 0._r8
     G_nc(:,:ncol)          = 0._r8

 ! Other

     qmin1(:,:ncol)         = 0._r8
     qmin2(:,:ncol)         = 0._r8
     qmin3(:,:ncol)         = 0._r8

     esat_a(:ncol)          = 0._r8     
     qsat_a(:ncol)          = 0._r8    
     Twb_aw(:,:ncol)        = 0._r8
     qvwb_aw(:,:ncol)       = 0._r8

   ! ---------------- !
   ! Main computation ! 
   ! ---------------- !

   ! ---------------------------------------------------------------------- !
   ! set to zero for levels above
   ! ---------------------------------------------------------------------- !
   if(top_lev .gt. 1)then
   ql0(:top_lev-1,:ncol) = 0._r8
   qi0(:top_lev-1,:ncol) = 0._r8
   nl0(:top_lev-1,:ncol) = 0._r8
   ni0(:top_lev-1,:ncol) = 0._r8
   end if
   
   ! ---------------------------------------------------------------------- !
   ! Check if input non-cumulus pixels satisfie a non-negative constraint.  !
   ! If not, force all water vapor substances to be positive in all layers. !
   ! We should use 'old' cumulus properties for this routine.               !                
   ! ---------------------------------------------------------------------- !

   T(:,:ncol)    =  T0(:,:ncol) 
   qv(:,:ncol)   = qv0(:,:ncol) 
   ql(:,:ncol)   = ql0(:,:ncol) 
   qi(:,:ncol)   = qi0(:,:ncol) 
   nl(:,:ncol)   = nl0(:,:ncol) 
   ni(:,:ncol)   = ni0(:,:ncol) 

   ! inquire cloud index
   qmin1(:,:ncol) = phy_tracer_info(1)%qmin
   do m = 2, pcnst
       if(phy_tracer_info(m)%longname .eq. 'cloud_liquid')then
           ixcldliq = m
           qmin2(:,:ncol) = phy_tracer_info(m)%qmin
       end if
       if(phy_tracer_info(m)%longname .eq. 'cloud_ice')then
           ixcldice = m
           qmin3(:,:ncol) = phy_tracer_info(m)%qmin
       end if
   end do

   call positive_moisture( ncol, dt, qmin1, qmin2, qmin3, dp, & 
                           qv, ql, qi, T, qvten_pwi1, qlten_pwi1, &
                           qiten_pwi1, Tten_pwi1, do_cldice)

   do k = top_lev, pver
   do i = 1, ncol
      if( ql(k,i) .lt. qsmall ) then
          nlten_pwi1(k,i) = -nl(k,i)/dt
          nl(k,i)        = 0._r8
      endif 
      if( qi(k,i) .lt. qsmall ) then
          niten_pwi1(k,i) = -ni(k,i)/dt
          ni(k,i)        = 0._r8
      endif 
   enddo
   enddo

   call subgrid_var(ncol, rpdel, top_lev, p, zm, zi, T, qv, ql,        &
                    lengi, shi, wstarPBL,                       &
                    qtu_shal, umf_shal, cnt_shal, cnb_shal, thlu_shal, &
                    sgm)

   call cldfrac_conden_single_only(ncol, top_lev, p, T, qv, ql, sgm, &
                                   al_st_nc ,conden_st_nc, G_nc,     &
                                   deltaq_sat, deltaq_uns, Q1_sat, Q1_uns, adjust_factor)
 
   do k = top_lev, pver

      ! "False" means that ice will not be taken into account.
      call findsp_vc(qv(k,:ncol),T(k,:ncol),p(k,:ncol), .false., &
           Twb_aw(k,:ncol),qvwb_aw(k,:ncol))

      call qsat_water(T(k,1:ncol), p(k,1:ncol), esat_a(1:ncol), qsat_a(1:ncol))

      !do i = 1, ncol
      !   U_nc(k,i)    =  qv(k,i)/qsat_a(i)
      !enddo

      !call aist_vector(qv(k,:ncol),T(k,:ncol),p(k,:ncol),qi(k,:ncol),landfrac(:ncol),snowh(:ncol),ai_st_nc(k,:ncol),ncol)

      ! use default cloud fraction scheme, for comparison?? LiXH
      !call astG_PDF(U_nc(k,:),p(k,:),qv(k,:),landfrac(:),snowh(:),alst_def(k,:),G_nc_def(k,:),ncol)

      do i = 1, ncol

         cond_gp(k,i) = adjust_factor(k,i)*(conden_st_nc(k,i)-ql(k,i))/dt
 
         QQ(k,i) = cond_gp(k,i)
 
         ! in order to diagnose the output cloud fraction and condensate
         alst_gp(k,i) = al_st_nc(k,i)
         !cond_gp(k,i) = conden_st_nc(k,i)
         delta_q(k,i) = G_nc(k,i)
 
         if(delta_q(k,i).ge.0._r8)then
             N1(k,i) = 1._r8
         else
             N1(k,i) = 0._r8
         endif
 
         ! adjust tendency
         if( QQ(k,i) .ge. 0._r8 ) then
             QQmax    = (qv(k,i) - phy_tracer_info(1)%qmin)/dt ! For ghost cumulus & semi-ghost ice stratus
             QQmax    = max(0._r8,QQmax)
             QQ(k,i)  = min(QQ(k,i),QQmax)
             QQw(k,i) = QQ(k,i)
             QQi(k,i) = 0._r8
         else
             QQmin  = 0._r8
             if( qv(k,i) .lt. qsat_a(i) ) QQmin = min(0._r8,cone*(qv(k,i)-qvwb_aw(k,i))/dt)
             QQ(k,i)  = max(QQ(k,i),QQmin)
             QQw(k,i) = QQ(k,i)
             QQi(k,i) = 0._r8
             QQwmin   = min(0._r8,-cone*ql(k,i)/dt)
             QQimin   = min(0._r8,-cone*qi(k,i)/dt)
             QQw(k,i) = min(0._r8,max(QQw(k,i),QQwmin))
             QQi(k,i) = min(0._r8,max(QQi(k,i),QQimin))
         endif

         if( QQw(k,i) .lt. 0._r8 ) then
             if( ql(k,i) .gt. qsmall ) then
                 QQnl(k,i) = QQw(k,i)*nl(k,i)/ql(k,i)
                 QQnl(k,i) = min(0._r8,cone*max(QQnl(k,i),-nl(k,i)/dt))
             else
                 QQnl(k,i) = 0._r8
             endif  
         endif 

         if( QQi(k,i) .lt. 0._r8 ) then
             if( qi(k,i) .gt. qsmall ) then
                 QQni(k,i) = QQi(k,i)*ni(k,i)/qi(k,i)
                 QQni(k,i) = min(0._r8,cone*max(QQni(k,i),-ni(k,i)/dt))
 
             else
                 QQni(k,i) = 0._r8
             endif  
         endif 

      enddo
   enddo

   do k = top_lev, pver
   do i = 1, ncol
       T_prime0(k,i)  = T(k,i)  + dt*( (latvap*QQw(k,i)+(latvap+latice)*QQi(k,i))/cp )
       qv_prime0(k,i) = qv(k,i) + dt*( 0._r8 - QQw(k,i) - QQi(k,i) )
       ql_prime0(k,i) = ql(k,i) + dt*( 0._r8 + QQw(k,i) )
       qi_prime0(k,i) = qi(k,i) + dt*( 0._r8 + QQi(k,i) )
       nl_prime0(k,i) = max(0._r8,nl(k,i) + dt*( 0._r8 + QQnl(k,i) ))
       ni_prime0(k,i) = max(0._r8,ni(k,i) + dt*( 0._r8 + QQni(k,i) ))
      if( ql_prime0(k,i) .lt. qsmall ) nl_prime0(k,i) = 0._r8
      if( qi_prime0(k,i) .lt. qsmall ) ni_prime0(k,i) = 0._r8
   enddo

   enddo

   ! -------------------------------------------------- !
   ! Perform diagnostic 'positive_moisture' constraint. !
   ! -------------------------------------------------- !
 
   call positive_moisture( ncol, dt, qmin1, qmin2, qmin3, dp,          & 
                           qv_prime0, ql_prime0, qi_prime0, T_prime0,  &
                           qvten_pwi2, qlten_pwi2, qiten_pwi2, Tten_pwi2, do_cldice)

   do k = top_lev, pver
   do i = 1, ncol
      if( ql_prime0(k,i) .lt. qsmall ) then
          nlten_pwi2(k,i) = -nl_prime0(k,i)/dt
          nl_prime0(k,i)   = 0._r8
      endif 
      if( qi_prime0(k,i) .lt. qsmall ) then
          niten_pwi2(k,i) = -ni_prime0(k,i)/dt
          ni_prime0(k,i)   = 0._r8
      endif 
   enddo
   enddo

   do k = top_lev, pver
   !Lixh Test: use updated state to diagnose aist
   call aist_vector(qv_prime0(k,:ncol),T_prime0(k,:ncol),p(k,:ncol),qi(k,:ncol),landfrac(:ncol),snowh(:ncol),ai_st_nc(k,:ncol),ncol)

   ! currently this subroutine is just used to diagnose the cloud fraction
      call instratus_condensate1( ncol              , k                  , ql_prime0(k,:ncol) , qi_prime0(k,:ncol), &
                                 a_cu0(k,:ncol)     , zeros(k,:ncol)     , zeros(k,:ncol)     ,                     & 
                                 zeros(k,:ncol)     , zeros(k,:ncol)     , zeros(k,:ncol)     ,                     &
                                 al_st_nc(k,:ncol)  , ai_st_nc(k,:ncol),                                            &
                                 al_st_star(k,:ncol), ai_st_star(k,:ncol), ql_st_star(k,:ncol), qi_st_star(k,:ncol) )

      a_st_star(k,:ncol)  = max(al_st_star(k,:ncol),ai_st_star(k,:ncol))

   enddo

   qme(top_lev:,:ncol)        = QQw(top_lev:,:ncol)+QQi(top_lev:,:ncol)
 
   s_tendout(top_lev:,:ncol)  = cp*( T_prime0(top_lev:,:ncol)  -  T0(top_lev:,:ncol) )/dt
   qv_tendout(top_lev:,:ncol) =    ( qv_prime0(top_lev:,:ncol) - qv0(top_lev:,:ncol) )/dt
   ql_tendout(top_lev:,:ncol) =    ( ql_prime0(top_lev:,:ncol) - ql0(top_lev:,:ncol) )/dt
   qi_tendout(top_lev:,:ncol) =    ( qi_prime0(top_lev:,:ncol) - qi0(top_lev:,:ncol) )/dt

   nl_tendout(top_lev:,:ncol) =    ( nl_prime0(top_lev:,:ncol) - nl0(top_lev:,:ncol) )/dt
   ni_tendout(top_lev:,:ncol) =    ( ni_prime0(top_lev:,:ncol) - ni0(top_lev:,:ncol) )/dt

   if (.not. do_cldice) then
      do k = top_lev, pver
         do i = 1, ncol

            ! Don't want either qi or ni tendencies, but the code above is somewhat convoluted and
            ! is trying to adjust both (small numbers). Just force it to zero here.
            qi_tendout(k,i) = 0._r8
            ni_tendout(k,i) = 0._r8
          end do
      end do
   end if

   ! ------------------ !
   ! Net cloud fraction !
   ! ------------------ !

   cld(top_lev:,:ncol) = a_st_star(top_lev:,:ncol) + a_cu0(top_lev:,:ncol)

   ! --------------------------------- !
   ! Updated grid-mean state variables !
   ! --------------------------------- !

   T0(top_lev:,:ncol)  = T_prime0(top_lev:,:ncol)
   qv0(top_lev:,:ncol) = qv_prime0(top_lev:,:ncol)
   ql0(top_lev:,:ncol) = ql_prime0(top_lev:,:ncol)
   qi0(top_lev:,:ncol) = qi_prime0(top_lev:,:ncol)
   nl0(top_lev:,:ncol) = nl_prime0(top_lev:,:ncol)
   ni0(top_lev:,:ncol) = ni_prime0(top_lev:,:ncol)

   return
   end subroutine mmacro_pcond_lin


! Lin_PDF scheme:
! Purpose: diagnose and check ql_st(qi_st) using grid-mean ql(qi) and alst(aist)
   subroutine instratus_condensate1( ncol, k, ql0_in, qi0_in,             & 
                                     a_dc_in, ql_dc_in, qi_dc_in,         &
                                     a_sc_in, ql_sc_in, qi_sc_in,         & 
                                     al0_st_nc_in, ai0_st_nc_in,          &
                                     al_st_out, ai_st_out, ql_st_out, qi_st_out)
   implicit none
! io
   integer,  intent(in)  :: ncol                 ! Number of atmospheric columns
   integer,  intent(in)  :: k                    ! Layer index
   real(r8), intent(in)  :: ql0_in(ncol)        ! Grid-mean LWC [kg/kg]
   real(r8), intent(in)  :: qi0_in(ncol)        ! Grid-mean IWC [kg/kg]
   real(r8), intent(in)  :: a_dc_in(ncol)       ! Deep cumulus cloud fraction
   real(r8), intent(in)  :: ql_dc_in(ncol)      ! In-deep cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_dc_in(ncol)      ! In-deep cumulus IWC [kg/kg]
   real(r8), intent(in)  :: a_sc_in(ncol)       ! Shallow cumulus cloud fraction
   real(r8), intent(in)  :: ql_sc_in(ncol)      ! In-shallow cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_sc_in(ncol)      ! In-shallow cumulus IWC [kg/kg]
   real(r8), intent(in)  :: al0_st_nc_in(ncol)
   real(r8), intent(in)  :: ai0_st_nc_in(ncol)

   real(r8), intent(out) :: al_st_out(ncol)     ! Liquid stratus fraction
   real(r8), intent(out) :: ai_st_out(ncol)     ! Ice stratus fraction
   real(r8), intent(out) :: ql_st_out(ncol)     ! In-stratus LWC [kg/kg]
   real(r8), intent(out) :: qi_st_out(ncol)     ! In-stratus IWC [kg/kg]

! local
   integer i                                     ! Column    index
   real(r8) a_dc   
   real(r8) ql_dc  
   real(r8) qi_dc  
   real(r8) a_sc   
   real(r8) ql_sc  
   real(r8) qi_sc  
   real(r8) ql_st
   real(r8) qi_st
   real(r8) al_st_nc
   real(r8) ai_st_nc
   real(r8) al_st
   real(r8) ai_st
 
   do i = 1, ncol

      a_dc  = a_dc_in(i)
      ql_dc = ql_dc_in(i)
      qi_dc = qi_dc_in(i)

      a_sc  = a_sc_in(i)
      ql_sc = ql_sc_in(i)
      qi_sc = qi_sc_in(i)


      ai_st  = (1._r8-a_dc-a_sc)*ai0_st_nc_in(i)
      al_st  = (1._r8-a_dc-a_sc)*al0_st_nc_in(i)

     if( al_st .eq. 0._r8 ) then
         if( ql0_in(i) .gt. 1.e-12_r8) then
!             print*,'Warning, dense liquid cloud!! Rank=',mpi_rank(),'i=',i,'k=',k,'ql=',ql0_in(i)
         end if   
         ql_st = 0._r8
     elseif(ql0_in(i) .le. 1.e-12_r8 )then
         if(al_st .gt. 1.e-4_r8)then
!             print*,'Warning, empty liquid cloud!! Rank=',mpi_rank(),'i=',i,'k=',k,'cld=',al_st
             al_st = 0._r8
         end if
         ql_st = 0._r8
     else
         ql_st = ql0_in(i)/al_st
     endif

     if( ai_st .eq. 0._r8) then
         if( qi0_in(i) .gt. 1.e-12_r8) then
             print*,'Warning, dense ice cloud!! Rank=',mpi_rank(),'i=',i,'k=',k,'qi=',qi0_in(i)
         end if
         qi_st = 0._r8
     elseif(qi0_in(i) .le. 1.e-12_r8)then
         if(ai_st .gt. 1.e-4_r8)then
             print*,'Warning, empty ice cloud!! Rank=',mpi_rank(),'i=',i,'k=',k,'cld=',ai_st
             ai_st = 0._r8
         end if
         qi_st = 0._r8
     else
         qi_st = qi0_in(i)/ai_st
     endif
   al_st_out(i) = al_st
   ai_st_out(i) = ai_st
   ql_st_out(i) = ql_st
   qi_st_out(i) = qi_st

   enddo 

   return
   end subroutine instratus_condensate1



!---------------------Lin Macrop-------------------------
! Purpose: to find saturation equilibrium state using a Newton 
!          iteration method, so that 'qc_st = qcst_crit' is satisfied.
!   subroutine instratus_core( icol, k, p,                             &
!                              T0, qv0, ql0, qi0,                      &
!                              a_dc, ql_dc, qi_dc,                     & 
!                              a_sc, ql_sc, qi_sc, ai_st,              &
!                              qcst_crit, Tmin, Tmax, landfrac, snowh, &
!                              T, qv, ql, qi )
!
!   implicit none
!! io
!   integer,  intent(in)  :: icol       ! Number of atmospheric columns
!   integer,  intent(in)  :: k          ! Layer index
!
!   real(r8), intent(in)  :: p          ! Pressure [Pa]
!   real(r8), intent(in)  :: T0         ! Temperature [K]
!   real(r8), intent(in)  :: qv0        ! Grid-mean water vapor [kg/kg]
!   real(r8), intent(in)  :: ql0        ! Grid-mean LWC [kg/kg]
!   real(r8), intent(in)  :: qi0        ! Grid-mean IWC [kg/kg]
!
!   real(r8), intent(in)  :: a_dc       ! Deep cumulus cloud fraction
!   real(r8), intent(in)  :: ql_dc      ! In-deep cumulus LWC [kg/kg]
!   real(r8), intent(in)  :: qi_dc      ! In-deep cumulus IWC [kg/kg]
!   real(r8), intent(in)  :: a_sc       ! Shallow cumulus cloud fraction
!   real(r8), intent(in)  :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
!   real(r8), intent(in)  :: qi_sc      ! In-shallow cumulus IWC [kg/kg]
!   real(r8), intent(in)  :: ai_st      ! Ice stratus fraction (fixed)
!
!   real(r8), intent(in)  :: Tmin       ! Minimum temperature system can have [K]
!   real(r8), intent(in)  :: Tmax       ! Maximum temperature system can have [K]
!   real(r8), intent(in)  :: qcst_crit  ! Critical in-stratus condensate [kg/kg]
!   real(r8), intent(in)  :: landfrac   ! Land fraction
!   real(r8), intent(in)  :: snowh      ! Snow depth (liquid water equivalent)
!
!   real(r8), intent(out) :: T          ! Temperature [K]
!   real(r8), intent(out) :: qv         ! Grid-mean water vapor [kg/kg]
!   real(r8), intent(out) :: ql         ! Grid-mean LWC [kg/kg]
!   real(r8), intent(out) :: qi         ! Grid-mean IWC [kg/kg]
!! local
!   integer i                           ! Iteration index
!   real(r8) muQ0, muQ
!   real(r8) ql_nc0, qi_nc0, qc_nc0, qc_nc    
!   real(r8) fice0, fice    
!   real(r8) ficeg0, ficeg   
!   real(r8) esat0
!   real(r8) qsat0
!   real(r8) dqcncdt, dastdt, dUdt
!   real(r8) alpha, beta
!   real(r8) U, U_nc
!   real(r8) al_st_nc, G_nc
!   real(r8) al_st
!
!   ! Variables for root-finding algorithm
!
!   integer j                          
!   real(r8)  x1, x2
!   real(r8)  rtsafe
!   real(r8)  df, dx, dxold, f, fh, fl, temp, xh, xl
!   real(r8), parameter :: xacc = 1.e-3_r8
!
!   ql_nc0 = max(0._r8,ql0-a_dc*ql_dc-a_sc*ql_sc)
!   qi_nc0 = max(0._r8,qi0-a_dc*qi_dc-a_sc*qi_sc)
!   qc_nc0 = max(0._r8,ql0+qi0-a_dc*(ql_dc+qi_dc)-a_sc*(ql_sc+qi_sc))
!   fice0  = 0._r8
!   ficeg0 = 0._r8
!   muQ0   = 1._r8
!
!   x1 = Tmin
!   x2 = Tmax
!   call funcd_instratus( x1, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
!                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
!                         qcst_crit, landfrac, snowh,                    &
!                         fl, df, qc_nc, fice, al_st )
!   call funcd_instratus( x2, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
!                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
!                         qcst_crit, landfrac, snowh,                    &
!                         fh, df, qc_nc, fice, al_st )
!    if((fl > 0._r8 .and. fh > 0._r8) .or. (fl < 0._r8 .and. fh < 0._r8)) then
!       call funcd_instratus( T0, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
!                             a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
!                             qcst_crit, landfrac, snowh,                    &
!                             fl, df, qc_nc, fice, al_st )
!       rtsafe = T0 
!       goto 10       
!   endif
!   if( fl == 0._r8) then
!           rtsafe = x1
!           goto 10
!   elseif( fh == 0._r8) then
!           rtsafe = x2
!           goto 10
!   elseif( fl < 0._r8) then
!           xl = x1
!           xh = x2
!   else
!           xh = x1
!           xl = x2
!   end if
!   rtsafe = 0.5_r8*(x1+x2)
!   dxold = abs(x2-x1)
!   dx = dxold
!   call funcd_instratus( rtsafe, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
!                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,     &
!                         qcst_crit, landfrac, snowh,                        &
!                         f, df, qc_nc, fice, al_st )
!   do j = 1, 20
!      if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) > 0._r8 .or. abs(2.0_r8*f) > abs(dxold*df) ) then
!           dxold = dx
!           dx = 0.5_r8*(xh-xl)
!           rtsafe = xl + dx
!           if(xl == rtsafe) goto 10
!      else
!           dxold = dx
!           dx = f/df
!           temp = rtsafe
!           rtsafe = rtsafe - dx
!           if (temp == rtsafe) goto 10
!      end if
!    ! if(abs(dx) < xacc) goto 10
!      call funcd_instratus( rtsafe, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
!                            a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,     &
!                            qcst_crit, landfrac, snowh,                        &
!                            f, df, qc_nc, fice, al_st )
! 
!    ! Sep.21.2010. Sungsu modified to enhance convergence and guarantee 'qlst_min <  qlst < qlst_max'.
!      if( qcst_crit < 0.5_r8 * ( qlst_min + qlst_max ) ) then
!          if( ( qc_nc*(1._r8-fice) .gt.          qlst_min*al_st .and. &
!                qc_nc*(1._r8-fice) .lt. 1.1_r8 * qlst_min*al_st ) ) goto 10
!      else
!          if( ( qc_nc*(1._r8-fice) .gt. 0.9_r8 * qlst_max*al_st .and. &
!                qc_nc*(1._r8-fice) .lt.          qlst_max*al_st ) ) goto 10
!      endif
!      if(f < 0._r8) then
!          xl = rtsafe
!      else
!          xh = rtsafe
!      endif
!
!   enddo
!
!10 continue
!
!   ! ------------------------------------------- !
!   ! Final safety check before sending to output !
!   ! ------------------------------------------- !
!
!   qc_nc = max(0._r8,qc_nc)
!
!   T  = rtsafe
!   ql = qc_nc*(1._r8-fice) + a_dc*ql_dc + a_sc*ql_sc
!   qi = qc_nc*fice + a_dc*qi_dc + a_sc*qi_sc
!   qv = qv0 + ql0 + qi0 - (qc_nc + a_dc*(ql_dc+qi_dc) + a_sc*(ql_sc+qi_sc))
!   qv = max(qv,1.e-12_r8) 
!
!   return
!   end subroutine instratus_core


! Purpose: to find function value and gradient at T
!   subroutine funcd_instratus( T, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0,   &
!                               a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,  &
!                               qcst_crit, landfrac, snowh,                     &
!                               f, fg, qc_nc, fice, al_st ) 
!
!   implicit none
!! io
!   real(r8), intent(in)  :: T          ! Iteration temperature [K]
!
!   real(r8), intent(in)  :: p          ! Pressure [Pa]
!   real(r8), intent(in)  :: T0         ! Initial temperature [K]
!   real(r8), intent(in)  :: qv0        ! Grid-mean water vapor [kg/kg]
!   real(r8), intent(in)  :: ql0        ! Grid-mean LWC [kg/kg]
!   real(r8), intent(in)  :: qi0        ! Grid-mean IWC [kg/kg]
!   real(r8), intent(in)  :: fice0      ! 
!   real(r8), intent(in)  :: muQ0       ! 
!   real(r8), intent(in)  :: qc_nc0     ! 
!
!   real(r8), intent(in)  :: a_dc       ! Deep cumulus cloud fraction
!   real(r8), intent(in)  :: ql_dc      ! In-deep cumulus LWC [kg/kg]
!   real(r8), intent(in)  :: qi_dc      ! In-deep cumulus IWC [kg/kg]
!   real(r8), intent(in)  :: a_sc       ! Shallow cumulus cloud fraction
!   real(r8), intent(in)  :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
!   real(r8), intent(in)  :: qi_sc      ! In-shallow cumulus IWC [kg/kg]
!
!   real(r8), intent(in)  :: ai_st      ! Ice stratus fraction (fixed)
!
!   real(r8), intent(in)  :: qcst_crit  ! Critical in-stratus condensate [kg/kg]
!   real(r8), intent(in)  :: landfrac   ! Land fraction
!   real(r8), intent(in)  :: snowh      ! Snow depth (liquid water equivalent)
!
!   real(r8), intent(out) :: f          ! Value of minimization function at T
!   real(r8), intent(out) :: fg         ! Gradient of minimization function 
!   real(r8), intent(out) :: qc_nc      !
!   real(r8), intent(out) :: al_st      !
!   real(r8), intent(out) :: fice       !
!! local
!   real(r8) es
!   real(r8) qs
!   real(r8) dqsdT
!   real(r8) dqcncdt
!   real(r8) alpha
!   real(r8) beta
!   real(r8) U
!   real(r8) U_nc
!   real(r8) al_st_nc
!   real(r8) G_nc
!   real(r8) dUdt
!   real(r8) dalstdt
!   real(r8) qv
!
!   call qsat_water(T, p, es, qs, dqsdt=dqsdT)
!
!   fice    = fice0 
!   qc_nc   = (cp/latvap)*(T-T0)+muQ0*qc_nc0       
!   dqcncdt = (cp/latvap) 
!   qv      = (qv0 + ql0 + qi0 - (qc_nc + a_dc*(ql_dc+qi_dc) + a_sc*(ql_sc+qi_sc)))
!   alpha   = (1._r8/qs)
!   beta    = (qv/qs**2._r8)*dqsdT 
!
!   U      =  (qv/qs)
!   U_nc   =   U
!   if( CAMstfrac ) then
!       call astG_RHU_single(U_nc,p,qv,landfrac,snowh,al_st_nc,G_nc)
!   else
!       call astG_PDF_single(U_nc,p,qv,landfrac,snowh,al_st_nc,G_nc)
!   endif
!   al_st   =  (1._r8-a_dc-a_sc)*al_st_nc 
!   dUdt    = -(alpha*dqcncdt+beta)
!   dalstdt =  (1._r8/G_nc)*dUdt
!   if( U_nc .eq. 1._r8 ) dalstdt = 0._r8
!
!   f  = qc_nc   - qcst_crit*al_st
!   fg = dqcncdt - qcst_crit*dalstdt
!
!   return
!   end subroutine funcd_instratus


! Purpose: to force grid-mean RH = 1 when RH > 1
!          This is condensation process similar to instratus_condensate.
!          During condensation, we assume 'fice' is maintained in this
!          verison for MG not for RK. 
!   subroutine gridmean_RH( icol, k, p, T, qv, ql, qi,               &
!                           a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc,  &
!                           landfrac, snowh )
!
!   implicit none
!
!   integer,  intent(in)    :: icol       ! Number of atmospheric columns
!   integer,  intent(in)    :: k          ! Layer index
!
!   real(r8), intent(in)    :: p          ! Pressure [Pa]
!   real(r8), intent(inout) :: T          ! Temperature [K]
!   real(r8), intent(inout) :: qv         ! Grid-mean water vapor [kg/kg]
!   real(r8), intent(inout) :: ql         ! Grid-mean LWC [kg/kg]
!   real(r8), intent(inout) :: qi         ! Grid-mean IWC [kg/kg]
!
!   real(r8), intent(in)    :: a_dc       ! Deep cumulus cloud fraction
!   real(r8), intent(in)    :: ql_dc      ! In-deep cumulus LWC [kg/kg]
!   real(r8), intent(in)    :: qi_dc      ! In-deep cumulus IWC [kg/kg]
!   real(r8), intent(in)    :: a_sc       ! Shallow cumulus cloud fraction
!   real(r8), intent(in)    :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
!   real(r8), intent(in)    :: qi_sc      ! In-shallow cumulus IWC [kg/kg]
!
!   real(r8), intent(in)    :: landfrac   ! Land fraction
!   real(r8), intent(in)    :: snowh      ! Snow depth (liquid water equivalent)
!
!   ! Local variables
!
!   integer m                             ! Iteration index
!
!   real(r8)  ql_nc0, qi_nc0, qc_nc0
!   real(r8)  Tscale
!   real(r8)  Tc, qt, qc, dqcdt, qc_nc    
!   real(r8)  es, qs, dqsdT
!   real(r8)  al_st_nc, G_nc
!   real(r8)  f, fg
!   real(r8), parameter :: xacc = 1.e-3_r8
!
!   ql_nc0 = max(0._r8,ql-a_dc*ql_dc-a_sc*ql_sc)
!   qi_nc0 = max(0._r8,qi-a_dc*qi_dc-a_sc*qi_sc)
!   qc_nc0 = max(0._r8,ql+qi-a_dc*(ql_dc+qi_dc)-a_sc*(ql_sc+qi_sc))
!   Tc    = T - (latvap/cp)*ql
!   qt    = qv + ql
!
!   do m = 1, 20
!      call qsat_water(T, p, es, qs, dqsdt=dqsdT)
!      Tscale = latvap/cp
!      qc     = (T-Tc)/Tscale
!      dqcdt  = 1._r8/Tscale
!      f      = qs + qc - qt 
!      fg     = dqsdT + dqcdt
!      fg     = sign(1._r8,fg)*max(1.e-10_r8,abs(fg))
!    ! Sungsu modified convergence criteria to speed up convergence and guarantee RH <= 1.
!      if( qc .ge. 0._r8 .and. ( qt - qc ) .ge. 0.999_r8*qs .and. ( qt - qc ) .le. 1._r8*qs ) then
!          goto 10
!      endif
!      T = T - f/fg
!   enddo
! ! write(iulog,*) 'Convergence in gridmean_RH is not reached. RH = ', ( qt - qc ) / qs
!10 continue
!
!   call qsat_water(T, p, es, qs)
! ! Sungsu modified 'qv = qs' in consistent with the modified convergence criteria above.
!   qv = min(qt,qs) ! Modified
!   ql = qt - qv
!   T  = Tc + (latvap/cp)*ql
!
!   return
!   end subroutine gridmean_RH


! Purpose: If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,
!          force them to be larger than minimum value by (1) condensating water vapor
!          into liquid or ice, and (2) by transporting water vapor from the very lower
!          layer. '2._r8' is multiplied to the minimum values for safety.
!          Update final state variables and tendencies associated with this correction.
!          If any condensation happens, update (s,t) too.
!          Note that (qv,ql,qi,t,s) are final state variables after applying corresponding
!          input tendencies.
!          Be careful the order of k : '1': top layer, 'pver' : near-surface layer
   subroutine positive_moisture( ncol, dt, qvmin, qlmin, qimin, dp, &
                                 qv,   ql, qi,    t,     qvten, &
                                 qlten,    qiten, tten,  do_cldice)
   implicit none
   integer,  intent(in)     :: ncol
   real(r8), intent(in)     :: dt
   real(r8), intent(in)     :: dp(pver, ncol), qvmin(pver, ncol), qlmin(pver, ncol), qimin(pver, ncol)
   real(r8), intent(inout)  :: qv(pver, ncol), ql(pver, ncol), qi(pver, ncol), t(pver, ncol)
   real(r8), intent(out)    :: qvten(pver, ncol), qlten(pver, ncol), qiten(pver, ncol), tten(pver, ncol)
   logical, intent(in)      :: do_cldice
   integer   i, k
   real(r8)  dql, dqi, dqv, sum, aa, dum 

   tten(:pver,:ncol)  = 0._r8
   qvten(:pver,:ncol) = 0._r8
   qlten(:pver,:ncol) = 0._r8
   qiten(:pver,:ncol) = 0._r8

   do i = 1, ncol
      do k = top_lev, pver
         if( qv(k,i) .lt. qvmin(k,i) .or. ql(k,i) .lt. qlmin(k,i) .or. qi(k,i) .lt. qimin(k,i) ) then
             goto 10
         endif
      enddo
      goto 11
   10 continue
      do k = top_lev, pver    ! From the top to the 1st (lowest) layer from the surface
         dql = max(0._r8,1._r8*qlmin(k,i)-ql(k,i))

         if (do_cldice) then
         dqi = max(0._r8,1._r8*qimin(k,i)-qi(k,i))
         else
           dqi = 0._r8
         end if

         qlten(k,i) = qlten(k,i) +  dql/dt
         qiten(k,i) = qiten(k,i) +  dqi/dt
         qvten(k,i) = qvten(k,i) - (dql+dqi)/dt
         tten(k,i)  = tten(k,i)  + (latvap/cp)*(dql/dt) + ((latvap+latice)/cp)*(dqi/dt)
         ql(k,i)    = ql(k,i) + dql
         qi(k,i)    = qi(k,i) + dqi
         qv(k,i)    = qv(k,i) - dql - dqi
         t(k,i)     = t(k,i)  + (latvap * dql + (latvap+latice) * dqi)/cp
         dqv        = max(0._r8,1._r8*qvmin(k,i)-qv(k,i))
         qvten(k,i) = qvten(k,i) + dqv/dt
         qv(k,i)    = qv(k,i)    + dqv
         if( k .ne. pver ) then 
             qv(k+1,i)    = qv(k+1,i)    - dqv*dp(k,i)/dp(k+1,i)
             qvten(k+1,i) = qvten(k+1,i) - dqv*dp(k,i)/dp(k+1,i)/dt
         endif
         qv(k,i) = max(qv(k,i),qvmin(k,i))
         ql(k,i) = max(ql(k,i),qlmin(k,i))
         qi(k,i) = max(qi(k,i),qimin(k,i))
      end do
      ! Extra moisture used to satisfy 'qv(pver,i)=qvmin' is proportionally 
      ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
      ! preserves column moisture. 
      if( dqv .gt. 1.e-20_r8 ) then
          sum = 0._r8
          do k = top_lev, pver
             if( qv(k,i) .gt. 2._r8*qvmin(k,i) ) sum = sum + qv(k,i)*dp(k,i)
          enddo
          aa = dqv*dp(pver,i)/max(1.e-20_r8,sum)
          if( aa .lt. 0.5_r8 ) then
              do k = top_lev, pver
                 if( qv(k,i) .gt. 2._r8*qvmin(k,i) ) then
                     dum        = aa*qv(k,i)
                     qv(k,i)    = qv(k,i) - dum
                     qvten(k,i) = qvten(k,i) - dum/dt
                 endif
              enddo 
          else 
              print*, 'Full positive_moisture is impossible in Park Macro'
              print*, 'rank = ',mpi_rank()
          endif
      endif 
11 continue
   enddo
   return

   end subroutine positive_moisture


!---------------------Lin Macrop-------------------------
! Purpose: Compute 'stratus fraction(a)' and Gs=(dU/da) from the
!          analytical formulation of triangular PDF.                 
!          Here, 'dV' is the ratio of 'half-width of PDF / qs(p,T)', 
!          so using constant 'dV' assume that width is proportional  
!          to the saturation specific humidity.                      
!             dV ~ 0.1.                                              
!             cldrh : RH of in-stratus( = 1 if no supersaturation)   
!          Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is 
!          G is discontinuous across U = 1.  In fact, it does not    
!          matter whether Ga = 1.e10 or 0 at a = 1: I derived that   
!          they will produce the same results.                       
!   subroutine astG_PDF_single( U, p, qv, landfrac, snowh, a, Ga, orhmin )
!
!   implicit none
!
!   real(r8), intent(in)  :: U                     ! Relative humidity
!   real(r8), intent(in)  :: p                     ! Pressure [Pa]
!   real(r8), intent(in)  :: qv                    ! Grid-mean water vapor specific humidity [kg/kg]
!   real(r8), intent(in)  :: landfrac              ! Land fraction
!   real(r8), intent(in)  :: snowh                 ! Snow depth (liquid water equivalent)
!
!   real(r8), intent(out) :: a                     ! Stratus fraction
!   real(r8), intent(out) :: Ga                    ! dU/da
!   real(r8), optional, intent(out) :: orhmin      ! Critical RH
!
!   ! Local variables
!   integer :: i                                   ! Loop indexes
!   real(r8) dV                                    ! Width of triangular PDF
!   real(r8) cldrh                                 ! RH of stratus cloud
!   real(r8) rhmin                                 ! Critical RH
!   real(r8) rhwght
!                            
!   ! Statement functions
!   logical land
!   land = nint(landfrac) == 1
!
!   cldrh  = 1.0_r8
!
!   if( p .ge. premib ) then
!
!       if( land .and. (snowh.le.0.000001_r8) ) then
!           rhmin = rhminl - rhminl_adj_land
!       else
!           rhmin = rhminl
!       endif
!
!       dV = cldrh - rhmin
!
!       if( U .ge. 1._r8 ) then
!           a  = 1._r8
!           Ga = 1.e10_r8
!       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
!           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
!           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
!       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
!           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
!                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
!           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
!       elseif( U .le. (cldrh-dV) ) then
!           a  = 0._r8
!           Ga = 1.e10_r8
!       endif
!
!       if( freeze_dry ) then
!           a  = a *max(0.15_r8,min(1.0_r8,qv/0.0030_r8))
!           Ga = Ga/max(0.15_r8,min(1.0_r8,qv/0.0030_r8)) 
!       endif
!
!   elseif( p .lt. premit ) then
!
!       rhmin = rhminh
!       dV    = cldrh - rhmin
!
!       if( U .ge. 1._r8 ) then
!           a  = 1._r8
!           Ga = 1.e10_r8
!       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
!           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
!           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
!       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
!           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
!                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
!           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
!       elseif( U .le. (cldrh-dV) ) then
!           a  = 0._r8
!           Ga = 1.e10_r8
!       endif
!   else
!
!       rhwght = (premib-(max(p,premit)))/(premib-premit)
!
!     ! if( land .and. (snowh.le.0.000001_r8) ) then
!     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
!     ! else
!           rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
!     ! endif
!
!       dV    = cldrh - rhmin
!
!       if( U .ge. 1._r8 ) then
!           a  = 1._r8
!           Ga = 1.e10_r8
!       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
!           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
!           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
!       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
!           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
!                         (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
!           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
!       elseif( U .le. (cldrh-dV) ) then
!           a  = 0._r8
!           Ga = 1.e10_r8
!       endif
!
!   endif
!
!   if (present(orhmin)) orhmin = rhmin
!
!   return
!   end subroutine astG_PDF_single
!---------------------Lin Macrop-------------------------


! Purpose: Compute 'stratus fraction(a)' and Gs=(dU/da) from the
!          analytical formulation of triangular PDF.
!          Here, 'dV' is the ratio of 'half-width of PDF / qs(p,T)',
!          so using constant 'dV' assume that width is proportional
!          to the saturation specific humidity.
!             dV ~ 0.1.
!             cldrh : RH of in-stratus( = 1 if no supersaturation) 
!          Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is
!          G is discontinuous across U = 1.  In fact, it does not
!          matter whether Ga = 1.e10 or 0 at a = 1: I derived that
!          they will produce the same results.
   subroutine astG_PDF( U_in, p_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, ncol )

   implicit none
   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: U_in(ncol)            ! Relative humidity
   real(r8), intent(in)  :: p_in(ncol)            ! Pressure [Pa]
   real(r8), intent(in)  :: qv_in(ncol)           ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(in)  :: landfrac_in(ncol)     ! Land fraction
   real(r8), intent(in)  :: snowh_in(ncol)        ! Snow depth (liquid water equivalent)
   real(r8), intent(out) :: a_out(ncol)           ! Stratus fraction
   real(r8), intent(out) :: Ga_out(ncol)          ! dU/da

   real(r8)              :: U                     ! Relative humidity
   real(r8)              :: p                     ! Pressure [Pa]
   real(r8)              :: qv                    ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8)              :: landfrac              ! Land fraction
   real(r8)              :: snowh                 ! Snow depth (liquid water equivalent)

   real(r8)              :: a                     ! Stratus fraction
   real(r8)              :: Ga                    ! dU/da

   ! Local variables
   integer :: i                                   ! Loop indexes
   real(r8) dV                                    ! Width of triangular PDF
   real(r8) cldrh                                 ! RH of stratus cloud
   real(r8) rhmin                                 ! Critical RH
   real(r8) rhwght
                            
   ! Statement functions
   logical land(ncol)

   do i = 1, ncol
       land(i) = nint(landfrac_in(i)) == 1
   end do

   cldrh  = 1.0_r8

   a_out(:)  = 0._r8
   Ga_out(:) = 0._r8

   do i = 1, ncol

   U        = U_in(i)      
   p        = p_in(i)        
   qv       = qv_in(i)       
   landfrac = landfrac_in(i) 
   snowh    = snowh_in(i)    

   if( p .ge. premib ) then

       if( land(i) .and. (snowh.le.0.000001_r8) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif

       dV = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

       if( freeze_dry ) then
           a  = a *max(0.15_r8,min(1.0_r8,qv/0.0030_r8))
           Ga = Ga/max(0.15_r8,min(1.0_r8,qv/0.0030_r8)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       dV    = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land(i) .and. (snowh.le.0.000001_r8) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
     ! endif

       dV    = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                         (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

   endif

   a_out(i)  = a
   Ga_out(i) = Ga 

   enddo

   return
   end subroutine astG_PDF


!---------------------Lin Macrop-------------------------
! Purpose: Compute 'stratus fraction(a)' and Gs=(dU/da) from the
!          Below is valid only for CAMUW at 1.9x2.5 fv dynamics core   
!          For the other cases, I should re-define 'rhminl,rhminh' & 
!          'premib,premit'.                                          
!          Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is 
!          G is discontinuous across U = 1.                          
!   subroutine astG_RHU_single( U, p, qv, landfrac, snowh, a, Ga, orhmin )
!
!   implicit none
!
!   real(r8), intent(in)  :: U               ! Relative humidity
!   real(r8), intent(in)  :: p               ! Pressure [Pa]
!   real(r8), intent(in)  :: qv              ! Grid-mean water vapor specific humidity [kg/kg]
!   real(r8), intent(in)  :: landfrac        ! Land fraction
!   real(r8), intent(in)  :: snowh           ! Snow depth (liquid water equivalent)
!
!   real(r8), intent(out) :: a               ! Stratus fraction
!   real(r8), intent(out) :: Ga              ! dU/da
!   real(r8), optional, intent(out) :: orhmin ! Critical RH
!
!   ! Local variables
!   real(r8) rhmin                                 ! Critical RH
!   real(r8) rhdif                                 ! Factor for stratus fraction
!   real(r8) rhwght
!
!   ! Statement functions
!   logical land
!   land = nint(landfrac) == 1
!
!   ! ---------------- !
!   ! Main computation !
!   ! ---------------- !
!
!   if( p .ge. premib ) then
!
!       if( land .and. (snowh.le.0.000001_r8) ) then
!           rhmin = rhminl - rhminl_adj_land
!       else
!           rhmin = rhminl
!       endif
!       rhdif = (U-rhmin)/(1.0_r8-rhmin)
!       a  = min(1._r8,(max(rhdif,0.0_r8))**2) 
!       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
!            Ga = 1.e20_r8
!       else          
!            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
!       endif
!       if( freeze_dry ) then
!           a  = a*max(0.15_r8,min(1.0_r8,qv/0.0030_r8))
!           Ga = Ga/max(0.15_r8,min(1.0_r8,qv/0.0030_r8)) 
!       endif
!
!   elseif( p .lt. premit ) then
!
!       rhmin = rhminh
!       rhdif = (U-rhmin)/(1.0_r8-rhmin)
!       a  = min(1._r8,(max(rhdif,0._r8))**2)
!       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
!            Ga = 1.e20_r8
!       else          
!            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
!       endif
!
!   else
!
!       rhwght = (premib-(max(p,premit)))/(premib-premit)
!
!     ! if( land .and. (snowh.le.0.000001_r8) ) then
!     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
!     ! else
!           rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
!     ! endif
!
!       rhdif = (U-rhmin)/(1.0_r8-rhmin)
!       a  = min(1._r8,(max(rhdif,0._r8))**2)
!       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
!            Ga = 1.e10_r8
!       else          
!            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
!       endif
!
!   endif
!
!   if (present(orhmin)) orhmin = rhmin
!
!   return
!   end subroutine astG_RHU_single


! Purpose: Compute 'stratus fraction(a)' and Gs=(dU/da) from the
!          CAM35 cloud fraction formula.
!          Below is valid only for CAMUW at 1.9x2.5 fv dynamics core
!          For the other cases, I should re-define 'rhminl,rhminh' &
!          'premib,premit'.
!          Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is
!          G is discontinuous across U = 1.
!   subroutine astG_RHU( U_in, p_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, ncol )
!
!   implicit none
!
!   integer,  intent(in)  :: ncol
!   real(r8), intent(in)  :: U_in(ncol)            ! Relative humidity
!   real(r8), intent(in)  :: p_in(ncol)            ! Pressure [Pa]
!   real(r8), intent(in)  :: qv_in(ncol)           ! Grid-mean water vapor specific humidity [kg/kg]
!   real(r8), intent(in)  :: landfrac_in(ncol)     ! Land fraction
!   real(r8), intent(in)  :: snowh_in(ncol)        ! Snow depth (liquid water equivalent)
!   real(r8), intent(out) :: a_out(ncol)           ! Stratus fraction
!   real(r8), intent(out) :: Ga_out(ncol)          ! dU/da
!
!   real(r8)              :: U                     ! Relative humidity
!   real(r8)              :: p                     ! Pressure [Pa]
!   real(r8)              :: qv                    ! Grid-mean water vapor specific humidity [kg/kg]
!   real(r8)              :: landfrac              ! Land fraction
!   real(r8)              :: snowh                 ! Snow depth (liquid water equivalent)
!   real(r8)              :: a                     ! Stratus fraction
!   real(r8)              :: Ga                    ! dU/da
!
!   ! Local variables
!   integer  i
!   real(r8) rhmin                                 ! Critical RH
!   real(r8) rhdif                                 ! Factor for stratus fraction
!   real(r8) rhwght
!
!   ! Statement functions
!   logical land(ncol)
!
!   do i = 1, ncol
!       land(i) = nint(landfrac_in(i)) == 1
!   end do
!
!
!   a_out(:) = 0._r8
!   Ga_out(:) = 0._r8
!
!   do i = 1, ncol
!
!   U        = U_in(i)      
!   p        = p_in(i)        
!   qv       = qv_in(i)       
!   landfrac = landfrac_in(i) 
!   snowh    = snowh_in(i)    
!
!   if( p .ge. premib ) then
!
!       if( land(i) .and. (snowh.le.0.000001_r8) ) then
!           rhmin = rhminl - rhminl_adj_land
!       else
!           rhmin = rhminl
!       endif
!       rhdif = (U-rhmin)/(1.0_r8-rhmin)
!       a  = min(1._r8,(max(rhdif,0.0_r8))**2) 
!       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
!            Ga = 1.e20_r8
!       else          
!            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
!       endif
!       if( freeze_dry ) then
!           a  = a*max(0.15_r8,min(1.0_r8,qv/0.0030_r8))
!           Ga = Ga/max(0.15_r8,min(1.0_r8,qv/0.0030_r8)) 
!       endif
!
!   elseif( p .lt. premit ) then
!
!       rhmin = rhminh
!       rhdif = (U-rhmin)/(1.0_r8-rhmin)
!       a  = min(1._r8,(max(rhdif,0._r8))**2)
!       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
!            Ga = 1.e20_r8
!       else          
!            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
!       endif
!
!   else
!
!       rhwght = (premib-(max(p,premit)))/(premib-premit)
!
!     ! if( land(i) .and. (snowh.le.0.000001_r8) ) then
!     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
!     ! else
!           rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
!     ! endif
!
!       rhdif = (U-rhmin)/(1.0_r8-rhmin)
!       a  = min(1._r8,(max(rhdif,0._r8))**2)
!       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
!            Ga = 1.e10_r8
!       else          
!            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
!       endif
!
!   endif
!
!   a_out(i)  = a
!   Ga_out(i) = Ga 
!
!   enddo
!
!   return
!   end subroutine astG_RHU
!---------------------Lin Macrop-------------------------


! Purpose: Compute non-physical ice stratus fraction
   subroutine aist_single( qv, T, p, qi, landfrac, snowh, aist )

   implicit none
  
   real(r8), intent(in)  :: qv              ! Grid-mean water vapor[kg/kg]
   real(r8), intent(in)  :: T               ! Temperature
   real(r8), intent(in)  :: p               ! Pressure [Pa]
   real(r8), intent(in)  :: qi              ! Grid-mean ice water content [kg/kg]
   real(r8), intent(in)  :: landfrac        ! Land fraction
   real(r8), intent(in)  :: snowh           ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: aist            ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   ! Local variables
   real(r8) rhmin                           ! Critical RH
   real(r8) rhwght

   real(r8) a,b,c,as,bs,cs                  ! Fit parameters
   real(r8) Kc                              ! Constant for ice cloud calc (wood & field)
   real(r8) ttmp                            ! Limited temperature
   real(r8) icicval                         ! Empirical IWC value [ kg/kg ]
   real(r8) rho                             ! Local air density
   real(r8) esl                             ! Liq sat vapor pressure
   real(r8) esi                             ! Ice sat vapor pressure
   real(r8) ncf,phi                         ! Wilson and Ballard parameters
   real(r8) es, qs

   real(r8) rhi                             ! grid box averaged relative humidity over ice
   real(r8) minice                          ! minimum grid box avg ice for having a 'cloud'
   real(r8) mincld                          ! minimum ice cloud fraction threshold
   real(r8) icimr                           ! in cloud ice mixing ratio
 ! real(r8) qist_min                        ! minimum in cloud ice mixing ratio
 ! real(r8) qist_max                        ! maximum in cloud ice mixing ratio                
   real(r8) rhdif                           ! working variable for slingo scheme


   ! Statement functions
   logical land
   land = nint(landfrac) == 1

   ! Wang and Sassen IWC paramters ( Option.1 )
     a = 26.87_r8
     b = 0.569_r8
     c = 0.002892_r8
   ! Schiller parameters ( Option.2 )
     as = -68.4202_r8
     bs = 0.983917_r8
     cs = 2.81795_r8
   ! Wood and Field parameters ( Option.3 )
     Kc = 75._r8
   ! Wilson & Ballard closure ( Option.4. smaller = more ice clouds)
   ! Slingo modified (option 5)
     minice = 1.e-12_r8
     mincld = 1.e-4_r8
   ! qist_min = 1.e-7_r8 
   ! qist_max = 5.e-3_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

     call qsat_water(T, p, es, qs)
     esl = svp_water(T)
     esi = svp_ice(T)
          
     if( iceopt.lt.3 ) then
         if( iceopt.eq.1 ) then
             ttmp = max(195._r8,min(T,253._r8)) - 273.16_r8
             icicval = a + b * ttmp + c * ttmp**2._r8
             rho = p/(rdry*T)
             icicval = icicval * 1.e-6_r8 / rho 
         else
             ttmp = max(190._r8,min(T,273.16_r8))
             icicval = 10._r8 **(as * bs**ttmp + cs)
             icicval = icicval * 1.e-6_r8 * 18._r8 / 28.97_r8
         endif
         aist =  max(0._r8,min(qi/icicval,1._r8)) 
     elseif( iceopt.eq.3 ) then
         aist = 1._r8 - exp(-Kc*qi/(qs*(esi/esl)))
         aist = max(0._r8,min(aist,1._r8))
     elseif( iceopt.eq.4) then
         if( p .ge. premib ) then
             if( land .and. (snowh.le.0.000001_r8) ) then
                 rhmin = rhminl - rhminl_adj_land
             else
                 rhmin = rhminl
             endif
         elseif( p .lt. premit ) then
             rhmin = rhminh
         else
             rhwght = (premib-(max(p,premit)))/(premib-premit)
           ! if( land .and. (snowh.le.0.000001_r8) ) then
           !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
           ! else
                 rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
           ! endif
         endif
         ncf = qi/((1._r8 - icecrit)*qs)
         if( ncf.le.0._r8 ) then 
             aist = 0._r8
         elseif( ncf.gt.0._r8 .and. ncf.le.1._r8/6._r8 ) then 
             aist = 0.5_r8*(6._r8 * ncf)**(2._r8/3._r8)
         elseif( ncf.gt.1._r8/6._r8 .and. ncf.lt.1._r8 ) then
             phi = (acos(3._r8*(1._r8-ncf)/2._r8**(3._r8/2._r8))+4._r8*3.1415927_r8)/3._r8
             aist = (1._r8 - 4._r8 * cos(phi) * cos(phi))
         else
             aist = 1._r8
         endif
             aist = max(0._r8,min(aist,1._r8))
     elseif (iceopt.eq.5) then 
! set rh ice cloud fraction
             rhi= (qv+qi)/qs * (esl/esi)
             rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
             aist = min(1.0_r8, max(rhdif,0._r8)**2)

! limiter to remove empty cloud and ice with no cloud
! and set icecld fraction to mincld if ice exists

             if (qi.lt.minice) then
                aist=0._r8
             else
                aist=max(mincld,aist)
             endif

! enforce limits on icimr
             if (qi.ge.minice) then
                icimr=qi/aist

!minimum
                if (icimr.lt.qist_min) then
                   aist = max(0._r8,min(1._r8,qi/qist_min))
                endif
!maximum
                if (icimr.gt.qist_max) then
                   aist = max(0._r8,min(1._r8,qi/qist_max))
                endif

             endif
     endif 

   ! 0.999_r8 is added to prevent infinite 'ql_st' at the end of instratus_condensate
   ! computed after updating 'qi_st'.  

     aist = max(0._r8,min(aist,0.999_r8))

   return
   end subroutine aist_single


! Purpose: Compute non-physical ice stratus fraction
   subroutine aist_vector( qv_in, T_in, p_in, qi_in, landfrac_in, snowh_in, aist_out, ncol )

   implicit none
  
   integer,  intent(in)  :: ncol 
   real(r8), intent(in)  :: qv_in(ncol)       ! Grid-mean water vapor[kg/kg]
   real(r8), intent(in)  :: T_in(ncol)        ! Temperature
   real(r8), intent(in)  :: p_in(ncol)        ! Pressure [Pa]
   real(r8), intent(in)  :: qi_in(ncol)       ! Grid-mean ice water content [kg/kg]
   real(r8), intent(in)  :: landfrac_in(ncol) ! Land fraction
   real(r8), intent(in)  :: snowh_in(ncol)    ! Snow depth (liquid water equivalent)
   real(r8), intent(out) :: aist_out(ncol)    ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   ! Local variables

   real(r8) qv                              ! Grid-mean water vapor[kg/kg]
   real(r8) T                               ! Temperature
   real(r8) p                               ! Pressure [Pa]
   real(r8) qi                              ! Grid-mean ice water content [kg/kg]
   real(r8) landfrac                        ! Land fraction
   real(r8) snowh                           ! Snow depth (liquid water equivalent)
   real(r8) aist                            ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   real(r8) rhmin                           ! Critical RH
   real(r8) rhwght

   real(r8) a,b,c,as,bs,cs                  ! Fit parameters
   real(r8) Kc                              ! Constant for ice cloud calc (wood & field)
   real(r8) ttmp                            ! Limited temperature
   real(r8) icicval                         ! Empirical IWC value [ kg/kg ]
   real(r8) rho                             ! Local air density
   real(r8) esl                             ! Liq sat vapor pressure
   real(r8) esi                             ! Ice sat vapor pressure
   real(r8) ncf,phi                         ! Wilson and Ballard parameters
   real(r8) qs
   real(r8) esat_in(ncol)
   real(r8) qsat_in(ncol)

   real(r8) rhi                             ! grid box averaged relative humidity over ice
   real(r8) minice                          ! minimum grid box avg ice for having a 'cloud'
   real(r8) mincld                          ! minimum ice cloud fraction threshold
   real(r8) icimr                           ! in cloud ice mixing ratio
 ! real(r8) qist_min                        ! minimum in cloud ice mixing ratio
 ! real(r8) qist_max                        ! maximum in cloud ice mixing ratio                
   real(r8) rhdif                           ! working variable for slingo scheme

   integer i


   ! Statement functions
   logical land(ncol)

   do i = 1, ncol
       land(i) = nint(landfrac_in(i)) == 1
   end do

   ! Wang and Sassen IWC paramters ( Option.1 )
     a = 26.87_r8
     b = 0.569_r8
     c = 0.002892_r8
   ! Schiller parameters ( Option.2 )
     as = -68.4202_r8
     bs = 0.983917_r8
     cs = 2.81795_r8
   ! Wood and Field parameters ( Option.3 )
     Kc = 75._r8
   ! Wilson & Ballard closure ( Option.4. smaller = more ice clouds)
   ! Slingo modified (option 5)
     minice = 1.e-12_r8
     mincld = 1.e-4_r8
   ! qist_min = 1.e-7_r8
   ! qist_max = 5.e-3_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

     aist_out(:) = 0._r8
     esat_in(:)  = 0._r8
     qsat_in(:)  = 0._r8

     call qsat_water(T_in(1:ncol), p_in(1:ncol), esat_in(1:ncol), qsat_in(1:ncol))

     do i = 1, ncol

     landfrac = landfrac_in(i)     
     snowh = snowh_in(i)   
     T = T_in(i)
     qv = qv_in(i)
     p = p_in(i)
     qi = qi_in(i)
     qs = qsat_in(i)
     esl = svp_water(T)
     esi = svp_ice(T)
          
     if( iceopt.lt.3 ) then
         if( iceopt.eq.1 ) then
             ttmp = max(195._r8,min(T,253._r8)) - 273.16_r8
             icicval = a + b * ttmp + c * ttmp**2._r8
             rho = p/(rdry*T)
             icicval = icicval * 1.e-6_r8 / rho 
         else
             ttmp = max(190._r8,min(T,273.16_r8))
             icicval = 10._r8 **(as * bs**ttmp + cs)
             icicval = icicval * 1.e-6_r8 * 18._r8 / 28.97_r8
         endif
         aist =  max(0._r8,min(qi/icicval,1._r8)) 
     elseif( iceopt.eq.3 ) then
         aist = 1._r8 - exp(-Kc*qi/(qs*(esi/esl)))
         aist = max(0._r8,min(aist,1._r8))
     elseif( iceopt.eq.4) then
         if( p .ge. premib ) then
             if( land(i) .and. (snowh.le.0.000001_r8) ) then
                 rhmin = rhminl - rhminl_adj_land
             else
                 rhmin = rhminl
             endif
         elseif( p .lt. premit ) then
             rhmin = rhminh
         else
             rhwght = (premib-(max(p,premit)))/(premib-premit)
           ! if( land(i) .and. (snowh.le.0.000001_r8) ) then
           !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
           ! else
                 rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
           ! endif
         endif
         ncf = qi/((1._r8 - icecrit)*qs)
         if( ncf.le.0._r8 ) then 
             aist = 0._r8
         elseif( ncf.gt.0._r8 .and. ncf.le.1._r8/6._r8 ) then 
             aist = 0.5_r8*(6._r8 * ncf)**(2._r8/3._r8)
         elseif( ncf.gt.1._r8/6._r8 .and. ncf.lt.1._r8 ) then
             phi = (acos(3._r8*(1._r8-ncf)/2._r8**(3._r8/2._r8))+4._r8*3.1415927_r8)/3._r8
             aist = (1._r8 - 4._r8 * cos(phi) * cos(phi))
         else
             aist = 1._r8
         endif
             aist = max(0._r8,min(aist,1._r8))
     elseif (iceopt.eq.5) then 
! set rh ice cloud fraction
             rhi= (qv+qi)/qs * (esl/esi)
             rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
             aist = min(1.0_r8, max(rhdif,0._r8)**2)

! limiter to remove empty cloud and ice with no cloud
! and set icecld fraction to mincld if ice exists

             if (qi.lt.minice) then
                aist=0._r8
             else
                aist=max(mincld,aist)
             endif

! enforce limits on icimr
             if (qi.ge.minice) then
                icimr=qi/aist

!minimum
                if (icimr.lt.qist_min) then
                   aist = max(0._r8,min(1._r8,qi/qist_min))
                endif
!maximum
                if (icimr.gt.qist_max) then
                   aist = max(0._r8,min(1._r8,qi/qist_max))
                endif

             endif
     endif 

   ! 0.999_r8 is added to prevent infinite 'ql_st' at the end of instratus_condensate
   ! computed after updating 'qi_st'.  

     aist = max(0._r8,min(aist,0.999_r8))

     aist_out(i) = aist

     enddo

   return
   end subroutine aist_vector



!---------------------Lin Macrop-------------------------
!    SUBROUTINE gaussj(a,n,np,b,m,mp)
!    INTEGER m,mp,n,np,NMAX
!    real(r8) a(np,np),b(np,mp)
!    real(r8) aa(np,np),bb(np,mp)
!    PARAMETER (NMAX=50)
!    INTEGER i,icol,irow,j,k,l,ll,ii,jj,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
!    real(r8) big,dum,pivinv
!
!    aa(:,:) = a(:,:)
!    bb(:,:) = b(:,:)
!
!    !---------------LiXH add----------------
!    ! icol should be initialized 
!    icol = 0
!    irow = 0
!    !---------------LiXH add----------------
!
!    do 11 j=1,n
!      ipiv(j)=0
!11  continue
!    do 22 i=1,n
!      big=0._r8
!      do 13 j=1,n
!        if(ipiv(j).ne.1)then
!          do 12 k=1,n
!            if (ipiv(k).eq.0) then
!              if (abs(a(j,k)).ge.big)then
!                big=abs(a(j,k))
!                irow=j
!                icol=k
!              endif
!            else if (ipiv(k).gt.1) then
!              do ii = 1, np
!              do jj = 1, np
!                 print*, "gassj", ii, jj, aa(ii,jj), bb(ii,1)
!              end do
!              end do   
!              call endrun('singular matrix in gaussj 1')
!            endif
!12        continue
!        endif
!13    continue
!
!    !---------------LiXH test----------------
!    if(icol .eq. 0 .and. irow .eq. 0)then
!        print*,'rank=',mpi_rank()
!        print*,'ERROR:',i,ipiv(1:n)
!        print*,'a:',abs(a(1,1)),abs(a(1,2)),abs(a(2,1)),abs(a(2,2))
!        call endrun('subroutine gaussj, in cldwat2m_macro.F90')
!    end if
!    !---------------LiXH test----------------
!
!      ipiv(icol)=ipiv(icol)+1
!
!      if (irow.ne.icol) then
!        do 14 l=1,n
!          dum=a(irow,l)
!          a(irow,l)=a(icol,l)
!          a(icol,l)=dum
!14      continue
!        do 15 l=1,m
!          dum=b(irow,l)
!          b(irow,l)=b(icol,l)
!          b(icol,l)=dum
!15      continue
!      endif
!      indxr(i)=irow
!      indxc(i)=icol
!      if (a(icol,icol).eq.0._r8) then
!          do ii = 1, np
!          do jj = 1, np
!             print*, "gaussj", ii, jj, aa(ii,jj), bb(ii,1)
!          end do
!          end do   
!          call endrun('singular matrix in gaussj 2')
!      endif 
!      pivinv=1._r8/a(icol,icol)
!      a(icol,icol)=1._r8
!      do 16 l=1,n
!        a(icol,l)=a(icol,l)*pivinv
!16    continue
!      do 17 l=1,m
!        b(icol,l)=b(icol,l)*pivinv
!17    continue
!      do 21 ll=1,n
!        if(ll.ne.icol)then
!          dum=a(ll,icol)
!          a(ll,icol)=0._r8
!          do 18 l=1,n
!            a(ll,l)=a(ll,l)-a(icol,l)*dum
!18        continue
!          do 19 l=1,m
!            b(ll,l)=b(ll,l)-b(icol,l)*dum
!19        continue
!        endif
!21      continue
!22    continue
!      do 24 l=n,1,-1
!        if(indxr(l).ne.indxc(l))then
!          do 23 k=1,n
!            dum=a(k,indxr(l))
!            a(k,indxr(l))=a(k,indxc(l))
!            a(k,indxc(l))=dum
!23        continue
!        endif
!24    continue
!
!    return
!    end subroutine gaussj
!---------------------Lin Macrop-------------------------

 end module cldwat2m_macro_lin
