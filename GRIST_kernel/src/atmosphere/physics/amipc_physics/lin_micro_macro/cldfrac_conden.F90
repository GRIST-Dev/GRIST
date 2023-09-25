!===================================================================================================
!
!  Created by LiXiaohan on 20/08/10
!
!  Condensation in cloud for Lin Macrophysics
!
!===================================================================================================

module cldfrac_conden_mod
    use grist_wv_saturation,            only: qsat_water
    use grist_constants,                only: r8, cp, latvap, latice, gravity, rdry, pi
    use grist_nml_module,               only: nlev, nlevp
    use grist_handle_error,             only: endrun

    ! the unit of mass flux from deep convection is mb/s
    implicit none
    private
    save
 
    public :: subgrid_var,cldfrac_conden_single_only

    contains

    subroutine subgrid_var(ncol,rpdel_in,top_lev,p_in,z_in,z_face,T_in,qv_in,ql_in, &
                           lengi_in,shi_in,wstarPBL_in,        &
                           qtu_shal, umf_shal, cnt_shal, cnb_shal, thlu_shal,       &
                           sgm_out)

    integer,intent(in) :: ncol, top_lev
    real(r8),intent(in) :: p_in(nlev, ncol) ! pressure
    real(r8),intent(in) :: z_in(nlev, ncol) ! height
    real(r8),intent(in) :: z_face(nlevp, ncol) ! height
    real(r8),intent(in) :: rpdel_in(nlev, ncol) ! pressure difference
    real(r8),intent(in) :: T_in(nlev, ncol) ! temperature
    real(r8),intent(in) :: qv_in(nlev, ncol)! water vapor mixing ratio
    real(r8),intent(in) :: ql_in(nlev, ncol) ! liquid water mixing ratio
    real(r8),intent(in) :: lengi_in(nlevp, ncol) ! master turbulent length scale
    real(r8),intent(in) :: shi_in(nlevp, ncol) ! stability function -- heat
    real(r8),intent(in) :: wstarPBL_in(ncol )  ! convective vertical velocity scale
    
    !shallow
    real(r8),intent(in) :: qtu_shal(nlevp, ncol) ! total specific humidity in updraft
    real(r8),intent(in) :: thlu_shal(nlevp, ncol) ! total liquid poten. temp. in updraft
    real(r8),intent(in) :: umf_shal(nlevp, ncol) ! shallow mass flux 
    real(r8),intent(in) :: cnt_shal(ncol ) ! cloud top
    real(r8),intent(in) :: cnb_shal(ncol ) ! cloud base 
    
    real(r8),intent(out) :: sgm_out(nlev, ncol) ! total subgrid scale variance
    
    
    real(r8),parameter :: c1_gp = 0.32_r8 ! c1=2*ck/cab
    real(r8),parameter :: kafa_gp = 7.44 ! consant for calculation of ld
    real(r8),parameter :: p00_gp = 100000._r8 ! reference surface pressure
    real(r8),parameter :: lkfix_gp=50._r8 ! assumed fixed turbulent length
    real(r8),parameter :: rcon22 = 3.57_r8 !default coefficient
    real(r8),parameter :: rcon23 = rcon22 !coefficient
    real(r8),parameter :: rcon33 = 1.0_r8 ! coefficient
    real(r8),parameter :: sgm_min = 0._r8
    
    !local variables
    integer :: i,j,k
    
    real(r8) :: exner_gp(nlev, ncol)
    real(r8) :: tl_gp(nlev, ncol) ! liquid water temperature
    real(r8) :: thl_gp(nlev, ncol) ! liquid potential temp.
    real(r8) :: thl_face(nlevp, ncol) ! liquid potential temp.
    real(r8) :: qw_gp(nlev, ncol) ! total specific humidity
    real(r8) :: qw_face(nlevp, ncol) ! total specific humidity
    real(r8) :: rho_gp(nlev, ncol) ! density
    
    
    real(r8) :: a_gp(nlev, ncol) ! coefficent "a"
    real(r8) :: b_gp(nlev, ncol)! coefficient "b"
    real(r8) :: beta_gp(nlev, ncol) ! coefficient "beta"
    real(r8) :: qs_gp(nlev, ncol) ! saturated specific humidity
    real(r8) :: es_gp(nlev, ncol) ! satureated water vapor pressure
    real(r8) :: qstl_gp(nlev, ncol) ! saturation specific humidity for liquid water temp.
    real(r8) :: estl_gp(nlev, ncol) ! saturation pressure for liquid water temp.
    
    real(r8) :: dqwdz(nlev, ncol) ! vertical gradient of total water
    real(r8) :: dthldz(nlev, ncol) ! vertical gradient of liquid potential temperature
    
    real(r8) :: cld_dep2(ncol ) ! cloud depth from shallow convection
    
    real(r8) :: qw2(nlev, ncol) ! qw*qw
    real(r8) :: thl2(nlev, ncol) ! thl*thl
    real(r8) :: qwthl(nlev, ncol) ! qu*thl
    
    integer :: cnt_int(ncol ) ! cloud top
    integer :: cnb_int(ncol ) ! cloud base
    
    real(r8) :: sgm_s1(nlev, ncol) ! variance from turbulent
    real(r8) :: sgm_s2(nlev, ncol) ! variance from shallow convection
    real(r8) :: sgm_s3(nlev, ncol) ! variance in free atmosphere

    real(r8) :: lengi, shi, qtu, thlu, umf
    !!!!!!!!!!!!
    ! intialize
    !!!!!!!!!!!!
    
    ! output variable
    !-------------LiXH's tuning-----------
    sgm_out(:nlev,:ncol) = 0._r8
    !-------------LiXH's tuning-----------
    
    ! local variables
    exner_gp(:nlev,:ncol) = 0._r8
    tl_gp(:nlev,:ncol) = 0._r8
    thl_gp(:nlev,:ncol) = 0._r8
    thl_face(:nlevp,:ncol) = 0._r8
    qw_gp(:nlev,:ncol) = 0._r8
    qw_face(:nlevp,:ncol) = 0._r8
    rho_gp(:nlev,:ncol) = 0._r8
    
    a_gp(:nlev,:ncol) = 0._r8
    b_gp(:nlev,:ncol)= 0._r8
    beta_gp(:nlev,:ncol)= 0._r8
    qs_gp(:nlev,:ncol)= 0._r8
    es_gp(:nlev,:ncol)= 0._r8
    qstl_gp(:nlev,:ncol)= 0._r8
    estl_gp(:nlev,:ncol)= 0._r8
    
    dqwdz(:nlev,:ncol) = 0._r8
    dthldz(:nlev,:ncol) = 0._r8
    
    qw2(:nlev,:ncol) = 0._r8
    thl2(:nlev,:ncol) = 0._r8
    qwthl(:nlev,:ncol) = 0._r8
    
    cnt_int(:ncol) = int(cnt_shal(:ncol))
    cnb_int(:ncol) = int(cnb_shal(:ncol))
    
    sgm_s1(:nlev,:ncol) = 0._r8
    sgm_s2(:nlev,:ncol) = 0._r8
    sgm_s3(:nlev,:ncol) = 0._r8

    !!!!!!!!!!!!!!!!!!!!!!
    ! calculation starts.
    !!!!!!!!!!!!!!!!!!!!!!
    
    ! full_level
    do i=1,ncol
        do k=1,nlev
    
            exner_gp(k,i) = (p_in(k,i)/p00_gp)**(rdry/cp)
            qw_gp(k,i) = qv_in(k,i)+ql_in(k,i)
            tl_gp(k,i) = T_in(k,i)-latvap/cp*ql_in(k,i)
            thl_gp(k,i) = tl_gp(k,i)/exner_gp(k,i)
            rho_gp(k,i) = p_in(k,i)/(rdry*T_in(k,i))
            
            
            call qsat_water(tl_gp(k,i),p_in(k,i),estl_gp(k,i),qstl_gp(k,i))
            call qsat_water(t_in(k,i),p_in(k,i),es_gp(k,i),qs_gp(k,i))
            
            beta_gp(k,i) = 0.622_r8*latvap*qstl_gp(k,i)/rdry/tl_gp(k,i)**2
            
            a_gp(k,i) = 1._r8/(1._r8+beta_gp(k,i)*latvap/cp)
            b_gp(k,i) = exner_gp(k,i)*beta_gp(k,i)/(1._r8+beta_gp(k,i)*latvap/cp)
    
        enddo
    enddo
    
    !Compute dqwdz,dthldz based on "slope" from uwshcu, LiXH
    dqwdz  = slope(nlev,ncol,qw_gp,z_in)
    dthldz = slope(nlev,ncol,thl_gp,z_in) 

    !================================================================!
    !=========================turbulent effect=======================!

    !--------------bk, delete---------------
    !do i=1,ncol
    !    do k = top_lev,nlev-1
    !        ! part I: caused by turbulent process
    !        !dqwdz(k,i) = -rho(k,i)*gravity*(qw(k+1,i)-qw(k,i))*rpdel_in(k,i)
    !        !dthldz(k,i) =-rho(k,i)*gravity*(thl(k+1,i)-thl(k,i))*rpdel_in(k,i) 
    !        
    !        dqwdz(k,i) = (qw_gp(k+1,i)-qw_gp(k,i))/(z_in(k+1,i)-z_in(k,i))
    !        dthldz(k,i) = (thl_gp(k+1,i)-thl_gp(k,i))/(z_in(k+1,i)-z_in(k,i))
    !    enddo
    !
    !    dqwdz(nlev,i) = 0._r8
    !    dthldz(nlev,i) = 0._r8
    !enddo
    !--------------bk, delete---------------


    do i=1,ncol
        do k=top_lev,nlev
          !  qw2(k,i) = 2.0_r8*rcon22*lengi_in(k+1,i)*lengi_in(k+1,i)*shi_in(k+1,i)*dqwdz(k,i)**2_r8
          !  thl2(k,i) = 2.0_r8*rcon22*lengi_in(k+1,i)*lengi_in(k+1,i)*shi_in(k+1,i)*dthldz(k,i)**2_r8
          !  qwthl(k,i) = 2.0_r8*rcon22*lengi_in(k+1,i)*lengi_in(k+1,i)*shi_in(k+1,i)*dqwdz(k,i)*dthldz(k,i)
          !  sgm_s1(k,i) = 0.25_r8*2.0_r8*rcon22*lengi_in(k+1,i)*lengi_in(k+1,i)*shi_in(k+1,i)*(a_gp(k,i)**2*dqwdz(k,i)**2 + &
          !                b_gp(k,i)**2*dthldz(k,i)**2 - 2._r8*a_gp(k,i)*b_gp(k,i)*dqwdz(k,i)*dthldz(k,i))
 
            lengi = (lengi_in(k+1,i)+lengi_in(k,i))*0.5_r8
            shi   = (shi_in(k+1,i)+shi_in(k,i))*0.5_r8
            sgm_s1(k,i) = rcon22*lengi*lengi*shi*(a_gp(k,i)**2*dqwdz(k,i)**2 + &
                          b_gp(k,i)**2*dthldz(k,i)**2 - 2._r8*a_gp(k,i)*b_gp(k,i)*dqwdz(k,i)*dthldz(k,i))


           ! LiXH add sgm for free atmosphere:
!           if(lengi .le. 0._r8 .or. shi .le. 0._r8)then
!               ! LiXH use 50 m for lengi 
!               lengi = 50._r8          
!               sgm_s3(k,i) = rcon23*lengi*lengi*(a_gp(k,i)**2*dqwdz(k,i)**2 + &
!                             b_gp(k,i)**2*dthldz(k,i)**2 - 2._r8*a_gp(k,i)*b_gp(k,i)*dqwdz(k,i)*dthldz(k,i))  
!           end if 
        enddo
    enddo 
    !=============================================================!
    !=====================shallow convection effect===============!
    
    do i=1,ncol
        do k=top_lev,nlev
    
            ! cloud depth
            cld_dep2(i) = z_in(cnt_int(i),i)-z_in(cnb_int(i),i) 
            
            if(cld_dep2(i).le.0)then
                    cld_dep2(i) =0._r8
            endif
            
            qtu = (qtu_shal(k,i)+qtu_shal(k+1,i))*0.5_r8
            thlu= (thlu_shal(k,i)+thlu_shal(k+1,i))*0.5_r8
            umf = (umf_shal(k,i)+umf_shal(k+1,i))*0.5_r8
            if((wstarPBL_in(i).ne.0._r8).and.(qtu.ne.0._r8).and.(thlu.ne.0._r8))then
             !   sgm_s2(k,i) = 0.5_r8*rcon33*umf_shal(k,i)/rho_gp(k,i)*cld_dep2(i)/wstarPBL_in(i)&
             !                *(a_gp(k,i)**2*(qtu_shal(k,i)-qw_gp(k,i))*dqwdz(k,i)&
             !                 +b_gp(k,i)**2*(thlu_shal(k,i)-thl_gp(k,i))*dthldz(k,i) &
             !                 -a_gp(k,i)*b_gp(k,i)*((thlu_shal(k,i)-thl_gp(k,i))*dqwdz(k,i)+(qtu_shal(k,i)-qw_gp(k,i))*dthldz(k,i)))

                sgm_s2(k,i) = 0.5_r8*rcon33*umf/rho_gp(k,i)*cld_dep2(i)/wstarPBL_in(i)&
                             *(a_gp(k,i)**2*(qtu-qw_gp(k,i))*dqwdz(k,i)&
                              +b_gp(k,i)**2*(thlu-thl_gp(k,i))*dthldz(k,i) &
                              -a_gp(k,i)*b_gp(k,i)*((thlu-thl_gp(k,i))*dqwdz(k,i)+(qtu-qw_gp(k,i))*dthldz(k,i)))
 
            endif
    
        enddo
    enddo
    
    ! calcuate the total subgrid scale variance
    do i=1,ncol
        do k=top_lev,nlev
        !--------------------LiXH's tuning------------------------
            sgm_out(k,i) = abs(sgm_s1(k,i))+abs(sgm_s2(k,i))
            if(sgm_out(k,i) .le. 0._r8)then
            sgm_out(k,i) = abs(sgm_s3(k,i))
            end if
        !--------------------LiXH's tuning------------------------

    !--------------LiXH Test-------------
    !if(sgm_out(k,i) .le. 0._r8)then
    !    print*,'LiXH:',k,dqwdz(k,i),dthldz(k,i),(lengi_in(k+1,i)+lengi_in(k,i))*0.5_r8,a_gp(k,i),b_gp(k,i)

    !end if
    !--------------LiXH Test-------------

        enddo
    enddo

    end subroutine subgrid_var


!    subroutine cldfrac_conden(p_in,T_in,qv_in,ql_in,sgm_in,cld_frac_out,conden_out,ncol,ilev)
!    
!    implicit none
!    
!    integer, intent(in) :: ncol, ilev
!    real(r8),intent(in) :: p_in(ncol ) ! pressure
!    real(r8),intent(in) :: T_in(ncol ) ! temperature
!    real(r8),intent(in) :: qv_in(ncol )! water vapor specific humidity
!    real(r8),intent(in) :: ql_in(ncol ) ! liquid water specific humidity
!    real(r8),intent(in) :: sgm_in(ncol ) ! total subgrid scale variance
!    
!    real(r8),intent(out) :: cld_frac_out(ncol ) ! diagnosed cloud fraction from Gauss-PDF
!    real(r8),intent(out) :: conden_out(ncol ) ! diagnosed cloud condensate
!    
!    real(r8),parameter:: sgm_min = 1.e-8_r8
!    real(r8),parameter:: p00 = 100000._r8
!    
!    !local variables
!    integer :: i,j,k
!    
!    real(r8) :: conden(ncol ) !ql/(2*sigma(s))
!    real(r8) :: cld_frac(ncol )
!    
!    real(r8) :: exner_gp(ncol )
!    real(r8) :: qw_gp(ncol )
!    real(r8) :: tl_gp(ncol )
!    real(r8) :: thl_gp(ncol )
!    real(r8) :: rho_gp(ncol )
!    
!    real(r8) :: qs_gp
!    real(r8) :: es_gp
!    real(r8) :: qstl_gp ! saturation specific humidity for liquid water temp.
!    real(r8) :: estl_gp ! saturation pressure for liquid water temp.
!    
!    real(r8) :: a_gp
!    real(r8) :: b_gp
!    real(r8) :: beta_gp
!    real(r8) :: deltaq_gp
!    real(r8) :: Q1_gp
!
!    exner_gp(:ncol) = 0._r8
!    qw_gp(:ncol) = 0._r8
!    tl_gp(:ncol)= 0._r8
!    thl_gp(:ncol) = 0._r8
!    rho_gp(:ncol) = 0._r8
!    
!    conden(:ncol) = 0._r8
!    cld_frac(:ncol) = 0._r8
!    
!    
!    do k=1,ncol
!    
!    qs_gp = 0._r8
!    es_gp = 0._r8
!    estl_gp = 0._r8
!    qstl_gp = 0._r8
!    
!    a_gp = 0._r8
!    b_gp = 0._r8
!    beta_gp = 0._r8
!    deltaq_gp = 0._r8
!    Q1_gp = 0._r8
!    
!    
!    exner_gp(k) = (p_in(k)/p00)**0.286_r8
!    qw_gp(k) = qv_in(k)+ql_in(k)
!    tl_gp(k) = T_in(k)-latvap/cp*ql_in(k)
!    thl_gp(k) = tl_gp(k)/exner_gp(k)
!    rho_gp(k) = p_in(k)/(rdry*T_in(k))
!    
!    call qsat_water(tl_gp(k),p_in(k),estl_gp,qstl_gp)
!    call qsat_water(t_in(k),p_in(k),es_gp,qs_gp)
!    
!    beta_gp = 0.622_r8*latvap*qstl_gp/rdry/tl_gp(k)**2
!    a_gp = 1._r8/(1._r8+beta_gp*latvap/cp)
!    b_gp = exner_gp(k)*beta_gp/(1._r8+beta_gp*latvap/cp)
!    
!    deltaq_gp = qw_gp(k)-qstl_gp
!    
!    if(sgm_in(k).ne.0._r8)then
!            Q1_gp = a_gp*deltaq_gp/(2._r8*sqrt(sgm_in(k)))
!    else
!            Q1_gp = 0._r8
!    endif
!
!    ! calculate cloud fraction and condensate
!    if(sqrt(sgm_in(k))>sgm_min)then
!        cld_frac(k) = 0.5_r8*(1._r8+erf(Q1_gp/sqrt(2._r8)))
!        conden(k) = cld_frac(k)*Q1_gp+exp(-0.5_r8*Q1_gp**2._r8)/sqrt(2._r8*pi)
!    
!        cld_frac_out(k) = min(1._r8,max(0._r8,cld_frac(k)))
!        conden_out(k) = conden(k)*2._r8*sqrt(sgm_in(k))
!        conden_out(k) = max(0._r8,conden_out(k))
!    
!    else
!        if(deltaq_gp.le.0._r8)then
!          cld_frac_out(k) = 0._r8
!          conden_out(k) = 0._r8
!    
!        else
!            !-----------LiXH test------------
!            print*,'advection of ql: k=',ilev,conden_out(k),qw_gp(k),qstl_gp,'We should check code!!!'
!            !-----------LiXH test------------
!          cld_frac_out(k) = 1._r8
!          conden_out(k) = max(0._r8,a_gp*deltaq_gp)
!    
!        endif
!    endif
!    
!    ! qinyi 2016-11-18 21:59:05
!    ! test by simple calculation indicates that there is unstable oscillation while
!    ! cldfrac is less than 1e-4.
!    ! the exact reason for this is the too small deltaq could cause this.
!    ! so in order to eliminate this inapproximate value, I set the threshold for
!    ! cloud fraction.
!    
!    if((conden_out(k).eq.0._r8).or.(cld_frac_out(k).lt.1.e-4_r8))then
!            conden_out(k) = 0._r8
!            cld_frac_out(k) = 0._r8
!    endif
!    
!    enddo !k=1,ncol
!    
!    end subroutine cldfrac_conden


    subroutine cldfrac_conden_single_only(ncol, top_lev, p_in,T_in,qv_in,ql_in,sgm_in, &
               cld_frac_out,conden_out,G,deltaq_sat,deltaq_uns,Q1_sat,Q1_uns,adjust_factor)
    
    integer,intent(in) :: ncol, top_lev
    real(r8),intent(in) :: p_in(nlev, ncol)
    real(r8),intent(in) :: T_in(nlev, ncol)
    real(r8),intent(in) :: qv_in(nlev, ncol)
    real(r8),intent(in) :: ql_in(nlev, ncol)
    real(r8),intent(in) :: sgm_in(nlev, ncol)
    
    real(r8),intent(out) :: cld_frac_out(nlev, ncol)
    real(r8),intent(out) :: conden_out(nlev, ncol)
    real(r8),intent(out) :: G(nlev, ncol)
    real(r8),intent(out) :: deltaq_sat(nlev, ncol)
    real(r8),intent(out) :: deltaq_uns(nlev, ncol)
    real(r8),intent(out) :: Q1_sat(nlev, ncol)
    real(r8),intent(out) :: Q1_uns(nlev, ncol)
    real(r8),intent(out) :: adjust_factor(nlev, ncol)
    
    real(r8),parameter:: sgm_min = 1.e-8_r8
    real(r8),parameter:: p00 = 100000._r8
    
    !local variables
    integer :: i,k
    
    real(r8) :: conden !ql/(2*sigma(s))
    real(r8) :: cld_frac
    
    real(r8) :: exner_gp
    real(r8) :: tl_gp
    real(r8) :: thl_gp
    real(r8) :: rho_gp
    real(r8) :: qw_gp
    real(r8) :: theta_gp
    
    real(r8) :: qstl_gp
    real(r8) :: estl_gp
    
    real(r8) :: a_gp
    real(r8) :: b_gp
    real(r8) :: beta_gp
    real(r8) :: deltaq_gp
    real(r8) :: Q1_gp
    
    real(r8) :: adjust_mid
    
    ! initialized
    ! output variables
    cld_frac_out = 0._r8
    conden_out = 0._r8
    G = 0._r8
    deltaq_sat = 0._r8
    deltaq_uns = 0._r8
    Q1_sat = 0._r8
    Q1_uns = 0._r8
    adjust_factor = 0._r8
    
    ! calculation starts
    do i =  1, ncol
    do k =  top_lev, nlev
        exner_gp = (p_in(k,i)/p00)**0.286_r8
        tl_gp = T_in(k,i)-latvap/cp*ql_in(k,i)
        thl_gp = tl_gp/exner_gp
        theta_gp = T_in(k,i)/exner_gp
        rho_gp = p_in(k,i)/(rdry*T_in(k,i))
        qw_gp = qv_in(k,i)+ql_in(k,i)
    
        call qsat_water(tl_gp,p_in(k,i),estl_gp,qstl_gp)
    
        beta_gp = 0.622_r8*latvap*qstl_gp/rdry/tl_gp**2
        
        a_gp = 1._r8/(1._r8+beta_gp*latvap/cp)
        b_gp = exner_gp*beta_gp/(1._r8+beta_gp*latvap/cp)

        deltaq_gp = qw_gp-qstl_gp
    
        if(sgm_in(k,i).ne.0_r8)then
                Q1_gp = a_gp*deltaq_gp/(2._r8*sqrt(sgm_in(k,i)))
        else
                Q1_gp = 0._r8
        endif
    
        ! calculate cloud fraction and condensate
        if(sqrt(sgm_in(k,i))>sgm_min)then
            cld_frac = 0.5_r8*(1._r8+erf(Q1_gp/sqrt(2._r8)))
            conden = cld_frac*Q1_gp+exp(-0.5_r8*Q1_gp**2._r8)/sqrt(2._r8*pi)
            
            cld_frac_out(k,i) = min(1._r8,max(0._r8,cld_frac))
            conden_out(k,i) = max(0._r8,conden*2._r8*sqrt(sgm_in(k,i)))
        else
            if(deltaq_gp.le.0._r8)then
              cld_frac_out(k,i) = 0._r8
              conden_out(k,i) = 0._r8
            else
              cld_frac_out(k,i) = 1._r8
              conden_out(k,i) = max(0._r8,a_gp*deltaq_gp)

            endif
        endif

        ! qinyi 2016-11-18 22:01:15
        if((conden_out(k,i).le.0._r8).or.(cld_frac_out(k,i).lt.1.e-4_r8))then
                conden_out(k,i) = 0._r8
                cld_frac_out(k,i) = 0._r8
        endif
    
        G(k,i) = deltaq_gp
    
        ! qinyi 2017-1-1 21:55:23
        ! this part is for the consistency b/t liquid water latent heating and temperature change
        adjust_mid = -1._r8*a_gp*(cld_frac_out(k,i)+1._r8/sqrt(2._r8*pi)*exp(-0.5_r8*Q1_gp**2)*(-1._r8)*Q1_gp)*0.622_r8*latvap*qstl_gp/(rdry*tl_gp**2*thl_gp**2)
        !adjust_factor(k,i) = cp*T_in(k,i)/(latvap*theta_gp)/(cp*T_in(k,i)/(latvap*theta_gp)-adjust_mid)
        adjust_factor(k,i) = max(0._r8,min(1._r8,cp*T_in(k,i)/(latvap*theta_gp)/(cp*T_in(k,i)/(latvap*theta_gp)-adjust_mid)) )     !LiXH, [0-1]

        !------------LiXH add------------
        !if(conden_out(k,i) .eq. 0._r8)adjust_factor(k,i) = 1._r8
        !------------LiXH add------------

        if(deltaq_gp.ge.0._r8) then
            deltaq_sat(k,i) = deltaq_gp
        
            if(sgm_in(k,i).ne.0_r8)then
            Q1_sat(k,i) = a_gp*deltaq_gp/(2._r8*sqrt(sgm_in(k,i)))
            else
            Q1_sat = 0._r8
            endif
        else
            deltaq_uns(k,i) = deltaq_gp
        
            if(sgm_in(k,i).ne.0_r8)then
            Q1_uns(k,i) = a_gp*deltaq_gp/(2._r8*sqrt(sgm_in(k,i)))
            else
            Q1_uns(k,i) = 0._r8
            endif
        endif
    enddo
    enddo

    end subroutine cldfrac_conden_single_only

    function slope(mkx,ncol,field,z0)
    ! Note that k is from top to bottom
    integer,  intent(in) :: mkx
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: field(mkx,ncol)
    real(r8), intent(in) :: z0(mkx,ncol)
    real(r8)             :: slope(mkx,ncol)

    real(r8)             :: below
    real(r8)             :: above
    integer              :: k, i

    do i = 1, ncol
        below = ( field(nlev-1,i) - field(nlev,i) ) / ( z0(nlev-1,i) - z0(nlev,i) )
        do k = nlev-1, 1, -1
            above = ( field(k,i) - field(k+1,i) ) / ( z0(k,i) -  z0(k+1,i) )
            if( above .gt. 0._r8 ) then
                slope(k+1,i) = max(0._r8,min(above,below))
            else
                slope(k+1,i) = min(0._r8,max(above,below))
            end if
            below = above
        end do
        slope(1,i) = slope(2,i)
    end do

    end function slope

!    subroutine cldfrac_conden_single(p_in,T_in,qv_in,ql_in,sgm_in,cld_frac_out,conden_out,G)
!    
!    real(r8),intent(in) :: p_in
!    real(r8),intent(in) :: T_in
!    real(r8),intent(in) :: qv_in
!    real(r8),intent(in) :: ql_in
!    real(r8),intent(in) :: sgm_in
!    
!    real(r8),intent(out) :: cld_frac_out
!    real(r8),intent(out) :: conden_out
!    real(r8),intent(out) :: G
!    
!    real(r8),parameter:: sgm_min = 1.e-8_r8 ! minimum value of SGV
!    real(r8),parameter:: p00 = 100000._r8
!    
!    !local variables
!    integer :: i,j,k
!    
!    real(r8) :: conden !ql/(2*sigma(s))
!    real(r8) :: cld_frac
!    
!    real(r8) :: exner_gp
!    real(r8) :: tl_gp
!    real(r8) :: thl_gp
!    real(r8) :: rho_gp
!    real(r8) :: qw_gp
!    
!    real(r8) :: qs_gp
!    real(r8) :: es_gp
!    real(r8) :: qstl_gp
!    real(r8) :: estl_gp
!    real(r8) :: Q1_gp
!    real(r8) :: deltaq_gp
!    real(r8) :: a_gp
!    real(r8) :: b_gp
!    real(r8) :: beta_gp
!    
!    !!!!!!!!!!!!!!!
!    ! initialized 
!    !!!!!!!!!!!!!!!
!    
!    ! output variables
!    cld_frac_out = 0._r8
!    conden_out = 0._r8
!    G = 0._r8
!    
!    ! local variables
!    
!    conden = 0._r8
!    cld_frac = 0._r8
!    
!    exner_gp = 0._r8
!    tl_gp = 0._r8
!    thl_gp = 0._r8
!    rho_gp = 0._r8
!    qw_gp = 0._r8
!    
!    qs_gp = 0._r8
!    es_gp = 0._r8
!    qstl_gp = 0._r8
!    estl_gp = 0._r8
!    Q1_gp = 0._r8
!    deltaq_gp = 0._r8
!    a_gp = 0._r8
!    b_gp = 0._r8
!    beta_gp = 0._r8
!    
!    !!!!!!!!!!!!!!!!!!!!!!!!!
!    !calculation starts here
!    !!!!!!!!!!!!!!!!!!!!!!!!!
!    
!    exner_gp = (p_in/p00)**0.286_r8
!    tl_gp = T_in-latvap/cp*ql_in
!    thl_gp = tl_gp/exner_gp
!    rho_gp = p_in/(rdry*T_in)
!    qw_gp = qv_in+ql_in
!    
!    call qsat_water(tl_gp,p_in,estl_gp,qstl_gp)
!    call qsat_water(t_in,p_in,es_gp,qs_gp)
!    
!    beta_gp = 0.622_r8*latvap*qstl_gp/rdry/tl_gp**2
!    
!    a_gp = 1._r8/(1._r8+beta_gp*latvap/cp)
!    b_gp = exner_gp*beta_gp/(1._r8+beta_gp*latvap/cp)
!    
!    deltaq_gp = qw_gp-qstl_gp
!    
!    if(sgm_in.ne.0_r8)then
!            Q1_gp = a_gp*deltaq_gp/(2._r8*sqrt(sgm_in))
!    else
!            Q1_gp = 0._r8
!    endif
!    
!    ! calculate cloud fraction and condensate
!    if(sqrt(sgm_in)>sgm_min)then
!        cld_frac = 0.5_r8*(1._r8+erf(Q1_gp/sqrt(2._r8)))
!        conden = cld_frac*Q1_gp+exp(-0.5_r8*Q1_gp**2._r8)/sqrt(2._r8*pi)
!        
!        cld_frac_out = min(1._r8,max(0._r8,cld_frac))
!        conden_out = conden*2._r8*sqrt(sgm_in)
!        conden_out = max(0._r8,conden_out)
!    
!    else
!        if(deltaq_gp.le.0._r8)then
!          cld_frac = 0._r8
!          conden = 0._r8
!          cld_frac_out = min(1._r8,max(0._r8,cld_frac))
!          conden_out = conden*2._r8*sqrt(sgm_in)
!          conden_out = max(0._r8,conden_out)
!    
!        else
!          cld_frac = 1._r8
!          conden_out = a_gp*deltaq_gp
!          cld_frac_out = min(1._r8,max(0._r8,cld_frac))
!          conden_out = max(0._r8,conden_out)
!    
!        endif
!    endif
!    
!    if((conden_out.eq.0._r8).or.(cld_frac_out.lt.1.e-4_r8))then
!    conden_out = 0._r8
!    cld_frac_out = 0._r8
!    endif
!    
!    G = deltaq_gp
!    
!    end subroutine cldfrac_conden_single

end module cldfrac_conden_mod
