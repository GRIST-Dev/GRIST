!-----------------------------------------------------------------------
! Created on 2016 
! Version 1.0
! Description: Time Manager
! Revision history:
!    This module should not use grist_nml_module, because nml will
!    need grist_time_manager; if you want to have info from nml, such as
!    start_ymd, start_tod, working_mode, use input stream. This has been
!    changed this to this present form, if something inconsistent occur, 
!    e.g., SCM, sligtly modify the code but follow this practice.
!    Revision: fix bug when using leap year, yizhang: 20211011
!-----------------------------------------------------------------------

 module grist_time_manager

   use grist_constants,      only: i4, r8, day2sec, i8
   use grist_handle_error,   only: endrun
 
   implicit none

   private

   public :: get_current_date,      &
             get_curr_calday,       &
             get_curr_date,         &
             get_curr_cdate,        &
             get_start_date,        &
             set_time_float_from_date,&
             set_date_from_time_float, &
             intToChar5 
!
! we have two logics here:
!(i) leap year matters; (ii) this year is leap
! whether "leap year matters" is controlled by compiling options
! whether "this year is leap: is controlled by internal flow
! if this year is leap but we donot consider leap-year, then counting leap-year must be
! same as counting normal years, to regress old results
!
   integer(i4) :: month(12)      =(/31,28,31,30,31,30,31,31,30,31,30,31/)
   integer(i4) :: allday         = 365

   integer(i4) :: norm_month(12) =(/31,28,31,30,31,30,31,31,30,31,30,31/)
   integer(i4) :: norm_allday    = 365
#ifdef USE_LEAP_YEAR
   integer(i4) :: leap_month(12) =(/31,29,31,30,31,30,31,31,30,31,30,31/)
   integer(i4) :: leap_allday    = 366
#else
   integer(i4) :: leap_month(12) =(/31,28,31,30,31,30,31,31,30,31,30,31/)
   integer(i4) :: leap_allday    = 365
#endif

 CONTAINS

!------------------------------------------------------------
! zy: if we use leap year, then set this; otherwise ignore it
!------------------------------------------------------------

  subroutine check_leap_year_present_step(start_ymd, start_tod, nstep, dtime)
! io
    integer(i4),  intent(in) :: start_ymd, start_tod, nstep
    real(r8), intent(in)     :: dtime
! local
    integer(i4)              :: yr, mn, dy, dayAfterIntegration
    real(r8)                 :: inc, sc
    integer(i4)              :: i

    yr  = start_ymd/10000
    mn  =(start_ymd-yr*10000)/100
    dy  = start_ymd-yr*10000-mn*100
    inc = nstep*dtime
    sc  = start_tod+inc
    dayAfterIntegration = 0
! check whether start year is leap
    call reset_month_allday(yr)

    do while(sc.ge.day2sec)
       dy = dy+1       ! total day since integration plus initial day of that month
       sc = sc-day2sec ! sc of this day
       dayAfterIntegration = dayAfterIntegration+1
    end do
    !print*,"total day after integration is", dayAfterIntegration

    do i = 1, mn-1
       dy = dy+month(i) ! total day since the start date of beg year
    end do

    do while(dy.gt.allday)
       dy = dy - allday
       yr = yr + 1
       call reset_month_allday(yr)
    end do

    return
  end subroutine check_leap_year_present_step

  subroutine get_current_date(itimestep,dtime,day,sec)
! io    
     integer,          intent(in)  :: itimestep
     real(r8),         intent(in)  :: dtime
     character(len=5), intent(out) :: day
     character(len=5), intent(out) :: sec
! local     
     character(len=5)              :: day_tmp
     character(len=5)              :: sec_tmp
     integer(i8)                   :: nday
     integer(i8)                   :: nsec
     real(r8)                      :: ntime

     ntime = dtime*itimestep     ! total time elapsed
     nday  = int(ntime/int(86400),i8)
     nsec  = int(mod(int(ntime,i8),int(86400)),i8)

     if(nday.lt.10)then
        write(day_tmp,'(i1)') nday
        day=trim('0000')//trim(day_tmp)
     else if(nday.ge.10 .and. nday.lt.100)then
        write(day_tmp,'(i2)') nday
        day=trim('000')//trim(day_tmp)
     else if(nday.ge.100.and. nday.lt.1000)then
        write(day_tmp,'(i3)') nday
        day=trim('00')//trim(day_tmp)
     else if(nday.ge.1000.and. nday.lt.10000)then
        write(day_tmp,'(i4)') nday
        day=trim('0')//trim(day_tmp)
     else
        write(day_tmp,'(i5)') nday
        day=trim(day_tmp)
     end if

     if(nsec.lt.10)then
        write(sec_tmp,'(i1)') nsec
        sec=trim('0000')//trim(sec_tmp)
     else if(nsec.ge.10 .and. nsec.lt.100)then
        write(sec_tmp,'(i2)') nsec
        sec=trim('000')//trim(sec_tmp)
     else if(nsec.ge.100.and. nsec.lt.1000)then
        write(sec_tmp,'(i3)') nsec
        sec=trim('00')//trim(sec_tmp)
     else if(nsec.ge.1000.and. nsec.lt.10000)then
        write(sec_tmp,'(i4)') nsec
        sec=trim('0')//trim(sec_tmp)
     else
        write(sec_tmp,'(i5)') nsec
        sec=trim(sec_tmp)
     end if

     return
   end subroutine get_current_date

! Purpose : calculate julian day.
! LiXH
   subroutine get_curr_calday(start_ymd, start_tod,nstep, dtime, calday)

   use grist_zenith,                       only: orb_params

   ! io
   integer,  intent(in)     :: start_ymd, start_tod
   integer,  intent(in)     :: nstep
   real(r8), intent(in)     :: dtime
   real(r8), intent(out)    :: calday
   ! local
   integer                  :: year, mon, day, i
   real(r8)                 :: inc, current_sec
   logical                  :: leap_year

       year        = start_ymd/10000
       mon         = (start_ymd-year*10000)/100
       day         = start_ymd-year*10000-mon*100

!================================================
!      inc   = nstep*nint(dtime)
! yizhang fix here for very large inc (cannot be int4), affect rad
! module, regression ok
!================================================

       inc         = nstep*dtime
       current_sec = start_tod+inc

!================================================
! yiz: correct here for leap year
!================================================
       call reset_month_allday(year)

       do while(current_sec .ge. day2sec)
          current_sec = current_sec-day2sec
          day = day+1
       end do

       calday = real(day,r8)+real(current_sec,r8)/day2sec

       do i = 1, mon-1
          calday = calday+real(month(i),r8)
       end do

       do while(calday.gt.allday)
          calday = calday - allday
          year   = year + 1
          call reset_month_allday(year)
           ! update orbital parameters, LiXH
          call orb_params(year, .false.)
       end do
   !else
   !    if(mpi_rank()==0)then
   !        print*,'ERROR! You should modify grist_time_manager!'
   !        print*,'start_ymd and start_tod should be set in grist.nml'
   !        call endrun('error in grist_time_manager')
   !    end if
   !end if

   end subroutine get_curr_calday

   subroutine get_curr_date(start_ymd, start_tod, nstep, dtime, yr, mn, dy, sc)
   ! io
   integer,  intent(in)     :: start_ymd, start_tod
   integer,  intent(in)     :: nstep
   real(r8), intent(in)     :: dtime
   integer,  intent(out)    :: yr, mn, dy, sc
   ! local
   real(r8)                 :: inc, current_sec
   integer                  :: i
 
       yr  = start_ymd/10000
       mn  =(start_ymd-yr*10000)/100
       dy  = start_ymd-yr*10000-mn*100
       inc = nstep*dtime
       current_sec = start_tod+inc

       call reset_month_allday(yr)

       do while(current_sec.ge.day2sec)
          current_sec = current_sec-day2sec
          dy = dy + 1
       end do
       sc = nint(current_sec)
!
! set dy to days since the start of startymd's 0101
!
       do i = 1, mn-1
          dy = dy+month(i) ! this month belones to startymd
       end do

       do while(dy.gt.allday)
          dy = dy - allday
          yr = yr + 1
          call reset_month_allday(yr)
       end do
!
! set the dy and mn of current year
!
       mn = 1
       do i = 1, 12
          if(dy.le.month(i)) exit ! this month belones to current year
          dy = dy - month(i)
          mn = mn + 1
       end do
   !else
   !    if(mpi_rank()==0)then
   !        print*,'ERROR! You should modify grist_time_manager!'
   !        print*,'start_ymd and start_tod should be set in grist.nml'
   !        call endrun('error in grist_time_manager')
   !    end if
   !end if
   end subroutine get_curr_date
!
! as get_curr_date but return character
!
   subroutine get_curr_cdate(start_ymd, start_tod, nstep, dtime, year, mon, day, sec)
   ! io
   integer,  intent(in)     :: start_ymd, start_tod
   integer,  intent(in)     :: nstep
   real(r8), intent(in)     :: dtime
   character(len=4), intent(out) :: year
   character(len=2), intent(out) :: mon
   character(len=2), intent(out) :: day
   character(len=5), intent(out) :: sec
   ! local
   integer                  :: yr, mn, dy, sc, i
   real(r8)                 :: inc, current_sec
   character(len=5)         :: tmp
 
       inc = nstep*dtime
       yr  = start_ymd/10000
       mn  =(start_ymd-yr*10000)/100
       dy  = start_ymd-yr*10000-mn*100
       current_sec  = start_tod+inc

       call reset_month_allday(yr)

       do while(current_sec.ge.day2sec)
          dy = dy+1
          current_sec = current_sec-day2sec
       end do
       sc = nint(current_sec)

! let dy be the total days since the 1st day of startymd
       do i = 1, mn-1
          dy = dy+month(i) ! this month belongs to startymd
       end do

       do while(dy.gt.allday)
          dy = dy - allday
          yr = yr + 1
          call reset_month_allday(yr)
       end do

       mn = 1
       do i = 1, 12
          if(dy .le. month(i)) exit ! this month belones to current year
          dy = dy - month(i)
          mn = mn + 1
       end do
! year
       write(year,'(i4)') yr
! month
       if(mn.lt.10)then
          write(tmp,'(i1)') mn
          mon=trim('0')//trim(tmp)
       else
          write(tmp,'(i2)') mn
          mon=trim(tmp)
       end if
! day
       if(dy.lt.10)then
          write(tmp,'(i1)') dy
          day=trim('0')//trim(tmp)
       else
          write(tmp,'(i2)') dy
          day=trim(tmp)
       end if
! sec
       call intToChar5(sc,sec)

      return
   end subroutine get_curr_cdate

   subroutine get_start_date(start_ymd, start_tod, yr, mn, dy, sc)
   ! io
   integer,  intent(in)     :: start_ymd, start_tod
   integer,  intent(out)    :: yr, mn, dy, sc
   ! local
   integer                  :: i
 
   !if(working_mode .eq. 'scm')then
       yr  = start_ymd/10000
       mn  = (start_ymd-yr*10000)/100
       dy  = start_ymd-yr*10000-mn*100
       sc  = start_tod
   !else
   !    if(mpi_rank()==0)then
   !        print*,'ERROR! You should modify grist_time_manager!'
   !        print*,'start_ymd and start_tod should be set in grist.nml'
   !        call endrun('error in grist_time_manager')
   !    end if
   !end if

   end subroutine get_start_date

! YiZ: if use below, please double check whether leap year works ok! I have not checked.
   subroutine set_time_float_from_date1(start_ymd, start_tod, nstep, dtime, time, yr, mn, dy, sc)
   ! io
   integer,  intent(in)  :: start_ymd, start_tod, nstep
   real(r8), intent(in)  :: dtime
   real(r8), intent(out) :: time
   integer,  intent(in)  :: yr, mn, dy, sc
   ! local
   integer               :: s_yr, s_mn, s_dy, s_sc, i
   integer               :: curr_dy

   call check_leap_year_present_step(start_ymd, start_tod, nstep, dtime)

   s_yr  = start_ymd/10000
   s_mn  = (start_ymd-s_yr*10000)/100
   s_dy  = start_ymd-s_yr*10000-s_mn*100
   s_sc  = start_tod

   do i = 1, s_mn-1
      s_dy = s_dy + month(i)
   end do 

   curr_dy = dy
   do i = 1, mn-1
      curr_dy = curr_dy + month(i)
   end do
 
   time = real(yr-s_yr,r8)*allday+real(curr_dy-s_dy,r8)+real(sc-s_sc,r8)/86400.

   end subroutine set_time_float_from_date1
   
   subroutine set_time_float_from_date(start_ymd, start_tod, time, yr, mn, dy, sc)
   ! io
   integer,  intent(in)  :: start_ymd, start_tod
   real(r8), intent(out) :: time
   integer,  intent(in)  :: yr, mn, dy, sc
   ! local
   integer               :: s_yr, s_mn, s_dy, s_sc, i
   integer               :: curr_dy

   s_yr  = start_ymd/10000
   s_mn  = (start_ymd-s_yr*10000)/100
   s_dy  = start_ymd-s_yr*10000-s_mn*100
   s_sc  = start_tod

   do i = 1, s_mn-1
      s_dy = s_dy + month(i)
   end do

   curr_dy = dy
   do i = 1, mn-1
      curr_dy = curr_dy + month(i)
   end do

   time = real(yr-s_yr,r8)*365._r8+real(curr_dy-s_dy,r8)+real(sc-s_sc,r8)/86400.

   end subroutine set_time_float_from_date


   subroutine set_date_from_time_float1(start_ymd, start_tod, nstep, dtime, curr_mod_time, yr, mn, day, nsec)
   !io 
   integer,  intent(in)      :: start_ymd, start_tod, nstep
   real(r8), intent(in)      :: dtime, curr_mod_time
   integer,  intent(out)     :: yr, mn, day, nsec
   ! local
   integer                   :: calday, i
   integer                   :: s_yr, s_mn, s_dy, s_sc
   real(r8)                  :: curr_time, s_time

   call check_leap_year_present_step(start_ymd, start_tod, nstep, dtime)

   s_yr  = start_ymd/10000
   s_mn  = (start_ymd-s_yr*10000)/100
   s_dy  = start_ymd-s_yr*10000-s_mn*100
   s_sc  = start_tod

   do i = 1, s_mn-1
      s_dy = s_dy + month(i)
   end do 

   s_time = real(s_yr*allday+s_dy+s_sc/86400._r8,r8)

   curr_time = curr_mod_time + s_time

   nsec = int((curr_time - floor(curr_time))*86400.)
   if(nsec .lt. 0)call endrun('Error in set_date_from_time_float')

   yr = floor(curr_time)/allday
   calday = mod(floor(curr_time), allday)

   mn = 1
   do i = 1, 12
       if(calday .gt. month(i))then
           calday = calday - month(i)
           mn = mn+1
       else
           exit
       end if
   end do

   day = calday

   end subroutine set_date_from_time_float1
   
   subroutine set_date_from_time_float(start_ymd, start_tod, curr_mod_time, yr, mn, day, nsec)
   !io 
   integer,  intent(in)      :: start_ymd, start_tod
   real(r8), intent(in)      :: curr_mod_time
   integer,  intent(out)     :: yr, mn, day, nsec
   ! local
   integer                   :: calday, i
   integer                   :: s_yr, s_mn, s_dy, s_sc
   real(r8)                  :: curr_time, s_time

   s_yr  = start_ymd/10000
   s_mn  = (start_ymd-s_yr*10000)/100
   s_dy  = start_ymd-s_yr*10000-s_mn*100
   s_sc  = start_tod

   do i = 1, s_mn-1
      s_dy = s_dy + month(i)
   end do

   s_time = real(s_yr*365._r8+s_dy+s_sc/86400._r8,r8)

   curr_time = curr_mod_time + s_time

   nsec = int((curr_time - floor(curr_time))*86400.)
   if(nsec .lt. 0)call endrun('Error in set_date_from_time_float')

   yr = floor(curr_time)/365
   calday = mod(floor(curr_time), 365)

   mn = 1
   do i = 1, 12
       if(calday .gt. month(i))then
           calday = calday - month(i)
           mn = mn+1
       else
           exit
       end if
   end do

   day = calday

   end subroutine set_date_from_time_float

   subroutine intToChar5(sc,sec)
    integer(i4), intent(in) :: sc
    character(len=5), intent(out) :: sec 
    character(len=5) :: tmp
    
       if(sc.lt.10)then
          write(tmp,'(i1)') sc
          sec=trim('0000')//trim(tmp)
       else if(sc.ge.10 .and. sc.lt.100)then
          write(tmp,'(i2)') sc
          sec=trim('000')//trim(tmp)
       else if(sc.ge.100.and. sc.lt.1000)then
          write(tmp,'(i3)') sc
          sec=trim('00')//trim(tmp)
       else if(sc.ge.1000.and. sc.lt.10000)then
          write(tmp,'(i4)') sc
          sec=trim('0')//trim(tmp)
       else
          write(tmp,'(i5)') sc
          sec=trim(tmp)
       end if
       return
   end subroutine intToChar5

   subroutine reset_month_allday(yr)
    integer(i4) :: yr

     if((mod(yr,4).eq.0.and.mod(yr,100).ne.0).or.(mod(yr,400).eq.0))then
         month(1:12)  = leap_month
         allday       = leap_allday
     else
         month(1:12)  = norm_month
         allday       = norm_allday
     end if
     return
   end subroutine reset_month_allday

end module grist_time_manager
