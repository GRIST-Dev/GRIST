 module grist_clocks
   use grist_mpi, only: grist_comm

   implicit none
   include 'mpif.h'

   private

   type :: clock
     character(len=32)  :: name
     integer(8)         :: tick
     integer(8)         :: total_ticks
     integer            :: nums
     logical            :: sync_on_begin !whether sync PEs before timing
   end type clock

   integer, parameter   :: MAX_CLOCKS=200
   integer              :: clock_num = 0
   integer              :: clock0 = 0
   type(clock),save     :: clocks(MAX_CLOCKS)

   integer(8)           :: ticks_per_sec
   real(8)              :: tick_rate

   interface clock_begin
       procedure :: clock_begin1
       procedure :: clock_begin2
   end interface clock_begin
   interface clock_end
       procedure :: clock_end1
       procedure :: clock_end2
   end interface clock_end

   public               :: clocks, clock_num, clock0, clock_id, &
                           clock_begin, clock_end, clock_summary
 contains

   function mpi_rank(comm) result(res)
     implicit none
     integer, optional :: comm
     integer :: ierr
     integer :: res 
     
     if(present(comm)) then
       call MPI_COMM_RANK(Comm, res, ierr)
     else
       call MPI_COMM_RANK(grist_comm(), res, ierr)
     end if
    
   end function

   function mpi_size(comm) result(res)
     implicit none
     integer, optional :: comm
     integer :: ierr
     integer :: res 
     
     if(present(comm)) then
        call MPI_COMM_SIZE(Comm, res, ierr)
     else
        call MPI_COMM_SIZE(grist_comm(), res, ierr)
     end if
     
   end function

   subroutine clock_init( id, name, flags)
     integer,           intent(in) :: id
     character(len=*),  intent(in) :: name
     integer, intent(in), optional :: flags

     clocks(id)%name = name
     clocks(id)%tick = 0
     clocks(id)%total_ticks = 0
     clocks(id)%nums = 0
     clocks(id)%sync_on_begin = .FALSE.
     if( PRESENT(flags) .and. flags .eq. 1 )then
        clocks(id)%sync_on_begin = .TRUE.
     end if

     return
   end subroutine clock_init

   !return an ID for a new or existing clock
   function clock_id(name, flags, comm)
     integer                       :: clock_id
     character(len=*),  intent(in) :: name
     integer, intent(in), optional :: flags
     integer, intent(in), optional :: comm
     integer                       :: comm_
     integer                       :: flags_

     clock_id = 1
     comm_ = grist_comm()
     if(present(comm)) comm_ = comm
     flags_ = 0
     if(present(flags)) flags_ = flags

     if (clock_num .eq. 0) then  !first
        clock_num = clock_id
        call clock_init(clock_id,name,flags_)
     else
        FIND_CLOCK: do while( trim(name) .ne. trim(clocks(clock_id)%name) )
           clock_id = clock_id + 1
           if (clock_id .gt. clock_num) then
              if (clock_id .gt. MAX_CLOCKS) then
                 if(mpi_rank(comm_) .eq. 0) print*, 'WARNING: CLOCK_ID: too many clock requests, this one is ignored.'
              else               !new clock: initialize
                 clock_num = clock_id
                 call clock_init(clock_id,name,flags_)
                 exit FIND_CLOCK
              end if
           end if
        end do FIND_CLOCK
     endif

     return
   end function clock_id

   subroutine clock_begin1(id, comm)
     integer, intent(in) :: id
     integer, optional   :: comm
     integer             :: ierr
     integer             :: comm_
     integer             :: errcode

     comm_ = grist_comm()
     if(present(comm)) comm_ = comm
     if (id .eq. 0) return
     if (id .lt. 0 .or. id .gt. clock_num) then
        if(mpi_rank(comm_) .eq. 0) print*,"FATAL: CLOCK_BEGIN: invalid id."
        call mpi_abort(comm_, errcode, ierr)
     end if

     if( clocks(id)%sync_on_begin )then
        !do an untimed sync at the beginning of the clock
        !to measure load imbalance for this
        call MPI_Barrier(comm_, ierr)
     end if
     call SYSTEM_CLOCK( clocks(id)%tick )
     clocks(id)%nums = clocks(id)%nums + 1
     return
   end subroutine clock_begin1

   subroutine clock_begin2(name, flags, comm)
     character(len=*),  intent(in) :: name
     integer, intent(in), optional :: flags
     integer, intent(in), optional :: comm
     integer                       :: flags_
     integer                       :: comm_
     integer                       :: id
     integer                       :: ierr
     integer                       :: errcode

     comm_ = grist_comm()
     if(present(comm)) comm_ = comm
     flags_ = 0
     if(present(flags)) flags_ = flags
     id = clock_id(name, flags_, comm_)

     if (id .eq. 0) return
     if (id .lt. 0 .or. id .gt. clock_num) then
        if(mpi_rank(comm_) .eq. 0) print*,"FATAL: CLOCK_BEGIN: invalid id."
        call mpi_abort(comm_, errcode, ierr)
     end if

     if( clocks(id)%sync_on_begin )then
        !do an untimed sync at the beginning of the clock
        !to measure load imbalance for this
        call MPI_Barrier(comm_, ierr)
     end if
     call SYSTEM_CLOCK( clocks(id)%tick )
     clocks(id)%nums = clocks(id)%nums + 1
     return
   end subroutine clock_begin2

   subroutine clock_end1(id, comm)
     integer, intent(in) :: id
     integer, optional   :: comm
     integer(8)          :: delta, end_tick
     integer             :: comm_
     integer             :: ierr
     integer             :: errcode

     comm_ = grist_comm()
     if(present(comm)) comm_ = comm

     if (id .eq. 0) return
     if (id .lt. 0 .or. id .gt. clock_num) then
        if(mpi_rank(comm_) .eq. 0) print*,"FATAL: CLOCK_BEGIN: invalid id."
        call mpi_abort(comm_, errcode, ierr)
     end if

     call SYSTEM_CLOCK(end_tick)
     delta = end_tick - clocks(id)%tick
     if (delta .lt. 0) then
        write(*,* )'pe, id, start_tick, end_tick, delta, max_ticks=', mpi_rank(comm_), id, clocks(id)%tick, end_tick, delta
        if(mpi_rank(comm_) .eq. 0) print*, 'WARNING: CLOCK_END: Clock rollover, assumed single roll.' 
     end if
     clocks(id)%total_ticks = clocks(id)%total_ticks + delta

     return
   end subroutine clock_end1

   subroutine clock_end2(name, flags, comm)
     character(len=*),  intent(in) :: name
     integer, intent(in), optional :: flags
     integer, intent(in), optional :: comm
     integer                       :: flags_
     integer                       :: comm_
     integer                       :: id
     integer                       :: ierr
     integer                       :: errcode
     integer(8)                    :: delta, end_tick

     comm_ = grist_comm()
     if(present(comm)) comm_ = comm
     flags_ = 0
     if(present(flags)) flags_ = flags
     id = clock_id(name, flags_, comm_)

     if (id .eq. 0) return
     if (id .lt. 0 .or. id .gt. clock_num) then
        if(mpi_rank(comm_) .eq. 0) print*,"FATAL: CLOCK_BEGIN: invalid id."
        call mpi_abort(comm_, errcode, ierr)
     end if

     call SYSTEM_CLOCK(end_tick)
     delta = end_tick - clocks(id)%tick
     if (delta .lt. 0) then
        write(*,* )'pe, id, start_tick, end_tick, delta, max_ticks=', mpi_rank(comm_), id, clocks(id)%tick, end_tick, delta
        if(mpi_rank(comm_) .eq. 0) print*, 'WARNING: CLOCK_END: Clock rollover, assumed single roll.' 
     end if
     clocks(id)%total_ticks = clocks(id)%total_ticks + delta

     return
   end subroutine clock_end2

   subroutine clock_summary(comm)
     integer, optional   :: comm
     integer             :: i, ierr
     real(8)             :: t, tmin, tmax, tavg, tstd, t_tmp, t_total
     integer             :: comm_

     comm_ = grist_comm()
     if(present(comm)) comm_ = comm

     if (clock_num .gt. 0) then
        if (mpi_rank(comm_) .eq. 0) then
           write(*,'(/a,i6,a)') 'Tabulating clock statistics across ', mpi_size(comm_), ' PEs...'
           write(*, '(/32x,a)') '   calls          tmin          tmax          tavg          tstd         tfrac'
        end if
        call MPI_Barrier(comm_, ierr)
        call SYSTEM_CLOCK(count_rate=ticks_per_sec)
        tick_rate = 1.d0/ticks_per_sec
        if (clock0 .ne. 0) then
           t_total = clocks(clock0)%total_ticks*tick_rate
        else
           t_total = 0.D0
        end if
        do i = 1, clock_num
           !times between mpp_clock ticks
           t = clocks(i)%total_ticks*tick_rate
           call MPI_Reduce(t, tmin, 1, MPI_REAL8, MPI_MIN, 0, comm_, ierr)
           call MPI_Reduce(t, tmax, 1, MPI_REAL8, MPI_MAX, 0, comm_, ierr)
           call MPI_Allreduce(t, tavg, 1, MPI_REAL8, MPI_SUM, comm_, ierr)
           tavg  = tavg/mpi_size(comm_)
           t_tmp = (t-tavg)**2
           call MPI_Reduce(t_tmp, tstd, 1, MPI_REAL8, MPI_SUM, 0, comm_, ierr)
           if (mpi_rank(comm_) .eq. 0) then
              tstd = sqrt( tstd/mpi_size(comm_) )
              if (clock0 .ne. 0) then
                 write(*,'(a32,i8,5f14.6)') & 
                   clocks(i)%name, clocks(i)%nums, tmin, tmax, tavg, tstd, tavg/t_total
              else
                 write(*,'(a32,i8,5f14.6)') & 
                   clocks(i)%name, clocks(i)%nums, tmin, tmax, tavg, tstd, 0.0!tavg/t_total
              end if
           end if
        end do
     end if

   end subroutine clock_summary
end module grist_clocks