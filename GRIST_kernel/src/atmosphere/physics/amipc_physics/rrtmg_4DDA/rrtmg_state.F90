!===================================================================
!
!  Created by LiXiaohan on 19/8/19.
!  Manages the absorber concentrations in the layers RRTMG operates
!  including an extra layer over the model if needed.
!
!===================================================================

 module rrtmg_state

    use grist_constants,                    only: i4, r8
    use grist_nml_module,                   only: nlev, nlevp
    use grist_handle_error,                 only: endrun
    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam
    use grist_mpi

    implicit none
    private
    save

    public :: rrtmg_state_init,     &
              rrtmg_state_create,   &
              rrtmg_state_update,   &
              rrtmg_state_t,        &
              rrtmg_state_destroy,  &
              num_rrtmg_levs

    type rrtmg_state_t

       real(r8), allocatable :: h2ovmr(:,:)      ! h2o volume mixing ratio
       real(r8), allocatable :: o3vmr(:,:)       ! o3 volume mixing ratio
       real(r8), allocatable :: co2vmr(:,:)      ! co2 volume mixing ratio 
       real(r8), allocatable :: ch4vmr(:,:)      ! ch4 volume mixing ratio 
       real(r8), allocatable :: o2vmr(:,:)       ! o2  volume mixing ratio 
       real(r8), allocatable :: n2ovmr(:,:)      ! n2o volume mixing ratio 
       real(r8), allocatable :: cfc11vmr(:,:)    ! cfc11 volume mixing ratio
       real(r8), allocatable :: cfc12vmr(:,:)    ! cfc12 volume mixing ratio
       real(r8), allocatable :: cfc22vmr(:,:)    ! cfc22 volume mixing ratio
       real(r8), allocatable :: ccl4vmr(:,:)     ! ccl4 volume mixing ratio

       real(r8), allocatable :: pmidmb(:,:)      ! Level pressure (hPa)
       real(r8), allocatable :: pintmb(:,:)      ! Model interface pressure (hPa)
       real(r8), allocatable :: tlay(:,:)        ! mid point temperature
       real(r8), allocatable :: tlev(:,:)        ! interface temperature

    end type rrtmg_state_t

    integer :: num_rrtmg_levs ! number of pressure levels greate than 1.e-4_r8 mbar

    real(r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
    real(r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
    real(r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
    real(r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
    real(r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
    real(r8), parameter :: amdo2 = 0.905140_r8   ! Molecular weight of dry air / oxygen
    real(r8), parameter :: amdc1 = 0.210852_r8   ! Molecular weight of dry air / CFC11
    real(r8), parameter :: amdc2 = 0.239546_r8   ! Molecular weight of dry air / CFC12

    integer :: o3_idx, co2_idx, ch4_idx, o2_idx, cfc11_idx, cfc12_idx, n2o_idx
 contains

! Purpose : sets the number of model levels RRTMG operates
    subroutine rrtmg_state_init
    use grist_constants,                    only: p00
    use grist_hpe_constants,                only: eta_face
 
    implicit none

    !local
    real(r8) :: pref_face(nlevp)
    integer  :: i

    pref_face(1:nlevp) = eta_face(1:nlevp)*p00 
    ! The following cuts off RRTMG at roughly the point where it becomes
    ! invalid due to low pressure.
    num_rrtmg_levs = count( pref_face(:) > 1._r8 ) ! pascals (1.e-2 mbar)

    !find ghg index, LiXH
    co2_idx = -1; o3_idx  = -1; ch4_idx=-1; o2_idx=-1; cfc11_idx=-1; cfc12_idx=-1; n2o_idx=-1
    do i = 1, pstate_cam%total_ghg_num
        if(trim(pstate_cam%ghg_at_pc_full_level(i)%name) .eq. 'CO2')co2_idx = i
        if(trim(pstate_cam%ghg_at_pc_full_level(i)%name) .eq. 'O2')o2_idx = i
        if(trim(pstate_cam%ghg_at_pc_full_level(i)%name) .eq. 'O3')o3_idx = i
        if(trim(pstate_cam%ghg_at_pc_full_level(i)%name) .eq. 'CH4')ch4_idx = i
        if(trim(pstate_cam%ghg_at_pc_full_level(i)%name) .eq. 'N2O')n2o_idx = i
        if(trim(pstate_cam%ghg_at_pc_full_level(i)%name) .eq. 'CFC11')cfc11_idx = i
        if(trim(pstate_cam%ghg_at_pc_full_level(i)%name) .eq. 'CFC12')cfc12_idx = i
    end do
    if(co2_idx .lt. 0 .or. o3_idx .lt. 0 .or. o2_idx .lt. 0  .or. ch4_idx .lt. 0 .or. &
      n2o_idx .lt. 0 .or. cfc11_idx .lt. 0 .or. cfc12_idx .lt. 0)then
        if(mpi_rank()==0)print*,'can not find index in ghg data!'
        call endrun('rrtmg_state_init')
    end if
 
    end subroutine rrtmg_state_init
 

    function rrtmg_state_create( ncol ) result( rstate )
    use grist_constants,              only: stebol

    implicit none

    integer, intent(in)           :: ncol
    type(rrtmg_state_t), pointer  :: rstate

    real(r8) :: dy                   ! Temporary layer pressure thickness
    real(r8) :: tint(nlevp, ncol)    ! Model interface temperature
    integer  :: i, kk, k

    allocate( rstate )

    allocate( rstate%h2ovmr(num_rrtmg_levs,ncol) )
    allocate( rstate%o3vmr(num_rrtmg_levs,ncol) )
    allocate( rstate%co2vmr(num_rrtmg_levs,ncol) )
    allocate( rstate%ch4vmr(num_rrtmg_levs,ncol) )
    allocate( rstate%o2vmr(num_rrtmg_levs,ncol) )
    allocate( rstate%n2ovmr(num_rrtmg_levs,ncol) )
    allocate( rstate%cfc11vmr(num_rrtmg_levs,ncol) )
    allocate( rstate%cfc12vmr(num_rrtmg_levs,ncol) )
    allocate( rstate%cfc22vmr(num_rrtmg_levs,ncol) )
    allocate( rstate%ccl4vmr(num_rrtmg_levs,ncol) )

    allocate( rstate%pmidmb(num_rrtmg_levs,ncol) )
    allocate( rstate%pintmb(num_rrtmg_levs+1,ncol) )
    allocate( rstate%tlay(num_rrtmg_levs,ncol) )
    allocate( rstate%tlev(num_rrtmg_levs+1,ncol) )

    ! Calculate interface temperatures (following method
    ! used in radtpl for the longwave), using surface upward flux and
    ! stebol constant in mks units
    do i = 1, ncol
       tint(1,i)     = pstate%temp_at_pc_full_level%f(1,i)
       tint(nlevp,i) = sqrt(sqrt(pstate%atm_in_lwup_at_pc_surface%f(i)/stebol))

       do k = 2, nlev
          dy = (log(pstate%pressure_at_pc_face_level%f(k,i))   - log(pstate%pressure_at_pc_full_level%f(k,i))) &
              /(log(pstate%pressure_at_pc_full_level%f(k-1,i)) - log(pstate%pressure_at_pc_full_level%f(k,i)))
          tint(k,i) = pstate%temp_at_pc_full_level%f(k,i)     &
                     -dy * (pstate%temp_at_pc_full_level%f(k,i) - pstate%temp_at_pc_full_level%f(k-1,i))
       end do
    end do

    do k = 1, num_rrtmg_levs
       kk = max(k + (nlevp-num_rrtmg_levs)-1,1)

       rstate%pmidmb(k,:ncol) = pstate%pressure_at_pc_full_level%f(kk,:ncol) * 1.e-2_r8
       rstate%pintmb(k,:ncol) = pstate%pressure_at_pc_face_level%f(kk,:ncol) * 1.e-2_r8

       rstate%tlay(k,:ncol) = pstate%temp_at_pc_full_level%f(kk,:ncol)
       rstate%tlev(k,:ncol) = tint(kk,:ncol)
    enddo

    ! bottom interface
    rstate%pintmb(num_rrtmg_levs+1,:ncol) = pstate%pressure_at_pc_face_level%f(nlevp,:ncol) * 1.e-2_r8 ! mbar
    rstate%tlev(num_rrtmg_levs+1,:ncol)   = tint(nlevp,:ncol)

    ! top layer thickness
    if (num_rrtmg_levs==nlevp) then
       rstate%pmidmb(1,:ncol) = 0.5_r8 * rstate%pintmb(2,:ncol) 
       rstate%pintmb(1,:ncol) = 1.e-4_r8 ! mbar
    endif

    endfunction rrtmg_state_create


! Purpose : updates the concentration fields
    subroutine rrtmg_state_update(ncol,rstate)
    implicit none

    integer,             intent(in) :: ncol
    type(rrtmg_state_t), pointer    :: rstate

    real(r8), dimension(:,:), allocatable :: sp_hum ! specific humidity
    real(r8), dimension(:,:), allocatable :: n2o    ! nitrous oxide mass mixing ratio
    real(r8), dimension(:,:), allocatable :: ch4    ! methane mass mixing ratio
    real(r8), dimension(:,:), allocatable :: o2     ! O2 mass mixing ratio
    real(r8), dimension(:,:), allocatable :: cfc11  ! cfc11 mass mixing ratio
    real(r8), dimension(:,:), allocatable :: cfc12  ! cfc12 mass mixing ratio
    real(r8), dimension(:,:), allocatable :: o3     ! Ozone mass mixing ratio
    real(r8), dimension(:,:), allocatable :: co2    ! co2   mass mixing ratio
    
    integer  ::  i, kk, k

    allocate(sp_hum(nlev,ncol));sp_hum=0._r8    !'H2O'
    allocate(o2(nlev,ncol));o2=0._r8            !'O2'
    allocate(o3(nlev,ncol));o3=0._r8            !'O3'
    allocate(co2(nlev,ncol));co2=0._r8          !'CO2'
    allocate(n2o(nlev,ncol));n2o=0._r8          !'N2O'
    allocate(ch4(nlev,ncol));ch4=0._r8          !'CH4'
    allocate(cfc11(nlev,ncol));cfc11=0._r8      !'CFC11'
    allocate(cfc12(nlev,ncol));cfc12=0._r8      !'CFC12'

    sp_hum(:,1:ncol) = pstate%tracer_mxrt_at_pc_full_level%f(1,:,1:ncol)
    o2(:,1:ncol)     = pstate_cam%ghg_at_pc_full_level(o2_idx)%f(:,1:ncol)
    o3(:,1:ncol)     = pstate_cam%ghg_at_pc_full_level(o3_idx)%f(:,1:ncol)
    co2(:,1:ncol)    = pstate_cam%ghg_at_pc_full_level(co2_idx)%f(:,1:ncol)
    n2o(:,1:ncol)    = pstate_cam%ghg_at_pc_full_level(n2o_idx)%f(:,1:ncol)
    ch4(:,1:ncol)    = pstate_cam%ghg_at_pc_full_level(ch4_idx)%f(:,1:ncol)
    cfc11(:,1:ncol)  = pstate_cam%ghg_at_pc_full_level(cfc11_idx)%f(:,1:ncol)
    cfc12(:,1:ncol)  = pstate_cam%ghg_at_pc_full_level(cfc12_idx)%f(:,1:ncol)

    do k = 1, num_rrtmg_levs

       kk = max(k + (nlevp-num_rrtmg_levs)-1,1)

       rstate%ch4vmr(k,:ncol)   = ch4(kk,:ncol) * amdm
       rstate%h2ovmr(k,:ncol)   = (sp_hum(kk,:ncol) / (1._r8 - sp_hum(kk,:ncol))) * amdw
       rstate%o3vmr(k,:ncol)    = o3(kk,:ncol) * amdo
       rstate%co2vmr(k,:ncol)   = co2(kk,:ncol) * amdc
       rstate%o2vmr(k,:ncol)    = o2(kk,:ncol) * amdo2
       rstate%n2ovmr(k,:ncol)   = n2o(kk,:ncol) * amdn
       rstate%cfc11vmr(k,:ncol) = cfc11(kk,:ncol) * amdc1
       rstate%cfc12vmr(k,:ncol) = cfc12(kk,:ncol) * amdc2
       rstate%cfc22vmr(k,:ncol) = 0._r8
       rstate%ccl4vmr(k,:ncol)  = 0._r8

    enddo

    deallocate(sp_hum)
    deallocate(o2)
    deallocate(o3)
    deallocate(co2)
    deallocate(n2o)
    deallocate(ch4)
    deallocate(cfc11)
    deallocate(cfc12)

!-------------LiXH Test--------------
!if(mpi_rank()==0)then
!print*,'start of test:'
!print*,'o2',rstate%o2vmr
!print*,'o3',rstate%o3vmr
!print*,'co2',rstate%co2vmr
!print*,'n2o',rstate%n2ovmr
!print*,'ch4',rstate%ch4vmr
!print*,'cfc11',rstate%cfc11vmr
!print*,'cfc12',rstate%cfc12vmr
!print*,'end of test:'
!end if
!-------------LiXH Test--------------
 
    end subroutine rrtmg_state_update


!--------------------------------------------------------------------------------
! de-allocates an rrtmg_state object
!--------------------------------------------------------------------------------
  subroutine rrtmg_state_destroy(rstate)

    implicit none

    type(rrtmg_state_t), pointer   :: rstate

    deallocate(rstate%h2ovmr)
    deallocate(rstate%o3vmr)
    deallocate(rstate%co2vmr)
    deallocate(rstate%ch4vmr)
    deallocate(rstate%o2vmr)
    deallocate(rstate%n2ovmr)
    deallocate(rstate%cfc11vmr)
    deallocate(rstate%cfc12vmr)
    deallocate(rstate%cfc22vmr)
    deallocate(rstate%ccl4vmr)

    deallocate(rstate%pmidmb)
    deallocate(rstate%pintmb)
    deallocate(rstate%tlay)
    deallocate(rstate%tlev)

    deallocate( rstate )
    nullify(rstate)

  end subroutine rrtmg_state_destroy



 end module rrtmg_state
