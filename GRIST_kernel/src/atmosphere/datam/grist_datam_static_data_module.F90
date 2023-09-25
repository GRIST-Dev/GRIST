 
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: Data Module for GRIST
!              currently contain surface static data
!              Name Convention:
!              static_data_${varname}_at_$(location)
! Revision history:
! (1) SST&SIC data from a sst-sic.nc
! (2) other land surface static data from a static.nc
! (3) SST/SIC daily-interpolation and topo-smooth also done here before they use
!----------------------------------------------------------------------------

  module grist_datam_static_data_module
!
! infrastructure
!
    use grist_domain_types, only: global_domain
    use grist_data_types,   only: scalar_1d_field, scalar_2d_field, exchange_field_list_1d, &
                                  wrap_allocate_data1d, wrap_deallocate_data1d, wrap_deallocate_data2d
    use grist_constants,    only: rearth, i4, i8, r8, pi, one, zero, t00, r4, gravity
    use grist_nml_module,   only: nlev, staticFilePath, sftopoFilePath, sstFilePath, aqua_sst_style, real_sst_style, numMonSST, test_real_case, topo_type, &
                                  start_ymd,start_tod, model_timestep, smooth_topo, smooth_type, nsmooth_topo, firstClipSic, doAquaPlanet, ko4_coef,topo_ko_coef
    use grist_fileio_list_1d_module_par, only: wrap_read_1d_group_rst
    use grist_fileio_list_2d_module_par, only: wrap_read_2d_group, wrap_read_2d_group_rst !---cheyz add
    use grist_config_partition,          only: exchange_data_1d_add, exchange_data_1d
    use grist_mpi
    use grist_datam_initial_data_module, only: initialData_soilh_at_pc_surface
    use grist_dycore_primal_flux_operators_2d, only: calc_primal_normal_flux_at_edge_1d

    implicit none

    private
    public :: grist_static_data_construct,&
              grist_static_data_destruct, &
              grist_static_fill_current_sst, &        ! real sst (sic) for ocean points
              grist_static_data_generate_analytic_sst ! full analytic sst

!================================================
! All data at primal cell
!================================================

! SST&SIC
   type(scalar_1d_field), public   :: staticData_sst_at_pc_surface      ! data directly used by other model components
   type(scalar_1d_field), public   :: staticData_sic_at_pc_surface      ! ...
   type(scalar_2d_field), public   :: staticData_sst_mon_at_pc_surface  ! assumed "monthly-mean mid-point" from sst data
   type(scalar_2d_field), public   :: staticData_sic_mon_at_pc_surface  ! assumed "monthly-mean mid-point" from sic data
! land-surface
   type(scalar_1d_field), public   :: staticData_phis_at_pc_surface     ! surface geopotential
   type(scalar_1d_field), public   :: staticData_phis2_at_pc_surface    ! for read topo from sftopo_file
   type(scalar_1d_field), public   :: staticData_landfrac_at_pc_surface ! land mask from 0->1
!--cheyz
   type(scalar_1d_field), public   :: staticData_luIndex_at_pc_surface
   type(scalar_1d_field), public   :: staticData_annual_deep_soilTemp_at_pc_surface
   type(scalar_1d_field), public   :: staticData_soilTypetop_at_pc_surface
   type(scalar_2d_field), public   :: staticData_albedo_at_pc_surface
   type(scalar_2d_field), public   :: staticData_greenfrac_at_pc_surface
   type(scalar_1d_field), public   :: staticData_snoalb_at_pc_surface
! GWDO-related
   type(scalar_1d_field), public   :: staticData_var2d_at_pc_surface
   type(scalar_1d_field), public   :: staticData_oc12d_at_pc_surface
   type(scalar_1d_field), public   :: staticData_oa2d1_at_pc_surface
   type(scalar_1d_field), public   :: staticData_oa2d2_at_pc_surface
   type(scalar_1d_field), public   :: staticData_oa2d3_at_pc_surface
   type(scalar_1d_field), public   :: staticData_oa2d4_at_pc_surface
   type(scalar_1d_field), public   :: staticData_ol2d1_at_pc_surface
   type(scalar_1d_field), public   :: staticData_ol2d2_at_pc_surface
   type(scalar_1d_field), public   :: staticData_ol2d3_at_pc_surface
   type(scalar_1d_field), public   :: staticData_ol2d4_at_pc_surface

   real(r4), parameter :: sstFillValue = -1e30 ! as in datafile, float
   real(r4), parameter :: sicFillValue = -1e30 ! as in datafile, float

  CONTAINS

    subroutine grist_static_data_construct(mesh)
! io
      type(global_domain), intent(inout), target :: mesh
! local
      integer(i4)  ::  it, ie, iv, inb,icell
      real(r8), allocatable     :: tmp_sum(:)
      integer(i8)  :: dim_2

!================================================
! primal cell, full level
!================================================

     call wrap_allocate_data1d(mesh%nv, staticData_phis_at_pc_surface)
     call wrap_allocate_data1d(mesh%nv, staticData_phis2_at_pc_surface)
     call wrap_allocate_data1d(mesh%nv, staticData_landfrac_at_pc_surface)
     call wrap_allocate_data1d(mesh%nv, staticData_sst_at_pc_surface)
     call wrap_allocate_data1d(mesh%nv, staticData_sic_at_pc_surface)

     if(numMonSST.eq.0)then
       if(mpi_rank().eq.0) print*,"Please set numMonSST in nml's data_para, model stops"
       call mpi_abort()
     end if

     if(.not.allocated(staticData_sst_mon_at_pc_surface%f))  allocate(staticData_sst_mon_at_pc_surface%f(mesh%nv,numMonSST),source=zero) ; staticData_sst_mon_at_pc_surface%pos = 0
     if(.not.allocated(staticData_sic_mon_at_pc_surface%f))  allocate(staticData_sic_mon_at_pc_surface%f(mesh%nv,numMonSST),source=zero) ; staticData_sic_mon_at_pc_surface%pos = 0

!--cheyz add  20200424
#ifdef USE_NOAHMP
     call wrap_allocate_data1d(mesh%nv, staticData_luIndex_at_pc_surface)
     call wrap_allocate_data1d(mesh%nv, staticData_annual_deep_soilTemp_at_pc_surface)
     call wrap_allocate_data1d(mesh%nv, staticData_soilTypetop_at_pc_surface)
     call wrap_allocate_data1d(mesh%nv, staticData_snoalb_at_pc_surface)
     if(.not.allocated(staticData_albedo_at_pc_surface%f)) allocate(staticData_albedo_at_pc_surface%f(12, mesh%nv),source=zero); staticData_albedo_at_pc_surface%pos   = 0
     if(.not.allocated(staticData_greenfrac_at_pc_surface%f)) allocate(staticData_greenfrac_at_pc_surface%f(12, mesh%nv),source=zero); staticData_greenfrac_at_pc_surface%pos   = 0
#endif

! GWDO-related, only for real-case
    if(.not.doAquaPlanet.and.test_real_case)then
       call wrap_allocate_data1d(mesh%nv, staticData_var2d_at_pc_surface)
       call wrap_allocate_data1d(mesh%nv, staticData_oc12d_at_pc_surface)
       call wrap_allocate_data1d(mesh%nv, staticData_oa2d1_at_pc_surface)
       call wrap_allocate_data1d(mesh%nv, staticData_oa2d2_at_pc_surface)
       call wrap_allocate_data1d(mesh%nv, staticData_oa2d3_at_pc_surface)
       call wrap_allocate_data1d(mesh%nv, staticData_oa2d4_at_pc_surface)
       call wrap_allocate_data1d(mesh%nv, staticData_ol2d1_at_pc_surface)
       call wrap_allocate_data1d(mesh%nv, staticData_ol2d2_at_pc_surface)
       call wrap_allocate_data1d(mesh%nv, staticData_ol2d3_at_pc_surface)
       call wrap_allocate_data1d(mesh%nv, staticData_ol2d4_at_pc_surface)
    end if

!------------
! SST/SIC
!------------

! all-points are first set to analytic sst, as in AquaPlanet
     call grist_static_data_generate_analytic_sst(mesh)
! realize ocean points with realistic SST/SIC if test_real_case is true
     if(test_real_case) call grist_static_fill_current_sst(mesh,0,mesh%nv_full,model_timestep)

!------------
! Static
!------------

#ifdef USE_NOAHMP
! 20210909, yizhang: seperate control surface topo file using sftopoFilePath
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(staticFilePath), '', 'PHIS',     0, staticData_phis_at_pc_surface%f)
     if(topo_type.eq.'sftopo')then ! unit uses gpm 
     ! means any 3rd-party (other than static and init) topo data for flexible testing
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(sftopoFilePath), '', 'PHIS',     0, staticData_phis2_at_pc_surface%f)
     end if
!---cheyz 20200424
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(staticFilePath), '', 'MASK',     0, staticData_landfrac_at_pc_surface%f)
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(staticFilePath), '', 'LU_INDEX', 0, staticData_luIndex_at_pc_surface%f)
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(staticFilePath), '', 'SOILTEMP', 0, staticData_annual_deep_soilTemp_at_pc_surface%f)
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(staticFilePath), '', 'SOILTOP' , 0, staticData_soilTypetop_at_pc_surface%f)
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(staticFilePath), '', 'SNOALB'  , 0, staticData_snoalb_at_pc_surface%f)
     dim_2  = 12
     call wrap_read_2d_group_rst(mesh%gcomm_read, trim(staticFilePath), '', 'ALBEDO'   , dim_2,  0, staticData_albedo_at_pc_surface%f)
     call wrap_read_2d_group_rst(mesh%gcomm_read, trim(staticFilePath), '', 'GREENFRAC', dim_2,  0, staticData_greenfrac_at_pc_surface%f)
#else
! old-name from NCAR-topo software
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(staticFilePath), '', 'PHIS',     0, staticData_phis_at_pc_surface%f)
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(staticFilePath), '', 'LANDFRAC' ,0, staticData_landfrac_at_pc_surface%f)
#endif
     staticData_landfrac_at_pc_surface%f = min(1._r8,staticData_landfrac_at_pc_surface%f)

!-----------------------------------------------------------------
! topo type
! for this, initial data must be allocated before static data
!-----------------------------------------------------------------

     select case(trim(topo_type))
     case('sftopo')
         where(staticData_landfrac_at_pc_surface%f(1:mesh%nv).gt.0.5_r8) staticData_phis_at_pc_surface%f(1:mesh%nv) = staticData_phis2_at_pc_surface%f(1:mesh%nv)
     case('init')
! set all land points to those from initial file
         where(staticData_landfrac_at_pc_surface%f(1:mesh%nv).gt.0.5_r8) staticData_phis_at_pc_surface%f(1:mesh%nv) = gravity*initialData_soilh_at_pc_surface%f(1:mesh%nv)
     case default
       !do nothing
     end select

     if(mpi_rank().eq.0) print*, "compute length is:", mesh%nv_compute, "halo1   length is:", mesh%nv_halo(1), "halo2   length is:", mesh%nv_halo(2), &
                                 "full    length is:", mesh%nv_full   , "dynamic length is:", mesh%nv
!---------------------
! topo smooth if any
!---------------------
    
    if(smooth_topo)then
       staticData_phis2_at_pc_surface%f(1:mesh%nv) = staticData_phis_at_pc_surface%f(1:mesh%nv)
       call grist_inline_topo_smooth(mesh,staticData_phis2_at_pc_surface)
       ! only repalce for land-topo
       where(staticData_landfrac_at_pc_surface%f(1:mesh%nv).gt.0.5_r8) staticData_phis_at_pc_surface%f(1:mesh%nv) = staticData_phis2_at_pc_surface%f(1:mesh%nv)
    end if

!--------------
! Static GWDO
!--------------
    if(.not.doAquaPlanet.and.test_real_case)then
       call wrap_read_1d_group_rst(mesh%gcomm_read,trim(staticFilePath),'','VAR2D',0,staticData_var2d_at_pc_surface%f)
       call wrap_read_1d_group_rst(mesh%gcomm_read,trim(staticFilePath),'','CON'  ,0,staticData_oc12d_at_pc_surface%f)
       call wrap_read_1d_group_rst(mesh%gcomm_read,trim(staticFilePath),'','OA1'  ,0,staticData_oa2d1_at_pc_surface%f)
       call wrap_read_1d_group_rst(mesh%gcomm_read,trim(staticFilePath),'','OA2'  ,0,staticData_oa2d2_at_pc_surface%f)
       call wrap_read_1d_group_rst(mesh%gcomm_read,trim(staticFilePath),'','OA3'  ,0,staticData_oa2d3_at_pc_surface%f)
       call wrap_read_1d_group_rst(mesh%gcomm_read,trim(staticFilePath),'','OA4'  ,0,staticData_oa2d4_at_pc_surface%f)
       call wrap_read_1d_group_rst(mesh%gcomm_read,trim(staticFilePath),'','OL1'  ,0,staticData_ol2d1_at_pc_surface%f)
       call wrap_read_1d_group_rst(mesh%gcomm_read,trim(staticFilePath),'','OL2'  ,0,staticData_ol2d2_at_pc_surface%f)
       call wrap_read_1d_group_rst(mesh%gcomm_read,trim(staticFilePath),'','OL3'  ,0,staticData_ol2d3_at_pc_surface%f)
       call wrap_read_1d_group_rst(mesh%gcomm_read,trim(staticFilePath),'','OL4'  ,0,staticData_ol2d4_at_pc_surface%f)
    end if

    return
  end subroutine grist_static_data_construct

  subroutine grist_static_data_destruct

! SST&SIC
      call wrap_deallocate_data1d(staticData_sst_at_pc_surface)
      call wrap_deallocate_data1d(staticData_sic_at_pc_surface)
      call wrap_deallocate_data2d(staticData_sst_mon_at_pc_surface)
      call wrap_deallocate_data2d(staticData_sic_mon_at_pc_surface)
! Land-surface static
      call wrap_deallocate_data1d(staticData_phis_at_pc_surface)
      call wrap_deallocate_data1d(staticData_phis2_at_pc_surface)
      call wrap_deallocate_data1d(staticData_landfrac_at_pc_surface)
      call wrap_deallocate_data1d(staticData_luIndex_at_pc_surface)
      call wrap_deallocate_data1d(staticData_annual_deep_soilTemp_at_pc_surface)
      call wrap_deallocate_data1d(staticData_soilTypetop_at_pc_surface)
      call wrap_deallocate_data1d(staticData_snoalb_at_pc_surface)
      call wrap_deallocate_data2d(staticData_albedo_at_pc_surface)
      call wrap_deallocate_data2d(staticData_greenfrac_at_pc_surface)
! static gwdo
      call wrap_deallocate_data1d(staticData_var2d_at_pc_surface)
      call wrap_deallocate_data1d(staticData_oc12d_at_pc_surface)
      call wrap_deallocate_data1d(staticData_oa2d1_at_pc_surface)
      call wrap_deallocate_data1d(staticData_oa2d2_at_pc_surface)
      call wrap_deallocate_data1d(staticData_oa2d3_at_pc_surface)
      call wrap_deallocate_data1d(staticData_oa2d4_at_pc_surface)
      call wrap_deallocate_data1d(staticData_ol2d1_at_pc_surface)
      call wrap_deallocate_data1d(staticData_ol2d2_at_pc_surface)
      call wrap_deallocate_data1d(staticData_ol2d3_at_pc_surface)
      call wrap_deallocate_data1d(staticData_ol2d4_at_pc_surface)

      return
   end subroutine grist_static_data_destruct

   subroutine grist_static_fill_current_sst(mesh,nstep,ncell,dtime,is_first_restart)
    use grist_time_manager,              only: get_curr_date
    use grist_nml_module,                only: sstFileNameHead, sstFileNameTail, sstFile_year_beg
!io
     type(global_domain), intent(in), target :: mesh
     integer(i4),  intent(in) :: nstep
     integer(i4),  intent(in) :: ncell
     real(r8),     intent(in) :: dtime
     logical, optional,      intent(in) :: is_first_restart
! local
     integer(i4)         :: year,mon,day,sec,mon_of_data
     character(len=1024) :: sstFileName
     character(len=4)    :: char_year, char_mon, char_day
     integer(i8)         :: dim_2, icell, it

     call get_curr_date(start_ymd, start_tod, nstep, dtime, year, mon, day, sec)

     !IF(CYCLE_SST)THEN ! 12-mon SST
    select case(trim(real_sst_style))
    case('CYCLE')

     mon_of_data = mon

     dim_2  = numMonSST ! must be 12
     if(numMonSST.ne.12)then
        print*,"for cycle sst, numMonSST must be 12"
        call mpi_abort()
     end if

     sstFileName = trim(sstFileNameHead)//'cycle'//trim(sstFileNameTail)
!
! only get sst data only for ONCE : 2d sst (12-mons)
!
     if(nstep.eq.0)then
        call wrap_read_2d_group(mesh%gcomm_read, trim(sstFilePath),trim(sstFileName), 'sst', dim_2, 0, staticData_sst_mon_at_pc_surface%f)
        call wrap_read_2d_group(mesh%gcomm_read, trim(sstFilePath),trim(sstFileName), 'sic', dim_2, 0, staticData_sic_mon_at_pc_surface%f)
        if(firstClipSic)then ! default no
           do icell = 1, ncell
              do it = 1, dim_2
                 staticData_sic_mon_at_pc_surface%f(icell,it) = max(min(staticData_sic_mon_at_pc_surface%f(icell,it),one),zero) 
              end do
           end do
        end if
     end if
!
! cycle sst here with multi-year runs
!
     where(staticData_sst_mon_at_pc_surface%f(:,mon_of_data).ne.sstFillValue) staticData_sst_at_pc_surface%f = staticData_sst_mon_at_pc_surface%f(:,mon_of_data)
     where(staticData_sic_mon_at_pc_surface%f(:,mon_of_data).ne.sicFillValue) staticData_sic_at_pc_surface%f = staticData_sic_mon_at_pc_surface%f(:,mon_of_data)
!
! linearly interpolate daily value in time
! 
     if(mon.ne.1.and.mon.ne.12)then
        call get_daily_sstsic_from_mon(year,mon,day,ncell, &
                                    staticData_sst_mon_at_pc_surface%f(1:ncell,mon_of_data)  ,&
                                    staticData_sst_mon_at_pc_surface%f(1:ncell,mon_of_data+1),&
                                    staticData_sst_mon_at_pc_surface%f(1:ncell,mon_of_data-1),&
                                    staticData_sst_at_pc_surface%f(1:ncell),&
                                    staticData_sic_mon_at_pc_surface%f(1:ncell,mon_of_data)  ,&
                                    staticData_sic_mon_at_pc_surface%f(1:ncell,mon_of_data+1),&
                                    staticData_sic_mon_at_pc_surface%f(1:ncell,mon_of_data-1),&
                                    staticData_sic_at_pc_surface%f(1:ncell))
     end if
     if(mon.eq.1)then
        call get_daily_sstsic_from_mon(year,mon,day,ncell, &
                                    staticData_sst_mon_at_pc_surface%f(1:ncell,mon_of_data)  ,&
                                    staticData_sst_mon_at_pc_surface%f(1:ncell,mon_of_data+1),&
                                    staticData_sst_mon_at_pc_surface%f(1:ncell,12)           ,&
                                    staticData_sst_at_pc_surface%f(1:ncell), &
                                    staticData_sic_mon_at_pc_surface%f(1:ncell,mon_of_data)  ,&
                                    staticData_sic_mon_at_pc_surface%f(1:ncell,mon_of_data+1),&
                                    staticData_sic_mon_at_pc_surface%f(1:ncell,12)           ,&
                                    staticData_sic_at_pc_surface%f(1:ncell) )
     end if
     if(mon.eq.12)then
        call get_daily_sstsic_from_mon(year,mon,day,ncell, &
                                    staticData_sst_mon_at_pc_surface%f(1:ncell,mon_of_data)  ,&
                                    staticData_sst_mon_at_pc_surface%f(1:ncell,1)            ,&
                                    staticData_sst_mon_at_pc_surface%f(1:ncell,mon_of_data-1),&
                                    staticData_sst_at_pc_surface%f(1:ncell),&
                                    staticData_sic_mon_at_pc_surface%f(1:ncell,mon_of_data)  ,&
                                    staticData_sic_mon_at_pc_surface%f(1:ncell,1)            ,&
                                    staticData_sic_mon_at_pc_surface%f(1:ncell,mon_of_data-1),&
                                    staticData_sic_at_pc_surface%f(1:ncell))
     end if

    case('AMIP')

!==================================================================================================
! AMIP SST
! we store AMIP-SST into a-decadal+2 mons (122-mon) data, e.g., 1970-1979, 1980-1989, 1990-1999, ...
! just the volume of ~4*30*nlev T data, not so big; with the 1st and last record for the last dec 
! and next jan; so only read data at the initial time or when a year reaches a 1970,1980-like
! number (mod(year,10)is 0)
!==================================================================================================

     if(nstep.eq.0.or.mod(year,10).eq.0.or.is_first_restart)then

        dim_2  = numMonSST  ! must be 122
        if(numMonSST.ne.122)then
           print*, "for AMIP (non-cycle) sst, numMonSST in nml must be 122=1+12*10+1"
           call mpi_abort()
        end if

        if(nstep.eq.0.or.is_first_restart)then
           write(char_year,'(i4)') sstFile_year_beg ! given in namelist
        else
           write(char_year,'(i4)') year
        end if
        sstFileName = trim(sstFileNameHead)//trim(char_year)//trim(sstFileNameTail)
        call wrap_read_2d_group(mesh%gcomm_read,trim(sstFilePath),trim(sstFileName),'sst', dim_2, 0, staticData_sst_mon_at_pc_surface%f)
        call wrap_read_2d_group(mesh%gcomm_read,trim(sstFilePath),trim(sstFileName),'sic', dim_2, 0, staticData_sic_mon_at_pc_surface%f)

        if(firstClipSic)then ! default no
           do icell = 1, ncell
              do it = 1, dim_2
                 staticData_sic_mon_at_pc_surface%f(icell,it) = max(min(staticData_sic_mon_at_pc_surface%f(icell,it),one),zero)
              end do
           end do
        end if

     end if

! 1mon+ 19X0-19X9 + 1mon
     mon_of_data = 1+mon+mod(year,10)*12

! raw monthly SST
     where(staticData_sst_mon_at_pc_surface%f(:,mon_of_data).ne.sstFillValue) staticData_sst_at_pc_surface%f = staticData_sst_mon_at_pc_surface%f(:,mon_of_data)
     where(staticData_sic_mon_at_pc_surface%f(:,mon_of_data).ne.sicFillValue) staticData_sic_at_pc_surface%f = staticData_sic_mon_at_pc_surface%f(:,mon_of_data)

     call get_daily_sstsic_from_mon(year,mon,day,ncell, &
                                    staticData_sst_mon_at_pc_surface%f(1:ncell,mon_of_data)  ,&
                                    staticData_sst_mon_at_pc_surface%f(1:ncell,mon_of_data+1),&
                                    staticData_sst_mon_at_pc_surface%f(1:ncell,mon_of_data-1),&
                                    staticData_sst_at_pc_surface%f    (1:ncell), &
                                    staticData_sic_mon_at_pc_surface%f(1:ncell,mon_of_data)  ,&
                                    staticData_sic_mon_at_pc_surface%f(1:ncell,mon_of_data+1),&
                                    staticData_sic_mon_at_pc_surface%f(1:ncell,mon_of_data-1),&
                                    staticData_sic_at_pc_surface%f    (1:ncell) )
     !where(staticData_sic_mon_at_pc_surface%f(:,mon_of_data).ne.sicFillValue) staticData_sic_at_pc_surface%f = staticData_sic_mon_at_pc_surface%f(:,mon_of_data)

    case('DAILY')
!
! in this case, sst/sic is one-day one value, so just read&fill based on current date,
! when call this grist_static_fill_current_sst routine
! same for initial and restart run because current data is always now
!
      dim_2  = numMonSST   !  must be 1
      if(numMonSST.ne.1)then
         print*, "for daily sst, numMonSST must be 1"
         call mpi_abort()
      end if

      write(char_year,'(i4)') year

      if(mon.ge.10)then
         write(char_mon,'(i2)') mon
      else
         write(char_mon,'(i1)') mon
         char_mon="0"//char_mon
      end if
      if(day.ge.10)then
         write(char_day,'(i2)') day
      else
         write(char_day,'(i1)') day
         char_day="0"//char_day
      end if

      sstFileName = trim(sstFileNameHead)//'daily'//trim(char_year)//trim(char_mon)//trim(char_day)//trim(sstFileNameTail)
      
      call wrap_read_2d_group(mesh%gcomm_read,trim(sstFilePath),trim(sstFileName),'sst', dim_2, 0, staticData_sst_mon_at_pc_surface%f)
      call wrap_read_2d_group(mesh%gcomm_read,trim(sstFilePath),trim(sstFileName),'sic', dim_2, 0, staticData_sic_mon_at_pc_surface%f)
      ! assign data
      where(staticData_sst_mon_at_pc_surface%f(:,numMonSST).ne.sstFillValue) staticData_sst_at_pc_surface%f = staticData_sst_mon_at_pc_surface%f(:,1)
      where(staticData_sic_mon_at_pc_surface%f(:,numMonSST).ne.sicFillValue) staticData_sic_at_pc_surface%f = staticData_sic_mon_at_pc_surface%f(:,1)
      mon_of_data = 1
    case default
       if(mpi_rank().eq.0) print*,"you must select a real_sst_style: 'CYCLE', 'AMIP', 'DAILY', model aborts"
  end select

     if(mpi_rank().eq.0) print*,"SST&SIC: current month is", mon, "mon_of_sst_data is", mon_of_data, "real_sst_style: ", trim(real_sst_style)
     return
   end subroutine grist_static_fill_current_sst

!--------------------------------------------------------------
! Generate sea surface temperature based on analytic function
! rather than data; seaice is zero; For Aqua Planet Run
!--------------------------------------------------------------

   subroutine grist_static_data_generate_analytic_sst(mesh)
! io
     type(global_domain), intent(in), target :: mesh
! local
     integer(i4) :: iv
     call wrap_allocate_data1d(mesh%nv_full, staticData_sst_at_pc_surface)
     call wrap_allocate_data1d(mesh%nv_full, staticData_sic_at_pc_surface)
     do iv = 1, mesh%nv_full
        staticData_sst_at_pc_surface%f(iv) = nh2000_aqua_planet_sst(real(mesh%vtx_lat(iv),r8),trim(aqua_sst_style))
        staticData_sic_at_pc_surface%f(iv) = zero
     end do
     return
   end subroutine grist_static_data_generate_analytic_sst

!-----------------------------------------------------
!  sst function following Neale and Hoskins (2000)
!-----------------------------------------------------

  real(r8) function nh2000_aqua_planet_sst(lat,flag)
     real(r8),         intent(in) :: lat
     character(len=*), intent(in) :: flag

     select case(flag)
       case('RJ')
           nh2000_aqua_planet_sst = 29._r8
       case('NK_CNTL')
           if(lat.gt.-pi/3._r8 .and. lat.lt.pi/3._r8)then
              nh2000_aqua_planet_sst = 27._r8*(one-(sin(1.5_r8*lat)**2))
           else
              nh2000_aqua_planet_sst = zero
           end if
       case('NK_FLAT')
           if(lat.gt.-pi/3._r8 .and. lat.lt.pi/3._r8)then
              nh2000_aqua_planet_sst = 27._r8*(one-(sin(1.5_r8*lat)**4))
           else
              nh2000_aqua_planet_sst = zero
           end if
       case('NK_QOBS')
           if(lat.gt.-pi/3._r8 .and. lat.lt.pi/3._r8)then
              nh2000_aqua_planet_sst = 13.5_r8*(2._r8-sin(1.5_r8*lat)**4-sin(1.5_r8*lat)**2)
           else
              nh2000_aqua_planet_sst = zero
           end if
       case default !  'RJ' 
           nh2000_aqua_planet_sst = 29._r8
     end select
     nh2000_aqua_planet_sst = nh2000_aqua_planet_sst + t00
     return
  end function nh2000_aqua_planet_sst

  subroutine get_daily_sstsic_from_mon(year,mon,day,ncell, &
                                       sst_iThis1,sst_iNext1,sst_iLast1,sst_day, &
                                       sic_iThis1,sic_iNext1,sic_iLast1,sic_day)

     integer(i4),      intent(in) :: year, mon, day  ! current date
     integer(i4),      intent(in) :: ncell
     real(r8),         intent(in) :: sst_iThis1(ncell)  ! this mon
     real(r8),         intent(in) :: sst_iNext1(ncell)  ! next mon
     real(r8),         intent(in) :: sst_iLast1(ncell)  ! last mon
     real(r8),         intent(out):: sst_day(ncell)     ! sst at target day
     real(r8),         intent(in) :: sic_iThis1(ncell)  ! this mon
     real(r8),         intent(in) :: sic_iNext1(ncell)  ! next mon
     real(r8),         intent(in) :: sic_iLast1(ncell)  ! last mon
     real(r8),         intent(out):: sic_day(ncell)     ! sic at target day
! local
     integer(i4)  :: d1, d2, icell
     integer(i4)  :: medInMon(14)=(/16, 16,15,16,15,16,15,16,16,15,16,15,16, 16/) ! no-leap, Appendix 3 of Taylor
    !integer(i4)  :: medInMon(14)=(/16, 16,15,16,16,16,16,16,16,16,16,16,16, 16/) ! no-leap
     integer(i4)  :: endInMon(14)=(/31, 31,28,31,30,31,30,31,31,30,31,30,31, 31/) ! Dec,1-12,Jan, 365 days
     character(len=20) :: intp_scheme

#ifdef USE_LEAP_YEAR
     if((mod(year,4).eq.0.and.mod(year,100).ne.0).or.(mod(year,400).eq.0))then
        endInMon=(/31,31,29,31,30,31,30,31,31,30,31,30,31,31/)
     else
        endInMon=(/31,31,28,31,30,31,30,31,31,30,31,30,31,31/)
     end if
#endif

     intp_scheme = 'raw'
!
! monthly-mean of this interpolated daily sst has at most <0.5 K difference
! (at certain areas 07-09, midhigh-lat) from the raw monthly midpoint sst
! if model monthly mean sst is compared against raw unprocessed monthly mean value (tos), nearly identical!!
!
     select case(trim(intp_scheme))
     case('raw')
!
! |med(lastMon)__________|end(lastMon)_______________|med(thisMon)..............
! if day is a middle-month day, just use sst of this mon as in the raw data
! if day lies in the second-half of this mon, interpolate using thisMon and nextMon sst
! if day lies in the first-half  of this mon, interpolate using lastMon and thisMon sst
!
     if(day.eq.medInMon(mon+1))then
        do icell = 1, ncell
           if(sst_iThis1(icell).ne.sstFillValue)then
              sst_day(icell) = sst_iThis1(icell)
              sst_day(icell) = max(sst_iThis1(icell),271.38_r8)
           end if
           if(sic_iThis1(icell).ne.sicFillValue)then
              sic_day(icell) = max(min(sic_iThis1(icell),one),zero)
           end if
        end do
      !where(sst_iThis1.ne.sstFillValue) sst_day = sst_iThis1 ! for entire day
     else if(day.gt.medInMon(mon+1).and.day.le.endInMon(mon+1))then ! use iThis1 and iNext1 as rear and forward
        d1 = day            -medInMon(mon+1) ! no +1
        d2 = endInMon(mon+1)-medInMon(mon+1)+medInMon(mon+2)
        do icell = 1, ncell
           if(sst_iNext1(icell).ne.sstFillValue.and.sst_iThis1(icell).ne.sstFillValue)then
              sst_day(icell) = (d1*sst_iNext1(icell)+(d2-d1)*sst_iThis1(icell))/d2
              sst_day(icell) = max(sst_day(icell),271.38_r8)
           end if
           if(sic_iNext1(icell).ne.sicFillValue.and.sic_iThis1(icell).ne.sicFillValue)then
              sic_day(icell) = (d1*sic_iNext1(icell)+(d2-d1)*sic_iThis1(icell))/d2
              sic_day(icell) = max(min(sic_day(icell),one),zero)
           end if
        end do
     else if(day.lt.medInMon(mon+1).and.day.ge.1)then ! use iLast1 and iThis1 as rear and forward
        d1 = endInMon(mon)  -medInMon(mon)+day
        d2 = medInMon(mon+1)+endInMon(mon)-medInMon(mon)
        do icell = 1, ncell
           if(sst_iThis1(icell).ne.sstFillValue.and.sst_iLast1(icell).ne.sstFillValue)then
              sst_day(icell) = (d1*sst_iThis1(icell)+(d2-d1)*sst_iLast1(icell))/d2
              sst_day(icell) = max(sst_day(icell),271.38_r8)
           end if
           if(sic_iThis1(icell).ne.sicFillValue.and.sic_iLast1(icell).ne.sicFillValue)then
              sic_day(icell) = (d1*sic_iThis1(icell)+(d2-d1)*sic_iLast1(icell))/d2
              sic_day(icell) = max(min(sic_day(icell),one),zero)
           end if
        end do
     else
        if(mpi_rank().eq.0) print*,"get_daily_sstsic_from_monMedian: ill condition, abort"
        call mpi_abort
     end if

    case default
        if(mpi_rank().eq.0) print*,"get_daily_sstsic_from_monMedian: no intp_scheme, abort"
        call mpi_abort
    end select

    return
  end subroutine get_daily_sstsic_from_mon

  subroutine grist_inline_topo_smooth(mesh, static_phis_local)
! io
    type(global_domain),    intent(inout) :: mesh
    type(scalar_1d_field),  intent(inout) :: static_phis_local
! local
    real(r8)          :: tmp_sum(mesh%nv_compute), area_sum(mesh%nv_compute)
    real(r8)          :: gradient_at_prime_edge(mesh%ne_full), tend_laplacian_4th_at_pc(mesh%nv_full)
    real(r8)          :: gradient_at_prime_edge_2(2,mesh%ne_full) ! store left and right gradients at two sides
    real(r8)          :: fake_wind(nlev,mesh%ne_full) ! store left and right gradients at two sides
    integer(i4)       :: itime, iv, inb, icell, v1, v2, ie
    real(r8)          :: v1v2(3), flag, div_sum,flag1
    !character(len=20) :: smooth_type
    type(exchange_field_list_1d),pointer :: field_head_1d
    type(scalar_1d_field) :: scalar_edge_value
    real(r8)          :: max_phis, min_phis


    field_head_1d =>null()

    !smooth_type='CellAvg'

    select case(trim(smooth_type))
    case('hyper4th')

    ! fourth-order hyperdiffusion
      DO itime = 1, nsmooth_topo ! how much pass
!
! 1st pass: gradient at edge, counter edge's normal direction
!
        do ie = 1, mesh%ne_halo(1)
           v1      = mesh%edt_v(1,ie)
           v2      = mesh%edt_v(2,ie)
           v1v2    = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
           flag    = sign(one,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))
           gradient_at_prime_edge(ie) = flag*(static_phis_local%f(v2)-static_phis_local%f(v1))/(rearth*mesh%edt_leng(ie))
        end do
!
! 1st pass: divergence at cell
!
        do iv = 1, mesh%nv_halo(1)
           tend_laplacian_4th_at_pc(iv) = zero
           div_sum = zero
           do inb = 1, mesh%vtx_nnb(iv)
              ie       = mesh%vtx_ed(inb,iv)
              div_sum  = div_sum+gradient_at_prime_edge(ie)*mesh%plg_nr(inb,iv)*mesh%edp_leng(ie)
           end do
           tend_laplacian_4th_at_pc(iv) = div_sum/(rearth*mesh%plg_areag(iv))
        end do
!
! 2nd pass: gradient at edge, counter edge's normal direction, untill ne_compute
!
        do ie = 1, mesh%ne_compute
           v1      = mesh%edt_v(1,ie)
           v2      = mesh%edt_v(2,ie)
           v1v2    = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
           flag    = sign(one,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))
           gradient_at_prime_edge(ie) = flag*(tend_laplacian_4th_at_pc(v2)-tend_laplacian_4th_at_pc(v1))/(rearth*mesh%edt_leng(ie))*topo_ko_coef*mesh%edt_scale_4th(ie)
        end do
!
! 2nd pass: divergence at cell
!
        do iv = 1, mesh%nv_compute
           tend_laplacian_4th_at_pc(iv) = zero
           div_sum = zero
           do inb = 1, mesh%vtx_nnb(iv)
              ie       = mesh%vtx_ed(inb,iv)
              div_sum  = div_sum+gradient_at_prime_edge(ie)*mesh%plg_nr(inb,iv)*mesh%edp_leng(ie)
           end do
           tend_laplacian_4th_at_pc(iv) = div_sum/(rearth*mesh%plg_areag(iv))
        end do

        static_phis_local%f(1:mesh%nv_compute) = static_phis_local%f(1:mesh%nv_compute)-tend_laplacian_4th_at_pc(1:mesh%nv_compute)*model_timestep
    ! renew halo region
        call exchange_data_1d_add(mesh,field_head_1d,static_phis_local)
        call exchange_data_1d(mesh%local_block,field_head_1d)

    END DO

    case('Laplace')
    ! normal Laplace

      DO itime = 1, nsmooth_topo ! how much pass
!
! 1st pass: gradient at edge, counter edge's normal direction
!
        do ie = 1, mesh%ne_halo(1)
           v1      = mesh%edt_v(1,ie)
           v2      = mesh%edt_v(2,ie)
           v1v2    = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
           flag    = sign(one,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))
           gradient_at_prime_edge(ie) = flag*(static_phis_local%f(v2)-static_phis_local%f(v1))/(rearth*mesh%edt_leng(ie))
           gradient_at_prime_edge(ie) = 0.5_r8*(mesh%vtxCellLeng(v1)+mesh%vtxCellLeng(v2))*rearth*2._r8*gradient_at_prime_edge(ie)
        end do
!
! 1st pass: divergence at cell
!
        do iv = 1, mesh%nv_halo(1)
           tend_laplacian_4th_at_pc(iv) = zero  ! actually 2nd but use 4th as name
           div_sum = zero
           do inb = 1, mesh%vtx_nnb(iv)
              ie       = mesh%vtx_ed(inb,iv)
              div_sum  = div_sum+gradient_at_prime_edge(ie)*mesh%plg_nr(inb,iv)*(mesh%edp_leng(ie))
           end do
           tend_laplacian_4th_at_pc(iv) = div_sum/(rearth*mesh%plg_areag(iv))
        end do

        static_phis_local%f(1:mesh%nv_compute) = static_phis_local%f(1:mesh%nv_compute)+model_timestep*tend_laplacian_4th_at_pc(1:mesh%nv_compute)
    ! renew halo region
        call exchange_data_1d_add(mesh,field_head_1d,static_phis_local)
        call exchange_data_1d(mesh%local_block,field_head_1d)

    END DO

    case('cellAvg')

    DO itime = 1, nsmooth_topo ! how much pass
       do iv = 1, mesh%nv_compute ! until compute layer
          tmp_sum(iv)  = static_phis_local%f(iv)*mesh%plg_areag(iv)
          area_sum(iv) = mesh%plg_areag(iv)
! only smooth topo higher than 10-m
          if(static_phis_local%f(iv).gt.gravity*10._r8)then
             do inb = 1, mesh%vtx_nnb(iv)
                icell       = mesh%vtx_nb(inb,iv)
                tmp_sum(iv) = tmp_sum(iv) +static_phis_local%f(icell)*mesh%plg_areag(icell)
                area_sum(iv)= area_sum(iv)+mesh%plg_areag(icell)
             end do
             tmp_sum(iv) = tmp_sum(iv)/area_sum(iv)
          end if
       end do
       static_phis_local%f(1:mesh%nv_compute) = tmp_sum(1:mesh%nv_compute)
! renew halo region
       call exchange_data_1d_add(mesh,field_head_1d,static_phis_local)
       call exchange_data_1d(mesh%local_block,field_head_1d)
    END DO

    case ('hiLaplace')
    ! Laplace with a high-order gradient

! Laplace operator with a 4th-order gradient operator
    fake_wind = 1._r8

    if(.not.allocated(scalar_edge_value%f)) allocate(scalar_edge_value%f(mesh%ne_full))
    scalar_edge_value%f    = zero
    scalar_edge_value%pos  = 6

    DO itime = 1, nsmooth_topo ! how much pass

    ! alter dynamica length temporarily
       mesh%ne = mesh%ne_compute
       call calc_primal_normal_flux_at_edge_1d(mesh, &
                                               fake_wind, & ! O4 never mind this
                                               fake_wind, & ! O4 never mind this
                                               static_phis_local%f , & ! raw data
                                               scalar_edge_value%f             , & ! at this time, not gradient but interpolated state
                                               4) ! never mind this

       call exchange_data_1d_add(mesh,field_head_1d,scalar_edge_value)
       call exchange_data_1d(mesh%local_block,field_head_1d)

       mesh%ne = mesh%ne_full
!
! 1st pass: gradient at edge, counter edge's normal direction
!
! V1----E-----V2
! V1->V2 is assumed positive direction
!
        do ie = 1, mesh%ne_halo(1)
           v1      = mesh%edt_v(1,ie)
           v2      = mesh%edt_v(2,ie)
           v1v2    = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
           flag    = sign(one,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))
! for v1->E
           flag1   = sign(one,scalar_edge_value%f(ie)-static_phis_local%f(v1))
           gradient_at_prime_edge_2(1,ie) = flag*flag1*abs(scalar_edge_value%f(ie)-static_phis_local%f(v1))/(rearth*mesh%edt_leng(ie)/2._r8)
! for E->v2
           flag1   = sign(one,static_phis_local%f(v2)-scalar_edge_value%f(ie))
           gradient_at_prime_edge_2(2,ie) = flag*flag1*abs(static_phis_local%f(v2)-scalar_edge_value%f(ie))/(rearth*mesh%edt_leng(ie)/2._r8)
           gradient_at_prime_edge_2(1:2,ie) = 0.5_r8*(mesh%vtxCellLeng(v1)+mesh%vtxCellLeng(v2))*rearth*2._r8*gradient_at_prime_edge_2(1:2,ie)
        end do
!
! 1st pass: divergence at cell
!
        Do iv = 1, mesh%nv_halo(1)
           tend_laplacian_4th_at_pc(iv) = zero  ! actually 2nd but use 4th as name
           div_sum = zero
           do inb = 1, mesh%vtx_nnb(iv)
              ie       = mesh%vtx_ed(inb,iv)
              v1       = mesh%edt_v(1,ie)
              v2       = mesh%edt_v(2,ie)
              if(iv.eq.v1)then
                 div_sum  = div_sum+gradient_at_prime_edge_2(1,ie)*mesh%plg_nr(inb,iv)*mesh%edp_leng(ie)
              elseif(iv.eq.v2)then
                 div_sum  = div_sum+gradient_at_prime_edge_2(2,ie)*mesh%plg_nr(inb,iv)*mesh%edp_leng(ie)
              else
                 if(mpi_rank().eq.0) print*,"neither v1 nor v2 conforms to iv, bad condition, model aborts"
                 call mpi_abort
              end if
           end do
           tend_laplacian_4th_at_pc(iv) = div_sum/(rearth*mesh%plg_areag(iv))
        End Do

        static_phis_local%f(1:mesh%nv_compute) = static_phis_local%f(1:mesh%nv_compute)+model_timestep*tend_laplacian_4th_at_pc(1:mesh%nv_compute)
    ! renew halo region for whole field
        call exchange_data_1d_add(mesh,field_head_1d,static_phis_local)
        call exchange_data_1d(mesh%local_block,field_head_1d)

    END DO

    if(allocated(scalar_edge_value%f)) deallocate(scalar_edge_value%f)

    case default
       if(mpi_rank().eq.0) print*,"no smooth_type is set when you activate smooth_topo"
       call mpi_abort
    end select
       if(mpi_rank().eq.0) print*,"topo smooth_type is:", trim(smooth_type)

    return
  end subroutine grist_inline_topo_smooth

  end module grist_datam_static_data_module
