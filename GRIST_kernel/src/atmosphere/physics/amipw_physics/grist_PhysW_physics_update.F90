
!--------------------------------------------------------------------------------------
! Created on 2020
! Author: Yi Zhang
! Version 1.0
! Description: Control how physics tendencies are organized for dynamics
! Revision history:
!     (1) P3&P2 are in active use and testing, others are for reference and document
!     (2) Filling logic:
!         Fill the raw tendency required in Ptend_f2 to Ptend_f3
!         Fill the raw tendency required in Ptend_rk to Ptend_f1
!         All these raw tendencies will be re-converted and re-organized
!         by DTP physics-dynamics coupling workflow and dynamics integration; whether
!         update the state (e.g. P2) or not (e.g., P3) only modify internal phys_state
!     (3) The coding can be flexible; but which should be used depend on our guiding 
!         principles and modeling experience; because f3&f1 are used to hold
!         phys tendencies, they must be initialized with nspecies size, even
!         thou some space may not be used
!--------------------------------------------------------------------------------------

 module grist_PhysW_physics_update

    use grist_constants,       only: i4, r8, zero
    use grist_nml_module,      only: nlev, physics_coupling
    use grist_physics_data_structure,      only: pstate, ptend_f1, ptend_f3
    use grist_PhysW_data_structure,        only: pstate_wrf, ptend_wrf, p_qv, p_qc, p_qr, p_qi, p_qs, p_qg, p_nr, p_ni, p_ns, p_ng
    use grist_mpi
 
    implicit none
    private
    save

    public  :: PhysW_update

contains
!
! for this update, the nlev dim directly comes from nml, which denotes
! the number of full level rather than face level as in WRF
!

  subroutine PhysW_update(dtime, ncell, process_name, scheme_name,first_called)

! io
    real(r8),              intent(in)    :: dtime
    integer(i4),           intent(in)    :: ncell
    character(len=*),      intent(in)    :: process_name
    character(len=*),      intent(in)    :: scheme_name
    logical,  optional,    intent(in)    :: first_called
! local
!
! zero tend, f3 only contains f2 tend, f1/rk contains all heating tend
!
     if(present(first_called).and.first_called)then
!
! f3, when converted, become f2
!
        ptend_f3%tend_u_wind_at_pc_full_level%f         = zero
        ptend_f3%tend_v_wind_at_pc_full_level%f         = zero
        ptend_f3%tend_potential_temp_at_pc_full_level%f = zero
        ptend_f3%tend_tracer_mxrt_at_pc_full_level%f    = zero
!
! f1, though it is a PDC operator by dynamical model, we temperatily use
! it to store the tendency that will be converted to ptend_RK, f1 is retained
! for tracer (latest 2021, f1 has be almost discarded as a compiling option)
!
        ptend_f1%tend_u_wind_at_pc_full_level%f         = zero
        ptend_f1%tend_v_wind_at_pc_full_level%f         = zero
        ptend_f1%tend_potential_temp_at_pc_full_level%f = zero
        ptend_f1%tend_tracer_mxrt_at_pc_full_level%f    = zero
     end if
!
! fill f3 with tendencies that want to be in f2
! fill f1 with tendencies that want to be in rk
!
! now, only P2 and P3 are supported by wrf physics
! other are retained for reference
!
    select case(trim(physics_coupling))
    case('P1')
       call phys_update_p1(dtime,  ncell, process_name, scheme_name)
    case('P1A')
       call phys_update_p1a(dtime, ncell, process_name, scheme_name)
    case('P2')
       call phys_update_p2(dtime,  ncell, process_name, scheme_name)
    case('P3')
       call phys_update_p3(dtime,  ncell, process_name, scheme_name)
    case('P3A')
       call phys_update_p3a(dtime, ncell, process_name, scheme_name)
    case('P5')
       call phys_update_p5(dtime,  ncell, process_name, scheme_name)
    case('P6')
       call phys_update_p6(dtime,  ncell, process_name, scheme_name)
    case default
       if(mpi_rank().eq.0) print*, "you must set physics_couplingL due to use of WRF physics"
       call mpi_abort()
    end select

    return
  end subroutine PhysW_update

!------------------------------------------------------------
! CU and RAD's uvpt tends are used as RK physics
! MP and PBL's uvpt tends are used as F2 operator spliting
! all q tends are F2 operator splitting
!------------------------------------------------------------

  subroutine phys_update_p1(dtime, ncell, process_name, scheme_name)
    real(r8),              intent(in)   :: dtime
    integer(i4),           intent(in)   :: ncell
    character(len=*),      intent(in)   :: process_name
    character(len=*),      intent(in)   :: scheme_name

    select case(trim(process_name))
    case('CU')
       
       if(trim(scheme_name).eq.'TDK')then

          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruucuten(1:ncell,1:nlev,1))

          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvcuten(1:ncell,1:nlev,1))

          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqccuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqicuten(1:ncell,1:nlev,1))

       end if

!================================================
! MP
! Kessler can only be run alone
! needs comment other processes
!================================================

    case('MP')
       if(trim(scheme_name).eq.'KES')then
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))
       end if
       if(trim(scheme_name).eq.'LIN')then
! store tendency
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(5 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(5 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqsmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(6 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(6 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqgmpten(1:ncell,1:nlev,1))
       end if

!================================================
! PBL
!================================================

    case('PBL')
       if(trim(scheme_name).eq.'YSU')then

          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruublten(1:ncell,1:nlev,1))

          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvblten(1:ncell,1:nlev,1))

          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqiblten(1:ncell,1:nlev,1))

       end if
    case('RAD')
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthraten(1:ncell,1:nlev,1))
    case default
       if(mpi_rank().eq.0) print*,"no process name is selected, please double check"
       call mpi_abort()
    end select
    return
  end subroutine phys_update_p1

!------------------------------------------------------------
! CU and RAD's uvpt tends are used as RK physics
! MP and PBL's uvpt tends are used as F2 operator spliting
! and PBL also updates state variables inside physics
! all q tends are F2 operator splitting
!------------------------------------------------------------

  subroutine phys_update_p1a(dtime, ncell, process_name, scheme_name)
    real(r8),              intent(in)   :: dtime
    integer(i4),           intent(in)   :: ncell
    character(len=*),      intent(in)   :: process_name
    character(len=*),      intent(in)   :: scheme_name

    select case(trim(process_name))
    case('CU')
       
       if(trim(scheme_name).eq.'TDK')then

          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruucuten(1:ncell,1:nlev,1))

          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvcuten(1:ncell,1:nlev,1))

          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqccuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqicuten(1:ncell,1:nlev,1))

       end if

!================================================
! MP
! Kessler can only be run alone
! needs comment other processes
!================================================

    case('MP')
       if(trim(scheme_name).eq.'KES')then
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))
       end if
       if(trim(scheme_name).eq.'LIN')then
! store tendency
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(5 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(5 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqsmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(6 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(6 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqgmpten(1:ncell,1:nlev,1))
!
! renew temperature
!
         ! pstate_wrf%t_phy(  1:ncell,1:nlev,1)      = pstate_wrf%th_phy( 1:ncell,1:nlev,1)*pstate_wrf%pi_phy( 1:ncell,1:nlev,1)
       end if

!================================================
! PBL, also update physics state in an OS style
!================================================

    case('PBL')

       if(trim(scheme_name).eq.'YSU')then

          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruublten(1:ncell,1:nlev,1))

          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvblten(1:ncell,1:nlev,1))

          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqiblten(1:ncell,1:nlev,1))

!
! update state here, as the original pbl does not update them
!
          pstate_wrf%u_phy(  1:ncell,1:nlev,1)      = pstate_wrf%u_phy(  1:ncell,1:nlev,1)+dtime*&
                                                          ptend_wrf%ruublten(1:ncell,1:nlev,1)
          pstate_wrf%v_phy(  1:ncell,1:nlev,1)      = pstate_wrf%v_phy(  1:ncell,1:nlev,1)+dtime*&
                                                        ptend_wrf%rvvblten(1:ncell,1:nlev,1)
          pstate_wrf%th_phy( 1:ncell,1:nlev,1)      = pstate_wrf%th_phy( 1:ncell,1:nlev,1)+dtime*&
                                                        ptend_wrf%rthblten(1:ncell,1:nlev,1)
          pstate_wrf%moist(  1:ncell,1:nlev,1,p_qv) = pstate_wrf%moist(  1:ncell,1:nlev,1,p_qv)+dtime*&
                                                        ptend_wrf%rqvblten(1:ncell,1:nlev,1)
          pstate_wrf%moist(  1:ncell,1:nlev,1,p_qc) = pstate_wrf%moist(  1:ncell,1:nlev,1,p_qc)+dtime*&
                                                        ptend_wrf%rqcblten(1:ncell,1:nlev,1)
          pstate_wrf%moist(  1:ncell,1:nlev,1,p_qi) = pstate_wrf%moist(  1:ncell,1:nlev,1,p_qi)+dtime*&
                                                          ptend_wrf%rqiblten(1:ncell,1:nlev,1)
!
! assume pressure remain unchanged during physics, and follow the convention in physics,
! also use mass pressure as air pressure when converting t and th; this is to say,
! we do not allow any acoustic effect inside model physics
!
          pstate_wrf%t_phy(  1:ncell,1:nlev,1)      = pstate_wrf%th_phy( 1:ncell,1:nlev,1)*pstate_wrf%pi_phy( 1:ncell,1:nlev,1)

       end if
    case('RAD')
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthraten(1:ncell,1:nlev,1))
    case default
       if(mpi_rank().eq.0) print*,"no process name is selected, please double check"
       call mpi_abort()
    end select
    return
  end subroutine phys_update_p1a

!------------------------------------------------------------
! All used as a pure F2 operator splitting approach;
! and internally update physics state, only devised for HDC
!------------------------------------------------------------

  subroutine phys_update_p2(dtime, ncell, process_name, scheme_name)

    use grist_hpe_constants,                only: deta_face, deta_full

    real(r8),              intent(in)   :: dtime
    integer(i4),           intent(in)   :: ncell
    character(len=*),      intent(in)   :: process_name
    character(len=*),      intent(in)   :: scheme_name
    !local
    real(r8) :: k

    SELECT CASE(trim(process_name))
    CASE('CU')
       select case(trim(scheme_name))
       case('NTDKV381','TDK')
          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruucuten(1:ncell,1:nlev,1))

          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv,1:nlev,1:ncell) + transpose(ptend_wrf%rqvcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc,1:nlev,1:ncell) + transpose(ptend_wrf%rqccuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi,1:nlev,1:ncell) + transpose(ptend_wrf%rqicuten(1:ncell,1:nlev,1))

! update due to CU
          pstate_wrf%u_phy(  1:ncell,1:nlev,1)      = pstate_wrf%u_phy(  1:ncell,1:nlev,1)     +dtime*ptend_wrf%ruucuten(1:ncell,1:nlev,1)
          pstate_wrf%v_phy(  1:ncell,1:nlev,1)      = pstate_wrf%v_phy(  1:ncell,1:nlev,1)     +dtime*ptend_wrf%rvvcuten(1:ncell,1:nlev,1)
          pstate_wrf%th_phy( 1:ncell,1:nlev,1)      = pstate_wrf%th_phy( 1:ncell,1:nlev,1)     +dtime*ptend_wrf%rthcuten(1:ncell,1:nlev,1)
          pstate_wrf%moist(  1:ncell,1:nlev,1,p_qv) = pstate_wrf%moist(  1:ncell,1:nlev,1,p_qv)+dtime*ptend_wrf%rqvcuten(1:ncell,1:nlev,1)
          pstate_wrf%moist(  1:ncell,1:nlev,1,p_qc) = pstate_wrf%moist(  1:ncell,1:nlev,1,p_qc)+dtime*ptend_wrf%rqccuten(1:ncell,1:nlev,1)
          pstate_wrf%moist(  1:ncell,1:nlev,1,p_qi) = pstate_wrf%moist(  1:ncell,1:nlev,1,p_qi)+dtime*ptend_wrf%rqicuten(1:ncell,1:nlev,1)
          pstate_wrf%t_phy(  1:ncell,1:nlev,1)      = pstate_wrf%th_phy( 1:ncell,1:nlev,1)*pstate_wrf%pi_phy( 1:ncell,1:nlev,1)

          do k = 2,nlev
            pstate_wrf%t8w(1:ncell,k,1) = 0.5*(deta_full(nlev+1-k)/deta_face(nlev+2-k)*pstate_wrf%t_phy(  1:ncell,k-1,1)+&
                                               deta_full(nlev+2-k)/deta_face(nlev+2-k)*pstate_wrf%t_phy(  1:ncell,k,1) )
          end do
          pstate_wrf%t8w(1:ncell,nlev+1,1) = pstate_wrf%t_phy(  1:ncell,nlev,1)
  
       case default
          if(mpi_rank().eq.0) print*,"no CU scheme name is selected, please double check"
          call mpi_abort()
      end select

!================================================
! MP
! Kessler can only be run along
! needs comment other processes
!================================================

    case('MP')
       select case(trim(scheme_name))
       case('KES')
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))
       case('LINV2','WSM6','WSM6V381') ! All six-classs 
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) + transpose(ptend_wrf%rqimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qs ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qs ,1:nlev,1:ncell) + transpose(ptend_wrf%rqsmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qg ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qg ,1:nlev,1:ncell) + transpose(ptend_wrf%rqgmpten(1:ncell,1:nlev,1))

       case('MORR_TWOM_V381','MORR_TWOM_V381_ACE') ! Double-moment
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) + transpose(ptend_wrf%rqimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qs ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qs ,1:nlev,1:ncell) + transpose(ptend_wrf%rqsmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qg ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qg ,1:nlev,1:ncell) + transpose(ptend_wrf%rqgmpten(1:ncell,1:nlev,1))
!
! number concentration, not necessarily for PDC
!
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_nr ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_nr ,1:nlev,1:ncell) + transpose(ptend_wrf%rnrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_ni ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_ni ,1:nlev,1:ncell) + transpose(ptend_wrf%rnimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_ns ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_ns ,1:nlev,1:ncell) + transpose(ptend_wrf%rnsmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_ng ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_ng ,1:nlev,1:ncell) + transpose(ptend_wrf%rngmpten(1:ncell,1:nlev,1))

       case('YLINV381','LINV381') ! five species
! store tendency
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) + transpose(ptend_wrf%rqimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qs ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qs ,1:nlev,1:ncell) + transpose(ptend_wrf%rqsmpten(1:ncell,1:nlev,1))
       case default
          if(mpi_rank().eq.0) print*,"no MP scheme name is selected, please double check"
          call mpi_abort()
      end select

!================================================
! PBL
!================================================

    case('PBL')
       select case(trim(scheme_name))
       case('YSU','YSUV341','YSUV381','S&HV381','CAMUW_V381')
          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruublten(1:ncell,1:nlev,1))

          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvblten(1:ncell,1:nlev,1))

          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) + transpose(ptend_wrf%rqiblten(1:ncell,1:nlev,1))

! update
          pstate_wrf%u_phy(  1:ncell,1:nlev,1)      = pstate_wrf%u_phy(  1:ncell,1:nlev,1)     +dtime*ptend_wrf%ruublten(1:ncell,1:nlev,1)
          pstate_wrf%v_phy(  1:ncell,1:nlev,1)      = pstate_wrf%v_phy(  1:ncell,1:nlev,1)     +dtime*ptend_wrf%rvvblten(1:ncell,1:nlev,1)
          pstate_wrf%th_phy( 1:ncell,1:nlev,1)      = pstate_wrf%th_phy( 1:ncell,1:nlev,1)     +dtime*ptend_wrf%rthblten(1:ncell,1:nlev,1)
          pstate_wrf%moist(  1:ncell,1:nlev,1,p_qv) = pstate_wrf%moist(  1:ncell,1:nlev,1,p_qv)+dtime*ptend_wrf%rqvblten(1:ncell,1:nlev,1)
          pstate_wrf%moist(  1:ncell,1:nlev,1,p_qc) = pstate_wrf%moist(  1:ncell,1:nlev,1,p_qc)+dtime*ptend_wrf%rqcblten(1:ncell,1:nlev,1)
          pstate_wrf%moist(  1:ncell,1:nlev,1,p_qi) = pstate_wrf%moist(  1:ncell,1:nlev,1,p_qi)+dtime*ptend_wrf%rqiblten(1:ncell,1:nlev,1)
          pstate_wrf%t_phy(  1:ncell,1:nlev,1)      = pstate_wrf%th_phy( 1:ncell,1:nlev,1)*pstate_wrf%pi_phy( 1:ncell,1:nlev,1)

          do k = 2,nlev
            pstate_wrf%t8w(1:ncell,k,1) = 0.5*(deta_full(nlev+1-k)/deta_face(nlev+2-k)*pstate_wrf%t_phy(  1:ncell,k-1,1)+&
                                               deta_full(nlev+2-k)/deta_face(nlev+2-k)*pstate_wrf%t_phy(  1:ncell,k,1) )
          end do
          pstate_wrf%t8w(1:ncell,nlev+1,1) = pstate_wrf%t_phy(  1:ncell,nlev,1)
 
       case default
          if(mpi_rank().eq.0) print*,"no PBL scheme name is selected, please double check"
          call mpi_abort()
      end select

    case('RAD')
!
! no selection here because any radiation module only modify heating
!
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthraten(1:ncell,1:nlev,1))
! update due to rad
          pstate_wrf%th_phy( 1:ncell,1:nlev,1)      = pstate_wrf%th_phy( 1:ncell,1:nlev,1)+dtime*ptend_wrf%rthraten(1:ncell,1:nlev,1)
          pstate_wrf%t_phy(  1:ncell,1:nlev,1)      = pstate_wrf%th_phy( 1:ncell,1:nlev,1)*pstate_wrf%pi_phy( 1:ncell,1:nlev,1)
          do k = 2,nlev
            pstate_wrf%t8w(1:ncell,k,1) = 0.5*(deta_full(nlev+1-k)/deta_face(nlev+2-k)*pstate_wrf%t_phy(  1:ncell,k-1,1)+&
                                               deta_full(nlev+2-k)/deta_face(nlev+2-k)*pstate_wrf%t_phy(  1:ncell,k,1) )
          end do
          pstate_wrf%t8w(1:ncell,nlev+1,1) = pstate_wrf%t_phy(  1:ncell,nlev,1)
 
    case default
       if(mpi_rank().eq.0) print*,"no process name is selected, please double check"
       call mpi_abort()
    end select
    return
  end subroutine phys_update_p2

!---------------------------------------------------------------------
! PBL, CU and RAD's uvpt tends are used as RK physics (filled in F1)
! MP 's uvpt tends are used as F2 operator spliting (filled in F3)
! All q tends are F2 operator splitting (filled in F3)
!---------------------------------------------------------------------

  subroutine phys_update_p3(dtime, ncell, process_name, scheme_name)
    real(r8),              intent(in)   :: dtime
    integer(i4),           intent(in)   :: ncell
    character(len=*),      intent(in)   :: process_name
    character(len=*),      intent(in)   :: scheme_name

    SELECT CASE(trim(process_name))
    CASE('CU')

      select case(trim(scheme_name))
        case('TDK','NTDKV381')
          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruucuten(1:ncell,1:nlev,1))

          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvcuten(1:ncell,1:nlev,1))

          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) + transpose(ptend_wrf%rqccuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) + transpose(ptend_wrf%rqicuten(1:ncell,1:nlev,1))
        case default
          if(mpi_rank().eq.0) print*,"no CU scheme name is selected, please double check"
          call mpi_abort()
      end select
!================================================
! MP
! Kessler can only be run along
! needs comment other processes
!================================================

    CASE('MP')
       select case(trim(scheme_name))
       case('KES')
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))
       case('LINV2','WSM6','WSM6V381') ! Six mxrt species
! store tendency
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) + transpose(ptend_wrf%rqimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qs ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qs ,1:nlev,1:ncell) + transpose(ptend_wrf%rqsmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qg ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qg ,1:nlev,1:ncell) + transpose(ptend_wrf%rqgmpten(1:ncell,1:nlev,1))

       case('MORR_TWOM_V381','MORR_TWOM_V381_ACE') ! Double-moment, use this as a template for other two-moment if species are different
! store tendency
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) + transpose(ptend_wrf%rqimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qs ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qs ,1:nlev,1:ncell) + transpose(ptend_wrf%rqsmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qg ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qg ,1:nlev,1:ncell) + transpose(ptend_wrf%rqgmpten(1:ncell,1:nlev,1))
!
! number concentration, not necessarily for PDC
!
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_nr ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_nr ,1:nlev,1:ncell) + transpose(ptend_wrf%rnrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_ni ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_ni ,1:nlev,1:ncell) + transpose(ptend_wrf%rnimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_ns ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_ns ,1:nlev,1:ncell) + transpose(ptend_wrf%rnsmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_ng ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_ng ,1:nlev,1:ncell) + transpose(ptend_wrf%rngmpten(1:ncell,1:nlev,1))

       case('YLINV381','LINV381') ! five species
! store tendency
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qr ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) + transpose(ptend_wrf%rqimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qs ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qs ,1:nlev,1:ncell) + transpose(ptend_wrf%rqsmpten(1:ncell,1:nlev,1))
        case default
          if(mpi_rank().eq.0) print*,"no MP scheme name is selected, please double check"
          call mpi_abort()
      end select
!================================================
! PBL
!================================================
    CASE('PBL')
       select case(trim(scheme_name))
       case('YSU','YSUV341','YSUV381','S&HV381','CAMUW_V381')
          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruublten(1:ncell,1:nlev,1))

          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvblten(1:ncell,1:nlev,1))

          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qv ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qc ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(p_qi ,1:nlev,1:ncell) + transpose(ptend_wrf%rqiblten(1:ncell,1:nlev,1))
       case default
          if(mpi_rank().eq.0) print*,"no PBL scheme name is selected, please double check"
          call mpi_abort()
      end select
    CASE('RAD')
! any radiation module will only update these
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthraten(1:ncell,1:nlev,1))
    CASE default
       if(mpi_rank().eq.0) print*,"no Process name is selected, please double check"
       call mpi_abort()
    END SELECT
    return
  end subroutine phys_update_p3

!------------------------------------------------------------

! PBL, CU, MP RAD's uvpt tends are used as RK physics
! all q tends are F2 operator splitting
!------------------------------------------------------------

  subroutine phys_update_p3a(dtime, ncell, process_name, scheme_name)
    real(r8),              intent(in)   :: dtime
    integer(i4),           intent(in)   :: ncell
    character(len=*),      intent(in)   :: process_name
    character(len=*),      intent(in)   :: scheme_name

    select case(trim(process_name))
    case('CU')
       
       if(trim(scheme_name).eq.'TDK')then

          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruucuten(1:ncell,1:nlev,1))

          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvcuten(1:ncell,1:nlev,1))

          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqccuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqicuten(1:ncell,1:nlev,1))

       end if

!================================================
! MP
! Kessler can only be run along
! needs comment other processes
!================================================

    case('MP')
       if(trim(scheme_name).eq.'KES')then
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))
       end if
       if(trim(scheme_name).eq.'LINV2'.or.trim(scheme_name).eq.'WSM6'.or.trim(scheme_name).eq.'WSM6V381')then ! add six species here
! store tendency
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(5 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(5 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqsmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(6 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(6 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqgmpten(1:ncell,1:nlev,1))
       end if

       if(trim(scheme_name).eq.'YLINV381'.or.trim(scheme_name).eq.'LINV381')then ! five species
! store tendency
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(5 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(5 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqsmpten(1:ncell,1:nlev,1))
       end if

!================================================
! PBL
!================================================
    case('PBL')
       if(trim(scheme_name).eq.'YSU'.or.trim(scheme_name).eq.'YSUV341'.or.trim(scheme_name).eq.'YSUV381'    &
         .or.trim(scheme_name).eq.'S&HV381' .or. trim(scheme_name).eq.'CAMUW_V381')then

          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruublten(1:ncell,1:nlev,1))

          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvblten(1:ncell,1:nlev,1))

          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqiblten(1:ncell,1:nlev,1))
          
       end if
    case('RAD')
    ! any radiation module will only update these
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthraten(1:ncell,1:nlev,1))
    case default
       if(mpi_rank().eq.0) print*,"no process name is selected, please double check"
       call mpi_abort()
    end select
    return
  end subroutine phys_update_p3a
!
! operator splitting: MP, SFC, PBL, CU
! tendency: RAD
!
  subroutine phys_update_p5(dtime, ncell, process_name, scheme_name)
    real(r8),              intent(in)   :: dtime
    integer(i4),           intent(in)   :: ncell
    character(len=*),      intent(in)   :: process_name
    character(len=*),      intent(in)   :: scheme_name

    select case(trim(process_name))

    case('CU')
       if(trim(scheme_name).eq.'TDK')then
          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruucuten(1:ncell,1:nlev,1))

          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvcuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqccuten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqicuten(1:ncell,1:nlev,1))
       end if

!================================================
! MP
! Kessler can only be run along
! needs comment other processes
!================================================

    case('MP')
       if(trim(scheme_name).eq.'KES')then
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))
       end if
       if(trim(scheme_name).eq.'LIN')then
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(5 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(5 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqsmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(6 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(6 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqgmpten(1:ncell,1:nlev,1))
       end if

!================================================
! PBL
!================================================

    case('PBL')
       if(trim(scheme_name).eq.'YSU')then
          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruublten(1:ncell,1:nlev,1))

          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvblten(1:ncell,1:nlev,1))

          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqiblten(1:ncell,1:nlev,1))
       end if
    case('RAD')
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthraten(1:ncell,1:nlev,1))
    case default
       if(mpi_rank().eq.0) print*,"no process name is selected, please double check"
       call mpi_abort()
    end select
    return
  end subroutine phys_update_p5

!
! F2 operator splitting: MP, BL
! RK-F1: RAD-CU
!
  subroutine phys_update_p6(dtime, ncell, process_name, scheme_name)
    real(r8),              intent(in)   :: dtime
    integer(i4),           intent(in)   :: ncell
    character(len=*),      intent(in)   :: process_name
    character(len=*),      intent(in)   :: scheme_name

    select case(trim(process_name))

    case('CU')
       if(trim(scheme_name).eq.'TDK')then
          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruucuten(1:ncell,1:nlev,1))

          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvcuten(1:ncell,1:nlev,1))

          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthcuten(1:ncell,1:nlev,1))

          ptend_f1%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f1%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvcuten(1:ncell,1:nlev,1))

          ptend_f1%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f1%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqccuten(1:ncell,1:nlev,1))

          ptend_f1%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f1%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqicuten(1:ncell,1:nlev,1))
       end if

!================================================
! MP
! Kessler can only be run along
! needs comment other processes
!================================================

    case('MP')
       if(trim(scheme_name).eq.'LIN')then
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(3 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqrmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqimpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(5 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(5 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqsmpten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(6 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(6 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqgmpten(1:ncell,1:nlev,1))
       end if

!================================================
! PBL
!================================================

    case('PBL')
       if(trim(scheme_name).eq.'YSU')then
          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%ruublten(1:ncell,1:nlev,1))

          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rvvblten(1:ncell,1:nlev,1))

          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f3%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqvblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(2 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqcblten(1:ncell,1:nlev,1))

          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) =  &
          ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(4 ,1:nlev,1:ncell) + transpose(ptend_wrf%rqiblten(1:ncell,1:nlev,1))
       end if
    case('RAD')
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) =  &
          ptend_f1%tend_potential_temp_at_pc_full_level%f(1:nlev,1:ncell) + transpose(ptend_wrf%rthraten(1:ncell,1:nlev,1))
    case default
       if(mpi_rank().eq.0) print*,"no process name is selected, please double check"
       call mpi_abort()
    end select
    return
  end subroutine phys_update_p6

 end module grist_PhysW_physics_update
