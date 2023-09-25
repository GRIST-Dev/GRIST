
!----------------------------------------------------------------------------
! Copyright: Copyright (c) 2016-2022
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!    http://www.apache.org/licenses/LICENSE-2.0
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
! Created on
! Version 1.0
! Description: Global-Regional Integrated forecast SysTem (GRIST)
! website: https://github.com/grist-dev/
! Revision history:
!----------------------------------------------------------------------------

 PROGRAM GRIST_atmos

! infrastructure
   use omp_lib
   use grist_constants, only: r8
   use grist_lib
   use grist_mpi
   use grist_destruct_module,only: destruct_global_mesh_var, destruct_local_mesh_element
   use grist_domain_types,   only: global_domain, block_structure,&
                                   global_domain_data, grist_domain_check_runtime_allocated
   use grist_nml_module,     only: nml_set_global_vars,lbdry_flag, model_instance
   use grist_nml_module,     only: use_tr_mbs, stencil_width, stencil_exchange_flag, pv_order
   use grist_grid_file_read, only: grist_read_grid_file
   use grist_config_mesh,    only: config_mesh_from_file
   use grist_geometric_info_icosh, only: init_geometric_info
   use grist_mesh_weight_icosh,only: calc_trsk_weights_plg,        &
                                     init_blocal_prime_cell,       &
                                     init_blocal_prime_cell_4th,   &
                                     init_blocal_dual_cell,        &
                                     init_blocal_dual_cell_mbs,    &
                                     init_blocal_dual_cell_4th,    &
                                     init_blocal_dual_cell_4th_mbs,&
                                     init_plg_sub_triangle_info,   &
                                     fill_key_compute_pattern2
! gcm model instance
   use grist_gcm_control_driver,        only: grist_gcm_init, grist_gcm_run, grist_gcm_final
#ifdef MISC
   use grist_gcm_dycore_control_driver, only: grist_gcm_dycore_init, grist_gcm_dycore_run, grist_gcm_dycore_final
   use grist_gcm_tracer_control_driver, only: grist_gcm_tracer_init, grist_gcm_tracer_run, grist_gcm_tracer_final
#endif
   use grist_config_partition
   use grist_datam_initial_data_module,  only: grist_initial_data_destruct
! timing
   use grist_clocks,             only: clock_id, clock_begin, clock_end, clock_summary, clock0
   
   implicit none

     type(block_structure) ,target     :: local_block
     integer                           :: comm,ierr
     logical                           :: flag_rank = .true.
     character(len=6)                  :: c_tmp
     type(global_domain_data)          :: mesh
     integer                           :: initclock, runclock, finalclock

     call grist_comm_set(MPI_COMM_WORLD)
     comm = grist_comm()

     call mpi_init(ierr)

     clock0     = clock_id('Total')
     initclock  = clock_id('Initialization')
     runclock   = clock_id('Main loop')
     finalclock = clock_id('Final')

     call clock_begin(clock0)
     call clock_begin(initclock)
!
! read mesh 
!
     call nml_set_global_vars()
     call grist_read_grid_file(mesh, comm)
     call config_mesh_from_file(mesh)
     call show_basic_info(mesh, flag_rank)
!
! mesh partition and init local mesh
!
     call domain_decompse(mesh, comm)
     if(mpi_rank() == 0) print *,"config_sub_domain ..."
     call config_sub_domain(mesh, local_block, comm)
     call destruct_global_mesh_var(mesh)
     local_block%full_domain%local_block => local_block
     if(mpi_rank() == 0) print *,"config_sub_domain end"

#ifdef USE_OMP
      if(mpi_rank().eq.0) print*, "thread number is:", omp_get_num_procs()
#endif

!
! init mesh weight
!
     call grist_all_mesh_weight_init(local_block)
     if(mpi_rank() == 0) print *,"grist_all_mesh_weight_init"
!
! init model instance here
!
     select case(trim(model_instance))
     case('gcm')
        call grist_gcm_init(local_block)
#ifdef MISC
     case('dycore')
        call grist_gcm_dycore_init(local_block)
     case('tracer')
        call grist_gcm_tracer_init(local_block)
#endif
     case default
        call grist_gcm_init(local_block)
     end select
!
! part2: added for runtime opt: when using any of CP2_XX (XX=VTX,PLG,EDT,EDP,TRI)
! it will completely delete mesh%XX element, and use all cp1 and cp2
! for runtime compute pattern
!
     call destruct_local_mesh_element(local_block%full_domain)
! clean initial data in memory
     call grist_initial_data_destruct
!
! print some info
!
     if(mpi_rank() == 0) then 
        print*, "MPI Pros Num:  ",mpi_size()
        print*, "low order bdry flag:  ",lbdry_flag
     end if
     if(stencil_width >= 4)then
        stencil_exchange_flag = .false.
        if(mpi_rank() == 0) print*,"stencil tpye:  local    halo:",stencil_width
     else
        stencil_exchange_flag = .true.
        if(mpi_rank() == 0) print*,"stencil tpye:  exchange    halo:",stencil_width
     end if
     call clock_end(initclock)
!
! major stepon
!
     call grist_domain_check_runtime_allocated(local_block%full_domain)

     call clock_begin(runclock)
     select case(trim(model_instance))
     case('gcm')
        call grist_gcm_run(local_block%full_domain)
#ifdef MISC
     case('dycore')   ! this should produce exactly the same results as called in gcm's dycore working_mode
        call grist_gcm_dycore_run(local_block%full_domain)
     case('tracer')   ! this should produce exactly the same results as called in gcm's tracer working_mode
        call grist_gcm_tracer_run(local_block%full_domain)
#endif
     case default
        call grist_gcm_run(local_block%full_domain)
     end select

     call clock_end(runclock)
!
! final
!
     call clock_begin(finalclock)
     select case(trim(model_instance))
     case('gcm')
        call grist_gcm_final(local_block%full_domain)
#ifdef MISC
     case('dycore')
        call grist_gcm_dycore_final(local_block%full_domain)
     case('tracer')
        call grist_gcm_tracer_final(local_block%full_domain)
#endif
     case default
        call grist_gcm_final(local_block%full_domain)
     end select

     call clock_end(finalclock)
     call clock_end(clock0)
     call clock_summary
     call mpi_finalize(ierr)

  contains

  subroutine grist_all_mesh_weight_init(local_block)
    type(block_structure) ,target, intent(inout)  :: local_block
    integer                           :: iii, jjj, kkk
!
! mesh weight
!
     call calc_trsk_weights_plg(local_block%full_domain)
     call init_blocal_prime_cell(local_block%full_domain)
     call init_plg_sub_triangle_info(local_block%full_domain)
#ifndef USE_HALO2
     call init_blocal_prime_cell_4th(local_block%full_domain)
#endif
     call init_geometric_info(local_block%full_domain)

     IF(pv_order(1).gt.2.or.pv_order(2).gt.2.or.pv_order(3).gt.2)then
        if(use_tr_mbs)then
           call init_blocal_dual_cell_mbs(local_block%full_domain)
        else
           call init_blocal_dual_cell(local_block%full_domain)
        end if
     END IF
!
! part1: added for runtime opt: when using any of CP2_XX (XX=VTX,PLG,EDT,EDP,TRI)
! it will completely delete mesh%XX element, and use all cp1 and cp2
! for runtime compute pattern
!
     call fill_key_compute_pattern2(local_block%full_domain)

!
! write some info
!
     write(c_tmp,'(I0)') mpi_rank()
     if(.false.)then
     open(1,file="stencil_number_4th-"//trim(adjustl(c_tmp))//".txt")
     open(2,file="stencil_index_4th-"//trim(adjustl(c_tmp))//".txt")
     open(3,file="my_edge_on_edge_num-"//trim(adjustl(c_tmp))//".txt")
     open(4,file="ur_cell_on_edge_num-"//trim(adjustl(c_tmp))//".txt")
     do kkk = 1,local_block%full_domain%nv_full
       !write(1,'(I0)') local_block%full_domain%plg(kkk)%stencil_number_4th
       !do jjj = 1,size(local_block%full_domain%plg(kkk)%stencil_index_4th)
       !   write(2,'(I0)') local_block%full_domain%plg(kkk)%stencil_index_4th(jjj)
       !end do
       !do iii = 1,2
         !write(3,'(I0)') local_block%full_domain%edt(k)%my_edge_on_edge_num(i)
         !write(4,'(I0)') local_block%full_domain%edt(k)%ur_cell_on_edge_num(i)
       !end do
     end do
     close(4)
     close(3)
     close(2)
     close(1)
     end if

     return
  end subroutine grist_all_mesh_weight_init

 END PROGRAM GRIST_atmos
