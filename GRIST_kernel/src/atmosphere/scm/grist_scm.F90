
!----------------------------------------------------------------------------
! Copyright: Copyright (c) 2016-2019
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
! Created on 2018
! Version 1.0
! Description: Global-Regional Integrated forecast SysTem-SCM (GRIST-SCM)
!          Major Program driver of GRIST-A, jointly developed by Chinese
!          Academy of Meteorological Sciences and National Supercomputing
!          Center in Wuxi. Unless otherwise mentioned, this model follows
!          the above License.
!
! Revision history:
!----------------------------------------------------------------------------

 PROGRAM GRIST_scm

! infrastructure
   use grist_constants, only: r8
   use grist_lib
   use grist_mpi
   use grist_destruct_module
   use grist_domain_types,   only: global_domain, &
                                   block_structure,&
                                   global_domain_data
   use grist_nml_module,     only: set_global_vars, lbdry_flag
   use grist_grid_file_read, only: grist_read_grid_file
   use grist_config_mesh,    only: config_mesh_from_file
   use grist_config_partition
! scm driver
   use grist_scm_control_driver, only: grist_scm_init, grist_scm_run, grist_scm_final
! clock
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
     call set_global_vars()
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
!
! init model instance here
!
     call grist_scm_init(local_block)
!
! print some info
!
     if(mpi_rank() == 0) then 
        print*,"MPI Pros Num:  ",mpi_size()
        print*,"low order bdry flag:  ",lbdry_flag
     end if
     call clock_end(initclock)
!
! major stepon
!
     call clock_begin(runclock)
     call grist_scm_run(local_block%full_domain)
     call clock_end(runclock)
!
!----final
!
     call clock_begin(finalclock)
     call grist_scm_final(local_block)
     call clock_end(finalclock)
     call clock_end(clock0)
     call clock_summary
     call mpi_finalize(ierr)

 END PROGRAM GRIST_scm
