
!-----------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: This routines contain some useful routines
!  for the GRIST model
!-----------------------------------------------------------

 module grist_util_module

  !use grist_constants,       only: i4, r8
  use grist_constants_dbl,   only: i4, r8
  use grist_domain_types,    only: global_domain, block_structure
  use grist_math_module,     only: norm, cross_product
  use grist_mpi
  use grist_lib

  implicit none

  private

  public ::  wind_edge_to_cell  , & ! deprecated version
             write_string       , &
             exchange_init_data  ! deprecated version

contains 

!================================================
!   For generating U&V winds
!    It returns the north and east unit vector
!    at a given point (edge point)
!    once called, stored in mesh
!================================================

  subroutine calc_north_east(mesh)

! io
   type(global_domain), intent(inout)  :: mesh
! local
   integer(i4)                          :: ie
   real(r8), dimension(1:3)             :: p1     ! point vector at the north
   real(r8), dimension(1:3)             :: p2     ! point vector at the north
   real(r8), dimension(1:3)             :: north  ! vector pointing the north
   real(r8), dimension(1:3)             :: east   ! vector pointing the east
   real(r8), dimension(1:3)             :: tmp    ! tmp vector
   real(r8), dimension(1:3)             :: z

   do ie = 1, mesh%ne
! north vector
       z(1) = 0._r8
       z(2) = 0._r8
       z(3) = 1._r8
       east = cross_product(z,mesh%edp(ie)%c%p)
       east = east/norm(east)
       north = cross_product(mesh%edp(ie)%c%p,east)

       mesh%edp(ie)%c%north = north
       mesh%edp(ie)%c%east  = east
   end do
 
  return

  end subroutine calc_north_east

! old and depreciated version

  subroutine wind_edge_to_cell(mesh,wind_at_edge,wind_at_cell)
! io
    type(global_domain)    , intent(in)    :: mesh
    real(r8), allocatable  , intent(in)    :: wind_at_edge(:)
    real(r8), allocatable  , intent(inout) :: wind_at_cell(:)
! local
    real(r8)                            :: tmp_sum
    integer(i4)                         :: iv
    integer(i4)                         :: ie
    integer(i4)                         :: global_index_edge

    do iv = 1, mesh%nv
          tmp_sum     = 0._r8
       do ie = 1, mesh%vtx(iv)%nnb
          global_index_edge       = mesh%vtx(iv)%ed(ie)
          tmp_sum                 = tmp_sum+wind_at_edge(global_index_edge)
       end do
          wind_at_cell(iv)      = tmp_sum/mesh%vtx(iv)%nnb
    end do

    return

   end subroutine wind_edge_to_cell

   subroutine write_string(num,c_num)
!io
   integer(i4),      intent(in)    :: num
   character(len=*), intent(inout) :: c_num

   if(num.gt.9)then
      write(c_num,'(i2)') num
   else
     write(c_num,'(i1)') num
   end if

   return
   end subroutine write_string
   
   subroutine exchange_init_data(local_block,tmp_data,type_data)
    !io
    type(block_structure), intent(in)  ,target          :: local_block
    !dim1 is data; dim2 is all data of nv/ne/nt
    real (r8), dimension(:,:),allocatable,intent(inout) :: tmp_data
    integer,intent(in)                                  :: type_data!0:nv;1:nt;6:ne

    !local
    integer, allocatable :: reqs_send(:), reqs_recv(:)
    integer :: dst, ierr, c, iv, ie, it
    integer :: i ,iblock,j,k,n
    integer :: tmp_n,n_v,n_t,n_e
    integer :: h_n,b_n 
    integer :: num_dim1,num_dim2,num_dim3,num_data_tmp
    INTEGER, allocatable  :: status(:,:) !STATUS(MPI_STATUS_SIZE)
    type(map), pointer    :: g2l_v, g2l_t, g2l_e
    character(len=128)    :: tmp_c

    real(r8),allocatable                    :: field_data_send(:,:,:)!bdry
    real(r8),allocatable                    :: field_data_recv(:,:,:)!halo
 
    h_n = 0
    b_n = 0
    iblock = mpi_rank()

    num_data_tmp = size(tmp_data,1)

    !print*,"num:",num_data_tmp,size(tmp_data,2)

    g2l_v => local_block%full_domain%map_g2l_v
    g2l_e => local_block%full_domain%map_g2l_e
    g2l_t => local_block%full_domain%map_g2l_t

    do c = 1, local_block%nnb
       select case(type_data)
       !nv
       case(0)
         if(h_n .lt. size(local_block%recv_sites(c)%v)) &
            h_n = size(local_block%recv_sites(c)%v)
         if(b_n .lt. size(local_block%send_sites(c)%v)) &
            b_n = size(local_block%send_sites(c)%v)
       !nt
       case(1)
         if(h_n .lt. size(local_block%recv_sites(c)%t)) &
            h_n = size(local_block%recv_sites(c)%t)
         if(b_n .lt. size(local_block%send_sites(c)%t)) &
            b_n = size(local_block%send_sites(c)%t)
       !ne
       case(6)
         if(h_n .lt. size(local_block%recv_sites(c)%e)) &
            h_n = size(local_block%recv_sites(c)%e)
         if(b_n .lt. size(local_block%send_sites(c)%e)) &
            b_n = size(local_block%send_sites(c)%e)
       case default
         print*,"exchange data error: Invalid scalar field data type"
         stop
       end select 
    end do

    allocate(field_data_send(num_data_tmp,b_n,local_block%nnb))
    allocate(field_data_recv(num_data_tmp,h_n,local_block%nnb))

    !assign send_data value
    do c = 1, local_block%nnb
       select case(type_data)
       !nv
       case(0)
         n_v = size(local_block%send_sites(c)%v)
         do iv=1, n_v
            tmp_n = g2l_v%find(local_block%send_sites(c)%v(iv))
            field_data_send(:,iv,c) = tmp_data(:,tmp_n)
         end do
       !nt
       case(1)
         n_t = size(local_block%send_sites(c)%t)
         do it=1, n_t
            tmp_n = g2l_t%find(local_block%send_sites(c)%t(it))
            field_data_send(:,it,c) = tmp_data(:,tmp_n)
         end do
       !ne
       case(6)
         n_e = size(local_block%send_sites(c)%e)
         do ie=1, n_e
            tmp_n = g2l_e%find(local_block%send_sites(c)%e(ie))
            field_data_send(:,ie,c) = tmp_data(:,tmp_n)
         end do
       case default
         print*,"Invalid scalar field data type"
       end select
    end do

    allocate(reqs_send(local_block%nnb))
    allocate(reqs_recv(local_block%nnb))
    !---------------------------------
    !exchange data
    do c = 1, local_block%nnb
       dst = local_block%send_sites(c)%cpuid

       select case(type_data)
       case(0)
         tmp_n = size(local_block%send_sites(c)%v)
       case(1)
         tmp_n = size(local_block%send_sites(c)%t)
       case(6)
         tmp_n = size(local_block%send_sites(c)%e)
       case default
         print*,"Invalid scalar field data type"
       end select

       if(r8 .eq. 8) then
         call MPI_Isend(field_data_send(:,:,c), &
                        tmp_n*num_data_tmp, MPI_REAL8, &
                        dst, 102, local_block%comm, reqs_send(c), ierr)
       elseif(r8 .eq. 4) then
         call MPI_Isend(field_data_send(:,:,c), &
                        tmp_n*num_data_tmp, MPI_REAL, &
                        dst, 102, local_block%comm, reqs_send(c), ierr)
       end if

       dst = local_block%recv_sites(c)%cpuid

       select case(type_data)
       case(0)
         tmp_n = size(local_block%recv_sites(c)%v)
       case(1)
         tmp_n = size(local_block%recv_sites(c)%t)
       case(6)
         tmp_n = size(local_block%recv_sites(c)%e)
       case default
         print*,"Invalid scalar field data type"
       end select

       if(r8 .eq. 8) then
         call MPI_Irecv(field_data_recv(:,:,c), &
                        tmp_n*num_data_tmp, MPI_REAL8, &
                        dst, 102, local_block%comm, reqs_recv(c), ierr)
       elseif(r8 .eq. 4) then
         call MPI_Irecv(field_data_recv(:,:,c), &
                        tmp_n*num_data_tmp, MPI_REAL, &
                        dst, 102, local_block%comm, reqs_recv(c), ierr)
       end if
   
    end do

    allocate(status(MPI_STATUS_SIZE, local_block%nnb))

    call MPI_WaitAll(local_block%nnb, reqs_send, status, ierr)
    call MPI_WaitAll(local_block%nnb, reqs_recv, status, ierr)

    !get recv_data value
    do c = 1, local_block%nnb
       select case(type_data)
       !nv
       case(0)
         n_v = size(local_block%recv_sites(c)%v)
         do iv=1, n_v
           tmp_n = g2l_v%find(local_block%recv_sites(c)%v(iv))
           tmp_data(:,tmp_n) = field_data_recv(:,iv,c)
         end do
       !nt
       case(1)
         n_t = size(local_block%recv_sites(c)%t)
         do it=1, n_t
           tmp_n = g2l_t%find(local_block%recv_sites(c)%t(it))
           tmp_data(:,tmp_n) = field_data_recv(:,it,c)
         end do
       !ne
       case(6)
         n_e = size(local_block%recv_sites(c)%e)
         do ie=1, n_e
           tmp_n = g2l_e%find(local_block%recv_sites(c)%e(ie))
           tmp_data(:,tmp_n) = field_data_recv(:,ie,c)
         end do
       case default
         print*,"Invalid scalar field data type"
       end select
    end do

!	call exchange_data_2d_clean(field_head)
!	print*,"scalar_height_at_prime_cell_send_data:",scalar_height_at_prime_cell_recv_data

    if(.false.)  then

      do c = 1, local_block%nnb

         select case(type_data)
         case(0)
           tmp_c = "cell"
         case(1)
           tmp_c = "tri"
         case(6)
           tmp_c = "edge"
         case default
           print*,"Invalid scalar field data type"
         end select

         open(101, file=i2s(iblock)//".send."//trim(tmp_c)//"."//i2s(local_block%send_sites(c)%cpuid)//".txt")
         write(101, "(G20.10)", advance='yes') field_data_send(:,:,c)
         close(101)
 
         open(2, file=i2s(iblock)//".recv."//trim(tmp_c)//"."//i2s(local_block%recv_sites(c)%cpuid)//".txt")
         write(2, "(G20.10)", advance='yes') field_data_recv(:,:,c)
         close(2)
       
      end do
    end if

    !clean  
    deallocate(field_data_send)
    deallocate(field_data_recv)
    deallocate(reqs_send)
    deallocate(reqs_recv)
    deallocate(status)
        
  end subroutine exchange_init_data

 end module grist_util_module
