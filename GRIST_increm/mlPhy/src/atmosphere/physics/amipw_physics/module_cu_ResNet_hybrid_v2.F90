
!----------------------------------------------------------------------------
! Yiming Wang 2023/05
! Description: Reconstruct ResNet in Fortran
! Revision history:
!----------------------------------------------------------------------------

module module_cu_ResNet_hybrid_v2

#ifdef ResNet
! external
  use, intrinsic :: iso_fortran_env, only : sp => real32
  use ftorch
! para 
  use grist_constants,                    only: r4, i4, r8, one, zero, latvap, cp, one, half,gravity,fillvalue
  use grist_math_module,                  only: lininterp
  use grist_hpe_constants,                only: eta_face_a, eta_face_b, eta_full_a, eta_full_b, p0
  use grist_nml_module,                   only: nlev, nlevp, nlev_inidata
  use grist_mps_nml_module,               only: set_mps_vars, CU_ml_model, mlFilePath
! data struct
  use grist_data_types,                   only: scalar_2d_field, wrap_deallocate_data2d
! function
  use grist_fileio_list_2d_module_par,    only: wrap_read_2d_group_mps

  implicit none

  private              ! Make default type private to the module
  save

  integer,parameter  :: nsize=3, in_channels=5, t0_channels=0
  integer,parameter  :: nkernels=128, out_channels =2
  integer,parameter  :: batch=3720,  nres=10
  integer,parameter  :: seq_len=1,begintime=4,casename=11
  integer,parameter  :: tin_channels=in_channels*seq_len+t0_channels
  integer, parameter :: wp = sp
  integer            :: length,lengtho

  real*4,parameter   :: ran=0.0d0, alpha=0.0d0
  type(torch_module) :: CU_model

  type(scalar_2d_field), public   :: mldata_imin
  type(scalar_2d_field), public   :: mldata_imax
  type(scalar_2d_field), public   :: mldata_omin
  type(scalar_2d_field), public   :: mldata_omax

  public :: ResNet_hybrid_init
  public :: ResNet_hybrid_tend
  public :: ResNet_hybrid_final

  contains
!====================================================================

  subroutine ResNet_hybrid_init(rthcuten,rqvcuten,         &
                     restart,                  &
                     ids, ide, jds, jde, kds, kde,                      &
                     ims, ime, jms, jme, kms, kme,                      &
                     its, ite, jts, jte, kts, kte)
!--------------------------------------------------------------------
   implicit none
!--------------------------------------------------------------------
   integer , intent(in)           ::  ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte

   real,     dimension( ims:ime , kms:kme , jms:jme ) , intent(out) ::  &
                                                              rthcuten, &
                                                              rqvcuten
   integer                :: i, j, k, itf, jtf, ktf
   logical, intent(in)    :: restart

   jtf = jme
   ktf = kme
   itf = ime

   call set_mps_vars()
   
   length  = nlev
   lengtho = nlev
   
   CU_model = torch_module_load(CU_ml_model)

   if(.not.restart)then
     !do j=jts,jtf
     !do k=kts,ktf
     !do i=its,itf
     !   rthcuten(i,k,j) = 0.
     !   rqvcuten(i,k,j) = 0.
     !enddo
     !enddo
     !enddo

        rthcuten = 0.
        rqvcuten = 0.

     !DO j=jts,jtf
     !DO k=kts,ktf
     !DO i=its,itf
     !ENDDO
     !ENDDO
     !ENDDO
   endif

   if(.not.allocated(mldata_imax%f)) allocate(mldata_imax%f(length ,5))
   if(.not.allocated(mldata_imin%f)) allocate(mldata_imin%f(length ,5))
   if(.not.allocated(mldata_omax%f)) allocate(mldata_omax%f(lengtho,2))
   if(.not.allocated(mldata_omin%f)) allocate(mldata_omin%f(lengtho,2))

   call wrap_read_2d_group_mps(trim(mlFilePath),'test_k2f_32.nc','inputsmax' ,length ,9,mldata_imax%f)
   call wrap_read_2d_group_mps(trim(mlFilePath),'test_k2f_32.nc','inputmin'  ,length ,9,mldata_imin%f)
   call wrap_read_2d_group_mps(trim(mlFilePath),'test_k2f_32.nc','outputsmax',lengtho,8,mldata_omax%f)
   call wrap_read_2d_group_mps(trim(mlFilePath),'test_k2f_32.nc','outputmin' ,lengtho,8,mldata_omin%f)

   return
end subroutine ResNet_hybrid_init

subroutine ResNet_hybrid_tend(dt,itimestep,stepcu               &
                ,u3d,v3d,th3d,qv3d                              &
                ,pcps,p8w,lat,qfx,raincv,cu_act_flag,dx         &
                ,ids,ide, jds,jde, kds,kde                      &
                ,ims,ime, jms,jme, kms,kme                      &
                ,its,ite, jts,jte, kts,kte                      &
                ,rthcuten,rqvcuten)
!inputs
!-- u3d         3d u-velocity interpolated to theta points (m/s)
!-- v3d         3d v-velocity interpolated to theta points (m/s)
!-- th3d        3d potential temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- therad3d        3d heating rate due to radiation process ratio (kg/kg)
!-- p8w         3d hydrostatic pressure at full levels (pa)
!-- pcps        3d hydrostatic pressure at half levels (pa)
!s-1)
!-- rthcuten          theta tendency due to 
!                 cumulus scheme precipitation (k/s)
!-- rqvcuten          qv tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- dz8w        dz between full levels (m)
!-- dt          time step (s)
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-------------------------------------------------------------------
      integer, intent(in) ::            ids,ide, jds,jde, kds,kde,      &
                                        ims,ime, jms,jme, kms,kme,      &
                                        its,ite, jts,jte, kts,kte,      &
                                        itimestep,                      &
                                        stepcu

      real,    intent(in) ::                                            &
                                        dt
      real,    dimension(ims:ime), intent(in) ::            dx


      real,    dimension(ims:ime, jms:jme), intent(out) ::               &
                                        raincv
      real,    dimension(ims:ime), intent(in) ::               &
                                        lat                                       
      real,    dimension(ims:ime, jms:jme) ::                           &
                                        qfx                            


      logical, dimension(ims:ime,jms:jme), intent(inout) ::             &
                                        cu_act_flag

      real,    dimension(ims:ime, kms:kme, jms:jme), intent(in) ::      &
                                        pcps,                           &
                                        qv3d,                           &
                                        p8w,                          &
                                        th3d,                            &
                                        u3d,                            &
                                        v3d

!--------------------------- optional vars ----------------------------
!--------------------------- optional vars ----------------------------

      real, dimension(ims:ime, kms:kme, jms:jme),                       &
               optional, intent(inout) ::                               &
                                        rqvcuten,                       &
                                        rthcuten
!--------------------------- local vars ------------------------------
      real      ::                                      &
                                        delt,                           &
                                        rdelt

      real     , dimension(its:ite) ::                  &
                                        rcs,                            &
                                        rn,                             &
                                        evap,                           &
                                        heatflux
      integer  , dimension(its:ite) ::  slimsk


      real     , dimension(its:ite, kts:kte+1) ::         &
                                        prsi,           &
                                        ghti,           &
                                          zl,           &
                                          zi
      integer, dimension(its:ite) ::                                    &
                                        kbot,                           &
                                        ktop

      integer ::                                                        &
                                        i,                              &
                                        im,                             &
                                        j,                              &
                                        k,                              &
                                        km,                             &
                                        kp,                             &
                                        kx,                             &
                                        kx1,                            &
                                        ilev
      real     , dimension(its:ite, 1:lengtho) ::  q1, q2
      real     , dimension(1:30) ::  q1w,q2w  

      integer                         :: nlev_eff,nlev_eff_t,nlev_eff_b
      !real                            :: pface(29+1), hpface(29+1)
      !real                            :: pfull(29),   hpfull(29)
      real                            :: dpfull(nlev) !,  dhpfull(29)
      real                            :: hpmfull(nlev), hpmface(nlevp)
      integer(i4)                     :: qlev
!-------other local variables----
      integer                         :: zz, kfly
!-----------------------------------------------------------------------
!
!
!***  check to see if this is a convection timestep
!

!-----------------------------------------------------------------------
      kfly = 18
      
      do j=jts,jte
         do i=its,ite
           cu_act_flag(i,j)=.true.
         enddo
      enddo

      im       = ite-its+1
      kx       = kte-kts+1
      kx1      = kx+1
      delt     = dt*stepcu
      rdelt    = 1./delt
      raincv   = 0

!-------------  j loop (outer)
!--------------------------------------------------

   do j=jts,jte
      
    call Resnet(q1(:,:),q2(:,:),u3d(:,:,j),v3d(:,:,j), &
                th3d(:,:,j),qv3d(:,:,j), &
                pcps(:,:,j),kms,kme, its, ite)
        q1=q1/86400   !K/s
        q2=q2/86400*cp/latvap !kg/kg/s

    do i = its, ite
!        if ( lat(i).lt. -70.0) then
!            q1(i,:)=zero
!            q2(i,:)=zero
!        endif

        do ilev=nlev-kfly-1,nlev
           dpfull(ilev)=p8w(i,31-ilev,j)-p8w(i,32-ilev,j)
           raincv(i,1)=raincv(i,1)+q2(i, ilev)/gravity*dpfull(ilev)
        enddo

#ifdef diag_qfx
        raincv(i,1)=raincv(i,1)+qfx(i,1)
#endif

        raincv(i,1)=amax1(0.0,raincv(i,1)) !kg/m^2/s ->mm/s

        rthcuten(i,1:kfly,1) =  q1(i, lengtho:lengtho-kfly+1:-1)
        rqvcuten(i,1:kfly,1) = -q2(i, lengtho:lengtho-kfly+1:-1)

!         rthcuten(i,1:30,1)=q1(i, 30:1:-1)
!         rqvcuten(i,1:30,1)=-q2(i, 30:1:-1)

      enddo
   enddo

end subroutine ResNet_hybrid_tend

subroutine Resnet(q1, q2, u3d, v3d,  &
                  th3d, qv3d, &
                  pcps, kms,kme,its, ite)
    implicit none
    integer , intent(in)           ::   its, ite
    integer , intent(in)           ::   kms, kme


    real,    dimension(its:ite, kms:kme), intent(in) ::      &
                                        pcps,                            &
                                        qv3d,                           &
                                        th3d,                            &
                                        u3d,                            &
                                        v3d                            

    real,     dimension(its:ite, 1:lengtho) , intent(out) ::  &
                                                              q1, &
                                                              q2
                                                             ! rthften,rqvften
    integer :: iname, ires
    integer :: iseq
    integer :: i, j, k
    integer :: ptop_res, ptop_reso

    real(r4)       :: inputs( length ,  in_channels*seq_len+t0_channels)
    real(r4)       :: outputs( lengtho,  out_channels)
    real(r4)       :: predictions(lengtho , out_channels)
    real(r4)       :: sinput(its:ite,length, in_channels*seq_len+t0_channels), stemp(length, nkernels )
    real(r4)       :: soutput(its:ite,lengtho, out_channels), spred(lengtho, out_channels)
    real(r4)       :: output_gap(nkernels),output_dense(lengtho*out_channels)
    real(r4)       :: inputs_max( length, in_channels*seq_len+t0_channels)
    real(r4)       :: inputs_min( length, in_channels*seq_len+t0_channels)
    real(r4)       :: outputs_max( lengtho, out_channels), outputs_min(lengtho,out_channels)
    real(r4)       :: weights0(nsize,tin_channels,nkernels), bias0(nkernels)
    real(r4)       :: weightsa(nsize,nkernels,nkernels),biasa(nkernels),weightsb(nsize,nkernels,nkernels),biasb(nkernels)
    real(r4)       :: weights1(nsize,nkernels,out_channels), bias1(out_channels)

    real           ::  qs,es,ilev,qtemp,qfix,tfix
    integer        :: status

    inputs_max  = mldata_imax%f
    inputs_min  = mldata_imin%f
    outputs_min = mldata_omin%f
    outputs_max = mldata_omax%f

    ptop_res    = nlev
    ptop_reso   = nlev

    sinput(:,1:ptop_res,1) = qv3d(:,ptop_res:1:-1)
    sinput(:,1:ptop_res,2) = th3d(:,ptop_res:1:-1)
    sinput(:,1:ptop_res,3) = u3d (:,ptop_res:1:-1)
    sinput(:,1:ptop_res,4) = v3d (:,ptop_res:1:-1)
    sinput(:,1:ptop_res,5) = pcps(:,ptop_res:1:-1)

    do i = its, ite
!        sinput(i,:,:)=amin1(sinput(i,:,:),inputs_max(:,:))
!        sinput(i,:,:)=amax1(sinput(i,:,:),inputs_min(:,:))
        sinput(i,:,:)=(sinput(i,:,:)-inputs_min(:,:))/(inputs_max(:,:)-inputs_min(:,:))*2.d0-1.d0
    end do

!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
!     if (.true.) then
!     do i=1, in_channels*seq_len+t0_channels
!     do j=1, length
!     do k=its,ite
!         if (sinput(k,j,i)>=1.d0+ran)then
!             sinput(k,j,i)=1.d0+ran
!         elseif (sinput(k,j,i)<=-1.d0-ran)then
!             sinput(k,j,i)=-1.d0-ran
!         endif
!     enddo
!     enddo
!     enddo
!     endif
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

!    if(.true.)then
       where(sinput>= 1.d0+ran)   sinput = 1.d0  + ran
       where(sinput<=-1.d0-ran)   sinput = -1.d0 - ran
!    end if

!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
!     if (.true.)then
!     do i=1,seq_len*in_channels+t0_channels
!         do j=1,length
!         do k=its,ite
!         if (inputs_max(j,i)<=1e-29)then
!             sinput(k,j,i)=0.d0
!         endif
!         enddo
!         enddo
!     enddo
!     endif
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

!    if(.true.)then
       do k=its,ite
          where(inputs_max(:,:)<=1e-29) sinput(k,:,:) = 0.d0
       end do
!    end if

   !!!!!
   !forward
   !!!!!
   call cu_dnn_forward(sinput, soutput, its, ite)

! ResNet complete
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
!     if (.true.)then
!     do i=1, out_channels
!         do j=1,lengtho
!         do k = its,ite
!         if (soutput(k,j,i)>=1.d0+ran)then
!             soutput(k,j,i)=1.d0+ran
!         elseif (soutput(k,j,i)<=-1.d0-ran)then
!             soutput(k,j,i)=-1.d0-ran
!         endif
!         enddo
!         enddo
!     enddo
!     endif
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

!    if (.true.)then
        where(soutput>=1.d0+ran)  soutput = 1.d0+ran
        where(soutput<=-1.d0-ran) soutput = -1.d0-ran
!    end if

    do i = its, ite
!        soutput(i,:,:)=amin1(soutput(i,:,:),outputs_max(:,:))
!        soutput(i,:,:)=amax1(soutput(i,:,:),outputs_min(:,:))
        soutput(i,:,:)=(soutput(i,:,:)+1.d0)*(outputs_max(:,:)-outputs_min(:,:))/2.d0+outputs_min(:,:)
    end do

!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
!     if (.true.)then
!     do i=1, out_channels
!         do j=1,lengtho
!         do k=its,ite
!         if (abs(soutput(k,j,i))<=1e-29)then
!             soutput(k,j,i)=0.d0
!         endif
!         enddo
!         enddo
!     enddo
!     endif
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

!    if (.true.)then
        where(abs(soutput)<=1e-29) soutput=0.d0    
!    end if

    q1 = 0.d0
    q2 = 0.d0

    q1(:,1:ptop_reso)=soutput(:,1:ptop_reso,1)!*3-sinput(1:ptop_res,4)*2
    q2(:,1:ptop_reso)=soutput(:,1:ptop_reso,2)!*3-sinput(1:ptop_res,3)*2

end subroutine Resnet

subroutine cu_dnn_forward(input, output, ims, ime)
    integer, intent(in) :: ims, ime

    real(r4), dimension(ims:ime, length, in_channels*seq_len+t0_channels), intent(in) :: input
    real(r4), dimension(ims:ime,lengtho, out_channels), intent(out) :: output

    real(wp), dimension(ims:ime,length, in_channels*seq_len+t0_channels), target :: in_data
    real(wp), dimension(ims:ime,lengtho, out_channels),  target :: out_data

    integer, parameter :: in_dims  = 3
    integer, parameter :: out_dims = 3
    integer, parameter :: n_inputs = 1

    integer :: in_layout(in_dims)   = [1,2,3]
    integer :: out_layout(out_dims) = [1,2,3]    

    type(torch_tensor), dimension(1) :: model_input
    type(torch_tensor)  :: model_output

    in_data(:,:,:)  = input(:,:,:)

    model_input(1)  = torch_tensor_from_array(in_data, in_layout, torch_kCPU)
    model_output    = torch_tensor_from_array(out_data, out_layout, torch_kCPU)
      ! Infer
    call torch_module_forward(CU_model, model_input, n_inputs, model_output)
    
    output(:,:,:)   = out_data(:,:,:)
 
    call torch_tensor_delete(model_input(1))
    call torch_tensor_delete(model_output)

end subroutine cu_dnn_forward

subroutine ResNet_hybrid_final

  if(allocated(mldata_imax%f)) deallocate(mldata_imax%f)
  if(allocated(mldata_imin%f)) deallocate(mldata_imin%f)
  if(allocated(mldata_omax%f)) deallocate(mldata_omax%f)
  if(allocated(mldata_omin%f)) deallocate(mldata_omin%f)

end subroutine ResNet_hybrid_final

#endif

end module module_cu_ResNet_hybrid_v2
