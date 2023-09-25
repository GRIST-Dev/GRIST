
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: Some utilities for tracer transport
! Revision history:
!                  1. fixer should not be used currently
!----------------------------------------------------------------------------

  module grist_tracer_transport_utils_module

    use grist_domain_types, only: global_domain
    use grist_data_types,   only: scalar_2d_field, scalar_3d_field
    use grist_constants,    only: i4, r8, zero, one, half
    use grist_nml_module,   only: nlev, ntracer, write_verbose, physpkg
    use grist_hpe_constants,only: deta_full, deta_face
    use grist_mpi
! data
    use grist_tracer_transport_vars_module, only: scalar_tracer_mxrt_min

    implicit none

    save

    PUBLIC :: tracer_transport_check_mxrt, &
              tracer_transport_fixer_mxrt, & ! low-efficient, needs optimization
              tracer_transport_qneg3_mxrt

  CONTAINS

   subroutine tracer_transport_check_mxrt(mesh, scalar_tracer_mxrt_at_pc_full_level_n,called_env)
! io
     use omp_lib
     type(global_domain),   intent(in) :: mesh
     type(scalar_3d_field), intent(in) :: scalar_tracer_mxrt_at_pc_full_level_n
     character(len=*),      intent(in) :: called_env
! local
     integer(i4)     :: itracer, iv, ilev
     integer(i4)     :: nsum, nsum_global, ierr
     integer(i4)     :: ii

        nsum = 0
!$omp parallel  private(ii,iv,ilev,itracer) 
!$omp do schedule(static,100) reduction(+:nsum)
     do ii = 1, mesh%nv_compute*nlev*ntracer,1
        iv=ceiling(ii/real(nlev*ntracer,r8))
        ilev=ceiling((ii-(iv-1)*nlev*ntracer)/real(ntracer,r8))
        itracer=ii-(iv-1)*nlev*ntracer-(ilev-1)*ntracer
!        do iv = 1, mesh%nv_compute
!           do itracer = 1, ntracer
!              do ilev = 1, nlev
                 if(scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev,iv).lt.scalar_tracer_mxrt_min(itracer))then
                    nsum = nsum + 1
                 end if
!              end do
!           end do
!        end do
     end do
!$omp end do nowait
!$omp end parallel
#ifndef SEQ_GRIST
        call reduce(nsum, nsum_global, 'sum')
#else
        nsum_global = nsum
#endif
        if(mpi_rank().eq.0.and.write_verbose) print*,"nsum_global=",nsum_global,"in violation of min",trim(called_env)

      return
    end subroutine tracer_transport_check_mxrt

!----------------------------------------------------------------------------------
! This filling procedure is similar to that used by CAM-FV; Although horizontal
! transport is monotonic and refuses undershooting to machine round off, certain
! small negative values may still be generated in some ill conditions
! (diabatic/hvsplit/boundary)
!----------------------------------------------------------------------------------

    subroutine tracer_transport_fixer_mxrt(ncell,scalar_tracer_mxrt_at_pc_full_level_n,&
                                                 scalar_delhp_at_pc_full_level_n)
! io
     integer(i4),           intent(in)    :: ncell
     type(scalar_3d_field), intent(inout) :: scalar_tracer_mxrt_at_pc_full_level_n
     type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_full_level_n
! local
     integer(i4)     :: itracer, iv, ilev
     integer(i4)     :: nsum, nsum_global, ierr
     real(r8)        :: qkold, qfill_up, qfill_lo, qfill

        do iv = 1, ncell
           do itracer = 1, ntracer
              do ilev = 1, nlev
                 if(scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev,iv).lt.scalar_tracer_mxrt_min(itracer))then
             ! as long as tracer_mxrt violates minimum constraint, it will be reset to min value,
             ! regardless whether we borrow from upper or lower levels; when not
             ! borrowing, this will violate mass conservation, but such computational
             ! sources are quite small.
                    qkold = scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev,iv)
                    scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev,iv) = scalar_tracer_mxrt_min(itracer)
             ! borrow from upper or lower levels
                    qfill_up = zero
                    qfill_lo = zero
                    if(ilev.ne.1)    qfill_up = (scalar_tracer_mxrt_min(itracer)-qkold)*&
                                scalar_delhp_at_pc_full_level_n%f(ilev,iv)/scalar_delhp_at_pc_full_level_n%f(ilev-1,iv)
                    if(ilev.ne.nlev) qfill_lo = (scalar_tracer_mxrt_min(itracer)-qkold)*&
                                scalar_delhp_at_pc_full_level_n%f(ilev,iv)/scalar_delhp_at_pc_full_level_n%f(ilev+1,iv)
             ! check how to borrow:
             ! 1) from upper level:
#ifndef NEWFIXER
                    if(ilev.ne.1.and.(scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev-1,iv)-qfill_up).gt.scalar_tracer_mxrt_min(itracer))then
#else
                    if(ilev.ne.1)then
                       if((scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev-1,iv)-qfill_up).gt.scalar_tracer_mxrt_min(itracer))then

#endif
                       scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev-1,iv) = scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev-1,iv)-qfill_up
                       cycle
#ifdef NEWFIXER
                       end if
#endif
             ! 2) from lower level:
#ifndef NEWFIXER
                    else if(ilev.ne.nlev.and.(scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev+1,iv)-qfill_lo).gt.scalar_tracer_mxrt_min(itracer))then
#else
                    else if(ilev.ne.nlev)then
                       if((scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev+1,iv)-qfill_lo).gt.scalar_tracer_mxrt_min(itracer))then
#endif
                       scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev+1,iv) = scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev+1,iv)-qfill_lo
                       cycle
#ifdef NEWFIXER
                       endif
#endif
             ! 3) not borrow:
                    else
                       cycle
                    end if
                 end if
              end do
           end do
        end do
       
      return
    end subroutine tracer_transport_fixer_mxrt
!
! this code is adapted from CAM's qneg3 code
!
 subroutine tracer_transport_qneg3_mxrt(ncell,nlev,ntracer,qqq,called_env)
! io
   integer(i4), intent(in)    :: ncell       ! number of atmospheric columns
   integer(i4), intent(in)    :: nlev        ! number of vertical levels in column
   integer(i4), intent(in)    :: ntracer     ! number of vertical levels in column
   real(r8),    intent(inout) :: qqq(:,:,:)  ! moisture/tracer field
   character*(*), intent(in)  :: called_env  ! name of calling routine
! local
   integer indx(nlev,ncell) ! array of indices of points < qmin
   integer nval(nlev)       ! number of points < qmin for 1 level
   integer nvals            ! number of values found < qmin
   integer nn
   integer iwtmp
   integer iv,ii,k          ! longitude, level indices
   integer m                ! constituent index
   integer iw,kw            ! i,k indices of worst violator
   logical found            ! true => at least 1 minimum violator found
   real(r8) worst           ! biggest violator


   do m = 1, ntracer
      nvals = 0
      found = .false.
      worst = 1.e35_r8
      iw    = -1
!
! Test all field values for being less than minimum value. Set q = qmin
! for all such points. Trace offenders and identify worst one.
!
      if(trim(physpkg).eq.'AMIPC_PHYSICS'  .and. (m == 4 .or. m == 5))then
          qqq(m,:,:) = max(1.e-12_r8,qqq(m,:,:))
          qqq(m,:,:) = min(1.e10_r8, qqq(m,:,:))
      else 
      do k = 1,nlev
         nval(k) = 0
         nn      = 0
         do iv = 1,ncell
            if(qqq(m,k,iv) < scalar_tracer_mxrt_min(m))then
               nn = nn + 1
               indx(k,nn) = iv
            end if
         end do
         nval(k) = nn ! how many for one level
      end do

      do k = 1,nlev
         if (nval(k) > 0) then
            found = .true.
            nvals = nvals + nval(k) ! total count
            iwtmp = -1
! count
            do ii = 1, nval(k)
               iv = indx(k,ii)
               if (qqq(m,k,iv) < worst) then
                  worst = qqq(m,k,iv)
                  iwtmp = ii
               end if
            end do

            if (iwtmp /= -1 ) kw = k
            if (iwtmp /= -1 ) iw = indx(k,iwtmp)
! modify in this way is good
!!!            do ii = 1, nval(k)
!!!               iv = indx(k,ii)
!!!               qqq(m,k,iv) = scalar_tracer_mxrt_min(m)
!!!            end do
            where(qqq(m,:,:) < scalar_tracer_mxrt_min(m)) qqq(m,:,:) = scalar_tracer_mxrt_min(m)

         end if
      end do
  !    if (found .and. abs(worst)>1.e-12_r8) write(6,9000) called_env, m,nvals,iwtmp,scalar_tracer_mxrt_min(m),worst,iw,kw
      end if
   end do

   return
9000 format(' QNEG3 from ',a,':m=',i3,' ncell=',i5, &
            ' Min. mixing ratio violated at ',i4,' points.  Reset to ', &
            1p,e8.1,' Worst =',e8.1,' at i,k=',i4,i3)

  end subroutine tracer_transport_qneg3_mxrt

  end module grist_tracer_transport_utils_module
