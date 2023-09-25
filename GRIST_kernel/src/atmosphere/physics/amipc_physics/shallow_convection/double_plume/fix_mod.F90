Module fix_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! chuwc 2020 !!!!!!!!!!!!!!!!!!!!!!!!!
!!!! kinetic energy fix !!!!!!!!!!!!!!!!!
!!!! Don't change !!!!!!!!!!!!!!!!!!!!!!!

use grist_constants, only: r8, gravity
use grist_mpi

implicit none

public :: fix

private
real(r8) :: g

contains

subroutine fix(uten_s,uten_d,vten_s,vten_d, &
               uflx_s,uflx_d,vflx_s,vflx_d,&
               cnt_sh, u0,v0,dt,dp0,t0,cp,    &
               iend,mkx,sten_s, isdeep)

integer, intent(in) :: iend
integer, intent(in) :: mkx
real(r8),intent(in) :: cnt_sh(iend)
real(r8),intent(in) :: dt
real(r8),intent(in) :: uten_s(mkx,iend)
real(r8),intent(in) :: uten_d(mkx,iend)
real(r8),intent(in) :: vten_s(mkx,iend)
real(r8),intent(in) :: vten_d(mkx,iend)
real(r8),intent(in) :: uflx_s(0:mkx,iend)
real(r8),intent(in) :: uflx_d(0:mkx,iend)
real(r8),intent(in) :: vflx_s(0:mkx,iend)
real(r8),intent(in) :: vflx_d(0:mkx,iend)
real(r8),intent(in) :: u0(mkx,iend)
real(r8),intent(in) :: v0(mkx,iend)
real(r8),intent(in) :: dp0(mkx,iend)
real(r8),intent(in) :: t0(mkx,iend)
real(r8),intent(in) :: cp
real(r8),intent(inout) :: sten_s(mkx,iend)
logical, intent(in) :: isdeep

integer i, k

real(r8) :: ufp,uf,ufm,vfp,vf,vfm
real(r8) :: origin,new
real(r8) :: uflx_tot(0:mkx,iend)
real(r8) :: vflx_tot(0:mkx,iend)
real(r8) :: uten_tot(mkx,iend)
real(r8) :: vten_tot(mkx,iend)
real(r8) :: uflx(0:mkx,iend)
real(r8) :: vflx(0:mkx,iend)
real(r8) :: uten(mkx,iend)
real(r8) :: vten(mkx,iend)

real(r8) :: uflx_sh(0:mkx,iend)
real(r8) :: vflx_sh(0:mkx,iend)
real(r8) :: uten_sh(mkx,iend)
real(r8) :: vten_sh(mkx,iend)
real(r8) :: uflx_dp(0:mkx,iend)
real(r8) :: vflx_dp(0:mkx,iend)
real(r8) :: uten_dp(mkx,iend)
real(r8) :: vten_dp(mkx,iend)
real(r8) :: sten_sh(mkx,iend)
real(r8) :: diff

g=gravity

uflx_dp(0:mkx,:)=uflx_d(mkx:0:-1,:)
vflx_dp(0:mkx,:)=vflx_d(mkx:0:-1,:)
uflx_sh(0:mkx,:)=uflx_s(mkx:0:-1,:)
vflx_sh(0:mkx,:)=vflx_s(mkx:0:-1,:)
uten_sh=uten_s
vten_sh=vten_s
uten_dp=uten_d
vten_dp=vten_d
sten_sh=sten_s


if (isdeep) then
	uflx=uflx_dp
	vflx=vflx_dp
	uten=uten_dp
	vten=vten_dp
else

    uflx=uflx_sh
	vflx=vflx_sh
	uten=uten_sh
	vten=vten_sh
endif

!!shallow

uflx_tot=uflx_sh+uflx_dp
vflx_tot=vflx_sh+vflx_dp
uten_tot=uten_sh+uten_dp
vten_tot=vten_sh+vten_dp

 do i=1,iend
     do k=1,int(cnt_sh(i))
       if (uten_sh(k,i) > 10._r8 .or. vten_sh(k,i) > 10._r8) then
         print*, "shallow error! Rank=", mpi_rank()
         print*, "k=",k,"uten_sh=",uten_sh(k,i),"vten_sh=",vten_sh(k,i)
         print*, "uflx_up=",uflx_sh(k,i),"uflx_dw=",uflx_sh(k-1,i)
         print*, "vflx_up=",vflx_sh(k,i),"vflx_dw=",vflx_sh(k-1,i)
       endif
     enddo
   enddo

do i=1,iend
	do k=1,int(cnt_sh(i))

		ufp=u0(k+1,i)+uten(k+1,i)*dt
		uf=u0(k,i)+uten(k,i)*dt
	!	ufm=u0(k-1,i)+uten(k-1,i)*dt
		vfp=v0(k+1,i)+vten(k+1,i)*dt
		vf=v0(k,i)+vten(k,i)*dt
	!	vfm=v0(k-1,i)+vten(k-1,i)*dt

	 if (k .eq. 1) then

		origin= g / 4._r8 / dp0(k,i) * (  &
                 uflx(k,i)*(ufp - uf + u0(k+1,i) - u0(k,i)) +   &
                 vflx(k,i)*(vfp - vf + v0(k+1,i) - v0(k,i)))

	 elseif( k .ge. 2 .and. k .le. int(cnt_sh(i))-1 ) then
                ufm=u0(k-1,i)+uten(k-1,i)*dt
                vfm=v0(k-1,i)+vten(k-1,i)*dt
	 	origin= g / 4._r8 / dp0(k,i) * (  &
                 uflx(k,i)*(ufp - uf + u0(k+1,i) - u0(k,i)) +     &
                 uflx(k-1,i)*(uf - ufm + u0(k,i) - u0(k-1,i)) +   &
                 vflx(k,i)*(vfp - vf + v0(k+1,i) - v0(k,i)) +     &
                 vflx(k-1,i)*(vf - vfm + v0(k,i) - v0(k-1,i)))
     
     elseif( k .eq. int(cnt_sh(i)) ) then
        ufm=u0(k-1,i)+uten(k-1,i)*dt
        vfm=v0(k-1,i)+vten(k-1,i)*dt
     	origin= g / 4._r8 / dp0(k,i) * (  &
     		     uflx(k-1,i)*(uf - ufm + u0(k,i) - u0(k-1,i)) +   &
                 vflx(k-1,i)*(vf - vfm + v0(k,i) - v0(k-1,i)))
     endif

	    ufp=u0(k+1,i)+uten_tot(k+1,i)*dt
            uf=u0(k,i)+uten_tot(k,i)*dt
	   ! ufm=u0(k-1,i)+uten_tot(k-1,i)*dt
	    vfp=v0(k+1,i)+vten_tot(k+1,i)*dt
	    vf=v0(k,i)+vten_tot(k,i)*dt
	   ! vfm=v0(k-1,i)+vten_tot(k-1,i)*dt

	if (k .eq. 1) then

		new= g / 4._r8 / dp0(k,i) * (  &
                 uflx(k,i)*(ufp - uf + u0(k+1,i) - u0(k,i)) +   &
                 vflx(k,i)*(vfp - vf + v0(k+1,i) - v0(k,i)))

	 elseif( k .ge. 2 .and. k .le. int(cnt_sh(i))-1 ) then
                ufm=u0(k-1,i)+uten_tot(k-1,i)*dt
                vfm=v0(k-1,i)+vten_tot(k-1,i)*dt
	 	new= g / 4._r8 / dp0(k,i) * (  &
                 uflx(k,i)*(ufp - uf + u0(k+1,i) - u0(k,i)) +     &
                 uflx(k-1,i)*(uf - ufm + u0(k,i) - u0(k-1,i)) +   &
                 vflx(k,i)*(vfp - vf + v0(k+1,i) - v0(k,i)) +     &
                 vflx(k-1,i)*(vf - vfm + v0(k,i) - v0(k-1,i)))
     
     elseif( k .eq. int(cnt_sh(i)) ) then
        ufm=u0(k-1,i)+uten_tot(k-1,i)*dt
        vfm=v0(k-1,i)+vten_tot(k-1,i)*dt
     	new= g / 4._r8 / dp0(k,i) * (  &
     		     uflx(k-1,i)*(uf - ufm + u0(k,i) - u0(k-1,i)) +   &
                 vflx(k-1,i)*(vf - vfm + v0(k,i) - v0(k-1,i)))
     endif

    diff=new-origin
    sten_sh(k,i)=sten_sh(k,i)-diff
enddo
enddo

sten_s(:mkx,:)=sten_sh(:mkx,:)
end subroutine fix

end module fix_mod
