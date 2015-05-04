#include <misc.h>
#include <params.h>
!
        subroutine semiu(uu,hh,p,du,dps
     _                  ,up,psp,dt,oux,dsig,c)
c
        implicit none

#include <PARADYN>

#if (defined SPMD)
#include <commpi.h>
#else
        integer jpole
#endif
        integer i,j,i1,i2,i0,k,np
	parameter(np=nx-1)
        real*8
     &   uu(nx,ny),hh(nx,ny),p(nx,ny),up(nx,ny),psp(nx,ny)
     &  ,du(nx,ny,nl),dps(nx,ny),dsig(nl)
     &  ,dt,oux(ny),c(nx,ny)
        real*8 ai(nx),bi(nx),di(nx)
	real*8 dx,us(nx),dus(nx)
	real*8 al(np-2,5),ak(np),wh(nx),byy,xx
c
c       Semi-implicit scheme
c
c
#if (defined SPMD)
        do j=jbeg,jend
#else
        do j=2,ny-1
#endif
           dx=0.25d0*oux(j)*dt
	     do i=2,np
	        dus(i)=0.0d0
	        do k=1,nl
		   dus(i)=dus(i)+du(i,j,k)*dsig(k)
		end do
	     enddo
!
           do i=2,np
              i1=i-1
              i2=i+1
              if (i1.eq.1) i1=np
              if (i2.eq.nx) i2=2
	        ai(i)=dx*c(i,j)*(p(i2,j)+p(i,j))
		bi(i)=dx*c(i,j)*(p(i,j)+p(i1,j))
		di(i)=dx*c(i1,j)*(p(i,j)+p(i1,j))
!
	        us(i)=uu(i,j)+dt*dus(i)
     _			 -bi(i)*hh(i,j)+di(i)*hh(i1,j)
		wh(i)=hh(i,j)+dt*dps(i,j)*c(i,j)
     _                   -ai(i)*uu(i2,j)+bi(i)*uu(i,j)
	     end do
c
           do i=2,np
              i1=i-1
              i2=i+1
              if (i1.eq.1) i1=np
              if (i2.eq.nx) i2=2
              if (i.lt.np) then
                 al(i-1,1)=-bi(i1)*di(i)
                 al(i-1,2)=1.0d0+bi(i)*bi(i)+ai(i1)*di(i)
                 al(i-1,3)=-ai(i)*bi(i)
                 al(i-1,4)=0.0
                 al(i-1,5)=us(i)-bi(i)*wh(i)+di(i)*wh(i1)
              else
                do i0=2,np
                 if (i0.eq.i ) then
                    ak(i0-1)= 1.0d0+bi(i)*bi(i)+ai(i1)*di(i)
                 else if (i0.eq.i1) then
                    ak(i0-1)=-bi(i1)*di(i)
                 else if (i0.eq.i2) then
                    ak(i0-1)=-ai(i)*bi(i)
                 else
                    ak(i0-1)=0.0d0
                 end if
                end do
                ak(np)=us(i)-bi(i)*wh(i)+di(i)*wh(i1)
              end if
           end do
!
           al(1,4)=al(1,1)
           al(1,1)=0.0
           al(np-2,4)=al(np-2,3)
           al(np-2,3)=0.0
           call gauss(al,ak,np-1)
           do i=2,np-1
              up(i,j)=al(i-1,5)
           end do
           i=np
           up(i,j)=ak(np)
	   up(1,j)=ak(np)
	   up(nx,j)=up(2,j)
c
	   do i=2,np
              i1=i-1
              i2=i+1
              if (i1.eq.1) i1=np
              if (i2.eq.nx) i2=2
		    psp(i,j)=wh(i)-ai(i)*up(i2,j)+bi(i)*up(i,j)
		    psp(i,j)=(psp(i,j)-hh(i,j))/(dt*c(i,j))-dps(i,j)
	   end do
	   psp(1,j)=psp(np,j)
	   psp(nx,j)=psp(2,j)
c
	   do i=2,np
              i1=i-1
              i2=i+1
              if (i1.eq.1) i1=np
              if (i2.eq.nx) i2=2
		    up(i,j)=(up(i,j)-uu(i,j))/dt-dus(i)
	   end do
	   up(1,j)=up(np,j)
	   up(nx,j)=up(2,j)
        end do
c
#if (defined SPMD)
      if ( myrank.eq.(nprocs-1) ) then
#else
        jpole = 1
#endif
        do i=1,nx
           up(i,jpole)=0.0
           psp(i,jpole)=0.0
        end do
#if (defined SPMD)
      endif
#endif
!
#if (defined SPMD)
      if ( myrank.eq.0 ) then
#else
        jpole = ny
#endif
        do i=1,nx
           up(i,jpole)=0.0
           psp(i,jpole)=0.0
        end do
#if (defined SPMD)
      endif
#endif
c
c       End of semi-implicit quadratic-conservation scheme
c
        return
        end
