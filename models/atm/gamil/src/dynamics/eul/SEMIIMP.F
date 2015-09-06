#include <misc.h>
#include <params.h>
!
        subroutine semiu(uuk,hhk,p,du,dps1
     _                  ,dus,dps0,dtdy,oux,dsig,c0)
c
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
        implicit none

#include <PARADYN>

        integer i,j,i1,i2,k
        real*8 uuk(beglonex:endlonex,beglatexdyn:endlatexdyn)
        real*8 hhk(beglonex:endlonex,beglatexdyn:endlatexdyn)
        real*8 p(beglonex:endlonex,beglatexdyn:endlatexdyn)
        real*8 du(beglonex:endlonex,beglatexdyn:endlatexdyn,nl)
        real*8 dps1(beglonex:endlonex,beglatexdyn:endlatexdyn)
        real*8 dus(beglonex:endlonex,beglatexdyn:endlatexdyn)
        real*8 dps0(beglonex:endlonex,beglatexdyn:endlatexdyn)
        real*8 c0(beglonex:endlonex,beglatexdyn:endlatexdyn)
        real*8 oux(beglatexdyn:endlatexdyn)
        real*8 dsig(nl)
        real*8 dtdy
        real*8 ai(beglonex:endlonex,jbeg0:jend0),
     _         bi(beglonex:endlonex,jbeg0:jend0),
     _         di(beglonex:endlonex,jbeg0:jend0)
        real*8 dx,us(beglonex:endlonex,jbeg0:jend0),duss(beglonex:endlonex,jbeg0:jend0)
        real*8 al(nx,jbeg0:jend0,5),ak(nx,jbeg0:jend0),wh(beglonex:endlonex,jbeg0:jend0)
        integer comm_request_id1, comm_request_id2
c
c       Semi-implicit scheme
c
c

        call register_comm_array(beglonex,endlonex,jbeg0,jend0,1,1,1,1,ai(:,jbeg0))
        call register_comm_array(beglonex,endlonex,jbeg0,jend0,1,1,1,1,bi(:,jbeg0))
        call register_comm_array(beglonex,endlonex,jbeg0,jend0,1,1,1,1,wh(:,jbeg0))
        call gamil_arrays_comm(COMM_TO_LEFT,1,p(:,beglatexdyn),uuk(:,beglatexdyn),request_id=comm_request_id1)
        call gamil_arrays_comm(COMM_TO_RIGHT,1,p(:,beglatexdyn),c0(:,beglatexdyn),uuk(:,beglatexdyn),request_id=comm_request_id2)

!$OMP PARALLEL DO PRIVATE (I,J,K)
        do j=jbeg1,jend1
          do i=ibeg1,iend1
            duss(i,j)=0.0d0
            do k=1,nl
	          duss(i,j)=duss(i,j)+du(i,j,k)*dsig(k)
            end do
          end do
        end do

        call wait_icomm_request(comm_request_id1)
        call wait_icomm_request(comm_request_id2)

!$OMP PARALLEL DO PRIVATE (I,J,dx)
        do j=jbeg1,jend1
          dx=0.25d0*oux(j)*dtdy
          do i=ibeg1,iend1
	      ai(i,j)=dx*c0(i,j)*(p(i+1,j)+p(i,j))
	      bi(i,j)=dx*c0(i,j)*(p(i,j)+p(i-1,j))
	      di(i,j)=dx*c0(i-1,j)*(p(i,j)+p(i-1,j))
              us(i,j)=uuk(i,j)+dtdy*duss(i,j)
     _                -bi(i,j)*hhk(i,j)+di(i,j)*hhk(i-1,j)
              wh(i,j)=hhk(i,j)+dtdy*dps1(i,j)*c0(i,j)
     _                -ai(i,j)*uuk(i+1,j)+bi(i,j)*uuk(i,j)
          end do
       end do
c
        call gamil_arrays_comm(COMM_TO_RIGHT,1,ai(:,jbeg0),bi(:,jbeg0),wh(:,jbeg0))

!$OMP PARALLEL DO PRIVATE (I,J,K,dx,i1,i2)
        do j=jbeg1,jend1
           do i=ibeg1,iend1
              ak(i,j) = 0.0
              al(i,j,:) = 0.0
              i1=i-1
              i2=i+1
              if (i1.eq.1) i1=nx-1
              if (i2.eq.nx) i2=2
              if (i.lt.nx-1) then
                 al(i,j,1)=-bi(i-1,j)*di(i,j)
                 al(i,j,2)=1.0d0+bi(i,j)*bi(i,j)+ai(i-1,j)*di(i,j)
                 al(i,j,3)=-ai(i,j)*bi(i,j)
                 al(i,j,4)=0.0
                 al(i,j,5)=us(i,j)-bi(i,j)*wh(i,j)+di(i,j)*wh(i-1,j)
              else
                 ak(i,j)= 1.0d0+bi(i,j)*bi(i,j)+ai(i-1,j)*di(i,j)
                 ak(i1,j)=-bi(i-1,j)*di(i,j)
                 ak(i2,j)=-ai(i,j)*bi(i,j)
                 ak(nx-3,j) = ak(2,j)
                 ak(nx,j)=us(i,j)-bi(i,j)*wh(i,j)+di(i,j)*wh(i-1,j)
              end if
           end do
!
           al(2,j,4)=al(2,j,1)
           al(2,j,1)=0.0
           al(nx-2,j,4)=al(nx-2,j,3)
           al(nx-2,j,3)=0.0
        end do

        call t_startf("GAUSS")
        call gauss(al,ak,nx)

        do j=jbeg1,jend1  
           al(nx-1,j,5)=ak(nx,j)
        enddo
        call t_stopf("GAUSS")
        
 !$OMP PARALLEL DO PRIVATE (I,J)
        do j=jbeg1,jend1  
	    do i=ibeg1,iend1
	     dus(i,j)=al(i,j,5)
	    end do
        end do

        call gamil_arrays_comm(COMM_TO_LEFT,1,dus(:,beglatexdyn))

 !$OMP PARALLEL DO PRIVATE (I,J)
        do j=jbeg1,jend1  
	   do i=ibeg1,iend1
	     dps0(i,j)=wh(i,j)-ai(i,j)*dus(i+1,j)+bi(i,j)*dus(i,j)
	     dps0(i,j)=(dps0(i,j)-hhk(i,j))/(dtdy*c0(i,j))-dps1(i,j)
	   end do
c
	   do i=ibeg1,iend1
	     dus(i,j)=(dus(i,j)-uuk(i,j))/dtdy-duss(i,j)
	   end do
        enddo

      do j=jbeg0,jend0
       if(j .eq. 1 .or. j .eq. plat) then
        do i=ibeg0,iend0
           dus(i,j)=0.0
           dps0(i,j)=0.0
        end do
       endif
      end do        

      call gamil_arrays_comm(COMM_TO_LEFT,1,dus(:,beglatexdyn),dps0(:,beglatexdyn))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,dus(:,beglatexdyn),dps0(:,beglatexdyn))

      call remove_comm_array(wh(:,jbeg0))
      call remove_comm_array(bi(:,jbeg0))
      call remove_comm_array(ai(:,jbeg0))

c
c       End of semi-implicit quadratic-conservation scheme
c
        return
        end
