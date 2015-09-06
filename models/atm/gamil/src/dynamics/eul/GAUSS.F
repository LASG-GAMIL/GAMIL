c
        subroutine gauss(a,b,nx)
      use pmgrid, only: plat
      use mpi_gamil

c
        implicit none
        integer nx,i,j,k
        real*8 a(nx,jbeg0:jend0,5),b(nx,jbeg0:jend0),y
        real*8 send_array(jbeg0:jend0,8), recv_array(jbeg0:jend0,8)
c

      if (iend0 .eq. nx) then
        do j=jbeg1,jend1
          send_array(j,1) = b(nx-3,j)
          send_array(j,2) = b(nx-1,j)
          send_array(j,3) = b(nx,j)
        end do
        call gamil_array_send_right(send_array,recv_array,(jend0-jbeg0+1)*3)
      end if
      
      if (ibeg0 .eq. 1) then
        call gamil_array_recv_left(send_array,recv_array,(jend0-jbeg0+1)*3)
        do j=jbeg1,jend1
          b(nx-3,j) = recv_array(j,1)
          b(nx-1,j) = recv_array(j,2)
          b(nx,j) = recv_array(j,3)
        end do
      end if

      do j=jbeg1,jend1
        b(2,j)=b(nx-3,j)
        b(nx-3,j)=0.0
      end do

      if (ibeg0 .ne. 1) then
        call gamil_array_recv_left(send_array,recv_array,(jend0-jbeg0+1)*8)
        do j=jbeg1,jend1
          do i = 1,5
            a(ibeg2-1,j,i) = recv_array(j,i)
          end do
          b(ibeg2-1,j) = recv_array(j,6)
          b(nx-1,j) = recv_array(j,7)
          b(nx,j) = recv_array(j,8)
        end do
      end if

!$OMP PARALLEL DO PRIVATE (j,i,y)
      do j=jbeg1,jend1
        do i=ibeg2,iend2
           y=a(i,j,1)/a(i-1,j,2)
           a(i,j,1)=0.0
           a(i,j,2)=a(i,j,2)-a(i-1,j,3)*y
           a(i,j,4)=a(i,j,4)-a(i-1,j,4)*y
           a(i,j,5)=a(i,j,5)-a(i-1,j,5)*y
c
           y=b(i-1,j)/a(i-1,j,2)
           b(i,j)=b(i,j)-a(i-1,j,3)*y
           b(nx-1,j)=b(nx-1,j)-a(i-1,j,4)*y
           b(nx,j)=b(nx,j)-a(i-1,j,5)*y
        end do
      end do

      if (iend0 .ne. nx) then
        do j=jbeg1,jend1
          do i = 1,5
            send_array(j,i) = a(iend2,j,i)
          end do
          send_array(j,6) = b(iend2,j)
          send_array(j,7) = b(nx-1,j)
          send_array(j,8) = b(nx,j)
        end do
        call gamil_array_send_right(send_array,recv_array,(jend0-jbeg0+1)*8)
      end if

      if (iend0 .eq. nx) then
!$OMP PARALLEL DO PRIVATE (j,i,y)
        do j=jbeg1,jend1
          i=nx-2
          y=b(i,j)/a(i,j,2)
          b(nx-1,j)=b(nx-1,j)-a(i,j,4)*y
          b(nx,j)=b(nx,j)-a(i,j,5)*y
          b(nx,j)=b(nx,j)/b(nx-1,j)
          a(nx-2,j,5)=(a(nx-2,j,5)-a(nx-2,j,4)*b(nx,j))/a(nx-2,j,2)
      end do
      end if

      if (iend0 .ne. nx) then
        call gamil_array_recv_right(send_array,recv_array,(jend0-jbeg0+1)*6)
        do j=jbeg1,jend1
          do i = 1,5
            a(iend2,j,i) = recv_array(j,i)
          end do
          b(nx,j) = recv_array(j,6)
        end do
      end if

!$OMP PARALLEL DO PRIVATE (j,i,y)
      do j=jbeg1,jend1
        do i=iend2,ibeg2,-1
           y=a(i-1,j,4)*b(nx,j)+a(i-1,j,3)*a(i,j,5)
           a(i-1,j,5)=(a(i-1,j,5)-y)/a(i-1,j,2)
        end do
      end do

      if (ibeg0 .ne. 1) then
        do j=jbeg1,jend1
          do i = 1,5
            send_array(j,i) = a(ibeg2-1,j,i)
          end do
          send_array(j,6) = b(nx,j)
        end do
        call gamil_array_send_left(send_array,recv_array,(jend0-jbeg0+1)*6)
      end if
c
      return
      end
