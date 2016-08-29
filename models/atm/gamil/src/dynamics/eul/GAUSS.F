
      SUBROUTINE GAUSS(A, B, NX)

      use mpi_gamil

      IMPLICIT NONE
      INTEGER NX, I, J, K
      REAL*8 A(NX,JBEG0:JEND0,5), B(NX,JBEG0:JEND0), Y
      REAL*8 SEND_ARRAY(JBEG0:JEND0,8), RECV_ARRAY(JBEG0:JEND0,8)

      IF (IEND0 .EQ. NX) THEN
        DO J = JBEG1, JEND1
          SEND_ARRAY(J,1) = B(NX-3,J)
          SEND_ARRAY(J,2) = B(NX-1,J)
          SEND_ARRAY(J,3) = B(NX,J)
        END DO
        CALL GAMIL_ARRAY_SEND_RIGHT(SEND_ARRAY, RECV_ARRAY, (JEND0-JBEG0+1)*3)
      END IF
      
      IF (IBEG0 .EQ. 1) THEN
        CALL GAMIL_ARRAY_RECV_LEFT(SEND_ARRAY, RECV_ARRAY, (JEND0-JBEG0+1)*3)
        DO J = JBEG1, JEND1
          B(NX-3,J) = RECV_ARRAY(J,1)
          B(NX-1,J) = RECV_ARRAY(J,2)
          B(NX,J) = RECV_ARRAY(J,3)
        END DO
      END IF

      DO J = JBEG1, JEND1
        B(2,   J) = B(NX-3,J)
        B(NX-3,J) = 0.0
      END DO

      IF (IBEG0 .NE. 1) THEN
        CALL GAMIL_ARRAY_RECV_LEFT(SEND_ARRAY,RECV_ARRAY,(JEND0-JBEG0+1)*8)
        DO J = JBEG1, JEND1
          DO I = 1, 5
            A(IBEG2-1,J,I) = recv_array(j,i)
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

      RETURN
      END
