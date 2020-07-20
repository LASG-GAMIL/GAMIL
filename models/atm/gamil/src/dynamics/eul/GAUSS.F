      SUBROUTINE GAUSS(A, B, NX)

      USE PMGRID, ONLY: PLAT
      USE MPI_GAMIL

      IMPLICIT NONE

      INTEGER NX, I, J, K
      REAL(8) A(NX,JBEG0:JEND0,5)
      REAL(8) B(NX,JBEG0:JEND0)
      REAL(8) Y
      REAL(8) SEND_ARRAY(JBEG0:JEND0,8)
      REAL(8) RECV_ARRAY(JBEG0:JEND0,8)

      IF (IEND0 == NX) THEN
        DO J = JBEG1, JEND1
          SEND_ARRAY(J,1) = B(NX-3,J)
          SEND_ARRAY(J,2) = B(NX-1,J)
          SEND_ARRAY(J,3) = B(NX  ,J)
        END DO
        CALL GAMIL_ARRAY_SEND_RIGHT(SEND_ARRAY, RECV_ARRAY, (JEND0 - JBEG0 + 1) * 3)
      END IF
      
      IF (IBEG0 == 1) THEN
        CALL GAMIL_ARRAY_RECV_LEFT(SEND_ARRAY, RECV_ARRAY, (JEND0 - JBEG0 + 1) * 3)
        DO J = JBEG1, JEND1
          B(NX-3,J) = RECV_ARRAY(J,1)
          B(NX-1,J) = RECV_ARRAY(J,2)
          B(NX  ,J) = RECV_ARRAY(J,3)
        END DO
      END IF

      DO J = JBEG1, JEND1
        B(2   ,J) = B(NX-3,J)
        B(NX-3,J) = 0.0
      END DO

      IF (IBEG0 /= 1) THEN
        CALL GAMIL_ARRAY_RECV_LEFT(SEND_ARRAY, RECV_ARRAY, (JEND0 - JBEG0 + 1) * 8)
        DO J = JBEG1, JEND1
          DO I = 1, 5
            A(IBEG2-1,J,I) = RECV_ARRAY(J,I)
          END DO
          B(IBEG2-1,J) = RECV_ARRAY(J,6)
          B(NX-1   ,J) = RECV_ARRAY(J,7)
          B(NX     ,J) = RECV_ARRAY(J,8)
        END DO
      END IF

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE (Y)
      DO J = JBEG1, JEND1
        DO I = IBEG2, IEND2
          Y = A(I,J,1) / A(I-1,J,2)
          A(I,J,1) = 0.0
          A(I,J,2) = A(I,J,2) - A(I-1,J,3) * Y
          A(I,J,4) = A(I,J,4) - A(I-1,J,4) * Y
          A(I,J,5) = A(I,J,5) - A(I-1,J,5) * Y

          Y = B(I-1,J) / A(I-1,J,2)
          B(I   ,J) = B(I   ,J) - A(I-1,J,3) * Y
          B(NX-1,J) = B(NX-1,J) - A(I-1,J,4) * Y
          B(NX  ,J) = B(NX  ,J) - A(I-1,J,5) * Y
        END DO
      END DO

      IF (IEND0 /= NX) THEN
        DO J = JBEG1, JEND1
          DO I = 1, 5
            SEND_ARRAY(J,I) = A(IEND2,J,I)
          END DO
          SEND_ARRAY(J,6) = B(IEND2,J)
          SEND_ARRAY(J,7) = B(NX-1 ,J)
          SEND_ARRAY(J,8) = B(NX   ,J)
        END DO
        CALL GAMIL_ARRAY_SEND_RIGHT(SEND_ARRAY, RECV_ARRAY, (JEND0 - JBEG0 + 1) * 8)
      END IF

      IF (IEND0 == NX) THEN
!$OMP PARALLEL DO PRIVATE (I, Y)
        DO J = JBEG1, JEND1
          I = NX - 2
          Y = B(I,J) / A(I,J,2)
          B(NX-1,J) = B(NX-1,J) - A(I,J,4) * Y
          B(NX  ,J) = B(NX  ,J) - A(I,J,5) * Y
          B(NX  ,J) = B(NX  ,J) / B(NX-1,J)
          A(NX-2,J,5) = (A(NX-2,J,5) - A(NX-2,J,4) * B(NX,J)) / A(NX-2,J,2)
        END DO
      END IF

      IF (IEND0 /= NX) THEN
        CALL GAMIL_ARRAY_RECV_RIGHT(SEND_ARRAY, RECV_ARRAY, (JEND0 - JBEG0 + 1) * 6)
        DO J = JBEG1, JEND1
          DO I = 1, 5
            A(IEND2,J,I) = RECV_ARRAY(J,I)
          END DO
          B(NX,J) = RECV_ARRAY(J,6)
        END DO
      END IF

!$OMP PARALLEL DO PRIVATE (I, Y)
      DO J = JBEG1, JEND1
        DO I = IEND2, IBEG2, -1
          Y = A(I-1,J,4) * B(NX,J) + A(I-1,J,3) * A(I,J,5)
          A(I-1,J,5) = (A(I-1,J,5) - Y) / A(I-1,J,2)
        END DO
      END DO

      IF (IBEG0 /= 1) THEN
        DO J = JBEG1, JEND1
          DO I = 1, 5
            SEND_ARRAY(J,I) = A(IBEG2-1,J,I)
          END DO
          SEND_ARRAY(J,6) = B(NX,J)
        END DO
        CALL GAMIL_ARRAY_SEND_LEFT(SEND_ARRAY, RECV_ARRAY, (JEND0 - JBEG0 + 1) * 6)
      END IF

      END
