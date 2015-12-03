#include <misc.h>
#include <params.h>

      SUBROUTINE AVNEGQ(Q, DSIG)

      use pmgrid, only: beglatexdyn, endlatexdyn, plat
      use mpi_gamil
      
C
C     AVOIDE   NEGATIVE MIXING RATIO
C
      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

      REAL*8  DSIG(NL)
      REAL*8  Q(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8  QR(NL), QI
      REAL*8  ZERO
      DATA    ZERO /0.0D0/
      INTEGER I, J, K

      DO J = jbeg0, jend0
        IF ( J .GE. 2 .AND. J .LE. PLAT-1) THEN
          DO I = ibeg1, iend1
            DO K = 1, NL
              QR(K) = Q(I,J,K)
            END DO
            DO K = 2, NL
              QI = QR(K-1)
              IF (QI .LT. ZERO) THEN
                QR(K-1) = ZERO
                QR(K )  = QR(K)+QI*DSIG(K-1)/DSIG(K)
              END IF
            END DO
            IF (QR(NL) .LT. ZERO) QR(NL) = ZERO
            DO K = 1, NL
              Q(I,J,K)  = QR(K)
            END DO
          END DO
        ELSE
          DO K = 1, NL
            QR(K) = Q(ibeg1,J,K)
          END DO
          DO K = 2, NL
            QI = QR(K-1)
            IF (QI .LT. ZERO) THEN
              QR(K-1) = ZERO
              QR(K )  = QR(K) + QI*DSIG(K-1)/DSIG(K)
            END IF
          END DO
          IF (QR(NL) .LT. ZERO ) QR(NL) = ZERO
          DO K = 1, NL
            DO I  = ibeg1,iend1
              Q(I,J,K)  = QR(K)
            END DO
          END DO
        END IF
      END DO

      call gamil_arrays_comm(COMM_TO_LEFT,  1, Q(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT, 1, Q(:,beglatexdyn,1))

      RETURN
      END
