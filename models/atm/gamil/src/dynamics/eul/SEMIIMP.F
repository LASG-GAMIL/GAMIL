#include <misc.h>
#include <params.h>

      SUBROUTINE SEMIU(UUK, HHK, P, DU, DPS1, DUS, DPS0, DTDY, OUX, DSIG, C0)

      USE PMGRID, ONLY: BEGLATEXDYN, ENDLATEXDYN
      USE MPI_GAMIL

      IMPLICIT NONE

#include <PARADYN>

      INTEGER I, J, I1, I2, K
      REAL*8 UUK(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 HHK(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 P(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 DU(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8 DPS1(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 DUS(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 DPS0(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 C0(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 OUX(BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 DSIG(NL)
      REAL*8 DTDY, DX
      REAL*8 AI(BEGLONEX:ENDLONEX,JBEG0:JEND0)
      REAL*8 BI(BEGLONEX:ENDLONEX,JBEG0:JEND0)
      REAL*8 DI(BEGLONEX:ENDLONEX,JBEG0:JEND0)
      REAL*8 US(BEGLONEX:ENDLONEX,JBEG0:JEND0)
      REAL*8 DUSS(BEGLONEX:ENDLONEX,JBEG0:JEND0)
      REAL*8 AL(NX,JBEG0:JEND0,5)
      REAL*8 AK(NX,JBEG0:JEND0)
      REAL*8 WH(BEGLONEX:ENDLONEX,JBEG0:JEND0)
      INTEGER COMM_REQUEST_ID1, COMM_REQUEST_ID2

      CALL REGISTER_COMM_ARRAY(BEGLONEX, ENDLONEX, JBEG0, JEND0, 1, 1, 1, 1, AI(:,JBEG0))
      CALL REGISTER_COMM_ARRAY(BEGLONEX, ENDLONEX, JBEG0, JEND0, 1, 1, 1, 1, BI(:,JBEG0))
      CALL REGISTER_COMM_ARRAY(BEGLONEX, ENDLONEX, JBEG0, JEND0, 1, 1, 1, 1, WH(:,JBEG0))
      CALL GAMIL_ARRAYS_COMM(COMM_TO_LEFT,  1, P(:,BEGLATEXDYN), UUK(:,BEGLATEXDYN), REQUEST_ID=COMM_REQUEST_ID1)
      CALL GAMIL_ARRAYS_COMM(COMM_TO_RIGHT, 1, P(:,BEGLATEXDYN), C0(:,BEGLATEXDYN), UUK(:,BEGLATEXDYN), REQUEST_ID=COMM_REQUEST_ID2)

!$OMP PARALLEL DO PRIVATE (I, J, K)
      DO J = JBEG1, JEND1
        DO I = IBEG1, IEND1
          DUSS(I,J) = 0.0D0
          DO K = 1, NL
            DUSS(I,J) = DUSS(I,J) + DU(I,J,K)*DSIG(K)
          END DO
        END DO
      END DO

      CALL WAIT_ICOMM_REQUEST(COMM_REQUEST_ID1)
      CALL WAIT_ICOMM_REQUEST(COMM_REQUEST_ID2)

!$OMP PARALLEL DO PRIVATE (I, J, DX)
      DO J = JBEG1, JEND1
        DX = 0.25D0*OUX(J)*DTDY
        DO I = IBEG1, IEND1
          AI(I,J) = DX*C0(I,J)*(P(I+1,J)+P(I,J))
          BI(I,J) = DX*C0(I,J)*(P(I,J)+P(I-1,J))
          DI(I,J) = DX*C0(I-1,J)*(P(I,J)+P(I-1,J))
          US(I,J) = UUK(I,J)+DTDY*DUSS(I,J)-BI(I,J)*HHK(I,J)+DI(I,J)*HHK(I-1,J)
          WH(I,J) = HHK(I,J)+DTDY*DPS1(I,J)*C0(I,J)-AI(I,J)*UUK(I+1,J)+BI(I,J)*UUK(I,J)
        END DO
      END DO

      CALL GAMIL_ARRAYS_COMM(COMM_TO_RIGHT, 1, AI(:,JBEG0), BI(:,JBEG0), WH(:,JBEG0))

!$OMP PARALLEL DO PRIVATE (I,J,K,DX,I1,I2)
      DO J = JBEG1, JEND1
        DO I = IBEG1, IEND1
          AK(I,J)   = 0.0
          AL(I,J,:) = 0.0
          I1 = I-1
          I2 = I+1
          IF (I1 .EQ. 1) I1 = NX-1
          IF (I2 .EQ. NX) I2 = 2
          IF (I .LT. NX-1) THEN
            AL(I,J,1) = -BI(I-1,J)*DI(I,J)
            AL(I,J,2) = 1.0D0+BI(I,J)*BI(I,J)+AI(I-1,J)*DI(I,J)
            AL(I,J,3) = -AI(I,J)*BI(I,J)
            AL(I,J,4) = 0.0
            AL(I,J,5) = US(I,J)-BI(I,J)*WH(I,J)+DI(I,J)*WH(I-1,J)
          ELSE
            AK(I,J)    = 1.0D0+BI(I,J)*BI(I,J)+AI(I-1,J)*DI(I,J)
            AK(I1,J)   = -BI(I-1,J)*DI(I,J)
            AK(I2,J)   = -AI(I,J)*BI(I,J)
            AK(NX-3,J) = AK(2,J)
            AK(NX,J)   = US(I,J)-BI(I,J)*WH(I,J)+DI(I,J)*WH(I-1,J)
          END IF
        END DO

        AL(2,J,4) = AL(2,J,1)
        AL(2,J,1) = 0.0
        AL(NX-2,J,4) = AL(NX-2,J,3)
        AL(NX-2,J,3) = 0.0
      END DO

      CALL T_STARTF("GAUSS")
      CALL GAUSS(AL,AK,NX)

      DO J = JBEG1,JEND1  
        AL(NX-1,J,5) = AK(NX,J)
      END DO
      CALL T_STOPF("GAUSS")
        
!$OMP PARALLEL DO PRIVATE (I,J)
      DO J = JBEG1,JEND1  
        DO I = IBEG1,IEND1
          DUS(I,J) = AL(I,J,5)
        END DO
      END DO

      CALL GAMIL_ARRAYS_COMM(COMM_TO_LEFT,1,DUS(:,BEGLATEXDYN))

!$OMP PARALLEL DO PRIVATE (I,J)
      DO J = JBEG1,JEND1  
        DO I = IBEG1,IEND1
          DPS0(I,J) = WH(I,J)-AI(I,J)*DUS(I+1,J)+BI(I,J)*DUS(I,J)
          DPS0(I,J) = (DPS0(I,J)-HHK(I,J))/(DTDY*C0(I,J))-DPS1(I,J)
        END DO

        DO I = IBEG1,IEND1
          DUS(I,J) = (DUS(I,J)-UUK(I,J))/DTDY-DUSS(I,J)
        END DO
      END DO

      DO J = JBEG0,JEND0
        IF(J .EQ. 1 .OR. J .EQ. PLAT) THEN
          DO I = IBEG0,IEND0
            DUS (I,J) = 0.0
            DPS0(I,J) = 0.0
          END DO
        END IF
      END DO        

      CALL GAMIL_ARRAYS_COMM(COMM_TO_LEFT,  1, DUS(:,BEGLATEXDYN), DPS0(:,BEGLATEXDYN))
      CALL GAMIL_ARRAYS_COMM(COMM_TO_RIGHT, 1, DUS(:,BEGLATEXDYN), DPS0(:,BEGLATEXDYN))

      CALL REMOVE_COMM_ARRAY(WH(:,JBEG0))
      CALL REMOVE_COMM_ARRAY(BI(:,JBEG0))
      CALL REMOVE_COMM_ARRAY(AI(:,JBEG0))

      RETURN
      END
