#include <misc.h>
#include <params.h>

      SUBROUTINE TSPASW(Q, W, DT, DSIG)

      use pmgrid, only: beglatexdyn, endlatexdyn, plat
      use mpi_gamil

      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

      REAL*8 Q(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 W(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 DT
      REAL*8 DSIG(NL)

      REAL*8 QWMIN, QWMAX

      REAL*8 WW(NL), FW(NZ), BETAW(NL), QW(NL), QWSTAR(NL)
      REAL*8 WSTAR(NL), AW(NL), HS(NL), HW(NL)
      REAL*8 north_pole_q_nx(NL), south_pole_q_1(NL)

      REAL*8 GAMA, CWSTAR, CW,TEMP1, TEMP2, TEMP3, TEMP4

      REAL*8  ZERO, HALF, FOURTH, EPSM
      DATA ZERO, HALF, FOURTH, EPSM /0.0D0,0.5D0,0.25D0,1.0D-80/
      INTEGER I, J, K, IS, IT

      DO K = 1, NL
        HS(K) = DSIG(K)/DT
      END DO

      DO K = 2, NL
        HW(K) = HALF * (HS(K) + HS(K-1))
      END DO

      DO J = jbeg0, jend0
        DO I = ibeg0, iend0
          DO K = 1, NL
            QW(K) = Q(I,J,K)
            WW(K) = W(I,J,K)
          END DO
          DO K = 2, NL
            FW(K) = HALF*(WW(K)         *(QW(K)+QW(K-1))
     $                   -WW(K)**2/HW(K)*(QW(K)-QW(K-1)))
          END DO
          DO K = 2, NM
            TEMP1 = ABS(WW(K  )/HW(K  ))*(1-ABS(WW(K  )/HW(K  )))
            TEMP2 = ABS(WW(K+1)/HW(K+1))*(1-ABS(WW(K+1)/HW(K+1)))
            GAMA  = MAX(TEMP1, TEMP2)
            BETAW(K) = 2.0D0/(2.0D0-GAMA)
            QWSTAR(K) = QW(K)-BETAW(K)*(FW(K+1)-FW(K))/HS(K)
          END DO
          QWSTAR(1 ) = QW(1 ) - BETAW(2 )*FW(2 )/HS(1 )
          QWSTAR(NL) = QW(NL) + BETAW(NM)*FW(NL)/HS(NL)

          DO K = 1, NL
            QWMIN =  1.0E15
            QWMAX = -1.0E15
            IF (K .EQ. 1) THEN
              QWMIN = MIN(QW(K), QW(K+1), QWMIN)
              QWMAX = MAX(QW(K), QW(K+1), QWMAX)
            ELSE IF (K .EQ. NL) THEN
              QWMIN = MIN(QW(K), QW(K-1), QWMIN)
              QWMAX = MAX(QW(K), QW(K-1), QWMAX)
            ELSE
              QWMIN = MIN(QW(K), QW(K-1), QW(K+1), QWMIN)
              QWMAX = MAX(QW(K), QW(K-1), QW(K+1), QWMAX)
            END IF
            AW(K) = (QWSTAR(K)-QWMIN)*(QWSTAR(K)-QWMAX)
          END DO
          DO K = 2, NL
            TEMP1 = (ABS(AW(K  ))+AW(K  ))/(AW(K  )+EPSM)
            TEMP2 = (ABS(AW(K-1))+AW(K-1))/(AW(K-1)+EPSM)
            TEMP3 = (ABS(AW(K))+AW(K))*(ABS(AW(K-1))+AW(K-1))
            TEMP4 = ABS(AW(K))*ABS(AW(K-1))+EPSM
            CWSTAR = HALF*(TEMP1+TEMP2)-FOURTH*TEMP3/TEMP4
            CW = CWSTAR+(1-CWSTAR)*ABS(WW(K)/HW(K))
            WSTAR(K) = CW*WW(K)
            FW(K) = HALF*(       WW(K) *(QW(K)+QW(K-1))
     $                   -ABS(CW*WW(K))*(QW(K)-QW(K-1)))
          END DO
          FW(1)  = ZERO
          FW(NZ) = ZERO

          DO K = 1, NL
            QW(K) = QW(K)-(FW(K+1)-FW(K))/HS(K)
          END DO

          DO K = 1, NL
            Q(I,J,K) = QW(K)
          END DO
        END DO
      END DO

      RETURN
      END
