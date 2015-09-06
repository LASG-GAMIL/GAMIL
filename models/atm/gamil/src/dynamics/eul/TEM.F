#include <misc.h>
#include <params.h>

      SUBROUTINE TEM(U, V, PS, P, TT, HS, CB, TB, DX, DY, DSIG,
     _     SINU, SINV, WTGU, WTGV, TE, TM)

      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil

      IMPLICIT NONE

#include <PARADYN>

      REAL*8 U(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 V(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 PS(beglonex:endlonex,beglatexdyn:endlatexdyn)
      REAL*8 P(beglonex:endlonex,beglatexdyn:endlatexdyn)
      REAL*8 TT(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)

      REAL*8 TB(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 CB(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 HS(beglonex:endlonex,beglatexdyn:endlatexdyn)
      REAL*8 WTGU(beglatexdyn:endlatexdyn)
      REAL*8 WTGV(beglatexdyn:endlatexdyn)
      REAL*8 SINU(beglatexdyn:endlatexdyn)
      REAL*8 SINV(beglatexdyn:endlatexdyn),DX,DY,DSIG(NL)

      INTEGER I, J, K
      REAL*8  DZ2, TE, DS0, DSU, DSV, DS,EE,EK,ES,TM
      REAL*16 TMJ(8,jbeg0:jend0),EUJ(8,jbeg0:jend0),EVJ(8,jbeg0:jend0)
      REAL*16 EEJ(8,jbeg0:jend0),ESJ(8,jbeg0:jend0)
      REAL*16 TMP_EE,TMP_EK,TMP_ES,TMP_TM

      INTEGER begj, endj

      DS0 = RAD*RAD*DX*DY

!$OMP PARALLEL DO PRIVATE (I,J,K,DZ2)
      DO J = jbeg0, jend0
         TMJ(1,J) = 0.0D0
         ESJ(1,J) = 0.0D0
         EUJ(1,J) = 0.0D0
         EVJ(1,J) = 0.0D0
         EEJ(1,J) = 0.0D0
         DO I = ibeg1, iend1
            TMJ(1,J) = TMJ(1,J)+PS(I,J)
            ESJ(1,J) = ESJ(1,J)+PS(I,J)*HS(I,J)
         END DO
         DO K = 1, NL
            DZ2 = DSIG(K)*0.5D0
            DO I = ibeg1, iend1
               EUJ(1,J) = EUJ(1,J)+U(I,J,K)*U(I,J,K)*DZ2
               EVJ(1,J) = EVJ(1,J)+V(I,J,K)*V(I,J,K)*DZ2
               EEJ(1,J) = EEJ(1,J)+(CPD*PS(I,J)*TB(I,J,K)
     &              +CB(I,J,K)/CAPA*P(I,J)*TT(I,J,K))*DSIG(K)
            END DO
         END DO
      END DO
     
      TMP_EK = 0.0D0
      TMP_EE = 0.0D0
      TMP_ES = 0.0D0
      TMP_TM = 0.0D0

      DO J = jbeg0, jend0
         IF (J .EQ. 1) THEN
            DS = 0.25D0*SINV(1)*DS0/WTGV(1)
            DSU = DS
            DSV = 4.0D0*DSU
         ELSE IF (J .EQ. PLAT) THEN
            DS = 0.25D0*SINV(J-1)*DS0/WTGV(J-1)
            DSU = DS
            DSV = 0.0D0
         ELSE
            DS = SINU(J)*DS0/WTGU(J)
            DSU = DS
            DSV = SINV(J)*DS0/WTGV(J)
         END IF

         TMP_EK = TMP_EK+DSU*EUJ(1,J)+DSV*EVJ(1,J)
         TMP_EE = TMP_EE+DSU*EEJ(1,J)
         TMP_ES = TMP_ES+DS*ESJ(1,J)
         TMP_TM = TMP_TM+DS*TMJ(1,J)
      END DO

      call gamil_all_reduce(TMP_EE,TMP_EK,TMP_ES,TMP_TM)

      EE = TMP_EE
      EK = TMP_EK
      ES = TMP_ES
      TM = TMP_TM
      
      TE = EK+EE+ES

      RETURN
      END

      SUBROUTINE QEM(QT, DX, DY, DSIG, SINU, SINV, WTGU, WTGV, QM)

      use pmgrid, only: beglatexdyn, endlatexdyn, plat
      use mpi_gamil

      IMPLICIT NONE

#include <PARADYN>

      REAL*8 QT(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 WTGU(beglatexdyn:endlatexdyn)
      REAL*8 WTGV(beglatexdyn:endlatexdyn)
      REAL*8 SINU(beglatexdyn:endlatexdyn)
      REAL*8 SINV(beglatexdyn:endlatexdyn)
      REAL*8 DX, DY, DSIG(NL)

      INTEGER I, J, K
      REAL*8  DS0, DS, QM
      REAL*16 QMJ(8,jbeg0:jend0)
      REAL*16 TMP_QM

      INTEGER begj, endj

      DS0 = RAD*RAD*DX*DY

!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO J = jbeg0, jend0
         QMJ(1,J) = 0.0D0
         DO K = 1, NL
            DO I = ibeg1, iend1
               QMJ(1,J) = QMJ(1,J)+QT(I,J,K)*DSIG(K)
            END DO
         END DO
      END DO
      TMP_QM = 0.0D0
      DO J = jbeg0, jend0
         IF (J .EQ. 1) THEN
            DS = 0.25D0*SINV(1)*DS0/WTGV(1)
         ELSE IF (J .EQ. PLAT) THEN
            DS = 0.25D0*SINV(J-1)*DS0/WTGV(J-1)
         ELSE
            DS = SINU(J)*DS0/WTGU(J)
         END IF
         TMP_QM = TMP_QM+DS*QMJ(1,J)
      END DO
      call gamil_all_reduce(TMP_QM)
      QM = TMP_QM

      RETURN
      END
