#include <misc.h>
#include <params.h>

      SUBROUTINE QPDATA(QT, U0, V0, W0, U, V, WS, DT, DX, DY,
     $                  SINU, SINV, WTGU, WTGV, DSIG)

      use pmgrid, only: beglatexdyn, endlatexdyn, plat
      use mpi_gamil

      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

      REAL*8 U0(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 V0(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 W0(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 U (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 V (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 WS(beglonex:endlonex,beglatexdyn:endlatexdyn,NZ)
      REAL*8 QT(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 DQ(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 UQ(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 VQ(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 WQ(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 DT, DX, DY
      REAL*8 SINU(beglatexdyn:endlatexdyn)
      REAL*8 SINV(beglatexdyn:endlatexdyn)
      REAL*8 WTGU(beglatexdyn:endlatexdyn)
      REAL*8 WTGV(beglatexdyn:endlatexdyn)
      REAL*8 DSIG(NL)
      INTEGER I, J, K

!     For debug.
      REAL*8 QM
!
!     GET THE ADVECTION VELOCITY
!
      DO K = 1, NL
        DO J = jbeg0, jend0
          DO I = beglonex, endlonex
            UQ(I,J,K) = 0.5D0*(U (I,J,K)+U0(I,J,K))
            VQ(I,J,K) = 0.5D0*(V (I,J,K)+V0(I,J,K))
            WQ(I,J,K) = 0.5D0*(WS(I,J,K)+W0(I,J,K))
            DQ(I,J,K) = QT(I,J,K)
          END DO
        END DO
      END DO
!
!     PERFORM HORIZONTAL ADVECTION IN SPHERICAL GEOMETRY
!
!     DO THE 2-D ADVECTION
!
      CALL TSPAS(QT, UQ, VQ, DT, DX, DY, SINU, SINV, WTGU, WTGV)
!
!     PERFORM THE VERTICAL ADVECTION
!       BY  R.C.Yu
      CALL TSPASW(QT, WQ, DT, DSIG)
!
!     PERFORM VERTICAL REDISTRIBUTION TO AVOID NEGATIVE Q-H2O
!
      CALL AVNEGQ(QT, DSIG)

!     For debug.
#ifdef DEBUG
      CALL QEM(QT, DX, DY, DSIG, SINU, SINV, WTGU, WTGV, QM)
      if (is_rootproc) write(6, *) "QM =", QM
#endif

      RETURN
      END
