#include <misc.h>
#include <params.h>

!!(2003.11.12)
!!--------------------

      SUBROUTINE AVNEGQ(Q,DSGHL)
C     **********************
C     **********************
C
C     AVOIDE   NEGATIVE MIXING RATIO
C
      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

#if (defined SPMD)
#include <commpi.h>
      character*50 filename
#endif

      REAL*8  DSGHL(NL)
      REAL*8  Q(NX,NY,NL),QR(NL),QI
      REAL*8 ZERO
      DATA    ZERO / 0.0E0 /
      INTEGER I,J,K

#if (defined SPMD)
      do 20  j  = jbeg,jend
#else
      DO 20  J  = JB,JE
#endif
      DO 20  I  = IB,IE
      DO 10  K  = 1 ,NL
10    QR(K)     = Q(I,J,K)
      DO 30  K  = 2 ,NL
      QI        = QR(K-1)
      IF( QI.LT.ZERO ) THEN
        QR(K-1)  = ZERO
        QR(K )  = QR(K ) + QI*DSGHL(K)
      ENDIF
30    CONTINUE
      IF( QR(NL).LT.ZERO ) QR(NL) = ZERO
      DO 20  K  = 1 ,NL
      Q(I,J,K)  = QR(K)
20    CONTINUE

#if (defined SPMD)
      if (inc_pole) then
         DO 40  K  = 1 ,NL
40       QR(K)     = Q(IB,Jpole,K)
         DO 60  K  = 2 ,NL
         QI        = QR(K-1)
         IF( QI.LT.ZERO ) THEN
           QR(K-1)  = ZERO
           QR(K )  = QR(K ) + QI*DSGHL(K)
         ENDIF
60       CONTINUE
         IF( QR(NL).LT.ZERO ) QR(NL) = ZERO
         DO 50  K  = 1 ,NL
         DO 50  I  = IB,IE
         Q(I,Jpole,K)  = QR(K)
50       CONTINUE
      endif
#else
      DO 50  J  = 1 ,NY,JE
      DO 40  K  = 1 ,NL
40    QR(K)     = Q(IB,J,K)
      DO 60  K  = 2 ,NL
      QI        = QR(K-1)
      IF( QI.LT.ZERO ) THEN
        QR(K-1)  = ZERO
        QR(K )  = QR(K ) + QI*DSGHL(K)
      ENDIF
60    CONTINUE
      IF( QR(NL).LT.ZERO ) QR(NL) = ZERO
      DO 50  K  = 1 ,NL
      DO 50  I  = IB,IE
      Q(I,J,K)  = QR(K)
50    CONTINUE
#endif
C
      DO 70  K  = 1 ,NL
      DO 70  J  = 1 ,NY
      Q(1 ,J,K) = Q(IE,J,K)
      Q(NX,J,K) = Q(IB,J,K)
70    CONTINUE

!- check ---------------------------------------------------------
!
!#if (defined SPMD)
!      write(filename,14) 'avnegq-p-',myrank,'.out'
!14    format(a9,i1,a4)
!
!      open (10,file=trim(filename))
!#else
!      open (10,file='avnegq-s.out')
!#endif
!
!      do j=1,ny
!        write(10,11) j,q(1,j,10),q(2,j,10)
!      enddo
!
!11    format(1x,i5,2e30.20)
!      close (10)
!
!#if (defined SPMD)
!      call mpi_finalize(j)
!#endif
!      stop 'avnegq'
!!--------------------------------------------------------------


      RETURN
      END

