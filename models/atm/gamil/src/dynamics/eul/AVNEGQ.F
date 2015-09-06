#include <misc.h>
#include <params.h>

!!(2003.11.12)
!!--------------------

      SUBROUTINE AVNEGQ(Q,DSGHL)

      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
      
C     **********************
C     **********************
C
C     AVOIDE   NEGATIVE MIXING RATIO
C
      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>


      REAL*8  DSGHL(NL)
      REAL*8  Q(beglonex:endlonex,beglatexdyn:endlatexdyn,NL),QR(NL),QI
      REAL*8 ZERO
      DATA    ZERO / 0.0E0 /
      INTEGER I,J,K


!$OMP PARALLEL DO PRIVATE (I,J,K,QR,QI)
      DO j = jbeg0, jend0
      if ( j .ge. 2 .and. j .le. PLAT-1) then
      DO I  = ibeg1, iend1
      DO K  = 1 ,NL
        QR(K)     = Q(I,J,K)
      ENDDO
      DO K  = 2 ,NL
      QI        = QR(K-1)
      IF( QI.LT.ZERO ) THEN
        QR(K-1)  = ZERO
        QR(K )  = QR(K ) + QI*DSGHL(K)
      ENDIF
      ENDDO
      IF( QR(NL).LT.ZERO ) QR(NL) = ZERO
      DO K  = 1 ,NL
      Q(I,J,K)  = QR(K)
      ENDDO
      ENDDO
      ELSE
      DO K  = 1 ,NL
         QR(K)     = Q(ibeg1,J,K)
      ENDDO
      DO K  = 2 ,NL
      QI        = QR(K-1)
      IF( QI.LT.ZERO ) THEN
        QR(K-1)  = ZERO
        QR(K )  = QR(K ) + QI*DSGHL(K)
      ENDIF
      ENDDO
      IF( QR(NL).LT.ZERO ) QR(NL) = ZERO
      DO K  = 1 ,NL
      DO I  = ibeg1,iend1
      Q(I,J,K)  = QR(K)
      ENDDO
      ENDDO
      ENDIF
      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,Q(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,Q(:,beglatexdyn,1))

      RETURN
      END

