#include <misc.h>
#include <params.h>

      SUBROUTINE SETMSA0(TBB,HBB,HS,P0,T0,PSB,TSB)
!
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
      IMPLICIT NONE
!
#include <PARADYN>
!
	REAL*8
     _       TBB (NA   )     ! INPUT , THE TEMPERATURE OF STANDARD ATMOSPHERE
     _      ,HBB (NA   )     ! INPUT , THE HEIGHT OF THE STANDARD ATMOSPHERE
     _      ,HS  (beglonex:endlonex,beglatexdyn:endlatexdyn)     ! INPUT , GEOPOTENTIAL OF THE GROUND ELEVATION
     _      ,P0              ! INPUT , PRESSURE AT THE SEA LEVEL
     _      ,T0              ! INPUT , TEMPERATURE AT THE SEA LEVEL
     _      ,PSB (beglonex:endlonex,beglatexdyn:endlatexdyn)     ! OUTPUT, THE SURFACE PRESSURE OF THE ST. ATMS.
     _      ,TSB (beglonex:endlonex,beglatexdyn:endlatexdyn)     ! OUTPUT, THE SURFACE TEMPERATURE OF THE STD. ATMS.
!
      INTEGER I,J,KPP
      REAL*8  XP,YP,WK5,WK6,HS1(1),EPSL
      REAL*8  SDH
      EXTERNAL SDH
!
	EPSL = 0.01D0
!
!     CALCULATING THE SURFACE PRESSURE OF STANDARD ATMOSPHERE PSB.
!
      DO J=jbeg0,jend0
      IF (J .ge. 2 .and. J .le. PLAT-1) then
	  DO I=ibeg1,iend1
          HS1(1)=HS(I,J)
          XP=P0*(1.0D0-HS1(1)/(RD*T0))
!
!         XP IS THE FIRST GUESS GIVEN BY USING PRESSURE-HEIGHT FORMULA.
!
          CALL NLS(XP,SDH,EPSL,HS1(1),HBB)
          PSB(I,J)=XP
        ENDDO
!
!     CALCULATING THE SURFACE TEMPERATURE OF STANDARD ATMOSPHERE TSB.
!
        DO I = ibeg1,iend1
          WK5=PSB(I,J)/DPALIB
          KPP=INT(WK5)
          WK6=WK5-DFLOAT(KPP)
          TSB(I,J)=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
        ENDDO
!
!	CALCULATE PSB,TSB AT THE PLOES
!
      ELSE
        if (beglonex .eq. 1) HS1(1)=HS(2,J)
        call broadcast_lon_data(2, J, HS1,1)
        
        XP=P0*(1.0D0-HS1(1)/(RD*T0))
        CALL NLS(XP,SDH,EPSL,HS1(1),HBB)
!
        WK5=XP/DPALIB
        KPP=INT(WK5)
        WK6=WK5-DFLOAT(KPP)
        YP=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
!
	  DO I=beglonex,endlonex
            PSB(I,J)=XP
	    TSB(I,J)=YP
	  END DO
	ENDIF
	END DO

        call gamil_arrays_comm(COMM_TO_LEFT,1,PSB(:,beglatexdyn),TSB(:,beglatexdyn))
        call gamil_arrays_comm(COMM_TO_RIGHT,1,PSB(:,beglatexdyn),TSB(:,beglatexdyn))

	RETURN
	END
!
      SUBROUTINE NLS(X,F,EPSL,HS1,HBB)
C     F(X) IS A KNOWN FUNCTIONS AND X IS THE ROOT OF F(X) TO BE FOUND.
C     EPSL IS AN ASSIGNNED ERROR BOUND.
C
      IMPLICIT NONE

#include <PARADYN>

      REAL*8 HBB(NA)
      REAL*8  F,ER,X,W,DX,DFDX,EPSL,HS1
      EXTERNAL F
C
C     JUDGING IF THE ITERATION HAS FINISHED.
   10 ER=F(X,HS1,HBB)
      IF(ABS(ER).LT.EPSL) RETURN
C     CALCULATING PARTIAL DERIVATIVES OF F AS WELL AS LAMDA(RLD).
      W=X
      DX=0.1
C     HERE, THE VALUE 0.1(MB) IS ONLY SUITABLE TO SURFACE PRESSURE.
      X=W+DX
      DFDX=(F(X,HS1,HBB)-ER)/DX
      X=W
C     FINDING A NEW TRIAL SOLUTION.
      IF(ABS(DFDX).GT.1.0D-19) GO TO 20
*     WRITE(6,15) X,ER
   15 FORMAT(/5X,'PSB1=',F8.1,3X,'HB(PSB1)-HS1=',E12.5
     *       /5X,'WARNING: X=PSB1 IS A TURNING POINT OF F(X)]'/)
   20 X=X-ER/DFDX
      GO TO 10
C
      END
!
      FUNCTION SDH(P,HS1,HBB)
!     SDH(P)=HB(P)-HS IS A FUNCTION TO BE USED IN FINDING PSB.
!
      IMPLICIT NONE

#include <PARADYN>

      REAL*8 HBB(NA)
      INTEGER KPP
      REAL*8  HS1,P,SDH,WK5,WK6,HB0
!
!     (HS1 IS THE GROUND ELEVATION AT A CERTAIN GRID POINT.)
!
      WK5=P/DPALIB
      KPP=INT(WK5)
      WK6=WK5-DFLOAT(KPP)
      HB0=(1.D0-WK6)*HBB(KPP)+WK6*HBB(KPP+1)
      SDH=HB0-HS1
!
      RETURN
      END

