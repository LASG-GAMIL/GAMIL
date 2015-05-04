#include <misc.h>
#include <params.h>

#if (defined SPMD)

!! (wh 2003.10.27)
!!---------------------

      SUBROUTINE SETMSA0(TBB,HBB,HS,P0,T0,PSB,TSB)

      IMPLICIT NONE

#include <PARADYN>

#include <mpif.h>
#include <commpi.h>

	REAL*8
     _       TBB (NA   )     ! INPUT , THE TEMPERATURE OF STANDARD ATMOSPHERE
     _      ,HBB (NA   )     ! INPUT , THE HEIGHT OF THE STANDARD ATMOSPHERE
     _      ,HS  (NX,NY)     ! INPUT , GEOPOTENTIAL OF THE GROUND ELEVATION
     _      ,P0              ! INPUT , PRESSURE AT THE SEA LEVEL
     _      ,T0              ! INPUT , TEMPERATURE AT THE SEA LEVEL
     _      ,PSB (NX,NY)     ! OUTPUT, THE SURFACE PRESSURE OF THE ST. ATMS.
     _      ,TSB (NX,NY)     ! OUTPUT, THE SURFACE TEMPERATURE OF THE STD. ATMS.

      character*50 filename
      INTEGER I,J,KPP
      REAL*8  XP,YP,WK5,WK6,HS1,EPSL
      REAL*8  SDH
      EXTERNAL SDH

	EPSL = 0.01D0
!
!     CALCULATING THE SURFACE PRESSURE OF STANDARD ATMOSPHERE PSB.
!

      DO J=jbeg,jend
	  DO I=2,NX-1
          HS1=HS(I,J)
          XP=P0*(1.0D0-HS1/(RD*T0))
!
!         XP IS THE FIRST GUESS GIVEN BY USING PRESSURE-HEIGHT FORMULA.
!
          CALL NLS(XP,SDH,EPSL,HS1,HBB)
          PSB(I,J)=XP
        ENDDO
	  PSB(1,J)=PSB(NX-1,J)
	  PSB(NX,J)=PSB(2,J)
!
!     CALCULATING THE SURFACE TEMPERATURE OF STANDARD ATMOSPHERE TSB.
!
        DO I = 2,NX-1
          WK5=PSB(I,J)/DPALIB
          KPP=INT(WK5)
          WK6=WK5-DFLOAT(KPP)
          TSB(I,J)=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
        ENDDO
	  TSB(1,J)=TSB(NX-1,J)
	  TSB(NX,J)=TSB(2,J)
      ENDDO
!
!	CALCULATE PSB,TSB AT THE PLOES
!
	if (inc_pole) then
         HS1=HS(2,jpole)
         XP=P0*(1.0D0-HS1/(RD*T0))
         CALL NLS(XP,SDH,EPSL,HS1,HBB)
!
         WK5=XP/DPALIB
         KPP=INT(WK5)
         WK6=WK5-DFLOAT(KPP)
         YP=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
!
	  DO I=1,NX
          PSB(I,jpole)=XP
	    TSB(I,jpole)=YP
	  END DO
       endif


!- check ---------------------------------------------------------------
!
!      write(filename,12) 'setmsa0-p-',myrank,'.out'
!12    format(a10,i1,a4)
!
!      open (10,file=trim(filename))
!
!        write(10,*) '----------------- psb -----------------'
!        write(10,11) (j,(psb(i,j),i=1,3),j=1,ny)
!
!        write(10,*) '----------------- tsb -----------------'
!        write(10,11) (j,(tsb(i,j),i=1,3),j=1,ny)
!
!11      format(1x,i5,3e30.20)
!      close(10)
!!      call mpi_finalize(j)
!!      stop 'setmsa0'
!
!!-----------------------------------------------------------------------
!
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
!
#include <PARADYN>
!
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

#else
!
      SUBROUTINE SETMSA0(TBB,HBB,HS,P0,T0,PSB,TSB)
!
      IMPLICIT NONE
!
#include <PARADYN>
!
	REAL*8
     _       TBB (NA   )     ! INPUT , THE TEMPERATURE OF STANDARD ATMOSPHERE
     _      ,HBB (NA   )     ! INPUT , THE HEIGHT OF THE STANDARD ATMOSPHERE
     _      ,HS  (NX,NY)     ! INPUT , GEOPOTENTIAL OF THE GROUND ELEVATION
     _      ,P0              ! INPUT , PRESSURE AT THE SEA LEVEL
     _      ,T0              ! INPUT , TEMPERATURE AT THE SEA LEVEL
     _      ,PSB (NX,NY)     ! OUTPUT, THE SURFACE PRESSURE OF THE ST. ATMS.
     _      ,TSB (NX,NY)     ! OUTPUT, THE SURFACE TEMPERATURE OF THE STD. ATMS.
!
      INTEGER I,J,KPP
      REAL*8  XP,YP,WK5,WK6,HS1,EPSL
      REAL*8  SDH
      EXTERNAL SDH
!
	EPSL = 0.01D0
!
!     CALCULATING THE SURFACE PRESSURE OF STANDARD ATMOSPHERE PSB.
!
      DO J=2,NY-1
	  DO I=2,NX-1
          HS1=HS(I,J)
          XP=P0*(1.0D0-HS1/(RD*T0))
!
!         XP IS THE FIRST GUESS GIVEN BY USING PRESSURE-HEIGHT FORMULA.
!
          CALL NLS(XP,SDH,EPSL,HS1,HBB)
          PSB(I,J)=XP
        ENDDO
	  PSB(1,J)=PSB(NX-1,J)
	  PSB(NX,J)=PSB(2,J)
!
!     CALCULATING THE SURFACE TEMPERATURE OF STANDARD ATMOSPHERE TSB.
!
        DO I = 2,NX-1
          WK5=PSB(I,J)/DPALIB
          KPP=INT(WK5)
          WK6=WK5-DFLOAT(KPP)
          TSB(I,J)=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
        ENDDO
	  TSB(1,J)=TSB(NX-1,J)
	  TSB(NX,J)=TSB(2,J)
      ENDDO
!
!	CALCULATE PSB,TSB AT THE PLOES
!
	DO J=1,NY,NY-1
	  HS1=HS(2,J)
        XP=P0*(1.0D0-HS1/(RD*T0))
        CALL NLS(XP,SDH,EPSL,HS1,HBB)
!
        WK5=XP/DPALIB
        KPP=INT(WK5)
        WK6=WK5-DFLOAT(KPP)
        YP=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
!
	  DO I=1,NX
          PSB(I,J)=XP
	    TSB(I,J)=YP
	  END DO
	END DO

!- check -------------------------------------------------------------
!
!       open (10,file='setmsa0-s.out')
!
!        write(10,*) '----------------- psb -----------------'
!        write(10,11) (j,(psb(i,j),i=1,3),j=1,ny)
!
!        write(10,*) '----------------- tsb -----------------'
!        write(10,11) (j,(tsb(i,j),i=1,3),j=1,ny)
!
!11      format(1x,i3,3e30.20)
!      close(10)
!
!--------------------------------------------------------------(wh)---
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

#endif
