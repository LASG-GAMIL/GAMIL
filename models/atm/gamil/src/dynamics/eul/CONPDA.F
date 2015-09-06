#include <misc.h>
#include <params.h>


C     ***************************
C     ***************************
      SUBROUTINE CONPDA(DTIMEQ,DX,DY,SINU,WTGV,DSIG,FIRST,ISOR
*    _                 ,IORD,IP,DSNP,DSSP,GC,DTDLT,DTDLN,DTDSG,DXYP)
     _                 ,IORD,IP,DSNP,DSSP,GC,DTDLT,DTDLN,DTDSG)

      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil, only: jbeg0, jend0, jbeg1, jend1

C     ***************************
C     ***************************
!
	IMPLICIT NONE
!
#include <PARADYN>
!
      REAL*8 ZERO,HALF,FOUR
      DATA   ZERO,HALF,FOUR/0.0D0,0.5D0,4.0D0/
!
	  REAL*8 WTGV(beglatexdyn:endlatexdyn),DX,DY,DTIMEQ,FIRST
     _      ,DXYP(beglatexdyn:endlatexdyn),SINU(beglatexdyn:endlatexdyn)
     _      ,DTDLT(beglatexdyn:endlatexdyn),DTDLN(beglatexdyn:endlatexdyn),DTDSG(NL)
     _      ,FIM,DSNP,DSSP,DSIG(NL),GC(beglatexdyn:endlatexdyn)
	INTEGER I,ISOR,IORD,J,JB,JE,K,IP(NX)
!
!
      DO J  = jbeg1,jend1
         DXYP(J) = HALF*(HALF*(RAD*DX)*(SINU(J-1)+SINU(J+1)+SINU(J)*2.0))*(RAD*DY)
      END DO
!
      if (jbeg0 .eq. 1) then
         DXYP(1 )  = (HALF*(RAD*DX)*(SINU(2)+SINU(1)))*(RAD*DY) / FOUR
      endif
      if (jend0 .eq. PLAT) then
         DXYP(PLAT)  = (HALF*(RAD*DX)*(SINU(PLAT-1)+SINU(PLAT)))*(RAD*DY)/FOUR
      endif
!
!     FOR THE USAGE OF SUB.MPDATA & VPDATA
!
      FIM       = DBLE( IM )
!
      IF( FIRST.LE.ZERO ) THEN
        IF( ISOR.EQ.3 ) IORD = MAX0(IORD,3)
        DO I= 1 ,NX
           IP(I)   = MOD( I+IM/2-1,IM ) + 1
	  END DO
	 if (jbeg0 .eq. 1) then
           DSNP    = DXYP(2) / (FIM*DXYP(1))
        endif
        if (jend0 .eq. PLAT) then
           DSSP    = DXYP(PLAT-1) / (FIM*DXYP(PLAT))
        endif
        DO J= jbeg0 ,jend0
           GC(J)   = SINU(J)
	  END DO
      ENDIF
C
      DO J  = jbeg0 ,jend1
         DTDLT(J)  = DTIMEQ / (RAD*DY)*WTGV(J)
	END DO
      DO J  = jbeg1,jend1
         DTDLN(J)  = DTIMEQ / (RAD*SINU(J)*DX)
	END DO
      DO K  = 1 ,NL
         DTDSG(K)  = DTIMEQ / DSIG(K)
      END DO


      RETURN
      END

