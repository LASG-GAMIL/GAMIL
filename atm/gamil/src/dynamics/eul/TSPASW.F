#include <misc.h>
#include <params.h>

C     =======================
      SUBROUTINE TSPASW(Q,W,DTDSG)
C     =======================
C
      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

#if (defined SPMD)
#include <commpi.h>
#endif
C
      REAL*8 Q(NX,NY,NL),W(NX,NY,NL)
      REAL*8 QWMIN(NX,NY,NL),QWMAX(NX,NY,NL)
!
      REAL*8 WW(NL),FW(NZ),BETAW(NL),QW(NL),QWSTAR(NL)
     _      ,WSTAR(NL),AW(NL),HS(NL),HW(NL),DTDSG(NL)
!
      REAL*8 GAMA,CWSTAR,CW,TEMP1,TEMP2,TEMP3,TEMP4
!
      REAL*8  ZERO,HALF,FOURTH,EPSM
      DATA ZERO,HALF,FOURTH,EPSM/ 0.0D0,0.5D0,0.25D0,
     $     1.0D-80/
*    $     1.0E-6/
      INTEGER I,J,K,IS,IT
      integer begj,endj

#if (defined SPMD)
      character*50 filename

      begj=2
      endj=ny-1
#else
      begj=1
      endj=ny
#endif
!
       DO I=1,NX
       DO J=begj,endj
       DO K=1,NL
          QWMIN(I,J,K)=1.0E15
          QWMAX(I,J,K)=-1.0E15
       ENDDO
       ENDDO
       ENDDO
!
      DO K=1,NL
	HS(K)=1.0D0/DTDSG(K)
      ENDDO
C
      DO K=2,NL
	HW(K)=HALF*(HS(K)+HS(K-1))
      ENDDO
C
      DO J=begj,endj

#if (defined SPMF)
!     (IF(J.EQ.1) THEN)
      if (myrank.eq.nprocs-1.and.j.eq.jpole) then
	IS=NX
	IT=NX
!     (ELSE IF(J.EQ.NY) THEN)
      else if (myrank.eq.0.and.j.eq.jpole) then
	IS=1
	IT=1
      ELSE
	IS=1
	IT=NX
      ENDIF
#else
      IF(J.EQ.1) THEN
	IS=NX
	IT=NX
      ELSE IF(J.EQ.NY) THEN
	IS=1
	IT=1
      ELSE
	IS=1
	IT=NX
      ENDIF
#endif
C
      DO I=IS,IT
C
      DO K=1,NL
	QW(K)=Q(I,J,K)
	WW(K)=W(I,J,K)
      ENDDO
C
      DO K=2,NL
	FW(K)=HALF*WW(K)*(QW(K)+QW(K-1))
     $ -HALF*WW(K)*WW(K)*(QW(K)-QW(K-1))/HW(K)
      ENDDO
C
      DO K=2,NM
 	TEMP1=ABS(WW(K)/HW(K))*(1-ABS(WW(K)/HW(K)))
	TEMP2=ABS(WW(K+1)/HW(K+1))*(1-ABS(WW(K+1)/HW(K+1)))
        GAMA=MAX(TEMP1,TEMP2)
	BETAW(K)=2.0D0/(2.0D0-GAMA)
	QWSTAR(K)=QW(K)-BETAW(K)*(FW(K+1)-FW(K))*DTDSG(K)
      ENDDO
C
      QWSTAR(1)=QW(1)-BETAW(2)*FW(2)*DTDSG(1)
      QWSTAR(NL)=QW(NL)+BETAW(NM)*FW(NL)*DTDSG(NL)
C
      DO K=1,NL
	IF(K.EQ.1) THEN
	  QWMIN(I,J,K)=MIN(QW(K),QW(K+1),QWMIN(I,J,K))
	  QWMAX(I,J,K)=MAX(QW(K),QW(K+1),QWMAX(I,J,K))
CCCC	ELSE IF(J.EQ.NL) THEN
	ELSE IF(K.EQ.NL) THEN
	  QWMIN(I,J,K)=MIN(QW(K),QW(K-1),QWMIN(I,J,K))
	  QWMAX(I,J,K)=MAX(QW(K),QW(K-1),QWMAX(I,J,K))
	ELSE
          QWMIN(I,J,K)=MIN(QW(K),QW(K-1),QW(K+1),QWMIN(I,J,K))
	  QWMAX(I,J,K)=MAX(QW(K),QW(K-1),QW(K+1),QWMAX(I,J,K))
	ENDIF
	  AW(K)=(QWSTAR(K)-QWMIN(I,J,K))*(QWSTAR(K)-QWMAX(I,J,K))
      ENDDO
C
      DO K=2,NL
	TEMP1=(ABS(AW(K))+AW(K))/(AW(K)+EPSM)
	TEMP2=(ABS(AW(K-1))+AW(K-1))/(AW(K-1)+EPSM)
	TEMP3=(ABS(AW(K))+AW(K))*(ABS(AW(K-1))+AW(K-1))
	TEMP4=ABS(AW(K))*ABS(AW(K-1))+EPSM
	CWSTAR=HALF*(TEMP1+TEMP2)-FOURTH*TEMP3/TEMP4
	CW=CWSTAR+(1-CWSTAR)*ABS(WW(K)/HW(K))
	WSTAR(K)=CW*WW(K)
      ENDDO
C
      DO K=2,NL
	FW(K)=HALF*WW(K)*(QW(K)+QW(K-1))
     $ -HALF*ABS(WSTAR(K))*(QW(K)-QW(K-1))
      ENDDO
C
      FW(1)=ZERO
      FW(NZ)=ZERO
C
      DO K=1,NL
	QW(K)=QW(K)-(FW(K+1)-FW(K))*DTDSG(K)
      ENDDO
C
      DO K=1,NL
	Q(I,J,K)=QW(K)
      ENDDO
C
      ENDDO

#if (defined SPMD)

!      IF(J.EQ.1) THEN
      if (myrank.eq.nprocs-1.and.j.eq.jpole) then
        DO K = 1 ,NL
          DO I = 1 ,NX
            Q(I,J,K)= Q(NX,J,K)
          ENDDO
        ENDDO
!      ELSE IF(J.EQ.NY) THEN
      else if (myrank.eq.0.and.j.eq.jpole) then
        DO K = 1 ,NL
          DO I = 1 ,NX
            Q(I,J,K)= Q(1,J,K)
          ENDDO
        ENDDO
      ENDIF
#else

      IF(J.EQ.1) THEN
        DO K = 1 ,NL
          DO I = 1 ,NX
            Q(I,J,K)= Q(NX,J,K)
          ENDDO
        ENDDO
      ELSE IF(J.EQ.NY) THEN
        DO K = 1 ,NL
          DO I = 1 ,NX
            Q(I,J,K)= Q(1,J,K)
          ENDDO
        ENDDO
      ENDIF
#endif
C
      ENDDO


!- check ---------------------------------------------------------
!#if (defined SPMD)
!
!      write(filename,14) 'tspasw-p-',myrank,'.out'
!14    format(a9,i1,a4)
!
!      open (10,file=trim(filename))
!
!#else
!      open (10,file='tspasw-s.out')
!#endif
!
!      write(10,*) 'qqqqqqqq--------------tspasw'
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
!      stop'tspasw'
!--------------------------------------------------------------

      RETURN
      END

