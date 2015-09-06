#include <misc.h>
#include <params.h>



      SUBROUTINE TSPAS(Q,U,V,P,DSNP,DSSP,GC,DTDLT,DTDLN)

      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil

      use comfm1, only: hu, hv, cu, cv, h, tuu, tvv, A, QH, QHSTAR, VSTAR, BETA, FX, FY
      
!
!     ==================
!     TWO-STEP SHAPE-PRESERVING ADVECTION SCHEME DEVELOPED
!     BY YU RUCONG
!     ==================
!
!     PERFORM 2-D ADVECTION IN GLOBAL SPHERICAL GEOMETRY
!     WITH THE UNIFORM LAT-LON C-GRID MESH
!
      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

!
      REAL*8 Q(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 U(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 V(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 P(beglonex:endlonex,beglatexdyn:endlatexdyn)
!
      REAL*8 QMIN,QMAX

      REAL*8 C(beglonex:endlonex,2,nl)
!
!!    REAL*8 GAMA,XNP,XSP,CXSTAR,CYSTAR,CX,CY,TEMP1,TEMP2,TEMP3,TEMP4
      REAL*8 GAMA,CXSTAR,CYSTAR,CX,CY,TEMP1,TEMP2,TEMP3,TEMP4
      REAL*16 XNP(nl),XNP_TMP(8,nl),XSP(nl),XSP_TMP(8,nl),TMP_SUM
!
      REAL*8  ZERO,HALF,FOURTH,EPSM
      REAL*8  QH_pole_1(NL), QH_pole_1_south(NL), QH_pole_1_north(NL)
      REAL*8  DSNP,DSSP,DTDLN(beglatexdyn:endlatexdyn)
      REAL*8 DTDLT(beglatexdyn:endlatexdyn),GC(beglatexdyn:endlatexdyn)
      DATA ZERO,HALF,FOURTH,EPSM/ 0.0D0,0.5D0,0.25D0,1.0D-80/
*     DATA ZERO,HALF,FOURTH,EPSM/ 0.0D0,0.5D0,0.25D0,1.0E-6/
      INTEGER I,J,K
      INTEGER comm_request1, comm_request2, comm_request3
!

C

!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = jbeg0 , jend0
        IF(J.GE.2.AND.J.LE.PLAT-1) THEN
          DO I = ibeg1,iend1
            H(I,J) = GC(J) * P(I,J)
          ENDDO
        ELSE
          DO I = ibeg0, iend0
            H(I,J) = ZERO
          ENDDO
        ENDIF
      ENDDO

      call gamil_arrays_comm(COMM_TO_RIGHT,1,h(:,beglatexdyn),Q(:,beglatexdyn,1),request_id=comm_request1)
      call gamil_arrays_comm(COMM_TO_TOP,1,h(:,beglatexdyn),Q(:,beglatexdyn,1),request_id=comm_request2)
      call wait_icomm_request(comm_request1)
      call wait_icomm_request(comm_request2)

!$OMP PARALLEL 
!$OMP DO PRIVATE (J,I)
      DO J = jbeg0 , jend0
        IF(J.GE.2.AND.J.LE.PLAT-1) THEN
          DO I = ibeg1,iend1
            HU(I,J)  = HALF*(H(I,J)+H(I-1,J))
          ENDDO
          DO I = beglonex ,endlonex
            CU(I,J) = DTDLN(J) * HU(I,J)
          ENDDO
        ELSE
          DO I = beglonex ,endlonex
            HU(I,J) = ZERO
            CU(I,J) = ZERO
          ENDDO
        ENDIF
      ENDDO
!$OMP END DO NOWAIT

!$OMP DO PRIVATE (I,J)      
      DO J = jbeg0,jend1
        DO I = ibeg1,iend1
          HV(I,J)  = HALF*(H(I,J)+H(I,J+1))
        ENDDO
        DO I = beglonex ,endlonex
          CV(I,J) = DTDLT(J) * HV(I,J)
        ENDDO
      ENDDO
!$OMP END PARALLEL


!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K=1,NL
      DO J=jbeg0,jend0
      DO I = beglonex ,endlonex
        TUU(I,J,k)=CU(I,J)*U(I,J,K)
        TVV(I,J,k)=CV(I,J)*V(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,hu(:,beglatexdyn),TUU(:,beglatexdyn,1),Q(:,beglatexdyn,1),request_id=comm_request1)
      call gamil_arrays_comm(COMM_TO_BOT,1,TVV(:,beglatexdyn,1),Q(:,beglatexdyn,1),hv(:,beglatexdyn),request_id=comm_request2)
C
C
!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K=1,NL
       DO J=jbeg1,jend1
	DO I=ibeg1,iend1
	  FX(I,J,k)=HALF*TUU(I,J,k)*(Q(I,J,k)+Q(I-1,J,k))
     $             -HALF*TUU(I,J,k)*TUU(I,J,k)*(Q(I,J,k)-Q(I-1,J,k))/HU(I,J)
	ENDDO
       ENDDO
       DO J=jbeg0,jend1
        DO I=ibeg1,iend1
	  FY(I,J,k)=HALF*TVV(I,J,k)*(Q(I,J+1,k)+Q(I,J,k))
     $             -HALF*TVV(I,J,k)*TVV(I,J,k)*(Q(I,J+1,k)-Q(I,J,k))/HV(I,J)
	ENDDO
       ENDDO
      ENDDO
C
      call wait_icomm_request(comm_request1)
      call wait_icomm_request(comm_request2)
      call gamil_arrays_comm(COMM_TO_LEFT,1,FX(:,beglatexdyn,1),request_id=comm_request1)
      call gamil_arrays_comm(COMM_TO_BOT,1,FY(:,beglatexdyn,1),request_id=comm_request2)

!$OMP PARALLEL DO PRIVATE (I,J,K,TEMP1,TEMP2,TEMP3,TEMP4,GAMA)
      DO K=1,NL
       DO J=jbeg1,jend1
	   DO I=ibeg1,iend1
          TEMP1=ABS(TUU(I,J,k)/HU(I,J))*(1-ABS(TUU(I,J,k)/HU(I,J)))
          TEMP2=ABS(TUU(I+1,J,k)/HU(I+1,J))*(1-ABS(TUU(I+1,J,k)/HU(I+1,J)))
          TEMP3=ABS(TVV(I,J-1,k)/HV(I,J-1))*(1-ABS(TVV(I,J-1,k)/HV(I,J-1)))
          TEMP4=ABS(TVV(I,J,k)/HV(I,J))*(1-ABS(TVV(I,J,k)/HV(I,J)))
          GAMA=MAX(TEMP1,TEMP2,TEMP3,TEMP4)
          BETA(I,J,k)=2.0D0/(2.0D0-2.0D0*GAMA)
        ENDDO
       ENDDO
      ENDDO

      call wait_icomm_request(comm_request1)
      call wait_icomm_request(comm_request2)

!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K=1,NL
       DO J=jbeg1,jend1
	   DO I=ibeg1,iend1
          QHSTAR(I,J,k)=Q(I,J,k)-BETA(I,J,k)*(FX(I+1,J,k)
     $                    -FX(I,J,k)+FY(I,J,k)-FY(I,J-1,k))/H(I,J)
        ENDDO
       ENDDO
      ENDDO


      IF (jbeg0 .eq. 1) THEN
!$OMP PARALLEL DO PRIVATE (I,J,K,TMP_SUM)
      DO K=1,NL
      TMP_SUM   = ZERO
      DO I = ibeg1,iend1
      TMP_SUM   = TMP_SUM + BETA(I,2,k)*FY(I,1,k)/H(I,2)
      ENDDO
      XNP_TMP(1,k)   = TMP_SUM
      ENDDO
      DO K=1,NL
        XNP(k) = XNP_TMP(1,k)
        QH_pole_1(k) = Q(ibeg1,jbeg0,k)
      ENDDO
      call gamil_sum_pole_data_phys(1,XNP,NL)
      call broadcast_lon_data(2,1,QH_pole_1,NL) 
!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K=1,NL
      XNP_TMP(1,k)   = QH_pole_1(k) - XNP(k)*DSNP
      DO I = ibeg1,iend1
      QHSTAR(I,1,k)   = XNP_TMP(1,k)
      ENDDO
      ENDDO
      ENDIF

      IF (jend0 .eq. PLAT) then
!$OMP PARALLEL DO PRIVATE (I,J,K,TMP_SUM)
      DO K=1,NL
      TMP_SUM   = ZERO
      DO I = ibeg1,iend1
      TMP_SUM   = TMP_SUM + BETA(I,PLAT-1,k)*FY(I,PLAT-1,k)/H(I,PLAT-1)
      ENDDO
      XSP_TMP(1,k)   = TMP_SUM
      ENDDO
      DO K=1,NL
        XSP(k) = XSP_TMP(1,k)
        QH_pole_1(k) = Q(ibeg1,jend0,k)
      ENDDO
      call gamil_sum_pole_data_phys(PLAT,XSP,NL)
      call broadcast_lon_data(2,PLAT,QH_pole_1,NL) 
!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K=1,NL
      XSP_TMP(1,k)   = QH_pole_1(k)+ XSP(k)*DSSP
      DO I = ibeg1,iend1
      QHSTAR(I,PLAT,k)   = XSP_TMP(1,k)
      ENDDO
      ENDDO
      ENDIF
C

!$OMP PARALLEL DO PRIVATE (I,J,K,QMIN,QMAX,TEMP1,TEMP2,TEMP3,TEMP4,CXSTAR,CX)
      DO K=1,NL
      DO J=jbeg0,jend0
        DO I=ibeg1,iend1
        QMIN=1.0E15
        QMAX=-1.0E15
	IF(J.EQ.1) THEN
	  QMIN=MIN(Q(I,J,k),Q(I,J+1,k),QMIN)
	  QMAX=MAX(Q(I,J,k),Q(I,J+1,k),QMAX)
	ELSE IF(J.EQ.PLAT) THEN
	  QMIN=MIN(Q(I,J,k),Q(I,J-1,k),QMIN)
	  QMAX=MAX(Q(I,J,k),Q(I,J-1,k),QMAX)
	ELSE
	  QMIN=MIN(Q(I,J,k),Q(I+1,J,k),Q(I-1,J,k),
     $                    Q(I,J+1,k),Q(I,J-1,k),QMIN)
          QMAX=MAX(Q(I,J,k),Q(I+1,J,k),Q(I-1,J,k),
     $                    Q(I,J+1,k),Q(I,J-1,k),QMAX)
	ENDIF
        A(I,J,k)=(QHSTAR(I,J,k)-QMAX)*(QHSTAR(I,J,k)-QMIN)
        ENDDO
      ENDDO
      ENDDO
C
      call gamil_arrays_comm(COMM_TO_RIGHT,1,A(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_TOP,1,A(:,beglatexdyn,1))


!$OMP PARALLEL DO PRIVATE (I,J,K,TEMP1,TEMP2,TEMP3,TEMP4,CXSTAR,CX,CYSTAR,CY)
      DO K=1,NL

      DO J=jbeg1,jend1
      DO I=ibeg1,iend1
	TEMP1=(ABS(A(I-1,J,k))+A(I-1,J,k))/(ABS(A(I-1,J,k))+EPSM)
	TEMP2=(ABS(A(I,J,k))+A(I,J,k))/(ABS(A(I,J,k))+EPSM)
	TEMP3=(ABS(A(I-1,J,k))+A(I-1,J,k))*(ABS(A(I,J,k))+A(I,J,k))
	TEMP4=ABS(A(I-1,J,k))*ABS(A(I,J,k))+EPSM
	CXSTAR=HALF*(TEMP1+TEMP2)-FOURTH*TEMP3/TEMP4
        CX=CXSTAR+(1-CXSTAR)*ABS(TUU(I,J,k)/HU(I,J))
        FX(I,J,k)=HALF*TUU(I,J,k)*(Q(I,J,k)+Q(I-1,J,k))
     $           -HALF*ABS(CX*TUU(I,J,k))*(Q(I,J,k)-Q(I-1,J,k))
	ENDDO
       ENDDO


      DO J=jbeg0,jend1
      DO I=ibeg1,iend1
	TEMP1=(ABS(A(I,J,k))+A(I,J,k))/(ABS(A(I,J,k))+EPSM)
	TEMP2=(ABS(A(I,J+1,k))+A(I,J+1,k))/(ABS(A(I,J+1,k))+EPSM)
	TEMP3=(ABS(A(I,J,k))+A(I,J,k))*(ABS(A(I,J+1,k))+A(I,J+1,k))
	TEMP4=ABS(A(I,J,k))*ABS(A(I,J+1,k))+EPSM
	CYSTAR=HALF*(TEMP1+TEMP2)-FOURTH*TEMP3/TEMP4
        CY=CYSTAR+(1-CYSTAR)*ABS(TVV(I,J,k)/HV(I,J))
        FY(I,J,k)=HALF*TVV(I,J,k)*(Q(I,J+1,k)+Q(I,J,k))
     $ -HALF*ABS(CY*TVV(I,J,k))*(Q(I,J+1,k)-Q(I,J,k))
      ENDDO
      ENDDO

      ENDDO

C
      call gamil_arrays_comm(COMM_TO_LEFT,1,FX(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_BOT,1,FY(:,beglatexdyn,1))

!$OMP PARALLEL DO PRIVATE (I,J,K,TMP_SUM)
      DO K=1,NL
      IF (jbeg0 .eq. 1) THEN
      TMP_SUM      = ZERO
      DO I = ibeg1,iend1
      TMP_SUM      = TMP_SUM + FY(I,1,k)/H(I,2)
      ENDDO
      XNP_TMP(1,k)      = TMP_SUM
      ENDIF
      IF (jend0 .eq. PLAT) THEN
      TMP_SUM      = ZERO
      DO I = ibeg1,iend1
      TMP_SUM      = TMP_SUM + FY(I,PLAT-1,k)/H(I,PLAT-1)
      ENDDO
      XSP_TMP(1,k)      = TMP_SUM
      ENDIF
      ENDDO

      DO K=1,NL
        XNP(k) = XNP_TMP(1,k)
        XSP(k) = XSP_TMP(1,k)
      ENDDO

      call gamil_sum_pole_data_phys(1,XNP,NL)
      call gamil_sum_pole_data_phys(PLAT,XSP,NL)

      IF (jbeg0 .eq. 1) THEN
      DO K=1,NL
         QH_pole_1_north(k) = Q(ibeg1,jbeg0,k)
      ENDDO
      call broadcast_lon_data(2,1,QH_pole_1_north,NL) 
      ENDIF
      IF (jend0 .eq. PLAT) then
      DO K=1,NL
         QH_pole_1_south(k) = Q(ibeg1,jend0,k)
      ENDDO
      call broadcast_lon_data(2,PLAT,QH_pole_1_south,NL) 
      ENDIF


!$OMP PARALLEL DO PRIVATE (I,J,K,TMP_SUM)
      DO K=1,NL
       DO J=jbeg1,jend1
        DO I=ibeg1,iend1
          QH(I,J,k)=Q(I,J,k)-(FX(I+1,J,k)-FX(I,J,k)
     $                        +FY(I,J,k)-FY(I,J-1,k))/H(I,J)
        ENDDO
       ENDDO
      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,QH(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,QH(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_TOP,1,QH(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_BOT,1,QH(:,beglatexdyn,1))

!$OMP PARALLEL DO PRIVATE (I,J,K,TEMP1,TEMP2,TEMP3)
      DO K=1,NL

      IF (jbeg0 .eq. 1) THEN
      XNP_TMP(1,k)      = QH_pole_1_north(k) - XNP(k)*DSNP
      DO I = beglonex ,endlonex
      QH(I,1,k)   = XNP_TMP(1,k)
      Q(I,1,K)=QH(I,1,k)
      ENDDO
      ENDIF

      IF (jend0 .eq. PLAT) THEN
      XSP_TMP(1,k)      = QH_pole_1_south(k) + XSP(k)*DSSP
      DO I = beglonex ,endlonex
      QH(I,PLAT,k)   = XSP_TMP(1,k)
      Q(I,PLAT,K)=QH(I,PLAT,k)
      ENDDO
      ENDIF


      DO J=jbeg1,jend1
      DO I=ibeg1,iend1
	TEMP1=0.5*(TUU(I,J,k)+TUU(I+1,J,k))
	TEMP2=0.5*(TVV(I,J-1,k)+TVV(I,J,k))
	TEMP3=0.25*(QH(I+1,J+1,k)+QH(I-1,J-1,k)-QH(I+1,J-1,k)-QH(I-1,J+1,k))
	Q(I,J,K)=QH(I,J,k)+TEMP1*TEMP2*TEMP3/(H(I,J)*H(I,J))
      ENDDO
      ENDDO
      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,Q(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,Q(:,beglatexdyn,1))

      RETURN
	END

