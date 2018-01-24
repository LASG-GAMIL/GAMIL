#include <misc.h>
#include <params.h>

      SUBROUTINE TSPAS(Q, U, V, DT, DX, DY, SINU, SINV, WTGU, WTGV)

      use pmgrid, only: beglatexdyn, endlatexdyn, plat
      use mpi_gamil
      use comfm1, only: CX, CY, A, QH, QHSTAR, VSTAR, BETA, FX, FY

      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

      REAL*8 Q(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 U(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 V(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 DT, DX, DY
      REAL*8 SINU(beglatexdyn:endlatexdyn)
      REAL*8 SINV(beglatexdyn:endlatexdyn)
      REAL*8 WTGU(beglatexdyn:endlatexdyn)
      REAL*8 WTGV(beglatexdyn:endlatexdyn)

      REAL*8 QMIN, QMAX, GAMA, CXSTAR, CYSTAR, TEMP1, TEMP2, TEMP3, TEMP4
      REAL*16 XNP(NL), XNP_TMP(8,NL), XSP(NL), XSP_TMP(8,NL), TMP_SUM
      REAL*8 QH_pole_1(NL), QH_pole_1_south(NL), QH_pole_1_north(NL)
      REAL*8 ZERO, HALF, FOURTH, EPSM
      DATA ZERO, HALF, FOURTH, EPSM /0.0D0,0.5D0,0.25D0,1.0D-80/

      INTEGER I, J, K
      INTEGER comm_request1, comm_request2

      call gamil_arrays_comm(COMM_TO_RIGHT, 1, Q(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_TOP,   1, Q(:,beglatexdyn,1))

      call gamil_arrays_comm(COMM_TO_LEFT, 1, Q(:,beglatexdyn,1), request_id=comm_request1)
      call gamil_arrays_comm(COMM_TO_BOT,  1, Q(:,beglatexdyn,1), request_id=comm_request2)

      DO K = 1, NL
        DO J = jbeg1, jend1
          DO I = ibeg1, iend1
            CX(I,J,K) = DT/(RAD*DX)*U(I,J,K)
          END DO
        END DO
        DO J = jbeg0, jend1
          DO I = ibeg1, iend1
            CY(I,J,K) = DT/(RAD*DY)*V(I,J,K)*SINV(J) ! 注意：缺少权重系数WTGV，在计算FY时添加上
          END DO
        END DO
      END DO

      call gamil_arrays_comm(COMM_TO_RIGHT, 1, CX(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_TOP,   1, CY(:,beglatexdyn,1)) ! TODO: 需要这个通信吗？

!     计算二阶Lax-Wendroff方案。
      DO K = 1, NL
        DO J = jbeg1, jend1
          DO I = ibeg1, iend1
            FX(I,J,K) = HALF*(CX(I,J,K)   *(Q(I,J,K)+Q(I-1,J,K))
     $                       -CX(I,J,K)**2*(Q(I,J,K)-Q(I-1,J,K))/SINU(J))
          END DO
        END DO
        DO J = jbeg0, jend1
          DO I = ibeg1, iend1
            FY(I,J,K) = HALF*(CY(I,J,K)   *(Q(I,J+1,K)+Q(I,J,K))
     $                       -CY(I,J,K)**2*(Q(I,J+1,K)-Q(I,J,K))/SINV(J)*WTGV(J))
          END DO
        END DO
      END DO

      call wait_icomm_request(comm_request1)
      call wait_icomm_request(comm_request2)
      call gamil_arrays_comm(COMM_TO_LEFT, 1, FX(:,beglatexdyn,1), request_id=comm_request1)
      call gamil_arrays_comm(COMM_TO_BOT,  1, FY(:,beglatexdyn,1), request_id=comm_request2)

      DO K = 1, NL
        DO J = jbeg1, jend1
          DO I = ibeg1, iend1
            TEMP1 = ABS(CX(I  ,J,K))*(1-ABS(CX(I  ,J,K)))
            TEMP2 = ABS(CX(I+1,J,K))*(1-ABS(CX(I+1,J,K)))
            TEMP3 = ABS(CY(I,J-1,K)/SINV(J-1)*WTGV(J-1))*(1-ABS(CY(I,J-1,K)/SINV(J-1)*WTGV(J-1)))
            TEMP4 = ABS(CY(I,J  ,K)/SINV(J  )*WTGV(J  ))*(1-ABS(CY(I,j  ,K)/SINV(J  )*WTGV(J  )))
            GAMA = MAX(TEMP1, TEMP2, TEMP3, TEMP4)
            BETA(I,J,K) = 2.0D0/(2.0D0-2.0D0*GAMA)
          END DO
        END DO
      END DO

      call wait_icomm_request(comm_request1)
      call wait_icomm_request(comm_request2)

      DO K = 1, NL
        DO J = jbeg1, jend1
          DO I = ibeg1, iend1
            QHSTAR(I,J,K) = Q(I,J,K) - BETA(I,J,K)*(
     $         FX(I+1,J,K) - FX(I,J  ,K) +
     $        (FY(I  ,J,K) - FY(I,J-1,K)) * WTGU(J)
     $      ) / SINU(J)
          END DO
        END DO
      END DO
!     处理北极点。
      IF (jbeg0 .eq. 1) THEN
        DO K = 1, NL
          TMP_SUM = ZERO
          DO I = ibeg1, iend1
            TMP_SUM = TMP_SUM + BETA(I,jbeg1,K)*FY(I,jbeg0,K)
          END DO
          XNP_TMP(1,K) = TMP_SUM
        END DO
        DO K = 1, NL
          XNP(K) = XNP_TMP(1,K)
        END DO
        call gamil_sum_pole_data_phys(1, XNP, NL)

        DO K = 1, NL
          XNP_TMP(1,K) = Q(ibeg1,jbeg0,K) - XNP(K) * 4.0d0/PLON * WTGV(jbeg0)/SINV(jbeg0)
          DO I = ibeg1, iend1
            QHSTAR(I,jbeg0,K) = XNP_TMP(1,K)
          END DO
        END DO
      END IF
!     处理南极点。
      IF (jend0 .eq. PLAT) then
        DO K = 1, NL
          TMP_SUM   = ZERO
          DO I = ibeg1, iend1
            TMP_SUM = TMP_SUM + BETA(I,jend1,K)*FY(I,jend1,K)
          END DO
          XSP_TMP(1,K) = TMP_SUM
        END DO
        DO K = 1, NL
          XSP(K) = XSP_TMP(1,K)
        END DO
        call gamil_sum_pole_data_phys(PLAT, XSP, NL)
        DO K = 1, NL
          XSP_TMP(1,K) = Q(ibeg1,jend0,K) + XSP(K) * 4.0d0/PLON * WTGV(jend1)/SINV(jend1)
          DO I = ibeg1, iend1
            QHSTAR(I,PLAT,K) = XSP_TMP(1,K)
          END DO
        END DO
      END IF

      DO K = 1, NL
        DO J = jbeg0, jend0
          DO I = ibeg1, iend1
            QMIN =  1.0E15
            QMAX = -1.0E15
            IF (J .EQ. 1) THEN
              QMIN = MIN(Q(I,J,K), Q(I,J+1,K), QMIN)
              QMAX = MAX(Q(I,J,K), Q(I,J+1,K), QMAX)
            ELSE IF (J .EQ. PLAT) THEN
              QMIN = MIN(Q(I,J,K), Q(I,J-1,K), QMIN)
              QMAX = MAX(Q(I,J,K), Q(I,J-1,K), QMAX)
            ELSE
              QMIN = MIN(Q(I,J  ,K), Q(I+1,J  ,K), Q(I-1,J,K),
     $                   Q(I,J+1,K), Q(I  ,J-1,K), QMIN)
              QMAX = MAX(Q(I,J  ,K), Q(I+1,J  ,K), Q(I-1,J,K),
     $                   Q(I,J+1,K), Q(I  ,J-1,K), QMAX)
            END IF
            A(I,J,K) = (QHSTAR(I,J,K)-QMAX)*(QHSTAR(I,J,K)-QMIN)
          END DO
        END DO
      END DO

      call gamil_arrays_comm(COMM_TO_RIGHT, 1, A(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_TOP,   1, A(:,beglatexdyn,1))

!     计算迎风格式。
      DO K = 1, NL
        DO J = jbeg1, jend1
          DO I = ibeg1, iend1
            TEMP1 = (ABS(A(I-1,J,K))+A(I-1,J,K))/(ABS(A(I-1,J,K))+EPSM)
            TEMP2 = (ABS(A(I  ,J,K))+A(I  ,J,K))/(ABS(A(I  ,J,K))+EPSM)
            TEMP3 = (ABS(A(I-1,J,K))+A(I-1,J,K))*(ABS(A(I,J,K))+A(I,J,K))
            TEMP4 = ABS(A(I-1,J,K))*ABS(A(I,J,K))+EPSM
            CXSTAR = HALF*(TEMP1+TEMP2)-FOURTH*TEMP3/TEMP4
            TEMP1 = CXSTAR+(1-CXSTAR)*ABS(CX(I,J,K))
            FX(I,J,K) = HALF*(          CX(I,J,K) *(Q(I,J,K)+Q(I-1,J,K))
     $                       -ABS(TEMP1*CX(I,J,K))*(Q(I,J,K)-Q(I-1,J,K)))
          END DO
        END DO

        DO J = jbeg0, jend1
          DO I = ibeg1, iend1
            TEMP1 = (ABS(A(I,J  ,K))+A(I,J  ,K))/(ABS(A(I,J  ,K))+EPSM)
            TEMP2 = (ABS(A(I,J+1,K))+A(I,J+1,K))/(ABS(A(I,J+1,K))+EPSM)
            TEMP3 = (ABS(A(I,J,K))+A(I,J,K))*(ABS(A(I,J+1,K))+A(I,J+1,K))
            TEMP4 = ABS(A(I,J,K))*ABS(A(I,J+1,K))+EPSM
            CYSTAR = HALF*(TEMP1+TEMP2)-FOURTH*TEMP3/TEMP4
            TEMP1 = CYSTAR+(1-CYSTAR)*ABS(CY(I,J,K)*WTGU(J))
            FY(I,J,K) = HALF*(          CY(I,J,K) *(Q(I,J+1,K)+Q(I,J,K))
     $                       -ABS(TEMP1*CY(I,J,K))*(Q(I,J+1,K)-Q(I,J,K)))
          END DO
        END DO
      END DO

      call gamil_arrays_comm(COMM_TO_LEFT, 1, FX(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_BOT,  1, FY(:,beglatexdyn,1))

      DO K = 1, NL
        IF (jbeg0 .eq. 1) THEN
          TMP_SUM = ZERO
          DO I = ibeg1, iend1
            TMP_SUM = TMP_SUM + FY(I,jbeg0,K)
          END DO
          XNP_TMP(1,K) = TMP_SUM
        END IF
        IF (jend0 .eq. PLAT) THEN
          TMP_SUM = ZERO
          DO I = ibeg1, iend1
            TMP_SUM = TMP_SUM + FY(I,jend1,K)
          END DO
          XSP_TMP(1,K) = TMP_SUM
        END IF
      END DO

      DO K = 1, NL
        XNP(K) = XNP_TMP(1,K)
        XSP(K) = XSP_TMP(1,K)
      END DO

      call gamil_sum_pole_data_phys(1,    XNP, NL)
      call gamil_sum_pole_data_phys(PLAT, XSP, NL)

      DO K = 1, NL
        DO J = jbeg1, jend1
          DO I = ibeg1, iend1
            QH(I,J,K) = Q(I,J,K) - (
     $        FX(I+1,J,K) - FX(I,J  ,K) +
     $       (FY(I  ,J,K) - FY(I,J-1,K)) * WTGU(J)
     $      ) / SINU(J)
          END DO
        END DO
      END DO

      call gamil_arrays_comm(COMM_TO_LEFT,  1, QH(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT, 1, QH(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_TOP,   1, QH(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_BOT,   1, QH(:,beglatexdyn,1))

      DO K = 1, NL
        IF (jbeg0 .eq. 1) THEN
          XNP_TMP(1,K) = Q(ibeg1,jbeg0,K) - XNP(K) * 4.0d0/PLON * WTGV(jbeg0)/SINV(jbeg0)
          DO I = beglonex, endlonex
            QH(I,1,K) = XNP_TMP(1,K)
            Q(I,1,K) = QH(I,1,K)
          END DO
        END IF

        IF (jend0 .eq. PLAT) THEN
          XSP_TMP(1,K) = Q(ibeg1,jend0,K) + XSP(K) * 4.0d0/PLON * WTGV(jend1)/SINV(jend1)
          DO I = beglonex, endlonex
            QH(I,PLAT,K) = XSP_TMP(1,K)
            Q(I,PLAT,K) = QH(I,PLAT,K)
          END DO
        END IF
!       计算空间交叉项。
        DO J = jbeg1, jend1
          DO I = ibeg1, iend1
!            TEMP1 = 0.5*(CX(I,J  ,K)+CX(I+1,J,K))/SINU(J)
!            TEMP2 = 0.5*(CY(I,J-1,K)+CY(I  ,J,K))*WTGU(J)
!            TEMP3 = 0.25D0*(QH(I+1,J+1,K) +
!     $                      QH(I-1,J-1,K) -
!     $                      QH(I+1,J-1,K) -
!     $                      QH(I-1,J+1,K))
!            Q(I,J,K) = QH(I,J,K)+TEMP1*TEMP2*TEMP3
            Q(I,J,K) = QH(I,J,K)
          END DO
        END DO
      END DO

      call gamil_arrays_comm(COMM_TO_LEFT, 1, Q(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1, Q(:,beglatexdyn,1))

      RETURN
      END
