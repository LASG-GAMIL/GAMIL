#include <misc.h>
#include <params.h>

      SUBROUTINE DYNAMICS(DTDY, ITIME, IHDIFUS, 
     _                    US, VS, WS, W0, PS, TS, H, HS, DP, SU, SV, ST, FAC, FBC, 
     _                    PMTOP, SIGL, DSIG, 
     _                    TBB, HBB, CBB, DCBB, PSB, TSB, 
     _                    DY, WTGU, WTGV, 
     _                    DX, SINU, SINV, OUX, OUY, OVX, OVY, FF, CUR, 
     _                    MM1, MP1, MM2, MP2, MM3, MP3, MDJ )

!---------------------------------------------------------------------------------
!
!	    This subroutine is the main part of the model dynamical
!	framework developed by Dr. Bin Wang from LASG/IAP in September
!	of 2001 recoded in May of 2002.
!         In this framework, some new numerical methods are used,
!     including:
!         1) the explicit difference scheme with exact linear and square
!            conservations, developed by Bin Wang, Zhongzhen Ji and
!            Qingcun Zeng, used to solve the atmopsheric equations. By
!            using this scheme, the model can conserve the total mass
!            exactly all the time and keep the conservation of the total
!            available energy when ignoring the outer forcing and friction,
!            and DLT1=DLT2=0.0
!         2) the weighted even-area coordinate along the latitude,
!            developed by Bin Wang, used to reduce the instability
!            of the model at the poles. The coordinate is produced by
!            the subroutine LATMESH
!         3) the flexible leaping-grid method, developed by Bin Wang to
!            further reduce the instability at the ploes. By using the
!            weighted even-area coordinate and the flexible leap-grid method,
!            no filter or smoother is needed in the model.
!         4) the reduction of the standard atmosphere, developed by
!            Qingcun Zeng, to improve the prediction accuracy.
!         This dynamica framework is normalized by Dr. Bin Wang and Rucong Yu,
!     Parallelized by Xin Zhang, translated to F90 by Pu Ye.
!
!---------------------------------------------------------------------------------
      USE COMFM1,    ONLY: P, C0, PLY, CB, DCB, TB, CB0, NIGW, UU, VV, TT, DPS0
      USE COMFM1,    ONLY: DPS1, HPS, PSK, DPS, UK, VK, TTK, DU, DV, DTT, DU1, DV1
      USE COMFM1,    ONLY: DTT1, UUK, HHK, DUS
      USE MPI_GAMIL, ONLY: BEGLONEX,ENDLONEX, IBEG0, JBEG0, IBEG1, IEND1, JBEG0, JEND0, JBEG1, JEND1
      USE MPI_GAMIL, ONLY: GAMIL_ARRAYS_COMM, GAMIL_MIN_LAT_ROW_DATA
      USE MPI_GAMIL, ONLY: COMM_TO_LEFT, COMM_TO_RIGHT
      USE PMGRID,    ONLY: BEGLATEXDYN, ENDLATEXDYN

      IMPLICIT NONE
                          
#include <PARADYN>

#if (defined SPMD)
#include <commpi.h>
#endif

#include <comfm2.h>

      REAL*8  PMTOP         !  MODEL TOP PRESSURE
      REAL*8  SIGL(NL)      !  HALF-MOVE VERTICAL LAYERS (SIGMA(K+1/2)) (UNIT: 1)
      REAL*8  DSIG(NL)      !  VERTICAL STEPSIZES (UNIT: 1)

      REAL*8  HS(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN) !  GEOPOTENTIAL OF THE GROUND ELEVATION  (UNIT: M^2/S^2)
      REAL*8  DTDY          !  TIME STEPSIZE (UNIT: S)

      INTEGER*4 ITIME    !  INTEGRATED TIME
      INTEGER   IHDIFUS  !  CONTROL PARAMETER RELATED TO THE HORIZONTAL

      REAL*8 SU(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN,NL)  !  THE FORCING ON U
      REAL*8 SV(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN,NL)  !  THE FORCING ON V
      REAL*8 ST(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN,NL)  !  THE FORCING ON T

      REAL*8 US(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN,NL)  !  ZONAL WIND        (UNIT: M/S)
      REAL*8 VS(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN,NL)  !  MERIDIONAL WIND   (UNIT: M/S)
      REAL*8 TS(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN,NL)  !  TEMPERATURE       (UNIT: K  )
      REAL*8 PS(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN   )  !  SURFACE PRESSURE  (UNIT: HPA)

      REAL*8 WS(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN,NZ)  !  DIAGNOSED VERTICAL VELOCITY AT THE PRESENT STEP (UNIT: 1/S)
      REAL*8 W0(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN,NL)  !  DIAGNOSED VERTICAL VELOCITY AT THE PREVIOUS STEP (UNIT: 1/S)
      REAL*8 H (BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN,NZ)  !  DIAGNOSED GEOPOTENTIAL HEIGHT (UNIT: M^2/S^2)
      REAL*8 DP(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN   )  !  THE TENDENCY OF THE SURFACE PRESSURE (UNIT: HPA/S)

      REAL*8 TBB (NA)       !  TEMPERATURE         (UNIT: K)
      REAL*8 HBB (NA)       !  GEOPOTENTIAL HEIGHT (UNIT: M^2/S^2)
      REAL*8 CBB (NA)       !
      REAL*8 DCBB(NA)       !
      REAL*8 PSB(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)    !  SURFACE PRESSURE OF THE STANDARD ATMOPSHERE  (UNIT: HPA)
      REAL*8 TSB(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)    !  SURFACE TEMPERATURE OF THE STANDARD ATMOPSHERE (UNIT: K  )

      INTEGER MM1(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)   !  MM1=I-MDJ      , IF MM1<2, MM1=NX-2+MM1
      INTEGER MP1(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)   !  MP1=I+MDJ      , IF MP1>NX-1, MP1=MP1-NX+2
      INTEGER MM2(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)   !  MM2=I-(MDJ-1)/2, IF MM2<2, MM2=NX-2+MM2
      INTEGER MP2(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)   !  MP2=I+(MDJ+1)/2, IF MP2>NX-1, MP2=MP2-NX+2
      INTEGER MM3(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)   !  MM3=I-(MDJ+1)/2, IF MM3<2, MM3=NX-2+MM3
      INTEGER MP3(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)   !  MP3=I+(MDJ-1)/2, IF MP3>NX-1, MP3=MP3-NX+2
      INTEGER MDJ(BEGLATEXDYN:ENDLATEXDYN)                     !  LEAPING SPAN OF THE DIFFERENCE

      REAL*8 WTGU(BEGLATEXDYN:ENDLATEXDYN)                     !  AREA WEIGHTING AT THE NORMAL GRID    (UNIT: 1)
      REAL*8 WTGV(BEGLATEXDYN:ENDLATEXDYN)                     !  AREA WEIGHTING AT THE HALF-MOVE GRID (UNIT: 1)
      REAL*8 DY                                                !  MERIDIONAL STEPSIZE IN COMPUTING MESH
      REAL*8 DX
      REAL*8 SINU(BEGLATEXDYN:ENDLATEXDYN)                     !  SIN(THETA) AT INTEGER GRID J
      REAL*8 SINV(BEGLATEXDYN:ENDLATEXDYN)                     !  SIN(THETA) AT HALF GRID J+1/2
      REAL*8 OUX (BEGLATEXDYN:ENDLATEXDYN)                     !  OUX=1/(RAD*SINU*DX*MDJ)
      REAL*8 OVX (BEGLATEXDYN:ENDLATEXDYN)                     !  OUX=1/(RAD*SINV*DX*MDJ)
      REAL*8 OUY (BEGLATEXDYN:ENDLATEXDYN)                     !  OUY=1/(RAD*SINU*DY*WTGU)
      REAL*8 OVY (BEGLATEXDYN:ENDLATEXDYN)                     !  OUY=1/(RAD*SINV*DY*WTGU)
      REAL*8 FF  (BEGLATEXDYN:ENDLATEXDYN)
      REAL*8 CUR (BEGLATEXDYN:ENDLATEXDYN)

      REAL*8  WK0, WK5, WK6, BYY1, BYY2, BYY3, DT2, TE, TM
      REAL*8  DPS2(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN)
      REAL*8  INNER
      REAL*8  FAC(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN,NZ)
      REAL*8  FBC(BEGLONEX:ENDLONEX,BEGLATEXDYN:ENDLATEXDYN,NL)
      REAL*8  TMS(BEGLATEXDYN:ENDLATEXDYN)
      REAL*8  TMS_BUF(16,BEGLATEXDYN:ENDLATEXDYN)
      INTEGER I, J, K, KPP, KWB

      EXTERNAL INNER

!--------------------------------------------------------------------------
! AS FOLLOWS, THE IMPROVED LEAP-FROG AND THE REGENERATION OF
! VARIABLES WILL BE FINISHED IN THE SAME CYCLES.
!--------------------------------------------------------------------------

      DT2 = 0.5D0*DTDY

      IF (ITIME .EQ. 0) THEN
!
!	FOR THE NORMAL RUN, DLT1 AND DTL2 MUST BE 1.0D0
!
        DLT1 = 1.0d0
        DLT2 = 1.0d0
!
!       DLT1=0.OD0 AND DLT2=0.0D0 ONLY FOR TESTING
!       THE CONSERVATION OF AVAILABLE ENERGY
!
        CALL GAMIL_ARRAYS_COMM(COMM_TO_LEFT,  1, PS(:,BEGLATEXDYN)) 
        CALL GAMIL_ARRAYS_COMM(COMM_TO_RIGHT, 1, PS(:,BEGLATEXDYN)) 

!$OMP PARALLEL DO PRIVATE (I, J, K, WK5, KPP, WK6)
        DO J = JBEG0, JEND0
          TMS_BUF(1,J)=0.0
          DO I = BEGLONEX, ENDLONEX
            P (I,J) = SQRT(PS(I,J))
            C0(I,J) = SQRT(RD*TSB(I,J)/PSB(I,J))
          END DO
          DO K = 1, NL
            DO I = BEGLONEX, ENDLONEX
              PLY(I,J,K) = PS(I,J)*SIGL(K) + PMTOP
              WK5 = PLY(I,J,K)/DPALIB
              KPP = INT(WK5)
              WK6 = WK5 - DFLOAT(KPP)
              CB (I,J,K) = (1.D0-WK6)*CBB (KPP) + WK6*CBB (KPP+1)
              DCB(I,J,K) = (1.D0-WK6)*DCBB(KPP) + WK6*DCBB(KPP+1)
              TB (I,J,K) = (1.D0-WK6)*TBB (KPP) + WK6*TBB (KPP+1)
              CB0(I,J,K) = CB(I,J,K)*P(I,J)/PLY(I,J,K)
              WK5 = CB0(I,J,K)*P(I,J)*DSIG(K)
              IF (TMS_BUF(1,J) .LT. WK5) TMS_BUF(1,J) = WK5
            END DO
          END DO
        END DO

        DO J = JBEG0, JEND0
          TMS(J) = TMS_BUF(1,J)
        END DO
        CALL GAMIL_MIN_LAT_ROW_DATA(TMS(JBEG0:), JEND0-JBEG0+1)

        DO J = jbeg0, jend0
          KPP     = DTDY*TMS(J)*OUX(J)*2.0 + 0.0001
          NIGW(J) = KPP + 1
        END DO

        CALL IAPTRSF(US, VS, WS, PS, TS, UU, VV, P, TT, TB, CB, 1)

        CALL DIFPS0(UU, P, DPS0, DSIG, OUX, MP1, MP2, MM1, MM2)
        CALL DIFPS1(UU, VV, P, PS, WS, DPS1, DPS0,
     _              DSIG, DY, OUX, OUY, SINV, MP1, MP2, MM1, MM2, WTGV)

        CALL DIAG(UU, VV, P, PS, PLY, TT, US, VS, TS, H,
     _            HPS, PMTOP, PSB, TSB, TB, CB, FAC, DSIG)

!$OMP PARALLEL DO PRIVATE(I, J, K)
        DO K = 1, NL
          DO J = JBEG0, JEND0
            DO I = BEGLONEX, ENDLONEX
              W0(I,J,K) = WS(I,J,K)
            END DO
          END DO
        END DO
      END IF

!$OMP PARALLEL DO PRIVATE(I, J, K)
      DO J = JBEG0, JEND0
        DO I = BEGLONEX, ENDLONEX
          PSK(I,J) = PS(I,J)
          DPS(I,J) = DPS0(I,J) + DPS1(I,J)
        END DO
        DO K = 1, NL
          DO I = BEGLONEX, ENDLONEX
            UK (I,J,K) = UU(I,J,K)
            VK (I,J,K) = VV(I,J,K)
            TTK(I,J,K) = TT(I,J,K)
          END DO
        END DO
        DO I = BEGLONEX,ENDLONEX
          HHK(I,J) = C0(I,J)*(PS(I,J)+PMTOP-PSB(I,J))
          UUK(I,J) = 0.0D0
          DO K = 1,NL
            UUK(I,J) = UUK(I,J)+UU(I,J,K)*DSIG(K)
          END DO
        END DO
      END DO

      DO KWB = 1, 2
!-------------------------------------------------------------------------------
! CALCULATING DU/DT, DV/DT AND DTT/DT.
!  FORMING THE FACTORS OUT OF ALL THE CYCLES AND CLEANING SOME SPECIAL ARRAIES.
!-------------------------------------------------------------------------------
        
        CALL DIFUVT(UU, US, VV, VS, WS, P, PS, PLY, DPS, TT, H, HPS, CB, DCB,
     _              SIGL, DSIG, DY, OUX, OVX, OUY, OVY, SINV, FF,
     _              CUR, DLT1, DLT2, MP1, MP2, MP3, MM1, MM2, MM3, WTGV,
     _              DU, DV, DTT, SU, SV, ST, FBC)

        CALL SEMIU(UUK, HHK, P, DU, DPS1, DUS, DPS0, DTDY, OUX, DSIG, C0)

!$OMP PARALLEL 
!$OMP DO PRIVATE (I, J, K, WK5, KPP, WK6)
        DO J = JBEG0, JEND0
          DO I = BEGLONEX, ENDLONEX
            DPS2(I,J) = DPS1(I,J) + DPS0(I,J)
            PS  (I,J) = PSK (I,J) + DT2*DPS2(I,J)
          END DO

          DO K = 1, NL
            DO I = BEGLONEX, ENDLONEX
              DU(I,J,K)  = DU (I,J,K) +     DUS(I,J  )
              UU(I,J,K)  = UK (I,J,K) + DT2*DU (I,J,K)
              TT(I,J,K)  = TTK(I,J,K) + DT2*DTT(I,J,K)
              VV(I,J,K)  = VK (I,J,K) + DT2*DV (I,J,K)
              PLY(I,J,K) = PS(I,J)*SIGL(K) + PMTOP
              WK5 = PLY(I,J,K)/DPALIB
              KPP = INT(WK5)
              WK6 = WK5-DFLOAT(KPP)
              CB (I,J,K) = (1.D0-WK6)*CBB (KPP)+WK6*CBB (KPP+1)
              DCB(I,J,K) = (1.D0-WK6)*DCBB(KPP)+WK6*DCBB(KPP+1)
              CB0(I,J,K) = CB(I,J,K)*P(I,J)/PLY(I,J,K)
            END DO
          END DO
        END DO
!$OMP END PARALLEL

        IF (KWB .LT. 2) THEN
!$OMP PARALLEL DO PRIVATE (I,J)
          DO J = JBEG0, JEND0
            DO I = BEGLONEX, ENDLONEX
              P(I,J) = SQRT(PS(I,J))
            END DO
          END DO
        END IF

        CALL DIFPS1(UU, VV, P, PS, WS, DPS1, DPS0,
     _              DSIG, DY, OUX, OUY, SINV, MP1, MP2, MM1, MM2, WTGV)

        CALL DIAG(UU, VV, P, PS, PLY, TT, US, VS, TS, H,
     _            HPS, PMTOP, PSB, TSB, TB, CB, FAC, DSIG)

!$OMP PARALLEL DO PRIVATE (I, J)
        DO J = JBEG0, JEND0
          DO I = BEGLONEX, ENDLONEX
            DPS (I,J) = DPS0(I,J)+DPS1(I,J)
          END DO
        END DO
      END DO

!     CALCULATING DU/DT, DV/DT AND DTT/DT.

!     FORMING THE FACTORS OUT OF ALL THE CYCLES AND CLEANING SOME
!     SPECIAL ARRAIES.

      CALL DIFUVT(UU, US, VV, VS, WS, P, PS, PLY, DPS, TT, H, HPS, CB,
     _            DCB, SIGL, DSIG, DY, OUX, OVX, OUY, OVY, SINV, FF,
     _            CUR, DLT1, DLT2, MP1, MP2, MP3, MM1, MM2, MM3, WTGV,
     _            DU1, DV1, DTT1, SU, SV, ST, FBC)

!$OMP PARALLEL 
!$OMP DO PRIVATE (I, J, K)
      DO K = 1, NL
        DO J = JBEG0, JEND0
          DO I = BEGLONEX,ENDLONEX
            DU1(I,J,K) = DU1(I,J,K)+DUS(I,J)
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT

!$OMP DO PRIVATE (I, J)
      DO J = JBEG0, JEND0
        DO I = BEGLONEX, ENDLONEX
          DPS1(I,J) = DPS (I,J)
          DPS (I,J) = DPS2(I,J)
        END DO
      END DO
!$OMP END PARALLEL

!     TO DEDUCT THE INNER GRAVITY WAVES FROM THE TENDENCES OF THE ZONAL WIND
!     AND THE TEMPERATURE DURING THE LONG-TIME-STEPSIZE INTEGRATIONS

      CALL MINUS_INGW(UU, P, TT, CB0, DSIG, OUX, DU1, DTT1, NIGW)

      BYY1 = INNER(DU1, DV1, DTT1, DPS1, 0, DU , DV , DTT , DPS , 0,
     _             DSIG, TSB, PSB, PMTOP, SINU, SINV, WTGU, WTGV)

!$OMP PARALLEL DO PRIVATE (I,J)
      DO J = JBEG0, JEND0
        DO I = BEGLONEX, ENDLONEX
          PS(I,J) = PSK(I,J)+DT2*DPS1(I,J)
          P (I,J) = SQRT(PS(I,J))
        END DO
      END DO
!
!       TO ADD THE INNER GRAVITY WAVES TO THE TENDENCES OF THE ZONAL WIND
!       AND THE TEMPERATURE BY THE METHOD OF SHORT-TIME-STEPSIZE INTEGRATIONS
!
      CALL PLUS_INGW(UK, P, TTK, CB0, DSIG, OUX, DTDY, SINU, WTGU,
     _               DU1, DTT1, BYY1, NIGW)

      BYY3 = INNER(DU1, DV1, DTT1, DPS1, 0, DU1, DV1, DTT1, DPS1, 0,
     _             DSIG, TSB, PSB, PMTOP, SINU, SINV, WTGU, WTGV)

      DT2 = DTDY*BYY1/BYY3

!$OMP PARALLEL DO PRIVATE (I, J, K, WK5, KPP, WK6)
      DO J = JBEG0, JEND0
        DO I = BEGLONEX, ENDLONEX
          PS(I,J) = PSK(I,J)+DT2*DPS1(I,J)
          P (I,J) = SQRT(PS(I,J))
        END DO

        TMS_BUF(1,J) = 0.0
        DO K = 1, NL
          DO I = BEGLONEX, ENDLONEX
            PLY(I,J,K) = PS(I,J)*SIGL(K) + PMTOP
            WK5 = PLY(I,J,K)/DPALIB
            KPP = INT(WK5)
            WK6 = WK5-DFLOAT(KPP)
            TB(I,J,K) = (1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
            CB(I,J,K) = (1.D0-WK6)*CBB(KPP)+WK6*CBB(KPP+1)
            DCB(I,J,K) = (1.D0-WK6)*DCBB(KPP)+WK6*DCBB(KPP+1)
            CB0(I,J,K) = CB(I,J,K)*P(I,J)/PLY(I,J,K)
            WK5 = CB0(I,J,K)*P(I,J)*DSIG(K)
            IF (TMS_BUF(1,J).LT.WK5) TMS_BUF(1,J) = WK5
            TT(I,J,K) = TTK(I,J,K)+DT2*DTT1(I,J,K)
            UU(I,J,K) = UK (I,J,K)+DT2*DU1 (I,J,K)
            VV(I,J,K) = VK (I,J,K)+DT2*DV1 (I,J,K)
          END DO
        END DO
      END DO

      DO J = JBEG0, JEND0
        TMS(J) = TMS_BUF(1,J)
      END DO
      CALL GAMIL_MIN_LAT_ROW_DATA(TMS(JBEG0:), JEND0-JBEG0+1)

      DO J = JBEG0, JEND0
        KPP = DTDY*TMS(J)*OUX(J)*2.0+0.0001
        NIGW(J) = KPP+1
      END DO

      I = ITIME
      K = I/86400
      IF (DLT1*DLT2 .EQ. 0) THEN
        BYY2 = INNER(UK, VK, TTK, PSK, 1, UK, VK, TTK, PSK, 1,
     _               DSIG, TSB, PSB, PMTOP, SINU, SINV, WTGU, WTGV)
        BYY1 = INNER(UU, VV, TT, PS, 1, UU, VV, TT, PS, 1,
     _               DSIG, TSB, PSB, PMTOP, SINU, SINV, WTGU, WTGV)
#IF (DEFINED SPMD)
        IF (MYRANK.EQ.0) THEN
#ENDIF
          WRITE(6, *) FLOAT(ITIME)+DTDY, BYY1, DT2
#IF (DEFINED SPMD)
        END IF
#ENDIF
      ELSE
        CALL TEM(UU, VV, PS, P, TT, HS, CB, TB, DX, DY, DSIG,
     _           SINU, SINV, WTGU, WTGV, TE, TM)

#IF (DEFINED SPMD)
        IF (MYRANK .EQ. 0) THEN
#ENDIF
          WRITE(6, '(1X,3E25.18)') TE,TM,DT2
#IF (DEFINED SPMD)
        END IF
#ENDIF
      END IF

      CALL DIFPS0(UU, P, DPS0, DSIG, OUX, MP1, MP2, MM1, MM2)
      CALL DIFPS1(UU, VV, P, PS, WS, DPS1, DPS0,
     _            DSIG, DY, OUX, OUY, SINV, MP1, MP2, MM1, MM2, WTGV)

      CALL DIAG(UU, VV, P, PS, PLY, TT, US, VS, TS, H,
     _          HPS, PMTOP, PSB, TSB, TB, CB, FAC, DSIG)

!$OMP PARALLEL DO PRIVATE (I, J)
      DO J = JBEG0, JEND0
        DO I = BEGLONEX, ENDLONEX
          DP(I,J) = DPS0(I,J)+DPS1(I,J)
        END DO
      END DO

      RETURN
      END
