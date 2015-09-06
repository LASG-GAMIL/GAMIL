#include <misc.h>
#include <params.h>

!!(2003.11.12-13)
!!(2004.03.30 wh)
!!-------------------


      SUBROUTINE HDIFUS(U,V,T,QV
     _                 ,FRDT,FRDS,FRDU,FRDV,FRDP,TB,PLY,DXVPN,DXVPS)
     
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
      use comfm1, only: D, DT, DS, DA, DB, VR, HQK, TK, HVK, HUK
      use comfm1, only: ROT, RLNT, RDLN, RDLT, TW
      
C     *****************
C     *****************
C
C     COMPUTE HORIZONTAL DIFFUSION
C     INCLOUDING: T & U & V & QV & QC & QR & QI & QS & QG
C             NON-LINEAR DIFFUSION OF           Smagorinsky (1963)
C             FORMULATED     AS W. M. Washington & D .L .Wiliamson
C             1977 :  A description of the NCAR global circulation
C             models, Methods in Computational Physics  17  113-73
C             DESCRETED IN C-GRID SYSTEM            BY X._Z. Liang
C             1986 :  The Design of IAP GCM and  the Simulation of
C                     Climate  and   Its Interseasonal Variability
C                     Ph.D Thesis                            250pp
C
      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

C
      REAL*8  FRDT(beglatexdyn:endlatexdyn,3)
     _           ,FRDS(beglatexdyn:endlatexdyn,3)
     _           ,FRDU(beglatexdyn:endlatexdyn,3)
     _           ,FRDV(beglatexdyn:endlatexdyn,3)
     _           ,FRDP(beglatexdyn:endlatexdyn,3)
C
      REAL*8  U(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
     _           ,V(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
     _            ,T(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 QV(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)


     
      REAL*8  TB(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
     _           ,PLY(beglonex:endlonex,beglatexdyn:endlatexdyn,NZ)
     _           ,DXVPN(beglatexdyn:endlatexdyn)
     _           ,DXVPS(beglatexdyn:endlatexdyn)
C
C
      REAL*8  ZERO,HALF,ONE
      DATA       ZERO,HALF,ONE / 0.0E0 , 0.5E0 , 1.0E0  /
      REAL*8  FRDTN,FRDTS,FRDSI,FRDSJ,TI,RIJ,RI,VRI,TDI,FT1,FT2,FT3
     &       ,FS1,FS2,FS3,U0,V0,DTJ,DSI,DIJ,R0,R1,D0,D1,RT
     &       ,RN,DIJ1,DIJ2,FU1,FU2,FU3,FV1,FV2,FV3,FA1,FA2,FA3,VR0
     &       ,DS0,DT0,DA0,DB0,RLNT0,ROT0,RTA0,RSB0,VRU,VRV,RTAI,RSBI
     &       ,RTAJ,RSBJ,FB1,FB2,FB3,FB4,TD0,Q0,Q0C,Q0R,Q0I,Q0S,Q0G
     &       ,UIJK,VIJK,TIJK,QIJK,QIJKC,QIJKR,QIJKI,QIJKS,QIJKG
     &       ,VRN,VRS,TDN,TDS,QKN,QKNC,QKNR,QKNI,QKNS,QKNG
     &       ,QKSV,QKSC,QKSR,QKSI,QKSS,QKSG,DQNC,DQNR,DQNI,DQNS,DQNG
     &       ,DQSC,DQSR,DQSI,DQSS,DQSG,FAN,FAS,EAN,EAS,T0N
     &       ,T0S,Q0N,Q0NC,Q0NR,Q0NI,Q0NS,Q0NG
     &       ,Q0S0,Q0SR,Q0SI,Q0SS,Q0SG,FVI,FVJ,DXVPNJ,DXVPSJ
      REAL*16 TMP_SUM1,TMP_SUM2,DTS_16(NL),DTN_16(NL),DQS_16(NL),DQN_16(NL)
      REAL*8 DTS_8,DTN_8,DQS_8,DQN_8,TMP_DATA1,TMP_DATA2
      REAL*8 T_pole_2(NL),PLY_pole_2(NL),TB_pole_2(NL),QV_pole_2(NL)
      INTEGER I,J,K,I1,II

CWB
C
!

      DO J  = jbeg0,jend0
      IF (J.EQ.1 .OR. J .EQ. PLAT) THEN
         IF(ibeg1 .EQ. 2) THEN
         DO K  = 1 ,NL
            T_pole_2(K) = T(ibeg1,J,K)
            PLY_pole_2(K) = PLY(ibeg1,J,K)
            TB_pole_2(K) = TB(ibeg1,J,K)
            QV_pole_2(K) = QV(ibeg1,J,K)
         ENDDO
         ENDIF
         call broadcast_lon_data(2,J,T_pole_2,NL,PLY_pole_2,NL,TB_pole_2,NL,QV_pole_2,NL)
         DO K  = 1 ,NL
            T(ibeg1,J,K) = T_pole_2(K)
            PLY(ibeg1,J,K) = PLY_pole_2(K)
            TB(ibeg1,J,K) = TB_pole_2(K)
            QV(ibeg1,J,K) = QV_pole_2(K)
         ENDDO       
      ENDIF
      ENDDO


!$OMP PARALLEL DO PRIVATE (I, J, K,TI,RI,RIJ,VRI,TDI)
      DO K  = 1 ,NL
C
C     CALCULATE DENSITY & SET 2-D FNX-1LDS TO BE DIFFUSED
C
      DO J = jbeg0 ,jend0
      DO I = beglonex ,endlonex
        HUK(I,J,K)  = U(I,J,K)
        HVK(I,J,K)  = V(I,J,K)
        HQK(I,J,K)  = QV(I,J,K)
        TW(I,J,K)  = T(I,J,K)
      ENDDO
      ENDDO
      
      DO J  = jbeg0,jend0
      IF (J .ge. 2 .and. J .le. PLAT-1) THEN
       DO I  = ibeg1,iend1
        TI        = TW(I,J,K)
        RIJ       = PLY(I,J,K) / (TI*RD)
        ROT(I,J,K)  = RIJ
        VR (I,J,K)  = ONE / RIJ
        TK (I,J,K)  = TI  - TB(I,J,K)
       ENDDO
      ELSE
        TI        = TW (ibeg1,J,K)
        RI        = PLY(ibeg1,J,K) / (TI*RD)
        TDI       = TI  - TB(ibeg1,J,K)
        VRI       = ONE / RI
        DO I  = ibeg1,iend1
          ROT(I,J,K)  = RI
          VR (I,J,K)  = VRI
          TK (I,J,K)  = TDI
        ENDDO
      ENDIF
      ENDDO
      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,TK(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,ROT(:,beglatexdyn,1),VR(:,beglatexdyn,1),TK(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_TOP,1,HUK(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_BOT,1,HVK(:,beglatexdyn,1))
C
C     CALCULATE DEFORMATION FNX-1LDS
C
!$OMP PARALLEL DO PRIVATE (I,K,TMP_SUM1)
      DO K  = 1 ,NL
      IF (jbeg0 .eq. 1) THEN
      TMP_SUM1   = ZERO
      DO I  = ibeg1,iend1
        TMP_SUM1 = TMP_SUM1 + HVK(I,1,K )
      ENDDO
      DTN_16(K) = TMP_SUM1
      ENDIF

      IF (jend0 .eq. PLAT) THEN
      TMP_SUM1   = ZERO
      DO I  = ibeg1,iend1
        TMP_SUM1 = TMP_SUM1 + HVK(I,PLAT-1,K)
      ENDDO
      DTS_16(K) = TMP_SUM1
      ENDIF
      ENDDO

      call gamil_sum_pole_data_phys(1,DTN_16,NL)
      call gamil_sum_pole_data_phys(PLAT,DTS_16,NL)


!$OMP PARALLEL DO PRIVATE (I,J,K,U0,V0,FRDTN,FRDSI,FRDSJ,FRDTS,DTS_8,
!$   &                                 DTN_8,FT1,FT2,FT3,FS1,FS2,FS3)
      DO K  = 1 ,NL
      DO J  = jbeg1,jend1
        FT1       = FRDT(J,1)
        FT2       = FRDT(J,2)
        FT3       = FRDT(J,3)
        FS1       = FRDS(J,1)
        FS2       = FRDS(J,2)
        FS3       = FRDS(J,3)
        DO I  = ibeg1,iend1
          U0        = HUK(I,J,K)
          V0        = HVK(I,J,K)
          DT(I,J,K)   = FT1*(HUK(I+1,J,K)-U0) - (FT2*V0 - FT3*HVK(I,J-1,K ))
          DS(I,J,K)   = FS1*(V0-HVK(I-1,J,K)) + (FS2*HUK(I,J+1,K) - FS3*U0)
        ENDDO
      ENDDO

      IF (jbeg0 .eq. 1) THEN
      FRDTN     = FRDT(1 ,1)      
      FRDSI     = FRDS(1 ,1)
      FRDSJ     = FRDS(1 ,2)
      DTN_8     = DTN_16(K)
      DTN_8     = DTN_8 * FRDTN
      DO I  = ibeg1,iend1
        DT(I,1,K )  = DTN_8
        DS(I,1,K )  = FRDSI*(HVK(I,1,K) - HVK(I-1,1,K)) + FRDSJ*HUK(I,2,K)
      ENDDO
      ENDIF

      IF (jend0 .eq. PLAT) THEN
      FRDTS     = FRDT(PLAT,1)
      DTS_8  = DTS_16(K)
      DTS_8  = DTS_8 * FRDTS
      DO I  = ibeg1,iend1
        DT(I,PLAT,K)  = DTS_8
        DS(I,PLAT,K)  = ZERO
      ENDDO
      ENDIF
      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,DS(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,DT(:,beglatexdyn,1),DS(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_TOP,1,DT(:,beglatexdyn,1),ROT(:,beglatexdyn,1))

      DO K  = 1 ,NL
      DO J  = jbeg0, jend1
        DO I  = ibeg1,iend1
          DTJ       = DT(I,J+1,K) + DT(I,J,K)
          DSI       = DS(I,J,K ) + DS(I-1,J,K)
          DIJ       = DTJ*DTJ  + DSI*DSI
          D(I,J,K)    = HALF * SQRT( DIJ )
        ENDDO
      ENDDO

      IF (jend0 .eq. PLAT) THEN
      DO I = beglonex ,endlonex
        D(I,PLAT,K) = ZERO
      ENDDO
      ENDIF
      ENDDO

      call gamil_arrays_comm(COMM_TO_RIGHT,1,D(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_BOT,1,D(:,beglatexdyn,1))

!$OMP PARALLEL DO PRIVATE (I,J,K,I1,R0,R1,D0,D1,RT,RN,DIJ,DIJ1,DIJ2)
      DO K  = 1 ,NL
      DO J  = jbeg0,jend0
      IF (J .ge. 2 .and. J .le. PLAT-1) THEN
      DO I  = ibeg1,iend1
        I1        = I - 1
        R0        = ROT(I,J,K)
        R1        = ROT(I,J+1,K)
        D0        = D(I,J,K)
        D1        = D(I,J-1,K)
        RT        = R0 + R1
        RN        = R0 + ROT(I1,J,K)
        RLNT(I,J,K) = RN + R1 + ROT(I1,J+1,K)
        DIJ1      = HALF * (D0 + D(I1,J,K))
        DIJ       = HALF * (D0 + D1)
        DIJ2      = HALF * (DIJ1 + HALF*(D1 + D(I1,J-1,K)) )
        RDLN(I,J,K) = RN   *  DIJ2
        RDLT(I,J,K) = RT   *  D0
        DA  (I,J,K) = DIJ
        DB  (I,J,K) = DIJ1
      ENDDO
      ELSE
      DO I = beglonex ,endlonex
        RDLN(I,J,K) = ZERO
        DA(I,J,K)   = ZERO
        RDLT(I,J,K) = ZERO
        RLNT(I,J,K) = ZERO
        DB(I,J,K)   = ZERO
      ENDDO      
      IF (J .eq. 1) THEN
      R0        = ROT(ibeg1,J,K)
      DO I  = ibeg1,iend1
        D0        = D(I,J,K)
        DB  (I,J,K) = (D0 + D(I-1,J,K )) * HALF
        RDLT(I,J,K) = (R0 + ROT(I,J+1,K)) * D0
        RLNT(I,J,K) =  R0 + R0 + ROT(I-1,J+1,K) + ROT(I,J+1,K)
      ENDDO
      ENDIF
      ENDIF
      ENDDO

      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,RLNT(:,beglatexdyn,1),RDLN(:,beglatexdyn,1),DB(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,DA(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_BOT,1,RLNT(:,beglatexdyn,1),DS(:,beglatexdyn,1),DB(:,beglatexdyn,1),
     &                       RDLT(:,beglatexdyn,1),TK(:,beglatexdyn,1),HQK(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_TOP,1,VR(:,beglatexdyn,1),TK(:,beglatexdyn,1),DA(:,beglatexdyn,1),
     &                       HQK(:,beglatexdyn,1))
C
C     UPDATE T & U & V & QV,QC,QI,QR,QS,QG DUE TO THE HORIZONTAL DIFFUSION
C
!$OMP PARALLEL DO PRIVATE (I,J,TDN,QKN,EAN,TDS,QKSV,EAS,TMP_SUM1,TMP_SUM2,TMP_DATA1,TMP_DATA2)
      DO K  = 1 ,NL
      IF (jbeg0 .eq. 1) THEN           ! North pole
      TDN       = TK(ibeg1,jbeg0,K )
      QKN       = HQK(ibeg1,jbeg0,K )
      TMP_SUM1  = ZERO
      TMP_SUM2  = ZERO
      DO I  = ibeg1,iend1
        EAN       = RDLT(I,1,K )
        TMP_DATA1  = EAN*(TK(I,2,K) - TDN)
        TMP_DATA2 = EAN*(HQK(I,2,K) - QKN)
        TMP_SUM1   = TMP_SUM1 + TMP_DATA1
        TMP_SUM2  = TMP_SUM2 + TMP_DATA2
      ENDDO
      DTN_16(K) = TMP_SUM1
      DQN_16(K) = TMP_SUM2
      ENDIF

      IF (jend0 .eq. PLAT) THEN   ! South pole
      TDS       = TK(ibeg1,jend0,K)
      QKSV      = HQK(ibeg1,jend0,K)
      TMP_SUM1  = ZERO
      TMP_SUM2  = ZERO
      DO I  = ibeg1,iend1
        EAS       = RDLT(I,PLAT-1,K)
        TMP_DATA1  = EAS*(TDS - TK(I,PLAT-1,K))
        TMP_DATA2 = EAS*(QKSV- HQK(I,PLAT-1,K))
        TMP_SUM1   = TMP_SUM1 + TMP_DATA1
        TMP_SUM2  = TMP_SUM2 + TMP_DATA2
      ENDDO
      DTS_16(K)    = TMP_SUM1
      DQS_16(K)    = TMP_SUM2
      ENDIF
      ENDDO

      call gamil_sum_pole_data_phys(1,DTN_16,NL,DQN_16,NL)
      call gamil_sum_pole_data_phys(PLAT,DTS_16,NL,DQS_16,NL)

!$OMP PARALLEL DO PRIVATE (I,J,K,FU1,FU2,FU3,FV1,FV2,FV3,FA1,FA2,FA3,II,I1,
!$   &                                     VR0,DS0,DT0,DA0,DB0,RLNT0,ROT0,RTA0,RSB0,
!$   &                                     VRU,VRV,RTAI,RSBI,RTAJ,RSBJ,FB1,FB2,FB3,FB4,
!$   &                                     TD0,Q0,UIJK,VIJK,TIJK,QIJK,VRN,TDN,QKN,
!$   &                                     FAN,EAN,T0N,Q0N,FVI,FVJ,VRS,TDS,QKSV,
!$   &                                     FAS,EAS,T0S,Q0S0,DTS_8,DTN_8,DQS_8,DQN_8)      
      DO K  = 1 ,NL
      DO J  = jbeg1,jend1
        FU1       = FRDU(J,1)
        FU2       = FRDU(J,2)
        FU3       = FRDU(J,3)
        FV1       = FRDV(J,1)
        FV2       = FRDV(J,2)
        FV3       = FRDV(J,3)
        FA1       = FRDP(J,1)
        FA2       = FRDP(J,2)
        FA3       = FRDP(J,3)
        DO I  = ibeg1,iend1
          II        = I + 1
          I1        = I - 1
          VR0       = VR  (I,J,K)
          DS0       = DS  (I,J,K)
          DT0       = DT  (I,J,K)
          DA0       = DA  (I,J,K)
          DB0       = DB  (I,J,K)
          RLNT0     = RLNT(I,J,K)
          ROT0      = ROT (I,J,K)
          RTA0      = DT0 * DA0 * ROT0
          RSB0      = DS0 * DB0 * RLNT0
          VRU       = VR0 + VR(I1,J,K)
          VRV       = VR0 + VR(I ,J+1,K)
          RTAI      = ROT (I1,J,K)*DT(I1,J,K)*DA(I1,J,K)
          RSBI      = RLNT(I,J-1,K)*DS(I,J-1,K)*DB(I,J-1,K)
          RTAJ      = ROT (I,J+1,K)*DT(I,J+1,K)*DA(I,J+1,K)
          RSBJ      = RLNT(II,J,K)*DS(II,J,K)*DB(II,J,K)
          FB1       = FA1 * RDLN(II,J,K)
          FB2       = FA1 * RDLN(I ,J,K)
          FB3       = FA2 * RDLT(I ,J,K)
          FB4       = FA3 * RDLT(I,J-1,K)
          TD0       = TK(I,J,K)
          Q0        = HQK(I,J,K)

          UIJK      = VRU * ( FU1*(RTA0-RTAI) + (FU2*RSB0 - FU3*RSBI) )
          VIJK      = VRV * ( FV1*(RSBJ-RSB0) - (FV2*RTAJ - FV3*RTA0) )
          TIJK      = VR0 * ( FB1*(TK(II,J,K)-TD0) - FB2*(TD0-TK(I1,J,K))
     &              +         FB3*(TK(I,J+1,K)-TD0) - FB4*(TD0-TK(I,J-1,K)) )
          QIJK      = VR0 * ( FB1*(HQK(II,J,K)-Q0 ) - FB2*(Q0 -HQK(I1,J,K))
     &              +         FB3*(HQK(I,J+1,K)-Q0 ) - FB4*(Q0 -HQK(I,J-1,K)) )
C
          U(I,J,K)  = HUK(I,J,K) + UIJK
          V(I,J,K)  = HVK(I,J,K) + VIJK
          T(I,J,K)  = TW(I,J,K) + TIJK
          QV(I,J,K)  = HQK (I,J,K) + QIJK
        ENDDO
      ENDDO


      IF (jbeg0 .eq. 1) THEN           ! North pole
      
      VRN       = VR(ibeg1,jbeg0,K )
      FAN       = FRDP(1 ,1) * VRN
      DTN_8     = DTN_16(K)
      T0N       = TW(ibeg1,jbeg0,K ) + DTN_8*FAN
      DQN_8     = DQN_16(K)
      Q0N       = HQK(ibeg1,jbeg0,K ) + DQN_8*FAN
      FVI       = FRDV(1,1)
      FVJ       = FRDV(1,2)      
      
      DO I  = ibeg1,iend1
        II        = I + 1
        T (I,1 ,K) = T0N
        QV(I,1 ,K) = Q0N
        U (I,1 ,K) = ZERO
        VRV       = VRN + VR(I,2,K)
        RSBI      = RLNT(II,1,K)*DS(II,1,K)*DB(II,1,K)
        RSB0      = RLNT(I ,1,K)*DS(I ,1,K)*DB(I ,1,K)
        RTA0      = ROT (I,2,K)*DT(I,2,K)*DA(I,2,K)
        VIJK      = VRV * ( FVI*(RSBI - RSB0) - FVJ*RTA0 )
        V(I,1 ,K) = HVK(I,1,K) + VIJK
      ENDDO

      ENDIF


      IF (jend0 .eq. PLAT) THEN   ! South pole
      VRS       = VR(ibeg1,jend0,K)
      FAS       = FRDP(PLAT,1) * VRS
      DTS_8     = DTS_16(K)
      T0S       = TW(ibeg1,jend0,K) + DTS_8*FAS
      DQS_8     = DQS_16(K)
      Q0S0      = HQK(ibeg1,jend0,K) + DQS_8*FAS
      DO I  = ibeg1,iend1
        T (I,PLAT,K) = T0S
        QV(I,PLAT,K) = Q0S0
        U (I,PLAT,K) = ZERO
        V (I,PLAT,K) = ZERO
      ENDDO
      
      ENDIF
C
C     ENERGY CONSERVATION DUE TO FRICTIONS
C
      DO J  = jbeg1,jend1
      DO I  = ibeg1,iend1
        DA(I,J,K)   = HUK(I,J,K) * (U(I,J,K) - HUK(I,J,K))
      ENDDO
      ENDDO

      DO J  = jbeg0 ,jend1
      DO I  = ibeg1,iend1
        DB(I,J,K)   = HVK(I,J,K) * (V(I,J,K) - HVK(I,J,K))
      ENDDO
      ENDDO
      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,DA(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_BOT,1,DB(:,beglatexdyn,1))

!$OMP PARALLEL DO PRIVATE (I,J,K,DXVPNJ,DXVPSJ)
      DO K  = 1 ,NL
      DO J  = jbeg1,jend1
        DXVPNJ    = DXVPN(J)
        DXVPSJ    = DXVPS(J)
        DO I  = ibeg1,iend1
          T(I,J,K)  = T(I,J,K) - (HALF * (DA(I,J,K) + DA(I+1,J,K))
     &              + DXVPNJ*DB(I,J-1,K) + DXVPSJ*DB(I,J,K))/CPD
        ENDDO
      ENDDO
      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,U(:,beglatexdyn,1),V(:,beglatexdyn,1),T(:,beglatexdyn,1),QV(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,U(:,beglatexdyn,1),V(:,beglatexdyn,1),T(:,beglatexdyn,1),QV(:,beglatexdyn,1))

      RETURN
      END

