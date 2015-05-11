#include <misc.h>
#include <params.h>

!!(2003.11.12-13)
!!(2004.03.30 wh)
!!-------------------

#if (defined SPMD)

*     SUBROUTINE HDIFUS(U,V,T,QV,QC,QR,QI,QS,QG
      SUBROUTINE HDIFUS(U,V,T,QV
     _                 ,FRDT,FRDS,FRDU,FRDV,FRDP,TB,PLY,DXVPN,DXVPS)
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
      use mpishorthand, only: mpicom

      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>
#include <commpi.h>
#include <mpif.h>

      real*8 workst1(  nx*nl),workrb1(  nx*nl),worksb1(3*nx*nl),workrt1(3*nx*nl)
!!    real workst2(6*nx*nl),workrb2(6*nx*nl),worksb2(6*nx*nl),workrt2(6*nx*nl) !!( 2004.03.30)
      real*8 workst2(  nx*nl),workrb2(  nx*nl),worksb2(  nx*nl),workrt2(  nx*nl) !!
      real*8 workst3(3*nx*nl),workrb3(3*nx*nl),worksb3(  nx*nl),workrt3(  nx*nl)
      real*8                                   worksb4(  nx*nl),workrt4(  nx*nl)

      integer isend11,irecv11
      integer isend12,irecv12
      integer isend21,irecv21
      integer isend22,irecv22
      integer isend31,irecv31
      integer isend32,irecv32
      integer isend42,irecv42

      integer istatus(mpi_status_size)
!      integer ii

      character*50 filename

C
      REAL*8 FRDT(NY,3),FRDS(NY,3),FRDU(NY,3),FRDV(NY,3),FRDP(NY,3)
C
      REAL*8 U(NX,NY,NL),V(NX,NY,NL),T(NX,NY,NL)

!!    REAL QV(NX,NY,NL),QC(NX,NY,NL),QR(NX,NY,NL)
!!   _    ,QI(NX,NY,NL),QS(NX,NY,NL),QG(NX,NY,NL)
      REAL*8 QV(NX,NY,NL)                              !!(wh 2004.03.30)

      REAL*8 TB(NX,NY,NL),PLY(NX,NY,NZ),DXVPN(NY),DXVPS(NY)
C
      REAL*8         D(NX,NY,nl),DT(NX,NY,nl),DS(NX,NY,nl),DA(NX,NY,nl),DB(NX,NY,nl)
     &       ,      VR(NX,NY,nl),QK(NX,NY,nl),TK(NX,NY,nl),VK(NX,NY,nl),UK(NX,NY,nl)
     &       ,      ROT(NX,NY,nl), RLNT(NX,NY,nl), RDLN(NX,NY,nl), RDLT(NX,NY,nl)
     &       ,      TW (NX,NY,nl)                    !!(wh 2004.03.30)
!!   &       ,      TW (NX,NY,nl),QKC(NX,NY,nl),QKR(NX,NY,nl),QKI(NX,NY,nl)
!!   &       ,      QKS(NX,NY,nl),QKG(NX,NY,nl)
C
      REAL*8  ZERO,HALF,ONE
      DATA       ZERO,HALF,ONE / 0.0E0 , 0.5E0 , 1.0E0  /
      REAL*8  FRDTN,FRDTS,FRDSI,FRDSJ,TI,RIJ,RI,VRI,TDI,FT1,FT2,FT3
      real*8    FS1,FS2,FS3,U0,V0,DTN,DTS,DTJ,DSI,DIJ,R0,R1,D0,D1,RT
      real*8    RN,DIJ1,DIJ2,FU1,FU2,FU3,FV1,FV2,FV3,FA1,FA2,FA3,VR0
      real*8    DS0,DT0,DA0,DB0,RLNT0,ROT0,RTA0,RSB0,VRU,VRV,RTAI,RSBI
!!   &       ,RTAJ,RSBJ,FB1,FB2,FB3,FB4,TD0,Q0,Q0C,Q0R,Q0I,Q0S,Q0G
!!   &       ,UIJK,VIJK,TIJK,QIJK,QIJKC,QIJKR,QIJKI,QIJKS,QIJKG
!!   &       ,VRN,VRS,TDN,TDS,QKN,QKNC,QKNR,QKNI,QKNS,QKNG
!!   &       ,QKSV,QKSC,QKSR,QKSI,QKSS,QKSG,DQN,DQNC,DQNR,DQNI,DQNS,DQNG
!!   &       ,DQS,DQSC,DQSR,DQSI,DQSS,DQSG,FAN,FAS,EAN,EAS,T0N
!!   &       ,T0S,Q0N,Q0NC,Q0NR,Q0NI,Q0NS,Q0NG
!!   &       ,Q0S0,Q0SC,Q0SR,Q0SI,Q0SS,Q0SG,FVI,FVJ,DXVPNJ,DXVPSJ
!!                                                                !!(wh 2004.03.30)
     &       ,RTAJ,RSBJ,FB1,FB2,FB3,FB4,TD0,Q0
     &       ,UIJK,VIJK,TIJK,QIJK
     &       ,VRN,VRS,TDN,TDS,QKN
     &       ,QKSV,DQN
     &       ,DQS,FAN,FAS,EAN,EAS,T0N
     &       ,T0S,Q0N
     &       ,Q0S0,FVI,FVJ,DXVPNJ,DXVPSJ

      INTEGER I,J,K,J1,JJ,I1,II
CWB
C
      if (myrank.eq.nprocs-1) then    ! the north pole
          FRDTN     = FRDT(jpole ,1)
          FRDSI     = FRDS(jpole ,1)
          FRDSJ     = FRDS(jpole ,2)
      endif

      if (myrank.eq.0) then           ! the south pole
          FRDTS     = FRDT(jpole,1)
      endif
!
C
C     CALCULATE DENSITY & SET 2-D FIELDS TO BE DIFFUSED
C
      do 300 k = 1 ,nl
      DO 300 J = 2 ,NY-1
      DO 300 I = 1 ,NX
      UK(I,J,k)  = U(I,J,K)
      VK(I,J,k)  = V(I,J,K)
      QK(I,J,k)  = QV(I,J,K)
*     QKC(I,J)  = QC(I,J,K)
*     QKR(I,J)  = QR(I,J,K)
*     QKI(I,J)  = QI(I,J,K)
*     QKS(I,J)  = QS(I,J,K)
*     QKG(I,J)  = QG(I,J,K)
!!      QKC(I,J,k)  = 0.0
!!      QKR(I,J,k)  = 0.0
!!      QKI(I,J,k)  = 0.0  !!(wh 2004.03.30)
!!      QKS(I,J,k)  = 0.0
!!      QKG(I,J,k)  = 0.0
      TW(I,J,k)  = T(I,J,K)
300   CONTINUE

      ii=1
      do k=1,nl
       do i=1,nx
          workst1(ii)=uk(i,2,   k)
          worksb1(ii)=vk(i,ny-1,k)
          ii=ii+1
       enddo
      enddo

      call mpi_isend( workst1,nx*nl,mpi_double_precision,itop,11,mpicom,isend11,ierr)
      call mpi_irecv( workrb1,nx*nl,mpi_double_precision,ibot,11,mpicom,irecv11,ierr)
      call mpi_isend( worksb1,nx*nl,mpi_double_precision,ibot,12,mpicom,isend12,ierr)
      call mpi_irecv( workrt1,nx*nl,mpi_double_precision,itop,12,mpicom,irecv12,ierr)

      ii=1
!--1
      do k=1,nl
       do i=1,nx
          workst2(ii)=qk(i,2,   k)
          worksb2(ii)=qk(i,ny-1,k)
          ii=ii+1
       enddo
      enddo
!!--2
!!      do k=1,nl
!!       do i=1,nx
!!          workst2(ii)=qkc(i,2,   k)
!!          worksb2(ii)=qkc(i,ny-1,k)
!!          ii=ii+1
!!       enddo
!!      enddo
!!!--3
!!      do k=1,nl
!!       do i=1,nx
!!          workst2(ii)=qkr(i,2,   k)
!!          worksb2(ii)=qkr(i,ny-1,k)
!!          ii=ii+1
!!       enddo
!!      enddo
!!!--4
!!      do k=1,nl
!!       do i=1,nx
!!          workst2(ii)=qki(i,2,   k)
!!          worksb2(ii)=qki(i,ny-1,k)
!!          ii=ii+1
!!       enddo
!!      enddo
!!!!--5
!!      do k=1,nl
!!       do i=1,nx
!!          workst2(ii)=qks(i,2,   k)
!!          worksb2(ii)=qks(i,ny-1,k)
!!          ii=ii+1
!!       enddo
!!      enddo
!!!--6
!!      do k=1,nl
!!       do i=1,nx
!!          workst2(ii)=qkg(i,2,   k)
!!          worksb2(ii)=qkg(i,ny-1,k)
!!          ii=ii+1
!!       enddo
!!      enddo
!!
!!      call mpi_isend( workst2,6*nx*nl,mpi_double_precision,itop,21,mpicom,isend21,ierr)
!!      call mpi_irecv( workrb2,6*nx*nl,mpi_double_precision,ibot,21,mpicom,irecv21,ierr)
!!      call mpi_isend( worksb2,6*nx*nl,mpi_double_precision,ibot,22,mpicom,isend22,ierr)
!!      call mpi_irecv( workrt2,6*nx*nl,mpi_double_precision,itop,22,mpicom,irecv22,ierr)

      call mpi_isend( workst2,nx*nl,mpi_double_precision,itop,21,mpicom,isend21,ierr)
      call mpi_irecv( workrb2,nx*nl,mpi_double_precision,ibot,21,mpicom,irecv21,ierr)
      call mpi_isend( worksb2,nx*nl,mpi_double_precision,ibot,22,mpicom,isend22,ierr)
      call mpi_irecv( workrt2,nx*nl,mpi_double_precision,itop,22,mpicom,irecv22,ierr)
!-----------------------------------

      do k = 1,nl

         DO 325 J  = JBeg,JEnd
         DO 325 I  = IB,IE
            TI        = TW(I,J,k)
            RIJ       = PLY(I,J,K) / (TI*RD)
            ROT(I,J,k)  = RIJ
            VR (I,J,k)  = ONE / RIJ
            TK (I,J,k)  = TI  - TB(I,J,K)
325      CONTINUE

         if (inc_pole) then
            TI        = TW (IB,Jpole,k)
            RI        = PLY(IB,Jpole,K) / (TI*RD)
            VRI       = ONE / RI
            TDI       = TI  - TB(IB,Jpole,K)
            DO 350 I  = IB,IE
            ROT(I,Jpole,k)  = RI
            VR (I,Jpole,k)  = VRI
            TK (I,Jpole,k)  = TDI
350      CONTINUE
         endif

         DO 375 J  = 2 ,NY-1
            ROT(1 ,J,k) = ROT(IE,J,k)
            ROT(NX,J,k) = ROT(IB,J,k)
            VR (1 ,J,k) = VR (IE,J,k)
            VR (NX,J,k) = VR (IB,J,k)
            TK (1 ,J,k) = TK (IE,J,k)
            TK (NX,J,k) = TK (IB,J,k)
375      CONTINUE
      enddo
!---
      ii=1
      do k=1,nl
       do i=1,nx
          workst3(ii)=tk(i,2,   k)
          worksb3(ii)=tk(i,ny-1,k)
          ii=ii+1
       enddo
      enddo
      do k=1,nl
       do i=1,nx
          workst3(ii)=vr(i,2,   k)
          ii=ii+1
       enddo
      enddo
      do k=1,nl
       do i=1,nx
          workst3(ii)=rot(i,2,   k)
          ii=ii+1
       enddo
      enddo

      call mpi_isend( workst3,3*nx*nl,mpi_double_precision,itop,31,mpicom,isend31,ierr)
      call mpi_irecv( workrb3,3*nx*nl,mpi_double_precision,ibot,31,mpicom,irecv31,ierr)
      call mpi_isend( worksb3,  nx*nl,mpi_double_precision,ibot,32,mpicom,isend32,ierr)
      call mpi_irecv( workrt3,  nx*nl,mpi_double_precision,itop,32,mpicom,irecv32,ierr)
!----------------------------------
C
C     CALCULATE DEFORMATION FIELDS
C

      call mpi_wait(isend11,istatus,ierr)
      call mpi_wait(irecv11,istatus,ierr)
      call mpi_wait(isend12,istatus,ierr)
      call mpi_wait(irecv12,istatus,ierr)
      ii=1
      do k=1,nl
       do i=1,nx
          uk(i,ny,k)=workrb1(ii)
          vk(i,1, k)=workrt1(ii)
          ii=ii+1
       enddo
      enddo
!--

      do k=1,nl

      DO 400 J  = JBeg,JEnd
      J1        = J - 1
      JJ        = J + 1
      FT1       = FRDT(J,1)
      FT2       = FRDT(J,2)
      FT3       = FRDT(J,3)
      FS1       = FRDS(J,1)
      FS2       = FRDS(J,2)
      FS3       = FRDS(J,3)
      DO 400 I  = IB,IE
      U0        = UK(I,J,k)
      V0        = VK(I,J,k)
      DT(I,J,k)   = FT1*(UK(I+1,J,k)-U0) - (FT2*V0 - FT3*VK(I,J1,k ))
      DS(I,J,k)   = FS1*(V0-VK(I-1,J,k)) + (FS2*UK(I,JJ,k) - FS3*U0)
400   CONTINUE

      if (myrank.eq.nprocs-1) then
          DTN       = ZERO
          do i=ib,ie
            DTN       = DTN + VK(I,jpole,k)
          enddo
            DTN       = DTN * FRDTN
          do i=ib,ie
            DT(I,jpole,k) = DTN
            DS(I,jpole,k) = FRDSI*(VK(I,jpole,k) -VK(I-1,jpole,k))+FRDSJ*UK(I,jpole+1,k)
          enddo
      endif

      if (myrank.eq.0) then
         DTS       = ZERO
         do i=ib,ie
            DTS       = DTS + VK(I,jpole-1,k)
         enddo
            DTS       = DTS * FRDTS
         DO  I  = IB,IE
            DT(I,jpole,k)  = DTS
            DS(I,jpole,k)  = ZERO
         enddo
      endif

      DO 475 J  = 2 ,NY-1
      DT(1 ,J,k)  = DT(IE,J,k)
      DT(NX,J,k)  = DT(IB,J,k)
      DS(1 ,J,k)  = DS(IE,J,k)
      DS(NX,J,k)  = DS(IB,J,k)
475   CONTINUE

      enddo
!------------------ send DS ----
      ii=1
      do k=1,nl
       do i=1,nx
          worksb4(ii)=ds(i,ny-1,k)
          ii=ii+1
       enddo
      enddo

      call mpi_isend( worksb4,nx*nl,mpi_double_precision,ibot,42,mpicom,isend42,ierr)
      call mpi_irecv( workrt4,nx*nl,mpi_double_precision,itop,42,mpicom,irecv42,ierr)

!----------------------send & recv DT-------
      ii=1
      do k=1,nl
       do i=1,nx
          workst1(ii)=dt(i,2,k)
          ii=ii+1
       enddo
      enddo

      call mpi_isend( workst1,nx*nl,mpi_double_precision,itop,11,mpicom,isend11,ierr)
      call mpi_irecv( workrb1,nx*nl,mpi_double_precision,ibot,11,mpicom,irecv11,ierr)
      call mpi_wait(isend11,istatus,ierr)
      call mpi_wait(irecv11,istatus,ierr)

      ii=1
      do k=1,nl
       do i=1,nx
          dt(i,ny,k)=workrb1(ii)
          ii=ii+1
       enddo
      enddo
!-----

      do k=1,nl

      DO 500 J  = 2 ,JEnd
      JJ        = J + 1
      DO 525 I  = IB,IE
      DTJ       = DT(I,JJ,k) + DT(I,J,k)
      DSI       = DS(I,J ,k) + DS(I-1,J,k)
      DIJ       = DTJ*DTJ  + DSI*DSI
      D(I,J,k)    = HALF * SQRT( DIJ )
525   CONTINUE
      D(1 ,J,k)   = D(IE,J,k)
      D(NX,J,k)   = D(IB,J,k)
500   CONTINUE

      if (myrank.eq.0) then
          DO I = 1 ,NX
            D(I,jpole,k) = ZERO
          ENDDO
      endif

      enddo

!------------------ send & recv D ----
      ii=1
      do k=1,nl
       do i=1,nx
          worksb1(ii)=d(i,ny-1,k)
          ii=ii+1
       enddo
      enddo

      call mpi_isend( worksb1,nx*nl,mpi_double_precision,ibot,12,mpicom,isend12,ierr)
      call mpi_irecv( workrt1,nx*nl,mpi_double_precision,itop,12,mpicom,irecv12,ierr)
      call mpi_wait( isend12,istatus,ierr)
      call mpi_wait( irecv12,istatus,ierr)

      ii=1
      do k=1,nl
       do i=1,nx
          d(i,1,k)=workrt1(ii)
          ii=ii+1
       enddo
      enddo
!------------------- recv TK,VR,ROT-------

      call mpi_wait(isend31,istatus,ierr)
      call mpi_wait(irecv31,istatus,ierr)
      call mpi_wait(isend32,istatus,ierr)
      call mpi_wait(irecv32,istatus,ierr)

      ii=1
      do k=1,nl
       do i=1,nx
          tk(i,ny,k)=workrb3(ii)
          tk(i,1, k)=workrt3(ii)
          ii=ii+1
       enddo
      enddo
      do k=1,nl
       do i=1,nx
          vr(i,ny,k)=workrb3(ii)
          ii=ii+1
       enddo
      enddo
      do k=1,nl
       do i=1,nx
          rot(i,ny,k)=workrb3(ii)
          ii=ii+1
       enddo
      enddo
!-----------------------------------------------------------

      do k=1,nl

      DO 550 J  = JBeg,JEnd
      JJ        = J + 1
      J1        = J - 1
      DO 550 I  = IB,IE
      I1        = I - 1
      R0        = ROT(I,J,k)
      R1        = ROT(I,JJ,k)
      D0        = D(I,J,k)
      D1        = D(I,J1,k)
      RT        = R0 + R1
      RN        = R0 + ROT(I1,J,k)
      RLNT(I,J,k) = RN + R1 + ROT(I1,JJ,k)
      DIJ1      = HALF * (D0 + D(I1,J,k))
      DIJ       = HALF * (D0 + D1)
      DIJ2      = HALF * (DIJ1 + HALF*(D1 + D(I1,J1,k)) )
      RDLN(I,J,k) = RN   *  DIJ2
      RDLT(I,J,k) = RT   *  D0
      DA  (I,J,k) = DIJ
      DB  (I,J,k) = DIJ1
550   CONTINUE

      if (inc_pole) then
        DO I = 1 ,NX
          RDLN(I,Jpole,k) = ZERO
          DA(I,Jpole,k)   = ZERO
        ENDDO
      endif

      if (myrank.eq.nprocs-1) then
         R0        = ROT (IB,jpole,k)
         DO 560 I  = IB,IE
            D0        = D(I,jpole,k)
            DB  (I,jpole,k) = (D0 + D(I-1,jpole ,k)) * HALF
            RDLT(I,jpole,k) = (R0 + ROT(I,Jpole+1,k)) * D0
            RLNT(I,jpole,k) =  R0 + R0 + ROT(I-1,Jpole+1,k) + ROT(I,Jpole+1,k)
560      CONTINUE
      endif

      if (myrank.eq.0) then
         DO I = 1 ,NX
           RDLT(I,jpole,k) = ZERO
           RLNT(I,jpole,k) = ZERO
           DB(I,jpole,k)   = ZERO
         ENDDO
      endif

      DO 580 J  = 2 ,JEnd
      RLNT(1 ,J,k)= RLNT(IE,J,k)
      RLNT(NX,J,k)= RLNT(IB,J,k)
      RDLN(1 ,J,k)= RDLN(IE,J,k)
      RDLN(NX,J,k)= RDLN(IB,J,k)
      RDLT(1 ,J,k)= RDLT(IE,J,k)
      RDLT(NX,J,k)= RDLT(IB,J,k)
      DA  (1 ,J,k)= DA  (IE,J,k)
      DA  (NX,J,k)= DA  (IB,J,k)
      DB  (1 ,J,k)= DB  (IE,J,k)
      DB  (NX,J,k)= DB  (IB,J,k)
580   CONTINUE

      enddo
!---------------------------------------------------------------------

C     UPDATE T & U & V & QV,QC,QI,QR,QS,QG DUE TO THE HORIZONTAL DIFFUSION
C

      ii=1
      do k=1,nl
       do i=1,nx
          workst1(ii)=da(i,2,k)
          worksb1(ii)=db(i,ny-1,k)
          ii=ii+1
       enddo
      enddo
      do k=1,nl
       do i=1,nx
          worksb1(ii)=rdlt(i,ny-1,k)
          ii=ii+1
       enddo
      enddo
      do k=1,nl
       do i=1,nx
          worksb1(ii)=rlnt(i,ny-1,k)
          ii=ii+1
       enddo
      enddo

      call mpi_isend( workst1,  nx*nl,mpi_double_precision,itop,11,mpicom,isend11,ierr)
      call mpi_irecv( workrb1,  nx*nl,mpi_double_precision,ibot,11,mpicom,irecv11,ierr)
      call mpi_isend( worksb1,3*nx*nl,mpi_double_precision,ibot,12,mpicom,isend12,ierr)
      call mpi_irecv( workrt1,3*nx*nl,mpi_double_precision,itop,12,mpicom,irecv12,ierr)
      call mpi_wait(isend11,istatus,ierr)
      call mpi_wait(irecv11,istatus,ierr)
      call mpi_wait(isend12,istatus,ierr)
      call mpi_wait(irecv12,istatus,ierr)

      ii=1
      do k=1,nl
       do i=1,nx
          da(i,ny,k)=workrb1(ii)
          db(i,1, k)=workrt1(ii)
          ii=ii+1
       enddo
      enddo
      do k=1,nl
       do i=1,nx
          rdlt(i,1,k)=workrt1(ii)
          ii=ii+1
       enddo
      enddo
      do k=1,nl
       do i=1,nx
          rlnt(i,1,k)=workrt1(ii)
          ii=ii+1
       enddo
      enddo
!-----
      call mpi_wait(isend42,istatus,ierr)
      call mpi_wait(irecv42,istatus,ierr)
      ii=1
      do k=1,nl
       do i=1,nx
          ds(i,1,k)=workrt4(ii)
          ii=ii+1
       enddo
      enddo

!-------------
      call mpi_wait(isend21,istatus,ierr)
      call mpi_wait(irecv21,istatus,ierr)
      call mpi_wait(isend22,istatus,ierr)
      call mpi_wait(irecv22,istatus,ierr)

      ii=1
!--1
      do k=1,nl
       do i=1,nx
          qk(i,ny,k)=workrb2(ii)
          qk(i,1, k)=workrt2(ii)
          ii=ii+1
       enddo
      enddo
!!--2
!!      do k=1,nl
!!       do i=1,nx
!!          qkc(i,ny,k)=workrb2(ii)
!!          qkc(i,1, k)=workrt2(ii)
!!          ii=ii+1
!!       enddo
!!      enddo
!!--3
!!      do k=1,nl
!!!       do i=1,nx
!!          qkr(i,ny,k)=workrb2(ii)
!!          qkr(i,1, k)=workrt2(ii)
!!          ii=ii+1
!!       enddo
!!      enddo
!!!--4
!!      do k=1,nl
!!       do i=1,nx
!!          qki(i,ny,k)=workrb2(ii)
!!          qki(i,1, k)=workrt2(ii)
!!          ii=ii+1
!!       enddo
!!      enddo
!!!--5
!!      do k=1,nl
!!       do i=1,nx
!!          qks(i,ny,k)=workrb2(ii)
!!          qks(i,1, k)=workrt2(ii)
!!          ii=ii+1
!!       enddo
!!      enddo
!!!--6
!!      do k=1,nl
!!       do i=1,nx
!!          qkg(i,ny,k)=workrb2(ii)
!!          qkg(i,1, k)=workrt2(ii)
!!          ii=ii+1
!!       enddo
!!      enddo                   !!(wh 2004.03.30)
!!!-------------------------------------------------
      do k=1,nl

      DO 600 J  = JBeg,JEnd
      FU1       = FRDU(J,1)
      FU2       = FRDU(J,2)
      FU3       = FRDU(J,3)
      FV1       = FRDV(J,1)
      FV2       = FRDV(J,2)
      FV3       = FRDV(J,3)
      FA1       = FRDP(J,1)
      FA2       = FRDP(J,2)
      FA3       = FRDP(J,3)
      JJ        = J + 1
      J1        = J - 1
      DO 600 I  = IB,IE
      II        = I + 1
      I1        = I - 1
      VR0       = VR  (I,J,k)
      DS0       = DS  (I,J,k)
      DT0       = DT  (I,J,k)
      DA0       = DA  (I,J,k)
      DB0       = DB  (I,J,k)
      RLNT0     = RLNT(I,J,k)
      ROT0      = ROT (I,J,k)
      RTA0      = DT0 * DA0 * ROT0
      RSB0      = DS0 * DB0 * RLNT0
      VRU       = VR0 + VR(I1,J,k)
      VRV       = VR0 + VR(I ,JJ,k)
      RTAI      = ROT (I1,J,k)*DT(I1,J,k)*DA(I1,J,k)
      RSBI      = RLNT(I,J1,k)*DS(I,J1,k)*DB(I,J1,k)
      RTAJ      = ROT (I,JJ,k)*DT(I,JJ,k)*DA(I,JJ,k)
      RSBJ      = RLNT(II,J,k)*DS(II,J,k)*DB(II,J,k)
      FB1       = FA1 * RDLN(II,J,k)
      FB2       = FA1 * RDLN(I ,J,k)
      FB3       = FA2 * RDLT(I ,J,k)
      FB4       = FA3 * RDLT(I,J1,k)
      TD0       = TK(I,J,k)
      Q0        = QK(I,J,k)
!!      Q0C       = QKC(I,J,k)
!!      Q0R       = QKR(I,J,k)
!!      Q0I       = QKI(I,J,k)  !!(wh 2004.03.30)
!!      Q0S       = QKS(I,J,k)
!!      Q0G       = QKG(I,J,k)

      UIJK      = VRU * ( FU1*(RTA0-RTAI) + (FU2*RSB0 - FU3*RSBI) )
      VIJK      = VRV * ( FV1*(RSBJ-RSB0) - (FV2*RTAJ - FV3*RTA0) )
      TIJK      = VR0 * ( FB1*(TK(II,J,k)-TD0) - FB2*(TD0-TK(I1,J,k))
     &          +         FB3*(TK(I,JJ,k)-TD0) - FB4*(TD0-TK(I,J1,k)) )
      QIJK      = VR0 * ( FB1*(QK(II,J,k)-Q0 ) - FB2*(Q0 -QK(I1,J,k))
     &          +         FB3*(QK(I,JJ,k)-Q0 ) - FB4*(Q0 -QK(I,J1,k)) )
!!      QIJKC     = VR0 * ( FB1*(QKC(II,J,k)-Q0C) - FB2*(Q0C-QKC(I1,J,k))      !!
!!     &          +         FB3*(QKC(I,JJ,k)-Q0C) - FB4*(Q0C-QKC(I,J1,k)) )    !!
!!      QIJKR     = VR0 * ( FB1*(QKR(II,J,k)-Q0R) - FB2*(Q0R-QKR(I1,J,k))      !!
!!     &          +         FB3*(QKR(I,JJ,k)-Q0R) - FB4*(Q0R-QKR(I,J1,k)) )    !!
!!      QIJKI     = VR0 * ( FB1*(QKI(II,J,k)-Q0I) - FB2*(Q0I-QKI(I1,J,k))      !!(wh 2004.03.30)
!!     &          +         FB3*(QKI(I,JJ,k)-Q0I) - FB4*(Q0I-QKI(I,J1,k)) )    !!
!!      QIJKS     = VR0 * ( FB1*(QKS(II,J,k)-Q0S) - FB2*(Q0S-QKS(I1,J,k))      !!
!!     &          +         FB3*(QKS(I,JJ,k)-Q0S) - FB4*(Q0S-QKS(I,J1,k)) )    !!
!!      QIJKG     = VR0 * ( FB1*(QKG(II,J,k)-Q0G) - FB2*(Q0G-QKG(I1,J,k))      !!
!!     &          +         FB3*(QKG(I,JJ,k)-Q0G) - FB4*(Q0G-QKG(I,J1,k)) )    !!
C
      U(I,J,K)  = UK(I,J,k) + UIJK
      V(I,J,K)  = VK(I,J,k) + VIJK
      T(I,J,K)  = TW(I,J,k) + TIJK
      QV(I,J,K)  = QK (I,J,k) + QIJK
!!      QC(I,J,K)  = QKC(I,J,k) + QIJKC !!
!!      QR(I,J,K)  = QKR(I,J,k) + QIJKR !!
!!      QI(I,J,K)  = QKI(I,J,k) + QIJKI !! (wh 2004.03.30)
!!      QS(I,J,K)  = QKS(I,J,k) + QIJKS !!
!!      QG(I,J,K)  = QKG(I,J,k) + QIJKG !!
600   CONTINUE

      if (myrank.eq.nprocs-1) then
         VRN       = VR(IB,jpole,k)
         TDN       = TK(IB,jpole,k)
         QKN       = QK(IB,jpole,k)
!!         QKNC      = QKC(IB,jpole,k)
!!         QKNR      = QKR(IB,jpole,k)
!!         QKNI      = QKI(IB,jpole,k)  !!(wh 2004.03.30)
!!         QKNS      = QKS(IB,jpole,k)
!!         QKNG      = QKG(IB,jpole,k)
         DQN       = ZERO
!!         DQNC      = ZERO
!!         DQNR      = ZERO
!!         DQNI      = ZERO             !!(wh 2004.03.30)
!!         DQNS      = ZERO
!!         DQNG      = ZERO
         DTN       = ZERO
         FAN       = FRDP(jpole ,1) * VRN
         DO  I  = IB,IE
            EAN       = RDLT(I,jpole,k)
            DTN       = DTN + EAN*(TK(I,Jpole+1,k) - TDN)
            DQN       = DQN + EAN*(QK(I,Jpole+1,k) - QKN)
!!            DQNC      = DQNC+ EAN*(QKC(I,Jpole+1,k) - QKNC)
!!            DQNR      = DQNR+ EAN*(QKR(I,Jpole+1,k) - QKNR)
!!            DQNI      = DQNI+ EAN*(QKI(I,Jpole+1,k) - QKNI)  !(wh 2004.03.30)
!!            DQNS      = DQNS+ EAN*(QKS(I,Jpole+1,k) - QKNS)
!!            DQNG      = DQNG+ EAN*(QKG(I,Jpole+1,k) - QKNG)
         enddo
      endif

      if (myrank.eq.0) then
         VRS       = VR(IB,jpole,k)
         TDS       = TK(IB,jpole,k)
         QKSV      = QK(IB,jpole,k)
!!         QKSC      = QKC(IB,jpole,k)
!!         QKSR      = QKR(IB,jpole,k)
!!         QKSI      = QKI(IB,jpole,k)  !!(wh 2004.03.30)
!!         QKSS      = QKS(IB,jpole,k)
!!         QKSG      = QKG(IB,jpole,k)
         DQS       = ZERO
!!         DQSC      = ZERO
!!         DQSR      = ZERO
!!         DQSI      = ZERO
!!         DQSS      = ZERO
!!         DQSG      = ZERO
         DTS       = ZERO
         FAS       = FRDP(jpole,1) * VRS
         DO  I  = IB,IE
            EAS       = RDLT(I,Jpole-1,k)
            DTS       = DTS + EAS*(TDS - TK(I,Jpole-1,k))
            DQS       = DQS + EAS*(QKSV- QK(I,Jpole-1,k))
!!            DQSC      = DQSC+ EAS*(QKSC - QKC(I,Jpole-1,k))
!!            DQSR      = DQSR+ EAS*(QKSR - QKR(I,Jpole-1,k))
!!            DQSI      = DQSI+ EAS*(QKSI - QKI(I,Jpole-1,k))
!!            DQSS      = DQSS+ EAS*(QKSS - QKS(I,Jpole-1,k))
!!            DQSG      = DQSG+ EAS*(QKSG - QKG(I,Jpole-1,k))
         enddo
      endif

      if (myrank.eq.nprocs-1) then
         T0N       = TW(IB,jpole,k) + DTN*FAN
         Q0N       = QK(IB,jpole,k) + DQN*FAN
!!         Q0NC      = QKC(IB,jpole,k) + DQNC*FAN
!!         Q0NR      = QKR(IB,jpole,k) + DQNR*FAN
!!         Q0NI      = QKI(IB,jpole,k) + DQNI*FAN
!!         Q0NS      = QKS(IB,jpole,k) + DQNS*FAN
!!         Q0NG      = QKG(IB,jpole,k) + DQNG*FAN
         FVI       = FRDV(jpole,1)
         FVJ       = FRDV(jpole,2)
         DO  I  = IB,IE
            II        = I + 1
            T (I,jpole,K) = T0N
            QV(I,jpole,K) = Q0N
!!            QC(I,jpole,K) = Q0NC
!!            QR(I,jpole,K) = Q0NR
!!            QI(I,jpole,K) = Q0NI
!!            QS(I,jpole,K) = Q0NS
!!            QG(I,jpole,K) = Q0NG
            U (I,jpole,K) = ZERO
            VRV       = VRN + VR(I,Jpole+1,k)
            RSBI      = RLNT(II,jpole,k)*DS(II,jpole,k)*DB(II,jpole,k)
            RSB0      = RLNT(I ,jpole,k)*DS(I ,jpole,k)*DB(I ,jpole,k)
            RTA0      = ROT (I,jpole+1,k)*DT(I,jpole+1,k)*DA(I,jpole+1,k)
            VIJK      = VRV * ( FVI*(RSBI - RSB0) - FVJ*RTA0 )
            V(I,jpole,K) = VK(I,jpole,k) + VIJK
         enddo
      endif

      if (myrank.eq.0) then
         T0S       = TW(IB,jpole,k) + DTS*FAS
         Q0S0      = QK(IB,jpole,k) + DQS*FAS
!!         Q0SC      = QKC(IB,jpole,k) + DQSC*FAS
!!         Q0SR      = QKR(IB,jpole,k) + DQSR*FAS
!!         Q0SI      = QKI(IB,jpole,k) + DQSI*FAS
!!         Q0SS      = QKS(IB,jpole,k) + DQSS*FAS
!!         Q0SG      = QKG(IB,jpole,k) + DQSG*FAS
         DO  I  = IB,IE
            II        = I + 1
            T (I,jpole,K) = T0S
            QV(I,jpole,K) = Q0S0
!!            QC(I,jpole,K) = Q0SC
!!            QR(I,jpole,K) = Q0SR
!!            QI(I,jpole,K) = Q0SI
!!            QS(I,jpole,K) = Q0SS
!!            QG(I,jpole,K) = Q0SG
            U (I,jpole,K) = ZERO
            V (I,jpole,K) = ZERO
         enddo
      endif

      enddo
!------
C
C     ENERGY CONSERVATION DUE TO FRICTIONS
C
      do k=1,nl
       DO 710 J  = JBeg,JEnd
         DO 700 I  = IB,IE
           DA(I,J,k)   = UK(I,J,k) * (U(I,J,K) - UK(I,J,k))
700      CONTINUE
CWB
           DA(1 ,J,k)  = DA(IE,J,k)
           DA(NX,J,k)  = DA(IB,J,k)
CWB
710    CONTINUE
      enddo
!----
      do 720 k  = 1,nl
      DO 720 J  = 2 ,JEnd
      DO 720 I  = IB,IE
      DB(I,J,k)   = VK(I,J,k) * (V(I,J,K) - VK(I,J,k))
720   CONTINUE

!------------------ send & recv DB ----
      ii=1
      do k=1,nl
       do i=1,nx
          worksb1(ii)=db(i,ny-1,k)
          ii=ii+1
       enddo
      enddo

      call mpi_isend( worksb1,nx*nl,mpi_double_precision,ibot,12,mpicom,isend12,ierr)
      call mpi_irecv( workrt1,nx*nl,mpi_double_precision,itop,12,mpicom,irecv12,ierr)
      call mpi_wait( isend12,istatus,ierr)
      call mpi_wait( irecv12,istatus,ierr)

      ii=1
      do k=1,nl
       do i=1,nx
          db(i,1,k)=workrt1(ii)
          ii=ii+1
       enddo
      enddo
!----
      do k=1,nl
      DO 750 J  = JBeg,JEnd
      JJ        = J - 1
      DXVPNJ    = DXVPN(J)
      DXVPSJ    = DXVPS(J)
      DO 750 I  = IB,IE
      D (I,J,k)   = HALF * (DA(I,J,k) + DA(I+1,J,k))
     &          + DXVPNJ*DB(I,JJ,k) + DXVPSJ*DB(I,J,k)
      TK(I,J,k)   = T(I,J,K)
750   CONTINUE
      DO 800 J  = JBeg,JEnd
      DO 800 I  = IB,IE
      T(I,J,K)  = TK(I,J,k) - D(I,J,k)/CPD
800   CONTINUE
      enddo
C
C     SPHERICAL BOUNDARY CONDITIONS
C
      DO 900 K  = 1 ,NL
      DO 900 J  = 2 ,NY-1
      U(1 ,J,K) = U(IE,J,K)
      U(NX,J,K) = U(IB,J,K)
      V(1 ,J,K) = V(IE,J,K)
      V(NX,J,K) = V(IB,J,K)
      T(1 ,J,K) = T(IE,J,K)
      T(NX,J,K) = T(IB,J,K)
      QV(1 ,J,K) = QV(IE,J,K)
      QV(NX,J,K) = QV(IB,J,K)
!!      QC(1 ,J,K) = QC(IE,J,K)
!!      QC(NX,J,K) = QC(IB,J,K)
!!      QR(1 ,J,K) = QR(IE,J,K)
!!      QR(NX,J,K) = QR(IB,J,K)
!!      QI(1 ,J,K) = QI(IE,J,K)
!!      QI(NX,J,K) = QI(IB,J,K)  !!(wh 2004.03.30)
!!      QS(1 ,J,K) = QS(IE,J,K)
!!      QS(NX,J,K) = QS(IB,J,K)
!!      QG(1 ,J,K) = QG(IE,J,K)
!!      QG(NX,J,K) = QG(IB,J,K)
900   CONTINUE
!
!- check ---------------------------------------------------------
!
!      write(filename,14) 'hdifus-p-',myrank,'.out'
!14    format(a9,i1,a4)
!
!      open (10,file=trim(filename))
!
!      write(10,*)'------------ u----------'
!      do j=1,ny
!        write(10,11) j,u(1,j,10),u(2,j,10)
!      enddo
!
!      write(10,*)'------------ v----------'
!      do j=1,ny
!!        write(10,11) j,v(1,j,10),v(2,j,10)
!      enddo
!
!      write(10,*)'------------ t----------'
!      do j=1,ny
!        write(10,11) j,t(1,j,10),t(2,j,10)
!      enddo
!
!      write(10,*)'------------ qv----------'
!      do j=1,ny
!        write(10,11) j,qv(1,j,10),qv(2,j,10)
!      enddo
!
!11    format(1x,i5,2e30.20)
!      close (10)
!
!!      call mpi_finalize(j)
!!      stop 'hdifus'
!--------------------------------------------------------------


      RETURN
      END

#else


*     SUBROUTINE HDIFUS(U,V,T,QV,QC,QR,QI,QS,QG
      SUBROUTINE HDIFUS(U,V,T,QV
     _                 ,FRDT,FRDS,FRDU,FRDV,FRDP,TB,PLY,DXVPN,DXVPS)
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
      REAL*8 FRDT(NY,3),FRDS(NY,3),FRDU(NY,3),FRDV(NY,3),FRDP(NY,3)
C
      REAL*8 U(NX,NY,NL),V(NX,NY,NL),T(NX,NY,NL)
      REAL*8 QV(NX,NY,NL),QC(NX,NY,NL),QR(NX,NY,NL)
     _    ,QI(NX,NY,NL),QS(NX,NY,NL),QG(NX,NY,NL)
      REAL*8 TB(NX,NY,NL),PLY(NX,NY,NZ),DXVPN(NY),DXVPS(NY)
C
      REAL*8         D(NX,NY),DT(NX,NY),DS(NX,NY),DA(NX,NY),DB(NX,NY)
     &       ,      VR(NX,NY),QK(NX,NY),TK(NX,NY),VK(NX,NY),UK(NX,NY)
     &       ,      ROT(NX,NY), RLNT(NX,NY), RDLN(NX,NY), RDLT(NX,NY)
     &       ,      TW (NX,NY),QKC(NX,NY),QKR(NX,NY),QKI(NX,NY)
     &       ,      QKS(NX,NY),QKG(NX,NY)
C
      REAL*8  ZERO,HALF,ONE
      DATA       ZERO,HALF,ONE / 0.0E0 , 0.5E0 , 1.0E0  /
      REAL*8  FRDTN,FRDTS,FRDSI,FRDSJ,TI,RIJ,RI,VRI,TDI,FT1,FT2,FT3
     &       ,FS1,FS2,FS3,U0,V0,DTN,DTS,DTJ,DSI,DIJ,R0,R1,D0,D1,RT
     &       ,RN,DIJ1,DIJ2,FU1,FU2,FU3,FV1,FV2,FV3,FA1,FA2,FA3,VR0
     &       ,DS0,DT0,DA0,DB0,RLNT0,ROT0,RTA0,RSB0,VRU,VRV,RTAI,RSBI
     &       ,RTAJ,RSBJ,FB1,FB2,FB3,FB4,TD0,Q0,Q0C,Q0R,Q0I,Q0S,Q0G
     &       ,UIJK,VIJK,TIJK,QIJK,QIJKC,QIJKR,QIJKI,QIJKS,QIJKG
     &       ,VRN,VRS,TDN,TDS,QKN,QKNC,QKNR,QKNI,QKNS,QKNG
     &       ,QKSV,QKSC,QKSR,QKSI,QKSS,QKSG,DQN,DQNC,DQNR,DQNI,DQNS,DQNG
     &       ,DQS,DQSC,DQSR,DQSI,DQSS,DQSG,FAN,FAS,EAN,EAS,T0N
     &       ,T0S,Q0N,Q0NC,Q0NR,Q0NI,Q0NS,Q0NG
     &       ,Q0S0,Q0SC,Q0SR,Q0SI,Q0SS,Q0SG,FVI,FVJ,DXVPNJ,DXVPSJ
      INTEGER I,J,K,J1,JJ,I1,II
CWB
C
      FRDTN     = FRDT(1 ,1)
      FRDTS     = FRDT(NY,1)
      FRDSI     = FRDS(1 ,1)
      FRDSJ     = FRDS(1 ,2)
!
      DO 800 K  = 1 ,NL
C
C     CALCULATE DENSITY & SET 2-D FIELDS TO BE DIFFUSED
C
      DO 300 J = 1 ,NY
      DO 300 I = 1 ,NX
      UK(I,J)  = U(I,J,K)
      VK(I,J)  = V(I,J,K)
      QK(I,J)  = QV(I,J,K)
*     QKC(I,J)  = QC(I,J,K)
*     QKR(I,J)  = QR(I,J,K)
*     QKI(I,J)  = QI(I,J,K)
*     QKS(I,J)  = QS(I,J,K)
*     QKG(I,J)  = QG(I,J,K)
      QKC(I,J)  = 0.0
      QKR(I,J)  = 0.0
      QKI(I,J)  = 0.0
      QKS(I,J)  = 0.0
      QKG(I,J)  = 0.0
      TW(I,J)  = T(I,J,K)
300   CONTINUE
      DO 325 J  = JB,JE
      DO 325 I  = IB,IE
      TI        = TW(I,J)
      RIJ       = PLY(I,J,K) / (TI*RD)
      ROT(I,J)  = RIJ
      VR (I,J)  = ONE / RIJ
      TK (I,J)  = TI  - TB(I,J,K)
325   CONTINUE
      DO 350 J  = 1 ,NY,JE
      TI        = TW (IB,J)
      RI        = PLY(IB,J,K) / (TI*RD)
      VRI       = ONE / RI
      TDI       = TI  - TB(IB,J,K)
      DO 350 I  = IB,IE
      ROT(I,J)  = RI
      VR (I,J)  = VRI
      TK (I,J)  = TDI
350   CONTINUE
      DO 375 J  = 1 ,NY
      ROT(1 ,J) = ROT(IE,J)
      ROT(NX,J) = ROT(IB,J)
      VR (1 ,J) = VR (IE,J)
      VR (NX,J) = VR (IB,J)
      TK (1 ,J) = TK (IE,J)
      TK (NX,J) = TK (IB,J)
375   CONTINUE
C
C     CALCULATE DEFORMATION FIELDS
C
      DO 400 J  = JB,JE
      J1        = J - 1
      JJ        = J + 1
      FT1       = FRDT(J,1)
      FT2       = FRDT(J,2)
      FT3       = FRDT(J,3)
      FS1       = FRDS(J,1)
      FS2       = FRDS(J,2)
      FS3       = FRDS(J,3)
      DO 400 I  = IB,IE
      U0        = UK(I,J)
      V0        = VK(I,J)
      DT(I,J)   = FT1*(UK(I+1,J)-U0) - (FT2*V0 - FT3*VK(I,J1 ))
      DS(I,J)   = FS1*(V0-VK(I-1,J)) + (FS2*UK(I,JJ) - FS3*U0)
400   CONTINUE
      DTN       = ZERO
      DTS       = ZERO
      DO 425 I  = IB,IE
      DTN       = DTN + VK(I,1 )
      DTS       = DTS + VK(I,JE)
425   CONTINUE
      DTN       = DTN * FRDTN
      DTS       = DTS * FRDTS
      DO 450 I  = IB,IE
      DT(I,1 )  = DTN
      DT(I,NY)  = DTS
      DS(I,NY)  = ZERO
      DS(I,1 )  = FRDSI*(VK(I,1) - VK(I-1,1)) + FRDSJ*UK(I,JB)
450   CONTINUE
      DO 475 J  = 1 ,NY
      DT(1 ,J)  = DT(IE,J)
      DT(NX,J)  = DT(IB,J)
      DS(1 ,J)  = DS(IE,J)
      DS(NX,J)  = DS(IB,J)
475   CONTINUE
C
      DO 500 J  = 1 ,JE
      JJ        = J + 1
      DO 525 I  = IB,IE
      DTJ       = DT(I,JJ) + DT(I,J)
      DSI       = DS(I,J ) + DS(I-1,J)
      DIJ       = DTJ*DTJ  + DSI*DSI
      D(I,J)    = HALF * SQRT( DIJ )
525   CONTINUE
      D(1 ,J)   = D(IE,J)
      D(NX,J)   = D(IB,J)
500   CONTINUE
      DO I = 1 ,NX
        D(I,NY) = ZERO
      ENDDO
      DO 550 J  = JB,JE
      JJ        = J + 1
      J1        = J - 1
      DO 550 I  = IB,IE
      I1        = I - 1
      R0        = ROT(I,J)
      R1        = ROT(I,JJ)
      D0        = D(I,J)
      D1        = D(I,J1)
      RT        = R0 + R1
      RN        = R0 + ROT(I1,J)
      RLNT(I,J) = RN + R1 + ROT(I1,JJ)
      DIJ1      = HALF * (D0 + D(I1,J))
      DIJ       = HALF * (D0 + D1)
      DIJ2      = HALF * (DIJ1 + HALF*(D1 + D(I1,J1)) )
      RDLN(I,J) = RN   *  DIJ2
      RDLT(I,J) = RT   *  D0
      DA  (I,J) = DIJ
      DB  (I,J) = DIJ1
550   CONTINUE
      DO J = 1 ,NY,JE
        DO I = 1 ,NX
          RDLN(I,J) = ZERO
          DA(I,J)   = ZERO
        ENDDO
      ENDDO
      R0        = ROT (IB,1)
      DO 560 I  = IB,IE
      D0        = D(I,1)
      DB  (I,1) = (D0 + D(I-1,1 )) * HALF
      RDLT(I,1) = (R0 + ROT(I,JB)) * D0
      RLNT(I,1) =  R0 + R0 + ROT(I-1,JB) + ROT(I,JB)
560   CONTINUE
      DO I = 1 ,NX
        RDLT(I,NY) = ZERO
        RLNT(I,NY) = ZERO
        DB(I,NY)   = ZERO
      ENDDO
      DO 580 J  = 1 ,JE
      RLNT(1 ,J)= RLNT(IE,J)
      RLNT(NX,J)= RLNT(IB,J)
      RDLN(1 ,J)= RDLN(IE,J)
      RDLN(NX,J)= RDLN(IB,J)
      RDLT(1 ,J)= RDLT(IE,J)
      RDLT(NX,J)= RDLT(IB,J)
      DA  (1 ,J)= DA  (IE,J)
      DA  (NX,J)= DA  (IB,J)
      DB  (1 ,J)= DB  (IE,J)
      DB  (NX,J)= DB  (IB,J)
580   CONTINUE
C
C     UPDATE T & U & V & QV,QC,QI,QR,QS,QG DUE TO THE HORIZONTAL DIFFUSION
C
      DO 600 J  = JB,JE
      FU1       = FRDU(J,1)
      FU2       = FRDU(J,2)
      FU3       = FRDU(J,3)
      FV1       = FRDV(J,1)
      FV2       = FRDV(J,2)
      FV3       = FRDV(J,3)
      FA1       = FRDP(J,1)
      FA2       = FRDP(J,2)
      FA3       = FRDP(J,3)
      JJ        = J + 1
      J1        = J - 1
      DO 600 I  = IB,IE
      II        = I + 1
      I1        = I - 1
      VR0       = VR  (I,J)
      DS0       = DS  (I,J)
      DT0       = DT  (I,J)
      DA0       = DA  (I,J)
      DB0       = DB  (I,J)
      RLNT0     = RLNT(I,J)
      ROT0      = ROT (I,J)
      RTA0      = DT0 * DA0 * ROT0
      RSB0      = DS0 * DB0 * RLNT0
      VRU       = VR0 + VR(I1,J)
      VRV       = VR0 + VR(I ,JJ)
      RTAI      = ROT (I1,J)*DT(I1,J)*DA(I1,J)
      RSBI      = RLNT(I,J1)*DS(I,J1)*DB(I,J1)
      RTAJ      = ROT (I,JJ)*DT(I,JJ)*DA(I,JJ)
      RSBJ      = RLNT(II,J)*DS(II,J)*DB(II,J)
      FB1       = FA1 * RDLN(II,J)
      FB2       = FA1 * RDLN(I ,J)
      FB3       = FA2 * RDLT(I ,J)
      FB4       = FA3 * RDLT(I,J1)
      TD0       = TK(I,J)
      Q0        = QK(I,J)
      Q0C       = QKC(I,J)
      Q0R       = QKR(I,J)
      Q0I       = QKI(I,J)
      Q0S       = QKS(I,J)
      Q0G       = QKG(I,J)

      UIJK      = VRU * ( FU1*(RTA0-RTAI) + (FU2*RSB0 - FU3*RSBI) )
      VIJK      = VRV * ( FV1*(RSBJ-RSB0) - (FV2*RTAJ - FV3*RTA0) )
      TIJK      = VR0 * ( FB1*(TK(II,J)-TD0) - FB2*(TD0-TK(I1,J))
     &          +         FB3*(TK(I,JJ)-TD0) - FB4*(TD0-TK(I,J1)) )
      QIJK      = VR0 * ( FB1*(QK(II,J)-Q0 ) - FB2*(Q0 -QK(I1,J))
     &          +         FB3*(QK(I,JJ)-Q0 ) - FB4*(Q0 -QK(I,J1)) )
      QIJKC     = VR0 * ( FB1*(QKC(II,J)-Q0C) - FB2*(Q0C-QKC(I1,J))
     &          +         FB3*(QKC(I,JJ)-Q0C) - FB4*(Q0C-QKC(I,J1)) )
      QIJKR     = VR0 * ( FB1*(QKR(II,J)-Q0R) - FB2*(Q0R-QKR(I1,J))
     &          +         FB3*(QKR(I,JJ)-Q0R) - FB4*(Q0R-QKR(I,J1)) )
      QIJKI     = VR0 * ( FB1*(QKI(II,J)-Q0I) - FB2*(Q0I-QKI(I1,J))
     &          +         FB3*(QKI(I,JJ)-Q0I) - FB4*(Q0I-QKI(I,J1)) )
      QIJKS     = VR0 * ( FB1*(QKS(II,J)-Q0S) - FB2*(Q0S-QKS(I1,J))
     &          +         FB3*(QKS(I,JJ)-Q0S) - FB4*(Q0S-QKS(I,J1)) )
      QIJKG     = VR0 * ( FB1*(QKG(II,J)-Q0G) - FB2*(Q0G-QKG(I1,J))
     &          +         FB3*(QKG(I,JJ)-Q0G) - FB4*(Q0G-QKG(I,J1)) )
C
      U(I,J,K)  = UK(I,J) + UIJK
      V(I,J,K)  = VK(I,J) + VIJK
      T(I,J,K)  = TW(I,J) + TIJK
      QV(I,J,K)  = QK (I,J) + QIJK
      QC(I,J,K)  = QKC(I,J) + QIJKC
      QR(I,J,K)  = QKR(I,J) + QIJKR
      QI(I,J,K)  = QKI(I,J) + QIJKI
      QS(I,J,K)  = QKS(I,J) + QIJKS
      QG(I,J,K)  = QKG(I,J) + QIJKG
600   CONTINUE

      VRN       = VR(IB,1 )
      VRS       = VR(IB,NY)
      TDN       = TK(IB,1 )
      TDS       = TK(IB,NY)
      QKN       = QK(IB,1 )
      QKSV      = QK(IB,NY)
      QKNC      = QKC(IB,1 )
      QKSC      = QKC(IB,NY)
      QKNR      = QKR(IB,1 )
      QKSR      = QKR(IB,NY)
      QKNI      = QKI(IB,1 )
      QKSI      = QKI(IB,NY)
      QKNS      = QKS(IB,1 )
      QKSS      = QKS(IB,NY)
      QKNG      = QKG(IB,1 )
      QKSG      = QKG(IB,NY)
      DQN       = ZERO
      DQS       = ZERO
      DQNC      = ZERO
      DQSC      = ZERO
      DQNR      = ZERO
      DQSR      = ZERO
      DQNI      = ZERO
      DQSI      = ZERO
      DQNS      = ZERO
      DQSS      = ZERO
      DQNG      = ZERO
      DQSG      = ZERO
      DTN       = ZERO
      DTS       = ZERO
      FAN       = FRDP(1 ,1) * VRN
      FAS       = FRDP(NY,1) * VRS
      DO 625 I  = IB,IE
      EAN       = RDLT(I,1 )
      EAS       = RDLT(I,JE)
      DTN       = DTN + EAN*(TK(I,JB) - TDN)
      DTS       = DTS + EAS*(TDS - TK(I,JE))
      DQN       = DQN + EAN*(QK(I,JB) - QKN)
      DQS       = DQS + EAS*(QKSV- QK(I,JE))
      DQNC      = DQNC+ EAN*(QKC(I,JB) - QKNC)
      DQSC      = DQSC+ EAS*(QKSC - QKC(I,JE))
      DQNR      = DQNR+ EAN*(QKR(I,JB) - QKNR)
      DQSR      = DQSR+ EAS*(QKSR - QKR(I,JE))
      DQNI      = DQNI+ EAN*(QKI(I,JB) - QKNI)
      DQSI      = DQSI+ EAS*(QKSI - QKI(I,JE))
      DQNS      = DQNS+ EAN*(QKS(I,JB) - QKNS)
      DQSS      = DQSS+ EAS*(QKSS - QKS(I,JE))
      DQNG      = DQNG+ EAN*(QKG(I,JB) - QKNG)
      DQSG      = DQSG+ EAS*(QKSG - QKG(I,JE))
625   CONTINUE
      T0N       = TW(IB,1 ) + DTN*FAN
      T0S       = TW(IB,NY) + DTS*FAS
      Q0N       = QK(IB,1 ) + DQN*FAN
      Q0S0       = QK(IB,NY) + DQS*FAS
      Q0NC      = QKC(IB,1 ) + DQNC*FAN
      Q0SC      = QKC(IB,NY) + DQSC*FAS
      Q0NR      = QKR(IB,1 ) + DQNR*FAN
      Q0SR      = QKR(IB,NY) + DQSR*FAS
      Q0NI      = QKI(IB,1 ) + DQNI*FAN
      Q0SI      = QKI(IB,NY) + DQSI*FAS
      Q0NS      = QKS(IB,1 ) + DQNS*FAN
      Q0SS      = QKS(IB,NY) + DQSS*FAS
      Q0NG      = QKG(IB,1 ) + DQNG*FAN
      Q0SG      = QKG(IB,NY) + DQSG*FAS
      FVI       = FRDV(1,1)
      FVJ       = FRDV(1,2)
      DO 650 I  = IB,IE
      II        = I + 1
      T (I,1 ,K) = T0N
      T (I,NY,K) = T0S
      QV(I,1 ,K) = Q0N
      QV(I,NY,K) = Q0S0
      QC(I,1 ,K) = Q0NC
      QC(I,NY,K) = Q0SC
      QR(I,1 ,K) = Q0NR
      QR(I,NY,K) = Q0SR
      QI(I,1 ,K) = Q0NI
      QI(I,NY,K) = Q0SI
      QS(I,1 ,K) = Q0NS
      QS(I,NY,K) = Q0SS
      QG(I,1 ,K) = Q0NG
      QG(I,NY,K) = Q0SG
      U (I,1 ,K) = ZERO
      U (I,NY,K) = ZERO
      V (I,NY,K) = ZERO
      VRV       = VRN + VR(I,JB)
      RSBI      = RLNT(II,1)*DS(II,1)*DB(II,1)
      RSB0      = RLNT(I ,1)*DS(I ,1)*DB(I ,1)
      RTA0      = ROT (I,JB)*DT(I,JB)*DA(I,JB)
      VIJK      = VRV * ( FVI*(RSBI - RSB0) - FVJ*RTA0 )
      V(I,1 ,K) = VK(I,1) + VIJK
650   CONTINUE
C
C     ENERGY CONSERVATION DUE TO FRICTIONS
C
      DO 700 J  = JB,JE
      DO 700 I  = IB,IE
      DA(I,J)   = UK(I,J) * (U(I,J,K) - UK(I,J))
700   CONTINUE
      DO 710 J  = JB,JE
*     DA(IE,J)  = DA(1,J)
CWB
      DA(1 ,J)  = DA(IE,J)
      DA(NX,J)  = DA(IB,J)
CWB
710   CONTINUE
      DO 720 J  = 1 ,JE
      DO 720 I  = IB,IE
      DB(I,J)   = VK(I,J) * (V(I,J,K) - VK(I,J))
720   CONTINUE
      DO 750 J  = JB,JE
      JJ        = J - 1
      DXVPNJ    = DXVPN(J)
      DXVPSJ    = DXVPS(J)
      DO 750 I  = IB,IE
      D (I,J)   = HALF * (DA(I,J) + DA(I+1,J))
     &          + DXVPNJ*DB(I,JJ) + DXVPSJ*DB(I,J)
      TK(I,J)   = T(I,J,K)
750   CONTINUE
      DO 800 J  = JB,JE
      DO 800 I  = IB,IE
      T(I,J,K)  = TK(I,J) - D(I,J)/CPD
800   CONTINUE
C
C     SPHERICAL BOUNDARY CONDITIONS
C
      DO 900 K  = 1 ,NL
      DO 900 J  = 1 ,NY
      U(1 ,J,K) = U(IE,J,K)
      U(NX,J,K) = U(IB,J,K)
      V(1 ,J,K) = V(IE,J,K)
      V(NX,J,K) = V(IB,J,K)
      T(1 ,J,K) = T(IE,J,K)
      T(NX,J,K) = T(IB,J,K)
      QV(1 ,J,K) = QV(IE,J,K)
      QV(NX,J,K) = QV(IB,J,K)
      QC(1 ,J,K) = QC(IE,J,K)
      QC(NX,J,K) = QC(IB,J,K)
      QR(1 ,J,K) = QR(IE,J,K)
      QR(NX,J,K) = QR(IB,J,K)
      QI(1 ,J,K) = QI(IE,J,K)
      QI(NX,J,K) = QI(IB,J,K)
      QS(1 ,J,K) = QS(IE,J,K)
      QS(NX,J,K) = QS(IB,J,K)
      QG(1 ,J,K) = QG(IE,J,K)
      QG(NX,J,K) = QG(IB,J,K)
900   CONTINUE
!
!- check ---------------------------------------------------------
!
!      open (10,file='hdifus-s.out')
!
!      write(10,*)'------------ u----------'
!      do j=1,ny
!        write(10,11) j,u(1,j,10),u(2,j,10)
!      enddo
!
!      write(10,*)'------------ v----------'
!      do j=1,ny
!        write(10,11) j,v(1,j,10),v(2,j,10)
!      enddo
!
!      write(10,*)'------------ t----------'
!      do j=1,ny
!        write(10,11) j,t(1,j,10),t(2,j,10)
!      enddo
!
!      write(10,*)'------------ qv----------'
!      do j=1,ny
!        write(10,11) j,qv(1,j,10),qv(2,j,10)
!      enddo
!
!11    format(1x,i5,2e30.20)
!      close (10)
!
!!      stop 'hdifus'
!--------------------------------------------------------------

      RETURN
      END

#endif
