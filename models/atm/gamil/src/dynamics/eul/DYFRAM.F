#include <misc.h>
#include <params.h>

      subroutine dyfram2( nseq,dtdy,itime
     _                   ,u,v,t,q,ws,ps,wpa,ghs,ghi,ply,tb
     _                   ,ut,vt,tt,qd
     _                   ,pmtop,sig,sigl,dsig
     _                   ,tbb,hbb,cbb,dcbb,psb,tsb
     _                   ,dy,wtgu,wtgv
     _                   ,dx,sinu,sinv,oux,ouy,ovx,ovy,ff,cur
     _                   ,mm1,mp1,mm2,mp2,mm3,mp3,mdj )

!     MAIN ROUTINE  OF THE DYNAMIC FRAME CALCULATION

      use pmgrid, only: beglatexdyn, endlatexdyn, plat
      use mpi_gamil
      use comfm1, only: P, QT, SU, SV, ST

      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

      real*8 dtdy
      INTEGER*4 ITIME

      real*8   pmtop,sig(nz),dsig(nl),sigl(nl)
      real*8   dy,wtgu(beglatexdyn:endlatexdyn),wtgv(beglatexdyn:endlatexdyn)
      real*8   dx,sinu(beglatexdyn:endlatexdyn),sinv(beglatexdyn:endlatexdyn)
      real*8   oux(beglatexdyn:endlatexdyn),ouy(beglatexdyn:endlatexdyn)
      real*8   ovx(beglatexdyn:endlatexdyn),ovy(beglatexdyn:endlatexdyn)
      real*8   ff(beglatexdyn:endlatexdyn),cur(beglatexdyn:endlatexdyn)

      integer  MM1(beglonex:endlonex,beglatexdyn:endlatexdyn)     !    MM1=i-MDJ      , if MM1<2, MM1=NX-2+MM1
      integer  MP1(beglonex:endlonex,beglatexdyn:endlatexdyn)     !    MP1=i+MDJ      , if MP1>NX-1, MP1=MP1-NX+2
      integer  MM2(beglonex:endlonex,beglatexdyn:endlatexdyn)     !    MM2=i-(MDJ-1)/2, if MM2<2, MM2=NX-2+MM2
      integer  MP2(beglonex:endlonex,beglatexdyn:endlatexdyn)     !    MP2=i+(MDJ+1)/2, if MP2>NX-1, MP2=MP2-NX+2
      integer  MM3(beglonex:endlonex,beglatexdyn:endlatexdyn)     !    MM3=i-(MDJ+1)/2, if MM3<2, MM3=NX-2+MM3
      integer  MP3(beglonex:endlonex,beglatexdyn:endlatexdyn)     !    MP3=i+(MDJ-1)/2, if MP3>NX-1, MP3=MP3-NX+2
      integer  MDJ(beglatexdyn:endlatexdyn)        !    leaping span of the difference

      REAL*8   CBB(NA),TBB(NA),hbb(na),dcbb(na),psb(beglonex:endlonex,beglatexdyn:endlatexdyn)
      REAL*8   tsb(beglonex:endlonex,beglatexdyn:endlatexdyn)

      REAL*8 U0(beglonex:endlonex,beglatexdyn:endlatexdyn,NL),V0(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 W0(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 CB(beglonex:endlonex,beglatexdyn:endlatexdyn,NL),TB(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 WS(beglonex:endlonex,beglatexdyn:endlatexdyn,NZ),PLY(beglonex:endlonex,beglatexdyn:endlatexdyn,NZ)
      REAL*8 GHI(beglonex:endlonex,beglatexdyn:endlatexdyn,NZ)
C
      REAL*8 GHS(beglonex:endlonex,beglatexdyn:endlatexdyn),DP(beglonex:endlonex,beglatexdyn:endlatexdyn)
C
      REAL*8 PS(beglonex:endlonex,beglatexdyn:endlatexdyn),U(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 V(beglonex:endlonex,beglatexdyn:endlatexdyn,NL),T(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 Q(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
C
C     SOURCE & SINK  OF THE MODEL PREDICTED VARIABLES   DYNAMICS/PHYS
      REAL*8 UT(beglonex:endlonex,beglatexdyn:endlatexdyn,NL),VT(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 TT(beglonex:endlonex,beglatexdyn:endlatexdyn,NL),QD(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
C
      REAL*8 WPA(beglonex:endlonex,beglatexdyn:endlatexdyn,NL),FAC(beglonex:endlonex,beglatexdyn:endlatexdyn,NZ)
      REAL*8 FBC(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
C     _______________________________________________________________
C
      REAL*8  WK5,WK6,PX,PY
      INTEGER NSEQ,I,J,K,NCYC,KPP
      integer begj, endj
!
!     START   THE DYNAMICAL INTEGRATING CYCLE
!
!$OMP PARALLEL DO PRIVATE (I,J)
      DO J = JBEG0, JEND0
         DO I = IBEG1, IEND1
            P(I,J) = SQRT(PS(I,J))
         END DO
      END DO

      call gamil_arrays_comm(COMM_TO_RIGHT,1,P(:,beglatexdyn))
      call gamil_arrays_comm(COMM_TO_TOP,1,P(:,beglatexdyn))

!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K=1,NL
      DO J=jbeg0,jend0
      DO I=beglonex,endlonex
      IF (I .eq. 1 .or. I .eq. NX) THEN
         WPA(I,J,K)=0.0
      ELSE
         QT(I,J,K)=Q(I,J,K)*PS(I,J)
         WPA(I,J,K)=0.0
      ENDIF
      END DO
      END DO
      END DO

!$OMP PARALLEL DO PRIVATE (I,J,K,PY,PX,WK5,KPP,WK6,cbt)
       DO K = 1, NL
         DO J = jbeg0, jend0

           DO I = ibeg1, iend1
!
! 	 CALCULATING THE AVERAGE VALUE OF pp AT V-GRID.
!
            IF (J .LT. plat) THEN
              PY = 0.5D0*(P(I,J)+P(I,J+1))
            ELSE
              PY = 0.0D0
            END IF
!
!	 CALCULATING THE AVERAGE VALUE OF pp AT U-GRID.
!
            PX=0.5D0*(P(I,J)+P(I-1,J))

            SU(I,J,K)=PX*UT(I,J,K)
            SV(I,J,K)=PY*VT(I,J,K)
!
!           CALCULATING STT=ST*P*RD/CB.
!
            WK5=(PS(I,J)*SIGL(K)+PMTOP)/DPALIB
            KPP=INT(WK5)
            WK6=WK5-DFLOAT(KPP)
            CB(I,J,K)=(1.D0-WK6)*CBB(KPP)+WK6*CBB(KPP+1)
            ST(I,J,K)=TT(I,J,K)*P(I,J)*RD/CB(I,J,K)
          END DO
        END DO
      END DO

      call gamil_arrays_comm(COMM_TO_LEFT,  1, QT(:,beglatexdyn,1), SU(:,beglatexdyn,1), SV(:,beglatexdyn,1), ST(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT, 1, QT(:,beglatexdyn,1), SU(:,beglatexdyn,1), SV(:,beglatexdyn,1), ST(:,beglatexdyn,1))
!
!     THE END OF THE TRANSFORMATION
!
      DO NCYC = 1, NSEQ
!
!     PREDICT DRY-ADIABATIC SYSTEM
!     ___________________________
!
!$OMP PARALLEL DO PRIVATE (I,J,K)
         DO K = 1, NL
            DO J = JBEG0, JEND0
               DO I = BEGLONEX, ENDLONEX
                  U0(I,J,K) = U(I,J,K)
                  V0(I,J,K) = V(I,J,K)
                  W0(I,J,K) = WS(I,J,K)
                  FAC(I,J,K) = 1.0
                  FBC(I,J,K) = 1.0
               END DO
            END DO
         END DO

!$OMP PARALLEL DO PRIVATE (I,J)
         DO J = JBEG0, JEND0
           DO I = BEGLONEX, ENDLONEX
              FAC(I,J,NZ) = 1.0
           END DO
         END DO

      call dynamics( dtdy,itime,ncyc,
     _               u,v,ws,w0,ps,t,ghi,ghs,dp,su,sv,st,fac,fbc,
     _               pmtop,sigl,dsig,
     _               tbb,hbb,cbb,dcbb,psb,tsb,
     _               dy,wtgu,wtgv,
     _               dx,sinu,sinv,oux,ouy,ovx,ovy,ff,cur,
     _               mm1,mp1,mm2,mp2,mm3,mp3,mdj )

      ITIME=ITIME+INT(DTDY+0.1)
!
!     PREDICT WATER VAPOR MIXING RATIO
!     BY POSITIVE DEFINITE ADVECTION TRANSPORT ALGORITHM
!
      CALL QPDATA(QT, U0, V0, W0, U, V, WS, DTDY, DX, DY,
     $            SINU, SINV, WTGU, WTGV, DSIG)
!
!     GET THE TIME AVERAGED P-SURFACE VERTICAL VELOCITY
!
!$OMP PARALLEL DO PRIVATE (I,J,K)
         DO K = 1, NL
            DO J = JBEG0, JEND0
               DO I = IBEG1, IEND1
                  WPA(I,J,K) = WPA(I,J,K)+(0.5*(WS(I,J,K+1)+WS(I,J,K))
     _                 *PS(I,J)+DP(I,J)*SIGL(K))/REAL(NSEQ)
               END DO
            END DO
         END DO

      END DO

      call gamil_arrays_comm(COMM_TO_LEFT,1,WPA(:,beglatexdyn,1)) 
      call gamil_arrays_comm(COMM_TO_RIGHT,1,WPA(:,beglatexdyn,1)) 

!     --------------------
!
!$OMP PARALLEL DO PRIVATE (I,J,K,WK5,KPP,WK6)
      DO J = jbeg0,jend0
      DO K=1,NL
      DO I=beglonex,endlonex
         Q(I,J,K)=QT(I,J,K)/PS(I,J)
         PLY(I,J,K)=PS(I,J)*SIGL(K)+PMTOP
         WK5=PLY(I,J,K)/DPALIB
         KPP=INT(WK5)
         WK6=WK5-DFLOAT(KPP)
         TB(I,J,K)=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
      END DO
      END DO
      DO I=beglonex,endlonex
      PLY(I,J,NZ)=PS(I,J)+PMTOP
      END DO
      END DO

!     ---------------------------

      RETURN
      END
