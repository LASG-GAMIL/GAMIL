#include <misc.h>
#include <params.h>


!!      SUBROUTINE DYFRAM2(NSEQ,DSGHL,IP
!!     _                 ,DSNP,DSSP,DTDLN,DTDLT,GC,DTDSG
!!     _                 ,PMTOP,GHS,SINU,SINV,WTGU,WTGV,SIG,SIGL,DSIG
!!     _                 ,CBB,TBB,DTDY,ITIME,UT,VT,TT,QD
!!     _                 ,PS,U,V,WS,WPA,T,Q,GHI,PLY,TB)

        subroutine dyfram2( nseq,dtdy,itime
     _                   ,u,v,t,q,ws,ps,wpa,ghs,ghi,ply,tb
     _                   ,ut,vt,tt,qd
     _                   ,nonos,iord,isor,ep,ip,dsnp,dssp,dtdln,dtdlt,gc,dtdsg,dsghl
     _                   ,pmtop,sig,sigl,dsig
     _                   ,tbb,hbb,cbb,dcbb,psb,tsb
     _                   ,dy,wtgu,wtgv
     _                   ,dx,sinu,sinv,oux,ouy,ovx,ovy,ff,cur
     _                   ,mm1,mp1,mm2,mp2,mm3,mp3,mdj )

!     **************************************************
!     **************************************************
!
!     MAIN ROUTINE  OF THE DYNAMIC FRAME CALCULATION

#if ( defined SPMD )
      use mpishorthand, only: mpicom
#endif


      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

#if (defined SPMD)
#include <mpif.h>
#include <commpi.h>

      integer   isend, irecv
      integer   istatus(mpi_status_size)
      character*50 filename
#endif

!
      real*8 dtdy
      INTEGER*4 ITIME

      real*8   pmtop,sig(nz),dsig(nl),sigl(nl)
      real*8   dy,wtgu(ny),wtgv(ny)
      real*8   dx,sinu(ny),sinv(ny),oux(ny),ouy(ny),ovx(ny),ovy(ny),ff(ny),cur(ny)

      integer  MM1(NX,NY)     !    MM1=i-MDJ      , if MM1<2, MM1=NX-2+MM1
      integer  MP1(NX,NY)     !    MP1=i+MDJ      , if MP1>NX-1, MP1=MP1-NX+2
      integer  MM2(NX,NY)     !    MM2=i-(MDJ-1)/2, if MM2<2, MM2=NX-2+MM2
      integer  MP2(NX,NY)     !    MP2=i+(MDJ+1)/2, if MP2>NX-1, MP2=MP2-NX+2
      integer  MM3(NX,NY)     !    MM3=i-(MDJ+1)/2, if MM3<2, MM3=NX-2+MM3
      integer  MP3(NX,NY)     !    MP3=i+(MDJ-1)/2, if MP3>NX-1, MP3=MP3-NX+2
      integer  MDJ(NY)        !    leaping span of the difference

      REAL*8   CBB(NA),TBB(NA),hbb(na),dcbb(na),psb(nx,ny),tsb(nx,ny)

      INTEGER NONOS,IORD,ISOR,IP(NX)
!
      REAL*8 EP,DSNP,DSSP,DTDLN(NY),DTDLT(NY),GC(NY),DTDSG(NL)
      REAL*8  DSGHL(NL)
C
      REAL*8 U0(NX,NY,NL),V0(NX,NY,NL),W0(NX,NY,NL)
      REAL*8 CB(NX,NY,NL),TB(NX,NY,NL)
      REAL*8 WS(NX,NY,NZ),PLY(NX,NY,NZ),GHI(NX,NY,NZ)
C
      REAL*8 GHS(NX,NY),DP(NX,NY)
C
      REAL*8 PS(NX,NY),U(NX,NY,NL),V(NX,NY,NL),T(NX,NY,NL)
      REAL*8 QT(NX,NY,NL),Q(NX,NY,NL)
C
C     SOURCE & SINK  OF THE MODEL PREDICTED VARIABLES   DYNAMICS/PHYS
      REAL*8 SU(NX,NY,NL),SV(NX,NY,NL),ST(NX,NY,NL)
      REAL*8 UT(NX,NY,NL),VT(NX,NY,NL),TT(NX,NY,NL),QD(NX,NY,NL)
C
      REAL*8 WPA(NX,NY,NL),FAC(NX,NY,NZ),FBC(NX,NY,NL)
      REAL*8 P(NX,NY)
!      COMMON /QTT/QT,QTC,QTR,QTI,QTS,QTG
      COMMON /QTT/QT
C     _______________________________________________________________
C
      REAL*8  WK5,WK6,PX,PY
      INTEGER NSEQ,I,J,K,NCYC,KPP
      integer begj, endj
!
!     START   THE DYNAMICAL INTEGRATING CYCLE
!

#if (defined SPMD)
      begj = 2
      endj = ny-1
#else
      begj = 1
      endj = ny
#endif

      DO J=begj,endj
      DO I=2,NX-1
	 P(I,J)=SQRT(PS(I,J))
      END DO
      P(1 ,J)=P(NX-1,J)
      P(NX,J)=P(2,J)
      END DO

#if (defined SPMD)
      call mpi_isend( p(1,2), nx, mpi_double_precision, itop, 1, mpicom, isend, ierr)
      call mpi_irecv( p(1,ny),nx, mpi_double_precision, ibot, 1, mpicom, irecv, ierr)
#endif


      DO J=begj,endj
      DO K=1,NL
      DO I=2,NX-1
         QT(I,J,K)=Q(I,J,K)*PS(I,J)
         WPA(I,J,K)=0.0
      END DO
      QT(1,J,K)=QT(NX-1,J,K)
      QT(NX,J,K)=QT(2,J,K)
      WPA(1,J,K)=0.0
      WPA(NX,J,K)=0.0
      END DO
      END DO
!
!     TRANSFORM THE FORCING TERMS OF WIND, TEMPERATURE AND MOISTURE
!     TO THE CONTRIBUTIONS TO THE TENDENCIES OF THE BAISC MODEL
!     PREDICTION VARIABLES
!
#if (defined SPMD)
      call mpi_wait(isend,istatus,ierr)
      call mpi_wait(irecv,istatus,ierr)
#endif

       DO J=begj,endj

         DO I=2,NX-1
!
! 	 CALCULATING THE AVERAGE VALUE OF P AT V-GRID.
!
#if (!defined SPMD)
         IF(J.LT.NY) THEN
 	    PY=0.5D0*(P(I,J)+P(I,J+1))
         ELSE
 	    PY=0.0D0
 	 ENDIF
#else
         if ( (myrank.eq.0).and.(j.eq.endj) ) then
 	     PY=0.0D0
         else
 	     PY=0.5D0*(P(I,J)+P(I,J+1))
         endif
#endif
!
!	 CALCULATING THE AVERAGE VALUE OF P AT U-GRID.
!
         PX=0.5D0*(P(I,J)+P(I-1,J))
!
         DO K=1,NL
!
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
!
      DO K = 1, NL
      DO J = begj,endj
	 SU(1,J,K)=SU(NX-1,J,K)
	 SU(NX,J,K)=SU(2,J,K)
	 SV(1,J,K)=SV(NX-1,J,K)
	 SV(NX,J,K)=SV(2,J,K)
	 ST(1,J,K)=ST(NX-1,J,K)
	 ST(NX,J,K)=ST(2,J,K)
      END DO
      END DO
!
!     THE END OF THE TRANSFORMATION
!
      DO 500 NCYC = 1 ,NSEQ
!
!     PREDICT DRY-ADIABATIC SYSTEM
!     ___________________________
!
      DO K=1,NL
	DO J=begj,endj
	   DO I=1,NX
	      U0(I,J,K)=U(I,J,K)
	      V0(I,J,K)=V(I,J,K)
	      W0(I,J,K)=WS(I,J,K)
              FAC(I,J,K)=1.0
              FBC(I,J,K)=1.0
	   END DO
	END DO
      END DO

        do j=begj,endj
          do i=1,nx
              FAC(I,J,NZ)=1.0
           enddo
        enddo

!!-----------------------------------------------------------------

!!    CALL DYNAMICS(U,V,WS,W0,PS,T,GHI,DP,SU,SV,ST,FAC,FBC     !!(wh 2003.11.04)
!!    _                  ,GHS          ,DTDY,ITIME,NCYC,contn) !!

       !write(*,*) "sxj--before dynamics"

      call dynamics( dtdy,itime,ncyc,
     _               u,v,ws,w0,ps,t,ghi,ghs,dp,su,sv,st,fac,fbc,
     _               pmtop,sigl,dsig,
     _               tbb,hbb,cbb,dcbb,psb,tsb,
     _               dy,wtgu,wtgv,
     _               dx,sinu,sinv,oux,ouy,ovx,ovy,ff,cur,
     _               mm1,mp1,mm2,mp2,mm3,mp3,mdj )

        !write(*,*) "sxj--after dynamics"

      ITIME=ITIME+INT(DTDY+0.1)

!
!     PREDICT WATER VAPOR MIXING RATIO
!     BY POSITIVE DEFINITE ADVECTION TRANSPORT ALGORITHM
!     ____________________________
!

!        CALL QPDATA2(QT,U0,V0,W0,DSGHL
!     _           ,U,V,WS,NONOS,IORD,ISOR,IP,EP
!     _           ,DSNP,DSSP,DTDLN,DTDLT,GC,DTDSG)

        CALL QPDATA1(QT,U0,V0,W0,DSGHL
     _           ,U,V,WS
     _           ,DSNP,DSSP,DTDLN,DTDLT,GC,DTDSG)


!
!     GET THE TIME AVERAGED P-SURFACE VERTICAL VELOCITY
!
      DO K = 1 ,NL
      DO J = begj,endj
         DO I = IB ,IE
            WPA(I,J,K) = WPA(I,J,K)+(0.5*(WS(I,J,K+1)+WS(I,J,K))
     +      *PS(I,J)+DP(I,J)*SIGL(K))/real(NSEQ)
	 END DO
      END DO
      END DO
!
 500  CONTINUE


      DO K = 1 ,NL
      DO J = begj,endj
         WPA(01,J,K) = WPA(IE,J,K)
         WPA(NX,J,K) = WPA(IB,J,K)
      END DO
      END DO
!
!     --------------------
!
      DO J = begj,endj
      DO I=1,NX
      DO K=1,NL
         Q(I,J,K)=QT(I,J,K)/PS(I,J)
         PLY(I,J,K)=PS(I,J)*SIGL(K)+PMTOP
         WK5=PLY(I,J,K)/DPALIB
         KPP=INT(WK5)
         WK6=WK5-DFLOAT(KPP)
         TB(I,J,K)=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
      END DO
      PLY(I,J,NZ)=PS(I,J)+PMTOP
      END DO
      END DO


!     ---------------------------

      RETURN
      END
