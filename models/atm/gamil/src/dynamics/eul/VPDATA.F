#include <misc.h>
#include <params.h>

!!(wh 2003.11.12)
!!-------------------

      SUBROUTINE VPDATA(Q,W,DTDSG,EP,NONOS,ISOR,IORD)
C
C     PERFORM 1-D ADVECTION IN THE VERTICAL DIRECTION
C             WITH A NON-UNIFORM SIGMA C-GRID MESH
C                             ________W________  K'
C
C     GET WORKSHOPS INDEPENDENT OF HORIZONTAL MESH
!
      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

#if (defined SPMD)
#include <commpi.h>
      character*50 filename
#endif
C
      REAL*8  Q(NX,NY,NL),W(NX,NY,NL)
      REAL*8  HW(NL   ),HS(NL   ),Z (NL   )
     &       ,VW(NZ   ),ZM(NL   ),ZN(NL   ),FW(NZ   )
     &       ,UW(NZ   ),BUW(  NL),BDW(  NL)
C
	REAL*8  EP,DTDSG(NL)
      REAL*8  ZERO,HALF,ONE
      DATA ZERO,HALF,ONE / 0.0D0,0.5D0,1.0D0 /
      INTEGER I,J,K,IO,IS,IT
      INTEGER NONOS,IORD,ISOR
C     =========================================================
      REAL*8  A,VDYF,R,X1,X2,YY,VC31,PP,Y,PN
      VDYF(R,A,X1,X2)= (ABS(A)-A**2/R)*(X2-X1)/(X2+X1+EP)
      VC31(R,A,YY)   =-(1.0-3.0*ABS(A)/R+2.0*(A/R)**2)*A*YY/3.0
      PP(Y)          = MAX(0.0,Y)
      PN(Y)          =-MIN(0.0,Y)
C     =========================================================
      DO 900 K = 1 ,NL
      HS(K)    = ONE  / DTDSG(K)
900   CONTINUE
      DO 905 K = 2 ,NL
      HW(K)    = HALF * ( HS(K)+HS(K-1) )
905   CONTINUE

!- check ---------------------------------------------------------
!
!#if (defined SPMD)
!      write(filename,14) 'qpdata-1-p-',myrank,'.out'
!14    format(a11,i1,a4)
!
!      open (10,file=trim(filename))
!#else
!      open (10,file='qpdata-1-s.out')
!#endif
!
!      write(10,*) '----------------- q -------------------'
!      do j=1,ny
!        write(10,11) j,q(1,j,10),q(2,j,10)
!      enddo
!
!      write(10,*) '----------------- w-------------------'
!      do j=1,ny
!        write(10,11) j,w(1,j,10),w(2,j,10)
!      enddo
!
!      write(10,*) '----------------- dtdsg -------------------'
!      write(10,12) (j,dtdsg(j),j=1,nl)
!
!      write(10,*) '----------------- hs-------------------'
!      write(10,12) (j,hs(j),j=1,nl)
!
!      write(10,*) '----------------- hw-------------------'
!      write(10,12) (j,hs(j),j=2,nl)
!
!11    format(1x,i5,2e30.20)
!12    format(1x,i5,e30.20)
!      close (10)
!
!#if (defined SPMD)
!!      call mpi_finalize(j)
!#endif
!!      stop'qpdata'
!--------------------------------------------------------------

C
C     START HORIZONTAL GRID DO LOOP
C     KEEP IN MIND THAT Z ARE THE SAME FOR ALL I AT POLE
C$DOACROSS LOCAL(J,IS,IT,I,K,Z,UW,VW,ZM,ZN,FW,BUW,BDW,IO)

!      write(999,*) 'vpdata -----------'

#if (defined SPMD)
      do 1000 j = 2 ,ny-1
#else
      DO 1000 J = 1 ,NY
#endif

!      write(999,*) 'j=',j

      FW(01)   = ZERO
      FW(NZ)   = ZERO
#if (defined SPMD)
      IF( (myrank.eq.nprocs-1).and.(J.EQ.jpole) ) THEN
        IS = NX
        IT = NX
      ELSE IF( (myrank.eq.0).and.(J.EQ.jpole) ) THEN
        IS = 1
        IT = 1
      ELSE
        IS = 1
        IT = NX
      ENDIF
#else
      IF(J.EQ.1) THEN
        IS = NX
        IT = NX
      ELSE IF(J.EQ.NY) THEN
        IS = 1
        IT = 1
      ELSE
        IS = 1
        IT = NX
      ENDIF
#endif

      DO 990 I = IS,IT
      DO 910 K = 1 ,NL
      Z (K)    = Q(I,J,K)
      UW(K)    = W(I,J,K)
      VW(K)    = UW(K)
910   CONTINUE
C
C     PREPARE FOR NON-OSSCILATORY OPTION
      IF( NONOS.EQ.1 ) THEN
      DO 920 K = 2 ,NM
      ZM(K)    = MAX( Z(K-1),Z(K),Z(K+1) )
      ZN(K)    = MIN( Z(K-1),Z(K),Z(K+1) )
920   CONTINUE
      ZM(1)    = MAX( Z(1),Z(2) )
      ZN(1)    = MIN( Z(1),Z(2) )
      ZM(NL)   = MAX( Z(NM),Z(NL) )
      ZN(NL)   = MIN( Z(NM),Z(NL) )
      ENDIF
C
      DO 980 IO= 1 ,IORD
C     ++++++++++++++++++++++++++++++++
C     PREDICTOR STEP : UPSTREAM SCHEME
C     ++++++++++++++++++++++++++++++++
      DO 930 K = 2 ,NL
Cb    FW(K)    = DONOR( Z(K-1),Z(K),VW(K) )
!      if ( k.eq.24) write(999,*) i,io,k,'vw(k)=',vw(k)

      IF(VW(K).GE.0.0D0) THEN
        FW(K) = Z(K-1)*VW(K)
      ELSE
        FW(K) = Z(K)*VW(K)
      ENDIF
930   CONTINUE
      DO 935 K = 1 ,NL
      Z(K)     = Z(K)-(FW(K+1)-FW(K))*DTDSG(K)
935   CONTINUE
C
      IF( IO.EQ.IORD ) GOTO 980
C     ++++++++++++++++++++++++++++++++++++++
C     CORRECTOR STEP : ANTI-DIFFUSION SCHEME
C     ++++++++++++++++++++++++++++++++++++++
      DO 940 K = 2 ,NL
      FW(K)    = VW(K)
940   CONTINUE
C
C     CALCULATE THE  PSEUDO VELOCITIES
      DO 945 K = 2 ,NL
      VW(K)    = VDYF( HW(K),FW(K),Z(K-1),Z(K) )
945   CONTINUE
C     ADD THE THIRD ORDER CORRECTION IF REQUESTED
      IF( ISOR.EQ.3 ) THEN
      DO 950 K = 3 ,NM
      VW(K)    = VW(K) + VC31( HW(K),FW(K)
     &         ,    (Z(K-2)+Z(K+1)-Z(K-1)-Z(K))
     &         / (EP+Z(K-2)+Z(K+1)+Z(K-1)+Z(K)) )
950   CONTINUE
C     ASSUME CONSTANT Z ABOVE K=1 & BELOW K=NL
      VW(2)    = VW(2) + VC31( HW(2),FW(2)
     &         ,         (Z(3)-Z(2))
     &         / (EP+Z(1)+Z(3)+Z(1)+Z(2)) )
      VW(NL)    = VW(NL) + VC31( HW(NL),FW(NL)
     &         ,    (Z(NL-2)-Z(NL-1))
     &         / (EP+Z(NL-2)+Z(NL)+Z(NL-1)+Z(NL)) )
      ENDIF
C
      DO 955 K = 2 ,NL
      VW(K)    = SIGN(ONE,VW(K))*MIN(ABS(UW(K)),ABS(VW(K)))
955   CONTINUE
C
C     PERFORM THE NON-OSSCILATORY OPTION
      IF( NONOS.EQ.1 ) THEN
      DO 960 K = 2 ,NM
      ZM(K)    = MAX( Z(K-1),Z(K),Z(K+1),ZM(K) )
      ZN(K)    = MIN( Z(K-1),Z(K),Z(K+1),ZN(K) )
960   CONTINUE
      ZM(1)    = MAX( Z(1),Z(2),ZM(1) )
      ZN(1)    = MIN( Z(1),Z(2),ZN(1) )
      ZM(NL)    = MAX( Z(NM),Z(NL),ZM(NL) )
      ZN(NL)    = MIN( Z(NM),Z(NL),ZN(NL) )
C
      DO 965 K = 2 ,NL
Cb    FW(K)    = DONOR( Z(K-1),Z(K),VW(K) )
      IF(VW(K).GE.0.0D0) THEN
        FW(K) = Z(K-1)*VW(K)
      ELSE
        FW(K) = Z(K)*VW(K)
      ENDIF
965   CONTINUE
      DO 970 K = 1 ,NL
      BUW(K)   = (ZM(K)-Z(K))*HS(K) / (PN(FW(K+1))+PP(FW(K))+EP)
      BDW(K)   = (Z(K)-ZN(K))*HS(K) / (PP(FW(K+1))+PN(FW(K))+EP)
970   CONTINUE
      DO 975 K = 2 ,NL
      VW(K)    = PP( VW(K) ) * MIN(ONE,BDW(K-1),BUW(K))
     &         - PN( VW(K) ) * MIN(ONE,BUW(K-1),BDW(K))
975   CONTINUE
      ENDIF
980   CONTINUE
C
C     UPDATE THE PREDICTED FIELD
      DO 990 K = 1 ,NL
      Q(I,J,K) = Z(K)
990   CONTINUE

#if (defined SPMD)
      IF( (myrank.eq.nprocs-1).and.(J.EQ.jpole) )THEN
        DO K = 1 ,NL
          DO I = 1 ,NX
            Q(I,J,K)= Q(NX,J,K)
          ENDDO
        ENDDO
      ELSE IF( (myrank.eq.nprocs-1).and.(J.EQ.jpole) .OR.
     _         (myrank.eq.0       ).and.(J.EQ.jpole)       ) THEN
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
      ELSE IF(J.EQ.1.OR.J.EQ.NY) THEN
        DO K = 1 ,NL
          DO I = 1 ,NX
            Q(I,J,K)= Q(1,J,K)
          ENDDO
        ENDDO
      ENDIF
#endif

1000  CONTINUE

!- check ---------------------------------------------------------
!#if (defined SPMD)
!      write(filename,18) 'qpdata-2-p-',myrank,'.out'
!18    format(a11,i1,a4)
!
!      open (10,file=trim(filename))
!#else
!      open (10,file='qpdata-2-s.out')
!#endif
!
!      do j=1,ny
!        write(10,19) j,q(1,j,10),q(2,j,10)
!      enddo
!
!19    format(1x,i5,2e30.20)
!      close (10)
!
!#if (defined SPMD)
!      call mpi_finalize(j)
!#endif
!      stop'qpdata-2'
!!--------------------------------------------------------------

      RETURN
      END
