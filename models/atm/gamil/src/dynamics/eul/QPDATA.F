#include <misc.h>
#include <params.h>

!!(2003.11.11)
!!-------------------

#if 0
       SUBROUTINE QPDATA2(QT,U0,V0,W0,DSGHL
     _                 ,U,V,WS,NONOS,IORD,ISOR,IP,EP
     _                 ,DSNP,DSSP,DTDLN,DTDLT,GC,DTDSG)
!     ********************************
!     ********************************
!
!     PREDICT POSITIVE DEFINITE FIELD Q        DUE TO 3-D ADVECTION
!           1)  BY USING THE SCHEME     PROPOSED BY R.C.Yu
!           2)  BY USING THE SCHEME     PROPOSED BY P.K.Smolarkiewicz
!
      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>
#include <commpi.h>

C
      REAL*8 U0(NX,NY,NL),V0(NX,NY,NL),W0(NX,NY,NL)
      REAL*8 U(NX,NY,NL),V(NX,NY,NL),WS(NX,NY,NZ)
      REAL*8 QT(NX,NY,NL),DQ(NX,NY,NL)
      REAL*8 DSGHL(NL)
!
      INTEGER NONOS,IORD,ISOR,IP(NX)
      REAL*8 EP,DSNP,DSSP,DTDLN(NY),DTDLT(NY),GC(NY),DTDSG(NL)
      REAL*8 UQ(NX,NY,NL),VQ(NX,NY,NL),WQ(NX,NY,NL),PQ(NX,NY)
      REAL*8 HALF,ONE
      DATA    HALF,ONE / 0.5E0,1.0E0 /
      INTEGER I,J,K
      integer begj,endj

#if (defined SPMD)
      begj = jbeg
      endj = jend
#else
      begj = 1
      endj = ny
#endif
!
!     GET THE ADVECTION VELOCITY
!
      DO K   = 1 ,NL
!        (DO J   = 1 ,NY)
          do j   = begj,endj
            DO I   = 1 ,NX
               UQ(I,J,K)  = HALF*(U (I,J,K)+U0(I,J,K))
               VQ(I,J,K)  = HALF*(V (I,J,K)+V0(I,J,K))
               WQ(I,J,K)  = HALF*(WS(I,J,K)+W0(I,J,K))
            ENDDO
         ENDDO
!
!     SAVE THE FIELD ON INPUT
!
!        (DO J   = 1 ,NY)
         do j   = begj,endj
            DO I   = 1 ,NX
               DQ(I,J,K)    = QT(I,J,K)
            ENDDO
         ENDDO
      ENDDO
!
!     PERFORM HORIZONTAL ADVECTION IN SPHERICAL GEOMETRY
!
!    (DO J   = 1 ,NY)
      do j   = begj,endj
         DO I = 1 ,NX
            PQ(I,J)    = ONE
         ENDDO
      ENDDO
!
!     DO THE 2-D ADVECTION BY MPDATA
!
      CALL MPDATA(QT,UQ,VQ,PQ,DSNP,DSSP,GC,DTDLT,DTDLN
     _           ,EP,NONOS,IORD,ISOR,IP)

!      write(*,*) 'mpdata called'
!
!     PERFORM THE VERTICAL ADVECTION
!       BY  P.K.Smolarkiewicz
      CALL VPDATA(QT,WQ,DTDSG,EP,NONOS,ISOR,IORD)

!      write(*,*) 'vpdata called'
!
!     PERFORM VERTICAL REDISTRIBUTION TO AVOID NEGATIVE Q-H2O
!
      CALL AVNEGQ(QT,DSGHL)

!      write(*,*) 'avnegq called'
      RETURN
      END
#endif

!##################################################################################
!! (2003.11.30)

      SUBROUTINE QPDATA1(QT,U0,V0,W0,DSGHL
     _                 ,U,V,WS
     _                 ,DSNP,DSSP,DTDLN,DTDLT,GC,DTDSG)

      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
!
!     PREDICT POSITIVE DEFINITE FIELD Q        DUE TO 3-D ADVECTION
!           1)  BY USING THE SCHEME     PROPOSED BY R.C.Yu
!           2)  BY USING THE SCHEME     PROPOSED BY P.K.Smolarkiewicz
!
      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

C
      REAL*8 U0(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 V0(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 W0(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 U(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 V(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 WS(beglonex:endlonex,beglatexdyn:endlatexdyn,NZ)
      REAL*8 QT(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      REAL*8 DSGHL(NL)
!
      INTEGER NONOS,IORD,ISOR
!
      REAL*8 DSNP,DSSP,DTDLN(beglatexdyn:endlatexdyn)
      REAL*8 DTDLT(beglatexdyn:endlatexdyn)
      REAL*8 GC(beglatexdyn:endlatexdyn),DTDSG(NL)
      REAL*8 HALF,ONE
      DATA    HALF,ONE / 0.5E0,1.0E0 /
      INTEGER I,J,K
      real*8 DQ(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      real*8 UQ(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      real*8 VQ(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      real*8 WQ(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)
      real*8 PQ(beglonex:endlonex,beglatexdyn:endlatexdyn)

!

!
!     GET THE ADVECTION VELOCITY
!
!$OMP PARALLEL DO PRIVATE (I,J,K)
      DO K   = 1 ,NL
         DO J   = jbeg0 , jend0
            DO I   = beglonex , endlonex
               UQ(I,J,K)  = HALF*(U (I,J,K)+U0(I,J,K))
               VQ(I,J,K)  = HALF*(V (I,J,K)+V0(I,J,K))
               WQ(I,J,K)  = HALF*(WS(I,J,K)+W0(I,J,K))
               DQ(I,J,K)    = QT(I,J,K)
            ENDDO
         ENDDO
      ENDDO
!
!     PERFORM HORIZONTAL ADVECTION IN SPHERICAL GEOMETRY
!

!$OMP PARALLEL DO PRIVATE (I,J)
      DO J = jbeg0 , jend0
         DO I = beglonex ,endlonex
            PQ(I,J)    = ONE
         ENDDO
      ENDDO
!
!     DO THE 2-D ADVECTION BY MPDATA
!
      CALL TSPAS(QT,UQ,VQ,PQ,DSNP,DSSP,GC,DTDLT,DTDLN)
!
!     PERFORM THE VERTICAL ADVECTION
!       BY  R.C.Yu
      CALL TSPASW(QT,WQ,DTDSG)
!
!     PERFORM VERTICAL REDISTRIBUTION TO AVOID NEGATIVE Q-H2O
!
      CALL AVNEGQ(QT,DSGHL)

      RETURN
      END
