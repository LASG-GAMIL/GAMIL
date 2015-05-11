#include <misc.h>
#include <params.h>

!!(wh 2003.11.11)
!!---------------------

#if (defined SPMD)

       SUBROUTINE MPDATA(Q,U,V,P,DSNP,DSSP,GC,DTDLT,DTDLN
     _                 ,EP,NONOS,IORD,ISOR,IP)
C     **************************
C     **************************
C
C     PERFORM 2-D ADVECTION IN GLOBAL SPHERICAL GEOMETRY
C             WITH THE UNIFORM LAT-LON C-GRID MESH
C

      use mpishorthand, only: mpicom

      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>
#include <commpi.h>
#include <mpif.h>

      character*50 filename

      integer isend1,irecv1
      integer isend2,irecv2
      integer istatus(mpi_status_size)

      real*8 workst(nx*nl),workrb(nx*nl)
      real*8 worksb(nx*nl),workrt(nx*nl)
      integer ii

      integer endj

!
      REAL*8  Q(NX,NY,NL),U(NX,NY,NL),V(NX,NY,NL)
     &       ,P(NX,NY)
      REAL*8  VU(NX,NY,nl),VV(NX,NY,nl),FU(NX,NY,nl),FV(NX,NY,nl)
     &       ,CU(NX,NY),CV(NX,NY),HU(NX,NY),HV(NX,NY)
     &       ,UU(NX,NY,nl),UV(NX,NY,nl),BU(NX,NY,nl),BD(NX,NY,nl)
     &       ,X (NX,NY+1,nl),H (NX,NY),XM(NX,NY,nl),XN(NX,NY,nl)
     &       ,C (NX,2 ,nl)
C
      REAL*8  EP,DSNP,DSSP,DTDLN(NY),DTDLT(NY),GC(NY)
      REAL*8  ZERO,HALF,ONE
      DATA ZERO,HALF,ONE / 0.0D0,0.5D0,1.0D0 /
      INTEGER I,J,K,IO
      INTEGER NONOS,IORD,ISOR,IP(NX)
      REAL*8  XNP(nl),XSP(nl)
C     =========================================================
      REAL*8  A,VDYF,R,X1,X2,VCOR,B,YY,VC31,VC32,PP,Y,PN
      VDYF(R,A,X1,X2)= (ABS(A)-A**2/R)*(X2-X1)/(X2+X1+EP)
      VCOR(R,A,B,YY) =-0.125*A*B*YY/R
      VC31(R,A,YY)   =-(1.0-3.0*ABS(A)/R+2.0*(A/R)**2)*A*YY/3.0
      VC32(R,A,B,YY) = 0.25*B/R*(ABS(A)-2.0*A**2/R)*YY
      PP(Y)          = MAX(0.0,Y)
      PN(Y)          =-MIN(0.0,Y)
C     =========================================================
C
C     GET WORKSHOPS INDEPENDENT OF VERTICAL COORDINATE
!-H ---
        do j=jbeg,jend
          DO I = IB,IE
            H(I,J) = GC(J) * P(I,J)
          ENDDO
          H(1 ,J)  = H(IE,J)
          H(NX,J)  = H(IB,J)
        enddo               !---------------
        if (inc_pole) then  ! for the poles
          DO I = 1 ,NX
            H(I,jpole) = ZERO
          ENDDO
        ENDIF

        call mpi_isend( h(1,2), nx,mpi_double_precision,itop,1,mpicom,isend1,ierr )
        call mpi_irecv( h(1,ny),nx,mpi_double_precision,ibot,1,mpicom,irecv1,ierr )

!-HU,CU --
        do j=jbeg,jend
          DO I = IB,IE
            HU(I,J)  = HALF*(H(I-1,J)+H(I,J))
          ENDDO
          HU(1 ,J) = HU(IE,J)
          HU(NX,J) = HU(IB,J)
          DO I = 1 ,NX
            CU(I,J) = DTDLN(J) * HU(I,J)
          ENDDO
        enddo                  !---------------
        if (inc_pole) then     ! for the poles
          DO I = 1 ,NX
            HU(I,jpole) = ZERO
            CU(I,jpole) = ZERO
          ENDDO
        ENDIF

!-HV,CV --
        call mpi_wait( isend1,istatus,ierr )
        call mpi_wait( irecv1,istatus,ierr )

        DO J = 2,jend
         DO I = IB,IE
            HV(I,J)  = HALF*(H(I,J)+H(I,J+1))
         ENDDO
            HV(1 ,J) = HV(IE,J)
            HV(NX,J) = HV(IB,J)
         DO I = 1 ,NX
            CV(I,J) = DTDLT(J) * HV(I,J)
         ENDDO
        ENDDO

!-------------------------------------------------------------
C
C     START VERTICAL LAYER DO LOOP
C$DOACROSS LOCAL(K,J,I,UU,UV,X,XM,XN,FU,FV,C,XNP,XSP,VU,VV,BU,BD,IO)

!-X ----
        do 250 k = 1 ,nl
        DO 250 J = 2 ,NY-1
        DO 250 I = 1 ,NX
           X (I,J,k)  = Q (I,J,K)
250     CONTINUE

        ii=1
        do k=1,nl
         do i=1,nx
            workst(ii)=x(i,2,   k)
            worksb(ii)=x(i,ny-1,k)
            ii=ii+1
         enddo
        enddo

        call mpi_isend(workst,nx*nl,mpi_double_precision,itop,2,mpicom,isend1,ierr)
        call mpi_isend(worksb,nx*nl,mpi_double_precision,ibot,3,mpicom,isend2,ierr)
        call mpi_irecv(workrb,nx*nl,mpi_double_precision,ibot,2,mpicom,irecv1,ierr)
        call mpi_irecv(workrt,nx*nl,mpi_double_precision,itop,3,mpicom,irecv2,ierr)
!-------
        if (inc_pole) then
          do k=1,nl
          DO I = 1 ,NX
            VU(I,Jpole,k) = ZERO
            FU(I,Jpole,k) = ZERO
          ENDDO
          enddo
        endif
!-------
C
C     CALCULATE COURANT NUMBERS MULTIPLIED BY H
!---
       do k=1,nl
        do j=jbeg,jend
         do i=ib,ie
            UU(I,J,k)  = CU(I,J) * U(I,J,K)
         enddo
            UU(1 ,J,k) = UU(IE,J,k)
            UU(NX,J,k) = UU(IB,J,k)
        enddo
       enddo
!---
       do k=1,nl
        do j=2,jend
         do i=ib,ie
            UV(I,J,k)  = CV(I,J) * V(I,J,K)
         enddo
            UV(1 ,J,k) = UV(IE,J,k)
            UV(NX,J,k) = UV(IB,J,k)
        enddo
       enddo
!---
C
       do 240 k = 1,nl
        DO 240 J = jbeg,jend
         DO 240 I = 1 ,NX
            VU(I,J,k)  = UU(I,J,k)
240    CONTINUE
!---
       do 245 k =1,nl
        DO 245 J = 2 ,jend
         DO 245 I = 1 ,NX
            VV(I,J,k)  = UV(I,J,k)
245    CONTINUE
!-----------------------------------------------------

C     PREPARE FOR NON-OSSCILATORY OPTION

        call mpi_wait( isend1,istatus,ierr )
        call mpi_wait( irecv1,istatus,ierr )
        call mpi_wait( isend2,istatus,ierr )
        call mpi_wait( irecv2,istatus,ierr )

        ii=1
        do k=1,nl
         do i=1,nx
            x(i,1, k)=workrt(ii)
            x(i,ny,k)=workrb(ii)
            ii=ii+1
         enddo
        enddo

      IF( NONOS.EQ.1 ) THEN
        do 300 k = 1,nl
        DO 300 J = JBeg,JEnd
        DO 300 I = IB,IE
          XM(I,J,k)  = MAX( X(I-1,J,k),X(I,J,k),X(I+1,J,k),X(I,J-1,k),X(I,J+1,k) )
          XN(I,J,k)  = MIN( X(I-1,J,k),X(I,J,k),X(I+1,J,k),X(I,J-1,k),X(I,J+1,k) )
300     CONTINUE
      ENDIF
C
      DO 700 IO= 1 ,IORD
C     ++++++++++++++++++++++++++++++++
C     PREDICTOR STEP : UPSTREAM SCHEME
C     ++++++++++++++++++++++++++++++++
      do 410 k = 1,nl
      DO 410 J = JBeg,JEnd
       DO 400 I = IB,IE
        IF(VU(I,J,k).GE.0.0E0) THEN
         FU(I,J,k)= X(I-1,J,k)*VU(I,J,k)
        ELSE
         FU(I,J,k)= X(I,J,k)*VU(I,J,k)
        ENDIF
400    CONTINUE
       FU(NX,J,k) = FU(IB,J,k)
410   CONTINUE

      do 420 k = 1, nl
      DO 420 J = 2 ,JEnd
      DO 420 I = IB,IE
        IF(VV(I,J,k).GE.0.0E0) THEN
          FV(I,J,k)= X(I,J,k)*VV(I,J,k)
        ELSE
          FV(I,J,k)= X(I,J+1,k)*VV(I,J,k)
        ENDIF
420   CONTINUE

      ii=1
      do k=1,nl
       do i=1,nx
          worksb(ii)=fv(i,ny-1,k)
          ii=ii+1
       enddo
      enddo
        call mpi_isend(worksb,nx*nl,mpi_double_precision,ibot,4,mpicom,isend1,ierr)
        call mpi_irecv(workrt,nx*nl,mpi_double_precision,itop,4,mpicom,irecv1,ierr)
        call mpi_wait(isend1,istatus,ierr)
        call mpi_wait(irecv1,istatus,ierr)
      ii=1
      do k=1,nl
       do i=1,nx
          fv(i,1,k)=workrt(ii)
          ii=ii+1
       enddo
      enddo

      do 430 k = 1,nl
      DO 430 J = JBeg,JEnd
      DO 430 I = IB,IE
      X(I,J,k)   = X(I,J,k)-(FU(I+1,J,k)-FU(I,J,k)+FV(I,J,k)-FV(I,J-1,k))/H(I,J)
430   CONTINUE


C     B.C. BY THE SPHERICAL CYCLICITY & MASS CONSERVATION RESTRICT

      if (nprocs.ge.32) then
         if (myrank.eq.1)
     _     call mpi_send( h(1,ny-1),nx,mpi_double_precision,0,5,mpicom,ierr)
         if (myrank.eq.0)
     _     call mpi_recv( h(1,1),   nx,mpi_double_precision,1,5,mpicom,ierr)
      endif

      if (myrank.eq.(nprocs-1)) then
         do  k = 1,nl
          XNP(k)      = ZERO
          do I = IB,IE
             C(I,1,k)= FV(I,jpole,k) / H(I,jpole+1)
             XNP(k)  = XNP(k) + C(I,1,k)
          enddo
             XNP(k)      = X(IB,jpole,k) - XNP(k)*DSNP
          do i=ib,ie
             X(I,jpole,k)   = XNP(k)
          enddo
         enddo
      endif

      if (myrank.eq.0) then
         do k=1,nl
          XSP(k)      = ZERO
          do I = IB,IE
             C(I,2,k) = FV(I,jpole-1,k) / H(I,jpole-1)
             XSP(k)   = XSP(k) + C(I,2,k)
          enddo
             XSP(k)   = X(IB,jpole,k) + XSP(k)*DSSP
          do i=ib,ie
             X(I,jpole,k)  = XSP(k)
          enddo
         enddo
      endif

        do 470 k = 1,nl
        DO 470 J = 2 ,NY-1
        X(1 ,J,k)  = X(IE,J,k)
        X(NX,J,k)  = X(IB,J,k)
470     CONTINUE
C
      IF( IO.EQ.IORD ) GOTO 700
C     ++++++++++++++++++++++++++++++++++++++
C     CORRECTOR STEP : ANTI-DIFFUSION SCHEME
C     ++++++++++++++++++++++++++++++++++++++

        ii=1
        do k=1,nl
         do i=1,nx
            workst(ii)=x(i,2,   k)
            worksb(ii)=x(i,ny-1,k)
            ii=ii+1
         enddo
        enddo
        call mpi_isend(workst,nx*nl,mpi_double_precision,itop,6,mpicom,isend1,ierr)
        call mpi_isend(worksb,nx*nl,mpi_double_precision,ibot,7,mpicom,isend2,ierr)
        call mpi_irecv(workrb,nx*nl,mpi_double_precision,ibot,6,mpicom,irecv1,ierr)
        call mpi_irecv(workrt,nx*nl,mpi_double_precision,itop,7,mpicom,irecv2,ierr)
!-----

      do 500 k = 1,nl
      DO 500 J = JBeg,JEnd
      DO 500 I = 1 ,NX
      FU(I,J,k)  = VU(I,J,k)
500   CONTINUE

      do 505 k = 1,nl
      DO 505 J = 2 ,JEnd
      DO 505 I = 1 ,NX
      FV(I,J,k)  = VV(I,J,k)
505   CONTINUE
!-----
        call mpi_wait( isend1,istatus,ierr )
        call mpi_wait( irecv1,istatus,ierr )
        call mpi_wait( isend2,istatus,ierr )
        call mpi_wait( irecv2,istatus,ierr )

        ii=1
        do k=1,nl
         do i=1,nx
            x(i,1, k)=workrt(ii)
            x(i,ny,k)=workrb(ii)
            ii=ii+1
         enddo
        enddo
!-----

C     CALCULATE THE  PSEUDO VELOCITIES
C               LONGITUDINAL DIRECTION

        ii=1
        do k=1,nl
         do i=1,nx
            workst(ii)=fu(i,2,   k)
            worksb(ii)=fv(i,ny-1,k)
            ii=ii+1
         enddo
        enddo

        call mpi_isend(workst,nx*nl,mpi_double_precision,itop,8,mpicom,isend1,ierr)
        call mpi_isend(worksb,nx*nl,mpi_double_precision,ibot,9,mpicom,isend2,ierr)
        call mpi_irecv(workrb,nx*nl,mpi_double_precision,ibot,8,mpicom,irecv1,ierr)
        call mpi_irecv(workrt,nx*nl,mpi_double_precision,itop,9,mpicom,irecv2,ierr)

        call mpi_wait( isend1,istatus,ierr )
        call mpi_wait( irecv1,istatus,ierr )
        call mpi_wait( isend2,istatus,ierr )
        call mpi_wait( irecv2,istatus,ierr )

        ii=1
        do k=1,nl
         do i=1,nx
            fv(i,1, k)=workrt(ii)
            fu(i,ny,k)=workrb(ii)
            ii=ii+1
         enddo
        enddo

      do 510 k = 1,nl
      DO 510 J = JBeg,JEnd
      DO 510 I = IB,IE
      VU(I,J,k)  = VDYF( HU(I,J),FU(I,J,k),X(I-1,J,k),X(I,J,k) )
     &         + VCOR( HU(I,J),FU(I,J,k)
     &         ,      (FV(I-1,J-1,k)+FV(I-1,J,k)+FV(I,J,k)+FV(I,J-1,k))
     &         ,      (X(I-1,J+1,k)+X(I,J+1,k)-X(I-1,J-1,k)-X(I,J-1,k))
     &         /   (EP+X(I-1,J+1,k)+X(I,J+1,k)+X(I-1,J-1,k)+X(I,J-1,k)) )
510   CONTINUE
C               LATITUDINAL  DIRECTION

      do 520 k = 1,nl
      DO 520 J = 2 ,JEnd
      DO 520 I = IB,IE
      VV(I,J,k)  = VDYF( HV(I,J),FV(I,J,k),X(I,J,k),X(I,J+1,k) )
     &         + VCOR( HV(I,J),FV(I,J,k)
     &         ,      (FU(I,J,k)+FU(I,J+1,k)+FU(I+1,J,k)+FU(I+1,J+1,k))
     &         ,      (X(I+1,J,k)+X(I+1,J+1,k)-X(I-1,J,k)-X(I-1,J+1,k))
     &         /   (EP+X(I+1,J,k)+X(I+1,J+1,k)+X(I-1,J,k)+X(I-1,J+1,k)) )
520   CONTINUE

!---------------------------------------------------------------
C     ADD THE THIRD ORDER CORRECTION IF REQUESTED
      IF( ISOR.EQ.3 ) THEN
C               LONGITUDINAL DIRECTION

      ii=1
      do k=1,nl
       do i=1,nx
          workst(ii)=x(i,3,k)
          ii=ii+1
       enddo
      enddo
      call mpi_isend(workst,nx*nl,mpi_double_precision,itop,10,mpicom,isend1,ierr)
      call mpi_irecv(workrb,nx*nl,mpi_double_precision,ibot,10,mpicom,irecv1,ierr)
!-----
      do 530 k = 1,nl
      DO 530 J = JBeg,JEnd
      DO 530 I = 3 ,IE
      VU(I,J,k)  = VU(I,J,k) + VC31( HU(I,J),FU(I,J,k)
     &         ,    (X(I-2,J,k)+X(I+1,J,k)-X(I-1,J,k)-X(I,J,k))
     &         / (EP+X(I-2,J,k)+X(I+1,J,k)+X(I-1,J,k)+X(I,J,k)) )
530   CONTINUE

      do 535 k=1,nl
      DO 535 J = JBeg,JEnd
      VU(2,J,k)  = VU(2,J,k) + VC31( HU(2,J),FU(2,J,k)
     &         ,    (X(IM,J,k)+X(3,J,k)-X(1,J,k)-X(2,J,k))
     &         / (EP+X(IM,J,k)+X(3,J,k)+X(1,J,k)+X(2,J,k)) )
535   CONTINUE

      do 540 k = 1,nl
      DO 540 J = JBeg,JEnd
      DO 540 I = IB,IE
      VU(I,J,k)  = VU(I,J,k) + VC32( HU(I,J),FU(I,J,k)
     &         ,    (FV(I-1,J-1,k)+FV(I-1,J,k)+FV(I,J-1,k)+FV(I,J,k))
     &         ,    (X(I,J+1,k)-X(I,J-1,k)-X(I-1,J+1,k)+X(I-1,J-1,k))
     &         / (EP+X(I,J+1,k)+X(I,J-1,k)+X(I-1,J+1,k)+X(I-1,J-1,k)) )
540   CONTINUE

C               LATITUDINAL  DIRECTION

      call mpi_wait(isend1,istatus,ierr)
      call mpi_wait(irecv1,istatus,ierr)

      ii=1
      do k=1,nl
       do i=1,nx
          x(i,ny+1,k)=workrb(ii)
          ii=ii+1
       enddo
      enddo

      endj=jend
      if (myrank.eq.0) endj=jend-1
      if ( (myrank.eq.1).and.(nprocs.ge.32) ) endj=jend-1

      do 550 k = 1,nl
      DO 550 J = JBeg,endj
      DO 550 I = IB,IE
      VV(I,J,k)  = VV(I,J,k) + VC31( HV(I,J),FV(I,J,k)
     &         ,    (X(I,J-1,k)+X(I,J+2,k)-X(I,J,k)-X(I,J+1,k))
     &         / (EP+X(I,J-1,k)+X(I,J+2,k)+X(I,J,k)+X(I,J+1,k)) )
550   CONTINUE

!--------------------------------------
C     B.C. BY THE SPHERICAL CYCLICITY

      if (myrank.eq.nprocs-1) then
        do k=1,nl
         do i=ib,ie
            C(I,1,k)  = X(IP(I),jpole+1,k)
           VV(I,jpole,k)  = VV(I,jpole,k) + VC31( HV(I,jpole),FV(I,jpole,k)
     &             ,    (C(I,1,k)+X(I,jpole+2,k)-X(I,jpole,k)-X(I,jpole+1,k))
     &             / (EP+C(I,1,k)+X(I,jpole+2,k)+X(I,jpole,k)+X(I,jpole+1,k)) )
         enddo
        enddo
      endif
!--
      if (nprocs.le.30) then
         if (myrank.eq.0) then
            do k=1,nl
            do i=ib,ie
               C(I,2,k)   = X(IP(I),jpole-1,k)
              VV(I,jpole-1,k)  = VV(I,jpole-1,k) + VC31( HV(I,jpole-1),FV(I,jpole-1,k)
     &           ,    (X(I,jpole-2,k)+C(I,2,k)-X(I,jpole-1,k)-X(I,jpole,k))
     &           / (EP+X(I,jpole-2,k)+C(I,2,k)+X(I,jpole-1,k)+X(I,jpole,k)) )
            enddo
            enddo
         endif
      else
         if (myrank.eq.1) then
            do k=1,nl
            do i=ib,ie
               C(I,2,k)   = X(IP(I),jend,k)
               VV(I,JEnd,k)  = VV(I,JEnd,k) + VC31( HV(I,JEnd),FV(I,JEnd,k)
     &         ,    (X(I,JEnd-1,k)+C(I,2,k)-X(I,JEnd,k)-X(I,NY,k))
     &         / (EP+X(I,JEnd-1,k)+C(I,2,k)+X(I,JEnd,k)+X(I,NY,k)) )
            enddo
            enddo
         endif
      endif
!--------------------


      do 565 k = 1,nl
      DO 565 J = 2 ,JEnd
      DO 565 I = IB,IE
      VV(I,J,k)  = VV(I,J,k) + VC32( HV(I,J),FV(I,J,k)
     &         ,    (FU(I,J,k)+FU(I+1,J,k)+FU(I+1,J+1,k)+FU(I,J+1,k))
     &         ,    (X(I+1,J+1,k)-X(I-1,J+1,k)-X(I+1,J,k)+X(I-1,J,k))
     &         / (EP+X(I+1,J+1,k)+X(I-1,J+1,k)+X(I+1,J,k)+X(I-1,J,k)) )
565   CONTINUE

      ENDIF
C
      DO 570 k = 1,nl
      DO 570 J = JBeg,JEnd
      DO 570 I = IB,IE
      VU(I,J,k)  = SIGN(ONE,VU(I,J,k))*MIN(ABS(UU(I,J,k)),ABS(VU(I,J,k)))
570   CONTINUE

      DO 580 k = 1 ,nl
      DO 580 J = 2 ,JEnd
      DO 580 I = IB,IE
      VV(I,J,k)  = SIGN(ONE,VV(I,J,k))*MIN(ABS(UV(I,J,k)),ABS(VV(I,J,k)))
580   CONTINUE

C     B.C. BY THE SPHERICAL CYCLICITY

      DO 590 k = 1,nl
      DO 590 J = JBeg,JEnd
      VU(1 ,J,k) = VU(IE,J,k)
      VU(NX,J,k) = VU(IB,J,k)
590   CONTINUE

      DO 595 k = 1 ,nl
      DO 595 J = 2 ,JEnd
      VV(1 ,J,k) = VV(IE,J,k)
      VV(NX,J,k) = VV(IB,J,k)
595   CONTINUE
C
C     PERFORM THE NON-OSSCILATORY OPTION
      IF( NONOS.EQ.1 ) THEN

      DO 600 k = 1,nl
      DO 600 J = JBeg,JEnd
      DO 600 I = IB,IE
      XM(I,J,k)  = MAX( X(I-1,J,k),X(I,J,k),X(I+1,J,k),X(I,J-1,k),X(I,J+1,k)
     &         ,        XM(I,J,k) )
      XN(I,J,k)  = MIN( X(I-1,J,k),X(I,J,k),X(I+1,J,k),X(I,J-1,k),X(I,J+1,k)
     &         ,        XN(I,J,k) )
600   CONTINUE
C
      DO 620 k = 1,nl
      DO 620 J = JBeg,JEnd
      DO 620 I = IB,IE
Cb    FU(I,J,k)  = DONOR( X(I-1,J,k),X(I,J,k),VU(I,J,k) )
      IF(VU(I,J,k).GE.0.0E0) THEN
        FU(I,J,k)= X(I-1,J,k)*VU(I,J,k)
      ELSE
        FU(I,J,k)= X(I,J,k)*VU(I,J,k)
      ENDIF
620   CONTINUE

      DO 625 k = 1,nl
      DO 625 J = JBeg,JEnd
      FU(NX,J,k) = FU(IB,J,k)
625   CONTINUE

      DO 630 k = 1 ,nl
      DO 630 J = 2 ,JEnd
      DO 630 I = IB,IE
      IF(VV(I,J,k).GE.0.0E0) THEN
        FV(I,J,k)= X(I,J,k)*VV(I,J,k)
      ELSE
        FV(I,J,k)= X(I,J+1,k)*VV(I,J,k)
      ENDIF
630   CONTINUE

!--------------
       ii=1
       do k=1,nl
        do i=1,nx
           worksb(ii)=fv(i,ny-1,k)
           ii=ii+1
        enddo
       enddo
        call mpi_isend(worksb,nx*nl,mpi_double_precision,ibot,12,mpicom,isend1,ierr)
        call mpi_irecv(workrt,nx*nl,mpi_double_precision,itop,12,mpicom,irecv1,ierr)
        call mpi_wait(isend1,istatus,ierr)
        call mpi_wait(irecv1,istatus,ierr)
       ii=1
       do k=1,nl
        do i=1,nx
           fv(i,1,k)=workrt(ii)
           ii=ii+1
        enddo
       enddo

      DO 640 k = 1,nl
      DO 640 J = JBeg,JEnd
      DO 640 I = IB,IE
      BU(I,J,k)  = (XM(I,J,k)-X(I,J,k))*H(I,J) /
     &  (PN(FU(I+1,J,k))+PP(FU(I,J,k))+PN(FV(I,J,k))+PP(FV(I,J-1,k))+EP)
      BD(I,J,k)  = (X(I,J,k)-XN(I,J,k))*H(I,J) /
     &  (PP(FU(I+1,J,k))+PN(FU(I,J,k))+PP(FV(I,J,k))+PN(FV(I,J-1,k))+EP)
640   CONTINUE

      DO 645 k = 1,nl
      DO 645 J = JBeg,JEnd
      BU(1,J,k)  = BU(IE,J,k)
      BD(1,J,k)  = BD(IE,J,k)
645   CONTINUE
C     KEEP IN MIND THAT H ARE ZERO AT POLE

      if (inc_pole) then
         do k=1,nl
          do i=1,ie
           BU(I,jpole,k) = ZERO
           BD(I,jpole,k) = ZERO
          enddo
         enddo
      endif
!----
      ii=1
      do k=1,nl
       do i=1,nx
          workst(ii)=bu(i,2,k)
          worksb(ii)=bd(i,2,k)
          ii=ii+1
       enddo
      enddo

        call mpi_isend(workst,nx*nl,mpi_double_precision,itop,13,mpicom,isend1,ierr)
        call mpi_isend(worksb,nx*nl,mpi_double_precision,itop,14,mpicom,isend2,ierr)
        call mpi_irecv(workrb,nx*nl,mpi_double_precision,ibot,13,mpicom,irecv1,ierr)
        call mpi_irecv(workrt,nx*nl,mpi_double_precision,ibot,14,mpicom,irecv2,ierr)
!----
      DO 660 k = 1,nl
      DO 660 J = JBeg,JEnd
      DO 660 I = IB,IE
      VU(I,J,k)  = PP( VU(I,J,k) ) * MIN(ONE,BD(I-1,J,k),BU(I,J,k))
     &         - PN( VU(I,J,k) ) * MIN(ONE,BU(I-1,J,k),BD(I,J,k))
660   CONTINUE
!----
        call mpi_wait( isend1,istatus,ierr )
        call mpi_wait( irecv1,istatus,ierr )
        call mpi_wait( isend2,istatus,ierr )
        call mpi_wait( irecv2,istatus,ierr )

      ii=1
      do k=1,nl
       do i=1,nx
          bu(i,ny,k)=workrb(ii)
          bd(i,ny,k)=workrt(ii)
          ii=ii+1
       enddo
      enddo
!-----
      DO 670 k = 1 ,nl
      DO 670 J = 2 ,JEnd
      DO 670 I = IB,IE
      VV(I,J,k)  = PP( VV(I,J,k) ) * MIN(ONE,BD(I,J,k),BU(I,J+1,k))
     &         - PN( VV(I,J,k) ) * MIN(ONE,BU(I,J,k),BD(I,J+1,k))
670   CONTINUE

C     B.C. BY THE SPHERICAL CYCLICITY

      DO 680 k = 1,nl
      DO 680 J = JBeg,JEnd
      VU(1 ,J,k) = VU(IE,J,k)
      VU(NX,J,k) = VU(IB,J,k)
680   CONTINUE

      do 685 k = 1,nl
      DO 685 J = 2 ,JEnd
      VV(1 ,J,k) = VV(IE,J,k)
      VV(NX,J,k) = VV(IB,J,k)
685   CONTINUE

      ENDIF

700   CONTINUE
C


C     UPDATE THE PREDICTED FIELD
      do 800 k = 1,nl
      DO 800 J = 2 ,NY-1
      DO 800 I = 1 ,NX
      Q(I,J,K) = X(I,J,k)
800   CONTINUE

!- check ---------------------------------------------------------
!
!      write(filename,14) 'mpdata-p-',myrank,'.out'
!14    format(a9,i1,a4)
!
!      open (10,file=trim(filename))
!
!      do j=1,ny
!        write(10,11) j,q(1,j,10),q(2,j,10)
!      enddo
!
!11    format(1x,i5,2e30.20)
!      close (10)
!
!      call mpi_finalize(j)
!      stop'mpdata'
!--------------------------------------------------------------

      RETURN
	END

#else


      SUBROUTINE MPDATA(Q,U,V,P,DSNP,DSSP,GC,DTDLT,DTDLN
     _                 ,EP,NONOS,IORD,ISOR,IP)
C     **************************
C     **************************
C
C     PERFORM 2-D ADVECTION IN GLOBAL SPHERICAL GEOMETRY
C             WITH THE UNIFORM LAT-LON C-GRID MESH
C
      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

!
      REAL*8  Q(NX,NY,NL),U(NX,NY,NL),V(NX,NY,NL)
     &       ,P(NX,NY)
      REAL*8  VU(NX,NY),VV(NX,NY),FU(NX,NY),FV(NX,NY)
     &       ,CU(NX,NY),CV(NX,NY),HU(NX,NY),HV(NX,NY)
     &       ,UU(NX,NY),UV(NX,NY),BU(NX,NY),BD(NX,NY)
     &       ,X (NX,NY),H (NX,NY),XM(NX,NY),XN(NX,NY)
     &       ,C (NX,2 )
C
      REAL*8  EP,DSNP,DSSP,DTDLN(NY),DTDLT(NY),GC(NY)
      REAL*8  ZERO,HALF,ONE
      DATA ZERO,HALF,ONE / 0.0D0,0.5D0,1.0D0 /
      INTEGER I,J,K,IO
      INTEGER NONOS,IORD,ISOR,IP(NX)
      REAL*8  XNP,XSP
C     =========================================================
      REAL*8  A,VDYF,R,X1,X2,VCOR,B,YY,VC31,VC32,PP,Y,PN
      VDYF(R,A,X1,X2)= (ABS(A)-A**2/R)*(X2-X1)/(X2+X1+EP)
      VCOR(R,A,B,YY) =-0.125*A*B*YY/R
      VC31(R,A,YY)   =-(1.0-3.0*ABS(A)/R+2.0*(A/R)**2)*A*YY/3.0
      VC32(R,A,B,YY) = 0.25*B/R*(ABS(A)-2.0*A**2/R)*YY
      PP(Y)          = MAX(0.0,Y)
      PN(Y)          =-MIN(0.0,Y)
C     =========================================================
C
C     GET WORKSHOPS INDEPENDENT OF VERTICAL COORDINATE
      DO J = 1 ,NY
        IF(J.GE.JB.AND.J.LE.JE) THEN
          DO I = IB,IE
            H(I,J) = GC(J) * P(I,J)
          ENDDO
          H(1 ,J)  = H(IE,J)
          H(NX,J)  = H(IB,J)
        ELSE
          DO I = 1 ,NX
            H(I,J) = ZERO
          ENDDO
        ENDIF
        IF(J.GE.JB.AND.J.LE.JE) THEN
          DO I = IB,IE
            HU(I,J)  = HALF*(H(I-1,J)+H(I,J))
          ENDDO
          H(1 ,J)  = H(IE,J)
          H(NX,J)  = H(IB,J)
          HU(1 ,J) = HU(IE,J)
          HU(NX,J) = HU(IB,J)
          DO I = 1 ,NX
            CU(I,J) = DTDLN(J) * HU(I,J)
          ENDDO
        ELSE
          DO I = 1 ,NX
            HU(I,J) = ZERO
            CU(I,J) = ZERO
          ENDDO
        ENDIF
      ENDDO
      DO J = 1 ,JE
        DO I = IB,IE
          HV(I,J)  = HALF*(H(I,J)+H(I,J+1))
        ENDDO
        HV(1 ,J) = HV(IE,J)
        HV(NX,J) = HV(IB,J)
        DO I = 1 ,NX
          CV(I,J) = DTDLT(J) * HV(I,J)
        ENDDO
      ENDDO
C
C     START VERTICAL LAYER DO LOOP
C$DOACROSS LOCAL(K,J,I,UU,UV,X,XM,XN,FU,FV,C,XNP,XSP,VU,VV,BU,BD,IO)
      DO 800 K = 1 ,NL
        DO J = 1 ,NY,JE
          DO I = 1 ,NX
            VU(I,J) = ZERO
            FU(I,J) = ZERO
          ENDDO
        ENDDO
C
C     CALCULATE COURANT NUMBERS MULTIPLIED BY H
      DO 200 J = JB,JE
      DO 200 I = 1 ,NX
      UU(I,J)  = CU(I,J) * U(I,J,K)
200   CONTINUE
      DO 210 J = 1 ,JE
      DO 210 I = 1 ,NX
      UV(I,J)  = CV(I,J) * V(I,J,K)
210   CONTINUE
      DO 220 J = JB,JE
      UU(1 ,J) = UU(IE,J)
      UU(NX,J) = UU(IB,J)
220   CONTINUE
      DO 230 J = 1 ,JE
      UV(1 ,J) = UV(IE,J)
      UV(NX,J) = UV(IB,J)
230   CONTINUE
C
      DO 240 J = JB,JE
      DO 240 I = 1 ,NX
      VU(I,J)  = UU(I,J)
240   CONTINUE
      DO 245 J = 1 ,JE
      DO 245 I = 1 ,NX
      VV(I,J)  = UV(I,J)
245   CONTINUE
      DO 250 J = 1 ,NY
      DO 250 I = 1 ,NX
      X (I,J)  = Q (I,J,K)
250   CONTINUE
C
C     PREPARE FOR NON-OSSCILATORY OPTION
      IF( NONOS.EQ.1 ) THEN
      DO 300 J = JB,JE
      DO 300 I = IB,IE
      XM(I,J)  = MAX( X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1) )
      XN(I,J)  = MIN( X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1) )
300   CONTINUE
      ENDIF
C
      DO 700 IO= 1 ,IORD
C     ++++++++++++++++++++++++++++++++
C     PREDICTOR STEP : UPSTREAM SCHEME
C     ++++++++++++++++++++++++++++++++
      DO 400 J = JB,JE
      DO 400 I = IB,IE
Cb    FU(I,J)  = DONOR( X(I-1,J),X(I,J),VU(I,J) )
      IF(VU(I,J).GE.0.0E0) THEN
        FU(I,J)= X(I-1,J)*VU(I,J)
      ELSE
        FU(I,J)= X(I,J)*VU(I,J)
      ENDIF
400   CONTINUE
      DO 410 J = JB,JE
      FU(NX,J) = FU(IB,J)
410   CONTINUE
      DO 420 J = 1 ,JE
      DO 420 I = IB,IE
Cb    FV(I,J)  = DONOR( X(I,J),X(I,J+1),VV(I,J) )
      IF(VV(I,J).GE.0.0E0) THEN
        FV(I,J)= X(I,J)*VV(I,J)
      ELSE
        FV(I,J)= X(I,J+1)*VV(I,J)
      ENDIF
420   CONTINUE

      DO 430 J = JB,JE
      DO 430 I = IB,IE
      X(I,J)   = X(I,J)-(FU(I+1,J)-FU(I,J)+FV(I,J)-FV(I,J-1))/H(I,J)
430   CONTINUE
C     B.C. BY THE SPHERICAL CYCLICITY & MASS CONSERVATION RESTRICT
      DO 440 I = IB,IE
      C(I,1)   = FV(I,1) / H(I,2)
      C(I,2)   = FV(I,JE) / H(I,JE)
440   CONTINUE
      XNP      = ZERO
      XSP      = ZERO
      DO 450 I = IB,IE
      XNP      = XNP     + C(I,1)
      XSP      = XSP     + C(I,2)
450   CONTINUE
      XNP      = X(IB,1) - XNP*DSNP
      XSP      = X(IB,NY) + XSP*DSSP
      DO 460 I = IB,IE
      X(I,1)   = XNP
      X(I,NY)   = XSP
460   CONTINUE
      DO 470 J = 1 ,NY
      X(1 ,J)  = X(IE,J)
      X(NX,J)  = X(IB,J)
470   CONTINUE
C
      IF( IO.EQ.IORD ) GOTO 700
C     ++++++++++++++++++++++++++++++++++++++
C     CORRECTOR STEP : ANTI-DIFFUSION SCHEME
C     ++++++++++++++++++++++++++++++++++++++
      DO 500 J = JB,JE
      DO 500 I = 1 ,NX
      FU(I,J)  = VU(I,J)
500   CONTINUE
      DO 505 J = 1 ,JE
      DO 505 I = 1 ,NX
      FV(I,J)  = VV(I,J)
505   CONTINUE
C
C     CALCULATE THE  PSEUDO VELOCITIES
C               LONGITUDINAL DIRECTION
      DO 510 J = JB,JE
      DO 510 I = IB,IE
      VU(I,J)  = VDYF( HU(I,J),FU(I,J),X(I-1,J),X(I,J) )
     &         + VCOR( HU(I,J),FU(I,J)
     &         ,      (FV(I-1,J-1)+FV(I-1,J)+FV(I,J)+FV(I,J-1))
     &         ,      (X(I-1,J+1)+X(I,J+1)-X(I-1,J-1)-X(I,J-1))
     &         /   (EP+X(I-1,J+1)+X(I,J+1)+X(I-1,J-1)+X(I,J-1)) )
510   CONTINUE
C               LATITUDINAL  DIRECTION
      DO 520 J = 1 ,JE
      DO 520 I = IB,IE
      VV(I,J)  = VDYF( HV(I,J),FV(I,J),X(I,J),X(I,J+1) )
     &         + VCOR( HV(I,J),FV(I,J)
     &         ,      (FU(I,J)+FU(I,J+1)+FU(I+1,J)+FU(I+1,J+1))
     &         ,      (X(I+1,J)+X(I+1,J+1)-X(I-1,J)-X(I-1,J+1))
     &         /   (EP+X(I+1,J)+X(I+1,J+1)+X(I-1,J)+X(I-1,J+1)) )
520   CONTINUE
C
C     ADD THE THIRD ORDER CORRECTION IF REQUESTED
      IF( ISOR.EQ.3 ) THEN
C               LONGITUDINAL DIRECTION
      DO 530 J = JB,JE
      DO 530 I = 3 ,IE
      VU(I,J)  = VU(I,J) + VC31( HU(I,J),FU(I,J)
     &         ,    (X(I-2,J)+X(I+1,J)-X(I-1,J)-X(I,J))
     &         / (EP+X(I-2,J)+X(I+1,J)+X(I-1,J)+X(I,J)) )
530   CONTINUE
      DO 535 J = JB,JE
      VU(2,J)  = VU(2,J) + VC31( HU(2,J),FU(2,J)
     &         ,    (X(IM,J)+X(3,J)-X(1,J)-X(2,J))
     &         / (EP+X(IM,J)+X(3,J)+X(1,J)+X(2,J)) )
535   CONTINUE
      DO 540 J = JB,JE
      DO 540 I = IB,IE
      VU(I,J)  = VU(I,J) + VC32( HU(I,J),FU(I,J)
     &         ,    (FV(I-1,J-1)+FV(I-1,J)+FV(I,J-1)+FV(I,J))
     &         ,    (X(I,J+1)-X(I,J-1)-X(I-1,J+1)+X(I-1,J-1))
     &         / (EP+X(I,J+1)+X(I,J-1)+X(I-1,J+1)+X(I-1,J-1)) )
540   CONTINUE
C               LATITUDINAL  DIRECTION
      DO 550 J = JB,JE-1
      DO 550 I = IB,IE
      VV(I,J)  = VV(I,J) + VC31( HV(I,J),FV(I,J)
     &         ,    (X(I,J-1)+X(I,J+2)-X(I,J)-X(I,J+1))
     &         / (EP+X(I,J-1)+X(I,J+2)+X(I,J)+X(I,J+1)) )
550   CONTINUE
C     B.C. BY THE SPHERICAL CYCLICITY
      DO 555 I = IB,IE
      C(I,1)   = X(IP(I),2)
      C(I,2)   = X(IP(I),JE)
555   CONTINUE
      DO 560 I = IB,IE
      VV(I,1)  = VV(I,1) + VC31( HV(I,1),FV(I,1)
     &         ,    (C(I,1)+X(I,3)-X(I,1)-X(I,2))
     &         / (EP+C(I,1)+X(I,3)+X(I,1)+X(I,2)) )
      VV(I,JE)  = VV(I,JE) + VC31( HV(I,JE),FV(I,JE)
     &         ,    (X(I,JE-1)+C(I,2)-X(I,JE)-X(I,NY))
     &         / (EP+X(I,JE-1)+C(I,2)+X(I,JE)+X(I,NY)) )
560   CONTINUE
      DO 565 J = 1 ,JE
      DO 565 I = IB,IE
      VV(I,J)  = VV(I,J) + VC32( HV(I,J),FV(I,J)
     &         ,    (FU(I,J)+FU(I+1,J)+FU(I+1,J+1)+FU(I,J+1))
     &         ,    (X(I+1,J+1)-X(I-1,J+1)-X(I+1,J)+X(I-1,J))
     &         / (EP+X(I+1,J+1)+X(I-1,J+1)+X(I+1,J)+X(I-1,J)) )
565   CONTINUE
      ENDIF
C
      DO 570 J = JB,JE
      DO 570 I = IB,IE
      VU(I,J)  = SIGN(ONE,VU(I,J))*MIN(ABS(UU(I,J)),ABS(VU(I,J)))
570   CONTINUE
      DO 580 J = 1 ,JE
      DO 580 I = IB,IE
      VV(I,J)  = SIGN(ONE,VV(I,J))*MIN(ABS(UV(I,J)),ABS(VV(I,J)))
580   CONTINUE
C     B.C. BY THE SPHERICAL CYCLICITY
      DO 590 J = JB,JE
      VU(1 ,J) = VU(IE,J)
      VU(NX,J) = VU(IB,J)
590   CONTINUE
      DO 595 J = 1 ,JE
      VV(1 ,J) = VV(IE,J)
      VV(NX,J) = VV(IB,J)
595   CONTINUE
C
C     PERFORM THE NON-OSSCILATORY OPTION
      IF( NONOS.EQ.1 ) THEN
      DO 600 J = JB,JE
      DO 600 I = IB,IE
      XM(I,J)  = MAX( X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1)
     &         ,        XM(I,J) )
      XN(I,J)  = MIN( X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1)
     &         ,        XN(I,J) )
600   CONTINUE
C
      DO 620 J = JB,JE
      DO 620 I = IB,IE
Cb    FU(I,J)  = DONOR( X(I-1,J),X(I,J),VU(I,J) )
      IF(VU(I,J).GE.0.0E0) THEN
        FU(I,J)= X(I-1,J)*VU(I,J)
      ELSE
        FU(I,J)= X(I,J)*VU(I,J)
      ENDIF
620   CONTINUE
      DO 625 J = JB,JE
      FU(NX,J) = FU(IB,J)
625   CONTINUE
      DO 630 J = 1 ,JE
      DO 630 I = IB,IE
Cb    FV(I,J)  = DONOR( X(I,J),X(I,J+1),VV(I,J) )
      IF(VV(I,J).GE.0.0E0) THEN
        FV(I,J)= X(I,J)*VV(I,J)
      ELSE
        FV(I,J)= X(I,J+1)*VV(I,J)
      ENDIF
630   CONTINUE
      DO 640 J = JB,JE
      DO 640 I = IB,IE
      BU(I,J)  = (XM(I,J)-X(I,J))*H(I,J) /
     &  (PN(FU(I+1,J))+PP(FU(I,J))+PN(FV(I,J))+PP(FV(I,J-1))+EP)
      BD(I,J)  = (X(I,J)-XN(I,J))*H(I,J) /
     &  (PP(FU(I+1,J))+PN(FU(I,J))+PP(FV(I,J))+PN(FV(I,J-1))+EP)
640   CONTINUE
      DO 645 J = JB,JE
      BU(1,J)  = BU(IE,J)
      BD(1,J)  = BD(IE,J)
645   CONTINUE
C     KEEP IN MIND THAT H ARE ZERO AT POLE
      DO 650 I = 1 ,IE
      BU(I,01) = ZERO
      BU(I,NY) = ZERO
      BD(I,01) = ZERO
      BD(I,NY) = ZERO
650   CONTINUE
      DO 660 J = JB,JE
      DO 660 I = IB,IE
      VU(I,J)  = PP( VU(I,J) ) * MIN(ONE,BD(I-1,J),BU(I,J))
     &         - PN( VU(I,J) ) * MIN(ONE,BU(I-1,J),BD(I,J))
660   CONTINUE
      DO 670 J = 1 ,JE
      DO 670 I = IB,IE
      VV(I,J)  = PP( VV(I,J) ) * MIN(ONE,BD(I,J),BU(I,J+1))
     &         - PN( VV(I,J) ) * MIN(ONE,BU(I,J),BD(I,J+1))
670   CONTINUE
C     B.C. BY THE SPHERICAL CYCLICITY
      DO 680 J = JB,JE
      VU(1 ,J) = VU(IE,J)
      VU(NX,J) = VU(IB,J)
680   CONTINUE
      DO 685 J = 1 ,JE
      VV(1 ,J) = VV(IE,J)
      VV(NX,J) = VV(IB,J)
685   CONTINUE
      ENDIF
700   CONTINUE
C
C     UPDATE THE PREDICTED FIELD
      DO 800 J = 1 ,NY
      DO 800 I = 1 ,NX
      Q(I,J,K) = X(I,J)
800   CONTINUE

!- check ---------------------------------------------------------
!
!      open (10,file='mpdats-s.out')
!
!      do j=1,ny
!        write(10,11) j,q(1,j,10),q(2,j,10)
!      enddo
!
!11    format(1x,i5,2e30.20)
!      close (10)
!
!      stop'mpdata'
!--------------------------------------------------------------


      RETURN
	END

#endif

