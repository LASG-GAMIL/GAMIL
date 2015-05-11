#include <misc.h>
#include <params.h>

#if (defined SPMD)


      SUBROUTINE TSPAS(Q,U,V,P,DSNP,DSSP,GC,DTDLT,DTDLN)

!     ==================
!     TWO-STEP SHAPE-PRESERVING ADVECTION SCHEME DEVELOPED
!     BY YU RUCONG
!     ==================

!     PERFORM 2-D ADVECTION IN GLOBAL SPHERICAL GEOMETRY
!     WITH THE UNIFORM LAT-LON C-GRID MESH
!
      use mpishorthand, only: mpicom

      IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

#include <commpi.h>
#include <mpif.h>

      character*50 filename

      integer isend11,irecv11
      integer isend12,irecv12
      integer isend22,irecv22
      integer istatus(mpi_status_size)

      real*8 workst1(nx*nl),workrb1(nx*nl)
      real*8 worksb1(nx*nl),workrt1(nx*nl)
      real*8 worksb2(nx*(nl+1)),workrt2(nx*(nl+1))
      integer ii


!
      REAL*8 Q(NX,NY,NL),U(NX,NY,NL),V(NX,NY,NL),P(NX,NY)
!
      REAL*8 QMIN(NX,NY,NL),QMAX(NX,NY,NL)
!
!!      REAL*8 HU(NX,NY),HV(NX,NY),CU(NX,NY),CV(NX,NY)
!!     _      ,UU(NX,NY),VV(NX,NY),H(NX,NY),A(NX,NY)
!!     _      ,QH(NX,NY),QHSTAR(NX,NY),USTAR(NX,NY),VSTAR(NX,NY)
!!      REAL*8 BETA(NX,NY),FX(NX,NY),FY(NX,NY),C(NX,2)

      REAL*8 HU(NX,NY),HV(NX,NY),CU(NX,NY),CV(NX,NY),H(NX,NY)

      real*8 UU(NX,NY,nl),VV(NX,NY,nl),A(NX,NY,nl)
      real*8 QH(NX,NY,nl),QHSTAR(NX,NY,nl),USTAR(NX,NY,nl),VSTAR(NX,NY,nl)
      REAL*8 BETA(NX,NY,nl),FX(NX,NY,nl),FY(NX,NY,nl),C(NX,2,nl)
!
!!    REAL*8 GAMA,XNP,XSP,CXSTAR,CYSTAR,CX,CY,TEMP1,TEMP2,TEMP3,TEMP4
      REAL*8 GAMA,XNP(nl),XSP(nl),CXSTAR,CYSTAR,CX,CY,TEMP1,TEMP2,TEMP3,TEMP4
!
      REAL*8  ZERO,HALF,FOURTH,EPSM
      REAL*8  DSNP,DSSP,DTDLN(NY),DTDLT(NY),GC(NY)
      DATA ZERO,HALF,FOURTH,EPSM/ 0.0D0,0.5D0,0.25D0,1.0D-80/
*     DATA ZERO,HALF,FOURTH,EPSM/ 0.0D0,0.5D0,0.25D0,1.0E-6/
      INTEGER I,J,K
      integer begj,endj

      begj = 2
      endj = ny-1
!
       DO K=1,NL
       DO J=begj,endj
       DO I=1,NX
         QMIN(I,J,K)=1.0E15
         QMAX(I,J,K)=-1.0E15
       ENDDO
       ENDDO
       ENDDO

!-H ---
      if (inc_pole) then
         DO I = 1 ,NX
            H(I,Jpole) = ZERO
         ENDDO
      endif

      do j = jbeg, jend
         DO I = IB,IE
            H(I,J) = GC(J) * P(I,J)
         ENDDO
            H(1 ,J)  = H(IE,J)
            H(NX,J)  = H(IB,J)
      enddo

      call mpi_isend( h(1,2), nx,mpi_double_precision,itop,1,mpicom,isend11,ierr )
      call mpi_irecv( h(1,ny),nx,mpi_double_precision,ibot,1,mpicom,irecv11,ierr )

!-HU,CU --

      if (inc_pole) then
         DO I = 1 ,NX
            HU(I,Jpole) = ZERO
            CU(I,Jpole) = ZERO
         ENDDO
      endif

      do j=jbeg, jend
          DO I = IB,IE
            HU(I,J)  = HALF*(H(I,J)+H(I-1,J))
          ENDDO
          HU(1 ,J) = HU(IE,J)
          HU(NX,J) = HU(IB,J)
          DO I = 1 ,NX
            CU(I,J) = DTDLN(J) * HU(I,J)
          ENDDO
      enddo

!-HV,CV --
      call mpi_wait( isend11,istatus,ierr )
      call mpi_wait( irecv11,istatus,ierr )

      DO J = 2 ,JEnd
        DO I = IB,IE
          HV(I,J)  = HALF*(H(I,J)+H(I,J+1))
        ENDDO
        HV(1 ,J) = HV(IE,J)
        HV(NX,J) = HV(IB,J)
        DO I = 1 ,NX
          CV(I,J) = DTDLT(J) * HV(I,J)
        ENDDO
      ENDDO
c----
      DO K=1,NL
      DO J=begj,endj
      DO I=1,NX
        QH(I,J,k)=Q(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      ii=1
      do k=1,nl
       do i=1,nx
          workst1(ii) = qh(i,2   ,k)
          worksb1(ii) = qh(i,ny-1,k)
          ii=ii+1
       enddo
      enddo

      call mpi_isend(workst1,nx*nl,mpi_double_precision,itop,2,mpicom,isend11,ierr)
      call mpi_irecv(workrb1,nx*nl,mpi_double_precision,ibot,2,mpicom,irecv11,ierr)
      call mpi_isend(worksb1,nx*nl,mpi_double_precision,ibot,3,mpicom,isend12,ierr)
      call mpi_irecv(workrt1,nx*nl,mpi_double_precision,itop,3,mpicom,irecv12,ierr)

C----
      DO K=1,NL
       DO J=JBeg,JEnd
        DO I=IB,IE
	  UU(I,J,k)=CU(I,J)*U(I,J,K)
        ENDDO
 	UU(1 ,J,k)=UU(IE,J,k)
        UU(NX,J,k)=UU(IB,J,k)
       ENDDO
      ENDDO
C----
      DO K=1,NL
       DO J=2,JEnd
        DO I=IB,IE
	  VV(I,J,k)=CV(I,J)*V(I,J,K)
        ENDDO
	VV(1 ,J,k)=VV(IE,J,k)
	VV(NX,J,k)=VV(IB,J,k)
       ENDDO
      ENDDO

      ii=1
       do i=1,nx
          worksb2(ii)=hv(i,ny-1)
          ii=ii+1
       enddo
      do k=1,nl
       do i=1,nx
          worksb2(ii)=vv(i,ny-1,k)
          ii=ii+1
       enddo
      enddo

      call mpi_isend(worksb2,nx*(nl+1),mpi_double_precision,ibot,4,mpicom,isend22,ierr)
      call mpi_irecv(workrt2,nx*(nl+1),mpi_double_precision,itop,4,mpicom,irecv22,ierr)

C----
      DO K=1,NL
       DO J=JBeg,JEnd
	DO I=IB,IE
	  FX(I,J,k)=HALF*UU(I,J,k)*(QH(I,J,k)+QH(I-1,J,k))
     $             -HALF*UU(I,J,k)*UU(I,J,k)*(QH(I,J,k)-QH(I-1,J,k))/HU(I,J)
	ENDDO
      	FX(NX,J,k)=FX(IB,J,k)
       ENDDO
      ENDDO
C----
      call mpi_wait( isend11,istatus,ierr )
      call mpi_wait( irecv11,istatus,ierr )
      call mpi_wait( isend12,istatus,ierr )
      call mpi_wait( irecv12,istatus,ierr )

      ii=1
      do k=1,nl
       do i=1,nx
          qh(i,ny,k) = workrb1(ii)
          qh(i,1 ,k) = workrt1(ii)
          ii=ii+1
       enddo
      enddo

      DO K=1,NL
       DO J=2,JEnd
        DO I=IB,IE
	  FY(I,J,k)=HALF*VV(I,J,k)*(QH(I,J+1,k)+QH(I,J,k))
     $             -HALF*VV(I,J,k)*VV(I,J,k)*(QH(I,J+1,k)-QH(I,J,k))/HV(I,J)
	ENDDO
       ENDDO
      ENDDO
C----
      ii=1
      do k=1,nl
       do i=1,nx
          worksb1(ii)=fy(i,ny-1,k)
          ii=ii+1
       enddo
      enddo

      call mpi_isend(worksb1,nx*nl,mpi_double_precision,ibot,5,mpicom,isend12,ierr)
      call mpi_irecv(workrt1,nx*nl,mpi_double_precision,itop,5,mpicom,irecv12,ierr)
!--
      call mpi_wait( isend22,istatus,ierr )
      call mpi_wait( irecv22,istatus,ierr )

      ii=1
       do i=1,nx
          hv(i,1)=workrt2(ii)
          ii=ii+1
       enddo
      do k=1,nl
       do i=1,nx
          vv(i,1,k)=workrt2(ii)
          ii=ii+1
       enddo
      enddo
!--
      call mpi_wait( isend12,istatus,ierr )
      call mpi_wait( irecv12,istatus,ierr )

      ii=1
      do k=1,nl
       do i=1,nx
          fy(i,1,k)=workrt1(ii)
          ii=ii+1
       enddo
      enddo

      DO K=1,NL
       DO J=JBeg,JEnd
	DO I=IB,IE
          TEMP1=ABS(UU(I,J,k)/HU(I,J))*(1-ABS(UU(I,J,k)/HU(I,J)))
          TEMP2=ABS(UU(I+1,J,k)/HU(I+1,J))*(1-ABS(UU(I+1,J,k)/HU(I+1,J)))
          TEMP3=ABS(VV(I,J-1,k)/HV(I,J-1))*(1-ABS(VV(I,J-1,k)/HV(I,J-1)))
	  TEMP4=ABS(VV(I,J,k)/HV(I,J))*(1-ABS(VV(I,J,k)/HV(I,J)))
	  GAMA=MAX(TEMP1,TEMP2,TEMP3,TEMP4)
          BETA(I,J,k)=2.0D0/(2.0D0-2.0D0*GAMA)
          QHSTAR(I,J,k)=QH(I,J,k)-BETA(I,J,k)*(FX(I+1,J,k)
     $                    -FX(I,J,k)+FY(I,J,k)-FY(I,J-1,k))/H(I,J)
        ENDDO
       ENDDO
      ENDDO
C------
      if (nprocs.ge.32) then

         if (myrank.eq.1) then
             call mpi_send( h(1,ny-1),nx,mpi_double_precision,0,61,mpicom,ierr)
             ii=1
             do k=1,nl
              do i=1,nx
                 worksb1(ii)=beta(i,ny-1,k)
                 ii=ii+1
              enddo
             enddo
             call mpi_isend(worksb1,nx*nl,mpi_double_precision,0,62,mpicom,isend12,ierr)
             call mpi_wait(isend12,istatus,ierr)
         endif

         if (myrank.eq.0 ) then
             call mpi_recv( h(1,1   ),nx,mpi_double_precision,1,61,mpicom,ierr)
             call mpi_irecv(workrt1,nx*nl,mpi_double_precision,1,62,mpicom,irecv12,ierr)
             call mpi_wait(irecv12,istatus,ierr)
             ii=1
             do k=1,nl
              do i=1,nx
                 beta(i,1,k)=workrt1(ii)
                 ii=ii+1
              enddo
             enddo
         endif

         if (myrank.eq.nprocs-2) then
             ii=1
             do k=1,nl
              do i=1,nx
                 workst1(ii)=beta(i,2,k)
                 ii=ii+1
              enddo
             enddo
             call mpi_isend(workst1,nx*nl,mpi_double_precision,nprocs-1,
     _                                 62,mpicom,isend12,ierr)
             call mpi_wait(isend12,istatus,ierr)
         endif

         if (myrank.eq.nprocs-1 ) then
             call mpi_irecv(workrb1,nx*nl,mpi_double_precision,nprocs-2,
     _                                 62,mpicom,irecv12,ierr)
             call mpi_wait(irecv12,istatus,ierr)
             ii=1
             do k=1,nl
              do i=1,nx
                 beta(i,ny,k)=workrb1(ii)
                 ii=ii+1
              enddo
             enddo
         endif

      endif
!--------------
      if (myrank.eq.nprocs-1) then
       DO K=1,NL
         DO 440 I = IB,IE
            C(I,1,k) =BETA(I,jpole+1,k)*FY(I,jpole,k)/H(I,jpole+1)
440      CONTINUE
            XNP(k) = ZERO
         DO 450 I = IB,IE
            XNP(k)   = XNP(k) + C(I,1,k)
450      CONTINUE
            XNP(k)   = QH(IB,jpole,k) - XNP(k)*DSNP
         DO 460 I = IB,IE
            QHSTAR(I,jpole,k)   = XNP(k)
460      CONTINUE
       ENDDO
      endif

      if (myrank.eq.0) then
       DO K=1,NL
         DO 441 I = IB,IE
            C(I,2,k) =BETA(I,Jpole-1,k)*FY(I,jpole-1,k)/H(I,jpole-1)
441      CONTINUE
            XSP(k)   = ZERO
         DO 451 I = IB,IE
            XSP(k)   = XSP(k) + C(I,2,k)
451      CONTINUE
            XSP(k)   = QH(IB,jpole,k)+ XSP(k)*DSSP
         DO 461 I = IB,IE
            QHSTAR(I,jpole,k)   = XSP(k)
461      CONTINUE
       ENDDO
      endif
C----
      DO K=1,NL
      DO J=begj,endj
	QHSTAR(1,J,k)=QHSTAR(IE,J,k)
	QHSTAR(NX,J,k)=QHSTAR(IB,J,k)
      ENDDO
      ENDDO
C----
      if (myrank.eq.nprocs-1) then
         do k=1,nl
          do i=ib,ie
	   QMIN(I,Jpole,K)=MIN(QH(I,Jpole,k),QH(I,Jpole+1,k),QMIN(I,Jpole,K))
	   QMAX(I,Jpole,K)=MAX(QH(I,Jpole,k),QH(I,Jpole+1,k),QMAX(I,Jpole,K))
           A(I,Jpole,k)=(QHSTAR(I,Jpole,k)-QMAX(I,Jpole,K))*(QHSTAR(I,Jpole,k)-QMIN(I,Jpole,K))
          enddo
	   A(1,Jpole,k)=A(IE,Jpole,k)
	   A(NX,Jpole,k)=A(IB,Jpole,k)
         enddo
      endif

      if (myrank.eq.0) then
         do k=1,nl
          do i=ib,ie
	    QMIN(I,Jpole,K)=MIN(QH(I,Jpole,k),QH(I,Jpole-1,k),QMIN(I,Jpole,K))
	    QMAX(I,Jpole,K)=MAX(QH(I,Jpole,k),QH(I,Jpole-1,k),QMAX(I,Jpole,K))
            A(I,Jpole,k)=(QHSTAR(I,Jpole,k)-QMAX(I,Jpole,K))*(QHSTAR(I,Jpole,k)-QMIN(I,Jpole,K))
          enddo
	    A(1,Jpole,k)=A(IE,Jpole,k)
	    A(NX,Jpole,k)=A(IB,Jpole,k)
         enddo
      endif

      DO K=1,NL
      DO J=jbeg,jend
        DO I=IB,IE
	  QMIN(I,J,K)=MIN(QH(I,J,k),QH(I+1,J,k),QH(I-1,J,k),
     $                    QH(I,J+1,k),QH(I,J-1,k),QMIN(I,J,K))
          QMAX(I,J,K)=MAX(QH(I,J,k),QH(I+1,J,k),QH(I-1,J,k),
     $                    QH(I,J+1,k),QH(I,J-1,k),QMAX(I,J,K))
          A(I,J,k)=(QHSTAR(I,J,k)-QMAX(I,J,K))*(QHSTAR(I,J,k)-QMIN(I,J,K))
        ENDDO
	A(1,J,k)=A(IE,J,k)
	A(NX,J,k)=A(IB,J,k)
      ENDDO
      ENDDO

      ii=1
      do k=1,nl
       do i=1,nx
          workst1(ii)=a(i,2,k)
          ii=ii+1
       enddo
      enddo

      call mpi_isend(workst1,nx*nl,mpi_double_precision,itop,7,mpicom,isend11,ierr)
      call mpi_irecv(workrb1,nx*nl,mpi_double_precision,ibot,7,mpicom,irecv11,ierr)

C----
      DO K=1,NL
      DO J=JBeg,JEnd
      DO I=IB,IE
	TEMP1=(ABS(A(I-1,J,k))+A(I-1,J,k))/(ABS(A(I-1,J,k))+EPSM)
	TEMP2=(ABS(A(I,J,k))+A(I,J,k))/(ABS(A(I,J,k))+EPSM)
	TEMP3=(ABS(A(I-1,J,k))+A(I-1,J,k))*(ABS(A(I,J,k))+A(I,J,k))
	TEMP4=ABS(A(I-1,J,k))*ABS(A(I,J,k))+EPSM
	CXSTAR=HALF*(TEMP1+TEMP2)-FOURTH*TEMP3/TEMP4
        CX=CXSTAR+(1-CXSTAR)*ABS(UU(I,J,k)/HU(I,J))
        USTAR(I,J,k)=CX*UU(I,J,k)
      ENDDO
      ENDDO
      ENDDO
C----
      call mpi_wait( isend11,istatus,ierr )
      call mpi_wait( irecv11,istatus,ierr )

      ii=1
      do k=1,nl
       do i=1,nx
          a(i,ny,k)=workrb1(ii)
          ii=ii+1
       enddo
      enddo

      DO K=1,NL
      DO J=2,JEnd
      DO I=IB,IE
	TEMP1=(ABS(A(I,J,k))+A(I,J,k))/(ABS(A(I,J,k))+EPSM)
	TEMP2=(ABS(A(I,J+1,k))+A(I,J+1,k))/(ABS(A(I,J+1,k))+EPSM)
	TEMP3=(ABS(A(I,J,k))+A(I,J,k))*(ABS(A(I,J+1,k))+A(I,J+1,k))
	TEMP4=ABS(A(I,J,k))*ABS(A(I,J+1,k))+EPSM
	CYSTAR=HALF*(TEMP1+TEMP2)-FOURTH*TEMP3/TEMP4
        CY=CYSTAR+(1-CYSTAR)*ABS(VV(I,J,k)/HV(I,J))
        VSTAR(I,J,k)=CY*VV(I,J,k)
      ENDDO
      ENDDO
      ENDDO
C----
      DO K=1,NL
      DO J=2,JEnd
      DO I=IB,IE
	  FY(I,J,k)=HALF*VV(I,J,k)*(QH(I,J+1,k)+QH(I,J,k))
     $             -HALF*ABS(VSTAR(I,J,k))*(QH(I,J+1,k)-QH(I,J,k))
      ENDDO
      ENDDO
      ENDDO

      ii=1
      do k=1,nl
       do i=1,nx
          worksb1(ii)=fy(i,ny-1,k)
          ii=ii+1
       enddo
      enddo

      call mpi_isend(worksb1,nx*nl,mpi_double_precision,ibot,8,mpicom,isend12,ierr)
      call mpi_irecv(workrt1,nx*nl,mpi_double_precision,itop,8,mpicom,irecv12,ierr)
C----
      DO K=1,NL
       DO J=JBeg,JEnd
	DO I=IB,IE
	  FX(I,J,k)=HALF*UU(I,J,k)*(QH(I,J,k)+QH(I-1,J,k))
     $           -HALF*ABS(USTAR(I,J,k))*(QH(I,J,k)-QH(I-1,J,k))
	ENDDO
	FX(NX,J,k)=FX(IB,J,k)
       ENDDO
      ENDDO
!--
      call mpi_wait( isend12,istatus,ierr )
      call mpi_wait( irecv12,istatus,ierr )

      ii=1
      do k=1,nl
       do i=1,nx
          fy(i,1,k)=workrt1(ii)
          ii=ii+1
       enddo
      enddo

C----
      DO K=1,NL
       DO J=JBeg,JEnd
        DO I=IB,IE
          QH(I,J,k)=QH(I,J,k)-(FX(I+1,J,k)-FX(I,J,k)
     $                        +FY(I,J,k)-FY(I,J-1,k))/H(I,J)
        ENDDO
       ENDDO
      ENDDO
C-----
       if (nprocs.ge.32) then
         if (myrank.eq.1) then
            ii=1
            do k=1,nl
             do i=1,nx
                worksb1(ii)=fy(i,ny-1,k)
                ii=ii+1
             enddo
            enddo
         call mpi_isend(worksb1,nx*nl,mpi_double_precision,0,9,mpicom,isend12,ierr)
         call mpi_wait(isend12,istatus,ierr)
         endif

         if (myrank.eq.0) then
         call mpi_irecv(workrt1,nx*nl,mpi_double_precision,1,9,mpicom,irecv12,ierr)
         call mpi_wait(irecv12,istatus,ierr)
            ii=1
            do k=1,nl
             do i=1,nx
                fy(i,1,k)=workrt1(ii)
                ii=ii+1
             enddo
            enddo
         endif

       endif
!---
      if (myrank.eq.nprocs-1) then
      DO K=1,NL
      DO 540 I = IB,IE
      C(I,1,k) = FY(I,jpole,k) / H(I,jpole+1)
540   CONTINUE
      XNP(k)   = ZERO
      DO 550 I = IB,IE
      XNP(k)   = XNP(k) + C(I,1,k)
550   CONTINUE
      XNP(k)   = QH(IB,jpole,k) - XNP(k)*DSNP
      DO 560 I = IB,IE
      QH(I,jpole,k)   = XNP(k)
560   CONTINUE
      ENDDO
      endif

      if (myrank.eq.0) then
      DO K=1,NL
      DO 541 I = IB,IE
      C(I,2,k) = FY(I,jpole-1,k) / H(I,jpole-1)
541   CONTINUE
      XSP(k)   = ZERO
      DO 551 I = IB,IE
      XSP(k)   = XSP(k)     + C(I,2,k)
551   CONTINUE
      XSP(k)   = QH(IB,jpole,k) + XSP(k)*DSSP
      DO 561 I = IB,IE
      QH(I,jpole,k)   = XSP(k)
561   CONTINUE
      ENDDO
      endif
C----
      DO K=1,NL
      DO J=begj,endj
        QH(1,J,k)=QH(IE,J,k)
        QH(NX,J,k)=QH(IB,J,k)
      ENDDO
      ENDDO
C----

      ii=1
      do k=1,nl
       do i=1,nx
          workst1(ii)=qh(i,2   ,k)
          worksb1(ii)=qh(i,ny-1,k)
          ii=ii+1
       enddo
      enddo

        call mpi_isend(workst1,nx*nl,mpi_double_precision,itop,10,mpicom,isend11,ierr)
        call mpi_isend(worksb1,nx*nl,mpi_double_precision,ibot,10,mpicom,isend12,ierr)
        call mpi_irecv(workrb1,nx*nl,mpi_double_precision,ibot,10,mpicom,irecv11,ierr)
        call mpi_irecv(workrt1,nx*nl,mpi_double_precision,itop,10,mpicom,irecv12,ierr)
!--
      DO K=1,NL
      DO J=begj,endj
      DO I=1,NX
	Q(I,J,K)=QH(I,J,k)
      ENDDO
      ENDDO
      ENDDO

!---
        call mpi_wait(isend11,istatus,ierr)
        call mpi_wait(isend12,istatus,ierr)
        call mpi_wait(irecv11,istatus,ierr)
        call mpi_wait(irecv12,istatus,ierr)

      ii=1
      do k=1,nl
       do i=1,nx
          qh(i,ny,k)=workrb1(ii)
          qh(i,1 ,k)=workrt1(ii)
          ii=ii+1
       enddo
      enddo

C
C*********BEGIN CROSS ITEM************
      DO K=1,NL
      DO J=JBeg,JEnd
      DO I=IB,IE
	TEMP1=0.5*(UU(I,J,k)+UU(I+1,J,k))
	TEMP2=0.5*(VV(I,J-1,k)+VV(I,J,k))
	TEMP3=0.25*(QH(I+1,J+1,k)+QH(I-1,J-1,k)-QH(I+1,J-1,k)-QH(I-1,J+1,k))
	Q(I,J,K)=QH(I,J,k)+TEMP1*TEMP2*TEMP3/(H(I,J)*H(I,J))
      ENDDO
      ENDDO
      ENDDO
C*********END CROSS ITEM*************
C
      DO K=1,NL
      DO J=begj,endj
        Q(1,J,K)=Q(IE,J,K)
        Q(NX,J,K)=Q(IB,J,K)
      ENDDO
      ENDDO

!- check ---------------------------------------------------------
!
!      write(filename,14) 'tspas-p-',myrank,'.out'
!14    format(a8,i1,a4)
!
!      open (10,file=trim(filename))
!
!      write(10,*) 'qqqqqqqqqqqqqqqqqqqqq'
!      do j=1,ny
!        write(10,11) j,q(1,j,10),q(2,j,10)
!      enddo
!
!11    format(1x,i5,2e30.20)
!      close (10)
!
!      call mpi_finalize(j)
!      stop'tspas'
!--------------------------------------------------------------

C
C
      RETURN
	END


#else


      SUBROUTINE TSPAS(Q,U,V,P,DSNP,DSSP,GC,DTDLT,DTDLN)
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
      REAL*8 Q(NX,NY,NL),U(NX,NY,NL),V(NX,NY,NL),P(NX,NY)
!
      REAL*8 QMIN(NX,NY,NL),QMAX(NX,NY,NL)
!
!!      REAL*8 HU(NX,NY),HV(NX,NY),CU(NX,NY),CV(NX,NY)
!!     _      ,UU(NX,NY),VV(NX,NY),H(NX,NY),A(NX,NY)
!!     _      ,QH(NX,NY),QHSTAR(NX,NY),USTAR(NX,NY),VSTAR(NX,NY)
!!      REAL*8 BETA(NX,NY),FX(NX,NY),FY(NX,NY),C(NX,2)

      REAL*8 HU(NX,NY),HV(NX,NY),CU(NX,NY),CV(NX,NY),H(NX,NY)

      real*8 UU(NX,NY,nl),VV(NX,NY,nl),A(NX,NY,nl)
      real*8 QH(NX,NY,nl),QHSTAR(NX,NY,nl),USTAR(NX,NY,nl),VSTAR(NX,NY,nl)
      REAL*8 BETA(NX,NY,nl),FX(NX,NY,nl),FY(NX,NY,nl),C(NX,2,nl)
!
!!    REAL*8 GAMA,XNP,XSP,CXSTAR,CYSTAR,CX,CY,TEMP1,TEMP2,TEMP3,TEMP4
      REAL*8 GAMA,XNP(nl),XSP(nl),CXSTAR,CYSTAR,CX,CY,TEMP1,TEMP2,TEMP3,TEMP4
!
      REAL*8  ZERO,HALF,FOURTH,EPSM
      REAL*8  DSNP,DSSP,DTDLN(NY),DTDLT(NY),GC(NY)
      DATA ZERO,HALF,FOURTH,EPSM/ 0.0D0,0.5D0,0.25D0,1.0D-80/
*     DATA ZERO,HALF,FOURTH,EPSM/ 0.0D0,0.5D0,0.25D0,1.0E-6/
      INTEGER I,J,K
!
       DO K=1,NL
       DO J=1,NY
       DO I=1,NX
         QMIN(I,J,K)=1.0E15
         QMAX(I,J,K)=-1.0E15
       ENDDO
       ENDDO
       ENDDO
C
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
            HU(I,J)  = HALF*(H(I,J)+H(I-1,J))
          ENDDO
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
C
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
      DO K=1,NL
       DO J=JB,JE
        DO I=IB,IE
	  UU(I,J,k)=CU(I,J)*U(I,J,K)
        ENDDO
 	UU(1 ,J,k)=UU(IE,J,k)
        UU(NX,J,k)=UU(IB,J,k)
       ENDDO
      ENDDO
C
      DO K=1,NL
       DO J=1,JE
        DO I=IB,IE
	  VV(I,J,k)=CV(I,J)*V(I,J,K)
        ENDDO
	VV(1 ,J,k)=VV(IE,J,k)
	VV(NX,J,k)=VV(IB,J,k)
       ENDDO
      ENDDO
C
      DO K=1,NL
      DO J=1,NY
      DO I=1,NX
        QH(I,J,k)=Q(I,J,K)
      ENDDO
      ENDDO
      ENDDO
C
      DO K=1,NL
       DO J=JB,JE
	DO I=IB,IE
	  FX(I,J,k)=HALF*UU(I,J,k)*(QH(I,J,k)+QH(I-1,J,k))
     $             -HALF*UU(I,J,k)*UU(I,J,k)*(QH(I,J,k)-QH(I-1,J,k))/HU(I,J)
	ENDDO
      	FX(NX,J,k)=FX(IB,J,k)
       ENDDO
      ENDDO
C
      DO K=1,NL
       DO J=1,JE
        DO I=IB,IE
	  FY(I,J,k)=HALF*VV(I,J,k)*(QH(I,J+1,k)+QH(I,J,k))
     $             -HALF*VV(I,J,k)*VV(I,J,k)*(QH(I,J+1,k)-QH(I,J,k))/HV(I,J)
	ENDDO
       ENDDO
      ENDDO
C
      DO K=1,NL
       DO J=JB,JE
	DO I=IB,IE
          TEMP1=ABS(UU(I,J,k)/HU(I,J))*(1-ABS(UU(I,J,k)/HU(I,J)))
          TEMP2=ABS(UU(I+1,J,k)/HU(I+1,J))*(1-ABS(UU(I+1,J,k)/HU(I+1,J)))
          TEMP3=ABS(VV(I,J-1,k)/HV(I,J-1))*(1-ABS(VV(I,J-1,k)/HV(I,J-1)))
	  TEMP4=ABS(VV(I,J,k)/HV(I,J))*(1-ABS(VV(I,J,k)/HV(I,J)))
	  GAMA=MAX(TEMP1,TEMP2,TEMP3,TEMP4)
          BETA(I,J,k)=2.0D0/(2.0D0-2.0D0*GAMA)
          QHSTAR(I,J,k)=QH(I,J,k)-BETA(I,J,k)*(FX(I+1,J,k)
     $                    -FX(I,J,k)+FY(I,J,k)-FY(I,J-1,k))/H(I,J)
        ENDDO
       ENDDO
      ENDDO
C
      DO K=1,NL
      DO 440 I = IB,IE
      C(I,1,k) =BETA(I,2,k)*FY(I,1,k)/H(I,2)
      C(I,2,k) =BETA(I,JE,k)*FY(I,JE,k)/H(I,JE)
440   CONTINUE
      XNP(k)   = ZERO
      XSP(k)   = ZERO
      DO 450 I = IB,IE
      XNP(k)   = XNP(k)     + C(I,1,k)
      XSP(k)   = XSP(k)     + C(I,2,k)
450   CONTINUE
      XNP(k)   = QH(IB,1,k) - XNP(k)*DSNP
      XSP(k)   = QH(IB,NY,k)+ XSP(k)*DSSP
      DO 460 I = IB,IE
      QHSTAR(I,1,k)   = XNP(k)
      QHSTAR(I,NY,k)   = XSP(k)
460   CONTINUE
      ENDDO
C
      DO K=1,NL
      DO J=1,NY
	QHSTAR(1,J,k)=QHSTAR(IE,J,k)
	QHSTAR(NX,J,k)=QHSTAR(IB,J,k)
      ENDDO
      ENDDO
C
      DO K=1,NL
      DO J=1,NY
        DO I=IB,IE
	IF(J.EQ.1) THEN
	  QMIN(I,J,K)=MIN(QH(I,J,k),QH(I,J+1,k),QMIN(I,J,K))
	  QMAX(I,J,K)=MAX(QH(I,J,k),QH(I,J+1,k),QMAX(I,J,K))
	ELSE IF(J.EQ.NY) THEN
	  QMIN(I,J,K)=MIN(QH(I,J,k),QH(I,J-1,k),QMIN(I,J,K))
	  QMAX(I,J,K)=MAX(QH(I,J,k),QH(I,J-1,k),QMAX(I,J,K))
	ELSE
	  QMIN(I,J,K)=MIN(QH(I,J,k),QH(I+1,J,k),QH(I-1,J,k),
     $                    QH(I,J+1,k),QH(I,J-1,k),QMIN(I,J,K))
          QMAX(I,J,K)=MAX(QH(I,J,k),QH(I+1,J,k),QH(I-1,J,k),
     $                    QH(I,J+1,k),QH(I,J-1,k),QMAX(I,J,K))
	ENDIF
        A(I,J,k)=(QHSTAR(I,J,k)-QMAX(I,J,K))*(QHSTAR(I,J,k)-QMIN(I,J,K))
        ENDDO
	A(1,J,k)=A(IE,J,k)
	A(NX,J,k)=A(IB,J,k)
      ENDDO
      ENDDO
C
      DO K=1,NL
      DO J=JB,JE
      DO I=IB,IE
	TEMP1=(ABS(A(I-1,J,k))+A(I-1,J,k))/(ABS(A(I-1,J,k))+EPSM)
	TEMP2=(ABS(A(I,J,k))+A(I,J,k))/(ABS(A(I,J,k))+EPSM)
	TEMP3=(ABS(A(I-1,J,k))+A(I-1,J,k))*(ABS(A(I,J,k))+A(I,J,k))
	TEMP4=ABS(A(I-1,J,k))*ABS(A(I,J,k))+EPSM
	CXSTAR=HALF*(TEMP1+TEMP2)-FOURTH*TEMP3/TEMP4
        CX=CXSTAR+(1-CXSTAR)*ABS(UU(I,J,k)/HU(I,J))
        USTAR(I,J,k)=CX*UU(I,J,k)
      ENDDO
      ENDDO
      ENDDO
C
      DO K=1,NL
      DO J=1,JE
      DO I=IB,IE
	TEMP1=(ABS(A(I,J,k))+A(I,J,k))/(ABS(A(I,J,k))+EPSM)
	TEMP2=(ABS(A(I,J+1,k))+A(I,J+1,k))/(ABS(A(I,J+1,k))+EPSM)
	TEMP3=(ABS(A(I,J,k))+A(I,J,k))*(ABS(A(I,J+1,k))+A(I,J+1,k))
	TEMP4=ABS(A(I,J,k))*ABS(A(I,J+1,k))+EPSM
	CYSTAR=HALF*(TEMP1+TEMP2)-FOURTH*TEMP3/TEMP4
        CY=CYSTAR+(1-CYSTAR)*ABS(VV(I,J,k)/HV(I,J))
        VSTAR(I,J,k)=CY*VV(I,J,k)
      ENDDO
      ENDDO
      ENDDO
C
      DO K=1,NL
       DO J=JB,JE
	DO I=IB,IE
	  FX(I,J,k)=HALF*UU(I,J,k)*(QH(I,J,k)+QH(I-1,J,k))
     $           -HALF*ABS(USTAR(I,J,k))*(QH(I,J,k)-QH(I-1,J,k))
	ENDDO
	FX(NX,J,k)=FX(IB,J,k)
       ENDDO
      ENDDO
C
      DO K=1,NL
      DO I=IB,IE
	DO J=1,JE
	  FY(I,J,k)=HALF*VV(I,J,k)*(QH(I,J+1,k)+QH(I,J,k))
     $ -HALF*ABS(VSTAR(I,J,k))*(QH(I,J+1,k)-QH(I,J,k))
        ENDDO
      ENDDO
      ENDDO

C
      DO K=1,NL
       DO J=JB,JE
        DO I=IB,IE
          QH(I,J,k)=QH(I,J,k)-(FX(I+1,J,k)-FX(I,J,k)
     $                        +FY(I,J,k)-FY(I,J-1,k))/H(I,J)
        ENDDO
       ENDDO
      ENDDO
C
      DO K=1,NL
      DO 540 I = IB,IE
      C(I,1,k)   = FY(I,1,k) / H(I,2)
      C(I,2,k)   = FY(I,JE,k) / H(I,JE)
540   CONTINUE
      XNP(k)      = ZERO
      XSP(k)      = ZERO
      DO 550 I = IB,IE
      XNP(k)      = XNP(k)     + C(I,1,k)
      XSP(k)      = XSP(k)     + C(I,2,k)
550   CONTINUE
      XNP(k)      = QH(IB,1,k) - XNP(k)*DSNP
      XSP(k)      = QH(IB,NY,k) + XSP(k)*DSSP
      DO 560 I = IB,IE
      QH(I,1,k)   = XNP(k)
      QH(I,NY,k)   = XSP(k)
560   CONTINUE
      ENDDO
C
      DO K=1,NL
      DO J=1,NY
        QH(1,J,k)=QH(IE,J,k)
        QH(NX,J,k)=QH(IB,J,k)
      ENDDO
      ENDDO
C
      DO K=1,NL
      DO I=1,NX
      DO J=1,NY
	Q(I,J,K)=QH(I,J,k)
      ENDDO
      ENDDO
      ENDDO
C
C*********BEGIN CROSS ITEM************
      DO K=1,NL
      DO I=IB,IE
      DO J=JB,JE
	TEMP1=0.5*(UU(I,J,k)+UU(I+1,J,k))
	TEMP2=0.5*(VV(I,J-1,k)+VV(I,J,k))
	TEMP3=0.25*(QH(I+1,J+1,k)+QH(I-1,J-1,k)-QH(I+1,J-1,k)-QH(I-1,J+1,k))
	Q(I,J,K)=QH(I,J,k)+TEMP1*TEMP2*TEMP3/(H(I,J)*H(I,J))
      ENDDO
      ENDDO
      ENDDO
C*********END CROSS ITEM*************
C
      DO K=1,NL
      DO J=1,NY
        Q(1,J,K)=Q(IE,J,K)
        Q(NX,J,K)=Q(IB,J,K)
      ENDDO
      ENDDO

!- check ---------------------------------------------------------
!
!      open (10,file='tspas-s.out')
!
!      write(10,*) 'qqqqqqqqqqqqqqqqqqqqq'
!      do j=1,ny
!        write(10,11) j,q(1,j,10),q(2,j,10)
!      enddo
!
!11    format(1x,i5,2e30.20)
!      close (10)
!
!      stop'tspas'
!--------------------------------------------------------------

C
C
      RETURN
	END

#endif
