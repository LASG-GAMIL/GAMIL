#include <misc.h>
#include <params.h>

!!(2003.11.10)
!!---------------------


      SUBROUTINE DIFUVT(U,US,V,VS,WS,P,PS,PLY,DPS,TT,H,HPS,CB,DCB
     _                ,SIGL,DSIG,DY,OUX,OVX,OUY,OVY,SINV,FF
     _                ,CUR,DLT1,DLT2,MP1,MP2,MP3,MM1,MM2,MM3,WTGV
     _                ,DU,DV,DTT,SU,SV,ST,FBC)

#if ( defined SPMD )
      use mpishorthand, only: mpicom
#endif

      IMPLICIT NONE
!
!	This subroutine is to calculate the tendency of the wind and the temperature:
!     DU, DV, DTT

#include <PARADYN>

#if (defined SPMD)
#include <commpi.h>
#include <mpif.h>

      character*50 filename

      integer num1,num2
      parameter(num1=nx*(nl*4+1),num2=nx*(8*nl+2*nz+1) )

      real*8  worksb( num1 )
      real*8  workrt( num1 )
      real*8  workst( num2 )
      real*8  workrb( num2 )

      integer isend1,irecv1
      integer isend2,irecv2
      integer istatus(mpi_status_size)
      integer ii
#else
      integer jpole
#endif
!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	  REAL*8
     _       U  (NX,NY,NL)  !  input variable, U  = u*sqrt(Ps)
     _	  ,US(NX,NY,NL )  !  input variable, US = u, zonal wind
     _      ,V  (NX,NY,NL)  !	 input variable, V  = v*sqrt(Ps)
     _      ,VS(NX,NY,NL )  !  input variable, VS = v, meridional wind
     _      ,WS(NX,NY,NZ )  !  input variable, WS = w, vertical velocity
     _      ,P  (NX,NY   )  !  input variable, P  = sqrt(Ps)
     _      ,PS (NX,NY   )  !	 input variable, Surface pressure
     _      ,PLY(NX,NY,NL)  !  input variable, PLY=p, Pressure in Sigma Layer
     _      ,DPS(NX,NY   )  !  input variable,
!                              the tendency of the surface pressure
     _      ,TT (NX,NY,NL)  !  input variable, TT=R*T'*Sqrt(Ps)/CB,
!                              where T'=T-TB, T is the temperatur,
!                              TBB	is Temperature of the standard atmosphere
     _	  ,H (NX,NY,NZ )  !  input variable,
     _	  ,HH(NX,NY,NZ )  !  input variable,
!  			     H=gz-HBB, gz is the geopotential height,
!                              HBB is the geopotential height of the standard atmopshere
     _      ,CB (NX,NY,NL)  !  input variable, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
!                              where, KK=R/Cp, R is a constant
     _      ,DCB(NX,NY,NL)  !  input variable,
     _      ,SIGL(NL     )  !  input constant, the vertical layers
     _      ,DSIG(NL     )  !  input constant, the vertical stepsizes
     _      ,DY             !  input constant,
!                              the horizontal stepsize in meridional direction
     _      ,SINV(NY)       !  input constant, sin(theta) at half grid j+1/2
     _      ,OUX(NY)        !  input constant, OUX=1/(RAD*SINU*DX*MDJ)
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
     _      ,OVX(NY)        !  input constant, OUX=1/(RAD*SINV*DX*MDJ)
     _      ,OUY(NY)        !  input constant, OUY=1/(RAD*SINU*DY*WTGU)
!                              where WTGU is the weighting at the integer grid j
     _      ,OVY(NY)        !  input constant, OUY=1/(RAD*SINV*DY*WTGU)
     _      ,WTGV(NY)       !	 input constant,
!                              the weighting at the half grid j+1/2,
     _	  ,DLT1           !  input constant
     _      ,DLT2           !	 input constant
     _      ,FF(NY)         !	 input constant
     _      ,CUR(NY)		  !	 input constant
!
	INTEGER
     _       MM1(NX,NY)     !	 input constant
     _      ,MP1(NX,NY)     !	 input constant
     _      ,MM2(NX,NY)     !	 input constant
     _      ,MP2(NX,NY)     !	 input constant
     _      ,MM3(NX,NY)     !	 input constant
     _      ,MP3(NX,NY)     !	 input constant
!
      REAL*8
     _       DU(NX,NY,NL)	  !  output variables
     _      ,DV(NX,NY,NL)   !  output variables
     _      ,DTT(NX,NY,NL)  !  output variables
     _      ,SU(NX,NY,NL)   !  input variables
     _      ,SV(NX,NY,NL)   !  input variables
     _      ,ST(NX,NY,NL)   !  input variables
!
      REAL*8  TTZ(NX,NY,NZ+1)
      REAL*8  UZ(NX,NY,NZ+1),VZ(NX,NY,NZ+1)
      INTEGER I,J,K
      INTEGER IM1,IP2,IP1,IM2,I1,J0,IM3,IP3
      REAL*8  OY1,R22,R21,R12,DYP0,DXP1,DXP0,DTTP,TO3P
      REAL*8  R11,OUX2,DYP2,O2,O1,TL3,TO3,TO2,TO1
      REAL*8  TL2,TL1,OPK0,OZ2
      REAL*8  OY4,RI0,TO1P,TL3P,TO2P
      REAL*8  O2P,TLP,O1P,WKP
      REAL*8  RV1,RV2,OVX4,OVY4,PDXP1,PDXP2,R14
      REAL*8  OUX4,R10,R20,R24,UL2,UL3
      REAL*8  H0,UL1,PX2,FSV,PX1,AXP0
      REAL*8  OPK1,PDYP2,AYP0,PDYP1
      REAL*8  VL2,VL3,OPK3,VL1,PY1
      REAL*8  FSU,PY2,FFYY,OZ4,OPK2
      REAL*8  FS0,FFXX,H1,OYY4,HPS(NX,NY)
      REAL*8  FBC(NX,NY,NL),TTV(NX,NY,NL)
      integer begj,endj

#if (defined SPMD)
      begj = 2
      endj = ny-1
#else
      begj = 1
      endj = ny
#endif

#if (defined SPMD)

      ii=1
      do i=1,nx
         worksb(ii)=p(i,ny-1)
         ii=ii+1
      enddo
!--1--
      do k=1,nl
       do i=1,nx
         worksb(ii)=u(i,ny-1,k)
         ii=ii+1
       enddo
      enddo
!--2--
      do k=1,nl
       do i=1,nx
         worksb(ii)=v(i,ny-1,k)
         ii=ii+1
       enddo
      enddo
!--3--
      do k=1,nl
       do i=1,nx
         worksb(ii)=tt(i,ny-1,k)
         ii=ii+1
       enddo
      enddo
!--4--
      do k=1,nl
       do i=1,nx
         worksb(ii)=vs(i,ny-1,k)
         ii=ii+1
       enddo
      enddo

      call mpi_isend( worksb,num1,mpi_double_precision,ibot,1,mpicom,isend1,ierr)
      call mpi_irecv( workrt,num1,mpi_double_precision,itop,1,mpicom,irecv1,ierr)
#endif

C     AS FOLLOWS, THE IMPROVED LEAP-FROG AND THE REGENERATION OF
C     VARIABLES WILL BE FINISHED IN THE SAME CYCLES.
C
      I1=NX-1
      J0=NY-1
      RI0=DFLOAT(NX-2)
      OY1=1.0D0/(RAD*DY)
      OY4=0.25D0*OY1

C
C     CALCULATING DU/DT, DV/DT AND DTT/DT.
C
C     FORMING THE FACTORS OUT OF ALL THE CYCLES AND CLEANING SOME
C     SPECIAL ARRAIES.
C
!     (DO J=1,NY)
      DO J=begj,endj
      DO I=1,NX
      DO K=1,NL
	 TTV(I,J,K)=TT(I,J,K)*FBC(I,J,K)
      END DO
      END DO
      END DO
c
C
        DO K=1,NZ
        DO J=begj,endj
        DO I=1,NX
           HH(I,J,K)=H(I,J,K)-HPS(I,J)
        ENDDO
        ENDDO
        ENDDO

!     (DO J=1,NY)
      do j=begj,endj
        DO K = 1,NZ+1
          IF(K.EQ.1.OR.K.EQ.NZ+1) THEN
            DO I = 1,NX
              TTZ(I,J,K) = 0.0D0
              UZ(I,J,K) = 0.0D0
              VZ(I,J,K) = 0.0D0
            ENDDO
          ELSE
            DO I = 1,NX
              TTZ(I,J,K) = TT(I,J,K-1)
              UZ(I,J,K) = U(I,J,K-1)
              VZ(I,J,K) = V(I,J,K-1)
            ENDDO
          ENDIF
        ENDDO
      enddo

#if (defined SPMD)
      ii=1
      do i=1,nx
         workst(ii)=p(i,2)
         ii=ii+1
      enddo
!--1--
      do k=1,nz
       do i=1,nx
         workst(ii)=h(i,2,k)
         ii=ii+1
       enddo
      enddo
!--2--
      do k=1,nz
       do i=1,nx
         workst(ii)=ws(i,2,k)
         ii=ii+1
       enddo
      enddo
!--1--
      do k=1,nl
       do i=1,nx
         workst(ii)=cb(i,2,k)
         ii=ii+1
       enddo
      enddo
!--2--
      do k=1,nl
       do i=1,nx
         workst(ii)=ply(i,2,k)
         ii=ii+1
       enddo
      enddo
!--3--
      do k=1,nl
       do i=1,nx
         workst(ii)=u(i,2,k)
         ii=ii+1
       enddo
      enddo
!--4--
      do k=1,nl
       do i=1,nx
         workst(ii)=v(i,2,k)
         ii=ii+1
       enddo
      enddo
!--5--
      do k=1,nl
       do i=1,nx
         workst(ii)=tt(i,2,k)
         ii=ii+1
       enddo
      enddo
!--6--
      do k=1,nl
       do i=1,nx
         workst(ii)=ttv(i,2,k)
         ii=ii+1
       enddo
      enddo
!--7--
      do k=1,nl
       do i=1,nx
         workst(ii)=us(i,2,k)
         ii=ii+1
       enddo
      enddo
!--8--
      do k=1,nl
       do i=1,nx
         workst(ii)=vs(i,2,k)
         ii=ii+1
       enddo
      enddo

      call mpi_isend( workst,num2,mpi_double_precision,itop,2,mpicom,isend2,ierr)
      call mpi_irecv( workrb,num2,mpi_double_precision,ibot,2,mpicom,irecv2,ierr)

!----------------------------------------------unpacking dataset 1----
      call mpi_wait( isend1,istatus,ierr )
      call mpi_wait( irecv1,istatus,ierr )

      ii=1
      do i=1,nx
         p(i,1)=workrt(ii)
         ii=ii+1
      enddo
!--1--
      do k=1,nl
       do i=1,nx
         u(i,1,k)=workrt(ii)
         ii=ii+1
       enddo
      enddo
!--2--
      do k=1,nl
       do i=1,nx
         v(i,1,k)=workrt(ii)
         ii=ii+1
       enddo
      enddo
!--3--
      do k=1,nl
       do i=1,nx
         tt(i,1,k)=workrt(ii)
         ii=ii+1
       enddo
      enddo
!--4--
      do k=1,nl
       do i=1,nx
         vs(i,1,k)=workrt(ii)
         ii=ii+1
       enddo
      enddo

!----------------------------------------------unpacking dataset 2----

      call mpi_wait( isend2,istatus,ierr )
      call mpi_wait( irecv2,istatus,ierr )

      ii=1
      do i=1,nx
         p(i,ny)=workrb(ii)
         ii=ii+1
      enddo
!--1--
      do k=1,nz
       do i=1,nx
         h(i,ny,k)=workrb(ii)
         ii=ii+1
       enddo
      enddo
!--2--
      do k=1,nz
       do i=1,nx
         ws(i,ny,k)=workrb(ii)
         ii=ii+1
       enddo
      enddo
!--1--
      do k=1,nl
       do i=1,nx
         cb(i,ny,k)=workrb(ii)
         ii=ii+1
       enddo
      enddo
!--2--
      do k=1,nl
       do i=1,nx
         ply(i,ny,k)=workrb(ii)
         ii=ii+1
       enddo
      enddo
!--3--
      do k=1,nl
       do i=1,nx
         u(i,ny,k)=workrb(ii)
         ii=ii+1
       enddo
      enddo
!--4--
      do k=1,nl
       do i=1,nx
         v(i,ny,k)=workrb(ii)
         ii=ii+1
       enddo
      enddo
!--5--
      do k=1,nl
       do i=1,nx
         tt(i,ny,k)=workrb(ii)
         ii=ii+1
       enddo
      enddo
!--6--
      do k=1,nl
       do i=1,nx
         ttv(i,ny,k)=workrb(ii)
         ii=ii+1
       enddo
      enddo
!--7--
      do k=1,nl
       do i=1,nx
         us(i,ny,k)=workrb(ii)
         ii=ii+1
       enddo
      enddo
!--8--
      do k=1,nl
       do i=1,nx
         vs(i,ny,k)=workrb(ii)
         ii=ii+1
       enddo
      enddo

#endif


C
C     CALCULATING DTT/DT AT POLES.
C
#if (!defined SPMD)
       jpole =1
#else
       if (myrank.eq.(nprocs-1)) then
#endif
!!        (IF(J.EQ.1) THEN)   ! for the north pole

          DO K=1,NL
            OZ2=0.5D0/DSIG(K)
            OZ4=0.5D0*OZ2
C
C       CALCULATING ADVECTION OF TT AND VERTICAL MOTION OVER PRESSURE
C       AT POLES.
C
            TLP=0.0D0
            O2P=0.0D0
            DO I=2,I1
              TLP=TLP+VS(I,jpole,K)*TT(I,jpole+1,K)
              O2P=O2P+(P(I,jpole+1)-P(1,jpole))*V(I,jpole,K)
            ENDDO
            WKP=2.0D0/RI0*OY1*WTGV(jpole)
            TLP=WKP*TLP
            O2P=2.0D0*WKP*O2P*SIGL(K)/PLY(1,jpole,K)
            O1P=(0.5D0*PS(1,jpole)*(WS(1,jpole,K+1)+WS(1,jpole,K))
     &        +DPS(1,jpole)*SIGL(K))/PLY(1,jpole,K)
            TL3P=OZ2*(WS(1,jpole,K+1)*TTZ(1,jpole,K+2)
     &               -WS(1,jpole,K)*TTZ(1,jpole,K))
C
            TO1P=CB(1,jpole,K)*P(1,jpole)*O1P
            TO2P=CB(1,jpole,K)*P(1,jpole)*O2P
            TO3P=(DLT1*CAPA-DLT2*DCB(1,jpole,K))*TTV(1,jpole,K)*(O1P+O2P)
C
C     CALCULATING DTT/DT AND FORMING BOTH POLAR BOUNDARIES.
            DTTP=-TLP-TL3P+TO1P+TO2P+TO3P
            DO I=1,NX
              DTT(I,jpole,K)=DTTP+ST(I,jpole,K)
              DU(I,jpole,K)=SU(I,jpole,K)
            ENDDO
            DO I=2,I1
              IP1=MP1(I,jpole)
              IM1=MM1(I,jpole)
              IP2=MP2(I,jpole)
              IM2=MM2(I,jpole)
C
              DYP0=P(I,jpole+1)-P(I,jpole)
              PDYP1=P(I,jpole+1)*DYP0
              PDYP2=P(I,jpole)*DYP0
              AYP0=0.25D0*(P(I,jpole)+P(I,jpole+1))
C
              OPK0=1.0D0/PLY(I,jpole,K)
              OPK3=1.0D0/PLY(I,jpole+1,K)
C
              VL1=0.25D0*OVX(jpole)*(US(IP2,jpole+1,K)*V(IP1,jpole,K)
     &                          -US(IM2,jpole+1,K)*V(IM1,jpole,K))
              VL2=0.25D0*OVY(jpole)*(SINV(jpole+1)*VS(I,jpole+1,K)
     &                          +SINV(jpole)*VS(I,jpole,K))*V(I,jpole+1,K)
              VL3=OZ4*((WS(I,jpole+1,K+1)+WS(I,jpole,K+1))*VZ(I,jpole,K+2)
     &                -(WS(I,jpole+1,K)+WS(I,jpole,K))*VZ(I,jpole,K))
C
              FSU=0.25D0*((FF(jpole+1)+CUR(jpole+1)*US(I,jpole+1,K))*U(I,jpole+1,K)
     &                   +(FF(jpole+1)+CUR(jpole+1)*US(I+1,jpole+1,K))*U(I+1,jpole+1,K))
              PY1=OY1*AYP0*(H(I,jpole+1,K+1)+H(I,jpole+1,K)-H(I,jpole,K+1)-H(I,jpole,K))
              PY2=OY1*(PDYP1*TTV(I,jpole+1,K)*OPK3*CB(I,jpole+1,K)
     &                +PDYP2*TTV(I,jpole,  K)*OPK0*CB(I,jpole,  K))*SIGL(K)
C
              FFYY=PY1+PY2
              DV(I,jpole,K)=-VL1-VL2-VL3+FSU-FFYY*WTGV(jpole)+SV(I,jpole,K)

C
            ENDDO
            DV(1,jpole,K)=DV(I1,jpole,K)
            DV(NX,jpole,K)=DV(2,jpole,K)
          ENDDO
#if (defined SPMD)
       ENDIF
#endif


#if (!defined SPMD)
       jpole =ny
#else
       if (myrank.eq.0) then
#endif
!!        (ELSE IF(J.EQ.NY) THEN)  ! for the south pole
          DO K=1,NL
            OZ2=0.5D0/DSIG(K)
            OZ4=0.5D0*OZ2
C
C       CALCULATING ADVECTION OF TT AND VERTICAL MOTION OVER PRESSURE
C       AT POLES.
C
            TLP=0.0D0
            O2P=0.0D0
            DO I=2,I1
              TLP=TLP+VS(I,jpole-1,K)*TT(I,jpole-1,K)
              O2P=O2P+(P(1,jpole)-P(I,jpole-1))*V(I,jpole-1,K)
            ENDDO
            WKP=2.0D0/RI0*OY1*WTGV(jpole-1)
            TLP=-WKP*TLP
            O2P=2.0D0*WKP*O2P*SIGL(K)/PLY(1,jpole,K)
            O1P=(0.5D0*PS(1,jpole)*(WS(1,jpole,K+1)+WS(1,jpole,K))
     &        +DPS(1,jpole)*SIGL(K))/PLY(1,jpole,K)
            TL3P=OZ2*(WS(1,jpole,K+1)*TTZ(1,jpole,K+2)
     &             -WS(1,jpole,K)*TTZ(1,jpole,K))
C
            TO1P=CB(1,jpole,K)*P(1,jpole)*O1P
            TO2P=CB(1,jpole,K)*P(1,jpole)*O2P
            TO3P=(DLT1*CAPA-DLT2*DCB(1,jpole,K))*TTV(1,jpole,K)*(O1P+O2P)
C
C     CALCULATING DTT/DT AND FORMING BOTH POLAR BOUNDARIES.
            DTTP=-TLP-TL3P+TO1P+TO2P+TO3P
            DO I=1,NX
              DTT(I,jpole,K)=DTTP+ST(I,jpole,K)
              DU(I,jpole,K)=SU(I,jpole,K)
              DV(I,jpole,K)=SV(I,jpole,K)
            ENDDO
          ENDDO
!!        ELSE
#if (defined SPMD)
       endif

       do j=jbeg,jend
#else
       do j=2,ny-1
#endif
C
C     CALCULATING DTT/DT FROM J=2 TO J=J0.
C
C     FORMING THE FACTORS INDEPENDENT ON I, K
C
          OUX2=0.5D0*OUX(J)
          R11=SINV(J)*OUY(J)
          R12=0.5D0*R11
          R21=SINV(J-1)*OUY(J)
          R22=0.5D0*R21
          OUX4=0.5D0*OUX2
          R14=0.5D0*R12
          R24=0.5D0*R22
          R10=R14/(OY1*WTGV(J))
          R20=R24/(OY1*WTGV(J-1))
          OVX4=0.25D0*OVX(J)
          OVY4=0.25D0*OVY(J)
          RV1=SINV(J+1)*OVY4
          RV2=SINV(J-1)*OVY4
	  OYY4=OY4*WTGV(J)
C
          DO K=1,NL
            OZ2=0.5D0/DSIG(K)
            OZ4=0.5D0*OZ2
C
            DO I=2,I1
C
              IP1=MP1(I,J)
              IM1=MM1(I,J)
              IP2=MP2(I,J)
              IM2=MM2(I,J)
C
C     FORMING THE FACTORS INDEPENDENT ON K.
C
              DXP0=P(I,J)-P(IM1,J)
              DXP1=P(IP1,J)-P(I,J)
              DYP0=P(I,J+1)-P(I,J)
              DYP2=P(I,J)-P(I,J-1)
C
C     TAKING THE ARRAY ELEMENTS APPEARING IN FOLLOWING FORMULAS
C     REPEATEDLY AND PLACING THEM INTO WORKING UNITS.
C
              OPK0=1.0D0/PLY(I,J,K)
C
C     CALCULATING DTT/DT FOR J=2--J0, I=2--I1 ANDK=1--NL.
C
              TL1=OUX2*(US(IP2,J,K)*TT(IP1,J,K)-US(IM2,J,K)*TT(IM1,J,K))
              TL2=R12*VS(I,J,K)*TT(I,J+1,K)-R22*VS(I,J-1,K)*TT(I,J-1,K)
              TL3=OZ2*(WS(I,J,K+1)*TTZ(I,J,K+2)-WS(I,J,K)*TTZ(I,J,K))
              O1=OPK0*(0.5D0*PS(I,J)*(WS(I,J,K+1)+WS(I,J,K))
     &          +DPS(I,J)*SIGL(K))
              O2=OPK0*SIGL(K)*(OUX(J)*(DXP1*U(IP2,J,K)+DXP0*U(IM2,J,K))
     &                  +R11*DYP0*V(I,J,K)+R21*DYP2*V(I,J-1,K))
              TO1=CB(I,J,K)*P(I,J)*O1
              TO2=CB(I,J,K)*P(I,J)*O2

              TO3=(DLT1*CAPA-DLT2*DCB(I,J,K))*TTV(I,J,K)*(O1+O2)
              DTT(I,J,K)=-TL1-TL2-TL3+TO1+TO2+TO3+ST(I,J,K)
            ENDDO
            DTT(1,J,K)=DTT(I1,J,K)
            DTT(NX,J,K)=DTT(2,J,K)
C
            DO I=2,I1
C
              IP1=MP1(I,J)
              IM1=MM1(I,J)
              IP2=MP2(I,J)
              IM2=MM2(I,J)
              IP3=MP3(I,J)
              IM3=MM3(I,J)
C
C     FORMING THE FACTORS INDEPENDENT ON K.
C
              DXP0=P(IP3,J)-P(IM3,J)
              PDXP1=P(IP3,J)*DXP0
              PDXP2=P(IM3,J)*DXP0
              DXP1=P(IP2,J)-P(IM2,J)
              AXP0=0.25D0*(P(IP3,J)+P(IM3,J))
              DYP0=P(I,J+1)-P(I,J)
              PDYP1=P(I,J+1)*DYP0
              PDYP2=P(I,J)*DYP0
              DYP2=P(I,J)-P(I,J-1)
              AYP0=0.25D0*(P(I,J+1)+P(I,J))
C
C     TAKING THE ARRAY ELEMENTS APPEARING IN FOLLOWING FORMULAS
C     REPEATEDLY AND PLACING THEM INTO WORKING UNITS.
C
              OPK0=1.0D0/PLY(I,J,K)
              OPK1=1.0D0/PLY(IP3,J,K)
              OPK2=1.0D0/PLY(IM3,J,K)
              OPK3=1.0D0/PLY(I,J+1,K)
              H1=H(I,J,K+1)+H(I,J,K)
!             H0=H(IP3,J,K+1)+H(IP3,J,K)
              H0=HH(IP3,J,K+1)+HH(IP3,J,K)
C
C     CALCULATING DU/DT, DV/DT J=2--J0, I=2--I1 AND K=1--NL.
C
              UL1=OUX4*((US(I,J,K)+US(IP1,J,K))*U(IP1,J,K)
     &                 -(US(I,J,K)+US(IM1,J,K))*U(IM1,J,K))
              UL2=R14*(VS(I,J,K)+VS(I-1,J,K))*U(I,J+1,K)
     &           -R24*(VS(I,J-1,K)+VS(I-1,J-1,K))*U(I,J-1,K)
              UL3=OZ4*((WS(I,J,K+1)+WS(I-1,J,K+1))*UZ(I,J,K+2)
     &                -(WS(I,J,K)+WS(I-1,J,K))*UZ(I,J,K))
              FS0=FF(J)+CUR(J)*US(I,J,K)
              FSV=FS0*(R10*(V(I,J,K)+V(I-1,J,K))
     &                +R20*(V(I,J-1,K)+V(I-1,J-1,K)))
!             PX1=OUX(J)*(AXP0*(H0-H(IM3,J,K+1)-H(IM3,J,K)))
              PX1=OUX(J)*(AXP0*(H0-HH(IM3,J,K+1)-HH(IM3,J,K)))
              PX2=OUX(J)*(PDXP1*TTV(IP3,J,K)*OPK1*CB(IP3,J,K)
     &                   +PDXP2*TTV(IM3,J,K)*OPK2*CB(IM3,J,K))*SIGL(K)
              FFXX=PX1+PX2
C
              VL1=OVX4*((US(IP2,J+1,K)+US(IP2,J,K))*V(IP1,J,K)
     &                 -(US(IM2,J+1,K)+US(IM2,J,K))*V(IM1,J,K))
              VL2=(RV1*VS(I,J+1,K)+OYY4*VS(I,J,K))*V(I,J+1,K)
     &           -(RV2*VS(I,J-1,K)+OYY4*VS(I,J,K))*V(I,J-1,K)
              VL3=OZ4*((WS(I,J+1,K+1)+WS(I,J,K+1))*VZ(I,J,K+2)
     &                -(WS(I,J+1,K)+WS(I,J,K))*VZ(I,J,K))
              FSU=0.25D0*(FS0*U(I,J,K)
     &                   +(FF(J)+CUR(J)*US(I+1,J,K))*U(I+1,J,K)
     &                   +(FF(J+1)+CUR(J+1)*US(I,J+1,K))*U(I,J+1,K)
     &                   +(FF(J+1)+CUR(J+1)*US(I+1,J+1,K))*U(I+1,J+1,K))
              PY1=OY1*(AYP0*(H(I,J+1,K+1)+H(I,J+1,K)-H1))
              PY2=OY1*(PDYP1*TTV(I,J+1,K)*OPK3*CB(I,J+1,K)
     &                +PDYP2*TTV(I,J,K)*OPK0*CB(I,J,K))*SIGL(K)
              FFYY=PY1+PY2
C
              DU(I,J,K)=-UL1-UL2-UL3-FSV-FFXX+SU(I,J,K)
              DV(I,J,K)=-VL1-VL2-VL3+FSU-FFYY*WTGV(J)+SV(I,J,K)
            ENDDO
            DU(1,J,K)=DU(I1,J,K)
            DU(NX,J,K)=DU(2,J,K)
            DV(1,J,K)=DV(I1,J,K)
            DV(NX,J,K)=DV(2,J,K)
          ENDDO
!!        (ENDIF)
	END DO
!- check ---------------------------------------------------------
!
!#if (defined SPMD)
!      write(filename,14) 'difuvt-p-',myrank,'.out'
!14    format(a9,i1,a4)
!      open (10,file=trim(filename))
!#else
!      open (10,file='difuvt-s.out')
!#endif
!
!      write(10,*) '------------------- du 26 -----------------'
!      write(10,11) (j,(du(i,j,26),i=1,2),j=1,ny)
!
!      write(10,*) '------------------- dv 26 -----------------'
!      write(10,11) (j,(dv(i,j,26),i=1,2),j=1,ny)
!
!      write(10,*) '------------------- dtt 26 -----------------'
!      write(10,11) (j,(dtt(i,j,26),i=1,2),j=1,ny)
!
!      write(10,*) '------------------- du 12 -----------------'
!      write(10,11) (j,(du(i,j,12),i=1,2),j=1,ny)
!
!      write(10,*) '------------------- dv 12 -----------------'
!      write(10,11) (j,(dv(i,j,12),i=1,2),j=1,ny)
!
!      write(10,*) '------------------- dtt 12 -----------------'
!      write(10,11) (j,(dtt(i,j,12),i=1,2),j=1,ny)
!
!11    format(1x,i5,2e30.20)
!
!
!      close (10)
!
!#if (defined SPMD)
!!      call mpi_finalize(j)
!#endif
!!      stop 'difuvt'
!
!-----------------------------------------------------------------


!
	RETURN
	END
