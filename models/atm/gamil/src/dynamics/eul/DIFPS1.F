#include <misc.h>
#include <params.h>

!!(2003.11.10)
!***********************************************************************
!
      SUBROUTINE DIFPS1(U,V,P,PS,WS,DPS,DPSU
     _                ,DSIG,DY,OUX,OUY,SINV,MP1,MP2,MM1,MM2,WTGV)
     _

#if ( defined SPMD )
       use mpishorthand, only: mpicom
#endif

	IMPLICIT NONE
!
!	This subroutine is to calculate the tendency of the surface prressure DPS,
!     the vertical velocity WS, the zonal wind US, meridional wind VS and the
!	departure of the geopotential height from the standard atmopshere

#include <PARADYN>

#if (defined SPMD)
#include <commpi.h>
#include <mpif.h>

      integer isend1,irecv1,isend2,irecv2
      integer istatus(mpi_status_size)
      integer ii
      character*50 filename

      real*8 worksb( (nx-2)*(nl+1) )
      real*8 workrt( (nx-2)*(nl+1) )
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
     _       U  (NX,NY,NL)  !  U=u*sqrt(Ps), input variable
     _      ,V  (NX,NY,NL)  !	 V=v*sqrt(Ps), input variable
     _      ,P  (NX,NY   )  !  P=sqrt(Ps)  , input variable
     _      ,PS (NX,NY   )  !	 Surface pressure, input variable
     _      ,DSIG(NL     )  !  The vertical stepsizes, input constant
     _      ,DY             !  The horizontal stepsize in meridional direction
!                              input constant
     _      ,SINV(NY)       !  sin(theta) at half grid j+1/2, input constant
     _      ,OUX(NY)        !  OUX=1/(RAD*SINU*DX*MDJ), input constant
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
     _      ,OUY(NY)        !  OUY=1/(RAD*SINU*DY*WTGU), input constant,
!                              where WTGU is the weighting at the integer grid j
     _      ,WTGV(NY)       !	 the weighting at the half grid j+1/2,
!                              input constant
     _      ,WS  (NX,NY,NZ) !  WS = w, vertical velocity, output variable
     _      ,DPS (NX,NY   ) !  the tendency of the surface pressure,
     _      ,DPSU(NX,NY   ) !  the tendency of the surface pressure,
!                              output variable
      REAL*8 WK0,WK1,WK2,WK3 !  working variables
     _      ,PXP,PXM,PYP	  !  working variables
     _      ,PYM,DPSP,WKQ   !  working variables
     _      ,D1(NZ),W1(NZ)  !  working variables
!
	INTEGER
     _       MM1(NX,NY)     !
     _      ,MP1(NX,NY)     !
     _      ,MM2(NX,NY)     !
     _      ,MP2(NX,NY)     !
!
      INTEGER I,I1,J,K,IP1,IM1,IP2,IM2 !  working variables
      integer begj,endj

      I1=NX-1

#if (defined SPMD)
      begj=2
      endj=ny-1
#else
      begj=1
      endj=ny
#endif


#if (defined SPMD)

      ii=1
      do i=2,nx-1
         worksb(ii)=p(i,ny-1)
         ii=ii+1
      enddo
      do k=1,nl
       do i=2,nx-1
         worksb(ii)=v(i,ny-1,k)
         ii=ii+1
       enddo
      enddo

      call mpi_isend( worksb,(nx-2)*(nl+1),mpi_double_precision,
     _                              ibot,1,mpicom,isend1,ierr)
      call mpi_isend( p(2,2), nx-2        ,mpi_double_precision,
     _                              itop,2,mpicom,isend2,ierr)
      call mpi_irecv( workrt,(nx-2)*(nl+1),mpi_double_precision,
     _                              itop,1,mpicom,irecv1,ierr)
      call mpi_irecv( p(2,ny),nx-2        ,mpi_double_precision,
     _                              ibot,2,mpicom,irecv2,ierr)

      call mpi_wait(isend1,istatus,ierr)
      call mpi_wait(isend2,istatus,ierr)
      call mpi_wait(irecv1,istatus,ierr)
      call mpi_wait(irecv2,istatus,ierr)

      ii=1
      do i=2,nx-1
         p(i,1)  =workrt(ii)
         ii=ii+1
      enddo
      do k=1,nl
       do i=2,nx-1
         v(i,1,k)=workrt(ii)
         ii=ii+1
       enddo
      enddo

#endif

C     DEALING WITH INTERNAL GRID POINTS.
C     CALCULATING FACTORS RELATIVE TO J IN DIVERGENCE FORMULA.
C
!     (DO J=1,NY)
      DO J=begj,endj
        DO I=1,NX
          WS(I,J,1)=0.0D0
          WS(I,J,NZ)=0.0D0
        ENDDO
      enddo

#if (!defined SPMD)
      do j=2,ny-1
#else
      do j=jbeg,jend
#endif
!!        IF(J.GE.2.AND.J.LE.NY-1) THEN
          WK1=0.5D0*OUX(J)
          WK2=0.5D0*OUY(J)*SINV(J)
          WK3=0.5D0*OUY(J)*SINV(J-1)
C     CALCULATING FACTORS RELATIVE TO I AND J IN DIVERGENCE FORMULA.
C
          DO I=2,I1
            IP1=MP1(I,J)
            IM1=MM1(I,J)
            IP2=MP2(I,J)
            IM2=MM2(I,J)
            PXP=WK1*(P(IP1,J)+P(I,J))
            PXM=WK1*(P(I,J)+P(IM1,J))
            PYP=WK2*(P(I,J+1)+P(I,J))
            PYM=WK3*(P(I,J)+P(I,J-1))
C     CALCULATING DIVERGENCES AS WELL AS SUM OF THEM.
            DPSP=0.0D0
            DO K=1,NL
              D1(K)=PYP*V(I,J,K)-PYM*V(I,J-1,K)
              DPSP=DPSP-DSIG(K)*D1(K)
              D1(K)=D1(K)+PXP*U(IP2,J,K)-PXM*U(IM2,J,K)
            ENDDO
C
C     CALCULATING DPS/DT, DPS/DT AND D(SIGMA)/DT.
C
            DPS(I,J)=DPSP
	    DPSP=DPSP+DPSU(I,J)
            WKQ=1.D0/PS(I,J)
            DO K=2,NL
              WS(I,J,K)=WS(I,J,K-1)-DSIG(K-1)*WKQ*(DPSP+D1(K-1))
            ENDDO
          ENDDO
C     CYCLIC EXTENSION TO I=1 AND I=NX.
          DPS(1,J)=DPS(I1,J)
          DPS(NX,J)=DPS(2,J)
          DO K=2,NL
            WS(1,J,K)=WS(I1,J,K)
            WS(NX,J,K)=WS(2,J,K)
          ENDDO
      enddo


!!        ELSE IF(J.EQ.1) THEN
C
C     FINDING DP/DT AND D(SIGMA)/DT AT POLES.
C     IN BELOW, SUBSCRIPTS 1 AND 2 REPRESENT J=1 AND NY RESPECTIVELY.

#if (defined SPMD)
      if ( myrank.eq.(nprocs-1) ) then
#else
         jpole = 1
#endif
          W1(1)=0.0D0
          WK0=2.0D0/(DFLOAT(NX-2)*RAD*DY)*WTGV(jpole)
C     CALCULATING DIVERGENCE AT POLES.
          DPSP=0.0D0
          DO K=1,NL
            D1(K)=0.0D0
            DO I=2,I1
              D1(K)=D1(K)+(P(1,jpole)+P(I,jpole+1))*V(I,jpole,K)
            ENDDO
            D1(K)=WK0*D1(K)
            DPSP=DPSP-DSIG(K)*D1(K)
          ENDDO
C
C     CALCULATING DPS/DT AND D(SIGMA)/DT AT POLES.
C
          DO K=2,NL
            W1(K)=W1(K-1)-DSIG(K-1)/PS(1,jpole)*(DPSP+D1(K-1))
          ENDDO
C
          DO I=1,NX
            DPS(I,jpole)=DPSP
          ENDDO
          DO K=2,NL
            DO I=1,NX
              WS(I,jpole,K)=W1(K)
            ENDDO
          ENDDO
#if (defined SPMD)
      endif
#endif

!!    ELSE IF(J.EQ.NY) THEN
C
C     FINDING DP/DT AND D(SIGMA)/DT AT POLES.
C     IN BELOW, SUBSCRIPTS 1 AND 2 REPRESENT J=1 AND NY RESPECTIVELY.

#if (defined SPMD)
      if ( myrank.eq.0 ) then
#else
         jpole = ny
#endif
          W1(1)=0.0D0
          WK0=2.0D0/(DFLOAT(NX-2)*RAD*DY)*WTGV(jpole-1)
C     CALCULATING DIVERGENCE AT POLES.
          DPSP=0.0D0
          DO K=1,NL
            D1(K)=0.0D0
            DO I=2,I1
              D1(K)=D1(K)+(P(1,jpole)+P(I,jpole-1))*V(I,jpole-1,K)
            ENDDO
            D1(K)=-WK0*D1(K)
            DPSP=DPSP-DSIG(K)*D1(K)
          ENDDO
C
C     CALCULATING DPS/DT AND D(SIGMA)/DT AT POLES.
C
          DO K=2,NL
            W1(K)=W1(K-1)-DSIG(K-1)/PS(1,jpole)*(DPSP+D1(K-1))
          ENDDO
C
          DO I=1,NX
            DPS(I,jpole)=DPSP
          ENDDO
          DO K=2,NL
            DO I=1,NX
              WS(I,jpole,K)=W1(K)
            ENDDO
          ENDDO
!!        ENDIF
#if (defined SPMD)
      endif
#endif

!- check ---------------------------------------------------------
!
!#if (defined SPMD)
!      write(filename,14) 'difps-p-',myrank,'.out'
!14    format(a8,i1,a4)
!      open (10,file=trim(filename))
!#else
!      open (10,file='difps-s.out')
!#endif
!
!      write(10,*) '------------------- ws -----------------'
!      write(10,11) (j,(ws(i,j,10),i=1,2),j=1,ny)
!
!      write(10,*) '------------------- dps -----------------'
!      write(10,11) (j,(dps(i,j),i=1,2),j=1,ny)
!11    format(1x,i5,2e30.20)
!
!      close (10)
!
!#if (defined SPMD)
!      call mpi_finalize(j)
!#endif
!      stop 'difps'
!
!!-----------------------------------------------------------------



      RETURN
      END
