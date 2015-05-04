#include <misc.h>
#include <params.h>

!!(2003.11.29)
!***********************************************************************
!
        SUBROUTINE DIFPS0(U,P,DPS,DSIG,OUX,MP1,MP2,MM1,MM2)
	IMPLICIT NONE
!
!	This subroutine is to calculate the tendency of the surface prressure DPS,
!     the vertical velocity WS, the zonal wind US, meridional wind VS and the
!	departure of the geopotential height from the standard atmopshere
!
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
     _      ,P  (NX,NY   )  !  P=sqrt(Ps)  , input variable
     _      ,DSIG(NL     )  !  The vertical stepsizes, input constant
     _      ,OUX(NY)        !  OUX=1/(RAD*SINU*DX*MDJ), input constant
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
     _      ,DPS(NX,NY   )  !  the tendency of the surface pressure,
!                              output variable
      REAL*8   WK1,WK2,WK3    !  working variables
     _      ,PXP,PXM,DPSP   !  working variables
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


C     DEALING WITH INTERNAL GRID POINTS.
C     CALCULATING FACTORS RELATIVE TO J IN DIVERGENCE FORMULA.
C

#if (!defined SPMD)
      do j=2,ny-1
#else
      do j=jbeg,jend
#endif
!!        IF(J.GE.2.AND.J.LE.NY-1) THEN
          WK1=0.5D0*OUX(J)
C     CALCULATING FACTORS RELATIVE TO I AND J IN DIVERGENCE FORMULA.
C
          DO I=2,I1
            IP1=MP1(I,J)
            IM1=MM1(I,J)
            IP2=MP2(I,J)
            IM2=MM2(I,J)
            PXP=WK1*(P(IP1,J)+P(I,J))
            PXM=WK1*(P(I,J)+P(IM1,J))
C     CALCULATING DIVERGENCES AS WELL AS SUM OF THEM.
            DPSP=0.0D0
            DO K=1,NL
              DPSP=DPSP-DSIG(K)*(PXP*U(IP2,J,K)-PXM*U(IM2,J,K))
            ENDDO
C
C     CALCULATING DPS/DT, DPS/DT AND D(SIGMA)/DT.
C
            DPS(I,J)=DPSP
          ENDDO
C     CYCLIC EXTENSION TO I=1 AND I=NX.
          DPS(1,J)=DPS(I1,J)
          DPS(NX,J)=DPS(2,J)
        END DO
!!        ELSE IF(J.EQ.1) THEN
C
C     FINDING DP/DT AND D(SIGMA)/DT AT POLES.
C     IN BELOW, SUBSCRIPTS 1 AND 2 REPRESENT J=1 AND NY RESPECTIVELY.

#if (defined SPMD)
      if ( myrank.eq.(nprocs-1) ) then
#else
         jpole = 1
#endif
          DO I=1,NX
            DPS(I,jpole)=0.0
          ENDDO
!
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
          DO I=1,NX
            DPS(I,jpole)=0.0
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
