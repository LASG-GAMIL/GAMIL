#include <misc.h>
#include <params.h>

!!(wh 2003.11.10)
!!-------------------


	FUNCTION INNER(U1,V1,T1,P1,L1,U2,V2,T2,P2,L2
     _                ,DSIG,TBS,PSB,PMTOP,SINU,SINV,WTGU,WTGV)

#if ( defined SPMD )
        use mpishorthand, only: mpicom
#endif

	IMPLICIT NONE
!
!	This subroutine is to define the inner product of the dynamical system
!     of the atmosphere. Here, INNER is the inner product of the vector
!     (U1,V1,T1,P1) and the vector (U2,V2,T2,P2)
!
#include <PARADYN>

#if (defined SPMD)
#include <commpi.h>
#include <mpif.h>

      integer nlat
      parameter(nlat=nprocessor*(ny-2))

      integer displs
      integer jj

      real*8 sinu_mpi(nlat)
      real*8 sinv_mpi(nlat)
      real*8 wtgu_mpi(nlat)
      real*8 wtgv_mpi(nlat)
      real*8 euj_mpi(nlat)
      real*8 evj_mpi(nlat)
      real*8 etj_mpi(nlat)
      real*8 epj_mpi(nlat)

      real*8 isend1,irecv1
      real*8 isend2,irecv2
      real*8 isend3,irecv3
      real*8 isend4,irecv4
      real*8 isend5,irecv5
      real*8 isend6,irecv6
      real*8 isend7,irecv7
      real*8 isend8,irecv8
      integer istatus(mpi_status_size)
#endif
!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	REAL*8
     _       U1(NX,NY,NL)  !	input variable
     _      ,V1(NX,NY,NL)  !	input variable
     _      ,T1(NX,NY,NL)  !	input variable
     _      ,P1(NX,NY)     !	input variable
     _      ,U2(NX,NY,NL)  !	input variable
     _      ,V2(NX,NY,NL)  !	input variable
     _      ,T2(NX,NY,NL)  !	input variable
     _      ,P2(NX,NY)     !	input variable
     _      ,TBS(NX,NY)    !	input variable
     _      ,PSB(NX,NY)    !	input variable
     _      ,SINU(NY)      !	input variable
     _      ,SINV(NY)      !	input variable
     _      ,WTGU(NY)		 !	input variable
     _      ,WTGV(NY)      !	input variable
     _      ,DSIG(NL)		 !	input variable
     _      ,PMTOP		 !	input variable
     _      ,INNER         !  output variable
!
	INTEGER
     _        L1			 !	input variable
     _       ,L2			 !	input variable
     _       ,I1,I,J,K     !  working variables
!
	REAL*8	TMPU(NX,NY)	 !  working variables
     _       ,TMPV(NX,NY)	 !  working variables
     _       ,TMPT(NX,NY)	 !  working variables
     _       ,EUJ(NY)      !  working variables
     _       ,EVJ(NY)      !  working variables
     _       ,ETJ(NY)      !  working variables
     _       ,EPJ(NY)      !  working variables
     _       ,DS           !  working variables
     _       ,DSU          !  working variables
     _       ,DSV          !  working variables
     _       ,DPS1         !  working variables
     _       ,DPS2         !  working variables

        integer begj,endj

#if (defined SPMD)
        begj = 2
        endj = ny-1
#else
        begj = 1
        endj = ny
#endif
!
	I1=NX-1
!
!      (DO J=1,NY)
       DO J=begj,endj
        DO I = 2,I1
          TMPU(I,J) = 0.0D0
          TMPV(I,J) = 0.0D0
          TMPT(I,J) = 0.0D0
        ENDDO
        DO K = 1,NL
          DO I = 2,I1
            TMPU(I,J) = TMPU(I,J)+U1(I,J,K)*U2(I,J,K)*DSIG(K)
            TMPV(I,J) = TMPV(I,J)+V1(I,J,K)*V2(I,J,K)*DSIG(K)
            TMPT(I,J) = TMPT(I,J)+T1(I,J,K)*T2(I,J,K)*DSIG(K)
          ENDDO
        ENDDO
        EUJ(J) = 0.0D0
        EVJ(J) = 0.0D0
        ETJ(J) = 0.0D0
        EPJ(J) = 0.0D0
        DO I = 2 ,I1
!
	    IF (L1.EQ.1) THEN
	       DPS1 = P1(I,J)+PMTOP-PSB(I,J)
	    ELSE
	       DPS1 = P1(I,J)
	    END IF
!
	    IF (L2.EQ.1) THEN
	       DPS2 = P2(I,J)+PMTOP-PSB(I,J)
	    ELSE
	       DPS2 = P2(I,J)
	    END IF
!
          EPJ(J) = EPJ(J) + RD*TBS(I,J)/PSB(I,J)*DPS1*DPS2
          EUJ(J) = EUJ(J) + TMPU(I,J)
          EVJ(J) = EVJ(J) + TMPV(I,J)
          ETJ(J) = ETJ(J) + TMPT(I,J)
        ENDDO
      ENDDO
!
      INNER=0.0D0
!
#if (!defined SPMD)
      DO J=1,NY
        IF(J.EQ.1) THEN
          DS=0.25D0*SINV(1)/WTGV(1)
          DSU=DS
          DSV=4.0D0*DSU
        ELSE IF(J.EQ.NY) THEN
          DS=0.25D0*SINV(NY-1)/WTGV(NY-1)
          DSU=DS
          DSV=0.0D0
        ELSE
          DS=SINU(J)/WTGU(J)
          DSU=DS
          DSV=SINV(J)/WTGV(J)
        ENDIF
        INNER=INNER+DSU*EUJ(J)+DSV*EVJ(J)+DSU*ETJ(J)+DS*EPJ(J)
      ENDDO

!      write(*,*) 'inner: inner=',inner
!      stop
#else

      if (myrank.ne.0) then
         call mpi_isend(sinu(2),ny-2,mpi_double_precision,0,1,mpicom,isend1,ierr)
         call mpi_isend(sinv(2),ny-2,mpi_double_precision,0,2,mpicom,isend2,ierr)
         call mpi_isend(wtgu(2),ny-2,mpi_double_precision,0,3,mpicom,isend3,ierr)
         call mpi_isend(wtgv(2),ny-2,mpi_double_precision,0,4,mpicom,isend4,ierr)
         call mpi_isend( euj(2),ny-2,mpi_double_precision,0,5,mpicom,isend5,ierr)
         call mpi_isend( evj(2),ny-2,mpi_double_precision,0,6,mpicom,isend6,ierr)
         call mpi_isend( etj(2),ny-2,mpi_double_precision,0,7,mpicom,isend7,ierr)
         call mpi_isend( epj(2),ny-2,mpi_double_precision,0,8,mpicom,isend8,ierr)
         call mpi_wait(isend1,istatus,ierr)
         call mpi_wait(isend2,istatus,ierr)
         call mpi_wait(isend3,istatus,ierr)
         call mpi_wait(isend4,istatus,ierr)
         call mpi_wait(isend5,istatus,ierr)
         call mpi_wait(isend6,istatus,ierr)
         call mpi_wait(isend7,istatus,ierr)
         call mpi_wait(isend8,istatus,ierr)
      else
       do jj=1,nprocessor-1
         displs=(nprocessor-1-jj)*(ny-2)
         call mpi_irecv(sinu_mpi(displs+1),ny-2,mpi_double_precision
     _                                    ,jj,1,mpicom,     irecv1,ierr)
         call mpi_irecv(sinv_mpi(displs+1),ny-2,mpi_double_precision
     _                                    ,jj,2,mpicom,     irecv2,ierr)
         call mpi_irecv(wtgu_mpi(displs+1),ny-2,mpi_double_precision
     _                                    ,jj,3,mpicom,     irecv3,ierr)
         call mpi_irecv(wtgv_mpi(displs+1),ny-2,mpi_double_precision
     _                                    ,jj,4,mpicom,     irecv4,ierr)
         call mpi_irecv( euj_mpi(displs+1),ny-2,mpi_double_precision
     _                                    ,jj,5,mpicom,     irecv5,ierr)
         call mpi_irecv( evj_mpi(displs+1),ny-2,mpi_double_precision
     _                                    ,jj,6,mpicom,     irecv6,ierr)
         call mpi_irecv( etj_mpi(displs+1),ny-2,mpi_double_precision
     _                                    ,jj,7,mpicom,     irecv7,ierr)
         call mpi_irecv( epj_mpi(displs+1),ny-2,mpi_double_precision
     _                                    ,jj,8,mpicom,     irecv8,ierr)
         call mpi_wait(irecv1,istatus,ierr)
         call mpi_wait(irecv2,istatus,ierr)
         call mpi_wait(irecv3,istatus,ierr)
         call mpi_wait(irecv4,istatus,ierr)
         call mpi_wait(irecv5,istatus,ierr)
         call mpi_wait(irecv6,istatus,ierr)
         call mpi_wait(irecv7,istatus,ierr)
         call mpi_wait(irecv8,istatus,ierr)
       enddo

       displs=(nprocessor-1)*(ny-2)
       do j=2,ny-1
          sinu_mpi(displs-1+j)=sinu(j)
       enddo
       do j=2,ny-1
          sinv_mpi(displs-1+j)=sinv(j)
       enddo
       do j=2,ny-1
          wtgu_mpi(displs-1+j)=wtgu(j)
       enddo
       do j=2,ny-1
          wtgv_mpi(displs-1+j)=wtgv(j)
       enddo
       do j=2,ny-1
           euj_mpi(displs-1+j)= euj(j)
       enddo
       do j=2,ny-1
           evj_mpi(displs-1+j)= evj(j)
       enddo
       do j=2,ny-1
           etj_mpi(displs-1+j)= etj(j)
       enddo
       do j=2,ny-1
           epj_mpi(displs-1+j)= epj(j)
       enddo
      endif

       if (myrank.eq.0) then
        DO J=1,nlat
          IF(J.EQ.1) THEN
            DS=0.25D0*SINV_mpi(1)/WTGV_mpi(1)
            DSU=DS
            DSV=4.0D0*DSU
          ELSE IF(J.EQ.nlat) THEN
            DS=0.25D0*SINV_mpi(nlat-1)/WTGV_mpi(nlat-1)
            DSU=DS
            DSV=0.0D0
          ELSE
            DS=SINU_mpi(J)/WTGU_mpi(J)
            DSU=DS
            DSV=SINV_mpi(J)/WTGV_mpi(J)
          ENDIF
          INNER=INNER+DSU*EUJ_mpi(J)+DSV*EVJ_mpi(J)+DSU*ETJ_mpi(J)+DS*EPJ_mpi(J)
        ENDDO
       endif

      call mpi_bcast(inner,1,mpi_double_precision,0,mpicom,ierr)

!      write(*,*)'inner: inner=',inner
!      call mpi_finalize(ierr)
!      stop 'inner'

#endif



!
	RETURN
	END
