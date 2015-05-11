#include <misc.h>
#include <params.h>

!!(wh 2003.11.10)
!!-------------------


      SUBROUTINE TEM(U,V,PS,P,TT,HS,CB,TB,DX,DY,DSIG
     _              ,SINU,SINV,WTGU,WTGV,TE,TM)

#if ( defined SPMD )
      use mpishorthand, only: mpicom
#endif

      IMPLICIT NONE

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
      real*8 eej_mpi(nlat)
      real*8 esj_mpi(nlat)
      real*8 tmj_mpi(nlat)

      real*8 isend1,irecv1
      real*8 isend2,irecv2
      real*8 isend3,irecv3
      real*8 isend4,irecv4
      real*8 isend5,irecv5
      real*8 isend6,irecv6
      real*8 isend7,irecv7
      real*8 isend8,irecv8
      real*8 isend9,irecv9
      integer istatus(mpi_status_size)
#endif

!
      REAL*8 U(NX,NY,NL),V(NX,NY,NL),PS(NX,NY),P(NX,NY),TT(NX,NY,NL)
!
      REAL*8 TB(NX,NY,NL),CB(NX,NY,NL),HS(NX,NY)
     _      ,WTGU(NY),WTGV(NY),SINU(NY),SINV(NY),DX,DY,DSIG(NL)
!
      INTEGER I,J,K,I1
      REAL*8  DZ2,TE,EE,EK,DS0,DSU,DSV,DS,ES,TM
      REAL*8  TMJ(NY),EUJ(NY),EVJ(NY),EEJ(NY),ESJ(NY)

      integer begj,endj

#if (defined SPMD)
      begj=2
      endj=ny-1
#else
      begj=1
      endj=ny
#endif

      I1=NX-1
!
      DS0=RAD*RAD*DX*DY
!
!     (DO J = 1 ,NY)
      DO J = begj ,endj
        TMJ(J) = 0.0D0
        ESJ(J) = 0.0D0
        EUJ(J) = 0.0D0
        EVJ(J) = 0.0D0
        EEJ(J) = 0.0D0
        DO I = 2 ,I1
          TMJ(J)=TMJ(J)+PS(I,J)
          ESJ(J)=ESJ(J)+PS(I,J)*HS(I,J)
        ENDDO
        DO K = 1 ,NL
          DZ2=DSIG(K)*0.5D0
          DO I = 2 ,I1
            EUJ(J)=EUJ(J)+U(I,J,K)*U(I,J,K)*DZ2
            EVJ(J)=EVJ(J)+V(I,J,K)*V(I,J,K)*DZ2
            EEJ(J)=EEJ(J)+(CPD*PS(I,J)*TB(I,J,K)
     &       +CB(I,J,K)/CAPA*P(I,J)*TT(I,J,K))*DSIG(K)
          ENDDO
        ENDDO
      ENDDO
C
      EK=0.0D0
      EE=0.0D0
      ES=0.0D0
      TM=0.0D0

#if (!defined SPMD)

      DO 60 J=1,NY
      IF(J.EQ.1) THEN
        DS=0.25D0*SINV(1)*DS0/WTGV(1)
        DSU=DS
        DSV=4.0D0*DSU
      ELSE IF(J.EQ.NY) THEN
        DS=0.25D0*SINV(NY-1)*DS0/WTGV(NY-1)
        DSU=DS
        DSV=0.0D0
      ELSE
        DS=SINU(J)*DS0/WTGU(J)
        DSU=DS
        DSV=SINV(J)*DS0/WTGV(J)
      ENDIF
C
      EK=EK+DSU*EUJ(J)+DSV*EVJ(J)
      EE=EE+DSU*EEJ(J)
      ES=ES+DS*ESJ(J)
      TM=TM+DS*TMJ(J)
*     TM=TM-DS*TMJ(J)
 60   CONTINUE

      TE=EK+EE+ES
#else
      if (myrank.ne.0) then
         call mpi_isend(sinu(2),ny-2,mpi_double_precision,0,1,mpicom,isend1,ierr)
         call mpi_isend(sinv(2),ny-2,mpi_double_precision,0,2,mpicom,isend2,ierr)
         call mpi_isend(wtgu(2),ny-2,mpi_double_precision,0,3,mpicom,isend3,ierr)
         call mpi_isend(wtgv(2),ny-2,mpi_double_precision,0,4,mpicom,isend4,ierr)
         call mpi_isend( euj(2),ny-2,mpi_double_precision,0,5,mpicom,isend5,ierr)
         call mpi_isend( evj(2),ny-2,mpi_double_precision,0,6,mpicom,isend6,ierr)
         call mpi_isend( eej(2),ny-2,mpi_double_precision,0,7,mpicom,isend7,ierr)
         call mpi_isend( esj(2),ny-2,mpi_double_precision,0,8,mpicom,isend8,ierr)
         call mpi_isend( tmj(2),ny-2,mpi_double_precision,0,9,mpicom,isend9,ierr)
         call mpi_wait(isend1,istatus,ierr)
         call mpi_wait(isend2,istatus,ierr)
         call mpi_wait(isend3,istatus,ierr)
         call mpi_wait(isend4,istatus,ierr)
         call mpi_wait(isend5,istatus,ierr)
         call mpi_wait(isend6,istatus,ierr)
         call mpi_wait(isend7,istatus,ierr)
         call mpi_wait(isend8,istatus,ierr)
         call mpi_wait(isend9,istatus,ierr)
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
         call mpi_irecv( eej_mpi(displs+1),ny-2,mpi_double_precision
     _                                    ,jj,7,mpicom,     irecv7,ierr)
         call mpi_irecv( esj_mpi(displs+1),ny-2,mpi_double_precision
     _                                    ,jj,8,mpicom,     irecv8,ierr)
         call mpi_irecv( tmj_mpi(displs+1),ny-2,mpi_double_precision
     _                                    ,jj,9,mpicom,     irecv9,ierr)
         call mpi_wait(irecv1,istatus,ierr)
         call mpi_wait(irecv2,istatus,ierr)
         call mpi_wait(irecv3,istatus,ierr)
         call mpi_wait(irecv4,istatus,ierr)
         call mpi_wait(irecv5,istatus,ierr)
         call mpi_wait(irecv6,istatus,ierr)
         call mpi_wait(irecv7,istatus,ierr)
         call mpi_wait(irecv8,istatus,ierr)
         call mpi_wait(irecv9,istatus,ierr)
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
           eej_mpi(displs-1+j)= eej(j)
       enddo
       do j=2,ny-1
           esj_mpi(displs-1+j)= esj(j)
       enddo
       do j=2,ny-1
           tmj_mpi(displs-1+j)= tmj(j)
       enddo

      endif

      if (myrank.eq.0) then
         DO  J=1,nlat
            IF(J.EQ.1) THEN
              DS=0.25D0*SINV_mpi(1)*DS0/WTGV_mpi(1)
              DSU=DS
              DSV=4.0D0*DSU
            ELSE IF(J.EQ.nlat) THEN
              DS=0.25D0*SINV_mpi(nlat-1)*DS0/WTGV_mpi(nlat-1)
              DSU=DS
              DSV=0.0D0
            ELSE
              DS=SINU_mpi(J)*DS0/WTGU_mpi(J)
              DSU=DS
              DSV=SINV_mpi(J)*DS0/WTGV_mpi(J)
            ENDIF
C
           EK=EK+DSU*EUJ_mpi(J)+DSV*EVJ_mpi(J)
           EE=EE+DSU*EEJ_mpi(J)
           ES=ES+DS*ESJ_mpi(J)
           TM=TM+DS*TMJ_mpi(J)
         enddo

         TE=EK+EE+ES
      endif

      call mpi_bcast(tm,1,mpi_double_precision,0,mpicom,ierr)
      call mpi_bcast(te,1,mpi_double_precision,0,mpicom,ierr)

!      call mpi_finalize(ierr)
#endif

!      write(*,*) tm,te
!      stop

      RETURN
      END
