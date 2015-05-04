#include <misc.h>
#include <params.h>

      SUBROUTINE IAPTRSF(US,VS,WS,PS,TS,U,V,P,TT,TB,CB,KP)
!
#if ( defined SPMD )
      use mpishorthand, only: mpicom
#endif

      IMPLICIT NONE

#include <PARADYN>

#if (defined SPMD)
#include <mpif.h>
#include <commpi.h>

      integer   isend, irecv
      integer   istatus(mpi_status_size)
      character*50 filename
#endif

!
!     The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	REAL*8
     _	     US(NX,NY,NL )  !  US = u, zonal wind,    input variable
     _      ,VS(NX,NY,NL )  !  VS = v, meridional wind, input variable
     _      ,WS(NX,NY,NZ )  !  vertical velocity, input variable
     _      ,PS (NX,NY   )  !  Surface pressure, input variable
     _      ,TS (NX,NY,NL)  !  TS = T, temperature, input variable
     _      ,TB (NX,NY,NL)  !  the temperature of the standard atmosphere a
!                           !  at the sigma layers, input constant
     _      ,CB (NX,NY,NL)  !  CB=Sqrt(R*(KK*TBB-dTB/dlnp)), input variable,
!                           !  where, KK=R/Cp, R is a constant
     _      ,U  (NX,NY,NL)  !  U=u*sqrt(Ps), output variable
     _      ,V  (NX,NY,NL)  !	 V=v*sqrt(Ps), output variable
     _      ,P  (NX,NY   )  !  P=sqrt(Ps)  , output variable
     _      ,TT (NX,NY,NL)  !  TT=R*T'*Sqrt(Ps)/CB, output variable
!                           !  where T'=T-TB, T is the temperatur,
     _      ,PX,PY          !  working variables
!
	INTEGER I,J,K,KP
        integer begj,endj

#if (defined SPMD)
      begj = 2
      endj = ny-1
#else
      begj = 1
      endj = ny
#endif

	IF (KP.EQ.1) THEN
!
        DO J=begj,endj
	  DO I=2,NX-1
	     P(I,J)=SQRT(PS(I,J))
	  END DO
	  PS (1 ,J) = PS (NX-1,J)
	  PS (NX,J) = PS (2   ,J)
!
	  P  (1 ,J) = P  (NX-1,J)
	  P  (NX,J) = P  (2   ,J)
        END DO

#if (defined SPMD)
      call mpi_isend( p(1,2), nx, mpi_double_precision, itop, 1, mpicom, isend, ierr)
      call mpi_irecv( p(1,ny),nx, mpi_double_precision, ibot, 1, mpicom, irecv, ierr)
#endif
!

        DO J=begj,endj
	  DO I=2,NX-1
             WS(I,J,1)=0.0D0
             WS(I,J,NZ)=0.0D0
	  END DO
	  WS(1 ,J,1) = WS(NX-1,J, 1)
	  WS(NX,J,1) = WS(2   ,J, 1)
	  WS(1,J,NZ) = WS(NX-1,J,NZ)
	  WS(NX,J,NZ)= WS(2   ,J,NZ)
	END DO

#if (defined SPMD)
      call mpi_wait(isend,istatus,ierr)
      call mpi_wait(irecv,istatus,ierr)
#endif

!
      DO J=begj,endj
        DO I=2,NX-1
!
!         CALCULATING THE AVERAGE VALUE OF P AT V-GRID.
!
#if (!defined SPMD)
          IF(J.LT.NY) THEN
            PY=0.5D0*(P(I,J)+P(I,J+1))
          ELSE
            PY=0.0D0
          ENDIF
#else
         if ((myrank.eq.0).and.(j.eq.endj)) then
            PY=0.0D0
         else
            PY=0.5D0*(P(I,J)+P(I,J+1))
         endif
#endif
!
!         CALCULATING THE AVERAGE VALUE OF P AT U-GRID.
!
          PX=0.5D0*(P(I,J)+P(I-1,J))
!
          DO K=1,NL
!
             U (I,J,K)=PX*US(I,J,K)
             V (I,J,K)=PY*VS(I,J,K)
!
!         CALCULATING TT=(TS-TB)*P*RD/CB.
!
	     TT (I,J,K)=(TS(I,J,K)-TB(I,J,K))*P(I,J)*RD/CB(I,J,K)
!
          ENDDO
        ENDDO
!
	  DO K=1,NL
!
	    US (1 ,J,K) = US (NX-1,J,K)
	    US (NX,J,K) = US (2   ,J,K)
!
	    VS (1 ,J,K) = VS (NX-1,J,K)
	    VS (NX,J,K) = VS (2   ,J,K)
!
	    TS (1 ,J,K) = TS (NX-1,J,K)
	    TS (NX,J,K) = TS (2   ,J,K)
!
	    U  (1 ,J,K) = U  (NX-1,J,K)
	    U  (NX,J,K) = U  (2   ,J,K)
!
	    V  (1 ,J,K) = V  (NX-1,J,K)
	    V  (NX,J,K) = V  (2   ,J,K)
!
	    TT (1 ,J,K) = TT (NX-1,J,K)
	    TT (NX,J,K) = TT (2   ,J,K)
!
	  END DO
      ENDDO
!
	ELSE
!
        DO J=begj,endj
	  DO I=2,NX-1
	     P(I,J)=SQRT(PS(I,J))
	  END DO
!
	  PS (1 ,J) = PS (NX-1,J)
	  PS (NX,J) = PS (2   ,J)
!
	  P  (1 ,J) = P  (NX-1,J)
	  P  (NX,J) = P  (2   ,J)
!
	END DO

#if (defined SPMD)
      call mpi_isend( p(1,2), nx, mpi_double_precision, itop, 1, mpicom, isend, ierr)
      call mpi_irecv( p(1,ny),nx, mpi_double_precision, ibot, 1, mpicom, irecv, ierr)
      call mpi_wait(isend,istatus,ierr)
      call mpi_wait(irecv,istatus,ierr)
#endif

!
      DO J=begj,endj
        DO I=2,NX-1
!
!         CALCULATING THE AVERAGE VALUE OF P AT V-GRID.
!
#if (!defined SPMD)
          IF(J.LT.NY) THEN
            PY=0.5D0*(P(I,J)+P(I,J+1))
          ELSE
            PY=0.0D0
          ENDIF
#else
         if ((myrank.eq.0).and.(j.eq.endj)) then
            PY=0.0D0
         else
            PY=0.5D0*(P(I,J)+P(I,J+1))
         endif
#endif
!
!         CALCULATING THE AVERAGE VALUE OF P AT U-GRID.
!
          PX=0.5D0*(P(I,J)+P(I-1,J))
!
          DO K=1,NL
!
          U (I,J,K)=PX*US(I,J,K)
          V (I,J,K)=PY*VS(I,J,K)
!
!         CALCULATING TT=(TS-TB)*P*RD/CB.
!
	    TT (I,J,K)=(TS(I,J,K)-TB(I,J,K))*P(I,J)*RD/CB(I,J,K)
!
          ENDDO
        ENDDO
!
	  DO K=1,NL
!
	    US (1 ,J,K) = US (NX-1,J,K)
	    US (NX,J,K) = US (2   ,J,K)
!
	    VS (1 ,J,K) = VS (NX-1,J,K)
	    VS (NX,J,K) = VS (2   ,J,K)
!
	    TS (1 ,J,K) = TS (NX-1,J,K)
	    TS (NX,J,K) = TS (2   ,J,K)
!
	    U  (1 ,J,K) = U  (NX-1,J,K)
	    U  (NX,J,K) = U  (2   ,J,K)
!
	    V  (1 ,J,K) = V  (NX-1,J,K)
	    V  (NX,J,K) = V  (2   ,J,K)
!
	    TT (1 ,J,K) = TT (NX-1,J,K)
	    TT (NX,J,K) = TT (2   ,J,K)
!
	  END DO
      ENDDO
!
	END IF
!- check ---------------------------------------------------------
!
!#if (defined SPMD)
!      write(filename,14) 'iaptrsf-p-',myrank,'.out'
!14    format(a10,i1,a4)
!      open (10,file=trim(filename))
!#else
!      open (10,file='iaptrsf-s.out')
!#endif
!
!      write(10,*) '------------------- u -----------------'
!      write(10,11) (j,(u(i,j,10),i=1,2),j=1,ny)
!
!      write(10,*) '------------------- v -----------------'
!      write(10,11) (j,(v(i,j,10),i=1,2),j=1,ny)
!
!      write(10,*) '------------------- tt -----------------'
!      write(10,11) (j,(tt(i,j,10),i=1,2),j=1,ny)
!
!      write(10,*) '------------------- p -----------------'
!      write(10,11) (j,(p(i,j),i=1,2),j=1,ny)
!
!      write(10,*) '------------------- pps -----------------'
!      write(10,11) (j,(pps(i,j),i=1,2),j=1,ny)
!
!11    format(1x,i5,2e30.20)
!
!
!      close (10)
!
!#if (defined SPMD)
!!      call mpi_finalize(j)
!#endif
!!      stop 'iaptrsf'
!
!!-----------------------------------------------------------------

!
	RETURN
	END
