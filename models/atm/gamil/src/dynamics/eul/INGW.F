#include <misc.h>
#include <params.h>

      SUBROUTINE INGW(UK,P,TTK,CB0,DSIG,OUX,DT,SINU,WTGU,DU0,DTT0,ITN,EE)
!
      IMPLICIT NONE
!
!	This subroutine is to calculate the tendency of the wind and the temperature:
!     DU, DV, DTT
!

#include <PARADYN>

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	  REAL*8
     _       U (NX,NL )  !  input variable, zonal wind velocity
     _      ,UP(NX,NL )  !
     _      ,UK(NX,NL )  !
     _      ,P  (NX   )  !  input variable, P  = sqrt(Ps)
     _      ,TT (NX,NL)  !  input variable, TT=R*T'*Sqrt(Ps)/CB,
!                              where T'=T-TB, T is the temperatur,
!                              TBB	is Temperature of the standard atmosphere
     _	  ,TTK(NX,NL )  !  input variable,
     _	  ,TTP(NX,NL )  !  input variable,
!							 H=gz-HBB, gz is the geopotential height,
!                              HBB is the geopotential height of the standard atmopshere
     _      ,CB0(NX,NL)     !  input variable, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
!                           !  CB0=CB*P/PLY
!                              where, KK=R/Cp, R is a constant
     _      ,DSIG(NL     )  !  input constant, the vertical stepsizes
     _      ,OUX            !  input constant, OUX=1/(RAD*SINU*DX*MDJ)
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
!
      INTEGER I,K,I1,KWB,JWB,ITN
!
      REAL*8
     _       DU0 (NX,NL)	!  input variables
     _      ,DTT0(NX,NL)  !  input variables
     _      ,DU1 (NX,NL)	!  output variables
     _      ,DTT1(NX,NL)  !  output variables
     _      ,DU2 (NX,NL)	!  output variables
     _      ,DTT2(NX,NL)  !  output variables
     _      ,DU  (NX,NL)	!  output variables
     _      ,DTT (NX,NL)  !  output variables
     _      ,DUX (NX,NL)	!  output variables
     _      ,DTTX(NX,NL)  !  output variables
!
	REAL*8 DTGW,DT2,DT,DJ,DS,EE,Y1,Y2,Y3,SINU,WTGU
!
	DTGW=DT/FLOAT(ITN)
	DT2=0.5*DTGW
	I1=NX-1
!
	DO K=1,NL
	DO I=1,NX
	   U (I,K)=UK (I,K)
	   TT(I,K)=TTK(I,K)
!
	   DU (I,K)=0.0
	   DTT(I,K)=0.0
!
	   DU1 (I,K)=0.0
	   DTT1(I,K)=0.0
!
	   DU2 (I,K)=0.0
	   DTT2(I,K)=0.0
	ENDDO
	ENDDO
!
	DO KWB=1,ITN
!
	DO K=1,NL
	DO I=1,NX
	   UP (I,K)=U (I,K)
	   TTP(I,K)=TT(I,K)
	ENDDO
	ENDDO
!
	DO JWB=1,2
!
	CALL DIFUTX(U,P,TT,CB0,DSIG,OUX,DUX,DTTX)
!
	DO K=1,NL
	DO I=1,NX
	   U (I,K)=UP (I,K)+DT2*(DUX (I,K)+DU0 (I,K))
	   TT(I,K)=TTP(I,K)+DT2*(DTTX(I,K)+DTT0(I,K))
	ENDDO
	ENDDO
!
	END DO
!
	DO K=1,NL
	DO I=1,NX
	   DU1 (I,K)=DU1 (I,K)+DUX (I,K)
	   DTT1(I,K)=DTT1(I,K)+DTTX(I,K)
	ENDDO
	ENDDO
!
	CALL DIFUTX(U,P,TT,CB0,DSIG,OUX,DUX,DTTX)
!
	DO K=1,NL
	DO I=1,NX
	   U (I,K)=UP (I,K)+DTGW*(DUX (I,K)+DU0 (I,K))
	   TT(I,K)=TTP(I,K)+DTGW*(DTTX(I,K)+DTT0(I,K))
!
	   DU2 (I,K)=DU2 (I,K)+DU (I,K)
	   DTT2(I,K)=DTT2(I,K)+DTT(I,K)
!
	   DU  (I,K)=DU  (I,K)+DUX (I,K)
	   DTT (I,K)=DTT (I,K)+DTTX(I,K)
	ENDDO
	ENDDO
!
	END DO
!
	Y1=0.0
	Y2=0.0
	Y3=0.0
        DJ=SINU/WTGU
	DO K=1,NL
	   DS=DJ*DSIG(K)
	   DO I=2,I1
	      DU  (I,K)=DU  (I,K)/FLOAT(ITN)
	      DTT (I,K)=DTT (I,K)/FLOAT(ITN)
!
		  Y1=Y1+(DU (I,K)*DU2 (I,K)+DTT(I,K)*DTT2(I,K))*DS
		  Y2=Y2+(DU (I,K)*DU1 (I,K)+DTT(I,K)*DTT1(I,K))*DS
		  Y3=Y3+(DU (I,K)*DU0 (I,K)+DTT(I,K)*DTT0(I,K))*DS
!
	      DU0 (I,K)=DU0 (I,K)+DU (I,K)
	      DTT0(I,K)=DTT0(I,K)+DTT(I,K)
	   ENDDO
!
	   DU0(1,K)=DU0(I1,K)
	   DU0(NX,K)=DU0(2,K)
!
	   DTT0(1,K)=DTT0(I1,K)
	   DTT0(NX,K)=DTT0(2,K)
!
	ENDDO
c	EE=EE+(Y1+Y1+Y2)/FLOAT(ITN*ITN)+Y3
	EE=(Y1+Y1+Y2)/FLOAT(ITN*ITN)+Y3
!
	RETURN
	END
!
!!!#include <misc.h>
!!!#include <params.h>

!!(2004.02.15)
!***********************************************************************
!
      SUBROUTINE MINUS_INGW(U,P,TT,CB0,DSIG,OUX,DU,DTT,NI)
!
#if ( defined SPMD )
       use mpishorthand, only: mpicom
#endif
!
      IMPLICIT NONE
!
!	This subroutine is to deduct the inner gravity waves from the tendency
!     of the zonal wind and the temperature: DU, DTT
!

#include <PARADYN>

#if (defined SPMD)
#include <commpi.h>
#endif

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	  REAL*8
     _       U (NX,NY,NL )  !  input array, zonal wind velocity
     _      ,UJ (NX,   NL) !  working array,
     _      ,P  (NX,NY   )  !  input array, P  = sqrt(Ps)
     _      ,TT (NX,NY,NL)  !  input array, TT=R*T'*Sqrt(Ps)/CB,
     _      ,TJ (NX,    NL)  !  input array, TT=R*T'*Sqrt(Ps)/CB,
!                              where T'=T-TB, T is the temperatur,
!                              TBB	is Temperature of the standard atmosphere
     _      ,CBJ(NX   ,NL)  !  input array, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
     _      ,CB0(NX,NY,NL)  !  input array, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
!                           !  CB0=CB*P/PLY
!                              where, KK=R/Cp, R is a constant
     _      ,DSIG(NL     )  !  input constant array, the vertical stepsizes
     _      ,OUX(NY)        !  input constant array, OUX=1/(RAD*SINU*DX*MDJ)
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
!
      REAL*8
     _       DU (NX,NY,NL)   !  input & output array
     _      ,DTT (NX,NY,NL)  !  input & output array
     _      ,DUX (NX   ,NL)  !  working array: inner gravity waves in DU
     _      ,DTTX(NX   ,NL)  !  working array: inner gravity waves in DTT
!
      INTEGER
     _       I               !  working variable
     _      ,J			   !  working variable
     _      ,K			   !  working variable
     _      ,NI(NY)		   !  times of moving length of inner gravity waves to
!						   !  the zonal gridsize
#if (defined SPMD)
	DO J=jbeg,jend
#else
	DO J=2,NY-1
#endif
!
	IF (NI(J).GT.1) THEN
!
	   DO K=1,NL
	   DO I=1,NX
	      UJ (I,K)=U  (I,J,K)
	      TJ (I,K)=TT (I,J,K)
	      CBJ(I,K)=CB0(I,J,K)
	   END DO
	   END DO
!
	   CALL DIFUTX(UJ,P(1,J),TJ,CBJ,DSIG,OUX(J),DUX,DTTX)
!
	   DO K=1,NL
	   DO I=1,NX
	      DU (I,J,K)=DU (I,J,K)-DUX (I,K)
	      DTT(I,J,K)=DTT(I,J,K)-DTTX(I,K)
	   END DO
	   END DO
!
	ENDIF
!
	END DO
!
	RETURN
	END
!
!!#include <misc.h>
!!#include <params.h>

!!(wb 2004.02.15)
!***********************************************************************
!
      SUBROUTINE PLUS_INGW(UK,P,TTK,CB0,DSIG,OUX,DT,SINU,WTGU,DU0,DTT0,EE,NI)
!
#if ( defined SPMD )
       use mpishorthand, only: mpicom
#endif
!
      IMPLICIT NONE
!
!	This subroutine is to calculate the tendency of the wind and the temperature:
!     DU, DV, DTT
!

#include <PARADYN>

#if (defined SPMD)
#include <commpi.h>
#include <mpif.h>

      integer isend1,irecv1,isend2,irecv2,jj
      integer istatus(mpi_status_size),nii(ny)
      real*8 EEJ(NY),FJ(NY)
#else
      real*8 EK
#endif

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	  REAL*8
     _       UK (NX,NY,NL )  !
     _      ,UKJ(NX   ,NL )  !
     _      ,P  (NX,NY   )  !  input variable, P  = sqrt(Ps)
     _	    ,TTK(NX,NY,NL )  !  input variable,
     _	    ,TKJ(NX   ,NL )  !  input variable,
     _      ,CBJ(NX   ,NL)  !  input variable, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
     _      ,CB0(NX,NY,NL)  !  input variable, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
!                           !  CB0=CB*P/PLY
!                              where, KK=R/Cp, R is a constant
     _      ,DSIG(NL     )  !  input constant, the vertical stepsizes
     _      ,OUX(NY)        !  input constant, OUX=1/(RAD*SINU*DX*MDJ)
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
!
      REAL*8
     _       DU0 (NX,NY,NL)	  !  input variables
     _      ,DTT0(NX,NY,NL)     !  input variables
     _      ,DUJ (NX   ,NL)	  !  working variables
     _      ,DTJ (NX   ,NL)     !  working variables
!
	REAL*8 DT,EE,EJ,SINU(NY),WTGU(NY)
!
	INTEGER I,J,K,MM,NI(NY)
!
	EJ=0.0
!
#if (!defined SPMD)
	DO J=2,NY-1
#else
        DO J=jbeg,jend
#endif
!
	MM=NI(J)
	IF (MM.GT.1) THEN
!
	DO I=1,NX
	   DO K=1,NL
	      UKJ(I,K)=UK (I,J,K)
	      TKJ(I,K)=TTK(I,J,K)
	      CBJ(I,K)=CB0(I,J,K)
	      DUJ(I,K)=DU0 (I,J,K)
	      DTJ(I,K)=DTT0(I,J,K)
	   ENDDO
	ENDDO

#if (defined SPMD)
	CALL INGW(UKJ,P(1,J),TKJ,CBJ,DSIG,OUX(J),DT
     _         ,SINU(J),WTGU(J),DUJ,DTJ,MM,FJ(J))
#else
	CALL INGW(UKJ,P(1,J),TKJ,CBJ,DSIG,OUX(J),DT
     _         ,SINU(J),WTGU(J),DUJ,DTJ,MM,EK)
	EJ=EJ+EK
#endif
!
	DO I=1,NX
	   DO K=1,NL
	      DU0 (I,J,K)=DUJ(I,K)
	      DTT0(I,J,K)=DTJ(I,K)
	   ENDDO
	ENDDO
!
	END IF
!
	ENDDO
!
#if (defined SPMD)
      if (myrank.ne.0) then
         call mpi_isend(FJ,ny,mpi_double_precision,0,1,mpicom,isend1,ierr)
         call mpi_isend(NI,ny,mpi_integer,0,2,mpicom,isend2,ierr)
         call mpi_wait(isend1,istatus,ierr)
         call mpi_wait(isend2,istatus,ierr)
      else
       DO J=2,NY-1
          IF(NI(J).GT.1) EJ=EJ+FJ(J)
       ENDDO
       do jj=1,nprocessor-1
         call mpi_irecv(EEJ,ny,mpi_double_precision,jj,1,mpicom,irecv1,ierr)
         call mpi_irecv(NII,ny,mpi_integer,jj,2,mpicom,irecv2,ierr)
         call mpi_wait(irecv1,istatus,ierr)
         call mpi_wait(irecv2,istatus,ierr)
	 DO J=2,NY-1
            IF(NII(J).GT.1) EJ=EJ+EEJ(J)
	 ENDDO
       enddo
      endif
      EE=EE+EJ
      call mpi_bcast(EE,1,mpi_double_precision,0,mpicom,ierr)
#else
	EE=EE+EJ
#endif
!
	RETURN
	END
