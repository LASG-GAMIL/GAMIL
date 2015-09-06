#include <misc.h>
#include <params.h>

      SUBROUTINE INGW(UK,P,TTK,CB0,DSIG,OUX,DT,SINU,WTGU,DU0,DTT0,ITN,EE)
      use mpi_gamil
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
     _       U (beglonex:endlonex,NL )  !  input variable, zonal wind velocity
     _      ,UP(beglonex:endlonex,NL )  !
     _      ,UK(beglonex:endlonex,NL )  !
     _      ,P  (beglonex:endlonex   )  !  input variable, P  = sqrt(Ps)
     _      ,TT (beglonex:endlonex,NL)  !  input variable, TT=R*T'*Sqrt(Ps)/CB,
!                              where T'=T-TB, T is the temperatur,
!                              TBB	is Temperature of the standard atmosphere
     _	  ,TTK(beglonex:endlonex,NL )  !  input variable,
     _	  ,TTP(beglonex:endlonex,NL )  !  input variable,
!							 H=gz-HBB, gz is the geopotential height,
!                              HBB is the geopotential height of the standard atmopshere
     _      ,CB0(beglonex:endlonex,NL)     !  input variable, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
!                           !  CB0=CB*P/PLY
!                              where, KK=R/Cp, R is a constant
     _      ,DSIG(NL     )  !  input constant, the vertical stepsizes
     _      ,OUX            !  input constant, OUX=1/(RAD*SINU*DX*MDJ)
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
!
      INTEGER I,K,KWB,JWB,ITN
!
      REAL*8
     _       DU0 (beglonex:endlonex,NL)	!  input variables
     _      ,DTT0(beglonex:endlonex,NL)  !  input variables
     _      ,DU1 (beglonex:endlonex,NL)	!  output variables
     _      ,DTT1(beglonex:endlonex,NL)  !  output variables
     _      ,DU2 (beglonex:endlonex,NL)	!  output variables
     _      ,DTT2(beglonex:endlonex,NL)  !  output variables
     _      ,DU  (beglonex:endlonex,NL)	!  output variables
     _      ,DTT (beglonex:endlonex,NL)  !  output variables
     _      ,DUX (beglonex:endlonex,NL)	!  output variables
     _      ,DTTX(beglonex:endlonex,NL)  !  output variables
!
	REAL*8 DTGW,DT2,DT,DJ,DS,SINU,WTGU
        REAL*16 Y1,Y2,Y3,EE,Y1_ARR(NL),Y2_ARR(NL),Y3_ARR(NL)
!
	DTGW=DT/FLOAT(ITN)
	DT2=0.5*DTGW
!
!$OMP PARALLEL DO PRIVATE (I,K)
	DO K=1,NL
	DO I=beglonex,endlonex
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
!$OMP PARALLEL DO PRIVATE (I,K)
	DO K=1,NL
	DO I=beglonex,endlonex
	   UP (I,K)=U (I,K)
	   TTP(I,K)=TT(I,K)
	ENDDO
	ENDDO
!
	DO JWB=1,2
!
	CALL DIFUTX(U,P,TT,CB0,DSIG,OUX,DUX,DTTX)
!
!$OMP PARALLEL DO PRIVATE (I,K)
	DO K=1,NL
	DO I=beglonex,endlonex
	   U (I,K)=UP (I,K)+DT2*(DUX (I,K)+DU0 (I,K))
	   TT(I,K)=TTP(I,K)+DT2*(DTTX(I,K)+DTT0(I,K))
	ENDDO
	ENDDO
!
	END DO
!
!$OMP PARALLEL DO PRIVATE (I,K)
	DO K=1,NL
	DO I=beglonex,endlonex
	   DU1 (I,K)=DU1 (I,K)+DUX (I,K)
	   DTT1(I,K)=DTT1(I,K)+DTTX(I,K)
	ENDDO
	ENDDO
!
	CALL DIFUTX(U,P,TT,CB0,DSIG,OUX,DUX,DTTX)
!
!$OMP PARALLEL DO PRIVATE (I,K)
	DO K=1,NL
	DO I=beglonex,endlonex
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
        DJ=SINU/WTGU
!$OMP PARALLEL DO PRIVATE (I,K,Y1,Y2,Y3,DS)
	DO K=1,NL
	   Y1=0.0
	   Y2=0.0
	   Y3=0.0
	   DS=DJ*DSIG(K)
	   DO I=ibeg1,iend1
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
           Y1_ARR(K)=Y1
           Y2_ARR(K)=Y2
           Y3_ARR(K)=Y3
	ENDDO

       call register_comm_array(beglonex,endlonex,1,1,1,NL,1,1,DU0(:,1)) 
       call register_comm_array(beglonex,endlonex,1,1,1,NL,1,1,DTT0(:,1)) 
       call gamil_arrays_comm(COMM_TO_LEFT,1,DU0(:,1),DTT0(:,1))
       call gamil_arrays_comm(COMM_TO_RIGHT,1,DU0(:,1),DTT0(:,1))
       call remove_comm_array(DTT0(:,1)) 
       call remove_comm_array(DU0(:,1)) 

	Y1=0.0
	Y2=0.0
	Y3=0.0
	DO K=1,NL
           Y1 = Y1+Y1_ARR(K)
           Y2 = Y2+Y2_ARR(K)
           Y3 = Y3+Y3_ARR(K)
	ENDDO

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
      SUBROUTINE MINUS_INGW(UU,P,TT,CB0,DSIG,OUX,DU,DTT,NIGW)
! 
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
!
      IMPLICIT NONE
!
!	This subroutine is to deduct the inner gravity waves from the tendency
!     of the zonal wind and the temperature: DU, DTT
!

#include <PARADYN>

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	  REAL*8
     _       UU (beglonex:endlonex,beglatexdyn:endlatexdyn,NL )  !  input array, zonal wind velocity
     _      ,UJ (beglonex:endlonex,   NL) !  working array,
     _      ,P  (beglonex:endlonex,beglatexdyn:endlatexdyn)  !  input array, P  = sqrt(Ps)
     _      ,TT (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  input array, TT=R*T'*Sqrt(Ps)/CB,
     _      ,TJ (beglonex:endlonex,    NL)  !  input array, TT=R*T'*Sqrt(Ps)/CB,
!                              where T'=T-TB, T is the temperatur,
!                              TBB	is Temperature of the standard atmosphere
     _      ,CBJ(beglonex:endlonex   ,NL)  !  input array, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
     _      ,CB0(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  input array, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
!                           !  CB0=CB*P/PLY
!                              where, KK=R/Cp, R is a constant
     _      ,DSIG(NL     )  !  input constant array, the vertical stepsizes
     _      ,OUX(beglatexdyn:endlatexdyn)        !  input constant array, OUX=1/(RAD*SINU*DX*MDJ)
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
!
      REAL*8
     _       DU (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)   !  input & output array
     _      ,DTT (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  input & output array
     _      ,DUX (beglonex:endlonex   ,NL)  !  working array: inner gravity waves in DU
     _      ,DTTX(beglonex:endlonex   ,NL)  !  working array: inner gravity waves in DTT
!
      INTEGER
     _       I               !  working variable
     _      ,J			   !  working variable
     _      ,K			   !  working variable
     _      ,NIGW(beglatexdyn:endlatexdyn)		   !  times of moving length of inner gravity waves to
!						   !  the zonal gridsize
!

       call gamil_arrays_comm(COMM_TO_LEFT,1,p(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT,1,p(:,beglatexdyn))

       DO j=jbeg1, jend1
	IF (NIGW(J).GT.1) THEN
!
!$OMP PARALLEL DO PRIVATE (I,K)
	   DO K=1,NL
	   DO I=beglonex,endlonex
	      UJ (I,K)=UU  (I,J,K)
	      TJ (I,K)=TT (I,J,K)
	      CBJ(I,K)=CB0(I,J,K)
	   END DO
	   END DO
!
	   CALL DIFUTX(UJ,P(beglonex,J),TJ,CBJ,DSIG,OUX(J),DUX,DTTX)
!
!$OMP PARALLEL DO PRIVATE (I,K)
	   DO K=1,NL
	   DO I=beglonex,endlonex
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
      SUBROUTINE PLUS_INGW(UK,P,TTK,CB0,DSIG,OUX,DT,SINU,WTGU,DU0,DTT0,EE,NIGW)
!
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
!
      IMPLICIT NONE

#include <PARADYN>
!
!	This subroutine is to calculate the tendency of the wind and the temperature:
!     DU, DV, DTT
!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	  REAL*8
     _       UK (beglonex:endlonex,beglatexdyn:endlatexdyn,NL )  !
     _      ,UKJ(beglonex:endlonex   ,NL )  !
     _      ,P  (beglonex:endlonex,beglatexdyn:endlatexdyn   )  !  input variable, P  = sqrt(Ps)
     _	    ,TTK(beglonex:endlonex,beglatexdyn:endlatexdyn,NL )  !  input variable,
     _	    ,TKJ(beglonex:endlonex   ,NL )  !  input variable,
     _      ,CBJ(beglonex:endlonex   ,NL)  !  input variable, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
     _      ,CB0(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  input variable, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
!                           !  CB0=CB*P/PLY
!                              where, KK=R/Cp, R is a constant
     _      ,DSIG(NL     )  !  input constant, the vertical stepsizes
     _      ,OUX(beglatexdyn:endlatexdyn)        !  input constant, OUX=1/(RAD*SINU*DX*MDJ)
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
!
      REAL*8
     _       DU0 (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)	  !  input variables
     _      ,DTT0(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)     !  input variables
     _      ,DUJ (beglonex:endlonex   ,NL)	  !  working variables
     _      ,DTJ (beglonex:endlonex   ,NL)     !  working variables
      real*16 EK(PLAT)
!
      REAL*8 DT,EE,SINU(beglatexdyn:endlatexdyn),WTGU(beglatexdyn:endlatexdyn)
      REAL*16 EJ
!
      INTEGER I,J,K,MM,NIGW(beglatexdyn:endlatexdyn)
!
!
       call gamil_arrays_comm(COMM_TO_LEFT,1,p(:,beglatexdyn))
       call gamil_arrays_comm(COMM_TO_RIGHT,1,p(:,beglatexdyn))

      DO J=jbeg1,jend1
!
      MM=NIGW(J)
      IF (MM.GT.1) THEN
!
!$OMP PARALLEL DO PRIVATE (I,K)
      DO K=1,NL
        DO I=beglonex,endlonex
          UKJ(I,K)=UK (I,J,K)
          TKJ(I,K)=TTK(I,J,K)
          CBJ(I,K)=CB0(I,J,K)
          DUJ(I,K)=DU0 (I,J,K)
          DTJ(I,K)=DTT0(I,J,K)
        ENDDO
      ENDDO


      CALL INGW(UKJ,P(beglonex,J),TKJ,CBJ,DSIG,OUX(J),DT
     _         ,SINU(J),WTGU(J),DUJ,DTJ,MM,EK(J))
!
!$OMP PARALLEL DO PRIVATE (I,K)
      DO K=1,NL
        DO I=beglonex,endlonex
          DU0 (I,J,K)=DUJ(I,K)
          DTT0(I,J,K)=DTJ(I,K)
        ENDDO
      ENDDO
!
      END IF
!
      ENDDO
!
      EJ=0.0
      DO J=jbeg1,jend1
      IF (NIGW(J) .GT. 1) THEN
         EJ=EJ+EK(J)
      ENDIF
      ENDDO
      call gamil_all_reduce(EJ)
    
      EE=EE+EJ
!
      RETURN
      END
