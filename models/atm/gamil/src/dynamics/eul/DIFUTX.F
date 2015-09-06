      SUBROUTINE DIFUTX(U,P,TT,CB0,DSIG,OUX,DU,DTT)
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
     _      ,UU(beglonex:endlonex,NL )  !
     _      ,UZ(beglonex:endlonex,NZ )  !
     _      ,P  (beglonex:endlonex   )  !  input variable, P  = sqrt(Ps)
     _      ,PP (beglonex:endlonex   )  !  input variable, P  = sqrt(Ps)
     _      ,TT (beglonex:endlonex,NL)  !  input variable, TT=R*T'*Sqrt(Ps)/CB,
!                              where T'=T-TB, T is the temperatur,
!                              TBB	is Temperature of the standard atmosphere
     _	  ,HZ(beglonex:endlonex,NZ )  !  input variable,
     _	  ,HH(beglonex:endlonex,NZ )  !  input variable,
!							 H=gz-HBB, gz is the geopotential height,
!                              HBB is the geopotential height of the standard atmopshere
     _      ,CB0(beglonex:endlonex,NL)  !  input variable, CB=Sqrt(R*(KK*TBB-dTB/dlnp)),
!                           !  CB0=CB*P/PLY
!                              where, KK=R/Cp, R is a constant
     _      ,DSIG(NL  )  !  input constant, the vertical stepsizes
     _      ,OUX         !  input constant, OUX=1/(RAD*SINU*DX*MDJ)
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
      REAL*8
     _       DU(beglonex:endlonex,NL)	  !  output variables
     _      ,DTT(beglonex:endlonex,NL)  !  output variables
!
      INTEGER I,K,K1
!
!
!     CALCULATING THE INNER GRAVITY WAVE PART OF DU/DT AND DTT/DT.
!
	DO I=ibeg1,iend1
	   UZ(I,1)=0.0
	   HZ(I,NZ)=0.0
       PP(I)=0.5*(P(I)+P(I-1))
	ENDDO
!
!
!$OMP PARALLEL DO PRIVATE (I,K,K1)
	DO I=ibeg1,iend1
	   DO K=1,NL
	      K1=K+1
	      UZ(I,K1)=UZ(I,K)+U(I,K)*DSIG(K)
	      UU(I,K)=0.5*(UZ(I,K)+UZ(I,K1))*PP(I)
	   ENDDO
	   DO K=NL,1,-1
	      K1=K+1
	      HZ(I,K)=HZ(I,K1)+TT(I,K)*CB0(I,K)*DSIG(K)
	      HH(I,K)=0.5*(HZ(I,K)+HZ(I,K1))
	   ENDDO
	ENDDO
!
      call register_comm_array(beglonex,endlonex,1,1,1,NL,1,1,UU(:,1)) 
      call register_comm_array(beglonex,endlonex,1,1,1,NL,1,1,HH(:,1)) 
      call register_comm_array(beglonex,endlonex,1,1,1,NL,1,1,DU(:,1)) 
      call register_comm_array(beglonex,endlonex,1,1,1,NL,1,1,DTT(:,1)) 
      call gamil_arrays_comm(COMM_TO_LEFT,1,UU(:,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,HH(:,1))
!
!$OMP PARALLEL DO PRIVATE (I,K,K1)
      DO K=1,NL
         DO I=ibeg1,iend1
            DTT(I,K)=-CB0(I,K)*OUX*(UU(I+1,K)-UU(I,K))
            DU(I,K)=-OUX*PP(I)*(HH(I,K)-HH(I-1,K))
         ENDDO
      ENDDO
!
       call gamil_arrays_comm(COMM_TO_LEFT,1,DU(:,1),DTT(:,1))
       call gamil_arrays_comm(COMM_TO_RIGHT,1,DU(:,1),DTT(:,1))

      call remove_comm_array(DTT(:,1)) 
      call remove_comm_array(DU(:,1)) 
      call remove_comm_array(HH(:,1)) 
      call remove_comm_array(UU(:,1)) 
	RETURN
	END
