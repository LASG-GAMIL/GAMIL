#include <misc.h>
#include <params.h>

      SUBROUTINE IAPTRSF(US,VS,WS,PS,TS,U,V,P,TT,TB,CB,KP)
!
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil

      IMPLICIT NONE

#include <PARADYN>

!
!     The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	REAL*8
     _	     US(beglonex:endlonex,beglatexdyn:endlatexdyn,NL )  !  US = U, zonal wind,    input variable
     _      ,VS(beglonex:endlonex,beglatexdyn:endlatexdyn,NL )  !  VS = V, meridional wind, input variable
     _      ,WS(beglonex:endlonex,beglatexdyn:endlatexdyn,NZ )  !  vertical velocity, input variable
     _      ,PS(beglonex:endlonex,beglatexdyn:endlatexdyn)  !  Surface pressure, input variable
     _      ,TS(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  temperature, input variable
     _      ,TB(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  the temperature of the standard atmosphere a
!                           !  at the sigma layers, input constant
     _      ,CB (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  CB=Sqrt(R*(KK*TBB-dTB/dlnp)), input variable,
!                           !  where, KK=R/Cp, R is a constant
     _      ,U  (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  U=U*sqrt(PS), output variable
     _      ,V  (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !	 V=V*sqrt(PS), output variable
     _      ,P  (beglonex:endlonex,beglatexdyn:endlatexdyn   )  !  P=sqrt(PS)  , output variable
     _      ,TT (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  TT=R*T'*Sqrt(PS)/CB, output variable
!                           !  where T'=T-TB, T is the temperatur,
     _      ,PX,PY          !  working variables
!
	INTEGER I,J,K,KP


!$OMP PARALLEL DO PRIVATE (I, J)
        DO J=jbeg0,jend0
	  DO I=ibeg1,iend1
	     P(I,J)=SQRT(PS(I,J))
	  END DO
        END DO

      call gamil_arrays_comm(COMM_TO_LEFT,1,p(:,beglatexdyn),PS(:,beglatexdyn))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,p(:,beglatexdyn),PS(:,beglatexdyn))
      call gamil_arrays_comm(COMM_TO_TOP,1,p(:,beglatexdyn))


      IF (KP.EQ.1) THEN

!$OMP PARALLEL DO PRIVATE (I, J)
        DO J=jbeg0,jend0
	  DO I=beglonex,endlonex
             WS(I,J,1)=0.0D0
             WS(I,J,NZ)=0.0D0
	  END DO
	END DO
      ENDIF
        

!$OMP PARALLEL DO PRIVATE (I, J, K, PY, PX)
      DO J=jbeg0,jend0
        DO I=ibeg1,iend1
!
!         CALCULATING THE AVERAGE VALUE OF P AT V-GRID.
!
          IF(J.LT.PLAT) THEN
            PY=0.5D0*(P(I,J)+P(I,J+1))
          ELSE
            PY=0.0D0
          ENDIF
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
!         CALCULATING TT=(T-TB)*P*RD/CB.
!
	     TT (I,J,K)=(TS(I,J,K)-TB(I,J,K))*P(I,J)*RD/CB(I,J,K)
!
          ENDDO
        ENDDO
      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,US(:,beglatexdyn,1),VS(:,beglatexdyn,1),
     &                       TS(:,beglatexdyn,1),U(:,beglatexdyn,1),V(:,beglatexdyn,1),TT(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,US(:,beglatexdyn,1),VS(:,beglatexdyn,1),
     &                       TS(:,beglatexdyn,1),U(:,beglatexdyn,1),V(:,beglatexdyn,1),TT(:,beglatexdyn,1))

!
	RETURN
	END
