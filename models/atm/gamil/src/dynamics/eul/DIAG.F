#include <misc.h>
#include <params.h>

!!(2003.11.29)
!!-------------------
	SUBROUTINE DIAG(U,V,P,PS,PLY,TT,US,VS,TS,H,HPS
     _               ,PMTOP,PSB,TSB,TB,CB,FAC,DSIG)
!
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
	IMPLICIT NONE

#include <PARADYN>

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	REAL*8
     _     U  (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  U=U*sqrt(PS), input variable
     _    ,V  (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !	 V=V*sqrt(PS), input variable
     _    ,P  (beglonex:endlonex,beglatexdyn:endlatexdyn   )  !  P=sqrt(PS)  , input variable
     _    ,PS (beglonex:endlonex,beglatexdyn:endlatexdyn   )  !	 Surface pressure, input variable
     _    ,PLY(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  PLY=p, Pressure in Sigma Layer, input variable
     _    ,TS (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  input variable, TS=T, TEMPERATURE
     _    ,TT (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  TT=R*T'*Sqrt(PS)/CB, input variable
!                            where T'=T-TB, T is the temperatur,
!                            TBB is Temperature of the standard atmosphere
     _    ,TB (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  input variable,
     _    ,CB (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  CB=Sqrt(R*(KK*TBB-dTB/dlnp)), input variable,
!                            where, KK=R/Cp, R is a constant
     _    ,PMTOP          !  PMTOP=10hPa
     _    ,TSB(beglonex:endlonex,beglatexdyn:endlatexdyn   )  !  TBB at the surface, input constant
     _    ,PSB(beglonex:endlonex,beglatexdyn:endlatexdyn   )  !  PSB is the surface pressure of the standard
!					 atmosphere, input constant
     _    ,DSIG(NL     )  !  The vertical stepsizes, input constant
     _	  ,US(beglonex:endlonex,beglatexdyn:endlatexdyn,NL )  !  US = U, zonal wind,    output variable
     _    ,VS(beglonex:endlonex,beglatexdyn:endlatexdyn,NL )  !  VS = V, meridional wind,output variable
     _    ,HPS(beglonex:endlonex,beglatexdyn:endlatexdyn   )  !  the surface geopotential height deviation
     _	  ,H (beglonex:endlonex,beglatexdyn:endlatexdyn,NZ )  !  H=gz-HBB, gz is the geopotential height,
!                              HBB is the geopotential height of the standard atmopshere
     _    ,FAC(beglonex:endlonex,beglatexdyn:endlatexdyn,NZ)
     _    ,WK1,WK2(beglonex:endlonex,beglatexdyn:endlatexdyn)        !  working variables
!
	INTEGER I,J,K

      call t_startf("DIAG COMP")

!$OMP PARALLEL DO PRIVATE(I,J,K,WK1)
      DO J=jbeg0,jend0
        DO I=ibeg1,iend1
	      HPS(I,J   )=RD*TSB(I,J)/PSB(I,J)*(PS(I,J)+PMTOP-PSB(I,J))
          H  (I,J,NZ)=HPS(I,J)
          IF (J .eq. PLAT) THEN
            WK2(I,J)=0.0D0
          ELSE
            WK2(I,J)=2.0D0/(P(I,J)+P(I,J+1))
          ENDIF
        ENDDO

        DO K=NL,1,-1
          DO I=ibeg1,iend1
            WK1=2.0D0/(P(I,J)+P(I-1,J))
            US(I,J,K)=WK1*U(I,J,K)
            VS(I,J,K)=WK2(I,J)*V(I,J,K)
            TS(I,J,K)=TT(I,J,K)*CB(I,J,K)/(P(I,J)*RD)+TB(I,J,K)
            H(I,J,K)=H(I,J,K+1)+DSIG(K)*P(I,J)*CB(I,J,K)
     &       /PLY(I,J,K)*TT(I,J,K)*.5*(FAC(I,J,K+1)+FAC(I,J,K))
          ENDDO
        ENDDO
      ENDDO

      call t_stopf("DIAG COMP")

      call gamil_arrays_comm(COMM_TO_LEFT,1,US(:,beglatexdyn,1),TS(:,beglatexdyn,1),
     &                       VS(:,beglatexdyn,1),HPS(:,beglatexdyn),H(:,beglatexdyn,1))
      call gamil_arrays_comm(COMM_TO_RIGHT,1,US(:,beglatexdyn,1),TS(:,beglatexdyn,1),
     &                       VS(:,beglatexdyn,1),HPS(:,beglatexdyn),H(:,beglatexdyn,1))


	RETURN
	END
