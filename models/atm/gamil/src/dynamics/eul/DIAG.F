#include <misc.h>
#include <params.h>

!!(2003.11.29)
!!-------------------

	SUBROUTINE DIAG(U,V,P,PS,PLY,TT,US,VS,TS,H,HPS
     _               ,PMTOP,PSB,TBS,TB,CB,FAC,DSIG)
!
	IMPLICIT NONE

#include <PARADYN>

#if (defined SPMD)
#include <commpi.h>
      character*50 filename
#endif
!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	REAL*8
     _     U  (NX,NY,NL)  !  U=u*sqrt(Ps), input variable
     _    ,V  (NX,NY,NL)  !	 V=v*sqrt(Ps), input variable
     _    ,P  (NX,NY   )  !  P=sqrt(Ps)  , input variable
     _    ,PS (NX,NY   )  !	 Surface pressure, input variable
     _    ,PLY(NX,NY,NL)  !  PLY=p, Pressure in Sigma Layer, input variable
     _    ,TS (NX,NY,NL)  !  input variable, TS=T, TEMPERATURE
     _    ,TT (NX,NY,NL)  !  TT=R*T'*Sqrt(Ps)/CB, input variable
!                            where T'=T-TB, T is the temperatur,
!                            TBB is Temperature of the standard atmosphere
     _    ,TB (NX,NY,NL)  !  input variable,
     _    ,CB (NX,NY,NL)  !  CB=Sqrt(R*(KK*TBB-dTB/dlnp)), input variable,
!                            where, KK=R/Cp, R is a constant
     _    ,PMTOP          !  PMTOP=10hPa
     _    ,TBS(NX,NY   )  !  TBB at the surface, input constant
     _    ,PSB(NX,NY   )  !  PSB is the surface pressure of the standard
!					 atmosphere, input constant
     _    ,DSIG(NL     )  !  The vertical stepsizes, input constant
     _	  ,US(NX,NY,NL )  !  US = u, zonal wind,    output variable
     _    ,VS(NX,NY,NL )  !  VS = v, meridional wind,output variable
     _    ,HPS(NX,NY   )  !  the surface geopotential height deviation
     _	  ,H (NX,NY,NZ )  !  H=gz-HBB, gz is the geopotential height,
!                              HBB is the geopotential height of the standard atmopshere
     _    ,FAC(NX,NY,NZ)
     _    ,WK1,WK2        !  working variables
!
	INTEGER I,J,K
        integer begj,endj

#if (defined SPMD)
      begj = 2
      endj = ny-1
#else
      begj = 1
      endj = ny
#endif

      DO J=begj,endj
        DO K=1,NL
          DO I=2,NX-1
            WK1=2.0D0/(P(I,J)+P(I-1,J))
#if (defined SPMD)
            if ( (myrank.eq.0).and.(j.eq.endj) ) then
              WK2=0.0D0
            else
              WK2=2.0D0/(P(I,J)+P(I,J+1))
            endif
#else
            IF(J.LT.NY) THEN
              WK2=2.0D0/(P(I,J)+P(I,J+1))
            ELSE
              WK2=0.0D0
            ENDIF
#endif
            US(I,J,K)=WK1*U(I,J,K)
            VS(I,J,K)=WK2*V(I,J,K)
	    TS(I,J,K)=TT(I,J,K)*CB(I,J,K)/(P(I,J)*RD)+TB(I,J,K)
          ENDDO
	  US(1,J,K)=US(NX-1,J,K)
	  US(NX,J,K)=US(2,J,K)
	  VS(1,J,K)=VS(NX-1,J,K)
	  VS(NX,J,K)=VS(2,J,K)
	  TS(1,J,K)=TS(NX-1,J,K)
	  TS(NX,J,K)=TS(2,J,K)
        ENDDO
C
C     CALCULATING H(K=NZ)=R*TB(PSB)/PSB*(PS-PSB) AS WELL AS H(K.LT.NZ).
        DO I=2,NX-1
	  HPS(I,J   )=RD*TBS(I,J)/PSB(I,J)*(PS(I,J)+PMTOP-PSB(I,J))
          H  (I,J,NZ)=HPS(I,J)
        ENDDO
	HPS(1,J   ) =HPS(NX-1,J)
	HPS(NX,J   )=HPS(2,J)
	H  (1,J,NZ) =H(NX-1,J,NZ)
	H  (NX,J,NZ)=H(2,J,NZ)
        DO K=NL,1,-1
          DO I=1,NX
           H(I,J,K)=H(I,J,K+1)+DSIG(K)*P(I,J)*CB(I,J,K)
     &       /PLY(I,J,K)*TT(I,J,K)*.5*(FAC(I,J,K+1)+FAC(I,J,K))
          ENDDO
	  H(1,J,K)=H(NX-1,J,K)
	  H(NX,J,K)=H(2,J,K)
        ENDDO
      ENDDO
c
!- check ---------------------------------------------------------
!
!#if (defined SPMD)
!      write(filename,14) 'diag-p-',myrank,'.out'
!14    format(a7,i1,a4)
!      open (10,file=trim(filename))
!#else
!      open (10,file='diag-s.out')
!#endif
!
!      write(10,*) '------------------- us -----------------'
!      write(10,11) (j,(us(i,j,10),i=1,2),j=1,ny)
!
!      write(10,*) '------------------- vs -----------------'
!      write(10,11) (j,(vs(i,j,10),i=1,2),j=1,ny)
!
!      write(10,*) '------------------- ts -----------------'
!      write(10,11) (j,(ts(i,j,10),i=1,2),j=1,ny)
!
!      write(10,*) '------------------- h -----------------'
!      write(10,11) (j,(h(i,j,10),i=1,2),j=1,ny)
!
!11    format(1x,i5,2e30.20)
!
!
!      close (10)
!
!#if (defined SPMD)
!!      call mpi_finalize(j)
!#endif
!!      stop 'diag'
!
!-----------------------------------------------------------------

	RETURN
	END
