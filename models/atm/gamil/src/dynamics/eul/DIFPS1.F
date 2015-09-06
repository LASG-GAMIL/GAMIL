#include <misc.h>
#include <params.h>

!!(2003.11.10)
!***********************************************************************
!
      SUBROUTINE DIFPS1(U,V,P,PS,WS,DPS,DPSU
     _                ,DSIG,DY,OUX,OUY,SINV,MP1,MP2,MM1,MM2,WTGV)
     _


      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
      
	IMPLICIT NONE
!
!	This subroutine is to calculate the tendency of the surface prressure DPS,
!     the vertical velocity WS, the zonal wind US, meridional wind VS and the
!	departure of the geopotential height from the standard atmopshere

#include <PARADYN>

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	REAL*8
     _       U  (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !  U=u*sqrt(Ps), input variable
     _      ,V  (beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !	 V=v*sqrt(Ps), input variable
     _      ,P  (beglonex:endlonex,beglatexdyn:endlatexdyn   )  !  P=sqrt(Ps)  , input variable
     _      ,PS (beglonex:endlonex,beglatexdyn:endlatexdyn   )  !	 Surface pressure, input variable
     _      ,DSIG(NL     )  !  The vertical stepsizes, input constant
     _      ,DY             !  The horizontal stepsize in meridional direction
!                              input constant
     _      ,SINV(beglatexdyn:endlatexdyn)       !  sin(theta) at half grid j+1/2, input constant
     _      ,OUX(beglatexdyn:endlatexdyn)        !  OUX=1/(RAD*SINU*DX*MDJ), input constant
!                              where, DX is the horizontal stepsize in zonal direction,
!                              MDJ is the leaping length of the central difference
!                              SINU is sin(theta) at intger grid j
     _      ,OUY(beglatexdyn:endlatexdyn)        !  OUY=1/(RAD*SINU*DY*WTGU), input constant,
!                              where WTGU is the weighting at the integer grid j
     _      ,WTGV(beglatexdyn:endlatexdyn)       !	 the weighting at the half grid j+1/2,
!                              input constant
     _      ,WS  (beglonex:endlonex,beglatexdyn:endlatexdyn,NZ) !  WS = w, vertical velocity, output variable
     _      ,DPS (beglonex:endlonex,beglatexdyn:endlatexdyn   ) !  the tendency of the surface pressure,
     _      ,DPSU(beglonex:endlonex,beglatexdyn:endlatexdyn   ) !  the tendency of the surface pressure,
!                              output variable
      REAL*8 WK0,WK1,WK2,WK3 !  working variables
     _      ,PXP,PXM,PYP	  !  working variables
     _      ,PYM,DPSP,WKQ   !  working variables
      REAL*8 W1(NZ),D1_8(NZ)
      REAL*16 D1_16(NZ),D1_16_TMP(8,NZ),TMP_SUM !  working variables
      REAL*8 pole_ps_1(1)
!
	INTEGER
     _       MM1(beglonex:endlonex,beglatexdyn:endlatexdyn)     !
     _      ,MP1(beglonex:endlonex,beglatexdyn:endlatexdyn)     !
     _      ,MM2(beglonex:endlonex,beglatexdyn:endlatexdyn)     !
     _      ,MP2(beglonex:endlonex,beglatexdyn:endlatexdyn)     !
!
      INTEGER I,J,K,IP1,IM1,IP2,IM2 !  working variables
      INTEGER J1,J2
      INTEGER icomm_request1, icomm_request2

      call gamil_arrays_comm(COMM_TO_TOP,1,p(:,beglatexdyn),request_id=icomm_request1)
      call gamil_arrays_comm(COMM_TO_BOT,1,p(:,beglatexdyn),v(:,beglatexdyn,1),request_id=icomm_request2)
      call wait_icomm_request(icomm_request1)
      call wait_icomm_request(icomm_request2)

C     DEALING WITH INTERNAL GRID POINTS.
C     CALCULATING FACTORS RELATIVE TO J IN DIVERGENCE FORMULA.
C

      call t_startf("DIFPS1 COMP")
!$OMP PARALLEL DO PRIVATE (I,J,K,WK1,WK2,WK3,IP1,IM1,IP2,IM2,PXP,PXM,
!$   &           PYP,PYM,DPSP,D1_8,W1,WK0,WKQ,J1)
      DO J=jbeg1,jend1
        DO I=ibeg1,iend1
          WS(I,J,1)=0.0D0
          WS(I,J,NZ)=0.0D0
        ENDDO
        
          WK1=0.5D0*OUX(J)
          WK2=0.5D0*OUY(J)*SINV(J)
          WK3=0.5D0*OUY(J)*SINV(J-1)
C     CALCULATING FACTORS RELATIVE TO I AND J IN DIVERGENCE FORMULA.
C
          DO I=ibeg1,iend1
            IP1=MP1(I,J)
            IM1=MM1(I,J)
            IP2=MP2(I,J)
            IM2=MM2(I,J)
            PXP=WK1*(P(IP1,J)+P(I,J))
            PXM=WK1*(P(I,J)+P(IM1,J))
            PYP=WK2*(P(I,J+1)+P(I,J))
            PYM=WK3*(P(I,J)+P(I,J-1))
C     CALCULATING DIVERGENCES AS WELL AS SUM OF THEM.
            DPSP=0.0D0
            DO K=1,NL
              D1_8(K)=PYP*V(I,J,K)-PYM*V(I,J-1,K)
              DPSP=DPSP-DSIG(K)*D1_8(K)
              D1_8(K)=D1_8(K)+PXP*U(IP2,J,K)-PXM*U(IM2,J,K)
            ENDDO
C
C     CALCULATING DPS/DT, DPS/DT AND D(SIGMA)/DT.
C
            DPS(I,J)=DPSP
	    DPSP=DPSP+DPSU(I,J)
            WKQ=1.D0/PS(I,J)
            DO K=2,NL
              WS(I,J,K)=WS(I,J,K-1)-DSIG(K-1)*WKQ*(DPSP+D1_8(K-1))
            ENDDO
          ENDDO
        ENDDO
      call t_stopf("DIFPS1 COMP")

      DO J=jbeg0,jend0 
      IF (J.EQ.1 .OR. J.EQ.PLAT) THEN
          call t_startf("DIFPS1 COMP")
          DO I=ibeg1,iend1
             WS(I,J,1)=0.0D0
             WS(I,J,NZ)=0.0D0
          ENDDO
 
          IF(J.EQ.1) THEN
             J1 = J+1
             J2 = J
          ELSE
             J1 = J-1
             J2 = J1
          ENDIF
C
C     FINDING DP/DT AND D(SIGMA)/DT AT POLES.
C     IN BELOW, SUBSCRIPTS 1 AND 2 REPRESENT J=1 AND PLAT RESPECTIVELY.

          W1(1)=0.0D0
C     CALCULATING DIVERGENCE AT POLES.
!$OMP PARALLEL DO PRIVATE (I,K,TMP_SUM)
          DO K=1,NL
            TMP_SUM=0.0D0
            DO I=ibeg1,iend1
              TMP_SUM=TMP_SUM+(P(I,J)+P(I,J1))*V(I,J2,K)
            ENDDO
            D1_16_TMP(1,K)=TMP_SUM
          ENDDO
          call t_stopf("DIFPS1 COMP")

          DO K=1,NL
            D1_16(K)=D1_16_TMP(1,K)
          ENDDO

          call gamil_sum_pole_data_phys(J,D1_16,NL)

          IF(J.EQ.1) THEN
             WK0=2.0D0/(DFLOAT(NX-2)*RAD*DY)*WTGV(J2)
          ELSE
             WK0=-2.0D0/(DFLOAT(NX-2)*RAD*DY)*WTGV(J2)
          ENDIF

          call t_startf("DIFPS1 COMP")
          DPSP=0.0D0
          DO K=1,NL
            D1_8(K) = D1_16(K)
            D1_8(K) = WK0*D1_8(K)
            DPSP=DPSP-DSIG(K)*D1_8(K)
          ENDDO
C
C     CALCULATING DPS/DT AND D(SIGMA)/DT AT POLES.
C
          IF (ibeg0 .eq. 1) pole_ps_1 = PS(1,J)
          call broadcast_lon_data(1,j,pole_ps_1,1)
          DO K=2,NL
            W1(K)=W1(K-1)-DSIG(K-1)/pole_ps_1(1)*(DPSP+D1_8(K-1))
          ENDDO
C
          DO I=ibeg1,iend1
            DPS(I,J)=DPSP
          ENDDO

!$OMP PARALLEL DO PRIVATE (I,K)
          DO K=2,NL
            DO I=ibeg1,iend1
              WS(I,J,K)=W1(K)
            ENDDO
          ENDDO

        call t_stopf("DIFPS1 COMP")
        ENDIF
      ENDDO

      call gamil_arrays_comm(COMM_TO_LEFT,1,DPS(:,beglatexdyn),WS(:,beglatexdyn,1),request_id=icomm_request1)
      call gamil_arrays_comm(COMM_TO_RIGHT,1,DPS(:,beglatexdyn),WS(:,beglatexdyn,1),request_id=icomm_request2)
      call wait_icomm_request(icomm_request1)
      call wait_icomm_request(icomm_request2)

      RETURN
      END
