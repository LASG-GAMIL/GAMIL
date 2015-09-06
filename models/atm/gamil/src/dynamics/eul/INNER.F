#include <misc.h>
#include <params.h>

!!(wh 2003.11.10)
!!-------------------


	FUNCTION INNER(U1,V1,T1,P1,L1,U2,V2,T2,P2,L2
     _                ,DSIG,TBS,PSB,PMTOP,SINU,SINV,WTGU,WTGV)

      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
      
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

!
!	The file PARA is to define the parameters related to the model resolution:
!     NX is the grid number in longitude
!     NY is the grid number in latitude
!     NL is the number of vertical layers
!
	REAL*8
     _       U1(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !	input variable
     _      ,V1(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !	input variable
     _      ,T1(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !	input variable
     _      ,P1(beglonex:endlonex,beglatexdyn:endlatexdyn)     !	input variable
     _      ,U2(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !	input variable
     _      ,V2(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !	input variable
     _      ,T2(beglonex:endlonex,beglatexdyn:endlatexdyn,NL)  !	input variable
     _      ,P2(beglonex:endlonex,beglatexdyn:endlatexdyn)     !	input variable
     _      ,TBS(beglonex:endlonex,beglatexdyn:endlatexdyn)    !	input variable
     _      ,PSB(beglonex:endlonex,beglatexdyn:endlatexdyn)    !	input variable
     _      ,SINU(beglatexdyn:endlatexdyn)      !	input variable
     _      ,SINV(beglatexdyn:endlatexdyn)      !	input variable
     _      ,WTGU(beglatexdyn:endlatexdyn)		 !	input variable
     _      ,WTGV(beglatexdyn:endlatexdyn)      !	input variable
     _      ,DSIG(NL)		 !	input variable
     _      ,PMTOP		 !	input variable
     _      ,INNER         !  output variable
!
	INTEGER
     _        L1			 !	input variable
     _       ,L2			 !	input variable
     _       ,I,J,K     !  working variables
!
      REAL*8  TMPU(beglonex:endlonex,jbeg0:jend0)	 !  working variables
     _       ,TMPV(beglonex:endlonex,jbeg0:jend0)	 !  working variables
     _       ,TMPT(beglonex:endlonex,jbeg0:jend0)	 !  working variables
     _       ,TMP_DATA

      REAL*16 EUJ(8,jbeg0:jend0)      !  working variables
     _       ,EVJ(8,jbeg0:jend0)      !  working variables
     _       ,ETJ(8,jbeg0:jend0)      !  working variables
     _       ,EPJ(8,jbeg0:jend0)      !  working variables
     _       ,DS           !  working variables
     _       ,DSU          !  working variables
     _       ,DSV          !  working variables
     _       ,DPS1         !  working variables
     _       ,DPS2         !  working variables

        REAL*16 TMP_INNER

!
!
!$OMP PARALLEL DO PRIVATE (I,J,K,DPS1,DPS2,TMP_DATA)
       DO J=jbeg0,jend0
        DO I = ibeg1,iend1
          TMPU(I,J) = 0.0D0
          TMPV(I,J) = 0.0D0
          TMPT(I,J) = 0.0D0
        ENDDO
        DO K = 1,NL
          DO I = ibeg1,iend1
            TMPU(I,J) = TMPU(I,J)+U1(I,J,K)*U2(I,J,K)*DSIG(K)
            TMPV(I,J) = TMPV(I,J)+V1(I,J,K)*V2(I,J,K)*DSIG(K)
            TMPT(I,J) = TMPT(I,J)+T1(I,J,K)*T2(I,J,K)*DSIG(K)
          ENDDO
        ENDDO
        EUJ(1,J) = 0.0D0
        EVJ(1,J) = 0.0D0
        ETJ(1,J) = 0.0D0
        EPJ(1,J) = 0.0D0
        DO I = ibeg1,iend1
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
          TMP_DATA = RD*TBS(I,J)/PSB(I,J)*DPS1*DPS2
          EPJ(1,J) = EPJ(1,J) + TMP_DATA
          EUJ(1,J) = EUJ(1,J) + TMPU(I,J)
          EVJ(1,J) = EVJ(1,J) + TMPV(I,J)
          ETJ(1,J) = ETJ(1,J) + TMPT(I,J)
        ENDDO
      ENDDO

      
!
      TMP_INNER=0.0D0
!

      DO J=jbeg0,jend0
        IF(J.EQ.1) THEN
          DS=0.25D0*SINV(1)/WTGV(1)
          DSU=DS
          DSV=4.0D0*DSU
        ELSE IF(J.EQ.PLAT) THEN
          DS=0.25D0*SINV(J-1)/WTGV(J-1)
          DSU=DS
          DSV=0.0D0
        ELSE
          DS=SINU(J)/WTGU(J)
          DSU=DS
          DSV=SINV(J)/WTGV(J)
        ENDIF
        TMP_INNER=TMP_INNER+DSU*EUJ(1,J)+DSV*EVJ(1,J)+DSU*ETJ(1,J)+DS*EPJ(1,J)
      ENDDO

      call gamil_all_reduce(TMP_INNER)
      INNER = TMP_INNER


!
	RETURN
	END
