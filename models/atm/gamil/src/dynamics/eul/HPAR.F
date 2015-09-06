#include <misc.h>
#include <params.h>


!-------------------------------------------------------
!
	SUBROUTINE HPAR(DX,DY,YTHU,YTHV,WTGU,WTGV,MDJ
     _               ,SINU,SINV,OUX,OUY,OVX,OVY,FF,CUR)
!
      use pmgrid, only: beglatexdyn,endlatexdyn, plat
      use mpi_gamil
      IMPLICIT NONE
!
#include <PARADYN>
!
	INTEGER J
!
!	INPUT VARIBALES:
!          1) HORIZONTAL RESOLUTION PARAMETERS (SPATIAL STEPSIZE): DX,DY
	REAL DX,DY
!
!		 2) VARIBLES RELATED TO ISO-AREA COORDINATE:
	REAL YTHU(beglatexdyn:endlatexdyn)
	REAL YTHV(beglatexdyn:endlatexdyn)
	REAL WTGU(beglatexdyn:endlatexdyn)
	REAL WTGV(beglatexdyn:endlatexdyn)
!
!          3) THE PARAMETER RELATED TO THE FLEXIBLE LEAPING-GRID METHOD
	INTEGER MDJ(beglatexdyn:endlatexdyn)
!
!     OUTPUT VARIABLES:
!
	REAL*8 OUX(beglatexdyn:endlatexdyn)
	REAL*8 OUY(beglatexdyn:endlatexdyn)
	REAL*8 OVX(beglatexdyn:endlatexdyn)
	REAL*8 OVY(beglatexdyn:endlatexdyn)
	REAL*8 SINU(beglatexdyn:endlatexdyn)
	REAL*8 SINV(beglatexdyn:endlatexdyn)
	REAL*8 FF(beglatexdyn:endlatexdyn)
	REAL*8 CUR(beglatexdyn:endlatexdyn)
!
!	WORKING VARIBALES:
!
	REAL YU,YV
!
      DO J=jbeg0,jend0
        IF(J.GE.2.AND.J.LE.PLAT-1) THEN
 	    YU=YTHU(J)
 	    YV=YTHV(J)
          FF(J)=2.0D0*OMGA*COS(YU)
          CUR(J)=COS(YU)/SIN(YU)/RAD
*         FF(J)=-2.0D0*OMGA*COS(YU)
*         CUR(J)=-COS(YU)/SIN(YU)/RAD
*         FF(NY-J+1)=2.0D0*OMGA*COS(YU)
*         CUR(NY-J+1)=COS(YU)/SIN(YU)/RAD
          SINU(J)=SIN(YU)
          SINV(J)=SIN(YV)
          OUX(J)=1.0D0/(RAD*SINU(J)*DX*MDJ(J))
          OUY(J)=1.0D0/(RAD*SINU(J)*DY)*WTGU(J)
          OVX(J)=1.0D0/(RAD*SINV(J)*DX*MDJ(J))
          OVY(J)=1.0D0/(RAD*SINV(J)*DY)*WTGV(J)
        ELSE IF(J.EQ.1) THEN
 		YV=YTHV(J)
 		SINV(1)=SIN(YV)
          OVX(1)=1.0D0/(RAD*SINV(1)*DX*MDJ(J))
          OVY(1)=1.0D0/(RAD*SINV(1)*DY)*WTGV(1)
          FF(J)=2.0D0*OMGA
          CUR(J)=0.0D0
*         FF(J)=-2.0D0*OMGA
*         CUR(J)=0.0D0
*         FF(NY)=2.0D0*OMGA
*         CUR(NY)=0.0D0
          SINU(J)=0.0D0
          OUX(J)=0.0D0
          OUY(J)=0.0D0
        ELSE
          FF(J)=-2.0D0*OMGA
          CUR(J)=0.0D0
*         FF(J)=2.0D0*OMGA
*         CUR(J)=0.0D0
*         FF(1)=-2.0D0*OMGA
*         CUR(1)=0.0D0
          SINU(J)=0.0D0
          SINV(J)=0.0D0
          OUX(J)=0.0D0
          OUY(J)=0.0D0
          OVX(J)=0.0D0
          OVY(J)=0.0D0
        ENDIF
      ENDDO

      call register_comm_array(1,1,beglatexdyn,endlatexdyn,1,1,1,1,SINV)
      call register_comm_array(1,1,beglatexdyn,endlatexdyn,1,1,1,1,SINU)
      call register_comm_array(1,1,beglatexdyn,endlatexdyn,1,1,1,1,FF)
      call register_comm_array(1,1,beglatexdyn,endlatexdyn,1,1,1,1,CUR)
      call register_comm_array(1,1,beglatexdyn,endlatexdyn,1,1,1,1,WTGV)
      call gamil_arrays_comm(COMM_TO_TOP,1,SINV,SINU,FF,CUR,WTGV)
      call gamil_arrays_comm(COMM_TO_BOT,1,SINV,SINU,FF,CUR,WTGV)
      call remove_comm_array(WTGV)
      call remove_comm_array(CUR)
      call remove_comm_array(FF)
      call remove_comm_array(SINU)
      call remove_comm_array(SINV)

	RETURN
	END

