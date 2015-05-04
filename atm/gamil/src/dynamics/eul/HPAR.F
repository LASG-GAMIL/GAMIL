#include <misc.h>
#include <params.h>


#if (defined SPMD)

!!(2003.10.27)
!!-------------------
!
	SUBROUTINE HPAR(DX,DY,YTHU,YTHV,WTGU,WTGV,MDJ
     _               ,SINU,SINV,OUX,OUY,OVX,OVY,FF,CUR)

      use mpishorthand, only: mpicom

      IMPLICIT NONE

#include <PARADYN>
#include <mpif.h>
#include <commpi.h>

!!      INTEGER J,J0
        integer j
        character*50 filename

!	INPUT VARIBALES:

	REAL*8 DX,DY
	REAL*8 YTHU(NY),YTHV(NY),WTGU(NY),WTGV(NY)
	INTEGER MDJ(NY)
!
!     OUTPUT VARIABLES:
!
	REAL*8 OUX(NY),OUY(NY),OVX(NY),OVY(NY)
     _      ,SINU(NY),SINV(NY),FF(NY),CUR(NY)
!
!	WORKING VARIBALES:
!
	REAL*8 YU,YV
        real*8 workst(4),workrb(4)
        real*8 worksb(4),workrt(4)

        integer isend1,isend2, irecv1,irecv2
        integer istatus(mpi_status_size)

!

      do j=jbeg,jend
         YU=YTHU(J)
         YV=YTHV(J)
         FF  (J) = 2.0D0*OMGA*COS(YU)
         CUR (J) = COS(YU)/SIN(YU)/RAD
         SINU(J) = SIN(YU)
         SINV(J) = SIN(YV)
         OUX (J) = 1.0D0/(RAD*SINU(J)*DX*MDJ(J))
         OUY (J) = 1.0D0/(RAD*SINU(J)*DY)*WTGU(J)
         OVX (J) = 1.0D0/(RAD*SINV(J)*DX*MDJ(J))
         OVY (J) = 1.0D0/(RAD*SINV(J)*DY)*WTGV(J)
      enddo

      if (myrank.eq.nprocs-1) then
          YV=YTHV(jpole)
          SINV(jpole) = SIN(YV)
          OVX (jpole) = 1.0D0/(RAD*SINV(jpole)*DX*MDJ(jpole))
          OVY (jpole) = 1.0D0/(RAD*SINV(jpole)*DY)*WTGV(jpole)
          FF  (jpole) = 2.0D0*OMGA
          CUR (jpole) = 0.0D0
          SINU(jpole) = 0.0D0
          OUX (jpole) = 0.0D0
          OUY (jpole) = 0.0D0
      endif

      if (myrank.eq.0) then
          FF  (jpole) = -2.0D0*OMGA
          CUR (jpole) = 0.0D0
          SINU(jpole) = 0.0D0
          SINV(jpole) = 0.0D0
          OUX (jpole) = 0.0D0
          OUY (jpole) = 0.0D0
          OVX (jpole) = 0.0D0
          OVY (jpole) = 0.0D0
      endif

      workst(1) = sinu(2)
      workst(2) = sinv(2)
      workst(3) = ff  (2)
      workst(4) = cur (2)
      worksb(1) = sinu(ny-1)
      worksb(2) = sinv(ny-1)
      worksb(3) = ff  (ny-1)
      worksb(4) = cur (ny-1)

      call mpi_isend( workst, 4, mpi_double_precision, itop, 1, mpicom, isend1, ierr)
      call mpi_isend( worksb, 4, mpi_double_precision, ibot, 2, mpicom, isend2, ierr)
      call mpi_irecv( workrb, 4, mpi_double_precision, ibot, 1, mpicom, irecv1, ierr)
      call mpi_irecv( workrt, 4, mpi_double_precision, itop, 2, mpicom, irecv2, ierr)

      call mpi_wait(isend1,istatus,ierr)
      call mpi_wait(isend2,istatus,ierr)
      call mpi_wait(irecv1,istatus,ierr)
      call mpi_wait(irecv2,istatus,ierr)

      sinu(1)  = workrt(1)
      sinv(1)  = workrt(2)
      ff  (1)  = workrt(3)
      cur (1)  = workrt(4)

      sinu(ny) = workrb(1)
      sinv(ny) = workrb(2)
      ff  (ny) = workrb(3)
      cur (ny) = workrb(4)


!-check-------------------------------------------------------
!      write(filename,12) 'hpar-p-',myrank,'.out'
!12    format(a7,i1,a4)
!
!      open (10,file=trim(filename))
!
!      write(10,*) '-------------------sinu-----------------'
!      write(10,11) (j,sinu(j),j=1,ny)
!
!      write(10,*) '-------------------sinv-----------------'
!      write(10,11) (j,sinv(j),j=1,ny)
!
!      write(10,*) '-------------------oux-----------------'
!      write(10,11) (j,oux(j),j=1,ny)
!
!      write(10,*) '-------------------ouy-----------------'
!!      write(10,11) (j,ouy(j),j=1,ny)
!
!      write(10,*) '-------------------ovx-----------------'
!      write(10,11) (j,ovx(j),j=1,ny)
!
!      write(10,*) '-------------------ovy-----------------'
!      write(10,11) (j,ovy(j),j=1,ny)
!
!      write(10,*) '-------------------ff-----------------'
!      write(10,11) (j,ff(j),j=1,ny)
!
!      write(10,*) '-------------------cur-----------------'
!      write(10,11) (j,cur(j),j=1,ny)
!
!11    format(1x,i5,f30.20)
!
!      close (10)
!
!!      call mpi_finalize(j)
!!      stop 'hpar'
!!-----------------------------------------------------------------


	RETURN
	END

#else
!
!-------------------------------------------------------
!
	SUBROUTINE HPAR(DX,DY,YTHU,YTHV,WTGU,WTGV,MDJ
     _               ,SINU,SINV,OUX,OUY,OVX,OVY,FF,CUR)
!
      IMPLICIT NONE
!
#include <PARADYN>
!
	INTEGER J,J0
!
!	INPUT VARIBALES:
!          1) HORIZONTAL RESOLUTION PARAMETERS (SPATIAL STEPSIZE): DX,DY
	REAL DX,DY
!
!		 2) VARIBLES RELATED TO ISO-AREA COORDINATE:
	REAL YTHU(NY),YTHV(NY),WTGU(NY),WTGV(NY)
!
!          3) THE PARAMETER RELATED TO THE FLEXIBLE LEAPING-GRID METHOD
	INTEGER MDJ(NY)
!
!     OUTPUT VARIABLES:
!
	REAL OUX(NY),OUY(NY),OVX(NY),OVY(NY)
     _      ,SINU(NY),SINV(NY),FF(NY),CUR(NY)
!
!	WORKING VARIBALES:
!
	REAL YU,YV
!
      J0=NY-1
      DO J=1,NY
        IF(J.GE.2.AND.J.LE.J0) THEN
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

*     DO J=1,NY
*	print *,sinv(j),WTGU(j),WTGV(J)
*	print *,OVX(j),OVY(j)
*     ENDDO
*	pause

!- check ----------------------------------------------------
!      open (10,file='hpar-s.out')
!
!      write(10,*) '-------------------sinu-----------------'
!      write(10,11) (j,sinu(j),j=1,ny)
!
!      write(10,*) '-------------------sinv-----------------'
!      write(10,11) (j,sinv(j),j=1,ny)
!
!      write(10,*) '-------------------oux-----------------'
!      write(10,11) (j,oux(j),j=1,ny)
!
!      write(10,*) '-------------------ouy-----------------'
!      write(10,11) (j,ouy(j),j=1,ny)
!
!      write(10,*) '-------------------ovx-----------------'
!      write(10,11) (j,ovx(j),j=1,ny)
!
!      write(10,*) '-------------------ovy-----------------'
!      write(10,11) (j,ovy(j),j=1,ny)
!
!      write(10,*) '-------------------ff-----------------'
!      write(10,11) (j,ff(j),j=1,ny)
!
!      write(10,*) '-------------------cur-----------------'
!      write(10,11) (j,cur(j),j=1,ny)
!
!11    format(1x,i10,f30.20)
!
!      close (10)
!!      stop
!-------------------------------------------------------------(wh)--------


	RETURN
	END

#endif
