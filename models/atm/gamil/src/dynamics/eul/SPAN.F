#include <misc.h>
#include <params.h>

!
!----------------------------------------------------------------
!
      SUBROUTINE SPAN(MM1,MM2,MM3,MP1,MP2,MP3,MDJ)
      use pmgrid, only: beglatexdyn, endlatexdyn
      use mpi_gamil
!
      IMPLICIT NONE
!
#include <PARADYN>
!
      INTEGER MM1(beglonex:endlonex,beglatexdyn:endlatexdyn)
      INTEGER MM2(beglonex:endlonex,beglatexdyn:endlatexdyn)
      INTEGER MM3(beglonex:endlonex,beglatexdyn:endlatexdyn)
      INTEGER MP1(beglonex:endlonex,beglatexdyn:endlatexdyn)
      INTEGER MP2(beglonex:endlonex,beglatexdyn:endlatexdyn)
      INTEGER MP3(beglonex:endlonex,beglatexdyn:endlatexdyn)
      INTEGER MDJ(beglatexdyn:endlatexdyn)
!
      INTEGER I,J
!

!************************************************************************
!               Set the values to MM1,MP1;MM2,MP2;MM3,MP3               *
!                    MM1=i-mdj,if MM1<2, MM1=NX-2+MM1                   *
!                  MP1=i+mdj,if MP1>NX-1, MP1=MP1-NX+2                  *
!                 MM2=i-(mdj-1)/2,if MM2<2, MM2=NX-2+MM2                *
!               MP2=i+(mdj+1)/2,if MP2>NX-1, MP2=MP2-NX+2               *
!                 MM3=i-(mdj+1)/2,if MM3<2, MM3=NX-2+MM3                *
!               MP3=i+(mdj-1)/2,if MP3>NX-1, MP3=MP3-NX+2               *
!************************************************************************
      DO J=jbeg0,jend0
         MDJ(J)=1
      END DO
!
      DO J=jbeg0,jend0
      DO I=beglonex,endlonex
	    MM1(I,J)=I-MDJ(J)
	    MP1(I,J)=I+MDJ(J)
	    MM2(I,J)=I-(MDJ(J)-1)/2
	    MP2(I,J)=I+(MDJ(J)+1)/2
	    MM3(I,J)=I-(MDJ(J)+1)/2
	    MP3(I,J)=I+(MDJ(J)-1)/2
      END DO
      END DO
!******************           Made by Wang Bin           *****************
!*************************************************************************

!- check -----------------------------------------------------------------
!      open (10,file='span-s.out')
!      write(10,11) mdj
!      write(10,*) '-------------------mp2-----------------'
!      do j=1,ny
!        write(10,*) 'j=',j,'----------------------'
!        write(10,11) (mp2(i,j),i=1,nx)
!11    format(1x,10i9)
!      enddo
!!      close (10)
!!      stop
!---------------------------------------------------------------(wh)------

	RETURN
	END

