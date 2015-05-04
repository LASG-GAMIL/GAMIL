#include <misc.h>
#include <params.h>

#if (defined SPMD)

!!(wh 2003.10.27)
!!-------------------

      SUBROUTINE SPAN(MM1,MM2,MM3,MP1,MP2,MP3,MDJ)

      IMPLICIT NONE

#include <PARADYN>

#include <commpi.h>

      INTEGER MM1(NX,NY),MM2(NX,NY),MM3(NX,NY)
      integer MP1(NX,NY),MP2(NX,NY),MP3(NX,NY),MDJ(NY)

      INTEGER I,J
      integer index   ! global index of each latitude
      integer max,min

      character*50 filename


!************************************************************************
!               Set the values to MM1,MP1;MM2,MP2;MM3,MP3               *
!                    MM1=i-mdj,if MM1<2, MM1=NX-2+MM1                   *
!                  MP1=i+mdj,if MP1>NX-1, MP1=MP1-NX+2                  *
!                 MM2=i-(mdj-1)/2,if MM2<2, MM2=NX-2+MM2                *
!               MP2=i+(mdj+1)/2,if MP2>NX-1, MP2=MP2-NX+2               *
!                 MM3=i-(mdj+1)/2,if MM3<2, MM3=NX-2+MM3                *
!               MP3=i+(mdj-1)/2,if MP3>NX-1, MP3=MP3-NX+2               *
!************************************************************************
!!      MDJ(1)=3
!!      MDJ(2)=3
!!      MDJ(3)=3
!!      MDJ(4)=3
!!      MDJ(NY)=3
!!      MDJ(NY-1)=3
!!      MDJ(NY-2)=3
!!      MDJ(NY-3)=3
!!      DO J=5,NY-4
!!         MDJ(J)=1
!!      END DO

       min = 4
       max = nprocs*(ny-2)-3

       do j=2,ny-1
          index = (nprocs-1-myrank)*(ny-2)+(j-1)
          if ( (index.le.min).or.(index.ge.max) ) then
!              mdj(j) = 3
               mdj(j) = 1
          else
	       mdj(j) = 1
          endif
       enddo

!!    DO J=1,NY
      do j=2,ny-1
      DO I=1,NX
	 MM1(I,J)=I-MDJ(J)
	 MP1(I,J)=I+MDJ(J)
	 MM2(I,J)=I-(MDJ(J)-1)/2
	 MP2(I,J)=I+(MDJ(J)+1)/2
	 MM3(I,J)=I-(MDJ(J)+1)/2
	 MP3(I,J)=I+(MDJ(J)-1)/2
	 IF (MM1(I,J).LT.2) MM1(I,J)=MM1(I,J)+NX-2
	 IF (MP1(I,J).GE.NX) MP1(I,J)=MP1(I,J)-NX+2
	 IF (MM2(I,J).LT.2) MM2(I,J)=MM2(I,J)+NX-2
	 IF (MP2(I,J).GE.NX) MP2(I,J)=MP2(I,J)-NX+2
	 IF (MM3(I,J).LT.2) MM3(I,J)=MM3(I,J)+NX-2
	 IF (MP3(I,J).GE.NX) MP3(I,J)=MP3(I,J)-NX+2
      END DO
      END DO
!******************           Made by Wang Bin           *****************
!*************************************************************************

!- check ----------------------------------------------------------------
!      write(filename,12) 'span-p-',myrank,'.out'
!12    format(a7,i1,a4)
!
!      open (10,file=trim(filename))
!      write(10,11) mdj
!      write(10,*) '-------------------mp2-----------------'
!      do j=1,ny
!        write(10,*) 'j=',j,'----------------------'
!        write(10,11) (mp2(i,j),i=1,nx)
!11    format(1x,10i9)
!      enddo
!      close (10)
!
!      call mpi_finalize(j)
!      stop 'span'
!----------------------------------------------------------------(wh)----


	RETURN
	END


#else
!
!----------------------------------------------------------------
!
	SUBROUTINE SPAN(MM1,MM2,MM3,MP1,MP2,MP3,MDJ)
!
      IMPLICIT NONE
!
#include <PARADYN>
!
      INTEGER MM1(NX,NY),MM2(NX,NY),MM3(NX,NY)
     &        ,MP1(NX,NY),MP2(NX,NY),MP3(NX,NY),MDJ(NY)
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
      MDJ(1)=1
      MDJ(2)=1
      MDJ(3)=1
      MDJ(4)=1
      MDJ(NY)=1
      MDJ(NY-1)=1
      MDJ(NY-2)=1
      MDJ(NY-3)=1
      DO J=5,NY-4
         MDJ(J)=1
      END DO
!
      DO J=1,NY
      DO I=1,NX
	 MM1(I,J)=I-MDJ(J)
	 MP1(I,J)=I+MDJ(J)
	 MM2(I,J)=I-(MDJ(J)-1)/2
	 MP2(I,J)=I+(MDJ(J)+1)/2
	 MM3(I,J)=I-(MDJ(J)+1)/2
	 MP3(I,J)=I+(MDJ(J)-1)/2
	 IF (MM1(I,J).LT.2) MM1(I,J)=MM1(I,J)+NX-2
	 IF (MP1(I,J).GE.NX) MP1(I,J)=MP1(I,J)-NX+2
	 IF (MM2(I,J).LT.2) MM2(I,J)=MM2(I,J)+NX-2
	 IF (MP2(I,J).GE.NX) MP2(I,J)=MP2(I,J)-NX+2
	 IF (MM3(I,J).LT.2) MM3(I,J)=MM3(I,J)+NX-2
	 IF (MP3(I,J).GE.NX) MP3(I,J)=MP3(I,J)-NX+2
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

#endif
