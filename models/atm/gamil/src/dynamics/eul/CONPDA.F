#include <misc.h>
#include <params.h>


#if (defined SPMD)

!!(2003.10.28)
!!(2003.12.12)
!!-------------------

      SUBROUTINE CONPDA(DTIMEQ,DX,DY,SINU,WTGV,DSIG,FIRST,ISOR
*    _                 ,IORD,IP,DSNP,DSSP,GC,DTDLT,DTDLN,DTDSG,DXYP)
     _                 ,IORD,IP,DSNP,DSSP,GC,DTDLT,DTDLN,DTDSG)

        use mpishorthand, only: mpicom

	IMPLICIT NONE

#include <PARADYN>
#include <mpif.h>
#include <commpi.h>

      REAL*8 ZERO,HALF,FOUR
      DATA ZERO,HALF,FOUR/0.0D0,0.5D0,4.0D0/

	  REAL*8 WTGV(NY),DX,DY,DTIMEQ,FIRST
     _      ,DXYP(NY),SINU(NY)
     _      ,DTDLT(NY),DTDLN(NY),DTDSG(NL)
     _      ,FIM,DSNP,DSSP,DSIG(NL),GC(NY)
	INTEGER I,ISOR,IORD,J,JB,JE,K,IP(NX)
        character*50 filename

!!	JB=2
!!	JE=NY-1

      do j=jbeg,jend
         DXYP(J) = HALF*(HALF*(RAD*DX)*(SINU(J-1)+SINU(J+1)+SINU(J)*2.0))*(RAD*DY)
      END DO

      if (myrank.eq.nprocs-1)
     _   DXYP(jpole) = (HALF*(RAD*DX)*(SINU(jpole)+SINU(jpole+1)))*(RAD*DY)/FOUR

      if (myrank.eq.0)
     _   DXYP(jpole) = (HALF*(RAD*DX)*(SINU(jpole)+SINU(jpole-1)))*(RAD*DY)/FOUR
!
!     FOR THE USAGE OF SUB.MPDATA & VPDATA
!
      FIM = DBLE( IM )
!
      IF( FIRST.LE.ZERO ) THEN
        IF( ISOR.EQ.3 ) IORD = MAX0(IORD,3)
        DO I= 1 ,NX
           IP(I)   = MOD( I+IM/2-1,IM ) + 1
	  END DO

	if (nprocs.lt.32) then
	    if (myrank.eq.nprocs-1) dsnp = dxyp(jpole+1)/(fim*dxyp(jpole))
            if (myrank.eq.0       ) dssp = dxyp(jpole-1)/(fim*dxyp(jpole))
	else
!!        -------------------------------for the north pole
            if (myrank.eq.nprocs-2) then
             call mpi_send( dxyp(2), 1, mpi_double_precision, nprocs-1,
     _                               1, mpicom,   ierr)
	    endif

            if (myrank.eq.nprocs-1) then
             call mpi_recv( dxyp(3), 1, mpi_double_precision, nprocs-2,
     _                               1, mpicom,   ierr)
	       dsnp=dxyp(3)/(fim*dxyp(2))
            endif
!!        --------------------------------for the south pole
            if (myrank.eq.1) then
             call mpi_send( dxyp(ny-1), 1, mpi_double_precision, 0,
     _                                  1, mpicom, ierr)
            endif

            if (myrank.eq.0) then
             call mpi_recv( dxyp(1),    1, mpi_double_precision, 1,
     _                                  1, mpicom, ierr)
	       dsnp=dxyp(1)/(fim*dxyp(2))
            endif
!!        ---------------------------------------------------
        endif

        DO J= 2,ny-1
           GC(J)   = SINU(J)
	  END DO
      ENDIF
C
      DO J  = 2 ,jend
         DTDLT(J)  = DTIMEQ / (RAD*DY)*WTGV(J)
	END DO
      DO J  = jbeg,jend
         DTDLN(J)  = DTIMEQ / (RAD*SINU(J)*DX)
      END DO
      DO K  = 1 ,NL
         DTDSG(K)  = DTIMEQ / DSIG(K)
      END DO

!- check ---------------------------------------------------------
!
!      write(filename,14) 'conpda-p-',myrank,'.out'
!14    format(a9,i1,a4)
!
!      open (10,file=trim(filename))
!
!      write(10,*) 'iord=',iord,'dsnp=',dsnp,'dssp=',dssp
!
!      write(10,*) '------------------- ip -----------------'
!      write(10,12) (ip(i),i=1,nx)
!12    format(1x,10i9)
!
!      write(10,*) '------------------- dtdsg -----------------'
!      write(10,13) (dtdsg(k),k=1,nl)
!13    format(1x,3e30.20)
!
!      write(10,*) '------------------------------------------------- gc ---------'
!      write(10,11) (j,gc(j),j=1,ny)
!
!      write(10,*) '------------------------------------------------dtdlt----------'
!      write(10,11) (j,dtdlt(j),j=1,ny)
!
!      write(10,*) '------------------------------------------------dtdln----------'
!      write(10,11) (j,dtdln(j),j=1,ny)
!!
!11    format(1x,i5,e30.20)
!
!      close (10)
!
!!      call mpi_finalize(j)
!!      stop 'conpda'
!
!!-----------------------------------------------------------------


      RETURN
      END
#else


C     ***************************
C     ***************************
      SUBROUTINE CONPDA(DTIMEQ,DX,DY,SINU,WTGV,DSIG,FIRST,ISOR
*    _                 ,IORD,IP,DSNP,DSSP,GC,DTDLT,DTDLN,DTDSG,DXYP)
     _                 ,IORD,IP,DSNP,DSSP,GC,DTDLT,DTDLN,DTDSG)

C     ***************************
C     ***************************
!
	IMPLICIT NONE
!
#include <PARADYN>
!
      REAL*8 ZERO,HALF,FOUR
      DATA   ZERO,HALF,FOUR/0.0D0,0.5D0,4.0D0/
!
	  REAL*8 WTGV(NY),DX,DY,DTIMEQ,FIRST
     _      ,DXYP(NY),SINU(NY)
     _      ,DTDLT(NY),DTDLN(NY),DTDSG(NL)
     _      ,FIM,DSNP,DSSP,DSIG(NL),GC(NY)
	INTEGER I,ISOR,IORD,J,JB,JE,K,IP(NX)
!
	JB=2
	JE=NY-1
!
      DO J  = JB,JE
      DXYP(J)   = HALF*(HALF*(RAD*DX)*(SINU(J-1)
     _           +SINU(J+1)+SINU(J)*2.0))*(RAD*DY)
	END DO
!
      DXYP(1 )  = (HALF*(RAD*DX)*(SINU(2)+SINU(1)))*(RAD*DY) / FOUR
      DXYP(NY)  = (HALF*(RAD*DX)*(SINU(JE)+SINU(NY)))*(RAD*DY)/FOUR
!
!     FOR THE USAGE OF SUB.MPDATA & VPDATA
!
      FIM       = DBLE( IM )
!
      IF( FIRST.LE.ZERO ) THEN
        IF( ISOR.EQ.3 ) IORD = MAX0(IORD,3)
        DO I= 1 ,NX
           IP(I)   = MOD( I+IM/2-1,IM ) + 1
	  END DO
        DSNP    = DXYP(JB) / (FIM*DXYP(01))
        DSSP    = DXYP(JE) / (FIM*DXYP(NY))
        DO J= 1 ,NY
           GC(J)   = SINU(J)
	  END DO
      ENDIF
C
      DO J  = 1 ,JE
         DTDLT(J)  = DTIMEQ / (RAD*DY)*WTGV(J)
	END DO
      DO J  = JB,JE
         DTDLN(J)  = DTIMEQ / (RAD*SINU(J)*DX)
	END DO
      DO K  = 1 ,NL
         DTDSG(K)  = DTIMEQ / DSIG(K)
      END DO


!- check ---------------------------------------------------------
!
!      open (10,file='conpda-s.out')
!
!      write(10,*) 'iord=',iord,'dsnp=',dsnp,'dssp=',dssp
!
!      write(10,*) '------------------- ip -----------------'
!      write(10,12) (ip(i),i=1,nx)
!12    format(1x,10i9)
!
!      write(10,*) '------------------- dtdsg -----------------'
!      write(10,13) (dtdsg(k),k=1,nl)
!13    format(1x,3e30.20)
!
!      write(10,*) '------------------------------------------------- gc ---------'
!      write(10,11) (j,gc(j),j=1,ny)
!
!      write(10,*) '------------------------------------------------dtdlt----------'
!      write(10,11) (j,dtdlt(j),j=1,ny)
!
!      write(10,*) '------------------------------------------------dtdln----------'
!      write(10,11) (j,dtdln(j),j=1,ny)
!
!11    format(1x,i10,e30.20)
!
!      close (10)
!!      stop
!!
!------------------------------------------------------------------------(wh)----

      RETURN
      END
#endif
