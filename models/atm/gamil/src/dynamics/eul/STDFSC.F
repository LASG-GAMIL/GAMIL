#include <misc.h>
#include <params.h>


#if (defined SPMD)

!!(wh 2003.10.28)
!!-------------------


      SUBROUTINE STDFSC(DFS0,DTHDFS,SINU,SINV,WTGU,WTGV,DX,DY
*    _                 ,FRDT,FRDS,FRDU,FRDV,FRDP,DXVPN,DXVPS,DXYV)
     _                 ,FRDT,FRDS,FRDU,FRDV,FRDP,DXVPN,DXVPS)

        use mpishorthand, only: mpicom

	IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

#include <mpif.h>
#include <commpi.h>
!
	REAL*8 DFS0,DTHDFS,SINU(NY),SINV(NY),WTGU(NY),WTGV(NY),DX,DY
     _      ,FRDT(NY,3),FRDS(NY,3),FRDU(NY,3),FRDV(NY,3),FRDP(NY,3)
     _      ,DXR,DYR,ALN,ALT,ONE,ANT,AVT,AUT,AUO,AVO,ANTS,DNSI,FRTP
     _      ,SNU,SNV,VST,SUI,SUJ,SVI,SVJ,FPO,CDFS,ZERO,HALF,FOUR,FIM
     _      ,DXYV(NY),DXVPN(NY),DXVPS(NY),DAN,DAS,DAP
	DATA ZERO,HALF,ONE,FOUR/0.0D0,0.5D0,1.0D0,4.0D0/
	INTEGER J,K
        character*50 filename

	integer isend,irecv
	integer istatus(mpi_status_size)
C
C     SET CONSTANTS     FOR COMPUTATION OF THE HORIZONTAL DIFFUSION
C
!!    DO J  = 1 ,JE
      do j=2,jend
         DXYV(J)   = HALF*(RAD*DX*(SINU(J+1)+SINU(J)))*(RAD*DY)
!!       DXYV(J)   = (RAD*DX)*(RAD*DY)
	END DO

!!    DXYV(NY)=0.D0
      if (myrank.eq.0) dxyv(jpole) = 0.0d0


      call mpi_isend( dxyv(ny-1), 1, mpi_double_precision, ibot, 1, mpicom, isend, ierr)
      call mpi_irecv( dxyv(1)   , 1, mpi_double_precision, itop, 1, mpicom, irecv, ierr)

      call mpi_wait(isend,istatus,ierr)
      call mpi_wait(irecv,istatus,ierr)

      DO 420 J  = jbeg,jend
      DAN       = DXYV(J-1)
      DAS       = DXYV(J  )
      DAP       = DAN + DAS
      DXVPN(J)  = DAN / DAP
      DXVPS(J)  = DAS / DAP
420   CONTINUE

      if (inc_pole) then
         DXVPN(jpole)  = ZERO
         DXVPS(jpole)  = ZERO
         DO 25  K  = 1 ,3
           FRDT(jpole,K) = ZERO
           FRDS(jpole,K) = ZERO
           FRDU(jpole,K) = ZERO
           FRDV(jpole,K) = ZERO
           FRDP(jpole,K) = ZERO
25       CONTINUE
      endif

 	FIM=DBLE(IM)
!
      DXR     = DX  * RAD
      DYR     = DY  * RAD
      ALN     = ONE   / DXR
      ALT     = ONE   / DYR
      ANT     = DX  / DY
      AVT     = DXR * ANT
      AUT     = AVT   / FOUR
      AVO     = DXR / FOUR
      AUO     = DXR
      ANTS    = ANT   * ANT
      DNSI    = DX*DX / FIM
      FRTP    = FOUR*ALT  / FIM
C
!!    DO 100 J  = JB,JE
      do 100 j=jbeg,jend
         SNU       = SINU(J)
         FRDT(J,1) = ALN   / SNU
         VST       = ALT   / SNU * WTGU(J)
         FRDT(J,2) = VST   * SINV(J)
         FRDT(J,3) = VST   * SINV(J-1)
100   CONTINUE
!!    FRDT(1 ,1)=-FRTP
      if (myrank.eq.nprocs-1) frdt(jpole,1)=-frtp
!!    FRDT(NY,1)= FRTP
      if (myrank.eq.0       ) frdt(jpole,1)= frtp

!!    DO 125 J  = 1 ,JE
      do 125 j  = 2, jend
         SNV       = SINV(J)
         VST       = ALT   / SNV * WTGV(J)
         FRDS(J,1) = ALN   / SNV
         FRDS(J,2) = VST   * SINU(J+1)
         FRDS(J,3) = VST   * SINU(J)
125   CONTINUE
C
!!    DO 150 J  = JB,JE
      do 150 j  = jbeg,jend
         SVI       = SINV(J)
         SVJ       = SINV(J-1)
         FRDU(J,1) = AUO   * SINU(J)
         FRDU(J,2) = AUT   * SVI*SVI*WTGU(J)
         FRDU(J,3) = AUT   * SVJ*SVJ*WTGU(J)
150   CONTINUE
C
!!    DO 175 J  = 1 ,JE
      do 175 j  = 2, jend
         SUI       = SINU(J)
         SUJ       = SINU(J+1)
         FRDV(J,1) = AVO   * SINV(J)
         FRDV(J,2) = AVT   * SUJ*SUJ*WTGV(J)
         FRDV(J,3) = AVT   * SUI*SUI*WTGV(J)
175   CONTINUE
C
!!    DO 200 J  = JB,JE
      do 200 j  = jbeg,jend
         SVI       = SINV(J)
         SVJ       = SINV(J-1)
         FPO       = ANTS  / SINU(J)
         FRDP(J,1) = ONE
         FRDP(J,2) = FPO   * SVI*SVI*SVI*WTGU(J)*WTGV(J)
         FRDP(J,3) = FPO   * SVJ*SVJ*SVJ*WTGU(J)*WTGV(J-1)
200   CONTINUE
!!    FRDP(1 ,1)= + DNSI
      if (myrank.eq.nprocs-1) frdp(jpole,1)= dnsi
!!    FRDP(NY,1)= - DNSI
      if (myrank.eq.0       ) frdp(jpole,1)=-dnsi

      CDFS      = DTHDFS * DFS0*DFS0
      DO 250 K  = 1 ,3
!!    DO 250 J  = 1 ,NY
      do 250 j  = 2 ,ny-1
         FRDU(J,K) = FRDU(J,K) * CDFS
         FRDV(J,K) = FRDV(J,K) * CDFS
         FRDP(J,K) = FRDP(J,K) * CDFS
250   CONTINUE

      DO 260 K  = 1 ,3
!!    DO 260 J  = 1 ,NY
      do 260 j  = 2 ,ny-1
         FRDP(J,K) = 3.0E0 * FRDP(J,K)
260   CONTINUE


!- check ------------------------------------------------------
!
!      write(filename,12) 'stdfsc-p-',myrank,'.out'
!12    format(a9,i1,a4)
!
!      open (10,file=trim(filename))
!
!      write(10,*) '----------------- dxvpn -----------------'
!      do j=1,ny
!        write(10,11) j,dxvpn(j)
!      enddo
!
!      write(10,*) '----------------- dxvps -----------------'
!      do j=1,ny
!        write(10,11) j,dxvps(j)
!      enddo
!!
!      write(10,*) '----------------- frdt -----------------'
!      do j=1,ny
!        write(10,11) (j,(frdt(j,k),k=1,3))
!      enddo
!!
!      write(10,*) '----------------- frds -----------------'
!      do j=1,ny
!        write(10,11) (j,(frds(j,k),k=1,3))
!      enddo
!
!      write(10,*) '----------------- frdu -----------------'
!      do j=1,ny
!        write(10,11) (j,(frdu(j,k),k=1,3))
!      enddo
!
!      write(10,*) '----------------- frdv -----------------'
!      do j=1,ny
!        write(10,11) (j,(frdv(j,k),k=1,3))
!      enddo
!
!      write(10,*) '----------------- frdp -----------------'
!      do j=1,ny
!        write(10,11) (j,(frdp(j,k),k=1,3))
!      enddo
!
!11    format(1x,i5,3e30.20)
!
!      close (10)
!
!!      call mpi_finalize(j)
!!      stop'stdfsc'
!!--------------------------------------------------------------

      RETURN
      END
#else
C     +++++++++++++++++++++++++++
C     :::::::::::::::::::::::::::
      SUBROUTINE STDFSC(DFS0,DTHDFS,SINU,SINV,WTGU,WTGV,DX,DY
*    _                 ,FRDT,FRDS,FRDU,FRDV,FRDP,DXVPN,DXVPS,DXYV)
     _                 ,FRDT,FRDS,FRDU,FRDV,FRDP,DXVPN,DXVPS)
C     :::::::::::::::::::::::::::
C     +++++++++++++++++++++++++++
	IMPLICIT NONE

#include <PARADYN>
#include <PARADD>

!
	REAL*8 DFS0,DTHDFS,SINU(NY),SINV(NY),WTGU(NY),WTGV(NY),DX,DY
     _      ,FRDT(NY,3),FRDS(NY,3),FRDU(NY,3),FRDV(NY,3),FRDP(NY,3)
     _      ,DXR,DYR,ALN,ALT,ONE,ANT,AVT,AUT,AUO,AVO,ANTS,DNSI,FRTP
     _      ,SNU,SNV,VST,SUI,SUJ,SVI,SVJ,FPO,CDFS,ZERO,HALF,FOUR,FIM
     _      ,DXYV(NY),DXVPN(NY),DXVPS(NY),DAN,DAS,DAP
	DATA ZERO,HALF,ONE,FOUR/0.0D0,0.5D0,1.0D0,4.0D0/
	INTEGER J,K
C
C     SET CONSTANTS     FOR COMPUTATION OF THE HORIZONTAL DIFFUSION
C
      DO J  = 1 ,JE
         DXYV(J)   = HALF*(RAD*DX*(SINU(J+1)+SINU(J)))*(RAD*DY)
!         DXYV(J)   = (RAD*DX)*(RAD*DY)
	END DO
      DXYV(NY)=0.D0
      DO 420 J  = JB,JE
      DAN       = DXYV(J-1)
      DAS       = DXYV(J  )
      DAP       = DAN + DAS
      DXVPN(J)  = DAN / DAP
      DXVPS(J)  = DAS / DAP
420   CONTINUE
      DXVPN(1)  = ZERO
      DXVPS(1)  = ZERO
      DXVPN(NY) = ZERO
      DXVPS(NY) = ZERO
!	JB=2
!	JE=NY-1
	FIM=DBLE(IM)
      DO 25  J  = 1 ,NY,JE
      DO 25  K  = 1 ,3
      FRDT(J,K) = ZERO
      FRDS(J,K) = ZERO
      FRDU(J,K) = ZERO
      FRDV(J,K) = ZERO
      FRDP(J,K) = ZERO
25    CONTINUE
!
      DXR     = DX  * RAD
      DYR     = DY  * RAD
      ALN       = ONE   / DXR
      ALT       = ONE   / DYR
      ANT       = DX  / DY
      AVT       = DXR * ANT
      AUT       = AVT   / FOUR
      AVO       = DXR / FOUR
      AUO       = DXR
      ANTS      = ANT   * ANT
      DNSI      = DX*DX / FIM
      FRTP      = FOUR*ALT  / FIM
C
      DO 100 J  = JB,JE
      SNU       = SINU(J)
      FRDT(J,1) = ALN   / SNU
      VST       = ALT   / SNU * WTGU(J)
      FRDT(J,2) = VST   * SINV(J)
      FRDT(J,3) = VST   * SINV(J-1)
100   CONTINUE
      FRDT(1 ,1)=-FRTP
      FRDT(NY,1)= FRTP
C
      DO 125 J  = 1 ,JE
      SNV       = SINV(J)
      VST       = ALT   / SNV * WTGV(J)
      FRDS(J,1) = ALN   / SNV
      FRDS(J,2) = VST   * SINU(J+1)
      FRDS(J,3) = VST   * SINU(J)
125   CONTINUE
C
      DO 150 J  = JB,JE
      SVI       = SINV(J)
      SVJ       = SINV(J-1)
      FRDU(J,1) = AUO   * SINU(J)
      FRDU(J,2) = AUT   * SVI*SVI*WTGU(J)
      FRDU(J,3) = AUT   * SVJ*SVJ*WTGU(J)
150   CONTINUE
C
      DO 175 J  = 1 ,JE
      SUI       = SINU(J)
      SUJ       = SINU(J+1)
      FRDV(J,1) = AVO   * SINV(J)
      FRDV(J,2) = AVT   * SUJ*SUJ*WTGV(J)
      FRDV(J,3) = AVT   * SUI*SUI*WTGV(J)
175   CONTINUE
C
      DO 200 J  = JB,JE
      SVI       = SINV(J)
      SVJ       = SINV(J-1)
      FPO       = ANTS  / SINU(J)
      FRDP(J,1) = ONE
      FRDP(J,2) = FPO   * SVI*SVI*SVI*WTGU(J)*WTGV(J)
      FRDP(J,3) = FPO   * SVJ*SVJ*SVJ*WTGU(J)*WTGV(J-1)
200   CONTINUE
      FRDP(1 ,1)= + DNSI
      FRDP(NY,1)= - DNSI
C
      CDFS      = DTHDFS * DFS0*DFS0
      DO 250 K  = 1 ,3
      DO 250 J  = 1 ,NY
         FRDU(J,K) = FRDU(J,K) * CDFS
         FRDV(J,K) = FRDV(J,K) * CDFS
         FRDP(J,K) = FRDP(J,K) * CDFS
250   CONTINUE
      DO 260 K  = 1 ,3
      DO 260 J  = 1 ,NY
      FRDP(J,K) = 3.0E0 * FRDP(J,K)
260   CONTINUE


!- check ----------------------------------------------------------
!
!      open (10,file='stdfsc-s.out')
!
!      write(10,*) '----------------- dxvpn -----------------'
!      do j=1,ny
!        write(10,11) j,dxvpn(j)
!      enddo
!
!      write(10,*) '----------------- dxvps -----------------'
!      do j=1,ny
!        write(10,11) j,dxvps(j)
!      enddo
!
!      write(10,*) '----------------- frdt -----------------'
!      write(10,11) (j,(frdt(j,k),k=1,3),j=1,ny)
!
!      write(10,*) '----------------- frds -----------------'
!      write(10,11) (j,(frds(j,k),k=1,3),j=1,ny)
!
!      write(10,*) '----------------- frdu -----------------'
!      write(10,11) (j,(frdu(j,k),k=1,3),j=1,ny)
!
!      write(10,*) '----------------- frdv -----------------'
!      write(10,11) (j,(frdv(j,k),k=1,3),j=1,ny)
!
!      write(10,*) '----------------- frdp -----------------'
!      write(10,11) (j,(frdp(j,k),k=1,3),j=1,ny)
!
!11    format(1x,i4,3e30.20)
!
!      close (10)
!!      stop
!--------------------------------------------------------------(wh)-----

      RETURN
      END

#endif
