! (wanhui 2003.04.01)
! -------------------
	SUBROUTINE  LATMESH(DY,YTHU,YTHV,WTGU,WTGV)

	IMPLICIT NONE

	INCLUDE 'PARADYN_serial'

	INTEGER NE,NG,I,J,NEE,NGG
	PARAMETER(NE=66,NG=64,NEE=6,NGG=9)

	REAL*8      YTHU1(NE),YTHV1(NE),WTGU1(NE),WTGV1(NE) !FOR E-GRID,LOCAL
	REAL*8      YTHU2(NG),  SG (NG),   W2(NG)           !FOR G-GRID,LOCAL
	REAL*8   DY,YTHU (NY),YTHV (NY),WTGU (NY),WTGV (NY) !FOR H-GRID,OUT
	REAL*8      B(NY,NY),M(NY)                          !FOR SPLINE,LOCAL

        write(6,*)
        write(6,*) '!!!  hybrid latmesh  !!!'
        write(6,*)

!-------------- EVEN-AREA GRID -----------------
	CALL LATMESH0(DY,YTHU1,YTHV1,WTGU1,WTGV1)

!-------------- GAUSSIAN GRID ------------------
	CALL GAUAW(SG,W2,NG)
	  DO J=1,NG/2
	    YTHU2(J)     =-ASIN(SG(J))
	    YTHU2(NG+1-J)= ASIN(SG(J))
	  ENDDO

!--------------- HYBRID GRID -------------------
	DO J=1,NEE
	  YTHU(J)     =YTHU1(J)     -PI/2.0d0
	  YTHU(NY+1-J)=YTHU1(NE+1-J)-PI/2.0d0
	ENDDO
	DO J=NEE+1,NY/2
	  YTHU(J)     =YTHU2(NGG-1+J-NEE)
	  YTHU(NY+1-J)=YTHU2(NG+1-(NGG-1+J-NEE))
	ENDDO


!--------- CUBIC SPLINE INTERPOLATION ----------
	DY=PI/(NE-1)
	DO J=1,NY
        DO I=1,NY
	  B(I,J)=0.0d0
	ENDDO
	ENDDO

	  B(1,1)    =2.0d0
	  B(2,1)    =1.0d0
	  B(NY-1,NY)=1.0d0
	  B(NY,NY)  =2.0d0
	  M(1)      =(YTHU(2)-YTHU(1)-DY/WTGU1(1))*6.0d0/DY/DY
	  M(NY)     =(DY/WTGU1(NE)-YTHU(NY)+YTHU(NY-1))*6.0d0/DY/DY

	DO J=2,NY-1
	  B(J-1,J) =1.0d0
	  B(  J,J) =4.0d0
	  B(J+1,J) =1.0d0
	  M(J)     =(YTHU(J+1)-2.0d0*YTHU(J)+YTHU(J-1))*6.0d0/DY/DY
	ENDDO

	CALL EQUATION_SOLVER(B,NY,M)

	  WTGU(1)=WTGU1(1)
	  YTHV(1)=(YTHU(1)+YTHU(2))/2.0d0-DY*DY/16.0d0*(M(1)+M(2))
	  WTGV(1)=1.0d0/((YTHU(2)-YTHU(1))/DY-DY/24.0d0*(M(2)-M(1)))

	DO J=2,NY-1
	  YTHV(J)=(YTHU(J)+YTHU(J+1))/2.0d0-DY*DY/16.0d0*(M(J)+M(J+1))
	  WTGU(J)=1.0/((YTHU(J)-YTHU(J-1))/DY+DY/6.0*(2.0*M(j)+M(J-1)))
	  WTGV(J)=1.0/((YTHU(J+1)-YTHU(J))/DY-DY/24.0*(M(J+1)-M(J)))
	ENDDO

	  WTGU(NY)=WTGU1(NE)
       	  YTHV(NY)=YTHU(NY)
	  WTGV(NY)=WTGU(NY)

!------------------------------------------
        DO J=1,NY
	  YTHU(J)=YTHU(J)+PI/2.0d0
	  YTHV(J)=YTHV(J)+PI/2.0d0
	ENDDO

!!	write(*,'(1x,i4,4f20.16)') (j,ythu(j),ythv(j),wtgu(j),wtgv(j),j=1,ny)

	RETURN
	END



!###############################################
	SUBROUTINE LATMESH0(DY,YTHU,YTHV,WTGU,WTGV)
!
!
        IMPLICIT NONE
        INCLUDE 'PARADYN_serial'
!!	REAL PI
!!	PARAMETER(PI=3.141592653589793)
	INTEGER NY0
	PARAMETER(NY0=66)
!
!     1) INPUT VARIABLE: NONE
!     2) INPUT CONSTANT: NY0,PI
!		   NY IS THE GRID NUMBER FOR LATITUDE, DEFINED IN 'PARADYN'
!		   PI IS THE RATIO OF THE CIRCUMFERENCE OF A CIRCLE TO
!			  ITS DIAMETER, DEFINED IN 'PARADYN'
!
!     3) OUTPUT VARIABLES: DY, YTHU, YTHV, WTGU, WTGV
!
	REAL*8
     _       DY		! MERIDIONAL STEPSIZE
     _      ,YTHU(NY0)	! LATITUDE AT THE NORNAL MERIDIONAL GRID Yj
     _      ,YTHV(NY0) ! LATITUDE AT THE HALF-MOVE MERIDIONAL GRID Yj+1/2
     _      ,WTGU(NY0)	! AREA-WEIGHTING AT THE NORNAL MERIDIONAL GRID Yj
     _      ,WTGV(NY0)	! AREA-WEIGHTING AT THE HALF-MOVE MERIDIONAL GRID Yj+1/2
!
!     4) WORKING VARIABLES: U,S,DA,A2,U2,A,B,J,M1
!
	REAL*8
     _       AS ! AREA SIZE WITH RESPECT TO LATITUDE
     _      ,B  !	A CONTROL PARAMETER RELATED TO THE MERIDIONAL RESOLUTION
!
!     The resolution formula:
!
!           DY = (180-10B)/(NY0-1) (in degree) PI(1-B/18)/(NY0-1) (in arc)
!
!           when B=2.0, NY0=41, DY=4 degree; when B=2.0, NY0=81, DY=2 degree
!
     _      ,A  ! A  = B/(0.5*PI)
!				A CONTROL PARAMETER RELATED TO THE AREA-SIZE COMPUTING
     _      ,A2 ! A2 = A*2
     _      ,DA ! DA = 1/A
     _      ,S  ! S  = AS(PI/2)	 THE TOTAL AREA-SIZE OF THE WHOLE REGION
     _      ,S1 ! S1 = AS(-PI/3)	 THE AREA-SIZE AT THE POINT THETA=-PI/3
     _      ,S2 ! S2 = AS(PI/3)	 THE AREA-SIZE AT THE POINT THETA=PI/3
     _      ,U1 ! U1 = 1+2B/3		 NEEDED WHEN CALCULATE YTHU, YTHV
     _      ,U2 ! U2 = (1-B/3)^2	 NEEDED WHEN CALCULATE YTHU, YTHV
	INTEGER
     _       J  ! AN INTEGER VARIABLE FOR LOOPS IN THIS SUBROUTINE
     _      ,M1 ! M1 = NY0-1
!
!	Description of the even-area method	developed by Bin Wang in 2000:
!
!     Under the weighting w(theta)=1.0-a(|theta|-PI/3) when |theta|>=PI/3
!                         w(theta)=1.0                 when |theta|< PI/3
!
!     The area size AS is:
!     AS(theta)=[(1+a*theta+a*PI/3)^2-(1-b/3)^2]/(2a)  when -PI/2<=theta<=-PI/3
!     AS(theta)=AS(-PI/3)+(theta+PI/3)                 when -PI/3< theta<= PI/3
!     AS(theta)=AS(PI/3)+[1-(1-a*theta+a*PI/3)^2]/(2a) when  PI/3< theta<= PI/2
!     where a= b/(0.5*Pi), 0<b<1, b=1.0-0.25*sqrt(25/(m1-20))
!
!	Suppose the total area-size of the whole region S is partitioned into NY0-1
!	equal small area: AS(theta(j+1))-AS(theta(j))=DY=constant, then theta(j)
!	can be calculated according to the formula of AS. Obviously,
!     theta(j+1)-theta(j) will not be a constant when they are not in the interval
!     [-PI/3, PI/3]. Especially, closer to poles theta(j) is, biger
!     theta(j+1)-theta(j) will become. In this way, the physical stepsizes in the
!	polar regions increase and the computational stability becomes better.
!	Note that: the physical mesh is not even, but the computing mesh is, which
!	makes the meridional discretization easy.
!
	M1= NY0-1
C
        B = 2.0
	A  = B*2.0D0/PI
	A2 = A*2.0D0
	DA = 1.0D0/A
        S  = PI*(1.0D0-B/18.0D0)
        S1 = PI*(1.0D0-B/6.0D0)/6.0D0
        S2 = PI*(5.0D0-B/6.0D0)/6.0D0
	U1 =  1.0D0+2.0D0*B/3.0D0
	U2 = (1.0D0-B/3.0D0)*(1.0D0-B/3.0D0)
	DY = S/DFLOAT(M1)
	DO J=0,M1
	   AS = DY*DFLOAT(J)
	   IF (AS.LE.S1) THEN
!!	      YTHU(J+1)=(DSQRT(AS*A2+U2)-U1)*DA+PI*0.5D0
	      YTHU(J+1)=( SQRT(AS*A2+U2)-U1)*DA+PI*0.5D0
              IF (YTHU(J+1).LT.0.0) YTHU(J+1)=0.0D0
	      WTGU(J+1)=1.0D0-A*(DABS(YTHU(J+1)-PI*0.5D0)-PI/3.0D0)
	   ELSE IF (AS.LE.S2) THEN
	      YTHU(J+1)=AS-S1-PI/3.0D0+PI*0.5D0
	      WTGU(J+1)=1.0D0
	   ELSE
!!	      YTHU(J+1)=(U1-DSQRT(1.0D0-(AS-S2)*A2))*DA+PI*0.5D0
	      YTHU(J+1)=(U1- SQRT(1.0D0-(AS-S2)*A2))*DA+PI*0.5D0
	      WTGU(J+1)=1.0D0-A*(DABS(YTHU(J+1)-PI*0.5D0)-PI/3.0D0)
	   END IF
	   AS = DY*(DFLOAT(J)+0.5D0)
	   IF (AS.LE.S1) THEN
!!	      YTHV(J+1)=(DSQRT(AS*A2+U2)-U1)*DA+PI*0.5D0
	      YTHV(J+1)=( SQRT(AS*A2+U2)-U1)*DA+PI*0.5D0
	      WTGV(J+1)=1.0D0-A*(DABS(YTHV(J+1)-PI*0.5D0)-PI/3.0D0)
	   ELSE IF (AS.LE.S2) THEN
	      YTHV(J+1)=AS-S1-PI/3.0+PI*0.5
	      WTGV(J+1)=1.0D0
	   ELSE
!!	      YTHV(J+1)=(U1-DSQRT(1.0D0-(AS-S2)*A2))*DA+PI*0.5D0
	      YTHV(J+1)=(U1- SQRT(1.0D0-(AS-S2)*A2))*DA+PI*0.5D0
	      WTGV(J+1)=1.0D0-A*(DABS(YTHV(J+1)-PI*0.5D0)-PI/3.0D0)
	   END IF
	END DO
!
	YTHV(NY0)=PI
	WTGV(NY0)=1.0D0-A*(PI*0.5-PI/3.0)
!
	RETURN
	END




!#############################################
      subroutine gauaw(a       ,w       ,k)
C-----------------------------------------------------------------------
C
C Calculate sine of latitudes a(k) and weights w(k) for the gaussian
C quadrature. The algorithm is described in Davis and Rabinowitz,
C Journal of Research of the NBS, V 56, Jan 1956.
C The zeros of the bessel function j0, which are obtained from bsslzr,
C are used as a first guess for the abscissa.
C
C Modified 1/23/97 by Jim Rosinski to use real*16 arithmetic in order to
C achieve (nearly) identical weights and latitudes on all machines.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      L. Bath, Jun 1992
C                    L. Buja, Feb 1996
C Reviewed:          D. Williamson, J. Hack, Aug 1992
C                    D. Williamson, J. Hack, Feb 1996
C
C-----------------------------------------------------------------------
c
c $Id: gauaw.F,v 1.5 1998/07/30 23:56:58 rosinski Exp $
c $Author: rosinski $
c
C-----------------------------------------------------------------------
      implicit none
C------------------------------Arguments--------------------------------
C
C Input argument
C
      integer k            ! number of latitudes pole to pole
C
C Output arguments
C
      real*8 a(k)            ! sine of latitudes
      real*8 w(k)            ! gaussian weights
C
C---------------------------Local workspace-----------------------------
C
      real*8 sinlat(k)    ! sine of latitudes
      real*8 wgt(k)       ! gaussian weights
      real*8 one          ! 1. in real*16.  Needed by atan
      real*8 eps          ! convergence criterion
      real*8 pi           ! value of pi
      real*8 c            ! constant combination
      real*8 fk           ! real k
      real*8 xz           ! abscissa estimate
      real*8 pkm1         ! |
      real*8 pkm2         ! |-polynomials
      real*8 pkmrk        ! |
      real*8 pk           ! |
      real*8 sp           ! current iteration latitude increment
      real*8 avsp         ! |sp|
      real*8 fn           ! real n
      parameter (one = 1.)
      parameter (eps = 1.D-10)

      integer kk           ! k/2 (number of latitudes in hemisphere)
      integer is           ! latitude index
      integer iter         ! iteration counter
      integer n,l          ! indices
C
C------------------------------Externals--------------------------------
C
      external bsslzr      ! provides zeroes of Bessel function or
C                          ! estimates thereof
C
C-----------------------------------------------------------------------
C
      pi  = 4.*atan(one)
C
C The value eps, used for convergence tests in the iterations,
C can be changed.  Newton iteration is used to find the abscissas.
C
      c = (1.-(2./pi)**2)*0.25
      fk = k
      kk = k/2
      call bsslzr(sinlat,kk)
      do is=1,kk
        xz = cos(sinlat(is)/sqrt((fk+0.5)**2+c))
C
C This is the first approximation to xz
C
        iter = 0
   10   pkm2 = 1.
        pkm1 = xz
        iter = iter + 1
        if (iter.gt.10) then
C
C Error exit
C
          write(6,*)'GAUAW:Error exit,no convergence in 10 iterations'
	    stop
        end if
C
C Computation of the legendre polynomial
C
        do n=2,k
          fn = n
          pk = ((2.*fn-1.)*xz*pkm1-(fn-1.)*pkm2)/fn
          pkm2 = pkm1
          pkm1 = pk
        enddo
        pkm1 = pkm2
        pkmrk = (fk*(pkm1-xz*pk))/(1.-xz**2)
        sp = pk/pkmrk
        xz = xz - sp
        avsp = abs(sp)
        if (avsp.gt.eps) go to 10
        sinlat(is) = xz
        wgt(is) = (2.*(1.-xz**2))/(fk*pkm1)**2
      end do
C
      if (k.ne.kk*2) then
C
C For odd k computation of weight at the equator
C
        sinlat(kk+1) = 0.
        pk = 2./fk**2
        do n=2,k,2
          fn = n
          pk = pk*fn**2/(fn-1.)**2
        end do
        wgt(kk+1) = pk
      end if
C
C Complete the sets of abscissas and weights, using the symmetry.
C Also note truncation from real*16 to real*8
C
      do n=1,kk
        l = k + 1 - n
        a(n) = sinlat(n)
        a(l) = -sinlat(n)

        w(n) = wgt(n)
        w(l) = wgt(n)
      end do
      return
      end



!################################
      subroutine bsslzr(bes,n)
C-----------------------------------------------------------------------
C
C Return n zeros (or if n>50, approximate zeros), of the Bessel function
C j0,in the array bes. The first 50 zeros will be given exactly, and the
C remaining zeros are computed by extrapolation,and therefore not exact.
C
C Modified 1/23/97 by Jim Rosinski to use real*16 arithmetic
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Hack, D. Williamson, August 1992
C Reviewed:          J. Hack, D. Williamson, April 1996
C
C-----------------------------------------------------------------------
      implicit none
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer n              ! Number of zeros to return
C
C Output arguments
C
      real*8 bes(n)         ! Array containing zeros of j0
C
C---------------------------Local workspace-----------------------------
C
      real*8 one            ! 1.
      real*8 pi             ! 3.14.......
      real*8 bz(50)         ! table of first 50 zeros
      save bz                ! ensure re-entrancy
      parameter (one = 1.)
      integer j,nn           ! loop indices
C
      data bz           / 2.4048255577,   5.5200781103,
     $    8.6537279129,  11.7915344391,  14.9309177086,  18.0710639679,
     $   21.2116366299,  24.3524715308,  27.4934791320,  30.6346064684,
     $   33.7758202136,  36.9170983537,  40.0584257646,  43.1997917132,
     $   46.3411883717,  49.4826098974,  52.6240518411,  55.7655107550,
     $   58.9069839261,  62.0484691902,  65.1899648002,  68.3314693299,
     $   71.4729816036,  74.6145006437,  77.7560256304,  80.8975558711,
     $   84.0390907769,  87.1806298436,  90.3221726372,  93.4637187819,
     $   96.6052679510,  99.7468198587, 102.8883742542, 106.0299309165,
     $  109.1714896498, 112.3130502805, 115.4546126537, 118.5961766309,
     $  121.7377420880, 124.8793089132, 128.0208770059, 131.1624462752,
     $  134.3040166383, 137.4455880203, 140.5871603528, 143.7287335737,
     $  146.8703076258, 150.0118824570, 153.1534580192, 156.2950342685/
C
      pi = 4.*atan(one)
      nn = n
      if (n.gt.50) then
        bes(50) = bz(50)
        do j=51,n
          bes(j) = bes(j-1) + pi
        end do
        nn = 49
      end if
      do j=1,nn
        bes(j) = bz(j)
      end do
      return
      end



!############################################
	SUBROUTINE EQUATION_SOLVER (A,N,B)
!
	IMPLICIT NONE
	INTEGER N,L,IS,JS(N),I,J,K
	REAL*8  A(N,N),B(N),T,D
!
	L=1
!
	DO 100 K=1,N
!
	   D=0.0
	   DO 10 I=K,N
	   DO 10 J=K,N
	     IF (ABS(A(I,J)).GT.D) THEN
	       D=ABS(A(I,J))
		   JS(K)=J
		   IS=I
	     ENDIF
10       CONTINUE
	   IF (D+1.0.EQ.1.0) THEN
	     WRITE(*,'(20A)')  'ATTENTION: FAIL!'
	     L=0
		 RETURN
         ENDIF
!
         DO 30 J=K,N
		 T=A(K,J)
	     A(K,J)=A(IS,J)
		 A(IS,J)=T
30	   CONTINUE
!
	   T=B(K)
	   B(K)=B(IS)
	   B(IS)=T
!
	   DO 50 I=1,N
		T=A(I,K)
	    A(I,K)=A(I,JS(K))
		A(I,JS(K))=T
50	   CONTINUE
!
	   T=A(K,K)
	   DO 60 J=K+1,N
	     IF (A(K,J).NE.0.0) A(K,J)=A(K,J)/T
60	   CONTINUE
!
	   B(K)=B(K)/T
	   DO 80 J=K+1,N
	     IF (A(K,J).NE.0.0) THEN
	        DO 70 I=1,N
	          IF ((I.NE.K).AND.(A(I,K).NE.0.0)) THEN
	          A(I,J)=A(I,J)-A(I,K)*A(K,J)
			  ENDIF
70	        CONTINUE
	     ENDIF
80	   CONTINUE
!
	   DO 90 I=1,N
	     IF ((I.NE.K).AND.(A(I,K).NE.0.0)) THEN
	       B(I)=B(I)-A(I,K)*B(K)
	     ENDIF
90	   CONTINUE
!
100	CONTINUE
!
!
      DO 110 K=N,1,-1
	  IF (K.NE.JS(K)) THEN
	    T=B(K)
	    B(K)=B(JS(K))
	    B(JS(K))=T
	  ENDIF
110	CONTINUE
!
	RETURN
	END
