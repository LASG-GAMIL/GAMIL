#include <misc.h>
module pft_module
!BOP
!
! !MODULE: dp_coupling --- dynamics-physics coupling module
!
!
! !PUBLIC MEMBER FUNCTIONS:
      public pft2d, pft_cf, rfftmlt, fftfax 
!
! !DESCRIPTION:
!
!      This module provides fast-Fourier transforms
!
!      \begin{tabular}{|l|l|} \hline \hline
!         pft2d     &  \\ \hline
!         pft\_cf   &  \\ \hline
!         rfftmlt   &  \\ \hline
!         fftfax    &  \\ \hline
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   01.01.30   Lin        Integrated into this module
!   01.03.26   Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pft2d --- Two-dimensional fast Fourier transform
!
! !INTERFACE: 
 subroutine pft2d(p, s, damp, im,  jp, ifax, trigs, q1, q2)

! !USES:
 use shr_kind_mod, only: r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
      integer im                   ! Total X dimension
      integer jp                   ! Total Y dimension
      integer ifax(13)             ! FFT internal information
      real(r8)   s(jp)             ! 3-point algebraic filter
      real(r8)  damp(im,jp)        ! FFT damping coefficients
      real(r8) trigs(3*im/2+1)

! !INPUT/OUTPUT PARAMETERS:
      real(r8) q1( im+2, *)        ! Work array
      real(r8) q2(*)               ! Work array
      real(r8)  p(im,jp)           ! Array to be polar filtered

! !DESCRIPTION:
!
!   Perform a two-dimensional fast Fourier transformation.
!
! !REVISION HISTORY:
!   01.01.30   Lin          Put into this module
!   01.03.26   Sawyer       Added ProTeX documentation
!   02.04.05   Sawyer       Integrated newest FVGCM version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real(r8) rsc, bt 
      integer  i, j, n, nj

!Local Auto arrays:
      real(r8) ptmp(0:im+1)
!!!      real(r8) q1(  im+2, jp)
!!!      real(r8) q2( (im+1)*jp )
      integer  jf(jp)

      nj = 0

      do 200 j=1,jp

#if defined ( ALT_PFT )
      if(s(j) > 1.01) then
#else
      if(s(j) > 1.01 .and. s(j) <= 4.) then

         rsc = 1./s(j)
         bt  = 0.5*(s(j)-1.)

         do i=1,im
            ptmp(i) = p(i,j)
         enddo
           ptmp(   0) = p(im,j)
           ptmp(im+1) = p(1 ,j)

         do i=1,im
            p(i,j) = rsc * ( ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)) )
         enddo

      elseif(s(j) > 4.) then
#endif

! Packing for FFT 
           nj  = nj + 1
         jf(nj) = j

         do i=1,im
            q1(i,nj) = p(i,j)
         enddo
            q1(im+1,nj) = 0.
            q1(im+2,nj) = 0.

      endif
200   continue

      if( nj == 0) return
      call rfftmlt(q1,  q2, trigs, ifax, 1, im+2, im, nj, -1)

      do n=1,nj
         do i=5,im+2
            q1(i,n) = q1(i,n) * damp(i-2,jf(n))
         enddo
      enddo

      call rfftmlt(q1, q2, trigs, ifax, 1, im+2, im, nj, 1)

      do n=1,nj
         do i=1,im
            p(i,jf(n)) = q1(i,n)
         enddo
      enddo
      return
!EOC
 end subroutine pft2d
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pft_cf --- Calculate algebraic and FFT polar filters
!
! !INTERFACE: 
 subroutine pft_cf(im, jm, js2g0, jn2g0, jn1g1, sc, se, dc, de,          &
                   cosp, cose, ycrit)

! !USES:
 use shr_kind_mod, only: r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
      integer im                      ! Total X dimension
      integer jm                      ! Total Y dimension
      integer js2g0                   ! j south limit ghosted 0 (SP: from 2)
      integer jn2g0                   ! j north limit ghosted 0 (NP: from jm-1)
      integer jn1g1                   ! j north limit ghosted 1 (starts jm)
      real (r8)   cosp(jm)            ! cosine array
      real (r8)   cose(jm)            ! cosine array
      real (r8)   ycrit               ! critical value

! !OUTPUT PARAMETERS:
      real (r8)   sc(js2g0:jn2g0)     ! Algebric filter at center
      real (r8)   se(js2g0:jn1g1)     ! Algebric filter at edge
      real (r8)   dc(im,js2g0:jn2g0)  ! FFT filter at center
      real (r8)   de(im,js2g0:jn1g1)  ! FFT filter at edge

! !DESCRIPTION:
!
!   Compute coefficients for the 3-point algebraic and the FFT
!   polar filters.
!
! !REVISION HISTORY:
!
!   99.01.01   Lin          Creation
!   99.08.20   Sawyer/Lin   Changes for SPMD mode
!   01.01.30   Lin          Put into this module
!   01.03.26   Sawyer       Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j
      real (r8)   dl, coszc, cutoff, phi, damp
      real (r8)  pi

      pi = 4.d0 * datan(1.d0)

      coszc = cos(ycrit*pi/180.)

! INIT fft polar coefficients:
      dl = pi/dble(im)
      cutoff = 1.e-20

      do j=js2g0,jn2g0
         do i=1,im
            dc(i,j) = 1.
         enddo
      enddo

      do j=js2g0,jn1g1
         do i=1,im
            de(i,j) = 1.
         enddo
      enddo

!     write(6,*) '3-point polar filter coefficients:'

!************
! Cell center
!************
      do j=js2g0,jn2g0
            sc(j) = (coszc/cosp(j))**2

#if defined ( ALT_PFT )
         if( sc(j) > 1. ) then
#else
         if(sc(j) > 1. .and. sc(j) <= 2.0) then
            sc(j) =  1. +  (sc(j)-1.)/(sc(j)+1.)
         elseif(sc(j) > 2. .and. sc(j) <= 4.) then
            sc(j) =  1. +  sc(j)/(8.-sc(j))
         elseif(sc(j) > 4. ) then
#endif

! FFT filter
            do i=1,im/2
               phi = dl * i
               damp = min((cosp(j)/coszc)/sin(phi), 1.)**2
               if(damp < cutoff) damp = 0.
               dc(2*i-1,j) = damp
               dc(2*i  ,j) = damp
            enddo

         endif
      enddo

!************
! Cell edges
!************
      do j=js2g0,jn1g1
            se(j) = (coszc/cose(j))**2

#if defined ( ALT_PFT )
         if( se(j) > 1. ) then
#else
         if(se(j) > 1. .and. se(j) <= 2.0 ) then
            se(j) =  1. +  (se(j)-1.)/(se(j)+1.)
         elseif(se(j) > 2. .and. se(j) <= 4.) then
            se(j) =  1. +  se(j)/(8.-se(j))
         elseif(se(j) > 4. ) then
#endif
! FFT
            do i=1,im/2
               phi = dl * i
               damp = min((cose(j)/coszc)/sin(phi), 1.)**2
               if(damp < cutoff) damp = 0.
               de(2*i-1,j) = damp
               de(2*i  ,j) = damp
            enddo
         endif
      enddo
      return
!EOC
 end subroutine pft_cf
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: rfftmlt --- Apply fast Fourier transform
!
! !INTERFACE: 
 subroutine rfftmlt(a,work,trigs,ifax,inc,jump,n,lot,isign)

! !USES:
 use shr_kind_mod, only: r8 => shr_kind_r8

! !DESCRIPTION:
!
!   Apply the fast Fourier transform.  If CPP token SGI_FFT is
!   set, SGI libraries will be used.  Otherwise the Fortran code
!   is inlined.
!
! !REVISION HISTORY:
!
!   99.11.24   Sawyer       Added wrappers for SGI
!   01.03.26   Sawyer       Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC


#if defined( SGI_FFT )
!
! WS 99.11.24 : Added SGI wrappers
!
      implicit none
      integer inc, jump, n, lot, isign
      real(r8)   a(jump,lot)
      real(r8)   trigs(1)
      real(r8)   work(1)                       ! Not used; here for plug reason
      integer ifax(*)
! Local
      integer*4 iisign,in,iinc,ijump,ilot
      integer i, j
      real(r8) scale

!-----convert to i4
      iisign = isign
      iinc = inc
      ijump = jump
      in = n
      ilot = lot

      if( iisign < 0 ) then
!-----forward
          call dzfftm1du (iisign,in,ilot,a,iinc,ijump,trigs)
       endif

      if( iisign > 0 ) then
!-----backward
          call zdfftm1du (iisign,in,ilot,a,iinc,ijump,trigs)

          scale = 1.0/float(n)
          do j=1,lot
             do i=1,jump
                a(i,j) = scale*a(i,j)
             enddo
          enddo
       endif
#else
!
! Default f77 version
!
      real(r8) a(jump*lot)
      real(r8) work((n+1)*lot)
      real(r8) trigs(3*n/2+1)
      integer ifax(13)
 
!     SUBROUTINE "FFT991" - MULTIPLE REAL/HALF-COMPLEX PERIODIC
!     FAST FOURIER TRANSFORM
 
!     SAME AS FFT99 EXCEPT THAT ORDERING OF DATA CORRESPONDS TO
!     THAT IN MRFFT2
 
!     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
!     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
!     (1970), 315-337)
 
!     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
!     WORK IS AN AREA OF SIZE (N+1)*LOT
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
 
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(0),X(1),X(2),...,X(N-1)
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!     PARALLEL
 
!     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
 
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
 
      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign==+1) go to 30
 
!     IF NECESSARY, TRANSFER DATA TO WORK AREA
      igo=50
      if (mod(nfax,2)==1) goto 40
      ibase=1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
!DIR$ IVDEP
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
 
      igo=60
      go to 40
 
!     PREPROCESSING (ISIGN=+1)
 
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
 
!     COMPLEX TRANSFORM
 
   40 continue
      ia=1
      la=1
      do 80 k=1,nfax
      if (igo==60) go to 60
   50 continue
      call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,     &
                  ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,     &
                  2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
 
      if (isign==-1) go to 130
 
!     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      if (mod(nfax,2)==1) go to 110
      ibase=1
      jbase=1
      do 100 l=1,lot
      i=ibase
      j=jbase
!DIR$ IVDEP
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
 
!     FILL IN ZEROS AT END
  110 continue
      ib=n*inc+1
!DIR$ IVDEP
      do 120 l=1,lot
      a(ib)=0.0
      a(ib+inc)=0.0
      ib=ib+jump
  120 continue
      go to 140
 
!     POSTPROCESSING (ISIGN=-1):
 
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
 
  140 continue
#endif
      return
!EOC
 end subroutine rfftmlt
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fftfax --- Initialize FFT
!
! !INTERFACE: 
 subroutine fftfax (n, ifax, trigs)

! !USES:

 use shr_kind_mod, only: r8 => shr_kind_r8

! !DESCRIPTION:
!
!   Initialize the fast Fourier transform.  If CPP token SGI_FFT is
!   set, SGI libraries will be used.  Otherwise the Fortran code
!   is inlined.
!
! !REVISION HISTORY:
!
!   99.11.24   Sawyer       Added wrappers for SGI
!   01.03.26   Sawyer       Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC

#if defined( SGI_FFT )
      implicit none
      real(r8)    trigs(1)
      integer ifax(*)
      integer n
! local
      integer*4 nn

      nn=n
      call dzfftm1dui (nn,trigs)
#else
       integer ifax(13)
       real(r8) trigs(3*n/2+1)
 
! MODE 3 IS USED FOR REAL/HALF-COMPLEX TRANSFORMS.  IT IS POSSIBLE
! TO DO COMPLEX/COMPLEX TRANSFORMS WITH OTHER VALUES OF MODE, BUT
! DOCUMENTATION OF THE DETAILS WERE NOT AVAILABLE WHEN THIS ROUTINE
! WAS WRITTEN.
 
      data mode /3/
      call fax (ifax, n, mode)
      i = ifax(1)
      if (ifax(i+1) > 5 .or. n <= 4) ifax(1) = -99
      if (ifax(1) <= 0 ) then
        write(6,*) ' set99 -- invalid n'
        stop 'set99'
      endif
      call fftrig (trigs, n, mode)
#endif
      return
!EOC
 end subroutine fftfax
!-----------------------------------------------------------------------

end module pft_module
