#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: geopk --- Calculate geopotential to the kappa
!
!-----------------------------------------------------------------------
! There are two versions of geopk below. The first is the standard
! version. The second (called geopk16) operates in yz space and performs
! semi-global communication in the z direction (to avoid transposes).
! It also uses 16-byte reals to preserve accuracy through round-off.
! Geopk16 can also be run with 8-byte reals (for debugging purposes)
! by setting Dsize equal to r8 (immediately below).
! Note that the interfaces to the two versions are slightly different.
! Also, geopk (the standard version with transposes) is called for
! the D-grid during the last small timestep in cd_core (id=1).

#define DSIZE 16
!#define DSIZE 8

#if (DSIZE == 16)
# define DTWO 2
#else
# define DTWO 1
#endif
!-----------------------------------------------------------------------
!
! !INTERFACE:
      subroutine geopk(ptop, pe, delp, pk, wz, hs, pt, ng_d, im, jm, km, &
                       jfirst, jlast, ifirst, ilast,              &
                       cp, akap, nx, id)

      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid, only: twod_decomp

      implicit none

! !INPUT PARAMETERS:

      integer im, jm, km, ng_d, jfirst, jlast, ifirst, ilast, id
      integer nx                        ! # of pieces in longitude direction
      real(r8)    akap, cp, ptop
      real(r8)    hs(ifirst:ilast,jfirst:jlast)
      real(r8)    pt(ifirst:ilast,jfirst-ng_d:jlast+ng_d,km)
      real(r8)  delp(ifirst:ilast,jfirst:jlast,km)

! !OUTPUT PARAMETERS
      real(r8) wz(ifirst:ilast,jfirst:jlast,km+1)  ! space N*1 S*1
      real(r8) pk(ifirst:ilast,jfirst:jlast,km+1)  ! space N*1 S*1
      real(r8) pe(ifirst:ilast,km+1,jfirst:jlast)  ! only if id .eq. 1

! !DESCRIPTION:
!     Calculates geopotential and pressure to the kappa.  This is an expensive
!     operation and several out arrays are kept around for future use.
!
! !REVISION HISTORY:
!
!  WS  99.10.22: MPIed SJ's original SMP version
!  SJL 00.01.01: Merged C-core and D-core computation
!                SMP "decmposition" in E-W by combining i and j loops
!  WS  00.12.01: Replaced MPI_ON with SPMD; hs now distributed
!  AAM 01.06.27: Generalize for 2D decomposition
!  AAM 01.07.24: Removed dpcheck
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local:
      integer i, j, k
      real(r8)    p1d(ifirst:ilast)
      integer ixj, jp, it, i1, i2, nxu, itot

      real(r8) ptk

      itot = ilast - ifirst + 1
      nxu = nx
      if (twod_decomp .eq. 1) nxu = 1
      it = itot / nxu
      jp = nxu * ( jlast - jfirst + 1 )

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i1, i2, ixj, i, j, k, p1d, ptk)

!     do 2000 j=jfirst,jlast
      do 2000 ixj=1, jp

         j  = jfirst + (ixj-1)/nxu
         i1 = ifirst + it * mod(ixj-1, nxu)
         i2 = i1 + it - 1

         ptk  = ptop ** akap
         do i=i1,i2
            p1d(i) = ptop
            pk(i,j,1) = ptk
            wz(i,j,km+1) = hs(i,j)
         enddo

         if(id .eq. 1) then
            do i=i1,i2
               pe(i,1,j) = ptop
            enddo
         endif

! Top down
         do k=2,km+1
            do i=i1,i2
               p1d(i)  = p1d(i) + delp(i,j,k-1)
#if !defined( VECTOR_MATH )
               pk(i,j,k) = p1d(i) ** akap
            enddo
#else
!
! Use a pipelined library for log and exp
!
            enddo
            call vlog( pk(i1,j,k),p1d(i1), it )
            do i=i1,i2
               pk(i,j,k) = akap * pk(i,j,k)
            enddo
            call vexp( pk(i1,j,k),pk(i1,j,k), it )
#endif
            if(id .eq. 1) then
               do i=i1,i2
                  pe(i,k,j) = p1d(i)
               enddo
            endif
         enddo

! Bottom up
         do k=km,1,-1
            do i=i1,i2
               wz(i,j,k) = wz(i,j,k+1)+cp*pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
            enddo
         enddo
2000  continue

      return
      end
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: geopk16 --- Calculate geopotential to the kappa
!
! !INTERFACE:
      subroutine geopk16(ptop, pe, delp, pk, wz, hs, pt, ng_d, im, jm, km,   &
                       jfirst, jlast, ifirst, ilast, cp, akap,         &
                       kfirst, klast)

      use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
      use decompmodule,only : decomptype
      use pmgrid,only : myid_z, npr_z
#if defined( SPMD )
      use spmd_dyn, only: comm_z
      use parutilitiesmodule, only : parexchangevector
#endif

      implicit none

! !INPUT PARAMETERS:

      integer im, jm, km, ng_d, jfirst, jlast, ifirst, ilast, kfirst, klast
      real(r8)    akap, cp, ptop
      real(r8)    hs(ifirst:ilast,jfirst:jlast)

! !INPUT PARAMETERS:
      real(r8)    pt(ifirst:ilast,jfirst-ng_d:jlast+ng_d,kfirst:klast) 
      real(r8)  delp(ifirst:ilast,jfirst:jlast,kfirst:klast)

! !OUTPUT PARAMETERS
      real(r8) wz(ifirst:ilast,jfirst:jlast,kfirst:klast+1)  ! space N*1 S*1
      real(r8) pk(ifirst:ilast,jfirst:jlast,kfirst:klast+1)  ! space N*1 S*1
      real(r8) pe(ifirst:ilast,kfirst:klast+1,jfirst:jlast)  ! temporary variable

! !DESCRIPTION:
!     Calculates geopotential and pressure to the kappa.  This is an expensive
!     operation and several out arrays are kept around for future use.
!     To preserve accuracy through round-off, 16-byte reals are used
!     for some intermediate data. The variable "id" is assumed to be
!     -1 or 0; the variable dp_check is assumed to be "false".
!
! !REVISION HISTORY:
!
!  AAM 00.12.18: Original version
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local:
      integer i, j, k, nk, ijtot
#if (DSIZE == 16)
      real*16  delp16(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*16  pe16(ifirst:ilast,jfirst:jlast,kfirst:klast+1)
      real*16  inbuf(ifirst:ilast,jfirst:jlast,0:npr_z-1)
      real*16  outbuf(ifirst:ilast,jfirst:jlast,0:npr_z-1)
#else
      real (r8) delp16(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real (r8) pe16(ifirst:ilast,jfirst:jlast,kfirst:klast+1)
      real (r8) inbuf(ifirst:ilast,jfirst:jlast,0:npr_z-1)
      real (r8) outbuf(ifirst:ilast,jfirst:jlast,0:npr_z-1)
#endif
      integer sendcount(0:npr_z-1), recvcount(0:npr_z-1)

! Arrays inbuf8 and outbuf8 are created to fool the compiler
!  into accepting them as calling arguments for parexchangevector.
!  The trickery below equivalences them to inbuf and outbuf
      real (r8) inbuf8(1), outbuf8(1)
      pointer (ptr_inbuf8, inbuf8)
      pointer (ptr_outbuf8, outbuf8)
      integer (i8) locinbuf, locoutbuf

      locinbuf = loc(inbuf)
      locoutbuf = loc(outbuf)
      ptr_inbuf8 = locinbuf
      ptr_outbuf8 = locoutbuf

      ijtot = (jlast-jfirst+1) * (ilast-ifirst+1)

! Top down

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k)
      do k = kfirst,klast
      do j = jfirst,jlast
      do i = ifirst,ilast
         delp16(i,j,k) = delp(i,j,k)
         pe16(i,j,k+1) = 0.
      enddo
      enddo
      enddo

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j)
      do j = jfirst,jlast
      do i = ifirst,ilast
        pe16(i,j,kfirst) = 0.
      enddo
      enddo

! compute partial sums

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k)
      do j = jfirst,jlast
        do k = kfirst+1,klast+1
        do i = ifirst,ilast
          pe16(i,j,k) = pe16(i,j,k-1) + delp16(i,j,k-1)
        enddo
        enddo
      enddo

#if defined( SPMD )
      if (npr_z .gt. 1) then

! communicate upward

        do nk = 0, npr_z-1
          sendcount(nk) = 0
          recvcount(nk) = 0
        enddo

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, nk)
        do nk = myid_z+1, npr_z-1
          do j = jfirst,jlast
          do i = ifirst,ilast
            inbuf(i,j,nk-myid_z-1) = pe16(i,j,klast+1)
          enddo
          enddo
! Double sendcount since quantities are 16-bytes long
          sendcount(nk) = DTWO*ijtot
        enddo

        call parexchangevector(comm_z, sendcount, inbuf8, recvcount, outbuf8)

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k, nk)
        do k = kfirst,klast+1
          do nk = 0, myid_z-1
          do j = jfirst,jlast
          do i = ifirst,ilast
             pe16(i,j,k) = pe16(i,j,k) + outbuf(i,j,nk)
          enddo
          enddo
          enddo
        enddo

      endif
#endif

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k)
      do k = kfirst,klast+1
      do j = jfirst,jlast
      do i = ifirst,ilast
        pe(i,k,j) = pe16(i,j,k) + ptop
        pk(i,j,k) = pe(i,k,j) ** akap
      enddo
      enddo
      enddo

! Bottom up

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k)
      do k = kfirst,klast
      do j = jfirst,jlast
      do i = ifirst,ilast
        delp16(i,j,k) = cp*pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
        pe16(i,j,k) = 0.
      enddo
      enddo
      enddo

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j)
      do j = jfirst,jlast
      do i = ifirst,ilast
        pe16(i,j,klast+1) = 0.
      enddo
      enddo

! compute partial sums

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k)
      do j = jfirst,jlast
        do k = klast,kfirst,-1
        do i = ifirst,ilast
          pe16(i,j,k) = pe16(i,j,k+1) + delp16(i,j,k)
        enddo
        enddo
      enddo

#if defined( SPMD )
      if (npr_z .gt. 1) then

! communicate downward

        do nk = 0, npr_z-1
          sendcount(nk) = 0
          recvcount(nk) = 0
        enddo

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, nk)
        do nk = 0, myid_z-1
          do j = jfirst,jlast
          do i = ifirst,ilast
            inbuf(i,j,nk) = pe16(i,j,kfirst)
          enddo
          enddo
! Double sendcount since quantities are 16-bytes long
          sendcount(nk) = DTWO*ijtot
        enddo

        call parexchangevector(comm_z, sendcount, inbuf8, recvcount, outbuf8)

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k, nk)
        do k = kfirst,klast+1
          do nk = myid_z+1, npr_z-1
          do j = jfirst,jlast
          do i = ifirst,ilast
            pe16(i,j,k) = pe16(i,j,k) + outbuf(i,j,nk-myid_z-1)
          enddo
          enddo
          enddo
        enddo

      endif
#endif

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k)
      do k = kfirst,klast+1
      do j = jfirst,jlast
      do i = ifirst,ilast
        wz(i,j,k) = pe16(i,j,k) + hs(i,j)
      enddo
      enddo
      enddo

      return
!EOC
      end
!-----------------------------------------------------------------------
