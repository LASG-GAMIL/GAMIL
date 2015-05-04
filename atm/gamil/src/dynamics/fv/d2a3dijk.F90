#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: d2a3ijk -- Generalized 2nd order D-to-A grid transform (3D)
!                      Output array is i,j,k
!
! !INTERFACE:

      subroutine d2a3dijk(u, v, ua, va, im, jm, ktot,                   & 
                          jfirst, jlast, ng_d, ngus, ngun, ngvs, ngvn,  &
                          ifirst, ilast, coslon, sinlon)

! !USES:

      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid, only : twod_decomp, myid_y, myidxy_y, myidxy_x, &
                         nprxy_x, iam

#if defined( SPMD )
      use spmd_dyn, only: comm_y, commxy_y, commxy_x
      use parutilitiesmodule, only : parcollective3d, sumop
      use mod_comm, only: mp_send3d, mp_recv3d
#endif

      implicit none
! !INPUT PARAMETERS:
      integer, intent(in)  :: im      ! Dimensions longitude (total)
      integer, intent(in)  :: jm      ! Dimensions latitude (total)
      integer, intent(in)  :: ktot    ! Dimensions vertical (strip)
      integer, intent(in)  :: jfirst  ! latitude strip start
      integer, intent(in)  :: jlast   ! latitude strip finish
      integer, intent(in)  :: ng_d    ! ghost latitudes on D grid
      integer, intent(in)  :: ngus    ! ghost latitudes on U south
      integer, intent(in)  :: ngun    ! ghost latitudes on U north
      integer, intent(in)  :: ngvs    ! ghost latitudes on V south
      integer, intent(in)  :: ngvn    ! ghost latitudes on V north
      integer, intent(in)  :: ifirst  ! longitude strip start
      integer, intent(in)  :: ilast   ! longitude strip finish
      real(r8), intent(in) :: u(ifirst:ilast,jfirst-ngus:jlast+ngun,ktot) ! U-Wind ghosted N1
      real(r8), intent(in) :: v(ifirst:ilast,jfirst-ngvs:jlast+ngvn,ktot) ! V-Wind
      real(r8) coslon(im),sinlon(im)  ! Sine and cosine in longitude

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: ua(ifirst:ilast,jfirst:jlast,ktot) ! U-Wind
      real(r8), intent(inout) :: va(ifirst:ilast,jfirst:jlast,ktot) ! V-Wind


! !DESCRIPTION:
!
!     This routine performs a second order 
!     interpolation of three-dimensional wind
!     fields on a D grid to an A grid.  !
!
! !REVISION HISTORY:
!     WS  00.12.22 : Creation from d2a3d
!     AAM 01.06.13 : Generalized to 2D decomposition
!     WS  02.04.25 : Newest mod_comm interfaces
!
!EOP
!-----------------------------------------------------------------------
!BOC

      integer imh, i, j, k, itot, jtot, ltot, lbegin, lend, ik
      real(r8) un(ktot), vn(ktot), us(ktot), vs(ktot)
      real(r8) veast(jfirst:jlast,ktot),unorth(ifirst:ilast,ktot)
      real(r8) uvaglob(im,ktot,4),uvaloc(ifirst:ilast,ktot,4)
      real(r8) uaglob(im),vaglob(im)

#if defined( SPMD )
      integer dest, src, incount, outcount
#endif

      itot = ilast-ifirst+1
      jtot = jlast-jfirst+1
      imh = im/2

#if defined( SPMD )
! Set ua on A-grid
      call mp_send3d( iam-nprxy_x, iam+nprxy_x, im, jm, ktot,          &
                      ifirst, ilast, jfirst-ngus, jlast+ngun, 1, ktot, &
                      ifirst, ilast, jfirst, jfirst, 1, ktot, u )
      call mp_recv3d( iam+nprxy_x, im, jm, ktot,                       &
                      ifirst, ilast, jlast+1, jlast+1, 1, ktot,        &
                      ifirst, ilast, jlast+1, jlast+1, 1, ktot, unorth )
                       
      if ( jlast .lt. jm ) then
!$omp  parallel do private(i, k)

         do k=1,ktot
            do i=ifirst,ilast
               ua(i,jlast,k) = 0.5 * ( u(i,jlast,k) + unorth(i,k) )
            enddo
         enddo
      endif
#endif

!$omp  parallel do private(i,j,k)

      do k=1,ktot
        do j=jfirst, jlast-1
          do i=ifirst,ilast
            ua(i,j,k) = 0.5*(u(i,j,k) + u(i,j+1,k))
          enddo
        enddo
      enddo

! Set va on A-grid

!$omp  parallel do private(j,k)

      do k = 1,ktot
         do j=jfirst,jlast
            veast(j,k) = v(ifirst,j,k)
         enddo
      enddo

#if defined( SPMD )
      if (itot .ne. im) then
         dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
         src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
         call mp_send3d( dest, src, im, jm, ktot,                           &
                         ifirst, ilast, jfirst-ngvs, jlast+ngvn, 1, ktot,   &
                         ifirst, ifirst, jfirst, jlast, 1, ktot, v )
         call mp_recv3d( src, im, jm, ktot,                                 & 
                         ilast+1, ilast+1, jfirst, jlast, 1, ktot,          &
                         ilast+1, ilast+1, jfirst, jlast, 1, ktot, veast )
      endif
#endif

!$omp  parallel do private(i,j,k)

      do k=1,ktot
         do j=jfirst, jlast
            do i=ifirst,ilast-1
               va(i,j,k) = 0.5*(v(i,j,k) + v(i+1,j,k))
            enddo
            va(ilast,j,k) = 0.5*(v(ilast,j,k) + veast(j,k))
         enddo
      enddo

!$omp  parallel do private(i,ik,k)

      do ik=1,4
         do k=1,ktot
            do i=1,im
               uvaglob(i,k,ik) = 0.
            enddo
         enddo
      enddo

      if (jfirst .eq. 1) then
!$omp  parallel do private(i,k)
         do k = 1,ktot
            do i=ifirst,ilast
               uvaloc(i,k,1) = ua(i,2,k)
               uvaloc(i,k,2) = va(i,2,k)
               uvaglob(i,k,1) = ua(i,2,k)
               uvaglob(i,k,2) = va(i,2,k)
            enddo
         enddo
         lbegin = 1
         lend = 2
      endif

      if (jlast .eq. jm) then
!$omp  parallel do private(i,k)
         do k = 1,ktot
            do i=ifirst,ilast
               uvaloc(i,k,3) = ua(i,jm-1,k)
               uvaloc(i,k,4) = va(i,jm-1,k)
               uvaglob(i,k,3) = ua(i,jm-1,k)
               uvaglob(i,k,4) = va(i,jm-1,k)
            enddo
         enddo
         lbegin = 3
         lend = 4
      endif
      if (jtot .eq. jm) lbegin=1

#if defined( SPMD )
      if (itot .ne. im) then
         ltot = lend-lbegin+1
         if (jfirst .eq. 1 .or. jlast .eq. jm) then
            call parcollective3d(commxy_x, sumop, im, ktot, ltot, uvaglob(1,1,lbegin))
         endif
      endif
#endif

      if ( jfirst .eq. 1 ) then
! Projection at SP
!$omp  parallel do private(i,k,uaglob,vaglob)
         do k=1,ktot
            us(k) = 0.
            vs(k) = 0.
            do i=1,imh
               us(k) = us(k) + (uvaglob(i+imh,k,1)-uvaglob(i,k,1))*sinlon(i)  &
                     + (uvaglob(i,k,2)-uvaglob(i+imh,k,2))*coslon(i)
               vs(k) = vs(k) + (uvaglob(i+imh,k,1)-uvaglob(i,k,1))*coslon(i)  & 
                     + (uvaglob(i+imh,k,2)-uvaglob(i,k,2))*sinlon(i)
            enddo

            us(k) = us(k)/im
            vs(k) = vs(k)/im
            do i=1,imh
               uaglob(i)   = -us(k)*sinlon(i) - vs(k)*coslon(i)
               vaglob(i)   =  us(k)*coslon(i) - vs(k)*sinlon(i)
               uaglob(i+imh) = -uaglob(i)
               vaglob(i+imh) = -vaglob(i)
            enddo
            do i=ifirst,ilast
               ua(i,1,k) = uaglob(i)
               va(i,1,k) = vaglob(i)
            enddo
         enddo
      endif

      if ( jlast .eq. jm ) then
! Projection at NP
!$omp  parallel do private(i,k,uaglob,vaglob)
         do k=1,ktot
            un(k) = 0.
            vn(k) = 0.
            do i=1,imh
               un(k) = un(k) + (uvaglob(i+imh,k,3)-uvaglob(i,k,3))*sinlon(i) &
                     + (uvaglob(i+imh,k,4)-uvaglob(i,k,4))*coslon(i)
               vn(k) = vn(k) + (uvaglob(i,k,3)-uvaglob(i+imh,k,3))*coslon(i) &
                     + (uvaglob(i+imh,k,4)-uvaglob(i,k,4))*sinlon(i)
            enddo

            un(k) = un(k)/im
            vn(k) = vn(k)/im
            do i=1,imh
               uaglob(i) = -un(k)*sinlon(i) + vn(k)*coslon(i)
               vaglob(i) = -un(k)*coslon(i) - vn(k)*sinlon(i)
               uaglob(i+imh) = -uaglob(i)
               vaglob(i+imh) = -vaglob(i)
            enddo
            do i=ifirst,ilast
               ua(i,jm,k) = uaglob(i)
               va(i,jm,k) = vaglob(i)
            enddo
         enddo
      endif

      return
!EOC
      end
!-----------------------------------------------------------------------
