#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: highp2 --- Compute finite-volume mean pressure forces
!                     using more accurate contour decomposition
!
! !INTERFACE:
      subroutine highp2(pk,   wz,    pzx,  pzy,        &
                       pzz,  im,  jm,   km,            &
                       jfirst,  jlast, kfirst, klast,  &
                       klastp, nx )

! !USES:
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid, only: npr_y, myid_z, npr_z, iam, strip3zatyj2,   &
                        strip3zatypj1, strip3zatypj2
#if defined (SPMD)
      use spmd_dyn, only: comm_z
      use parutilitiesmodule, only: pargatherreal, parscatterreal
      use mod_comm, only: mp_send3d, mp_recv3d
#endif
      implicit none

! !INPUT PARAMETERS:

      integer im, jm, km, jfirst, jlast, kfirst, klast, klastp
      integer nx                        ! # of pieces in longitude direction

      real(r8) pk(im,jfirst-1:jlast+1,kfirst:klast+1)
      real(r8) wz(im,jfirst-1:jlast+1,kfirst:klast+1) 

! !OUTPUT PARAMETERS
      real(r8) pzx(im,jfirst-1:jlast  ,kfirst:klast+1)  
      real(r8) pzy(im,jfirst  :jlast+1,kfirst:klast+1) 
      real(r8) pzz(im,jfirst-1:jlast+1,kfirst:klast)

! ! !DESCRIPTION:
! Note: this routine needs to be further optimized for speed.
!
! !REVISION HISTORY:
!  AAM    01.06.27     Institute 2D decomposition
!  WS     02.01.15     Transitioned to mod_comm
!  WS     02.04.24     New mod_comm interfaces
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local:
      integer i, j, k
      integer ixj, jp, it, i1, i2

      real(r8) pka(km+1)
      real(r8) pkb(km+1)
      real(r8) wza(km+1)
      real(r8) wzb(km+1)
      real(r8) pf(km+1)
      real(r8) pza(km)
      real(r8) pzb(km)
      real(r8) dka(km)
      real(r8) dkb(km)

      real(r8)    tmp
      integer im1

      integer js1g1
      integer js2g0
!     integer js2g1
      integer jn2g0
      integer jn1g1

!
! Additional variables for yz decomposition
!
      real (r8), allocatable:: pkz(:,:,:), wzz(:,:,:), pzxz(:,:,:),  &
            pzyz(:,:,:), pzzz(:,:,:)
      real (r8), allocatable:: pky(:,:,:), wzy(:,:,:), pzxy(:,:,:),  &
            pzyy(:,:,:)

#if defined(SPMD)
      integer incount, outcount
#endif

      allocate(pkz(im,jfirst-1:jlast+1,km+1))
      allocate(wzz(im,jfirst-1:jlast+1,km+1))
      allocate(pzxz(im,jfirst-1:jlast,km+1))
      allocate(pzyz(im,jfirst:jlast+1,km+1))
      allocate(pzzz(im,jfirst-1:jlast+1,km))

#if defined(SPMD)
      allocate(pky(im,jfirst-1:jlast+1,kfirst:klastp))
      allocate(wzy(im,jfirst-1:jlast+1,kfirst:klastp))
      allocate(pzxy(im,jfirst-1:jlast,kfirst:klastp))
      allocate(pzyy(im,jfirst:jlast+1,kfirst:klastp))
#endif

!
! Gather pk,wz into full z arrays
!
#if defined(SPMD)
      if (myid_z .eq. 0) then
#endif

!$omp  parallel do     &
!$omp  default(shared) &
!$omp  private(i,j,k)

        do k = kfirst, klast+1
        do j = jfirst-1, jlast+1
        do i = 1, im
          pkz(i,j,k) = pk(i,j,k)
          wzz(i,j,k) = wz(i,j,k)
        enddo
        enddo
        enddo
#if defined(SPMD)
      endif
#endif
#if defined(SPMD)
      if (npr_z .gt. 1) then
        do k = kfirst, klastp
        do j = jfirst-1, jlast+1
        do i = 1, im
          pky(i,j,k) = pk(i,j,k)
          wzy(i,j,k) = wz(i,j,k)
        enddo
        enddo
        enddo
        call pargatherreal( comm_z, 0, pky, strip3zatypj2, pkz)
        call pargatherreal( comm_z, 0, wzy, strip3zatypj2, wzz)
      endif
#endif

      js1g1  = max(1,jfirst-1)

! E-W PG
!     js2g1  = max(2,jfirst-1)
      jn2g0  = min(jm-1,jlast)

! N-S PG
      js2g0  = max(2,jfirst)
      jn1g1  = min(jm,jlast+1)

      it = im / nx
      jp = nx * ( jn1g1 - js1g1  + 1 )

!$omp  parallel do           &
!$omp  default(shared)       &
!$omp  private(i1, i2, ixj, i, j, k, tmp, pka, pkb, pf) &
!$omp  private(pza, pzb, wza, wzb, dka, dkb, im1)

      do 1000 ixj=1, jp
! ***
! E-W
! ***
         j  = js1g1 + (ixj-1)/nx
         i1 = 1 + it * mod(ixj-1, nx)
         i2 = i1 + it - 1

        do i=i1,i2

        if( i .eq. i1 ) then

           if ( i .eq. 1) then
              im1 = im
           else
              im1 = i-1
           endif 

        do k=1,km+1
           pka(k) =  pkz(im1,j,k)
           wza(k) =  wzz(im1,j,k)
        enddo

        do k=1,km
           dka(k) =  pka(k+1) - pka(k)
           pza(k) = 0.5*(wza(k)-wza(k+1))*(pka(k+1)+pka(k))
        enddo

        else

        do k=1,km+1
           pka(k) = pkb(k)
           wza(k) = wzb(k)
        enddo

        do k=1,km
           pza(k) = pzb(k)
           dka(k) = dkb(k)
        enddo

        endif

        do k=1,km+1
           pkb(k) = pkz(i,j,k)
           wzb(k) = wzz(i,j,k)
        enddo

        do k=1,km
           dkb(k) = pkb(k+1) - pkb(k)
           tmp    = (wzb(k)-wzb(k+1))*(pkb(k+1)+pkb(k))
           pzb(k) = 0.5*tmp
           pzzz(i,j,k) = tmp
        enddo

        if ( j .le. jn2g0 ) then
        call pxt3 (km,   pf,   pka,  pkb,  wza,  wzb, &
                   pza,  pzb,  dka,  dkb )

        do k=1,km+1
           pzxz(i,j,k) = pf(k)
        enddo
        endif

! ***
! N-S
! ***

! N <-- E at (i,j)

        if ( j .ge. js2g0 ) then
        do k=1,km+1
           pka(k) = pkz(i,j-1,k)
           wza(k) = wzz(i,j-1,k)
        enddo

        do k=1,km
           dka(k) = pka(k+1) - pka(k)
           pza(k) = 0.5*(wza(k)-wza(k+1))*(pka(k+1)+pka(k)) 
        enddo

! need pzy(js2g0:jn2g1)
        call pxt3 (km,    pf,   pka,  pkb,  wza,  wzb,  &
                   pza,  pzb,   dka,  dkb )

        do k=1,km+1
           pzyz(i,j,k) = pf(k)
        enddo
        endif

        enddo                  ! i-loop
1000  continue

!
! Scatter pz[xyz]z into local arrays
!
#if defined(SPMD)
      if (myid_z .eq. 0) then
#endif

!$omp  parallel do     &
!$omp  default(shared) &
!$omp  private(i,j,k)

        do k = kfirst, klast+1
        do j = jfirst-1, jlast
        do i = 1, im
          pzx(i,j,k) = pzxz(i,j,k)
          pzy(i,j+1,k) = pzyz(i,j+1,k)
        enddo
        enddo
        enddo

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i,j,k)

        do k = kfirst, klast
        do j = jfirst-1, jlast+1
        do i = 1, im
          pzz(i,j,k) = pzzz(i,j,k)
        enddo
        enddo
        enddo
#if defined(SPMD)
      endif
#endif

#if defined(SPMD)
      if (npr_z .gt. 1) then
        call parscatterreal(comm_z, 0, pzxz, strip3zatypj1, pzxy)
        call parscatterreal(comm_z, 0, pzyz, strip3zatypj1, pzyy)
        call parscatterreal(comm_z, 0, pzzz, strip3zatyj2, pzz)

        if (myid_z .ne. 0) then

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i,j,k)
          do k = kfirst, klastp
          do j = jfirst-1, jlast
          do i = 1, im
            pzx(i,j,k) = pzxy(i,j,k)
            pzy(i,j+1,k) = pzyy(i,j+1,k)
          enddo
          enddo
          enddo
        endif

! Fill in klast+1
        call mp_send3d( iam-npr_y, iam+npr_y, im, jm, km+1,             &
                        1, im, jfirst-1, jlast, kfirst, klast+1,        &
                        1, im, jfirst-1, jlast, kfirst, kfirst, pzx )
        call mp_send3d( iam-npr_y, iam+npr_y, im, jm, km+1,             &
                        1, im, jfirst, jlast+1, kfirst, klast+1,        &
                        1, im, jfirst, jlast+1, kfirst, kfirst, pzy )
        call mp_recv3d( iam+npr_y, im, jm, km+1,                        &
                        1, im, jfirst-1, jlast, kfirst, klast+1,        &
                        1, im, jfirst-1, jlast, klast+1, klast+1, pzx )
        call mp_recv3d( iam+npr_y, im, jm, km+1,                        &
                        1, im, jfirst, jlast+1, kfirst, klast+1,        &
                        1, im, jfirst, jlast+1, klast+1, klast+1, pzy )

      endif
#endif

#if defined(SPMD)
      deallocate(pky)
      deallocate(wzy)
      deallocate(pzxy)
      deallocate(pzyy)
#endif

      deallocate(pkz)
      deallocate(wzz)
      deallocate(pzxz)
      deallocate(pzyz)
      deallocate(pzzz)

      return
!EOC
      end

      subroutine pxt3(km,  pf,  pka, pkb, wza, wzb, pza, pzb, &
                      dka, dkb )
 
! Compute pressure forcing along "horizontal" coordinate surface
! for a single column

! !USES:
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none

      integer km
 
! Input
      real(r8) pka(km+1)     ! p**kappa at the west side of the box
      real(r8) pkb(km+1)

      real(r8) wza(km+1)
      real(r8) wzb(km+1)

      real(r8) pza(km)
      real(r8) pzb(km)

      real(r8) dka(km)
      real(r8) dkb(km)

! Output
      real(r8)  pf(km+1)

! Local
      integer k
      integer ka
      integer kb
      integer kk
      integer kt
      real(r8) wa
      real(r8) wb
      real(r8) tmp
      logical found

      pf(1) = -(pka(1) + pkb(1))*(wzb(1) - wza(1))

         ka = km/2
         kb = km/2

      do 1000 k=2,km+1

          wa = 0.
          wb = 0.

      if( wzb(k) .gt. wza(k) ) then

!
!                     * ---->  wzb(k)
!                    /
!                   /
!                  /
!                 /
!                /
!               /
!              /
!             /
!            /
!           /
!          /
!         /
!        /
!       /
!      /
!     /
!    * ----> wza(k)
!
! Compute wa (along left edge)
!
! Top
      if(wzb(k) .gt. wza(1)) then 
         wa = 0.5*(pkb(k)+pka(k))*(wzb(k)-wza(k))
! DEBUG
!        write(*,*) 'pxt3: ', wa
!        stop
! DEBUG
      else

            kk = k
            found = .false.
         do while (.not. found )
            if(wzb(k) .le. wza(kk-1) .and.       &
               wzb(k) .ge. wza(kk)       )  then
!
! find p**cappa at left edge (side A) at height wzb(k)
!
               tmp = pka(kk-1) + dka(kk-1) *     &
                    (wza(kk-1)-wzb(k))/(wza(kk-1)-wza(kk))
               wa = wa + 0.5*(wzb(k)-wza(kk))*(pka(kk)+tmp)
               found = .true.
            else
               wa  = wa + pza(kk-1)
               kk = kk - 1
            endif
         enddo
      endif

! Compute wb

         if( k .ne. (km+1)) then
            do kk=k,km
               if( wza(k) .le. wzb(kk)   .and.    &
                   wza(k) .ge. wzb(kk+1)      ) then
                   tmp = pkb(kk) + dkb(kk)*(wzb(kk)-wza(k)) &
                                 / (wzb(kk)-wzb(kk+1))   
                   wb = wb + 0.5*(tmp+pkb(kk))*(wzb(kk)-wza(k))
                   goto 222
               else
                  wb = wb + pzb(kk)
               endif
            enddo
222      continue
         endif        ! k safety check

         if( wza(k) .lt. wzb(km+1)) then
! Integrate down along the left edge!!

            do kk=kb,k-1
            if( wzb(km+1) .le. wza(kk)   .and.    &
                wzb(km+1) .ge. wza(kk+1)      ) then
                tmp = pka(kk) + dka(kk)*(wza(kk)-wzb(km+1)) &
                              / (wza(kk)-wza(kk+1))
                wb = wb + 0.5*(pka(kk+1)+tmp)*(wzb(km+1)-wza(kk+1))
                if( kk .le. k-2) then
                    do kt=kk+1, k-1
                       wb = wb + pza(kt)
                    enddo
                endif
                kb = kk
                goto 444
            endif
            enddo
444      continue
         endif

      elseif( wzb(k) .lt. wza(k) ) then

!============================
! for wzb(k) .lt. wza(k) 
!============================

! Compute wa

         if(k .ne. km+1) then
            do kk=k,km
               if( wzb(k) .le. wza(kk)   .and.   &
                   wzb(k) .ge. wza(kk+1)      ) then
                   tmp = pka(kk) + dka(kk)*(wza(kk)-wzb(k)) &
                                         / (wza(kk)-wza(kk+1))
                   wa = wa - 0.5*(tmp+pka(kk))*(wza(kk)-wzb(k))
                  goto 666
               else
                  wa  = wa - pza(kk)
               endif
            enddo
666      continue
         endif        ! k safety check

         if( wzb(k) .lt. wza(km+1)) then
! Integrate  down along the right edge!!
            do kk=ka,k-1
              if( wza(km+1) .le. wzb(kk)   .and.    &
                  wza(km+1) .ge. wzb(kk+1)   ) then
                  tmp = pkb(kk) + dkb(kk)*(wzb(kk)-wza(km+1)) &
                                        / (wzb(kk)-wzb(kk+1))
                  wa = wa - 0.5*(tmp+pkb(kk+1))*(wza(km+1)-wzb(kk+1))
                  if(kk .le. k-2) then
                     do kt=kk+1, k-1
                        wa  = wa - pzb(kt)
                     enddo
                  endif
                  ka = kk
                 goto 777
              endif
            enddo
777      continue
         endif

! Top
       if(wza(k) .gt. wzb(1)) then 
          wb = 0.5*(pkb(k)+pka(k))*(wzb(k)-wza(k))
! DEBUG
!         write(*,*) 'pxt3: ', wb
!         stop
! DEBUG
       else

          do kk=k,2,-1
             if(wza(k) .le. wzb(kk-1) .and.    &
                wza(k) .ge. wzb(kk)       )  then
                tmp = pkb(kk-1)+dkb(kk-1)*(wzb(kk-1)-wza(k)) &
                                       / (wzb(kk-1)-wzb(kk))
                wb = wb + 0.5*(tmp+pkb(kk))*(wzb(kk)-wza(k))
                goto 999
             else
                wb  = wb - pzb(kk-1)
             endif
         enddo
999     continue
       endif   ! end top checck
      endif

      pf(k) = -(wa + wb) 

1000  continue

      return
      end
