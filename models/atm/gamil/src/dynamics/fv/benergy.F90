#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: benergy --- Calculate the total energy
!
! !INTERFACE:

      subroutine benergy(im, jm, km, u, v, pt, ng_d, ng_s, delp, pe, pk, pkz, hs, &
                        cp,  te0, te, dz, jfirst, jlast,              &
                        kfirst, klast, klastp)

! !USES:

      use shr_kind_mod, only: r8 => shr_kind_r8
      use dynamics_vars, only : cosp, acap
      use pmgrid, only : npr_z, myid_z, npr_y, myid_y, strip3zaty, strip3zatypt, iam
#if defined( SPMD )
      use spmd_dyn, only : comm_z, comm_y
      use mod_comm, only: mp_send3d, mp_recv3d
      use parutilitiesmodule, only : sumop, pargatherreal, parcollective
#endif
      implicit none

! !INPUT PARAMETERS:
      integer im, jm, km, jfirst, jlast                  ! Dimensions
      integer kfirst, klast, klastp
      integer ng_d                                       ! Ghosting zone of pt
      integer ng_s                                       ! max(ng_c+1,ng_d)
      real (r8)  u(im, jfirst-ng_d:jlast+ng_s, kfirst:klast)   ! Winds x
      real (r8)  v(im, jfirst-ng_s:jlast+ng_d, kfirst:klast)   ! Winds y
      real (r8) pt(im, jfirst-ng_d:jlast+ng_d, kfirst:klast)   ! Potential temperature
      real (r8) delp(im, jfirst:jlast, kfirst:klast)           ! Delta pressure
      real (r8) pkz(im, jfirst:jlast, kfirst:klast)
      real (r8) pe(im, kfirst:klast+1, jfirst:jlast)           ! Edge pressure
      real (r8) pk(im, jfirst:jlast, kfirst:klast+1)
      real (r8) hs(im,jfirst:jlast)
      real (r8) cp

! !INPUT work arrays:
      real (r8) te(im, jfirst:jlast, kfirst:klast)             ! Work array (cache perf.)
      real (r8) dz(im, jfirst:jlast, kfirst:klast)             ! Work array (cache perf.)

! !OUTPUT PARAMETERS:
      real (r8) te0              ! globally integrated total energy

! !DESCRIPTION:
!    Determines the globally integrated total energy, if jfirst == 1
!    and jlast == jnp, otherwise it calculates the total energy in
!    the latitude slice jfirst:jlast. 
!
! !REVISION HISTORY:
!
! SJL 99.04.13 : Delivered as release 0.9.8
! WS  99.05.18 : Added im, jm, km, te, dz as arguments
! WS  99.05.25 : Replaced IMR by IM, JMR by JM-1; removed fvcore.h
! WS  99.10.11 : Ghosted U, now fully limited to jfirst:jlast
! WS  99.11.23 : Pruned te, additional cleaning
! WS  00.05.14 : Renamed ghost indices as per Kevin's definitions
! WS  00.07.13 : Changed PILGRIM API
! WS  00.08.28 : Cosmetic changes
! AAM 00.08.28 : Added kfirst,klast
! WS  00.12.01 : Replaced MPI_ON with SPMD; hs now distributed
! AAM 01.06.15 : Changes for zero diff
! WS  01.12.10 : Ghosted PT
! WS  01.12.31 : Ghosted U,V
! WS  02.07.04 : Fixed 2D decomposition bug dest/src for mp_send3d
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local
      real (r8) u2(im,jfirst:jlast+1)
      real (r8) v2(im,jfirst:jlast)
      real (r8) tte(jfirst:jlast)
      real (r8) bte(im)
      real (r8) tteglob(jm)
      real (r8) te_sp
      real (r8) te_np
      real (r8) gztop(im,jfirst:jlast)

      integer i, j, k, js1g0, js2g0, jn1g0, jn1g1, jn2g0, ktot, jtot
      integer dest, src
      real (r8) xsum

! Geometric arrays

#if !defined( SPMD )
      integer comm_y
#endif
      real (r8), allocatable :: dzte_tmp(:,:,:)
      real (r8), allocatable :: pez(:,:,:), pey(:,:,:)


! WS 99.07.27 :  Set loop limits appropriately

      js1g0  = max(1,jfirst)
      js2g0  = max(2,jfirst)
      jn2g0  = min(jm-1,jlast)
      jn1g0  = min(jm,jlast)
      jn1g1  = min(jm,jlast+1)
      ktot = klast-kfirst+1
      jtot = jlast-jfirst+1
 
#if defined(SPMD)
!
! Send one latitude to south; ghost U one latitude at jlast+1
!
      dest = iam-1
      src  = iam+1
      if ( mod(iam,npr_y) == 0 ) dest = -1
      if ( mod(iam+1,npr_y) == 0 ) src = -1
      call mp_send3d( dest, src, im, jm, km,                      &
                      1, im, jfirst-ng_d, jlast+ng_s, kfirst, klast, &
                      1, im, jfirst, jfirst, kfirst, klast, u )
      call mp_recv3d( src, im, jm, km,                             &
                      1, im, jfirst-ng_d, jlast+ng_s, kfirst, klast, &
                      1, im, jlast+1, jlast+1, kfirst, klast, u )
#endif

!$omp  parallel do         &
!$omp  default(shared)     &
!$omp  private(i,j,k, u2, v2,te_sp,te_np)

      do k=kfirst,klast
        do j=js2g0,jn1g1
           do i=1,im
              u2(i,j) = u(i,j,k)**2
           enddo
        enddo

        do j=js2g0,jn2g0
           do i=1,im
              v2(i,j) = v(i,j,k)**2
           enddo
        enddo

        do j=js2g0,jn2g0
           do i=1,im-1
           te(i,j,k) = 0.25 * ( u2(i,j) + u2(i,j+1) +         &
                                v2(i,j) + v2(i+1,j)  )
           enddo
! i=im
           te(im,j,k) = 0.25 * ( u2(im,j) + u2(im,j+1) +      &
                                 v2(im,j) + v2(1,j)  )
        enddo


        do j=js2g0,jn2g0
           do i=1,im
              te(i,j,k) = delp(i,j,k) * ( te(i,j,k) +         &
                             cp*pt(i,j,k)*pkz(i,j,k)  ) 
           enddo
        enddo


! poles

        if ( jfirst .EQ. 1 ) then
          te_sp = 0.
          do i=1,im
            te_sp = te_sp + u2(i,  2) + v2(i,  2)
          enddo
          te_sp =   delp(1,  1,k) * (0.5*te_sp/float(im) +    &
                    cp*pt(1, 1,k)*pkz(1,1,k)      )
          do i=1,im
            te(i,  1,k) = te_sp
          enddo
        endif
  
        if ( jlast .EQ. jm ) then
          te_np = 0.
          do i=1,im
            te_np = te_np + u2(i,jm) + v2(i,jm-1)
          enddo
          te_np =   delp(1,jm,k) * (0.5*te_np/float(im) +     &
                    cp*pt(1,jm,k)*pkz(1,jm,k)   )
          do i=1,im
            te(i,jm,k) = te_np
          enddo
        endif

! Compute dz
        do j=js1g0,jn1g0
           do i=1,im
              dz(i,j,k) = cp*pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
           enddo
        enddo

      enddo

! Perform vertical integration

      allocate(dzte_tmp(im,jfirst:jlast,km))
      allocate(pez(im,km+1,jfirst:jlast))
      allocate(pey(im,kfirst:klastp,jfirst:jlast))

! Compute gztop

#if defined( SPMD )
      if (myid_z .eq. 0) then
#endif
!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i,j,k)
        do k = kfirst, klast
         do j = jfirst, jlast
          do i = 1, im
           dzte_tmp(i,j,k) = dz(i,j,k)
          enddo
         enddo
        enddo
#if defined( SPMD )
      endif
#endif
#if defined( SPMD )
      if (npr_z .gt. 1) then
        call pargatherreal(comm_z, 0, dz, strip3zaty, dzte_tmp)
      endif
#endif

#if defined( SPMD )
      if (myid_z .eq. 0) then
#endif

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i,j,k)
        do j = jfirst, jlast
          do i=1,im
            gztop(i,j) = hs(i,j)
          enddo
          do k=1,km
            do i=1,im
              gztop(i,j) = gztop(i,j) + dzte_tmp(i,j,k)
            enddo
          enddo
        enddo

#if defined( SPMD )
      endif
#endif

#if defined( SPMD )
      if (myid_z .eq. 0) then
#endif
!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i,j,k)
        do j = jfirst, jlast
         do k = kfirst, klastp
          do i = 1, im
           pez(i,k,j) = pe(i,k,j)
          enddo
         enddo
        enddo
#if defined( SPMD )
      endif
#endif
#if defined( SPMD )
      if (npr_z .gt. 1) then
!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i,j,k)
        do j = jfirst, jlast
         do k = kfirst, klastp
          do i = 1, im
           pey(i,k,j) = pe(i,k,j)
          enddo
         enddo
        enddo
        call pargatherreal(comm_z, 0, pey, strip3zatypt, pez)
      endif
#endif

#if defined( SPMD )
      if (myid_z .eq. 0) then
#endif
!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i,j,k)
        do k = kfirst, klast
         do j = jfirst, jlast
          do i = 1, im
           dzte_tmp(i,j,k) = te(i,j,k)
          enddo
         enddo
        enddo
#if defined( SPMD )
      endif
#endif
#if defined( SPMD )
      if (npr_z .gt. 1) then
        call pargatherreal(comm_z, 0, te, strip3zaty, dzte_tmp)
      endif
#endif

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(j)

      do j=js1g0, jn1g0
         tte(j) = 0.
      enddo

#if defined( SPMD )
      if (myid_z .eq. 0) then
#endif

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i,j,k,bte)
         do j=js1g0, jn1g0
            tte(j) = 0.
            if  (j .eq. 1) then
               tte(1) = tte(1) - pez(1,1,1)*gztop(1,1)
               tte(1) = tte(1) + pez(1,km+1,1)*hs(1,1)
               do k=1,km
                  tte(1) = tte(1) + dzte_tmp(1,1,k)
               enddo
               tte(1) = acap * tte(1)
            elseif (j .eq. jm) then
               tte(jm) = tte(jm) - pez(1,1,jm)*gztop(1,jm)
               tte(jm) = tte(jm) + pez(1,km+1,jm)*hs(1,jm)
               do k=1,km
                  tte(jm) = tte(jm) + dzte_tmp(1,jm,k)
               enddo
               tte(jm) = acap * tte(jm)
            else
               do i = 1,im
                  bte(i) = 0.
               enddo
               do i = 1,im
                  bte(i) = bte(i) - pez(i,1,j)*gztop(i,j)
                  bte(i) = bte(i) + pez(i,km+1,j)*hs(i,j)
                  do k=1,km
                     bte(i) = bte(i) + dzte_tmp(i,j,k)
                  enddo
               enddo
               do i = 1,im
                  tte(j) = tte(j) + bte(i)
               enddo
               tte(j) = tte(j) * cosp(j)
            endif
         enddo
#if defined( SPMD )
      endif
#endif

#if defined( SPMD )
      if (npr_z .gt. 1) then
        call parcollective(comm_z, sumop, jtot, tte)
      endif
#endif

      deallocate(dzte_tmp)
      deallocate(pez)
      deallocate(pey)

! Sum over j

      te0=0.
      do j = jfirst, jlast
        tteglob(j) = tte(j)
      enddo
      call par_vecsum(jm, jfirst, jlast, tteglob, te0, comm_y, npr_y)

      return
!EOC
      end
!-----------------------------------------------------------------------
