#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: hswf --- Held-Suarez forcing
!
! !INTERFACE:

      subroutine hswf(im,     jm,     km,     jfirst,   jlast,          &
                      ifirst, ilast,  u,      v,        pt,             &
                      ng_d,   ngus,   ngun,   ngvs,     ngvn,           &
                      pe,     pkz,    pdt,    akap,     grav,           &
                      rg,     dcaf,   strat,  rayf,     sinp,           &
                      cosp,   sine,   cose,   coslon,   sinlon)

! !USES:

      use shr_kind_mod, only: r8 => shr_kind_r8
      use dynamics_vars, only : sinp2, cosp2, cosp4, rf, ks
      use pmgrid, only: iam, myidxy_x, myidxy_y, nprxy_x, twod_decomp

#if defined( SPMD )
      use parutilitiesmodule, only : parcollective, sumop
      use spmd_dyn, only : commxy_x, commxy_y, comm_y
      use mod_comm, only : mp_send3d, mp_recv3d
#endif
      implicit none

! !INPUT PARAMETERS:
      integer im, jm, km
      integer jfirst, jlast
      integer ifirst, ilast
      integer pdt
      integer ng_d               ! number of ghosted latitudes on D grid
      integer ngus               ! number of ghosted latitudes U south
      integer ngun               ! number of ghosted latitudes U north
      integer ngvs               ! number of ghosted latitudes V south
      integer ngvn               ! number of ghosted latitudes V north
      real (r8) akap, grav, rg
      logical strat
      logical rayf
      logical dcaf
      real (r8) cosp(jm),sinp(jm),cose(jm),sine(jm)
      real (r8) coslon(im), sinlon(im)
      real (r8), intent(in):: pe(ifirst:ilast,km+1,jfirst:jlast)
      real (r8), intent(in):: pkz(ifirst:ilast,jfirst:jlast,km)

! !INPUT/OUTPUT PARAMETERS:
      real (r8), intent(inout):: u(ifirst:ilast,jfirst-ngus:jlast+ngun,km)
      real (r8), intent(inout):: v(ifirst:ilast,jfirst-ngvs:jlast+ngvn,km)
      real (r8), intent(inout):: pt(ifirst:ilast,jfirst-ng_d:jlast+ng_d,km)

! !DESCRIPTION:
!    Author: Shian-Jiann Lin, NASA/GSFC
!
! !REVISION HISTORY:
!   SJL 99.09.30:  Delivery
!   WS  99.10.28:  jfirst:jlast; documentation
!   WS  99.11.07:  pruned arrays
!   WS  00.07.10:  Incorporated simplfications to PILGRIM
!   WS  00.08.28:  Change to ParEndTransfer, SPMD instead of MPI_ON
!   WS  00.12.22:  Replaced d2a3d call with d2a3dijk
!   AAM 01.06.13:  2-D decomposition
!   WS  01.12.30:  Ghosted U,V
!   WS  02.01.15:  Transitioned to mod_comm
!   WS  02.04.25:  New mod_comm interfaces
!
!EOP
!-----------------------------------------------------------------------
!BOC
!

      real (r8) p0, t0, sday, rkv, rka, rks, rkt, sigb, rsgb
      real (r8) tmp
      real (r8) ap0k, algpk, rfc
      real (r8) tey, tez, fac, pw, sigl, tmin
      real (r8) pl(ifirst:ilast,km)
      real (r8) pl1(km,jfirst:jlast)
      real (r8) frac(ifirst-1:ilast,jfirst-1:jlast,km)
      real (r8) teq(ifirst:ilast,km)
      real (r8) h0, dz
      real (r8) dt_tropic
      real (r8) rmr, rms
      real (r8) relx, tau
      real (r8) t_st, t_ms
      real (r8) f1
      real (r8) pc, c1
      real (r8) t2(ifirst:ilast, km)
      real (r8) dp(ifirst:ilast, km)

      real (r8) ua(ifirst:ilast,jfirst:jlast,km)
      real (r8) va(ifirst:ilast,jfirst:jlast,km)
      real (r8) u2(ifirst:ilast,km)
      real (r8) v2(ifirst:ilast,km)
      real (r8) fu(ifirst:ilast,km)
      real (r8) fv(ifirst:ilast,km,jfirst:jlast)
      real (r8) fvwest(km,jfirst:jlast)

#if defined( SPMD )
      real (r8) , allocatable :: uasouth(:,:), pesouth(:,:)
      integer dest, src, incount, outcount
#endif

      integer i, j, k, js2g0, js2gm1, jn2g0, count, ierror, itot, jtot

      real (r8) p_drycon
      parameter (p_drycon  = 100000.)

      logical first
      data first /.true./

      js2g0  = max(2,jfirst)
      js2gm1 = max(2,jfirst+1)
      jn2g0  = min(jm-1,jlast)
      itot = ilast - ifirst + 1
      jtot = jlast - jfirst + 1

      p0 = 100000.
      t0 = 200.
      h0 = 7.
      sday = 24*3600
      rkv = 0.5*pdt/sday
      rka = pdt/ (40.*sday)      ! was 40 days
      rfc = 1./(1.+rka)
      rks = pdt/ (4.*sday)       ! was 4 days

! For strat-mesosphere
      t_ms = 10.
      t_st = 40.
      tau = (t_st - t_ms) / log(100.)
      rms = pdt/(t_ms*sday)
      rmr =  1./(1.+rms)

      sigb = 0.7
      rsgb = 1./(1.-sigb)
      ap0k = 1./p0**akap
      algpk = log(ap0k)

#if defined (SPMD)
      allocate( uasouth( ifirst:ilast,km ) )
      allocate( pesouth( ifirst:ilast,km+1 ) )
#endif

      if (dcaf) then
!
! WS 99.10.26 : Replaced loop 500 with one call to d2a3d
! WS 00.12.22 : Replaced d2a3d with d2a3dijk (d2a2 inlined)
!
        call d2a3dijk(u,v,ua,va,im,jm,km,jfirst,jlast,ng_d, &
                      ngus,ngun,ngvs,ngvn,ifirst,ilast,     &
                      coslon,sinlon)
      endif

      if (myidxy_x .eq. 0) then
         do j = jfirst,jlast
            do k = 1,km
               pl1(k,j) = 0.5*(pe(1,k,j)+pe(1,k+1,j))
            enddo
         enddo
      else
         do j=jfirst,jlast
            do k=1,km
               pl1(k,j) = 0.
            enddo
         enddo
      endif
#if defined (SPMD)
      if (itot .ne. im) then
         call parcollective(commxy_x, sumop, km, jtot, pl1)
      endif
#endif

!$omp  parallel do     &
!$omp  default(shared) &
!$omp  private(i,j,k,pl,tey,tez,dz,relx,dt_tropic)   &
!$omp  private(teq,tmin,sigl,f1,rkt, tmp, u2, v2, fu, t2, dp)

      do 1000 j=jfirst,jlast

        tey = ap0k*(315.-60.*sinp2(j))
        tez = ap0k*10./akap*cosp2(j)

        do k=1,km
          do i=ifirst,ilast
            pl(i,k) = 0.5*(pe(i,k,j)+pe(i,k+1,j))
          enddo
        enddo

        do k=km,1,-1
          do i=ifirst,ilast
            if (strat .and. pl(i,k) .lt. 10000.  &
                      .and. pl(i,k) .gt. 100.  )  then
              dz = h0 * log(pl(i,k+1)/pl(i,k))
!
! Lapse rate above tropic stratopause is 2.25 deg/km
! Relaxation time is t_st days at 100 mb (as H-S) and gradually
! decreases to t_ms Days at and above the stratopause
!
              relx =  t_ms + tau*log(0.01*pl(i,k))
              relx = pdt/(relx*sday)
              dt_tropic = 2.25*cosp(j) * dz
              teq(i,k) = (teq(i,k+1)*pkz(i,j,k+1) +  &
                          dt_tropic)/pkz(i,j,k)
              pt(i,j,k) = (pt(i,j,k)+relx*teq(i,k))/(1.+relx)
!!!              pt(i,j,k) = teq(i,k)
            elseif (strat .and. pl(i,k) .le. 100.)  then
!
! Mesosphere
!
              dz = h0 * log(pl(i,k+1)/pl(i,k))
              dt_tropic = -2.25*cosp(j) * dz
              tmp = teq(i,k+1)*pkz(i,j,k+1) + dt_tropic
              teq(i,k) =  tmp / pkz(i,j,k)
!!!              teq(i,k) = max(200., tmp) / pkz(i,j,k)
              pt(i,j,k) = (pt(i,j,k)+rms*teq(i,k))*rmr
!!!              pt(i,j,k) = teq(i,k)
            else
!
! Trop:  strictly Held-Suarez
!
              sigl = pl(i,k)/pe(i,km+1,j)
              f1 = max(0., (sigl-sigb) * rsgb )
              tmin = t0/pkz(i,j,k)
              teq(i,k) = tey - tez*(log(pkz(i,j,k))+algpk)
              teq(i,k) = max(tmin, teq(i,k))
              rkt = rka + (rks-rka)*f1*cosp4(j)
              pt(i,j,k) = (pt(i,j,k)+rkt*teq(i,k))/(1.+rkt)
            endif
          enddo     !i-loop
        enddo     !k-loop

! Do dry_con
        if (dcaf) then
          do k=1, km
            do i=ifirst,ilast
              dp(i,k) = pe(i,k+1,j) - pe(i,k,j)
              fu(i,k) = 0.
              fv(i,k,j) = 0.
              u2(i,k) = ua(i,j,k)
              v2(i,k) = va(i,j,k)
              t2(i,k) = pt(i,j,k)
            enddo
          enddo

          call dry_adj(itot, km, 1., p_drycon, pl1(1,j), fu, fv(ifirst,1,j), t2,  &
                       u2, v2, dp )

          do k=1, km
            do i=ifirst,ilast
              ua(i,j,k) = fu(i,k)
              pt(i,j,k) = t2(i,k)
            enddo
          enddo
        endif

1000  continue

! Adjust D-grid v-winds

      if (dcaf) then

!$omp  parallel do     &
!$omp  default(shared) &
!$omp  private(j,k)

         do k=1,km
            do j=jfirst,jlast
               fvwest(k,j) = fv(ilast,k,j)
            enddo
         enddo

#if defined( SPMD )
! communicate to get fvwest (AAM)
        if (itot .ne. im) then
           dest = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
           src  = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)

           call mp_send3d( dest, src, im, km, jm,                     &
                           ifirst, ilast, 1, km, jfirst, jlast,       &
                           ilast, ilast, 1, km, jfirst, jlast, fv)
           call mp_recv3d( src, im, km, jm,                           &
                           ifirst-1, ifirst-1, 1, km, jfirst, jlast,  &
                           ifirst-1, ifirst-1, 1, km, jfirst, jlast, fvwest )

        endif
#endif

!$omp  parallel do     &
!$omp  default(shared) &
!$omp  private(i,j,k)

        do k=1, km
           do j=jfirst, jlast
              v(ifirst,j,k) = v(ifirst,j,k) + 0.5*(fv(ifirst,k,j)+fvwest(k,j))
              do i=ifirst+1, ilast
                 v(i,j,k) = v(i,j,k) + 0.5*(fv(i,k,j)+fv(i-1,k,j))
              enddo
           enddo
        enddo
      endif

#if defined( SPMD )
!
! Communication might include ua and/or pe on the south only
!
      call mp_send3d( iam+nprxy_x, iam-nprxy_x, im, jm, km,            &
                      ifirst, ilast, jfirst, jlast, 1, km,             &
                      ifirst, ilast, jlast, jlast, 1, km, ua )
      call mp_send3d( iam+nprxy_x, iam-nprxy_x, im, km+1, jm,          &
                      ifirst, ilast, 1, km+1, jfirst, jlast,           &
                      ifirst, ilast, 1, km+1, jlast, jlast, pe )
      call mp_recv3d( iam-nprxy_x, im, jm, km,                         &
                      ifirst, ilast, jfirst-1, jfirst-1, 1, km,        &
                      ifirst, ilast, jfirst-1, jfirst-1, 1, km, uasouth ) 
      call mp_recv3d( iam-nprxy_x, im, km+1, jm,                       &
                      ifirst, ilast, 1, km+1, jfirst-1, jfirst-1,      &
                      ifirst, ilast, 1, km+1, jfirst-1, jfirst-1, pesouth )
!
! Communication finished
!
#endif

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i,j,k,sigl,fac)

      do 2000 k=1,km

        if (dcaf) then
          do j=js2gm1,jlast
            do i=ifirst,ilast
              u(i,j,k) = u(i,j,k) + 0.5*(ua(i,j,k)+ua(i,j-1,k))
            enddo
          enddo
#if defined( SPMD )
          if ( jfirst .gt. 1 ) then
            do i=ifirst,ilast
              u(i,jfirst,k) = u(i,jfirst,k)     &
                              + 0.5*(ua(i,jfirst,k)+uasouth(i,k))
            enddo
          endif
#endif
        endif

        if (rayf .and. k.le. ks) then
! Apply Rayleigh friction
          do j=js2g0,jlast
            do i=ifirst,ilast
              u(i,j,k) = u(i,j,k)*rf(k)
            enddo
          enddo

          do j=js2g0,jn2g0
            do i=ifirst,ilast
              v(i,j,k) = v(i,j,k)*rf(k)
            enddo
          enddo
        else
! Surface Rayleigh friction according to Held-Suarez
          do j=jfirst,jlast
            do i=ifirst,ilast
              sigl = 0.5*(pe(i,k,j)+pe(i,k+1,j)) / pe(i,km+1,j)
              frac(i,j,k) = max(0., (sigl-sigb)*rsgb )
            enddo
          enddo
#if defined( SPMD )
          if ( jfirst .gt. 1 ) then
            do i=ifirst,ilast
              sigl = 0.5*(pesouth(i,k)+pesouth(i,k+1)) / pesouth(i,km+1)
              frac(i,jfirst-1,k) = max(0., (sigl-sigb)*rsgb )
            enddo
          endif
#endif

! Backward adjustment
          do j=js2g0,jlast
            do i=ifirst,ilast
              fac = frac(i,j,k)+frac(i,j-1,k)
              if (fac .gt. 0.) then
                u(i,j,k) = u(i,j,k)/(1.+rkv*fac)
              endif
            enddo
          enddo
        endif
2000  continue

!$omp  parallel do         &
!$omp  default(shared)     &
!$omp  private(j,k)

      do k=1,km
         do j=jfirst,jlast
            frac(ifirst-1,j,k) = frac(ilast,j,k)
         enddo
      enddo

#if defined( SPMD )
! communicate to get frac(ifirst-1,j,k) (AAM)
      if (itot .ne. im) then
         dest = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
         src  = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
         call mp_send3d( dest, src, im, jm, km,                         &
                         ifirst-1, ilast, jfirst-1, jlast, 1, km,       &
                         ilast, ilast, jfirst, jlast, 1, km, frac )
         call mp_recv3d( src, im, jm, km,                               &
                         ifirst-1, ilast, jfirst-1, jlast, 1, km,       &
                         ifirst-1, ifirst-1, jfirst, jlast, 1, km, frac )
      endif
#endif

!$omp  parallel do         &
!$omp  default(shared)     &
!$omp  private(i,j,k,fac)

      do k = 1, km
         if (.not. rayf .or. k .gt. ks) then
            do j=js2g0,jn2g0
               do i=ifirst,ilast
                  fac = frac(i,j,k)+frac(i-1,j,k)
                  if (fac .gt. 0.) then
                     v(i,j,k) = v(i,j,k)/(1.+rkv*fac)
                  endif
               enddo
            enddo
         endif
      enddo

#if defined( SPMD )
      deallocate( pesouth )
      deallocate( uasouth )
#endif

      return
!EOP
      end
!-----------------------------------------------------------------------

