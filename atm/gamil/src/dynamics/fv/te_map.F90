#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: te_map --- Map vertical Lagrangian coordinates to normal grid
!
! !INTERFACE:

   subroutine te_map(consv,   convt,  ps,      omga,    pe,         &
                     delp,    pkz,    pk,      mdt,     im,         &
                     jm,      km,     nx,      jfirst,  jlast,      &
                     ng_d,    ngus,   ngun,    ngvs,    ngvn,       &
                     ifirst,  ilast,  nq,                           &
                     u,       v,      pt,      q,       hs,         &
                     cp,      akap,   kord,    peln,    te0,        &
                     te,      dz,     nc )
!
! !USES:

   use shr_kind_mod, only: r8 => shr_kind_r8
   use dynamics_vars, only : acap, ks, ptop, cosp, ak, bk
   use mapz_module, only : map1_ppm, mapn_ppm
!  use fill_module
   use pmgrid, only: npr_y, nprxy_y, nprxy_x, myid_y, myidxy_y,      &
                     myidxy_x, twod_decomp, iam

#if defined( SPMD )
   use parutilitiesmodule, only: sumop,  parcollective
   use spmd_dyn, only: comm_y, commxy_x, commxy_y
   use mod_comm, only: mp_send3d, mp_recv3d
#endif

   implicit none

#if defined( SPMD )
#define CPP_PRT_PREFIX  if(iam==0)
#else
#define CPP_PRT_PREFIX
#endif

! !INPUT PARAMETERS:
   logical consv                 ! flag to force TE conservation
   logical convt                 ! flag to control pt output (see below)
   integer mdt                   ! mapping time step (same as phys)
   integer im, jm, km            ! x, y, z dimensions
   integer nq                    ! number of tracers (including h2o)
   integer nc                    ! 
   integer ng_d                  ! number of ghost latitudes on D grid
   integer ngus                  ! number of ghost latitudes for U south
   integer ngun                  ! number of ghost latitudes for U north
   integer ngvs                  ! number of ghost latitudes for V south
   integer ngvn                  ! number of ghost latitudes for V north
   integer nx                    ! number of SMP "decomposition" in x
   integer jfirst, jlast         ! starting & ending latitude index
   integer ifirst, ilast         ! starting & ending longitude index
   real(r8) hs(ifirst:ilast,jfirst:jlast) ! surface geopotential
   real(r8) cp
   real(r8) te0
 
! !INPUT/OUTPUT PARAMETERS:
   real(r8) pk(ifirst:ilast,jfirst:jlast,km+1) ! pe to the kappa
   real(r8) u(ifirst:ilast,jfirst-ngus:jlast+ngun,km)    ! u-wind (m/s)
   real(r8) v(ifirst:ilast,jfirst-ngvs:jlast+ngvn,km)    ! v-wind (m/s)
   real(r8) q(ifirst:ilast,jfirst-ng_d:jlast+ng_d,km,nc)! tracers including specific humidity
   real(r8) pe(ifirst:ilast,km+1,jfirst:jlast) ! pressure at layer edges
   real(r8) ps(ifirst:ilast,jfirst:jlast)      ! surface pressure
   real(r8) pt(ifirst:ilast,jfirst-ng_d:jlast+ng_d,km)   ! virtual potential temperature as input
                                     ! Output: virtual temperature if convt is true
                                     ! false: output is (virtual) potential temperature 
   real(r8)  te(ifirst:ilast,jfirst:jlast,km)  ! Work array (cache performance)
   real(r8)  dz(ifirst:ilast,jfirst:jlast,km)     ! Work array (cache performance)
 
! !OUTPUT PARAMETERS:
   real(r8) delp(ifirst:ilast,jfirst:jlast,km) ! pressure thickness
   real(r8) omga(ifirst:ilast,km,jfirst:jlast)    ! vertical press. velocity (pascal/sec)
   real(r8) peln(ifirst:ilast,km+1,jfirst:jlast)  ! log(pe)
   real(r8) pkz(ifirst:ilast,jfirst:jlast,km)     ! layer-mean pk for converting t to pt

! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! WS 99.05.19 : Replaced IMR, JMR, JNP and NL with IM, JM-1, JM and KM
! WS 99.05.25 : Revised conversions with IMR and JMR; removed fvcore.h
! WS 99.07.29 : Reduced computation region to jfirst:jlast
! WS 99.07.30 : Tested concept with message concatenation of te_map calls
! WS 99.10.01 : Documentation; indentation; cleaning
! SJL 99.12.31: SMP "decomposition" in-E-W direction
! WS 00.05.14 : Renamed ghost indices as per Kevin's definitions
! WS 00.07.13 : Changed PILGRIM API
! AM 00.08.29 : Variables in this routine will ultimately be decomposed in (i,j).
! AM 01.06.13 : 2-D decomposition; reordering summation causes roundoff difference.
! WS 01.06.10 : Removed "if(first)" section in favor of a variable module
! AM 01.06.27 : Merged yz decomposition technology into ccm code.
! WS 02.01.14 : Upgraded to mod_comm
! WS 02.04.22 : Use mapz_module from FVGCM
! WS 02.04.25 : New mod_comm interfaces
!
!EOP
!-----------------------------------------------------------------------
!BOC
! Local arrays:
      real(r8) rmin(nx*jm), rmax(nx*jm)
      real(r8) tte(jm)
! x-y
      real(r8)  u2(ifirst:ilast,jfirst:jlast+1)
      real(r8)  v2(ifirst:ilast+1,jfirst:jlast)
      real(r8)  t2(ifirst:ilast,jfirst:jlast)
      real(r8)  veast(jfirst:jlast,km)
      real(r8)  pewest(km+1,jfirst:jlast)
! y-z
      real(r8)  pe0(ifirst:ilast,km+1)
      real(r8)  pe1(ifirst:ilast,km+1)
      real(r8)  pe2(ifirst:ilast,km+1)
      real(r8)  pe3(ifirst:ilast,km+1)
      real(r8) phis(ifirst:ilast,km+1)
      real(r8) u2_sp(ifirst:ilast,km)
      real(r8) v2_sp(ifirst:ilast,km)
      real(r8) t2_sp(ifirst:ilast,km)
      real(r8) u2_np(ifirst:ilast,km)
      real(r8) v2_np(ifirst:ilast,km)
      real(r8) t2_np(ifirst:ilast,km)

! x
      real(r8)     gz(ifirst:ilast)
      real(r8)  ratio(ifirst:ilast)
      real(r8)    bte(ifirst:ilast)
! z
      real(r8) pe1w(km+1)
      real(r8) pe2w(km+1)

      integer i1w, nxu
      integer i, j, k, ic, js2g0, jn2g0, jn1g1
      integer kord
      integer krd

      real(r8) akap, dak, bkh, qmax, qmin
      real(r8) te_sp(km), te_np(km)
      real(r8) xysum(jfirst:jlast,2)
      real(r8) tmpik(ifirst:ilast,km)
      real(r8) tmpij(ifirst:ilast,jfirst:jlast,2)
      real(r8) dtmp
      real(r8) sum
      real(r8) rdt5
      real(r8) rg
      real(r8) te1
      real(r8) dlnp
      real(r8) tvm

      integer ixj, jp, it, i1, i2

#if defined( SPMD )
      integer :: dest, src
      real(r8), allocatable, save :: pesouth(:,:)
#endif

      integer comm_use, npry_use, itot

      logical diag
      logical first

      data diag    /.false./
      data first   /.true./

! WS 99.07.27 :  Set loop limits appropriately
      js2g0  = max(2,jfirst)
      jn1g1  = min(jm,jlast+1)
      jn2g0  = min(jm-1,jlast)
      do j=jfirst,jlast
         xysum(j,1) = 0.
         xysum(j,2) = 0.
      enddo
      do j=jfirst,jlast
         do i=ifirst,ilast
            tmpij(i,j,1) = 0.
            tmpij(i,j,2) = 0.
         enddo
      enddo
      do k=1,km
         do i=ifirst,ilast
            tmpik(i,k) = 0.
         enddo
      enddo

      itot = ilast-ifirst+1
      nxu = 1
      if (itot == im) nxu = nx

      if (first) then
#if defined( SPMD )
        if (.not. allocated(pesouth)) then
          allocate( pesouth(ifirst:ilast,km+1) )
        endif
#endif
        first = .false.
      endif

#if defined( SPMD )
      if (twod_decomp == 1) then
         comm_use = commxy_y
         npry_use = nprxy_y
      else
         comm_use = comm_y
         npry_use = npr_y
      endif

      call mp_send3d( iam-nprxy_x, iam+nprxy_x, im, jm, km,              &
                      ifirst, ilast, jfirst-ngus, jlast+ngun, 1, km,     &
                      ifirst, ilast, jfirst, jfirst, 1, km, u )
! Nontrivial x decomposition
      if (itot /= im) then
        dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
        src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
        call mp_send3d( dest, src, im, jm, km,                       &
                        ifirst, ilast, jfirst-ngvs, jlast+ngvn, 1,km,&
                        ifirst, ifirst, jfirst, jlast, 1, km, v )
      endif
#endif
      call pkez(nxu, im, km, jfirst, jlast, 1, km, ifirst, ilast,      &
                pe, pk, akap, ks, peln, pkz, .false.)
 
! Single subdomain case (periodic)
      do k=1,km
         do j=jfirst,jlast
            veast(j,k) = v(ifirst,j,k)
         enddo
      enddo
#if defined( SPMD ) 
      call mp_recv3d( iam+nprxy_x, im, jm, km,                           &
                      ifirst, ilast, jfirst-ngus, jlast+ngun, 1, km,     &
                      ifirst, ilast, jlast+1, jlast+1, 1, km, u )
! Nontrivial x decomposition
      if (itot /= im) then
        call mp_recv3d( src, im, jm, km,                                 &
                        ilast+1, ilast+1, jfirst, jlast, 1, km,          &
                        ilast+1, ilast+1, jfirst, jlast, 1, km, veast )
        dest = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
        src  = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
        call mp_send3d( dest, src, im, km+1, jm,                    &
                        ifirst, ilast, 1, km+1, jfirst, jlast,      &
                        ilast, ilast, 1, km+1, jfirst, jlast, pe )
      endif
      call mp_send3d( iam+nprxy_x, iam-nprxy_x, im, km+1, jm,            &
                      ifirst, ilast, 1, km+1, jfirst, jlast,             &
                      ifirst, ilast, 1, km+1, jlast, jlast, pe )
#endif

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(i,j, k, u2, v2, t2)

! Compute cp*T + KE

      do 1000 k=1,km

        do j=js2g0,jn1g1
           do i=ifirst,ilast
              u2(i,j) = u(i,j,k)**2
           enddo
        enddo

        do j=js2g0,jn2g0
          do i=ifirst,ilast
             v2(i,j) = v(i,j,k)**2
          enddo
          v2(ilast+1,j) = veast(j,k)**2
        enddo

        do j=jfirst,jlast
           do i=ifirst,ilast
              t2(i,j) = cp*pt(i,j,k)
           enddo
        enddo

        do j=js2g0,jn2g0
          do i=ifirst,ilast
            te(i,j,k) = 0.25 * ( u2(i,j) + u2(i,j+1) +        &
                                 v2(i,j) + v2(i+1,j) ) +      &
                        t2(i,j)*pkz(i,j,k)
          enddo
        enddo

! WS 99.07.29 : Restructuring creates small round-off.  Not clear why...

! Do collective Mpisum (in i) for te_sp and te_np below (AAM)
!
        if ( jfirst == 1 ) then
! South pole
          do i=ifirst,ilast
             u2_sp(i,k) = u2(i,2)
             v2_sp(i,k) = v2(i,2)
             t2_sp(i,k) = t2(i,1)
          enddo
        endif

        if ( jlast == jm ) then
! North pole
          do i=ifirst,ilast
             u2_np(i,k) = u2(i,jm)
             v2_np(i,k) = v2(i,jm-1)
             t2_np(i,k) = t2(i,jm)
          enddo
        endif

! Compute dz; geo-potential increments
        do j=jfirst,jlast
           do i=ifirst,ilast
              dz(i,j,k) = t2(i,j)*(pk(i,j,k+1)-pk(i,j,k))
           enddo
        enddo
1000  continue

#if defined( SPMD )
      if (itot /= im) then
        call mp_recv3d( src, im, km+1, jm,                          &
                        ifirst-1, ifirst-1, 1, km+1, jfirst, jlast, &
                        ifirst-1, ifirst-1, 1, km+1, jfirst, jlast, pewest )   
      endif
      call mp_recv3d( iam-nprxy_x, im, km+1, jm,                         &
                      ifirst, ilast, 1, km+1, jfirst-1, jfirst-1,        &
                      ifirst, ilast, 1, km+1, jfirst-1, jfirst-1, pesouth )
#endif

      if ( jfirst == 1 ) then

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(i, k)

         do k = 1, km
            te_sp(k) = 0.
            do i=ifirst,ilast
              tmpik(i,k) = u2_sp(i,k) + v2_sp(i,k)
              te_sp(k) = te_sp(k) + tmpik(i,k)
            enddo
         enddo

#if defined( SPMD )
         if (nprxy_x > 1) then
            call par_xsum(tmpik, ifirst, ilast, im, km, te_sp)
         endif
#endif

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(i, k)

         do k = 1, km
            te_sp(k) = 0.5*te_sp(k)/float(im) + t2_sp(ifirst,k)*pkz(ifirst,1,k)
            do i=ifirst,ilast
              te(i,  1,k) = te_sp(k)
            enddo
         enddo
      endif

      if ( jlast == jm ) then

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i, k)

         do k = 1, km
            te_np(k) = 0.
            do i=ifirst,ilast
              tmpik(i,k) = u2_np(i,k) + v2_np(i,k)
              te_np(k) = te_np(k) + tmpik(i,k)
            enddo
         enddo

#if defined( SPMD )
         if (nprxy_x > 1) then
            call par_xsum(tmpik, ifirst, ilast, im, km, te_np)
         endif
#endif

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i, k)

         do k = 1, km
            te_np(k) = 0.5*te_np(k)/float(im) + t2_np(ifirst,k)*pkz(ifirst,jm,k)
            do i=ifirst,ilast
              te(i,jm,k) = te_np(k)
            enddo
         enddo
      endif

      it = itot / nxu
      jp = nxu * ( jlast - jfirst + 1 )

!$omp  parallel do           &
!$omp  default(shared)       &
!$omp  private(i,j,k,ic,i1w,pe0,pe1,pe2,pe3,ratio)   &
!$omp  private(dak,bkh,rdt5,phis,krd, ixj,i1,i2) &
!$omp  private(pe1w, pe2w )

!     do 2000 j=jfirst,jlast
      do 2000 ixj=1,jp

        j  = jfirst + (ixj-1) / nxu
        i1 = ifirst + it * mod(ixj-1, nxu)
        i2 = i1 + it - 1

! Copy data to local 2D arrays.
        i1w = i1-1
        if (i1 == 1) i1w = im
        do k=1,km+1
           do i=i1,i2
              pe1(i,k) = pe(i,k,j)
           enddo
           if( itot == im ) then
               pe1w(k) = pe(i1w,k,j)
           else
               pe1w(k) = pewest(k,j)
           endif
        enddo

        do k=1,ks+1
           do i=i1,i2
              pe0(i,k) = ak(k)
              pe2(i,k) = ak(k)
              pe3(i,k) = ak(k)
            enddo
        enddo

        do k=ks+2,km
           do i=i1,i2
              pe0(i,k) = ak(k) + bk(k)* ps(i,j)
              pe2(i,k) = ak(k) + bk(k)*pe1(i,km+1)
           enddo
        enddo

        do i=i1,i2
           pe0(i,km+1) =  ps(i,j)
           pe2(i,km+1) = pe1(i,km+1)
        enddo

! Ghosting for v mapping
        do k=ks+2,km
           pe2w(k) = ak(k) + bk(k)*pe1w(km+1)
        enddo
        pe2w(km+1) = pe1w(km+1)

! Compute omga (dp/dt)
        rdt5 = 0.5 / float(mdt)
        do k=2,km+1
           do i=i1,i2
              pe0(i,k) = pe1(i,k) - pe0(i,k)
           enddo
        enddo

        do i=i1,i2
! update ps
          ps(i,j)   = pe1(i,km+1)
          omga(i,1,j) = rdt5 * pe0(i,2)
        enddo

        do k=2,km
          do i=i1,i2
             omga(i,k,j) = rdt5 * ( pe0(i,k) + pe0(i,k+1) )
          enddo
        enddo

        if(ks /= 0) then
           do k=1,ks
             dak = ak(k+1) - ak(k)
             do i=i1,i2
                delp(i,j,k) = dak
             enddo
           enddo
        endif

        do k=ks+1,km
          do i=i1,i2
             delp(i,j,k) = pe2(i,k+1) - pe2(i,k)
          enddo
        enddo

! Compute correction terms to Total Energy
        do i=i1,i2
           phis(i,km+1) = hs(i,j)      
        enddo

        do k=km,1,-1
          do i=i1,i2
             phis(i,k) = phis(i,k+1) + dz(i,j,k)   
          enddo
        enddo

        do k=1,km+1
          do i=i1,i2
             phis(i,k) = phis(i,k) * pe1(i,k)
          enddo
        enddo

! <<< Compute Total Energy >>>
        do k=1,km
          do i=i1,i2
            te(i,j,k) =  te(i,j,k) + (phis(i,k+1) - phis(i,k)) /   &
                         (pe1(i,k+1) - pe1(i,k) )
          enddo
        enddo

! Map Total Energy
        call map1_ppm ( km,   pe1,   te,                           &
                        km,   pe2,   te,  0,  0,                   &
                        itot, i1-ifirst+1, i2-ifirst+1,            &
                        j, jfirst, jlast, 1, kord)

! Map constituents

       if( nq /= 0 ) then
          if(kord == 8) then
             krd = 8
          else
             krd = 7
          endif

          call mapn_ppm ( km,   pe1,   q, nq,                      &
                          km,   pe2,   q, ng_d, ng_d,              &
                          itot, i1-ifirst+1, i2-ifirst+1,          &
                          j, jfirst, jlast, 0, krd)
       endif

! map u
        if(j /= 1) then

! WS 99.07.29 : protect j==jfirst case
          if (j > jfirst) then
            do k=2,km+1
              do i=i1,i2
                pe0(i,k) = 0.5*(pe1(i,k)+pe(i,k,j-1))
              enddo
            enddo

            do k=ks+2,km+1
              bkh = 0.5*bk(k)
              do i=i1,i2
                pe3(i,k) = ak(k) + bkh*(pe1(i,km+1)+pe(i,km+1,j-1))
              enddo
            enddo

#if defined( SPMD )
          else
!  WS 99.10.01 : Read in pe(:,:,jfirst-1) from the pesouth buffer
            do k=2,km+1
              do i=i1,i2
                pe0(i,k) = 0.5*(pe1(i,k)+pesouth(i,k))
              enddo
            enddo

            do k=ks+2,km+1
              bkh = 0.5*bk(k)
              do i=i1,i2
                pe3(i,k) = ak(k) + bkh*(pe1(i,km+1)+pesouth(i,km+1))
              enddo
            enddo
#endif
          endif

          call map1_ppm ( km,   pe0,    u,                                 &
                          km,   pe3,    u,                                 &
                          ngus, ngun, itot, i1-ifirst+1, i2-ifirst+1,      &
                          j,    jfirst, jlast,  -1,    kord)

        endif

! map v
        if(j /= 1 .and. j /= jm) then
          do k=2,km+1
! pe1(i1-1,1:km+1) must be ghosted
            pe0(i1,k) = 0.5*(pe1(i1,k)+pe1w(k))
            do i=i1+1,i2
               pe0(i ,k) = 0.5*(pe1(i,k)+pe1(i-1,k))
            enddo
          enddo

          do k=ks+2,km+1
! pe2(i1-1,ks+2:km+1) must be ghosted
            pe3(i1,k) = 0.5*(pe2(i1,k)+pe2w(k))
            do i=i1+1,i2
               pe3(i,k) = 0.5*(pe2(i,k)+pe2(i-1,k))
            enddo
          enddo

          call map1_ppm ( km,   pe0,    v,                             &
                          km,   pe3,    v,                             &
                          ngvs, ngvn, itot, i1-ifirst+1, i2-ifirst+1,  &
                          j,    jfirst, jlast, -1, kord)

        endif

! Save new PE to temp storage peln
        do k=2,km
          do i=i1,i2
             peln(i,k,j) = pe2(i,k)
          enddo
        enddo

! Check deformation.
       if( diag ) then
          rmax(ixj) = 0.
          rmin(ixj) = 1.
          do k=1,km
             do i=i1,i2
              ratio(i) = (pe1(i,k+1)-pe1(i,k)) / (pe2(i,k+1)-pe2(i,k))
             enddo

             do i=i1,i2
              if(ratio(i) > rmax(ixj)) then
                 rmax(ixj) = ratio(i)
              elseif(ratio(i) < rmin(ixj)) then
                 rmin(ixj) = ratio(i)
              endif
            enddo
          enddo
       endif
2000  continue


#if defined( SPMD )
! Send u southward
      call mp_send3d( iam-nprxy_x, iam+nprxy_x, im, jm, km,                  &
                      ifirst, ilast, jfirst-ngus, jlast+ngun, 1, km,         &
                      ifirst, ilast, jfirst, jfirst, 1, km, u )
      if (itot /= im) then
        dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
        src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
        call mp_send3d( dest, src, im, jm, km,                          &
                        ifirst, ilast, jfirst-ngvs, jlast+ngvn, 1, km,  &
                        ifirst, ifirst, jfirst, jlast, 1, km, v )
      endif
#endif

      if( diag ) then
        qmin = rmin(1)
        do ixj=2, jp
          if(rmin(ixj) < qmin) then
            qmin = rmin(ixj)
          endif
        enddo
        CPP_PRT_PREFIX write(6,*) 'rmin=', qmin

        qmax = rmax(1)
        do ixj=2, jp
          if(rmax(ixj) > qmax) then
            qmax = rmax(ixj)
          endif
        enddo
        CPP_PRT_PREFIX write(6,*) 'rmax=', qmax
      endif

!$omp  parallel do          &
!$omp  default(shared)      &
!$omp  private(i,j,k)

      do j=jfirst,jlast
        do k=2,km
          do i=ifirst,ilast
            pe(i,k,j) = peln(i,k,j)
          enddo
        enddo
      enddo

      call pkez(nxu, im, km, jfirst, jlast, 1, km, ifirst, ilast,  &
                pe, pk, akap, ks, peln, pkz, .true.)



! Single x-subdomain case (periodic)
      do k = 1, km
      do j = jfirst, jlast
        veast(j,k) = v(ifirst,j,k)
      enddo
      enddo

#if defined( SPMD )
! Recv u from north
      call mp_recv3d( iam+nprxy_x, im, jm, km,                               &
                      ifirst, ilast, jfirst-ngus, jlast+ngun, 1, km,         &
                      ifirst, ilast, jlast+1, jlast+1, 1, km, u )
      if (itot /= im) then
        call mp_recv3d( src, im, jm, km,                                &
                        ilast+1, ilast+1, jfirst, jlast, 1, km,         &
                        ilast+1, ilast+1, jfirst, jlast, 1, km, veast )
      endif
#endif

! ((((((((((((((((( compute globally integrated TE >>>>>>>>>>>>>>>>

      if( consv ) then

!$omp  parallel do         &
!$omp  default(shared)     &
!$omp  private(i,j,k)

        do k=1,km
          do j=jfirst,jlast
             do i=ifirst,ilast
                dz(i,j,k) = te(i,j,k) * delp(i,j,k)
             enddo
          enddo
        enddo

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(i,j,k,bte)

! Perform vertical integration

        do 4000 j=jfirst,jlast

          if ( j == 1 ) then
! SP
            tte(1) = 0.
  
            do k=1,km
              tte(1) = tte(1) + dz(ifirst,1,k)
            enddo

          elseif ( j == jm) then
! NP
            tte(jm) = 0.

            do k=1,km
              tte(jm) = tte(jm) + dz(ifirst,jm,k)
            enddo

          else
! Interior
            do i=ifirst,ilast
              bte(i) = 0.
            enddo

            do k=1,km
              do i=ifirst,ilast
                bte(i) = bte(i) + dz(i,j,k)
              enddo
            enddo

            xysum(j,1) = 0.
            do i=ifirst,ilast
              xysum(j,1) = xysum(j,1) + bte(i)
              tmpij(i,j,1) = bte(i)
            enddo

          endif
4000    continue

#if defined (SPMD)
        if (nprxy_x > 1) then
          call par_xsum(tmpij, ifirst, ilast, im, jlast-jfirst+1, xysum)
        endif
#endif

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(j)

        do j = max(jfirst,2), min(jlast,jm-1)
           tte(j) = xysum(j,1)*cosp(j)
        enddo

        if ( jfirst == 1 ) tte(1)  = acap * tte(1)
        if ( jlast == jm ) tte(jm) = acap * tte(jm)

        te1 = 0.
        call par_vecsum(jm, jfirst, jlast, tte, te1, comm_use, npry_use)

      endif   ! consv

      if( consv ) then

!$omp  parallel do       &
!$omp& default(shared)   &
!$omp& private(i,j)
 
       do j=js2g0, jn2g0
          xysum(j,1) = 0.
          xysum(j,2) = 0.
        do i=ifirst,ilast
          xysum(j,1) = xysum(j,1) + ps(i,j)
          xysum(j,2) = xysum(j,2) + peln(i,km+1,j)
          tmpij(i,j,1) = ps(i,j)
          tmpij(i,j,2) = peln(i,km+1,j) 
        enddo
       enddo

#if defined( SPMD )
       if (nprxy_x > 1) then
          call par_xsum(tmpij, ifirst, ilast, im, 2*(jlast-jfirst+1), xysum)
       endif
#endif

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(j)
 
       do j=js2g0, jn2g0
        tte(j) = cp*cosp(j)*(xysum(j,1) - ptop*float(im) -           &
                 akap*ptop*(xysum(j,2) - peln(ifirst,1,j)*float(im)) )
! peln(i,1,j) should be independent of i (AAM)
       enddo

       if ( jfirst == 1 ) tte(1) = acap*cp * (ps(ifirst,1) - ptop -    &
                akap*ptop*(peln(ifirst,km+1,1) - peln(ifirst,1,1) ) )
       if ( jlast == jm ) tte(jm)= acap*cp * (ps(ifirst,jm) - ptop -   &
                akap*ptop*(peln(ifirst,km+1,jm) - peln(ifirst,1,jm) ) )
      endif ! consv

      if (consv) then

       sum=0.
       call par_vecsum(jm, jfirst, jlast, tte, sum, comm_use, npry_use)

       dtmp = (te0 - te1) / sum
       if( diag ) then
         CPP_PRT_PREFIX write(6,*) 'te=',te0, ' Energy deficit in T = ', dtmp
       endif

      endif              ! end consv check

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i,j,k, u2, v2)

      do 8000 k=1,km
! Compute KE
        do j=js2g0,jn1g1
          do i=ifirst,ilast
            u2(i,j) = u(i,j,k)**2
          enddo
        enddo

        do j=js2g0,jn2g0
          do i=ifirst,ilast
            v2(i,j) = v(i,j,k)**2
          enddo
          v2(ilast+1,j) = veast(j,k)**2
        enddo

        do j=js2g0,jn2g0
          do i=ifirst,ilast
            te(i,j,k) = te(i,j,k) - 0.25 * ( u2(i,j) + u2(i,j+1)     &
                                            +v2(i,j) + v2(i+1,j) )
          enddo
        enddo

        if ( jfirst == 1 ) then
! South pole
          do i=ifirst,ilast
            u2_sp(i,k) = u2(i,2)
            v2_sp(i,k) = v2(i,2)
          enddo
        endif

        if ( jlast == jm ) then
! North pole
          do i=ifirst,ilast
            u2_np(i,k) = u2(i,jm)
            v2_np(i,k) = v2(i,jm-1)
          enddo
        endif

8000  continue

      if ( jfirst == 1 ) then

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i, k)

         do k = 1, km
            te_sp(k) = 0.
            do i=ifirst,ilast
              tmpik(i,k) = u2_sp(i,k) + v2_sp(i,k)
              te_sp(k) = te_sp(k) + tmpik(i,k)
            enddo
         enddo

#if defined( SPMD )
         if (nprxy_x > 1) then
            call par_xsum(tmpik, ifirst, ilast, im, km, te_sp)
         endif
#endif

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i, k)

         do k = 1, km
            te_sp(k) = te(ifirst,1,k) - 0.5*te_sp(k)/float(im)
            do i=ifirst,ilast
              te(i,  1,k) = te_sp(k)
            enddo
         enddo
      endif

      if ( jlast == jm ) then

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i, k)

         do k = 1, km
            te_np(k) = 0.
            do i=ifirst,ilast
              tmpik(i,k) = u2_np(i,k) + v2_np(i,k)
              te_np(k) = te_np(k) + tmpik(i,k)
            enddo
         enddo

#if defined( SPMD )
         if (nprxy_x > 1) then
            call par_xsum(tmpik, ifirst, ilast, im, km, te_np)
         endif
#endif

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i, k)

         do k = 1, km
            te_np(k) = te(ifirst,jm,k) - 0.5*te_np(k)/float(im)
            do i=ifirst,ilast
              te(i,jm,k) = te_np(k)
            enddo
         enddo
      endif

! Recover (virtual) temperature
!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(ixj, i1, i2, i, j, k, rg, gz, tvm, dlnp)

!     do 9000 j=jfirst,jlast
      do 9000 ixj=1,jp

         j  = jfirst + (ixj-1) / nxu
         i1 = ifirst + it * mod(ixj-1, nxu)
         i2 = i1 + it - 1

         rg = akap * cp
         do i=i1,i2
            gz(i) = hs(i,j)      
         enddo

        do k=km,1,-1
          do i=i1,i2
            dlnp  = rg*(peln(i,k+1,j) - peln(i,k,j))
            tvm   = delp(i,j,k)*(te(i,j,k) - gz(i)) /     &
                     ( cp*delp(i,j,k) - pe(i,k,j)*dlnp )
! Update phis
            gz(i) = gz(i) + dlnp*tvm
            pt(i,j,k) = tvm         ! pt is now (virtual) temperature
          enddo

          if( consv ) then
              do i=i1,i2
                 pt(i,j,k) = pt(i,j,k) + dtmp
              enddo
          endif

          if( .not. convt ) then
              do i=i1,i2
                 pt(i,j,k) = pt(i,j,k) / pkz(i,j,k)
              enddo
          endif
        enddo           ! end k-loop
9000  continue

      return
!EOC
      end subroutine te_map
!-----------------------------------------------------------------------
