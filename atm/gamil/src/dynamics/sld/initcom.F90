#include <misc.h>
#include <params.h>

subroutine initcom

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize Model commons, including COMCON, COMHYB, COMMAP, COMSPE,
! and COMTRCNM
! 
! Method: 
! 
! Author: 
! Original version:  CCM1
! Standardized:      L. Bath, Jun 1992
!                    L. Buja, Feb 1996
!
!-----------------------------------------------------------------------
!
! $Id: initcom.F90,v 1.18.2.2 2002/06/15 13:48:26 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use rgrid
   use gauaw_mod, only: gauaw
   use commap
   use dynconst, only: rearth, ra, dynconsti
   use comslt
   use physconst, only: rair, rga
   use constituents, only: ppcnst, qmin, qmincg
   use time_manager, only: get_step_size
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comfft.h>
!-----------------------------------------------------------------------
#include <comhyb.h>
!-----------------------------------------------------------------------
!
! Local workspace
!
   real(r8) zsi(plat)      ! sine of latitudes
   real(r8) zw(plat)       ! Gaussian weights
   real(r8) zra2           ! ra squared
   real(r8) zalp(2*pspt)   ! Legendre function array
   real(r8) zdalp(2*pspt)  ! Derivative array
   real(r8) zslat          ! sin of lat  and cosine of colatitude
   real(r8) tmp            ! tmp variable
   real(r8) tmp1           ! tmp variable
   real(r8) psref          ! reference surf pres in Hortal calc.
   real(r8) pref           ! reference pres at level k in Hortal calc
   real(r8) tnot           ! reference surf temp in Hortal calc.
   real(r8) tref           ! reference temp at level k in Hortal calc.
   real(r8) lapse          ! standard atmosphere lapse rate
   real(r8) tcut           ! cutoff temperature below which "Hortal"

   integer i           ! longitude index
   integer j           ! Latitude index
   integer k           ! Level index
   integer kk          ! Level index
   integer kkk         ! Level index
   integer m           ! Index for legendre array
#ifdef PVP
   integer n           ! Index for legendre array
#endif
   integer nkk         ! Print control variables
   integer ik1         ! Print index temporary variable
   integer ik2         ! Print index temporary variable
   integer itmp        ! Dimension of polynomial arrays temporary.
   integer iter        ! Iteration index
   real(r8)    zdt         ! Time step for settau

   logical lprint      ! Debug print flag
   integer irow        ! Latitude pair index
   integer lat         ! Latitude index

   real(r8) xlat           ! Latitude (radians)
   real(r8) pi             ! Mathematical pi (3.14...)
   real(r8) dtime          ! timestep size [seconds]
!
!-----------------------------------------------------------------------
   call dynconsti
!
   lprint = masterproc .and. .FALSE.

   dtime = get_step_size()

   call hdinti  (rearth  ,dtime   )

   epssld = 0.2
!
! Initialize commap.  Set hybrid level dependent arrays
!
   call hycoef
!
! NMAX dependent arrays
!
   if (pmmax.gt.plon/2) then
      write(6,*)'INITCOM:mmax=ptrm+1 .gt. plon/2'
      call endrun
   end if
   zra2 = ra*ra
   do j=2,pnmax
      sq(j)  = j*(j-1)*zra2
      rsq(j) = 1./sq(j)
   end do
   sq(1)  = 0.
   rsq(1) = 0.
!
! MMAX dependent arrays
!
   do j=1,pmmax
      xm(j) = j-1
   end do
!
! Gaussian latitude dependent arrays
!
   call gauaw(zsi     ,zw      ,plat    )
   do irow=1,plat/2
      slat(irow) = zsi(irow)
      w(irow)              = zw(irow)
      w(plat - irow + 1)   = zw(irow)
      cs(irow)  = 1. - zsi(irow)*zsi(irow)
      xlat = asin(slat(irow))
      clat(irow) = -xlat
      clat(plat - irow + 1) = xlat
   end do

   do lat=1,plat
      latdeg(lat) = clat(lat)*45./atan(1._r8)
   end do
!
! Integration matrices of hydrostatic equation(href) and conversion
! term(ecref).  ecref calculated to be consistent with continuity Eq.;
! href calculated to conserve energy.
!
   do k=1,plev
      do kk=1,plev
         href(kk,k) = 0.
         ecref(kk,k) = 0.
      end do
   end do
!
! Mean atmosphere energy conversion term is consistent with continiuty
! Eq.  In ecref, 1st index = column; 2nd index = row of matrix.
! Mean atmosphere energy conversion term is energy conserving
!
   do k=1,plev
      ecref(k,k) = 0.5/hypm(k) * hypd(k)
      do kk=1,k-1
         ecref(kk,k) = 1./hypm(k) * hypd(kk)
      end do
   end do
!
! Reference hydrostatic integration matrix consistent with conversion
! term for energy conservation.  In href, 1st index = column; 
! 2nd index = row of matrix.
!
   do k = 1,plev
      do kk = k,plev
         href(kk,k) = ecref(k,kk)*hypd(kk)/hypd(k)
      end do
   end do
!
! Print statements
!
   if (lprint) then
      nkk = plev/13
      if (mod(plev,13).ne.0) nkk = nkk + 1
      write(6,*)' '
      write(6,*)'INITCOM: Hydrostatic matrix href'
      do kk=1,nkk
         ik1 = 1 + (kk-1)*13
         ik2 = min0( ik1+12, plev )
         write(6,9920) (k,k=ik1,ik2)
         do kkk=1,plev
            write(6,9910) kkk,(href(kkk,k),k=ik1,ik2)
         end do
      end do
      write(6,*)' '
      write(6,*)'INITCOM: Thermodynamic matrix ecref'
      do kk=1,nkk
         ik1 = 1 + (kk-1)*13
         ik2 = min0( ik1+12, plev )
         write(6,9920) (k,k=ik1,ik2)
         do kkk=1,plev
            write(6,9910) kkk,(ecref(kkk,k),k=ik1,ik2)
         end do
      end do
   end if
!
! Multiply href by r
!
   do k=1,plev
      do kk=1,plev
         href(kk,k) = href(kk,k)*rair
      end do
   end do
!
! Compute truncation parameters
!
   if (masterproc) then
      write(6,9950) ptrm,ptrn,ptrk
   end if
!
! Compute semi-implicit timestep constants (COMSPE)
!
   zdt = 0.5*dtime
!
! The CMIC$ DO ALL ... construct is a "phony loop" to fool the low level
! Cray matrix library utilities into *not* multitasking, since these 
! utilities give DIFFERENT answers for different values of $NCPUS.  Useful 
! work is done only for iter=1.
!
!MIC$ DO ALL PRIVATE (ZDT, ITER)
   do iter=1,2
      call settau(zdt, iter)
   end do
!
! Set minimum mixing ratio for moisture and advected tracers
!
   qmin(1) = 1.e-12          ! Minimum mixing ratio for moisture
   do m=2,ppcnst
      qmin(m) = 0.0
   end do
!
! Set the minimum mixing ratio for the counter-gradient term.  
! Normally this should be the same as qmin, but in order to 
! match control case 414 use zero for water vapor.
!
   qmincg(1) = 0.
   do m=2,ppcnst
      qmincg(m) = qmin(m)
   end do
!
! Compute constants related to Legendre transforms
! Compute and reorder ALP and DALP
!
   do j=1,plat/2
      zslat = slat(j)
      itmp = 2*pspt - 1
      call phcs  (zalp    ,zdalp   ,itmp    ,zslat    )
      call reordp(j       ,itmp    ,zalp    ,zdalp   )
   end do
!
! Initialize array of constants for Hortal temperature correction
!
   psref = 101325.
   tnot  = 288.
   lapse = 0.0065
   tcut  = 216.5
   tmp   = rair*lapse*rga
!
   do k = 1,plev
      pref = hyam(k)*ps0 + hybm(k)*psref
      tmp1 = pref/psref
      tref = tnot*tmp1**tmp
      if (tref .gt. tcut) then
         hortalc(k) = hybm(k)*tmp*tref/tmp1
      else
         hortalc(k) = 0.
      endif
   end do
!
! Compute vertical difference
!
   hdel(0) = hortalc(1)
   do k = 1,plev-1
      hdel(k) = hortalc(k+1) - hortalc(k)
   end do
!
! Determine whether full or reduced grid
!
   fullgrid = .true.
   do j=1,plat
      if (masterproc) then
         write(6,*)'nlon(',j,')=',nlon(j),' wnummax(',j,')=',wnummax(j)
      end if
      if (nlon(j).lt.plon) fullgrid = .false.
   end do
!
! Mirror latitudes south of south pole
!
   lat = 1
   do j=j1-2,1,-1
      nlonex(j) = nlon(lat)
      lat = lat + 1
   end do
   nlonex(j1-1) = nlon(1)     ! south pole
!
! Real latitudes
!
   j = j1
   do lat=1,plat
      nlonex(j) = nlon(lat)
      j = j + 1
   end do
   nlonex(j1+plat) = nlon(plat)  ! north pole
!
! Mirror latitudes north of north pole
!
   lat = plat
   do j=j1+plat+1,platd
      nlonex(j) = nlon(lat)
      lat = lat - 1
   end do
!
! Longitude array
!
   pi = 4.0*atan(1.0)
   do lat=1,plat
      do i=1,nlon(lat)
         londeg(i,lat) = (i-1)*360./nlon(lat)
         clon(i,lat)   = (i-1)*2.0*pi/nlon(lat)
      end do
   end do

   do j=1,plat/2
      nmmax(j) = wnummax(j) + 1
   end do
#ifdef PVP
   do irow=1,plat/2
      do n=1,pmax
         nmreduced(n,irow) = min(nm(n),nmmax(irow))
      end do
   end do
#endif
   do m=1,pmmax
      do irow=1,plat/2
         if (nmmax(irow) .ge. m) then
            beglatpair(m) = irow
            goto 10
         end if
      end do
      write(6,*)'INITCOM: Should not ever get here'
      call endrun
10    continue
   end do
!
! Set up trigonometric tables for fft
!
   do j=1,plat
      call set99(trig(1,j),ifax(1,j),nlon(j))
   end do
!
! Set flag indicating dynamics grid is now defined.
! NOTE: this ASSUMES initcom is called after spmdinit.  The setting of nlon done here completes
! the definition of the dynamics grid.
!
   dyngrid_set = .true.

   return

9910 format( 1x,i3,13f9.5)
9920 format(/,      13i9)
9950 format(/,'     Truncation Parameters',/,'     NTRM = ',i4,/, &
      '     NTRN = ',i4,/,'     NTRK = ',i4,/)

end subroutine initcom
