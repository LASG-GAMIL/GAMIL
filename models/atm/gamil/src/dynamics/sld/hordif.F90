#include <misc.h>
#include <params.h>
subroutine hordif(k,ztdt)
!-----------------------------------------------------------------------
!
! Purpose:
! Horizontal diffusion of z,d,t,q
! 1. implicit del**2 form above level kmnhd4
! 2. implicit del**4 form at level kmnhd4 and below
! 3. courant number based truncation at level kmxhdc and above
! 4. increased del**2 coefficient at level kmxhd2 and above
!
! Computational note: this routine is multitasked by level, hence it 
! is called once for each k
!
! Note: Most of the "ifdef" constructs employed in this routine are related
! to the fact that storage order for spectral coefficients is different
! depending upon whether the target architecture is PVP or not.  The token
! SPMD has an "ifdef" test associated with it since the message-passing 
! implementation of CAM distributes subregions of Fourier wavenumber space
! to individual processes.
!
! Original version:  CCM1
!
!-----------------------------------------------------------------------
!
! $Id: hordif.F90,v 1.4.2.1 2002/06/15 13:48:25 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  use time_manager, only: get_step_size, get_nstep
  implicit none
#include <comhd.h>
#include <comctl.h>

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: k    ! level index
  real(r8), intent(in)   :: ztdt ! 2 times time step unless nstep=0
!
!---------------------------Local workspace-----------------------------
!
#if ( defined PVP )
  real(r8) tqfac(2*pnmax)        ! time-split implicit diffusion factors (t,q)
  real(r8) zdfac(2*pnmax)        ! time-split implicit diffusion factors (z,d)
  integer isp,ne                 ! spectral indices
#else
  integer ir,ii                  ! spectral indices       
  integer mr,mc                  ! spectral indices
#endif
  real(r8) dfac                  ! large coefficient on del^n multipliers to
!                                ! strongly damp waves req'd by Courant limiter
  integer m,n                    ! spectral indices
  real(r8) zdt                   ! model time step
  real(r8) dmpini                ! used to compute divergence damp rate
  real(r8) dmptim                ! used to compute divergence damp rate
  real(r8) dmprat                ! divergence damping rate
  real(r8) coef                  ! coeff. used to apply damping rate to divergence
!
!-----------------------------------------------------------------------
!
! Set the horizontal diffusion factors for each wavenumer at this level
! depending on: whether del^2 or del^4 diffusion is to be applied; and 
! whether the courant number limit is to be applied.
!
  if (k .ge. kmnhd4) then        ! Del^4 diffusion factors
     do n=1,pnmax
        hdiftq(n,k) = hdfst4(n)
        hdifzd(n,k) = hdfsd4(n)
     end do
!
! Spectrally truncate selected levels (if courant number too large)
!
     if (k.le. kmxhdc .and. nindex(k).le.pnmax) then
        dfac = 1000.
        do n=nindex(k),pnmax
           hdiftq(n,k) = dfac*hdfst4(n)
           hdifzd(n,k) = dfac*hdfsd4(n)
        end do
     end if
  else                      ! Del^2 diffusion factors
     if (k.le.kmxhd2) then
        dfac = sqrt(2.)**(real(kmxhd2-k+1))
     else
        dfac = 1.0
     end if
     do n=1,pnmax
        hdiftq(n,k) = dfac*hdfst2(n)
        hdifzd(n,k) = dfac*hdfsd2(n)
     end do
!
! Spectrally truncate selected levels (if courant number too large)
!
     if ((k.le.kmxhdc).and.(nindex(k).le.pnmax)) then
        dfac = 1000.
        do n=nindex(k),pnmax
           hdiftq(n,k) = dfac*hdfst2(n)
           hdifzd(n,k) = dfac*hdfsd2(n)
        end do
     end if
  end if
!
! Define damping rate for divergence damper
!
  zdt = get_step_size()
!
! Initial damping rate (e-folding time = zdt) and then linearly decrease
! to 0. over number of days specified by "divdampn".
!
  coef = 1.
  if (divdampn .gt. 0.0) then
     dmpini = 1./(zdt)
     dmptim = divdampn*86400.
     dmprat = dmpini * (dmptim - float(get_nstep())*zdt) / dmptim
     if (dmprat .gt. 0.0) coef = 1.0 / (1.0+zdt*dmprat)
  endif
!
! Compute time-split implicit factors for this level
!
#if ( defined PVP )
  do n = 1, pnmax
     tqfac(2*n)   = 1./(1. + ztdt*hdiftq(n,k))
     tqfac(2*n-1) = 1./(1. + ztdt*hdiftq(n,k))
     zdfac(2*n)   = 1./(1. + ztdt*hdifzd(n,k))
     zdfac(2*n-1) = 1./(1. + ztdt*hdifzd(n,k))
  end do
!
! Apply the horizontal diffusion for this level
!
  do n=1,pmax
     isp = nco2(n) - 2
     ne = 2*(n-1)
     do m=1,2*nm(n)
        t(isp+m,k)  =  t(isp+m,k)*tqfac(ne+m)
        d(isp+m,k)  =  d(isp+m,k)*zdfac(ne+m)*coef
        vz(isp+m,k) = vz(isp+m,k)*zdfac(ne+m)
     end do
  end do
#else
  do m=begm(iam),endm(iam)
     mr = nstart(m)
     mc = 2*mr
     do n=1,nlen(m)
        ir = mc + 2*n - 1
        ii = ir + 1
!
! time-split implicit factors
!
        t(ir,k) = t(ir,k)/(1. + ztdt*hdiftq(n+m-1,k))
        t(ii,k) = t(ii,k)/(1. + ztdt*hdiftq(n+m-1,k))
!
        d(ir,k) = d(ir,k)*coef/(1. + ztdt*hdifzd(n+m-1,k))
        d(ii,k) = d(ii,k)*coef/(1. + ztdt*hdifzd(n+m-1,k))
!
        vz(ir,k) = vz(ir,k)/(1. + ztdt*hdifzd(n+m-1,k))
        vz(ii,k) = vz(ii,k)/(1. + ztdt*hdifzd(n+m-1,k))
     end do
  end do
#endif
!
  return
end subroutine hordif

