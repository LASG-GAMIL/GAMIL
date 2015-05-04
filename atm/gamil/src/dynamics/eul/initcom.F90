#include <misc.h>
#include <params.h>

subroutine initcom

!! (wanhui 2003.04.30)
!! (wanhui 2003.07.10)
!! (wanhui 2004.04.14)
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
! $Id: initcom.F90,v 1.16.2.2 2002/06/15 13:47:48 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
!! use pspect
!! use comspe
   use rgrid
!! use gauaw_mod, only: gauaw
   use commap
!! use prognostics, only: phis
   use comfm1, only: ghs
!! use dynconst, only: rearth, ra, dynconsti
!! use physconst, only: rair
   use constituents, only: ppcnst, qmin, qmincg               !!(wh 2003.10.24)
!! use constituents, only: ppcnst, qmin, qmincg, &
!!                           isor, iord, ipq, dsnp,dssp, gc, dtdlt,dtdln,dtdsg
   use qadv, only:  isor, iord, ipq, dsnp,dssp, gc, dtdlt,dtdln,dtdsg,initialize_qadv !!(wh)
!! use time_manager, only: get_step_size
   use time_manager, only: get_step_size,dtdy  !! (wh 2004.04.14)
   use stdatm                           !!(wh 2003.10.23)
   use comhd                            !!(wh 2003.10.23)
   use fspan                            !!(wh 2003.10.24)
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
!!#include <comfft.h>
!-----------------------------------------------------------------------
#include <comhyb.h>
!-----------------------------------------------------------------------
!!#include <comhd.h>
!!----------------------------------------------------------------------
! Local workspace
!

   integer i           ! longitude index
   integer j           ! Latitude index
   integer begj

   integer m           ! lengendre array index

   logical lprint      ! Debug print flag
   integer lat         ! Latitude index

   real(r8) pi             ! Mathematical pi (3.14...)
   real(r8) north,south

   real(r8) dtime          ! timestep size [seconds]
!
!-----------------------------------------------------------------------
     pi = 4.0*atan(1.0)

     begj=beglatexdyn

!!
!! Set the fm2003 std.atm.
!!
     call initialize_stdatm                 !!(wh 2003.10.23)
     call stdatm0(TBB,CBB,DCBB,HBB,P00,T00 )
     call setmsa0(TBB,HBB,ghs ,P00,T00,PSB,TSB)

!
     lprint = masterproc .and. .FALSE.

     dtime = get_step_size()
!!
!! Set vertical layers
!!
     call vpar (pmtop,p00,sig,sigl,dsig )

     pmtop = pmtop * 100.0d0
     p00   = p00   * 100.0d0

!!
!! calculate hypi & hypm for 'inti'
!!
     call hycoef
    
!!
!! Initialize commap.
!!
     call initialize_fspan                   !!(wh 2003.10.24)
     call span (mm1,mm2,mm3,mp1,mp2,mp3,mdj)

!!     
     call latmesh (dy,ythu(1),ythv(1),wtgu(1),wtgv(1))


      w(1)    = 1-cos(0.5*ythu(2))
      w(plat) = w(1) 

      do j=2,plat/2
         north = 0.5*( ythu(j-1)+ythu(j) )
         south = 0.5*( ythu(j+1)+ythu(j) )
         w(j) = cos(north)-cos(south)
         w(plat+1-j) = w(j)
      enddo

      do j=1,plat
         clat(j) = ythu(j)-0.5d0*pi     
      enddo
 
      do lat=1,plat
         latdeg(lat) = clat(lat)*45./atan(1._r8)
      end do
!!
     dx = pi*2.0/dble(plon)

     call initialize_hpar                   !!(wh 2003.10.24)
     call hpar (dx,dy,ythu(begj),ythv(begj),wtgu(begj),wtgv(begj),mdj            &
                     ,sinu,sinv,oux,ouy,ovx,ovy,ff,cur)
!!
!!
!! Set parameters for horizontal diffusion
!!
     call initialize_comhd             !!(wh 2003.10.23)

!!   dfs0 = 0.02d0                     !!(namelist variable.  wh 2004.04.14)
     dthdfs = dtime 
     
     call stdfsc (dfs0,dthdfs,sinu,sinv,wtgu(begj),wtgv(begj),dx,dy   &
                   ,frdt,frds,frdu,frdv,frdp,dxvpn,dxvps)

!!
!! Set parameters for the advection integration of q-H2O
!!
     call initialize_qadv               !!(wh 2003.10.24)

     call conpda (dtdy,dx,dy,sinu,wtgv(begj),dsig,-1.0d0,isor,iord   &
                          ,ipq,dsnp,dssp,gc,dtdlt,dtdln,dtdsg) 

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
!
! Determine whether full or reduced grid
!
   fullgrid = .true.
   do j=1,plat
      if (masterproc) then
         write(6,*)'nlon(',j,')=',nlon(j)
      end if
      if (nlon(j).lt.plon) fullgrid = .false.
   end do
!
!
! Longitude array
!
   do lat=1,plat
      do i=1,nlon(lat)
         londeg(i,lat) = (i-1)*360./nlon(lat)
         clon(i,lat)   = (i-1)*2.0*pi/nlon(lat)
      end do
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

end subroutine initcom


