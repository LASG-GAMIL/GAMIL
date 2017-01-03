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
  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use rgrid
  use commap
  use comfm1,       only: ghs
  use constituents, only: ppcnst, qmin, qmincg
  use time_manager, only: get_step_size, dtdy
  use stdatm
  use comhd
  use fspan

  implicit none

#include <comctl.h>
#include <comhyb.h>

  integer i, j, begj, m

  real(r8) pi             ! Mathematical pi (3.14...)
  real(r8) north,south

  real(r8) dtime          ! timestep size [seconds]

  pi = 4.0*atan(1.0)

  begj = beglatexdyn

  CALL INITIALIZE_STDATM
  CALL STDATM0(TBB, CBB, DCBB, HBB, P00, T00)
  CALL SETMSA0(TBB, HBB, GHS , P00, T00, PSB, TSB)

  dtime = get_step_size()

  P00   = P00  * 100.0D0

!
! calculate hypi & hypm for 'inti'
!
  call hycoef
!
! Initialize commap.
!
  call initialize_fspan
  CALL SPAN(MM1, MM2, MM3, MP1, MP2, MP3, MDJ)

  CALL LATMESH(DY, YTHU(1), YTHV(1), WTGU(1), WTGV(1))

  w(1)    = 1-cos(0.5*ythu(2))
  w(plat) = w(1) 

  do j = 2, plat/2
    north = 0.5*( ythu(j-1)+ythu(j) )
    south = 0.5*( ythu(j+1)+ythu(j) )
    w(j) = cos(north)-cos(south)
    w(plat+1-j) = w(j)
  end do

  do j = 1, plat
    clat(j) = ythu(j)-0.5d0*pi     
  end do
 
  do j = 1, plat
    latdeg(j) = clat(j)*45./atan(1._r8)
  end do

  dx = pi*2.0/dble(plon)

  call initialize_hpar
  CALL HPAR(DX, DY, YTHU(BEGJ), YTHV(BEGJ), WTGU(BEGJ), WTGV(BEGJ), MDJ, &
            SINU, SINV, OUX, OUY, OVX, OVY, FF, CUR)
!
! Set parameters for horizontal diffusion
!
  call initialize_comhd             !!(wh 2003.10.23)

  dthdfs = dtime 
     
  CALL STDFSC(DFS0, DTHDFS, SINU, SINV, WTGU(BEGJ), WTGV(BEGJ), DX, DY, &
              FRDT, FRDS, FRDU, FRDV, FRDP, DXVPN, DXVPS)
!
! Set minimum mixing ratio for moisture and advected tracers
!
  qmin(1) = 1.e-12          ! Minimum mixing ratio for moisture
  do m = 2, ppcnst
    qmin(m) = 0.0
  end do
!
! Set the minimum mixing ratio for the counter-gradient term.  
! Normally this should be the same as qmin, but in order to 
! match control case 414 use zero for water vapor.
!
  qmincg(1) = 0.
  do m = 2, ppcnst
    qmincg(m) = qmin(m)
  end do
!
!
! Determine whether full or reduced grid
!
  fullgrid = .true.
  do j = 1, plat
    if (nlon(j) < plon) fullgrid = .false.
  end do
!
!
! Longitude array
!
  do j = 1, plat
    do i = 1, nlon(j)
      londeg(i,j) = (i-1)*360./nlon(j)
      clon(i,j)   = (i-1)*2.0*pi/nlon(j)
    end do
  end do

!
! Set flag indicating dynamics grid is now defined.
! NOTE: this ASSUMES initcom is called after spmdinit.  The setting of nlon done here completes
! the definition of the dynamics grid.
!
  dyngrid_set = .true.
end subroutine initcom


