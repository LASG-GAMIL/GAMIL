#include <misc.h>
#include <params.h>
subroutine radini(gravx   ,cpairx  ,epsilox ,stebolx, pstdx )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize various constants for radiation scheme; note that
! the radiation scheme uses cgs units.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: W. Collins (H2O parameterization) and J. Kiehl
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pver, pverp
   use comozp,       only: cplos, cplol
   use pmgrid,       only: masterproc, plev, plevp, plevstd,plond,plat
   use radae,        only: radaeini
   use physconst,    only: mwdry, mwco2
   use constituents, only: co2vmr
#if ( defined SPMD )
    use mpishorthand
#endif
   implicit none

#include <crdcon.h>
#include <comhyb.h>
#include <comctl.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: gravx      ! Acceleration of gravity (MKS)
   real(r8), intent(in) :: cpairx     ! Specific heat of dry air (MKS)
   real(r8), intent(in) :: epsilox    ! Ratio of mol. wght of H2O to dry air
   real(r8), intent(in) :: stebolx    ! Stefan-Boltzmann's constant (MKS)
   real(r8), intent(in) :: pstdx      ! Standard pressure (Pascals)
!
!---------------------------Local variables-----------------------------

   real(r8), parameter :: rmwco2 = mwco2/mwdry   ! ratio of molecular weights of co2 to dry air
!
   integer k       ! Loop variable

   real(r8) v0         ! Volume of a gas at stp (m**3/kmol)
   real(r8) p0         ! Standard pressure (pascals)
   real(r8) amd        ! Effective molecular weight of dry air (kg/kmol)
   real(r8) goz        ! Acceleration of gravity (m/s**2)
!
!-----------------------------------------------------------------------
!
! Set general radiation consts; convert to cgs units where appropriate:
!
   gravit  =  100.*gravx
   rga     =  1./gravit
   cpair   =  1.e4*cpairx
   epsilo  =  epsilox
   sslp    =  1.013250e6
   stebol  =  1.e3*stebolx
   rgsslp  =  0.5/(gravit*sslp)
   dpfo3   =  2.5e-3
   dpfco2  =  5.0e-3
   dayspy  =  365.
   pie     =  4.*atan(1.)
!
! Determine co2 global mass mixing ratio - only if not ramping co2vmr
! Note: when ramping co2vmr, co2mmr is computed in subroutine ramp_ghg
!
   if (.not. doRamp_ghg) then
      co2mmr = rmwco2 * co2vmr
   endif
!
! Initialize ozone data.
!
   v0  = 22.4136         ! Volume of a gas at stp (m**3/kmol)
   p0  = 0.1*sslp        ! Standard pressure (pascals)
   amd = 28.9644         ! Molecular weight of dry air (kg/kmol)
   goz = gravx           ! Acceleration of gravity (m/s**2)
!
! Constants for ozone path integrals (multiplication by 100 for unit
! conversion to cgs from mks):
!
   cplos = v0/(amd*goz)       *100.0
   cplol = v0/(amd*goz*p0)*0.5*100.0
!
! Derived constants
! mxaerl = max number of levels (from bottom) for background aerosol
! Limit background aerosol height to regions below 900 mb
!
   mxaerl = 0
   do k=pver,1,-1
      if (hypm(k) >= 9.e4) mxaerl = mxaerl + 1
   end do
   mxaerl = max(mxaerl,1)
   if (masterproc) then
      write(6,*)'RADINI:  Background aerosol will be limited to ', &
                'bottom ',mxaerl,' model interfaces. Top interface is ', &
                hypi(pverp-mxaerl),' pascals'
   end if

! If the top model level is above ~90 km (0.1 Pa), set the top level to compute
! longwave cooling to about 80 km (1 Pa)
   if (hypm(1) .lt. 0.1) then
      do k = 1, pver
         if (hypm(k) .lt. 1.) ntoplw  = k
      end do
   else
      ntoplw = 1
   end if
   if (masterproc) then
      write (6,*) 'RADINI: ntoplw =',ntoplw, ' pressure:',hypm(ntoplw)
   endif

   call radaeini( pstdx, mwdry, mwco2 )
   return
end subroutine radini
