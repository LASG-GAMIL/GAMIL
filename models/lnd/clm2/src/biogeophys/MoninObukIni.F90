#include <misc.h>
#include <preproc.h>

subroutine MoninObukIni(ur, thv, dthv, zldis, z0m, &
                        um, obu  )

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely
!  L                           M  available land surface process model.
!  M --COMMUNITY LAND MODEL--  C
!  C                           L
!  LMCLMCLMCLMCLMCLMCLMCLMCLMCLM
!
!-----------------------------------------------------------------------
! Purpose:
! Initialization of the Monin-Obukhov length.
!
! Method:
! The scheme is based on the work of Zeng et al. (1998): 
! Intercomparison of bulk aerodynamic algorithms for the computation 
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate, 
! Vol. 11, 2628-2644.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: MoninObukIni.F90,v 1.2.10.3 2002/06/15 13:50:16 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varcon, only : grav
  implicit none

!----Arguments----------------------------------------------------------

  real(r8), intent(in) :: ur    ! wind speed at reference height [m/s]
  real(r8), intent(in) :: thv   ! virtual potential temperature [K]
  real(r8), intent(in) :: dthv  ! diff of thv between ref. height and surface [K]
  real(r8), intent(in) :: zldis ! reference height "minus" zero displacement height [m]
  real(r8), intent(in) :: z0m   ! roughness length, momentum [m]

  real(r8), intent(out) :: um   ! wind speed including the stability effect [m/s]
  real(r8), intent(out) :: obu  ! monin-obukhov length [m]

!----Local Variables----------------------------------------------------

  real(r8)  wc    ! convective velocity [m/s]
  real(r8)  rib   ! bulk Richardson number
  real(r8)  zeta  ! dimensionless height used in Monin-Obukhov theory [-]
  real(r8)  ustar ! friction velocity [m/s]     

!----End Variable List--------------------------------------------------

!
! Initial values of u* and convective velocity
!

  ustar=0.06
  wc=0.5
  if (dthv >= 0.) then
     um=max(ur,0.1_r8)
  else
     um=sqrt(ur*ur+wc*wc)
  endif

  rib=grav*zldis*dthv/(thv*um*um)
#if (defined PERGRO)
  rib = 0.
#endif

  if (rib >= 0.) then      ! neutral or stable
     zeta = rib*log(zldis/z0m)/(1.-5.*min(rib,0.19_r8))
     zeta = min(2._r8,max(zeta,0.01_r8 ))
  else                    !unstable
     zeta=rib*log(zldis/z0m)
     zeta = max(-100._r8,min(zeta,-0.01_r8 ))
  endif

  obu=zldis/zeta

end subroutine MoninObukIni
