#include <misc.h>
#include <preproc.h>

subroutine FrictionVelocity (displa, z0m,   z0h,   z0q,   obu, &
                             iter, ur, um, ustar, temp1, temp2, clm)

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
! Calculation of the friction velocity, relation for potential 
! temperature and humidity profiles of surface boundary layer. 
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
! $Id: FrictionVelocity.F90,v 1.3.6.5.6.1 2002/10/03 20:07:27 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varcon, only : vkc
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm	 !CLM 1-D Module

  real(r8), intent(in) :: displa ! displacement height [m]
  real(r8), intent(in) :: z0m    ! roughness length, momentum [m]
  real(r8), intent(in) :: z0h    ! roughness length, sensible heat [m]
  real(r8), intent(in) :: z0q    ! roughness length, latent heat [m]
  real(r8), intent(in) :: obu    ! monin-obukhov length [m]
  real(r8), intent(in) :: um     ! wind speed including the stablity effect [m/s]
  real(r8), intent(in) :: ur     ! wind speed [m/s]
  integer , intent(in) :: iter   ! interation number [-]

  real(r8), intent(out) :: ustar ! friction velocity [m/s]
  real(r8), intent(out) :: temp1 ! relation for potential temperature profile
  real(r8), intent(out) :: temp2 ! relation for specific humidity profile

!----Local Variables----------------------------------------------------

  real(r8) zldis   ! reference height "minus" zero displacement height [m]
  real(r8) StabilityFunc ! stability function for unstable case
  real(r8) zetam   ! transition point of flux-gradient relation (wind profile) [-]
  real(r8) zetat   ! transition point of flux-gradient relation (temp. profile) [-]
  real(r8) zeta    ! dimensionless height used in Monin-Obukhov theory [-]

!----End Variable List--------------------------------------------------

!
! Adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.
! Wind profile
!

  zldis=clm%forc_hgt_u-displa
  zeta=zldis/obu
  zetam=1.574

  if (zeta < -zetam) then           ! zeta < -1
     ustar=vkc*um/(log(-zetam*obu/z0m)- &
          StabilityFunc(1,-zetam) +StabilityFunc(1,z0m/obu) &
          +1.14*((-zeta)**0.333-(zetam)**0.333))
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     ustar=vkc*um/(log(zldis/z0m)- &
          StabilityFunc(1,zeta)+StabilityFunc(1,z0m/obu))
  else if (zeta <= 1.) then        !  0 <= zeta <= 1
     ustar=vkc*um/(log(zldis/z0m) + &
          5.*zeta -5.*z0m/obu)
  else                             !  1 < zeta, phi=5+zeta
     ustar=vkc*um/(log(obu/z0m)+5.-5.*z0m/obu &
          +(5.*log(zeta)+zeta-1.))
  endif

#if (defined PERGRO)
  if (zeta < -zetam) then           ! zeta < -1
     ustar=vkc*um/log(-zetam*obu/z0m)
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     ustar=vkc*um/log(zldis/z0m)
  else if (zeta <= 1.) then        !  0 <= zeta <= 1
     ustar=vkc*um/log(zldis/z0m)
  else                             !  1 < zeta, phi=5+zeta
     ustar=vkc*um/log(obu/z0m)
  endif
#endif


!
! Temperature profile
!

  zldis=clm%forc_hgt_t-displa
  zeta=zldis/obu
  zetat=0.465
  if (zeta < -zetat) then           ! zeta < -1
     temp1=vkc/(log(-zetat*obu/z0h)-StabilityFunc(2,-zetat) &
          + StabilityFunc(2,z0h/obu) &
          + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     temp1=vkc/(log(zldis/z0h) - StabilityFunc(2,zeta) + &
          StabilityFunc(2,z0h/obu))
  else if (zeta <= 1.) then        !  0 <= zeta <= 1
     temp1=vkc/(log(zldis/z0h) + 5.*zeta - 5.*z0h/obu)
  else                             !  1 < zeta, phi=5+zeta
     temp1=vkc/(log(obu/z0h) + 5. - 5.*z0h/obu &
          + (5.*log(zeta)+zeta-1.))
  endif

#if (defined PERGRO)
  if (zeta < -zetat) then           ! zeta < -1
     temp1=vkc/log(-zetat*obu/z0h)
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     temp1=vkc/log(zldis/z0h)
  else if (zeta <= 1.) then        !  0 <= zeta <= 1
     temp1=vkc/log(zldis/z0h)
  else                             !  1 < zeta, phi=5+zeta
     temp1=vkc/log(obu/z0h)
  endif
#endif

!
! Humidity profile
!

  zldis=clm%forc_hgt_q-displa
  zeta=zldis/obu
  zetat=0.465
  if (zeta < -zetat) then          ! zeta < -1
     temp2=vkc/(log(-zetat*obu/z0q) - &
          StabilityFunc(2,-zetat) + StabilityFunc(2,z0q/obu) &
          + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     temp2=vkc/(log(zldis/z0q) - &
          StabilityFunc(2,zeta)+StabilityFunc(2,z0q/obu))
  else if (zeta <= 1.) then        !  0 <= zeta <= 1
     temp2=vkc/(log(zldis/z0q)+5.*zeta-5.*z0q/obu)
  else                             !  1 < zeta, phi=5+zeta
     temp2=vkc/(log(obu/z0q) + 5. - 5.*z0q/obu &
          + (5.*log(zeta)+zeta-1.))
  endif

#if (defined PERGRO)
  if (zeta < -zetat) then          ! zeta < -1
     temp2=vkc/log(-zetat*obu/z0q)
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     temp2=vkc/log(zldis/z0q)
  else if (zeta <= 1.) then        !  0 <= zeta <= 1
     temp2=vkc/log(zldis/z0q)
  else                             !  1 < zeta, phi=5+zeta
     temp2=vkc/log(obu/z0q)
  endif
#endif

end subroutine FrictionVelocity
