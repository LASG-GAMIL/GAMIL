#include <misc.h>
#include <params.h>

subroutine cpslec (ncol, pmid, phis, ps, t, psl, gravit, rair)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Hybrid coord version:  Compute sea level pressure for a latitude line
! 
! Method: 
! CCM2 hybrid coord version using ECMWF formulation
! Algorithm: See section 3.1.b in NCAR NT-396 "Vertical 
! Interpolation and Truncation of Model-Coordinate Data
!
! Author: Stolen from the Processor by Erik Kluzek
! 
!-----------------------------------------------------------------------
!
! $Id: cpslec.F90,v 1.1.28.1 2002/06/15 13:49:14 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: pcols, pver

  implicit none

!-----------------------------Arguments---------------------------------
  integer , intent(in) :: ncol             ! longitude dimension

  real(r8), intent(in) :: pmid(pcols,pver) ! Atmospheric pressure (pascals)
  real(r8), intent(in) :: phis(pcols)      ! Surface geopotential (m**2/sec**2)
  real(r8), intent(in) :: ps(pcols)        ! Surface pressure (pascals)
  real(r8), intent(in) :: T(pcols,pver)    ! Vertical slice of temperature (top to bot)
  real(r8), intent(in) :: gravit           ! Gravitational acceleration
  real(r8), intent(in) :: rair             ! gas constant for dry air

  real(r8), intent(out):: psl(pcols)       ! Sea level pressures (pascals)
!-----------------------------------------------------------------------

!-----------------------------Parameters--------------------------------
  real(r8), parameter :: xlapse = 6.5e-3   ! Temperature lapse rate (K/m)
!-----------------------------------------------------------------------

!-----------------------------Local Variables---------------------------
  integer i              ! Loop index
  real(r8) alpha         ! Temperature lapse rate in terms of pressure ratio (unitless)
  real(r8) Tstar         ! Computed surface temperature
  real(r8) TT0           ! Computed temperature at sea-level
  real(r8) alph          ! Power to raise P/Ps to get rate of increase of T with pressure
  real(r8) beta          ! alpha*phis/(R*T) term used in approximation of PSL
!-----------------------------------------------------------------------
!
  alpha = rair*xlapse/gravit
  do i=1,ncol
     if ( abs(phis(i)/gravit) < 1.e-4 )then
        psl(i)=ps(i)
     else
        Tstar=T(i,pver)*(1.+alpha*(ps(i)/pmid(i,pver)-1.)) ! pg 7 eq 5

        TT0=Tstar + xlapse*phis(i)/gravit                  ! pg 8 eq 13

        if ( Tstar<=290.5 .and. TT0>290.5 ) then           ! pg 8 eq 14.1
           alph=rair/phis(i)*(290.5-Tstar)  
        else if (Tstar>290.5  .and. TT0>290.5) then        ! pg 8 eq 14.2
           alph=0.
           Tstar= 0.5 * (290.5 + Tstar)  
        else  
           alph=alpha  
           if (Tstar<255.) then  
              Tstar= 0.5 * (255. + Tstar)                  ! pg 8 eq 14.3
           endif
        endif

        beta = phis(i)/(rair*Tstar)
        psl(i)=ps(i)*exp( beta*(1.-alph*beta/2.+((alph*beta)**2)/3.))
     end if
  enddo

  return
end subroutine cpslec
