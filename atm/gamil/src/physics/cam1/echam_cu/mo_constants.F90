MODULE mo_constants

  use shr_kind_mod, only: r8 => shr_kind_r8

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_constants* basic universal constants  and derived constants.
  !
  ! ----------------------------------------------------------------

  REAL(r8), PARAMETER :: dayl = 86400.  ! length of the day (in seconds).

  REAL(r8) :: api     ! 2.*arcsin(1.).
  REAL(r8) :: a       ! radius of the earth.
  REAL(r8) :: omega   ! solid rotation velocity of the earth.
  REAL(r8) :: g       ! gravity acceleration.
  REAL(r8) :: cpd     ! specific heat at constant pressure (dry air).
  REAL(r8) :: cpv     !            idem               (water vapour).
  REAL(r8) :: rd      ! gas constant for dry air.
  REAL(r8) :: rv      !    idem      for water vapour.
  REAL(r8) :: rcpd    ! rcpd=1./cpd.
  REAL(r8) :: vtmpc1  ! vtmpc1=rv/rd-1.
  REAL(r8) :: vtmpc2  ! vtmpc2=cpv/cpd-1.
  REAL(r8) :: rhoh2o  ! density of liquid water.
  REAL(r8) :: alv     ! latent heat for vaporisation.
  REAL(r8) :: als     ! latent heat for sublimation.
  REAL(r8) :: alf     ! latent heat for fusion.
  REAL(r8) :: clw     ! specific heat for liquid water.
  REAL(r8) :: tmelt   ! temperature of fusion of ice.
  REAL(r8) :: solc    ! solar constant.
  REAL(r8) :: stbo    ! stephan boltzmann constant.
  REAL(r8) :: yearl   ! length of the year (in days).

  ! constants used for computation of saturation mixing ratio
  !   over liquid water(*c_les*) or ice(*c_ies*).

  REAL(r8) :: c1es    ! 610.78
  REAL(r8) :: c2es    ! 1es*rd/rv
  REAL(r8) :: c3les   ! 17.269
  REAL(r8) :: c3ies   ! 21.875
  REAL(r8) :: c4les   ! 35.86
  REAL(r8) :: c4ies   !  7.66
  REAL(r8) :: c5les   ! c3les*(tmelt-c4les)
  REAL(r8) :: c5ies   ! c3ies*(tmelt-c4ies)
  REAL(r8) :: c5alvcp ! c5les*alv/cpd
  REAL(r8) :: c5alscp ! c5ies*als/cpd
  REAL(r8) :: alvdcp  ! alv/cpd
  REAL(r8) :: alsdcp  ! als/cpd

CONTAINS

  SUBROUTINE inicon

    ! Description:
    ! Preset constants in mo_constants.
    !
    ! Method:
    !
    ! *inicon* is called from *setdyn*.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, December 1982, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! H.-S. Bauer, MPI, Jul 1998, changed
    ! A. Rhodin, MPI, Jan 1999, subroutine inicon put into module mo_constants
    !
    ! for more details see file AUTHORS
    !

    IMPLICIT NONE

    !  Intrinsic functions 
    INTRINSIC ASIN

    !  Executable statements 

    !-- 1. Preset constants

    api    = 2.*ASIN(1.)
    a      = 6371000.
    omega  = .7292E-4
    g      = 9.80665
    cpd    = 1005.46
    cpv    = 1869.46
    rd     = 287.05
    rv     = 461.51

    rcpd   = 1./cpd
    vtmpc1 = rv/rd - 1.
    vtmpc2 = cpv/cpd - 1.

    rhoh2o = 1000.
    alv    = 2.5008E6
    als    = 2.8345E6
    alf    = als - alv

    clw    = 4186.84
    tmelt  = 273.16

    solc   = 1365.
    stbo   = 5.67E-8

    c1es    = 610.78
    c2es    = c1es*rd/rv
    c3les   = 17.269
    c3ies   = 21.875
    c4les   = 35.86
    c4ies   =  7.66
    c5les   = c3les*(tmelt-c4les)
    c5ies   = c3ies*(tmelt-c4ies)
    c5alvcp = c5les*alv/cpd
    c5alscp = c5ies*als/cpd
    alvdcp  = alv/cpd
    alsdcp  = als/cpd

  END SUBROUTINE inicon

END MODULE mo_constants
