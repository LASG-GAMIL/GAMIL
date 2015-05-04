MODULE mo_cumulus_flux

  use shr_kind_mod, only: r8 => shr_kind_r8

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_cumulus_flux* - parameters for cumulus massflux scheme
  !
  ! ----------------------------------------------------------------

  REAL(r8) :: entrpen      !    entrainment rate for penetrative convection
  REAL(r8) :: entrscv      !    entrainment rate for shallow convection
  REAL(r8) :: entrmid      !    entrainment rate for midlevel convection
  REAL(r8) :: entrdd       !    entrainment rate for cumulus downdrafts
  REAL(r8) :: cmfctop      !    relat. cloud massflux at level above nonbuoyanc
  REAL(r8) :: cmfcmax      !    maximum massflux value allowed for
  REAL(r8) :: cmfcmin      !    minimum massflux value (for safety)
  REAL(r8) :: cmfdeps      !    fractional massflux for downdrafts at lfs
  REAL(r8) :: rhcdd        !    relative saturation in downdrafts
  REAL(r8) :: cprcon       !    coefficients for determining conversion
                       !    from cloud water to rain
  LOGICAL :: lmfpen    !    true if penetrative convection is switched on
  LOGICAL :: lmfscv    !    true if shallow     convection is switched on
  LOGICAL :: lmfmid    !    true if midlevel    convection is switched on
  LOGICAL :: lmfdd     !    true if cumulus downdraft      is switched on
  LOGICAL :: lmfdudv   !    true if cumulus friction       is switched on

CONTAINS

SUBROUTINE cuparam

  ! Description:
  !
  ! Defines disposable parameters for massflux scheme
  !
  ! Method:
  !
  ! This routine is called from *iniphy*
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, February 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, Jan 1999, subroutine cuparam -> module mo_cumulus_flux
  ! 
  ! for more details see file AUTHORS
  ! 

! USE mo_control, ONLY: lamip2

  IMPLICIT NONE


  !  Executable Statements 
!-- 0. For cam run
  lmfpen = .TRUE.
  lmfscv = .TRUE.
  lmfmid = .TRUE.
  lmfdd = .TRUE.
  lmfdudv = .TRUE.

!-- 1. Specify parameters for massflux-scheme

  entrpen = 1.0E-4 ! Average entrainment rate for penetrative convection

!!(ljli)  entrscv = 3.0E-4 ! Average entrainment rate for shallow convection
  entrscv = 1.2E-3 ! Average entrainment rate for shallow convection

  entrmid = 1.0E-4 ! Average entrainment rate for midlevel convection

  entrdd  = 2.0E-4 ! Average entrainment rate for downdrafts

! IF (lamip2) THEN
!!(ljli)    cmfctop = 0.33 ! Relative cloud massflux at level above nonbuoyancy level
    cmfctop = 0.28 ! Relative cloud massflux at level above nonbuoyancy level
! ELSE
!   cmfctop = 0.1
! END IF

  cmfcmax = 1.0    ! Maximum massflux value allowed for updrafts etc

  cmfcmin = 1.E-10 ! Minimum massflux value (for safety)

  cmfdeps = 0.3    ! Fractional massflux for downdrafts at lfs
!!(ljli)  cmfdeps = 0.2    ! Fractional massflux for downdrafts at lfs

!!(ljli)  cprcon  = 6.E-4  ! Coefficients for determining conversion from cloud water
  cprcon  = 2.E-4  ! Coefficients for determining conversion from cloud water

  ! Next value is relative saturation in downdrafrs
  ! but is no longer used ( formulation implies saturation)

  rhcdd = 1.

  RETURN
END SUBROUTINE cuparam

END MODULE mo_cumulus_flux
