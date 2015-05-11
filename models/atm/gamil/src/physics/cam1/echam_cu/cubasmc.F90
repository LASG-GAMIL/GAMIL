!+  Calculates cloud base values for midlevel convection
!+ $Id: cubasmc.f90,v 1.5 1999/09/17 11:18:36 m214089 Exp $

SUBROUTINE cubasmc(klon,klev,kk,pten,pqen,pqsen,puen,pven,&
!ktrac, &
!&      pxten,pxtu,pmfuxt,
   pverv,pgeo,pgeoh,ldcum,ktype,klab,pmfu,pmfub,pentr, &
&      kcbot,ptu,pqu,plu,puu,pvu,pmfus,pmfuq,pmful,pdmfup,pmfuu,pmfuv)

  ! Description:
  !
  ! This routine calculates cloud base values for midlevel convection
  !
  ! Method:
  !
  ! This routine is called from *cuasc*.
  ! Input are environmental values t,q etc.
  ! It returns cloudbase values for midlevel convection
  !
  ! Authors:
  !
  ! S. Tiedtke, ECMWF, in 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !
  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE mo_constants,    ONLY: cpd,     &! specific heat at constant pressure
                             rcpd,    &! rcpd=1./cpd
                             g         ! gravity acceleration
  USE mo_cumulus_flux, ONLY: cmfcmin, &! minimum massflux value (for safety)
                             cmfcmax, &! maximum massflux value allowed for
                             entrmid, &! entrainment rate for midlevel convect.
                             lmfdudv   ! true if cumulus friction is switched on

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: klon, kk, klev
!, ktrac

  !  Array arguments with intent(In):
  REAL(r8), INTENT (IN) :: pgeo(klon,klev), pgeoh(klon,klev), pqen(klon,klev),  &
&                      pqsen(klon,klev), pten(klon,klev), puen(klon,klev), &
&                      pven(klon,klev), pverv(klon,klev)
!, pxten(klon,klev,ktrac)
  LOGICAL, INTENT (IN) :: ldcum(klon)

  !  Array arguments with intent(InOut):
  REAL(r8), INTENT (INOUT) :: pdmfup(klon,klev), pentr(klon), plu(klon,klev), &
&                         pmfu(klon,klev), pmfub(klon), pmful(klon,klev), &
&                         pmfuq(klon,klev), pmfus(klon,klev), pmfuu(klon), &
&                         pmfuv(klon), &
!pmfuxt(klon,klev,ktrac),  &
&                         pqu(klon,klev), ptu(klon,klev), puu(klon,klev), &
&                         pvu(klon,klev)
!, pxtu(klon,klev,ktrac)
  INTEGER, INTENT (INOUT) :: kcbot(klon), klab(klon,klev), ktype(klon)

  !  Local scalars: 
  REAL(r8) :: zzzmb
  INTEGER :: jl, jt

  !  Local arrays: 
  LOGICAL :: llo3(klon)

  !  Intrinsic functions 
  INTRINSIC MAX, MIN


  !  Executable statements 

!-- 1. CALCULATE ENTRAINMENT AND DETRAINMENT RATES

!DIR$ IVDEP
!OCL NOVREC
  DO jl = 1, klon
    llo3(jl) = .FALSE.
    IF ( .NOT. ldcum(jl) .AND. klab(jl,kk+1)==0 .AND. &
&          pqen(jl,kk)>0.90*pqsen(jl,kk)) THEN
      llo3(jl) = .TRUE.
      ptu(jl,kk+1) = (cpd*pten(jl,kk)+pgeo(jl,kk)-pgeoh(jl,kk+1))*rcpd
      pqu(jl,kk+1) = pqen(jl,kk)
      plu(jl,kk+1) = 0.
      zzzmb = MAX(cmfcmin,-pverv(jl,kk)/g)
      zzzmb = MIN(zzzmb,cmfcmax)
      pmfub(jl) = zzzmb
      pmfu(jl,kk+1) = pmfub(jl)
      pmfus(jl,kk+1) = pmfub(jl)*(cpd*ptu(jl,kk+1)+pgeoh(jl,kk+1))
      pmfuq(jl,kk+1) = pmfub(jl)*pqu(jl,kk+1)
      pmful(jl,kk+1) = 0.
      pdmfup(jl,kk+1) = 0.
      kcbot(jl) = kk
      klab(jl,kk+1) = 1
      ktype(jl) = 3
      pentr(jl) = entrmid
      IF (lmfdudv) THEN
        puu(jl,kk+1) = puen(jl,kk)
        pvu(jl,kk+1) = pven(jl,kk)
        pmfuu(jl) = pmfub(jl)*puu(jl,kk+1)
        pmfuv(jl) = pmfub(jl)*pvu(jl,kk+1)
      END IF
    END IF
  END DO

!DIR$ IVDEP
!OCL NOVREC
! DO jt = 1, ktrac
!   DO jl = 1, klon
!     IF (llo3(jl)) THEN
!       pxtu(jl,kk+1,jt) = pxten(jl,kk,jt)
!       pmfuxt(jl,kk+1,jt) = pmfub(jl)*pxtu(jl,kk+1,jt)
!     END IF
!   END DO
! END DO

  RETURN
END SUBROUTINE cubasmc
