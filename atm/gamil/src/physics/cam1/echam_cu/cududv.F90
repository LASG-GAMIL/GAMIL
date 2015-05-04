!+ updates u and v tendencies, does global diagnostic of dissipation
!+ $Id: cududv.f90,v 1.3 1998/10/28 12:28:40 m214003 Exp $

SUBROUTINE cududv(klon,klev,klevp1,ktopm2,ktype,kcbot,paphp1, &
&      ldcum,puen,pven,pvom,pvol,puu,pud,pvu,pvd,pmfu,pmfd)

  ! Description:
  !
  ! Updates u and v tendencies, does global diagnostic of dissipation
  !
  ! Method:
  !
  ! *cududv* is called from *cumastr*
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, in July 1986, and December 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  !
  ! for more details see file AUTHORS
  !

  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE mo_constants, ONLY: g ! gravity acceleration

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: klon, klev, klevp1, ktopm2

  !  Array arguments with intent(In):
  REAL(r8), INTENT (IN) :: paphp1(klon,klevp1), pmfd(klon,klev), pmfu(klon,klev), &
&           pud(klon,klev), puen(klon,klev), puu(klon,klev), pvd(klon,klev), &
&           pven(klon,klev), pvu(klon,klev)
  INTEGER, INTENT (IN) :: kcbot(klon), ktype(klon)
  LOGICAL(r8), INTENT (IN) :: ldcum(klon)

  !  Scalar arguments with intent(InOut):

  !  Array arguments with intent(InOut):
  REAL(r8), INTENT (INOUT) :: pvol(klon,klev), pvom(klon,klev)

  !  Local scalars: 
  REAL(r8) :: zdudt, zdvdt, zsum, zzp
  INTEGER :: ik, ikb, jk, jl

  !  Local arrays: 
  REAL(r8) :: zdiss(klon), zmfdu(klon,klev), zmfdv(klon,klev), zmfuu(klon,klev), &
&      zmfuv(klon,klev)

  !  Intrinsic functions 
  INTRINSIC MERGE, SUM


  !  Executable statements 

!-- 1. Calculate fluxes and update u and v tendencies

  DO jk = ktopm2, klev
    ik = jk - 1
    DO jl = 1, klon
      IF (ldcum(jl)) THEN
        zmfuu(jl,jk) = pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
        zmfuv(jl,jk) = pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
        zmfdu(jl,jk) = pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
        zmfdv(jl,jk) = pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
      END IF
    END DO
  END DO
  DO jk = ktopm2, klev
!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      IF (ldcum(jl) .AND. jk>kcbot(jl)) THEN
        ikb = kcbot(jl)
        zzp = ((paphp1(jl,klevp1)-paphp1(jl,jk))/(paphp1(jl, &
&            klevp1)-paphp1(jl,ikb)))
        zzp = MERGE(zzp**2,zzp,ktype(jl)==3)
        zmfuu(jl,jk) = zmfuu(jl,ikb)*zzp
        zmfuv(jl,jk) = zmfuv(jl,ikb)*zzp
        zmfdu(jl,jk) = zmfdu(jl,ikb)*zzp
        zmfdv(jl,jk) = zmfdv(jl,ikb)*zzp
      END IF
    END DO
  END DO

  DO jl = 1, klon
    zdiss(jl) = 0.
  END DO

  DO jk = ktopm2, klev

    IF (jk<klev) THEN
      DO jl = 1, klon
        IF (ldcum(jl)) THEN
          zdudt = (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*(zmfuu(jl,jk+1)-zmfuu( &
&              jl,jk)+zmfdu(jl,jk+1)-zmfdu(jl,jk))
          zdvdt = (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*(zmfuv(jl,jk+1)-zmfuv( &
&              jl,jk)+zmfdv(jl,jk+1)-zmfdv(jl,jk))
          zdiss(jl) = zdiss(jl) + puen(jl,jk)*(zmfuu(jl,jk+1)-zmfuu(jl,jk)+ &
&              zmfdu(jl,jk+1)-zmfdu(jl,jk)) + pven(jl,jk)*(zmfuv(jl,jk+1)- &
&              zmfuv(jl,jk)+zmfdv(jl,jk+1)-zmfdv(jl,jk))
          pvom(jl,jk) = pvom(jl,jk) + zdudt
          pvol(jl,jk) = pvol(jl,jk) + zdvdt
        END IF
      END DO

    ELSE
      DO jl = 1, klon
        IF (ldcum(jl)) THEN
          zdudt = -(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*(zmfuu(jl,jk)+zmfdu(jl &
&              ,jk))
          zdvdt = -(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*(zmfuv(jl,jk)+zmfdv(jl &
&              ,jk))
          zdiss(jl) = zdiss(jl) - (puen(jl,jk)*(zmfuu(jl,jk)+zmfdu(jl, &
&              jk))+pven(jl,jk)*(zmfuv(jl,jk)+zmfdv(jl,jk)))
          pvom(jl,jk) = pvom(jl,jk) + zdudt
          pvol(jl,jk) = pvol(jl,jk) + zdvdt
        END IF
      END DO

    END IF
  END DO

  zsum = SUM(zdiss(1:klon))

  RETURN
END SUBROUTINE cududv
