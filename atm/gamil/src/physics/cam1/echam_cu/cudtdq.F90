!+ updates t and q tendencies, precipitation rates does global diagnostics
!+ $Id: cudtdq.f90,v 1.4 1998/10/28 12:28:37 m214003 Exp $

SUBROUTINE cudtdq(klon,klev,klevp1,ktopm2,paphp1,ldcum,pten,ptte,pqte, &
&      pxtec,pmfus,pmfds,pmfuq,pmfdq,pmful,pdmfup, &
&      pdmfdp,plude,pdpmel,prain,prfl,psfl,prsfc, &
&      pssfc,paprc,paprs,ztodt)

  ! Description:
  !
  ! Updates t and q tendencies, precipitation rates does global diagnostics
  !
  ! Method:
  !
  ! *cudtdq* is called from *cumastr*
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, in July 1986 and December 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  !
  ! for more details see file AUTHORS
  !

  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE mo_constants, ONLY: rhoh2o, &! density of liquid water
                          tmelt,  &! temperature of fusion of ice
                          alv,    &! latent heat for vaporisation
                          g,      &! gravity acceleration
                          alf,    &! latent heat for fusion
                          als,    &! latent heat for sublimation
                          rcpd     ! rcpd=1./cpd

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: klev, klevp1, klon, ktopm2

  !  Array arguments with intent(In):
  REAL(r8), INTENT (IN) :: paphp1(klon,klevp1), pdmfdp(klon,klev), &
&      pdmfup(klon,klev), pdpmel(klon,klev), plude(klon,klev), &
&      pmfdq(klon,klev), pmfds(klon,klev), &
&      pmful(klon,klev), pmfuq(klon,klev), pmfus(klon,klev), &
&       prain(klon), prfl(klon), &
&      psfl(klon), pten(klon,klev),ztodt
  LOGICAL, INTENT (IN) :: ldcum(klon)

  !  Scalar arguments with intent(InOut):

  !  Array arguments with intent(InOut):
  REAL(r8), INTENT (INOUT) :: paprc(klon), paprs(klon), pqte(klon,klev), &
&                         prsfc(klon), pssfc(klon), ptte(klon,klev), &
&                         pxtec(klon,klev)

  !  Local scalars: 
  REAL(r8) :: zalv, zdiagt, zdiagw, zdqdt, zdtdt, zdxtdt
  INTEGER :: jk, jl, jt
  LOGICAL :: llo1

  !  Intrinsic functions 
  INTRINSIC MERGE


  !  Executable statements 

!-- 1. Specify parameters

! zdiagt = 0.5*ztodt
  zdiagt = ztodt
  zdiagw = zdiagt/rhoh2o

!-- 2. Incrementation of t and q tendencies

  DO jk = ktopm2, klev

    IF (jk<klev) THEN
      DO jl = 1, klon
        IF (ldcum(jl)) THEN
          llo1 = (pten(jl,jk)-tmelt) > 0.
          zalv = MERGE(alv,als,llo1)
          zdtdt = (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*rcpd*(pmfus(jl,jk+1)- &
&              pmfus(jl,jk)+pmfds(jl,jk+1)-pmfds(jl,jk)-alf*pdpmel(jl,jk)- &
&              zalv*(pmful(jl,jk+1)-pmful(jl,jk)-plude(jl, &
&              jk)-(pdmfup(jl,jk)+pdmfdp(jl,jk))))
          ptte(jl,jk) = ptte(jl,jk) + zdtdt
          zdqdt = (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*(pmfuq(jl,jk+1)-pmfuq( &
&              jl,jk)+pmfdq(jl,jk+1)-pmfdq(jl,jk)+pmful(jl,jk+1)-pmful(jl,jk)- &
&              plude(jl,jk)-(pdmfup(jl,jk)+pdmfdp(jl,jk)))
          pqte(jl,jk) = pqte(jl,jk) + zdqdt
          pxtec(jl,jk) = (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*plude(jl,jk)
        END IF
      END DO

    ELSE
      DO jl = 1, klon
        IF (ldcum(jl)) THEN
          llo1 = (pten(jl,jk)-tmelt) > 0.
          zalv = MERGE(alv,als,llo1)
          zdtdt = -(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*rcpd*(pmfus(jl,jk)+ &
&              pmfds(jl,jk)+alf*pdpmel(jl,jk)-zalv*(pmful(jl,jk)+pdmfup(jl, &
&              jk)+pdmfdp(jl,jk)+plude(jl,jk)))
          ptte(jl,jk) = ptte(jl,jk) + zdtdt
          zdqdt = -(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*(pmfuq(jl,jk)+pmfdq(jl &
&              ,jk)+plude(jl,jk)+(pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
          pqte(jl,jk) = pqte(jl,jk) + zdqdt
          pxtec(jl,jk) = (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*plude(jl,jk)
        END IF
      END DO

    END IF
  END DO
!-- 3. Update surface fields and do global budgets
  DO jl = 1, klon
    prsfc(jl) = prfl(jl)
    pssfc(jl) = psfl(jl)
    paprc(jl) = paprc(jl) + zdiagw*(prfl(jl)+psfl(jl))
    paprs(jl) = paprs(jl) + zdiagw*psfl(jl)
  END DO

  RETURN
END SUBROUTINE cudtdq
