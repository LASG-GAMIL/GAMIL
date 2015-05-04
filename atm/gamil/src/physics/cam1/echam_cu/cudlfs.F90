!+ calculates level of free sinking for cumulus downdrafts and specifies
!  t,q,u and v values
!+ $Id: cudlfs.f90,v 1.5 1998/11/03 14:32:02 m214003 Exp $

SUBROUTINE cudlfs(klon,klev,klevp1,ptenh,pqenh,puen,pven,&
& pgeoh,paphp1,ptu,pqu,puu,pvu,ldcum,kcbot,kctop, &
&      pmfub,prfl,ptd,pqd,pud,pvd,pmfd,pmfds,pmfdq,pdmfdp,kdtop,lddraf)

  ! Description:
  !
  ! Calculates level of free sinking for cumulus downdrafts and specifies
  ! t,q,u and v values
  !
  ! Method:
  !
  ! Produce lfs-values for cumulus downdrafts for massflux cumulus
  ! parameterization
  !
  ! This routine is called from *cumastr*.
  ! Input are environmental values of t,q,u,v,p,phi and updraft values
  ! t,q,u and v and also cloud base massflux and cu-precipitation rate.
  ! It returns t,q,u and v values and massflux at lfs.
  !
  ! Check for negative buoyancy of air of equal parts of moist environmental 
  ! air and cloud air.
  !
  ! *cuadjtq* for calculating wet bulb t and q at lfs
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, in December 1986 and December 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE mo_constants,    ONLY: vtmpc1,  &! vtmpc1=rv/rd-1 
                             cpd       ! specific heat at constant pressure
  USE mo_cumulus_flux, ONLY: lmfdd,   &! true if cumulus downdr. is switched on
                             cmfdeps, &! fractional massflux for downdr. at lfs
                             lmfdudv   ! true if cumulus friction is switched on

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: klon, klev, klevp1

  !  Array arguments with intent(In):
  REAL(r8), INTENT (IN) :: paphp1(klon,klevp1), pgeoh(klon,klev), pmfub(klon), &
&                      pqenh(klon,klev), pqu(klon,klev), ptenh(klon,klev), &
&                      ptu(klon,klev), puen(klon,klev), puu(klon,klev), &
&                      pven(klon,klev), pvu(klon,klev)
  INTEGER, INTENT (IN) :: kcbot(klon), kctop(klon)
  LOGICAL, INTENT (IN) :: ldcum(klon)

  !  Array arguments with intent(InOut):
  REAL(r8), INTENT (INOUT) :: pdmfdp(klon,klev), pmfd(klon,klev), pmfdq(klon,klev), &
&                         pmfds(klon,klev), &
&                         pqd(klon,klev), prfl(klon), ptd(klon,klev), &
&                         pud(klon,klev), pvd(klon,klev)
  INTEGER, INTENT (INOUT) :: kdtop(klon)
  LOGICAL, INTENT (INOUT) :: lddraf(klon)

  !  Local scalars: 
  REAL(r8) :: zbuo, zmftop, zqtest, zttest
  INTEGER :: icall, ik, is, jk, jl, jt, ke

  !  Local arrays: 
  REAL(r8) :: zcond(klon), zph(klon), zqenwb(klon,klev), ztenwb(klon,klev)
  LOGICAL :: llo2(klon), llo3(klon)

  !  External subroutines 
  EXTERNAL cuadjtq

  !  Intrinsic functions 
  INTRINSIC MERGE


  !  Executable statements 

!-- 1. Set default values for downdrafts

  DO jl = 1, klon
    lddraf(jl) = .FALSE.
    kdtop(jl) = klevp1
  END DO

  IF ( .NOT. lmfdd) RETURN

!-- 2. Determine level of free sinking by doing a scan from top to base
!      of cumulus clouds for every point and proceed as follows:
!      (1) detemine wet bulb environmental t and q
!      (2) do mixing with cumulus cloud air
!      (3) check for negative buoyancy the assumption is that air of
!          downdrafts is mixture of 50% cloud air + 50% environmental air 
!          at wet bulb temperature (i.e. which became saturated due to
!          evaporation of rain and cloud water)

  ke = klev - 3

!-- 2.1  Calculate wet-bulb temperature and moisture
!        for environmental air in *cuadjtq*

  DO jk = 3, ke
    is = 0
    DO jl = 1, klon
      ztenwb(jl,jk) = ptenh(jl,jk)
      zqenwb(jl,jk) = pqenh(jl,jk)
      zph(jl) = paphp1(jl,jk)
      llo2(jl) = ldcum(jl) .AND. prfl(jl) > 0. .AND. .NOT. lddraf(jl) .AND. &
&          (jk<kcbot(jl) .AND. jk>kctop(jl))
      is = is + MERGE(1,0,llo2(jl))
    END DO
    IF (is==0) CYCLE

    ik = jk
    icall = 2
    CALL cuadjtq(klon,klev,ik,zph,ztenwb,zqenwb,llo2,icall)

!-- 2.2  Do mixing of cumulus and environmental air and check for negative
!        buoyancy.Then set values for downdraft at lfs.

!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      llo3(jl) = .FALSE.
      IF (llo2(jl)) THEN
        zttest = 0.5*(ptu(jl,jk)+ztenwb(jl,jk))
        zqtest = 0.5*(pqu(jl,jk)+zqenwb(jl,jk))
        zbuo = zttest*(1.+vtmpc1*zqtest) - ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl, &
&            jk))
        zcond(jl) = pqenh(jl,jk) - zqenwb(jl,jk)
        zmftop = -cmfdeps*pmfub(jl)
        IF (zbuo<0. .AND. prfl(jl)>10.*zmftop*zcond(jl)) THEN
          llo3(jl) = .TRUE.
          kdtop(jl) = jk
          lddraf(jl) = .TRUE.
          ptd(jl,jk) = zttest
          pqd(jl,jk) = zqtest
          pmfd(jl,jk) = zmftop
          pmfds(jl,jk) = pmfd(jl,jk)*(cpd*ptd(jl,jk)+pgeoh(jl,jk))
          pmfdq(jl,jk) = pmfd(jl,jk)*pqd(jl,jk)
          pdmfdp(jl,jk-1) = -0.5*pmfd(jl,jk)*zcond(jl)
          prfl(jl) = prfl(jl) + pdmfdp(jl,jk-1)
        END IF
      END IF
    END DO

    IF (lmfdudv) THEN
      DO jl = 1, klon
        IF (pmfd(jl,jk)<0.) THEN
          pud(jl,jk) = 0.5*(puu(jl,jk)+puen(jl,jk-1))
          pvd(jl,jk) = 0.5*(pvu(jl,jk)+pven(jl,jk-1))
        END IF
      END DO

    END IF

  END DO

  RETURN
END SUBROUTINE cudlfs
