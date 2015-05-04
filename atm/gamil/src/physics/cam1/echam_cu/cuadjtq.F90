!+ produce t,q and l values for cloud ascent
!+ $Id: cuadjtq.f90,v 1.14 2000/09/14 13:14:12 m214003 Exp $

SUBROUTINE cuadjtq(klon,klev,kk,pp,pt,pq,ldflag,kcall)

  ! Description:
  !
  ! Produce t,q and l values for cloud ascent
  !
  ! Method:
  !
  ! This routine is called from subroutines:
  !    *cubase*   (t and q at condenstion level)
  !    *cuasc*    (t and q at cloud levels)
  !    *cuini*    (environmental t and qs values at half levels)
  ! input are unadjusted t and q values,
  ! it returns adjusted values of t and q
  ! note: input parameter kcall defines calculation as
  !    kcall=0    env. t and qs in*cuini*
  !    kcall=1    condensation in updrafts  (e.g. cubase, cuasc)
  !    kcall=2    evaporation in downdrafts (e.g. cudlfs,cuddraf)
  !
  ! Externals:
  ! 3 lookup tables ( tlucua, tlucub, tlucuc )
  ! for condensation calculations.
  ! the tables are initialised in *setphys*.
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, December 1989, original source
  ! D. Salmond, CRAY(UK), August 1991, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, December 1998, lookup tables removed
  ! 
  ! for more details see file AUTHORS
  !

  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE mo_constants,      ONLY: vtmpc1
  USE mo_convect_tables, ONLY: tlucua,   & ! table a
                               tlucub,   & ! table b
                               tlucuc,   & ! table c
                               jptlucu1, jptlucu2, &
                               lookuperror, lookupoverflow

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: kcall, kk, klev, klon

  !  Array arguments with intent(In):
  REAL(r8), INTENT (IN) :: pp(klon)
  LOGICAL, INTENT (IN) :: ldflag(klon)

  !  Array arguments with intent(InOut):
  REAL(r8), INTENT (INOUT) :: pq(klon,klev), pt(klon,klev)

  !  Local scalars: 
  REAL(r8)    :: zcond1, zcor, zqsat
  INTEGER :: isum, jl
  INTEGER :: it

  !  Local arrays: 
  REAL(r8) :: zcond(klon), zqp(klon)

  !  Intrinsic functions 
  INTRINSIC MAX, MIN
  INTRINSIC INT

  !  Executable statements 

  ! 1.  Calculate condensation and adjust t and q accordingly

  lookupoverflow = .FALSE.

  zcond(:) = 0.

  IF (kcall==1) THEN

    isum = 0
!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      IF (ldflag(jl)) THEN
        zqp(jl) = 1./pp(jl)
        it = INT(pt(jl,kk)*1000.)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsat = tlucua(it)*zqp(jl)
        zqsat = MIN(0.5,zqsat)
        zcor  = 1./(1.-vtmpc1*zqsat)
        zqsat = zqsat*zcor
        zcond(jl) = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
        zcond(jl) = MAX(zcond(jl),0.)
        pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond(jl)
        pq(jl,kk) = pq(jl,kk) - zcond(jl)
        IF (ABS(zcond(jl)) > 0.) isum = isum + 1
      END IF
    END DO
    IF (lookupoverflow) CALL lookuperror ('cuadjtq (1) ')

    IF (isum/=0) THEN
!DIR$ IVDEP
!OCL NOVREC
      DO jl = 1, klon
        IF (ldflag(jl) .AND. ABS(zcond(jl)) > 0.) THEN
          it = INT(pt(jl,kk)*1000.)
          IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
          it = MAX(MIN(it,jptlucu2),jptlucu1)
          zqsat = tlucua(it)*zqp(jl)
          zqsat = MIN(0.5,zqsat)
          zcor = 1./(1.-vtmpc1*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
          pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond1
          pq(jl,kk) = pq(jl,kk) - zcond1
        END IF
      END DO
      IF (lookupoverflow) CALL lookuperror ('cuadjtq (2) ')
    END IF

  END IF

  IF (kcall==2) THEN

    isum = 0
!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      IF (ldflag(jl)) THEN
        zqp(jl) = 1./pp(jl)
        it = INT(pt(jl,kk)*1000.)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsat = tlucua(it)*zqp(jl)
        zqsat = MIN(0.5,zqsat)
        zcor = 1./(1.-vtmpc1*zqsat)
        zqsat = zqsat*zcor
        zcond(jl) = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
        zcond(jl) = MIN(zcond(jl),0.)
        pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond(jl)
        pq(jl,kk) = pq(jl,kk) - zcond(jl)
        IF (ABS(zcond(jl)) > 0.) isum = isum + 1
      END IF
    END DO
    IF (lookupoverflow) CALL lookuperror ('cuadjtq (3) ')

    IF (isum/=0) THEN
!DIR$ IVDEP
!OCL NOVREC
      DO jl = 1, klon
        IF (ldflag(jl) .AND. ABS(zcond(jl)) > 0.) THEN
          it = INT(pt(jl,kk)*1000.)
          IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
          it = MAX(MIN(it,jptlucu2),jptlucu1)
          zqsat = tlucua(it)*zqp(jl)
          zqsat = MIN(0.5,zqsat)
          zcor = 1./(1.-vtmpc1*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
          pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond1
          pq(jl,kk) = pq(jl,kk) - zcond1
        END IF
      END DO
      IF (lookupoverflow) CALL lookuperror ('cuadjtq (4) ')

    END IF

  END IF

  IF (kcall==0) THEN

    isum = 0
!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      zqp(jl) = 1./pp(jl)
      it = INT(pt(jl,kk)*1000.)
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zqsat = tlucua(it)*zqp(jl)
      zqsat = MIN(0.5,zqsat)
      zcor = 1./(1.-vtmpc1*zqsat)
      zqsat = zqsat*zcor
      zcond(jl) = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
      pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond(jl)
      pq(jl,kk) = pq(jl,kk) - zcond(jl)
      IF (ABS(zcond(jl)) > 0.) isum = isum + 1
    END DO
    IF (lookupoverflow) CALL lookuperror ('cuadjtq (5) ')

    IF (isum/=0) THEN
!DIR$ IVDEP
!OCL NOVREC
      DO jl = 1, klon
        it = INT(pt(jl,kk)*1000.)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsat = tlucua(it)*zqp(jl)
        zqsat = MIN(0.5,zqsat)
        zcor = 1./(1.-vtmpc1*zqsat)
        zqsat = zqsat*zcor
        zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
        pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond1
        pq(jl,kk) = pq(jl,kk) - zcond1
      END DO
      IF (lookupoverflow) CALL lookuperror ('cuadjtq (6) ')

    END IF

  END IF

  IF (kcall==4) THEN

!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      zqp(jl) = 1./pp(jl)
      it = INT(pt(jl,kk)*1000.)
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zqsat = tlucua(it)*zqp(jl)
      zqsat = MIN(0.5,zqsat)
      zcor = 1./(1.-vtmpc1*zqsat)
      zqsat = zqsat*zcor
      zcond(jl) = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
      pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond(jl)
      pq(jl,kk) = pq(jl,kk) - zcond(jl)
    END DO
    IF (lookupoverflow) CALL lookuperror ('cuadjtq (7) ')

!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      it = INT(pt(jl,kk)*1000.)
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zqsat = tlucua(it)*zqp(jl)
      zqsat = MIN(0.5,zqsat)
      zcor = 1./(1.-vtmpc1*zqsat)
      zqsat = zqsat*zcor
      zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
      pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond1
      pq(jl,kk) = pq(jl,kk) - zcond1
    END DO

    IF (lookupoverflow) CALL lookuperror ('cuadjtq (8) ')

  END IF

  RETURN
END SUBROUTINE cuadjtq
