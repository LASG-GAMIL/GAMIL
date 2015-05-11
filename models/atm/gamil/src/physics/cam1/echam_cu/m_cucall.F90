!+ master routine - provides interface for: cumastr (cumulus parameterization)

MODULE m_cucall

CONTAINS

SUBROUTINE cucall(klon,klev,klevp1,klevm1, &
                  ptm1,pqm1,pum1,pvm1,pxm1,ptte,pqte, &
                  pvom,pvol,pxte,pverv,pqhfl,pxtec,   &
                  papp1,paphp1,pgeo, paprc,&
                  ktype,ldland,ptopmax,ptopmaxm, &
                  cnb,cmfmc,zlude,ztodt,mu2,md2,du2,eu2,ed2)

  ! ------------------> ldland included for AMIP2

  ! Description:
  !
  ! Master routine - provides interface for: cumastr (cumulus parameterization)
  !
  ! Method:
  !
  ! *cucall* - interface for *cumastr*:
  ! Provides input for cumastr.
  ! Receives updated tendencies, precipitation.
  !
  ! *cucall* is called from *physc*
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, December 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, December 1998, lookup tables removed
  ! 
  ! for more details see file AUTHORS
  !

  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE mo_cumulus_flux,   ONLY: cuparam
  USE time_manager,      ONLY: is_first_step,         &! initial run
                               is_first_restart_step   ! restart run 
  USE mo_constants,      ONLY: vtmpc1,  &! vtmpc1=rv/rd-1
                               inicon
  USE mo_convect_tables, ONLY: tlucua, jptlucu1, jptlucu2, &
                               lookuperror, lookupoverflow, &
                               set_lookup_tables
  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: klev, klevm1, klevp1, klon

  !  Array arguments with intent(In):
  REAL(r8), INTENT (IN) :: paphp1(klon,klevp1), papp1(klon,klev), pgeo(klon,klev), &
                       pqhfl(klon), pqm1(klon,klev), ptm1(klon,klev),        &
                       ptopmaxm(klon), pum1(klon,klev), pverv(klon,klev),     &
                       pvm1(klon,klev), pxm1(klon,klev), pxte(klon,klev),     &
                       ztodt
  LOGICAL, INTENT (IN) :: ldland(klon)    ! INCLUDED FOR AMIP2 !!

  !  Scalar arguments with intent(InOut):

  !  Array arguments with intent(InOut):
  REAL(r8), INTENT (INOUT) :: paprc(klon), &
   pqte(klon,klev),           &
                          ptopmax(klon),             &
                          ptte(klon,klev), pvol(klon,klev), pvom(klon,klev), &
                          pxtec(klon,klev),&
                          cmfmc(klon,klev),cnb(klon),mu2(klon,klev),md2(klon,klev),&
                          du2(klon,klev),eu2(klon,klev),ed2(klon,klev)
! INTEGER, INTENT (INOUT) :: ilab(klon,klev), ktype(klon)
  INTEGER, INTENT (INOUT) :: ktype(klon)

  !  Local scalars: 
  REAL(r8)    :: ztmst
  INTEGER :: ilevmin, jk, jl, jt
  INTEGER :: it

  !  Local arrays: 
  INTEGER :: ilab(klon,klev)
  REAL(r8) :: zlu(klon,klev), zlude(klon,klev), zmfd(klon,klev), zmfu(klon,klev), &
          zqp1(klon,klev), zqsat(klon,klev), zqu(klon,klev), zrain(klon), &
          ztp1(klon,klev), ztu(klon,klev), zup1(klon,klev), zvp1(klon,klev), &
          zxp1(klon,klev), prsfc(klon),pssfc(klon),paprs(klon)
  INTEGER :: icbot(klon), ictop(klon), itopec2(klon)
  LOGICAL :: locum(klon)

  !  External subroutines 
  EXTERNAL cumastr

  !  Intrinsic functions 
  INTRINSIC MIN

  !  Executable statements 

  lookupoverflow = .FALSE.

!-- 1. Calculate t,q and qs at main levels

  ztmst = ztodt
  IF (is_first_step())  ztmst = 0.5*ztodt
  IF (is_first_step().or.is_first_restart_step()) then
    CALL inicon()
    CALL cuparam()
    CALL set_lookup_tables
  END IF
  DO jk = 1, klev
    DO jl = 1, klon
      ztp1(jl,jk) = ptm1(jl,jk) + ptte(jl,jk)*ztmst
      zqp1(jl,jk) = pqm1(jl,jk) + pqte(jl,jk)*ztmst
      zxp1(jl,jk) = pxm1(jl,jk) + pxte(jl,jk)*ztmst
      zup1(jl,jk) = pum1(jl,jk) + pvom(jl,jk)*ztmst
      zvp1(jl,jk) = pvm1(jl,jk) + pvol(jl,jk)*ztmst
      it = INT(ztp1(jl,jk)*1000.)
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zqsat(jl,jk) = tlucua(it)/papp1(jl,jk)
      zqsat(jl,jk) = MIN(0.5,zqsat(jl,jk))
      zqsat(jl,jk) = zqsat(jl,jk)/(1.-vtmpc1*zqsat(jl,jk))
    END DO
    IF (lookupoverflow) CALL lookuperror ('cucall      ')

  END DO

  DO jl = 1, klon
    zrain(jl) = 0.
    locum(jl) = .FALSE.
  END DO
!-- 2. Call 'cumastr'(master-routine for cumulus parameterization)
!      -----------------------------------------------------------
  CALL cumastr(klon,klev,klevp1,klevm1,ilab,ztp1,zqp1,zxp1,zup1,          &
               zvp1,ldland,pverv,zqsat,pqhfl,papp1,paphp1,pgeo,ptte, &
               pqte,pvom,pvol,prsfc,pssfc,paprc,paprs,pxtec,locum,ktype,icbot,ictop, &
               ztu,zqu,zlu,zlude,zmfu,zmfd,zrain,ztodt,du2,eu2,ed2)

!-- 3. Pressure altitude of convective cloud tops
  ilevmin = klev - 4

  DO jl = 1, klon
    itopec2(jl) = klevp1
  END DO

  DO jk = 1, ilevmin
    DO jl = 1, klon
      IF (ilab(jl,jk)==2 .AND. itopec2(jl)==klevp1) THEN
        itopec2(jl) = jk
      END IF
    END DO
  END DO

  DO jl = 1, klon
    IF (itopec2(jl)==1) THEN
      ptopmax(jl) = papp1(jl,1)
    ELSE IF (itopec2(jl)/=klevp1) THEN
      ptopmax(jl) = paphp1(jl,itopec2(jl))
    ELSE
      ptopmax(jl) = 99999.
    END IF
    ptopmax(jl) = MIN(ptopmax(jl),ptopmaxm(jl))
  END DO

  cnb(:)=icbot(:)
! ptopmax(:)=itopec2(:)
  ptopmax(:)=ictop(:)
  cmfmc(:,:)=zmfu(:,:)+zmfd(:,:)
  mu2(:,:)=zmfu(:,:)
  md2(:,:)=zmfd(:,:)

  RETURN
END SUBROUTINE cucall

END MODULE m_cucall
