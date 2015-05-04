!+ master routine for cumulus massfluxscheme

SUBROUTINE cumastr(klon,klev,klevp1,klevm1,ilab,pten,pqen,pxen, puen,pven,   &
                   ldland,pverv,pqsen,pqhfl,papp1,paphp1,pgeo,   &
                   ptte,pqte,pvom,pvol,prsfc,pssfc,paprc,paprs,pxtec,ldcum,ktype,kcbot, &
                   kctop,ptu,pqu,plu,plude,pmfu,pmfd,prain,&
                   ztodt,du,eu,ed)

  ! -------------> ldland included for AMIP2 !!

  ! Description:
  !
  ! Master routine for cumulus massfluxscheme
  !
  ! Method:
  !
  ! This routine computes the physical tendencies of the
  ! prognostic variables t,q,u and v due to convective processes.
  ! Processes considered are: 
  !    Convective fluxes, 
  !    formation of precipitation,
  !    evaporation of falling rain below cloud base,
  !    saturated cumulus downdrafts.
  !
  ! *cumastr* is called from *cucall*.
  ! The routine takes its input from the long-term storage
  ! t,q,u,v,phi and p and moisture tendencies.
  ! It returns its output to the same space
  !   1. modified tendencies of model variables
  !   2. rates of convective precipitation (used in subroutine surf)
  !
  ! Parameterization is done using a massflux-scheme.
  ! (1) Define constants and parameters.
  ! (2) Specify values (t,q,qs.) at half levels and
  !     initialize updraft- and downdraft-values in 'cuini'.
  ! (3) Calculate cloud base in 'cubase'.
  !     and specify cloud base massflux from pbl moisture budget
  ! (4) Do cloud ascent in 'cuasc' in absence of downdrafts.
  ! (5) Do downdraft calculations:
  !     (a) Determine values at lfs in 'cudlfs'.
  !     (b) Determine moist descent in 'cuddraf'.
  !     (c) Decalculate cloud base massflux considering the
  !         effect of cu-downdrafts.
  ! (6) Do final cloud ascent in 'cuasc'.
  ! (7) Do final adjusments to convective fluxes in 'cuflx',
  !     do evaporation in subcloud layer.
  ! (8) Calculate increments of t and q in 'cudtdq'.
  ! (9) Calculate increments of u and v in 'cududv'.
  !
  ! Externals:
  !
  ! cuini:  initializes values at vertical grid used in cu-parametr.
  ! cubase: cloud base calculation for penetr.and shallow convection
  ! cuasc:  cloud ascent for entraining plume
  ! cudlfs: determines values at lfs for downdrafts
  ! cuddraf:does moist descent for cumulus downdrafts
  ! cuflx:  final adjustments to convective fluxes (also in pbl)
  ! cudqdt: updates tendencies for t and q
  ! cududv: updates tendencies for u and v
  !
  ! Switches:
  ! lmfpen=.t.   penetrative convection is switched on
  ! lmfscv=.t.   shallow convection is switched on
  ! lmfmid=.t.   midlevel convection is switched on
  ! lmfdd=.t.    cumulus downdrafts switched on
  ! lmfdudv=.t.  cumulus friction switched on
  !
  ! model parameters (defined in subroutine cuparam):
  ! entrpen    entrainment rate for penetrative convection
  ! entrscv    entrainment rate for shallow convection
  ! entrmid    entrainment rate for midlevel convection
  ! entrdd     entrainment rate for cumulus downdrafts
  ! cmfctop    relative cloud massflux at level above nonbuoyancy leve
  ! cmfcmax    maximum massflux value allowed for
  ! cmfcmin    minimum massflux value (for safety)
  ! cmfdeps    fractional massflux for downdrafts at lfs
  ! cprcon     coefficient for conversion from cloud water to rain
  !
  ! Paper on massflux scheme (tiedtke,1989)
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, in 1986, 1987, and 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE time_manager,     only: is_first_step          ! initial run
  USE mo_constants,     only: rd,      &! gas constant for dry air
                              c4les,   &!
                              c5les,   &!
                              cpd,     &! specific heat at constant pressure
                              g,       &! gravity acceleration
                              vtmpc1,  &! vtmpc1=rv/rd-1
                              alv       ! latent heat for vaporisation
  USE mo_cumulus_flux,  only: lmfdudv, &! true if cum. friction is switched on
                              cmfdeps, &! fractional massflux for downdr. at lfs
                              lmfdd,   &! true if cum. downdraft is switched on
                              entrpen, &! entrainment rate for penetrative conv.
                              entrscv   ! entrainment rate for shallow conv.

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: klev, klevm1, klevp1, klon

  !  Array arguments with intent(In):
  REAL(r8), INTENT (IN) :: paphp1(klon,klevp1), papp1(klon,klev), pgeo(klon,klev),  &
                       pqen(klon,klev), pqhfl(klon), pqsen(klon,klev),          &
                       pten(klon,klev), puen(klon,klev),                       &
                       pven(klon,klev), pverv(klon,klev),                      &
                       pxen(klon,klev),      &
                       ztodt
!llj(2005-07-25)
  REAL(r8)             :: dept,rhcr0,zrelh(klon)  !(llj05-07-25)
  INTEGER              :: ktop0         !(llj05-07-25)
  LOGICAL, INTENT (IN) :: ldland(klon)  !  included for AMIP2 !!

  !  Scalar arguments with intent(InOut):

  !  Array arguments with intent(InOut):
  REAL(r8), INTENT (INOUT) :: paprc(klon), paprs(klon), plu(klon,klev),            &
                          plude(klon,klev), pmfd(klon,klev), pmfu(klon,klev),  &
                          pqte(klon,klev), pqu(klon,klev), prain(klon),        &
                          prsfc(klon), pssfc(klon), ptte(klon,klev),           &
                          ptu(klon,klev), pvol(klon,klev), pvom(klon,klev),  &
                          pxtec(klon,klev),du(klon,klev),eu(klon,klev),        &
                          ed(klon,klev)
  INTEGER, INTENT (INOUT) :: ilab(klon,klev), kcbot(klon), kctop(klon),  &
                             ktype(klon)
  LOGICAL, INTENT (INOUT) :: ldcum(klon)

  !  Local scalars: 
  REAL(r8) :: zalvdcp, zb, zbi, zcons2, zdepth, zdhdz, zdqmin, zdz, zeps, zfac,  &
          zgam, zhhat, zhsat, zmfmax, zpbmpt, zqalv, zqumqe, zrh, zro, ztau, &
          ztmst, zzz
  INTEGER :: icum, ikb, itopm2, jk, jl, jt
  LOGICAL :: llo1

  !  Local arrays: 
  REAL(r8) :: zcape(klon), zdmfdp(klon,klev), zdmfup(klon,klev),                    &
          zdpmel(klon,klev), zdqcv(klon), zdqpbl(klon), zentr(klon),            &
          zgeoh(klon,klev), zhcbase(klon), zheat(klon), zhhatt(klon,klev),      &
          zhmin(klon), zmfdq(klon,klev), zmfds(klon,klev),                      &
          zmfub(klon), zmfub1(klon), zmful(klon,klev), &
          zmfuq(klon,klev), zmfus(klon,klev), &
          zqd(klon,klev), zqenh(klon,klev), zqsenh(klon,klev), zrfl(klon),      &
          zsfl(klon), ztd(klon,klev), ztenh(klon,klev), zud(klon,klev),         &
          zuu(klon,klev), zvd(klon,klev), zvu(klon,klev), zxenh(klon,klev)
  INTEGER :: ictop0(klon), idtop(klon), ihmin(klon), ilwmin(klon)
  LOGICAL :: loddraf(klon)


  !  External subroutines 
  EXTERNAL cuasc, cubase, cuddraf, cudlfs, cudtdq, cududv, cuflx, cuini

  !  Intrinsic functions 
  INTRINSIC ABS, MAX, MERGE, MIN, SQRT


  !  Executable statements 

!-- 1. Specify constants and parameters

  ztmst = ztodt
! IF (nstep==nstart) ztmst = 0.5*ztodt
  IF (is_first_step()) ztmst = 0.5*ztodt
  zcons2 = 1./(g*ztmst)

!-- 2. Initialize values at vertical grid points in 'cuini'

  CALL cuini(klon,klev,klevp1,klevm1,pten,pqen,pqsen,pxen,puen,          &
             pven, &
     pverv,pgeo,paphp1,     &
             zgeoh,ztenh,zqenh,zqsenh,zxenh,ilwmin,ptu,pqu,ztd,zqd,zuu,zvu,zud,zvd, &
             pmfu,pmfd,zmfus,zmfds,zmfuq,zmfdq,zdmfup,zdmfdp,zdpmel,plu,plude,ilab)
!-- 3. Cloud base calculations

  ! (A) Determine cloud base values in 'cubase'

  CALL cubase(klon,klev,klevp1,klevm1,ztenh,zqenh,zgeoh,paphp1, &
              ptu,pqu,plu,puen,pven,zuu,zvu,ldcum,kcbot,ilab)

  ! (B) Determine total moisture convergence and
  !     then decide on type of cumulus convection

  jk = 1
  DO jl = 1, klon
    zdqcv(jl) = pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
    zdqpbl(jl) = 0.0
    idtop(jl) = 0
  END DO
  DO jk = 2, klev
    DO jl = 1, klon
      zdqcv(jl) = zdqcv(jl) + pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
      IF (jk>=kcbot(jl)) zdqpbl(jl) = zdqpbl(jl) + pqte(jl,jk)*(paphp1(jl,jk+1) &
             -paphp1(jl,jk))
    END DO
  END DO
  DO JL = 1, klon
!   if(ldland(jl))then
!     ktype(jl) = MERGE(1,2,zdqcv(jl)>MAX(0.,-1.05*pqhfl(jl)*g))
!   else
!     ktype(jl) = MERGE(1,2,zdqcv(jl)>MAX(0.,-1.3*pqhfl(jl)*g))
!   endif
    ktype(jl) = MERGE(1,2,zdqcv(jl)>MAX(0.,-1.1*pqhfl(jl)*g))
!   ktype(jl) = MERGE(1,2,zdqcv(jl)>MAX(0.,-1.1*pqhfl(jl)*g))
!   ktype(jl) = MERGE(2,0,pqhfl(jl)<0.)
!   IF (pqhfl(jl)<0.)then
!     ktype(jl) = MERGE(1,2,(zdqcv(jl)>-1.1*pqhfl(jl)*g) )
!     ktype(jl) = MERGE(1,2,(zdqcv(jl)>-0.5*pqhfl(jl)*g) )
!   END IF
  END DO
  ! (C) Determine moisture supply for boundary layer
  !     and determine cloud base massflux ignoring
  !     the effects of downdrafts at this stage

!DIR$ IVDEP
!OCL NOVREC
  DO jl = 1, klon
    ikb = kcbot(jl)
    zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zqenh(jl,ikb)
    zdqmin = MAX(0.01*zqenh(jl,ikb),1.E-10)
    llo1 = zdqpbl(jl) > 0. .AND. zqumqe > zdqmin .AND. ldcum(jl)
    zmfub(jl) = MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),0.01,llo1)
    zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
    zmfub(jl) = MIN(zmfub(jl),zmfmax)
    IF ( .NOT. llo1) ldcum(jl) = .FALSE.
    zcape(jl) = 0.
    zheat(jl) = 0.
    zentr(jl) = MERGE(entrpen,entrscv,ktype(jl)==1)
  END DO

!-- 4. Determine cloud ascent for entraining plume

  ! (a) Estimate cloud height for entrainment/detrainment
  !     calculations in cuasc (max.possible cloud height
  !     for non-entraining plume, following a.-s.,1974)

!DIR$ IVDEP
!OCL NOVREC
  DO jl = 1, klon
    ikb = kcbot(jl)
    zhcbase(jl) = cpd*ptu(jl,ikb) + zgeoh(jl,ikb) + alv*pqu(jl,ikb)
    ictop0(jl) = kcbot(jl) - 1
  END DO
  zalvdcp = alv/cpd
  zqalv = 1./alv
  DO jk = klevm1, 3, -1
!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      zhsat = cpd*ztenh(jl,jk) + zgeoh(jl,jk) + alv*zqsenh(jl,jk)
      zgam = c5les*zalvdcp*zqsenh(jl,jk)/((1.-vtmpc1*zqsenh(jl,jk)) &
          *(ztenh(jl,jk)-c4les)**2)
      zzz = cpd*ztenh(jl,jk)*0.608
      zhhat = zhsat - (zzz+zgam*zzz)/(1.+zgam*zzz*zqalv)*MAX(zqsenh(jl,jk)- &
              zqenh(jl,jk),0.)
      zhhatt(jl,jk) = zhhat
      IF (jk<ictop0(jl) .AND. zhcbase(jl)>zhhat) ictop0(jl) = jk
    END DO
  END DO

  DO jl = 1, klon
    jk = kcbot(jl)
    zhsat = cpd*ztenh(jl,jk) + zgeoh(jl,jk) + alv*zqsenh(jl,jk)
    zgam = c5les*zalvdcp*zqsenh(jl,jk)/((1.-vtmpc1*zqsenh(jl,jk))* &
           (ztenh(jl,jk)-c4les)**2)
    zzz = cpd*ztenh(jl,jk)*0.608
    zhhat = zhsat - (zzz+zgam*zzz)/(1.+zgam*zzz*zqalv)*MAX(zqsenh(jl,jk)- &
            zqenh(jl,jk),0.)
    zhhatt(jl,jk) = zhhat
  END DO

  ! Find lowest possible org. detrainment level

  DO jl = 1, klon
    llo1 = ldcum(jl) .AND. ktype(jl) == 1
    IF (llo1) THEN
      ikb = kcbot(jl)
      zhmin(jl) = 0.
      ihmin(jl) = ikb
    ELSE
      ihmin(jl) = -1
    END IF
  END DO

  zb = 25.
  zbi = 1./(zb*g)
  DO jk = klev, 1, -1
    DO jl = 1, klon
      llo1 = ldcum(jl) .AND. ktype(jl) == 1 .AND. ihmin(jl) == kcbot(jl)
      IF (llo1 .AND. jk<kcbot(jl) .AND. jk>=ictop0(jl)) THEN
        ikb = kcbot(jl)
        zro = paphp1(jl,jk)/(rd*ztenh(jl,jk))
        zdz = (paphp1(jl,jk)-paphp1(jl,jk-1))/(g*zro)
        zdhdz = (cpd*(pten(jl,jk-1)-pten(jl,jk))+alv*(pqen(jl,jk-1)-pqen(jl,jk)) &
              +(pgeo(jl,jk-1)-pgeo(jl,jk)))*g/(pgeo(jl,jk-1)-pgeo(jl,jk))
        zdepth = zgeoh(jl,jk) - zgeoh(jl,ikb)
        zfac = SQRT(1.+zdepth*zbi)
        zhmin(jl) = zhmin(jl) + zdhdz*zfac*zdz
        zrh = -alv*(zqsenh(jl,jk)-zqenh(jl,jk))*zfac
        IF (zhmin(jl)>zrh) ihmin(jl) = jk
      END IF
    END DO
  END DO

  DO jl = 1, klon
    IF (ldcum(jl) .AND. ktype(jl)==1) THEN
      IF (ihmin(jl)<ictop0(jl)) ihmin(jl) = ictop0(jl)
    END IF
  END DO

  ! (B) Do ascent in 'cuasc'in absence of downdrafts

  CALL cuasc(klon,klev,klevp1,klevm1,ztenh,zqenh,zxenh,puen,pven,            &
     pten,pqen,pqsen,pgeo,zgeoh,papp1,paphp1,    &
             pverv,ilwmin,ldcum,ldland,ktype,ilab,ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr, &
             zmfus,zmfuq,zmful,plude,zdmfup,ihmin,zhhatt,zhcbase,zqsenh,kcbot,kctop,    &
             ictop0,icum,ztodt,eu,du)

  ! ----------> ldland included for AMIP2 !!

  ! (C) Check cloud depth and change entrainment rate accordingly
  !     calculate precipitation rate (for downdraft calculation)

!DIR$ IVDEP
!OCL NOVREC
  DO jl = 1, klon
    zpbmpt = paphp1(jl,kcbot(jl)) - paphp1(jl,kctop(jl))
    IF (ldcum(jl) .AND. ktype(jl)==1 .AND. zpbmpt<2.E4) ktype(jl) = 2

    IF (ktype(jl)==2) zentr(jl) = entrscv
    zrfl(jl) = zdmfup(jl,1)
  END DO
  DO jk = 2, klev
    DO jl = 1, klon
      zrfl(jl) = zrfl(jl) + zdmfup(jl,jk)
    END DO
  END DO

!-- 5. Cumulus downdraft calculations

  IF (lmfdd) THEN

    ! (A) Determine lfs in 'cudlfs'

    CALL cudlfs(klon,klev,klevp1,ztenh,zqenh,puen,pven,&
   zgeoh,paphp1,ptu,pqu,zuu,zvu,ldcum,kcbot, &
                kctop,zmfub,zrfl,ztd,zqd,zud,zvd,pmfd,zmfds,zmfdq,zdmfdp,idtop,   &
               loddraf)

    ! (B)  Determine downdraft t,q and fluxes in 'cuddraf'

    CALL cuddraf(klon,klev,klevp1,ztenh,zqenh,puen,pven,&
    zgeoh,paphp1,zrfl,ztd,zqd,zud,zvd,pmfd,zmfds, &
                 zmfdq,zdmfdp,loddraf,ed)
!   ed(:,:)=pmfd(:,:)
  END IF

!-- 5.1 Recalculate cloud base massflux from a
!       cape closure for deep convection (ktype=1)
!       and by pbl equilibrum taking downdrafts into
!       account for shallow convection (ktype=2)

  DO jl = 1, klon
    zheat(jl) = 0.
    zcape(jl) = 0.
    zrelh(jl) = 0.
    zmfub1(jl) = zmfub(jl)
  END DO

  DO jl = 1, klon
    llo1 = ldcum(jl) .AND. ktype(jl) == 1
    ktop0=max(17,kctop(jl))
    DO jk = 1, klev
      IF (jk<=kcbot(jl) .AND. jk>kctop(jl)) THEN
        ikb = kcbot(jl)
        zro = paphp1(jl,jk)/(rd*ztenh(jl,jk))
        zdz = (paphp1(jl,jk)-paphp1(jl,jk-1))/(g*zro)
        zheat(jl) = zheat(jl) + ((pten(jl,jk-1)-pten(jl,jk)               &
              +g*zdz/cpd)/ztenh(jl,jk)+0.608*(pqen(jl,jk-1)-pqen(jl,jk))) &
              *(g*(pmfu(jl,jk)+pmfd(jl,jk)))/zro
        zcape(jl) = zcape(jl) + (g*(ptu(jl,jk)-ztenh(jl,jk))  &
              /ztenh(jl,jk)+g*0.608*(pqu(jl,jk)-zqenh(jl,jk)) &
              -g*plu(jl,jk))*zdz
      END IF
      IF (jk<=kcbot(jl) .AND. jk>ktop0) THEN
        dept=(PAPHP1(JL,JK)-PAPHP1(JL,JK-1))/(PAPHP1(JL,KCBOT(JL))-  &
            PAPHP1(JL,KTOP0))
        ZRELH(JL)=ZRELH(JL)+dept*PQEN(JL,JK)/PQSEN(JL,JK)
      END IF
     END DO
!llj    END IF
  END DO

!llj
  ztau = 3.*3600.
!!  ztau = 2.*3600.
!
  rhcr0 = 0.8

  DO jl = 1, klon
    IF (ldcum(jl) .AND. ktype(jl)==1) THEN
      IF(zrelh(jl).ge.rhcr0) THEN
      ikb = kcbot(jl)
      zmfub1(jl) = (zcape(jl)*zmfub(jl))/(zheat(jl)*ztau)
      zmfub1(jl) = MAX(zmfub1(jl),0.001)
      zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
      zmfub1(jl) = MIN(zmfub1(jl),zmfmax)
      ELSE
      zmfub1(jl) = 0.01
      zmfub(jl) = 0.01
      ldcum(jl)=.false.
      ENDIF
    END IF
  END DO

  ! Recalculate convective fluxes due to effect of
  ! downdrafts on boundary layer moisture budget
  ! for shallow convection (ktype=2)

!DIR$ IVDEP
!OCL NOVREC
  DO jl = 1, klon
    IF (ktype(jl)/=1) THEN
      ikb = kcbot(jl)
      llo1 = pmfd(jl,ikb) < 0. .AND. loddraf(jl)
      zeps = MERGE(cmfdeps,0.,llo1)
      zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zeps*zqd(jl,ikb) - &
               (1.-zeps)*zqenh(jl,ikb)
      zdqmin = MAX(0.01*zqenh(jl,ikb),1.E-10)
      zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
      llo1 = zdqpbl(jl) > 0. .AND. zqumqe > zdqmin .AND. ldcum(jl) .AND. &
             zmfub(jl) < zmfmax
      zmfub1(jl) = MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),zmfub(jl),llo1)
      zmfub1(jl) = MERGE(zmfub1(jl),zmfub(jl),(ktype(jl)==2) &
                   .AND. ABS(zmfub1(jl)-zmfub(jl))<0.2*zmfub(jl))
    END IF
  END DO
  DO jk = 1, klev
    DO jl = 1, klon
      IF (ldcum(jl)) THEN
        zfac = zmfub1(jl)/MAX(zmfub(jl),1.E-10)
        pmfd(jl,jk) = pmfd(jl,jk)*zfac
        zmfds(jl,jk) = zmfds(jl,jk)*zfac
        zmfdq(jl,jk) = zmfdq(jl,jk)*zfac
        zdmfdp(jl,jk) = zdmfdp(jl,jk)*zfac
      ELSE
        pmfd(jl,jk) = 0.
        zmfds(jl,jk) = 0.
        zmfdq(jl,jk) = 0.
        zdmfdp(jl,jk) = 0.
      END IF
    END DO

  END DO

  ! New values of cloud base mass flux

  DO jl = 1, klon
    IF (ldcum(jl)) zmfub(jl) = zmfub1(jl)
  END DO

!-- 6. Determine final cloud ascent for entraining plume
!      for penetrative convection (type=1),
!      for shallow to medium convection (type=2)
!      and for mid-level convection (type=3).

  CALL cuasc(klon,klev,klevp1,klevm1,ztenh,zqenh,zxenh,puen,pven,            &
             pten,pqen,pqsen,pgeo,zgeoh,papp1,paphp1,    &
             pverv,ilwmin,ldcum,ldland,ktype,ilab,ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr, &
             zmfus,zmfuq,zmful,plude,zdmfup,ihmin,zhhatt,zhcbase,zqsenh,kcbot,kctop,    &
             ictop0,icum,ztodt,eu,du)

!-- 7. Determine final convective fluxes in 'cuflx'

  CALL cuflx(klon,klev,klevp1,pqen,pqsen,ztenh,zqenh,&
    paphp1,zgeoh,kcbot,kctop,idtop,ktype,loddraf,ldcum,pmfu,  &
             pmfd,zmfus,zmfds,zmfuq,zmfdq,zmful,plude,zdmfup,zdmfdp,zrfl,prain,pten, &
             zsfl,zdpmel,itopm2,ztodt)

!-- 8. Update tendencies for t and q in subroutine cudtdq

  CALL cudtdq(klon,klev,klevp1,itopm2,paphp1,ldcum,pten,ptte,pqte,&
   pxtec,zmfus,zmfds,zmfuq,zmfdq,zmful,zdmfup,zdmfdp, &
             plude,zdpmel,prain,zrfl,zsfl,prsfc,pssfc,  &
             paprc,paprs,ztodt)

!-- 9. Update tendencies for u and u in subroutine cududv
  IF (lmfdudv) THEN
    CALL cududv(klon,klev,klevp1,itopm2,ktype,kcbot,paphp1,ldcum, &
                puen,pven,pvom,pvol,zuu,zud,zvu,zvd,pmfu,pmfd)

  END IF

  RETURN
END SUBROUTINE cumastr
