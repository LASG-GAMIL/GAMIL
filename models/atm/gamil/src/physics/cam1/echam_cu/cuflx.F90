!+ does the final calculation of convective fluxes

SUBROUTINE cuflx(klon,klev,klevp1,pqen,pqsen,ptenh,pqenh,&
&   paphp1,pgeoh,kcbot,kctop,kdtop,ktype,lddraf,ldcum,pmfu,   &
                 pmfd,pmfus,pmfds,pmfuq,pmfdq,pmful,plude,pdmfup,pdmfdp,prfl,prain,pten, &
                 psfl,pdpmel,ktopm2,ztodt)

  ! Description:
  !
  ! Does the final calculation of convective fluxes.
  !
  ! Method:
  !
  ! This routine does the final calculation of convective
  ! fluxes in the cloud layer and in the subcloud layer
  !
  ! This routine is called from *cumastr*.
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, July 1989, original source
  ! M. Tiedtke, ECMWF, December 1989, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE mo_constants,     ONLY: cpd,   &! specific heat at constant pressure
                              alf,   &! latent heat for vaporisation
                              g,     &! gravity acceleration
                              tmelt   ! temperature of fusion of ice
  USE time_manager,     ONLY: is_first_step,         &! initial run
                              is_first_restart_step, &! restart run
                              get_step_size           ! step size
   use pmgrid,             only: plev, plevp !Needed for hypm passed to vd_inti

  IMPLICIT NONE

#include <comhyb.h>

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: klev, klevp1, klon

  !  Array arguments with intent(In):
  REAL(r8), INTENT (IN) :: paphp1(klon,klevp1), pgeoh(klon,klev), pqen(klon,klev), &
                       pqenh(klon,klev), pqsen(klon,klev), pten(klon,klev),    &
                       ptenh(klon,klev),ztodt
  INTEGER, INTENT (IN) :: kcbot(klon), kctop(klon), kdtop(klon)
  LOGICAL, INTENT (IN) :: ldcum(klon)

  !  Array arguments with intent(InOut):
  REAL(r8), INTENT (INOUT) :: pdmfdp(klon,klev), pdmfup(klon,klev), pdpmel(klon,klev), &
                          plude(klon,klev), pmfd(klon,klev), pmfdq(klon,klev),     &
                          pmfds(klon,klev), &
                          pmfu(klon,klev), pmful(klon,klev), pmfuq(klon,klev),     &
                          pmfus(klon,klev), &
!pmfuxt(klon,klev,ktrac),               &
                          prain(klon), prfl(klon), psfl(klon)
  INTEGER, INTENT (INOUT) :: ktype(klon)
  LOGICAL, INTENT (INOUT) :: lddraf(klon)

  !  Scalar arguments with intent(Out):
  INTEGER, INTENT (OUT) :: ktopm2

  !  Local scalars: 
  REAL(r8) :: zcons1, zcons2, zcucov, zdpevap, zdrfl, zfac, zrfl, zrfln, zrmin, &
          zrnew, zrsum, zsnmlt, ztmelp2, ztmst, zzp  ! zrsum new for AMIP2
  INTEGER :: ikb, jk, jl, jt

  !  Local arrays: 
  REAL(r8) :: zpsubcl(klon)
  REAL(r8) :: cevapcu(klev)

  !  Intrinsic functions 
  INTRINSIC MAX, MERGE, MIN, SQRT


  !  Executable statements

  !  eta level for CAM
  cevapcu(:)=1.93E-6*261.*SQRT(1.E3/(38.3*0.293)*SQRT(hyam(:)+hybm(:)))*0.5/g 

  !  Specify constants
  ztmst = ztodt
  IF (is_first_step()) ztmst = 0.5*ztodt
  zcons1 = cpd/(alf*g*ztmst)
  zcons2 = 1./(g*ztmst)
  zcucov = 0.05
  ztmelp2 = tmelt + 2.

!-- 1. Determine final convective fluxes

!  itop = klev
  DO jl = 1, klon
!    itop = MIN(itop,kctop(jl))
    IF ( .NOT. ldcum(jl) .OR. kdtop(jl)<kctop(jl)) lddraf(jl) = .FALSE.
    IF ( .NOT. ldcum(jl)) ktype(jl) = 0
  END DO
!  ktopm2 = itop - 2
  ktopm2 = 2
  DO jk = ktopm2, klev
!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      IF (ldcum(jl) .AND. jk>=kctop(jl)-1) THEN
        pmfus(jl,jk) = pmfus(jl,jk) - pmfu(jl,jk)*(cpd*ptenh(jl,jk)+pgeoh(jl,jk))
        pmfuq(jl,jk) = pmfuq(jl,jk) - pmfu(jl,jk)*pqenh(jl,jk)
        IF (lddraf(jl) .AND. jk>=kdtop(jl)) THEN
          pmfds(jl,jk) = pmfds(jl,jk) - pmfd(jl,jk)*(cpd*ptenh(jl,jk)+pgeoh(jl,jk))
          pmfdq(jl,jk) = pmfdq(jl,jk) - pmfd(jl,jk)*pqenh(jl,jk)
        ELSE
          pmfd(jl,jk) = 0.
          pmfds(jl,jk) = 0.
          pmfdq(jl,jk) = 0.
          pdmfdp(jl,jk-1) = 0.
        END IF
      ELSE
        pmfu(jl,jk) = 0.
        pmfd(jl,jk) = 0.
        pmfus(jl,jk) = 0.
        pmfds(jl,jk) = 0.
        pmfuq(jl,jk) = 0.
        pmfdq(jl,jk) = 0.
        pmful(jl,jk) = 0.
        pdmfup(jl,jk-1) = 0.
        pdmfdp(jl,jk-1) = 0.
        plude(jl,jk-1) = 0.
      END IF
    END DO

  END DO
  DO jk = ktopm2, klev
!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      IF (ldcum(jl) .AND. jk>kcbot(jl)) THEN
        ikb = kcbot(jl)
        zzp = ((paphp1(jl,klevp1)-paphp1(jl,jk))/(paphp1(jl,klevp1)-paphp1(jl,ikb)))
        zzp = MERGE(zzp**2,zzp,ktype(jl)==3)
        pmfu(jl,jk) = pmfu(jl,ikb)*zzp
        pmfus(jl,jk) = pmfus(jl,ikb)*zzp
        pmfuq(jl,jk) = pmfuq(jl,ikb)*zzp
        pmful(jl,jk) = pmful(jl,ikb)*zzp
      END IF
    END DO

  END DO

!-- 2. Calculate rain/snow fall rates
!      Calculate melting of snow
!      Calculate evaporation of precip

  DO jl = 1, klon
    prfl(jl) = 0.
    psfl(jl) = 0.
    prain(jl) = 0.
  END DO
  DO jk = ktopm2, klev
    DO jl = 1, klon
      IF (ldcum(jl)) THEN
        prain(jl) = prain(jl) + pdmfup(jl,jk)
        IF (pten(jl,jk)>tmelt) THEN
          prfl(jl) = prfl(jl) + pdmfup(jl,jk) + pdmfdp(jl,jk)
          IF (psfl(jl)>0. .AND. pten(jl,jk)>ztmelp2) THEN
            zfac = zcons1*(paphp1(jl,jk+1)-paphp1(jl,jk))
            zsnmlt = MIN(psfl(jl),zfac*(pten(jl,jk)-ztmelp2))
            pdpmel(jl,jk) = zsnmlt
            psfl(jl) = psfl(jl) - zsnmlt
            prfl(jl) = prfl(jl) + zsnmlt
          END IF
        ELSE
          psfl(jl) = psfl(jl) + pdmfup(jl,jk) + pdmfdp(jl,jk)
        END IF
      END IF
    END DO
  END DO
  DO jl = 1, klon
    prfl(jl) = MAX(prfl(jl),0.)
    psfl(jl) = MAX(psfl(jl),0.)
    zpsubcl(jl) = prfl(jl) + psfl(jl)
  END DO
  DO jk = ktopm2, klev
    DO jl = 1, klon
      IF (ldcum(jl) .AND. jk>=kcbot(jl) .AND. zpsubcl(jl)>1.E-20) THEN
        zrfl = zpsubcl(jl)
        zrnew = (MAX(0.,SQRT(zrfl/zcucov)-cevapcu(jk)*(paphp1(jl,jk+1)- &
                paphp1(jl,jk))*MAX(0.,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
        zrmin = zrfl - zcucov*MAX(0.,0.8*pqsen(jl,jk)-pqen(jl,jk))*zcons2* &
                (paphp1(jl,jk+1)-paphp1(jl,jk))
        zrnew = MAX(zrnew,zrmin)
        zrfln = MAX(zrnew,0.)
        zdrfl = MIN(0.,zrfln-zrfl)
        pdmfup(jl,jk) = pdmfup(jl,jk) + zdrfl
        zpsubcl(jl) = zrfln
      END IF
    END DO
  END DO

    ! ------ CORRECTION OF PRECIPITATION EVAPORATION ERROR (AMIP2) --------
    !                (Found by M. Werner)
    DO jl = 1, klon
      zrsum    = prfl(jl)+psfl(jl)
      zdpevap  = zpsubcl(jl) - zrsum
      prfl(jl) = prfl(jl) + zdpevap*prfl(jl)*(1./MAX(1.E-20,zrsum))
      psfl(jl) = psfl(jl) + zdpevap*psfl(jl)*(1./MAX(1.E-20,zrsum))
    END DO

  RETURN
END SUBROUTINE cuflx
