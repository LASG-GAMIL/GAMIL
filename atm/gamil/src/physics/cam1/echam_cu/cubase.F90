!+ calculates cloud base values (t and q) for cumulus parameterization
!+ $Id: cubase.f90,v 1.9 1999/07/19 19:07:06 m214030 Exp $

SUBROUTINE cubase(klon,klev,klevp1,klevm1,ptenh,pqenh,pgeoh,paph, &
&      ptu,pqu,plu,puen,pven,puu,pvu,ldcum,kcbot,klab)

  ! Description:
  !
  ! Calculates cloud base values (t and q) for cumulus parameterization
  !
  ! Method:
  !
  ! Produce cloud base values for cu-parametrization
  !
  ! This routine is called from *cumastr*.
  ! Input are environm. values of t,q,p,phi at half levels.
  ! It returns cloud base values and flags as follows;
  !    klab=1 for subcloud levels
  !    klab=2 for condensation level
  !
  ! Lift surface air dry-adiabatically to cloud base
  ! (non entraining plume,i.e.constant massflux)
  !
  ! *cuadjtq* for adjusting t and q due to condensation in ascent
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
  USE mo_constants,    ONLY: cpd,   &! specific heat at constant pressure
                             rcpd,  &! rcpd=1./cpd
                             vtmpc1  ! vtmpc1=rv/rd-1
  USE mo_cumulus_flux, ONLY: lmfdudv ! true if cumulus friction is switched on

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: klon, klev, klevm1, klevp1

  !  Array arguments with intent(In):
  REAL(r8), INTENT (IN) :: paph(klon,klevp1), pgeoh(klon,klev), pqenh(klon,klev),  &
&                      ptenh(klon,klev), puen(klon,klev), pven(klon,klev)

  !  Array arguments with intent(InOut):
  REAL(r8), INTENT (INOUT) :: plu(klon,klev), pqu(klon,klev), ptu(klon,klev), &
&                         puu(klon,klev), pvu(klon,klev)
  INTEGER, INTENT (INOUT) :: kcbot(klon), klab(klon,klev)
  LOGICAL, INTENT (INOUT) :: ldcum(klon)

  !  Local scalars: 
  REAL(r8) :: zbuo, zz
  INTEGER :: icall, ik, ikb, is, jk, jl

  !  Local arrays: 
  REAL(r8) :: zph(klon), zqold(klon)
  LOGICAL :: loflag(klon)

  !  External subroutines 
  EXTERNAL cuadjtq

  !  Executable statements 

!-- 1. Initialize values at lifting level

  DO jl = 1, klon
    klab(jl,klev) = 1
    kcbot(jl) = klevm1
    ldcum(jl) = .FALSE.
    puu(jl,klev) = puen(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))
    pvu(jl,klev) = pven(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))

!-- 2. Do ascent in subcloud layer,
!      check for existence of condensation level,
!      adjust t,q and l accordingly in *cuadjtq*,
!      check for buoyancy and set flags

  END DO
  DO jk = klevm1, 2, -1
    is = 0
    DO jl = 1, klon
      IF (klab(jl,jk+1)==1) is = is +1
      loflag(jl) = klab(jl,jk+1) == 1
      zph(jl) = paph(jl,jk)
    END DO
    IF (is==0) CYCLE
    DO jl = 1, klon
      IF (loflag(jl)) THEN
        pqu(jl,jk) = pqu(jl,jk+1)
        ptu(jl,jk) = (cpd*ptu(jl,jk+1)+pgeoh(jl,jk+1)-pgeoh(jl,jk))*rcpd
        zbuo = ptu(jl,jk)*(1.+vtmpc1*pqu(jl,jk)) - ptenh(jl,jk)*(1.+vtmpc1* &
&            pqenh(jl,jk)) + 0.5
        IF (zbuo>0.) klab(jl,jk) = 1
        zqold(jl) = pqu(jl,jk)
      END IF
    END DO

    ik = jk
    icall = 1
    CALL cuadjtq(klon,klev,ik,zph,ptu,pqu,loflag,icall)

!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      IF (loflag(jl) .AND. ABS(pqu(jl,jk)-zqold(jl))>0.) THEN
        klab(jl,jk) = 2
        plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
        zbuo = ptu(jl,jk)*(1.+vtmpc1*pqu(jl,jk)) - ptenh(jl,jk)*(1.+vtmpc1* &
&            pqenh(jl,jk)) + 0.5
        IF (zbuo>0.) THEN
          kcbot(jl) = jk
          ldcum(jl) = .TRUE.
        END IF
      END IF
    END DO

    ! Calculate averages of u and v for subcloud ara,.
    ! The values will be used to define cloud base values.

    IF (lmfdudv) THEN
      DO jl = 1, klon
        IF (jk>=kcbot(jl)) THEN
          puu(jl,klev) = puu(jl,klev) + puen(jl,jk)*(paph(jl,jk+1)-paph(jl,jk &
&              ))
          pvu(jl,klev) = pvu(jl,klev) + pven(jl,jk)*(paph(jl,jk+1)-paph(jl,jk &
&              ))
        END IF
      END DO

    END IF
  END DO

  IF (lmfdudv) THEN
    DO jl = 1, klon
      IF (ldcum(jl)) THEN
        ikb = kcbot(jl)
        zz = 1./(paph(jl,klevp1)-paph(jl,ikb))
        puu(jl,klev) = puu(jl,klev)*zz
        pvu(jl,klev) = pvu(jl,klev)*zz
      ELSE
        puu(jl,klev) = puen(jl,klevm1)
        pvu(jl,klev) = pven(jl,klevm1)
      END IF
    END DO
  END IF

  RETURN
END SUBROUTINE cubase
