#include <misc.h>
#include <params.h>

module zm_conv

!---------------------------------------------------------------------------------
! Purpose:
!
! Interface from Zhang-McFarlane convection scheme, includes evaporation of convective
! precip from the ZM scheme
!
! Author: Byron Boville, from code in tphysbc
!
!---------------------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,    only: pcols, pver, pverp
  use moistconvection,  only: cp, grav, rgrav, rgas, limcnv

  implicit none

  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public zm_convi                 ! ZM schemea
  public zm_convr                 ! ZM schemea
  public zm_conv_evap             ! evaporation of precip from ZM schemea
!
! Private data
!
   public rl, cpres, capelmt
   real(r8) qmin       ! Lower bound on specific humidity
   real(r8) rl         ! wg latent heat of vaporization.
   real(r8) cpres      ! specific heat at constant pressure in j/kg-degk.
   real(r8), parameter :: capelmt = 70.  ! threshold value for cape for deep convection.
!   real(r8), parameter :: capelmt = 80.  ! threshold value for cape for deep convection. Modified by LI Lijuan

contains

    subroutine zm_convi(tmelt, epsilo, latvap, cpair)
!
! Initialization of ZM constants
!
        real(r8), intent(in) :: tmelt   ! Freezing temperature
        real(r8), intent(in) :: epsilo  ! Freezing temperature
        real(r8), intent(in) :: latvap  ! Latent heat of vaporization
        real(r8), intent(in) :: cpair   ! Specific heat of dry air (J/K/kg)

#include <guang.h>

        cpres   = cpair
        a       = 21.656
        b       = 5418.
        c1      = 6.112
        c2      = 17.67
        c3      = 243.5
        eps1    = epsilo
        qmin    = 1.E-20
        tfreez  = tmelt
        rl      = latvap

    end subroutine zm_convi

subroutine zm_convr(lchnk   ,ncol    , &
                    t       ,qh      ,pcpc    ,jctop   ,jcbot   , &
                    pblh    ,zm      ,geos    ,zi      ,qtg     , &
                    ttg     ,pap     ,paph    ,dpp     ,ts      , &
                    delt    ,mcon    ,cme     ,nstep   ,          &
                    tpert   ,dlf     ,pflx    ,zdu     ,cmfdqr  , &
                    mu2     ,md2     ,du2     ,eu2     ,ed2     , &
                    dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
                    lengath ,ql      ,ptendlaq  ,ptendlat  )
!-----------------------------------------------------------------------
!
! Purpose:
! Main driver for zhang-mcfarlane convection scheme
!
! Method:
! performs deep convective adjustment based on mass-flux closure
! algorithm.
!
! Author:guang jun zhang, m.lazare, n.mcfarlane. CAM Contact: P. Rasch
!
! This is contributed code not fully standardized by the CAM core group.
! All variables have been typed, where most are identified in comments
! The current procedure will be reimplemented in a subsequent version
! of the CAM where it will include a more straightforward formulation
! and will make use of the standard CAM nomenclature
!
!-----------------------------------------------------------------------
   use tracers, only: pcnst, pnats

!
! ************************ index of variables **********************
!
!  wg * alpha    array of vertical differencing used (=1. for upstream).
!  w  * cape     convective available potential energy.
!  wg * capeg    gathered convective available potential energy.
!  c  * capelmt  threshold value for cape for deep convection.
!  ic  * cpres    specific heat at constant pressure in j/kg-degk.
!  i  * dpp
!  ic  * delt     length of model time-step in seconds.
!  wg * dp       layer thickness in mbs (between upper/lower interface).
!  wg * dqdt     mixing ratio tendency at gathered points.
!  wg * dsdt     dry static energy ("temp") tendency at gathered points.
!  wg * dudt     u-wind tendency at gathered points.
!  wg * dvdt     v-wind tendency at gathered points.
!  wg * dsubcld  layer thickness in mbs between lcl and maxi.
!  ic  * grav     acceleration due to gravity in m/sec2.
!  wg * du       detrainment in updraft. specified in mid-layer
!  wg * ed       entrainment in downdraft.
!  wg * eu       entrainment in updraft.
!  wg * hmn      moist static energy.
!  wg * hsat     saturated moist static energy.
!  w  * ideep    holds position of gathered points vs longitude index.
!  ic  * pver     number of model levels.
!  wg * j0       detrainment initiation level index.
!  wg * jd       downdraft   initiation level index.
!  ic  * jlatpr   gaussian latitude index for printing grids (if needed).
!  wg * jt       top  level index of deep cumulus convection.
!  w  * lcl      base level index of deep cumulus convection.
!  wg * lclg     gathered values of lcl.
!  w  * lel      index of highest theoretical convective plume.
!  wg * lelg     gathered values of lel.
!  w  * lon      index of onset level for deep convection.
!  w  * maxi     index of level with largest moist static energy.
!  wg * maxg     gathered values of maxi.
!  wg * mb       cloud base mass flux.
!  wg * mc       net upward (scaled by mb) cloud mass flux.
!  wg * md       downward cloud mass flux (positive up).
!  wg * mu       upward   cloud mass flux (positive up). specified
!                at interface
!  ic  * msg      number of missing moisture levels at the top of model.
!  w  * p        grid slice of ambient mid-layer pressure in mbs.
!  i  * pblt     row of pbl top indices.
!  w  * pcpdh    scaled surface pressure.
!  w  * pf       grid slice of ambient interface pressure in mbs.
!  wg * pg       grid slice of gathered values of p.
!  w  * q        grid slice of mixing ratio.
!  wg * qd       grid slice of mixing ratio in downdraft.
!  wg * qg       grid slice of gathered values of q.
!  i/o * qh       grid slice of specific humidity.
!  w  * qh0      grid slice of initial specific humidity.
!  wg * qhat     grid slice of upper interface mixing ratio.
!  wg * ql       grid slice of cloud liquid water.
!  wg * qs       grid slice of saturation mixing ratio.
!  w  * qstp     grid slice of parcel temp. saturation mixing ratio.
!  wg * qstpg    grid slice of gathered values of qstp.
!  wg * qu       grid slice of mixing ratio in updraft.
!  ic  * rgas     dry air gas constant.
!  wg * rl       latent heat of vaporization.
!  w  * s        grid slice of scaled dry static energy (t+gz/cp).
!  wg * sd       grid slice of dry static energy in downdraft.
!  wg * sg       grid slice of gathered values of s.
!  wg * shat     grid slice of upper interface dry static energy.
!  wg * su       grid slice of dry static energy in updraft.
!  wg * sumde    row of vertically-integrated moist static energy
!                change.
!  wg * sumdq    row of vertically-integrated scaled mixing ratio
!                change.
!  wg * sumdt    row of vertically-integrated dry static energy change.
!  wg * sumq     row of vertically-integrated mixing ratio change.
!  i/o * t
!  o  * jctop    row of top-of-deep-convection indices passed out.
!  o  * jcbot    row of base of cloud indices passed out.
!  wg * tg       grid slice of gathered values of t.
!  w  * tl       row of parcel temperature at lcl.
!  wg * tlg      grid slice of gathered values of tl.
!  w  * tp       grid slice of parcel temperatures.
!  wg * tpg      grid slice of gathered values of tp.
!  i/o * u        grid slice of u-wind (real).
!  wg * ug       grid slice of gathered values of u.
!  i/o * utg      grid slice of u-wind tendency (real).
!  i/o * v        grid slice of v-wind (real).
!  w  * va       work array re-used by called subroutines.
!  wg * vg       grid slice of gathered values of v.
!  i/o * vtg      grid slice of v-wind tendency (real).
!  i  * w        grid slice of diagnosed large-scale vertical velocity.
!  w  * z        grid slice of ambient mid-layer height in metres.
!  w  * zf       grid slice of ambient interface height in metres.
!  wg * zfg      grid slice of gathered values of zf.
!  wg * zg       grid slice of gathered values of z.
!
!-----------------------------------------------------------------------
!
! multi-level i/o fields:
!  i      => input arrays.
!  i/o    => input/output arrays.
!  w      => work arrays.
!  wg     => work arrays operating only on gathered points.
!  ic     => input data constants.
!  c      => data constants pertaining to subroutine itself.
!
! input/output arguments:
!
   real(r8),intent(inout) :: t(pcols,pver)              ! grid slice of temperature at mid-layer.
   real(r8),intent(inout) :: qh(pcols,pver,pcnst+pnats) ! grid slice of specific humidity.
   real(r8) u(pcols,pver)              ! grid slice of U-wind (real).
   real(r8) v(pcols,pver)              ! grid slice of U-wind (real).
   real(r8),intent(inout) :: qtg(pcols,pver)
   real(r8),intent(inout) :: ttg(pcols,pver)
!      real(r8),intent(inout) :: utg(pcols,pver)        ! grid slice of u-wind tendency (real).
!      real(r8),intent(inout) :: vtg(pcols,pver)        ! grid slice of v-wind tendency (real).
!
! input arguments
!
   integer, intent(in) :: lchnk                   ! chunk identifier
   integer, intent(in) :: ncol                    ! number of atmospheric columns

   real(r8), intent(in) :: pap(pcols,pver)
   real(r8), intent(in) :: paph(pcols,pver+1)
   real(r8), intent(in) :: dpp(pcols,pver)        ! local sigma half-level thickness (i.e. dshj).
   real(r8), intent(in) :: zm(pcols,pver)
   real(r8), intent(in) :: ptendlaq(pcols,pver)
   real(r8), intent(in) :: ptendlat(pcols,pver)
   real(r8), intent(in) :: geos(pcols)
   real(r8), intent(in) :: zi(pcols,pver+1)
   real(r8), intent(in) :: pblh(pcols)
   real(r8) zs(pcols)
   real(r8), intent(in) :: tpert(pcols)
!
! output arguments
!
   real(r8) pcpck(pcols,pver)
   real(r8), intent(out) :: mcon(pcols,pverp)
   real(r8) dlg(pcols,pver)    ! gathrd version of the detraining cld h2o tend
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8) pflxg(pcols,pverp) ! gather precip flux at each level
   real(r8) cug(pcols,pver)    ! gathered condensation rate
   real(r8) evpg(pcols,pver)   ! gathered evap rate of rain in downdraft
   real(r8) mumax(pcols)
   real(r8), intent(out) :: cme(pcols,pver)
   real(r8), intent(out) :: zdu(pcols,pver)
   real(r8), intent(out) :: cmfdqr(pcols,pver)
! move these vars from local storage to output so that convective
! transports can be done in outside of conv_cam.
   real(r8), intent(out) :: mu2(pcols,pver)
   real(r8), intent(out) :: eu2(pcols,pver)
   real(r8), intent(out) :: du2(pcols,pver)
   real(r8), intent(out) :: md2(pcols,pver)
   real(r8), intent(out) :: ed2(pcols,pver)
   real(r8), intent(out) :: dp(pcols,pver)       ! wg layer thickness in mbs (between upper/lower interface).
   real(r8), intent(out) :: dsubcld(pcols)       ! wg layer thickness in mbs between lcl and maxi.
   integer jt(pcols)                          ! wg top  level index of deep cumulus convection.
   integer maxg(pcols)                        ! wg gathered values of maxi.
   integer ideep(pcols)                       ! w holds position of gathered points vs longitude index.
   integer lengath
!     diagnostic field used by chem/wetdep codes
   real(r8) ql(pcols,pver)                    ! wg grid slice of cloud liquid water.
!
! single-level i/o fields:
!
!       input arguments
!
   real(r8), intent(in) :: ts(pcols)
   real(r8) pblt(pcols)           ! i row of pbl top indices.
!
!       input/output arguments:
!
   real(r8) paprc(pcols)
   real(r8) paprs(pcols)
!
!       output arguments:
!
   real(r8), intent(out) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8), intent(out) :: jcbot(pcols)  ! o row of base of cloud indices passed out.
   real(r8) pcpr(pcols)
   real(r8) pcps(pcols)
   real(r8), intent(out) :: pcpc(pcols)
!
!-----------------------------------------------------------------------
!
! general work fields (local variables):
!
   real(r8) q(pcols,pver)              ! w  grid slice of mixing ratio.
   real(r8) p(pcols,pver)              ! w  grid slice of ambient mid-layer pressure in mbs.
   real(r8) z(pcols,pver)              ! w  grid slice of ambient mid-layer height in metres.
   real(r8) s(pcols,pver)              ! w  grid slice of scaled dry static energy (t+gz/cp).
   real(r8) tp(pcols,pver)             ! w  grid slice of parcel temperatures.
   real(r8) zf(pcols,pver+1)           ! w  grid slice of ambient interface height in metres.
   real(r8) pf(pcols,pver+1)           ! w  grid slice of ambient interface pressure in mbs.
   real(r8) qstp(pcols,pver)           ! w  grid slice of parcel temp. saturation mixing ratio.

   real(r8) cape(pcols)                ! w  convective available potential energy.
   real(r8) negcape(pcols)             ! w  negative convective available potential energy.
   real(r8) tl(pcols)                  ! w  row of parcel temperature at lcl.
   real(r8) sumq(pcols)                ! wg row of vertically-integrated mixing ratio change.
   real(r8) pcpdh(pcols)               ! w  scaled surface pressure.
!      real(r8) sumdt(pcols)           ! wg row of vertically-integrated dry static energy change.
!      real(r8) sumdq(pcols)           ! wg row of vertically-integrated scaled mixing ratio change
!      real(r8) sumde(pcols)           ! wg row of vertically-integrated moist static energy change

   integer lcl(pcols)                  ! w  base level index of deep cumulus convection.
   integer lel(pcols)                  ! w  index of highest theoretical convective plume.
   integer lon(pcols)                  ! w  index of onset level for deep convection.
   integer maxi(pcols)                 ! w  index of level with largest moist static energy.
   integer index(pcols)
   real(r8) precip
!
! gathered work fields:
!
   real(r8) qg(pcols,pver)             ! wg grid slice of gathered values of q.
   real(r8) tg(pcols,pver)             ! w  grid slice of temperature at interface.
   real(r8) ptendtg(pcols,pver)        ! w  grid slice of temp. tend. at interface.
   real(r8) ptendqg(pcols,pver)        ! w  grid slice of moist. tend. at interface.
   real(r8) pg(pcols,pver)             ! wg grid slice of gathered values of p.
   real(r8) zg(pcols,pver)             ! wg grid slice of gathered values of z.
   real(r8) sg(pcols,pver)             ! wg grid slice of gathered values of s.
   real(r8) tpg(pcols,pver)            ! wg grid slice of gathered values of tp.
   real(r8) zfg(pcols,pver+1)          ! wg grid slice of gathered values of zf.
   real(r8) qstpg(pcols,pver)          ! wg grid slice of gathered values of qstp.
   real(r8) ug(pcols,pver)             ! wg grid slice of gathered values of u.
   real(r8) vg(pcols,pver)             ! wg grid slice of gathered values of v.
   real(r8) cmeg(pcols,pver)

   real(r8) cmfdqrg(pcols,pver)
   real(r8) capeg(pcols)               ! wg gathered convective available potential energy.
   real(r8) negcapeg(pcols)            ! DONG Li: gathered negative convective available potential energy
   real(r8) tlg(pcols)                 ! wg grid slice of gathered values of tl.
   integer lclg(pcols)                 ! wg gathered values of lcl.
   integer lelg(pcols)
!
! work fields arising from gathered calculations.
!
   real(r8) mu(pcols,pver)             ! wg upward cloud mass flux (positive up). specified at interface
   real(r8) eu(pcols,pver)             ! wg entrainment in updraft.
   real(r8) dqdt(pcols,pver)           ! wg mixing ratio tendency at gathered points.
   real(r8) dsdt(pcols,pver)           ! wg dry static energy ("temp") tendency at gathered points.
   real(r8) du(pcols,pver)             ! wg detrainment in updraft. specified in mid-layer
   real(r8) md(pcols,pver)             ! wg downward cloud mass flux (positive up).
   real(r8) ed(pcols,pver)             ! wg entrainment in downdraft.
!      real(r8) alpha(pcols,pver)      ! array of vertical differencing used (=1. for upstream).
   real(r8) sd(pcols,pver)             ! wg grid slice of dry static energy in downdraft.
   real(r8) qd(pcols,pver)             ! wg grid slice of mixing ratio in downdraft.
   real(r8) mc(pcols,pver)             ! wg net upward (scaled by mb) cloud mass flux.
   real(r8) qhat(pcols,pver)           ! wg grid slice of upper interface mixing ratio.
   real(r8) qu(pcols,pver)             ! wg grid slice of mixing ratio in updraft.
   real(r8) su(pcols,pver)             ! wg grid slice of dry static energy in updraft.
   real(r8) qs(pcols,pver)             ! wg grid slice of saturation mixing ratio.
   real(r8) shat(pcols,pver)           ! wg grid slice of upper interface dry static energy.
   real(r8) hmn(pcols,pver)            ! wg moist static energy.
   real(r8) hsat(pcols,pver)           ! wg saturated moist static energy.
   real(r8) qlg(pcols,pver)
   real(r8) dudt(pcols,pver)           ! wg u-wind tendency at gathered points.
   real(r8) dvdt(pcols,pver)           ! wg v-wind tendency at gathered points.
!      real(r8) ud(pcols,pver)
!      real(r8) vd(pcols,pver)

!      real(r8) deltat(pcols,pver)
!      real(r8) deltaq(pcols,pver)

   real(r8) mb(pcols)                  ! wg cloud base mass flux.




   real(r8) totpcp(pcols)
   real(r8) totevp(pcols)

   integer jlcl(pcols)
   integer j0(pcols)                 ! wg detrainment initiation level index.
   integer jd(pcols)                 ! wg downdraft initiation level index.

   real(r8) delt                     ! length of model time-step in seconds.

   integer i
   integer ii
   integer k
   integer msg                      !  ic number of missing moisture levels at the top of model.
   integer nstep
   real(r8) psdiss
   real(r8) psevap
   real(r8) psheat
   real(r8) psrain
   real(r8) qdifr
   real(r8) qeff
   real(r8) sdifr
!
#include <guang.h>
!--------------------------Data statements------------------------------
!
   logical momentm
   data momentm/.FALSE./
!
! Set internal variable "msg" (convection limit) to "limcnv-1"
!
   msg = limcnv - 1
!
! initialize necessary arrays.
! zero out variables not used in cam
!
   qtg(:,:) = 0.
   ttg(:,:) = 0.
   mcon(:,:) = 0.
   do i = 1,ncol
      paprc(i) = 0.
      paprs(i) = 0.
      pcpc(i) = 0
   end do
   psdiss = 0.
   psheat = 0.
   psevap = 0.
   psrain = 0.
!      jlatpr = 32
!
! initialize convective tendencies
!
   do k = 1,pver
      do i = 1,ncol
         dqdt(i,k) = 0.
         dsdt(i,k) = 0.
         dudt(i,k) = 0.
         dvdt(i,k) = 0.
!            deltaq(i,k) = qh(i,k,1)
!            deltat(i,k) = t(i,k)
         pcpck(i,k) = 0.
         pflx(i,k) = 0.
         pflxg(i,k) = 0.
         cme(i,k) = 0.
!++bee
         cmfdqr(i,k) = 0.
!--bee
         zdu(i,k) = 0.
         ql(i,k) = 0.
         qlg(i,k) = 0.
      end do
   end do
   do i = 1,ncol
      pflx(i,pverp) = 0
      pflxg(i,pverp) = 0
   end do
   if (.not.momentm) then
      do k = 1,pver
         do i = 1,ncol
            u(i,k) = 0.
            v(i,k) = 0.
         end do
      end do
   end if
!
   do i = 1,ncol
      pblt(i) = pver
      pcpr(i) = 0.
      pcps(i) = 0.
      dsubcld(i) = 0.
      sumq(i) = 0.
!         sumdt(i) = 0.
!         sumdq(i) = 0.
      pcpdh(i) = rgrav
      jctop(i) = pver
      jcbot(i) = 1
   end do
!
! calculate local pressure (mbs) and height (m) for both interface
! and mid-layer locations.
!
   do i = 1,ncol
      zs(i) = geos(i)*rgrav
      pf(i,pver+1) = paph(i,pver+1)*0.01
      zf(i,pver+1) = zi(i,pver+1) + zs(i)
   end do
   do k = 1,pver
      do i = 1,ncol
         p(i,k) = pap(i,k)*0.01
         pf(i,k) = paph(i,k)*0.01
         z(i,k) = zm(i,k) + zs(i)
         zf(i,k) = zi(i,k) + zs(i)
      end do
   end do
!
   do k = pver - 1,msg + 1,-1
      do i = 1,ncol
         if (abs(z(i,k)-zs(i)-pblh(i)) < (zf(i,k)-zf(i,k+1))*0.5) pblt(i) = k
      end do
   end do
!
! store incoming specific humidity field for subsequent calculation
! of precipitation (through change in storage).
! convert from specific humidity (bounded by qmin) to mixing ratio.
! define dry static energy (normalized by cp).
!
   do k = 1,pver
      do i = 1,ncol
         qeff = max(qh(i,k,1),qmin)
         q(i,k) = qeff
         s(i,k) = t(i,k) + (grav/cpres)*z(i,k)
         tp(i,k)=0.0
         shat(i,k) = s(i,k)
         qhat(i,k) = q(i,k)
         dp(i,k) = dpp(i,k)*0.01
         qg(i,k) = q(i,k)
         tg(i,k) = t(i,k)
         pg(i,k) = p(i,k)
         zg(i,k) = z(i,k)
         sg(i,k) = s(i,k)
         tpg(i,k) = tp(i,k)
         zfg(i,k) = zf(i,k)
         qstpg(i,k) = q(i,k)
         ug(i,k) = u(i,k)
         vg(i,k) = v(i,k)
         dlg(i,k) = 0
         dlf(i,k) = 0
      end do
   end do
   do i = 1,ncol
      zfg(i,pver+1) = zf(i,pver+1)
      capeg(i) = 0.
      lclg(i) = 1
      lelg(i) = pver
      maxg(i) = 1
      tlg(i) = 400.
      dsubcld(i) = 0.
   end do

    ! evaluate covective available potential energy (CAPE).
    call buoyan(lchnk,  ncol,   &
                q,      t,      p,      z,    pf,    &
                tp,     qstp,   tl,     rl,   cape,  &
                pblt,   lcl,    lel,    lon,  maxi,  &
                rgas,   grav,   cpres,  msg,  nstep, &
                tpert,  negcape)

    ! determine whether grid points will undergo some deep convection
    ! (ideep=1) or not (ideep=0), based on values of cape, lcl, lel
    ! (require cape > 0 and lel < lcl as minimum conditions).
    lengath = 0
    do i = 1, ncol
        if (cape(i) > capelmt) then
            lengath = lengath+1
            index(lengath) = i
        end if
    end do

    if (lengath == 0) return
!CDIR$ IVDEP
        do ii = 1, lengath
            i = index(ii)
            ideep(ii) = i
        end do
!c        jyes = 0
!c        jno = nlon - 1 + 2
!c        do il = 1,nlon
!c           if (cape(il).gt.capelmt) then
!c              jyes = jyes + 1
!c              ideep(jyes) = il
!c           else
!c              jno = jno - 1
!c              ideep(jno) = il
!c           end if
!c        end do
!c        lengath = jyes
!c        if (lengath.eq.0) return
!
! obtain gathered arrays necessary for ensuing calculations.
!
        do k = 1, pver
            do i = 1, lengath
                dp(i,k)      = 0.01*dpp(ideep(i),k)
                qg(i,k)      = q(       ideep(i),k)
                tg(i,k)      = t(       ideep(i),k)
                ptendtg(i,k) = ptendlat(ideep(i),k)
                ptendqg(i,k) = ptendlaq(ideep(i),k)
                pg(i,k)      = p(       ideep(i),k)
                zg(i,k)      = z(       ideep(i),k)
                sg(i,k)      = s(       ideep(i),k)
                tpg(i,k)     = tp(      ideep(i),k)
                zfg(i,k)     = zf(      ideep(i),k)
                qstpg(i,k)   = qstp(    ideep(i),k)
                ug(i,k)      = u(       ideep(i),k)
                vg(i,k)      = v(       ideep(i),k)
            end do
        end do

        do i = 1, lengath
            zfg(i,pver+1) = zf(         ideep(i),pver+1)
        end do
        do i = 1, lengath
            capeg(i)    = cape(   ideep(i))
            negcapeg(i) = negcape(ideep(i))
            lclg(i)     = lcl(    ideep(i))
            lelg(i)     = lel(    ideep(i))
            maxg(i)     = maxi(   ideep(i))
            tlg(i)      = tl(     ideep(i))
        end do
!
! calculate sub-cloud layer pressure "thickness" for use in
! closure and tendency routines.
!
   do k = msg + 1,pver
      do i = 1,lengath
         if (k >= maxg(i)) then
            dsubcld(i) = dsubcld(i) + dp(i,k)
         end if
      end do
   end do
!
! define array of factors (alpha) which defines interfacial
! values, as well as interfacial values for (q,s) used in
! subsequent routines.
!
   do k = msg + 2,pver
      do i = 1,lengath
!            alpha(i,k) = 0.5
         sdifr = 0.
         qdifr = 0.
         if (sg(i,k) > 0. .or. sg(i,k-1) > 0.) &
            sdifr = abs((sg(i,k)-sg(i,k-1))/max(sg(i,k-1),sg(i,k)))
         if (qg(i,k) > 0. .or. qg(i,k-1) > 0.) &
            qdifr = abs((qg(i,k)-qg(i,k-1))/max(qg(i,k-1),qg(i,k)))
         if (sdifr > 1.E-6) then
            shat(i,k) = log(sg(i,k-1)/sg(i,k))*sg(i,k-1)*sg(i,k)/(sg(i,k-1)-sg(i,k))
         else
            shat(i,k) = 0.5* (sg(i,k)+sg(i,k-1))
         end if
         if (qdifr > 1.E-6) then
            qhat(i,k) = log(qg(i,k-1)/qg(i,k))*qg(i,k-1)*qg(i,k)/(qg(i,k-1)-qg(i,k))
         else
            qhat(i,k) = 0.5* (qg(i,k)+qg(i,k-1))
         end if
      end do
   end do
!
! obtain cloud properties.
!
   call cldprp(lchnk   , &
               qg      ,tg      ,ug      ,vg      ,pg      , &
               zg      ,sg      ,mu      ,eu      ,du      , &
               md      ,ed      ,sd      ,qd      ,mc      , &
               qu      ,su      ,zfg     ,qs      ,hmn     , &
               hsat    ,shat    ,qlg     ,totpcp  ,totevp  , &
               cmeg    ,maxg    ,lelg    ,jt      ,jlcl    , &
               maxg    ,j0      ,jd      ,rl      ,lengath , &
               rgas    ,grav    ,cpres   ,msg     ,nstep   , &
                        pflxg   ,evpg    ,cug     ,mu2     , &
               eu2     ,du2     ,md2     ,ed2     ,cmfdqrg , &
               limcnv  )
!
! convert detrainment from units of "1/m" to "1/mb".
!
   do k = msg + 1,pver
      do i = 1,lengath
         du(i,k) = du(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         eu(i,k) = eu(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         ed(i,k) = ed(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cug(i,k) = cug(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
!++bee, bug fix from /home/pjr/ccm2/omega0.10.1/fspj01/.
         cmeg(i,k) = cmeg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cmfdqrg(i,k) = cmfdqrg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
!--bee
         evpg(i,k) = evpg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         du2(i,k) = du2(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         eu2(i,k) = eu2(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         ed2(i,k) = ed2(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
      end do
   end do

   call closure(lchnk   , &
                qg      ,tg      ,pg      ,zg      ,sg      , &
                tpg     ,qs      ,qu      ,su      ,mc      , &
                du      ,mu      ,md      ,qd      ,sd      , &
                qhat    ,shat    ,dp      ,qstpg   ,zfg     , &
                qlg     ,dsubcld ,mb      ,capeg   ,tlg     , &
                lclg    ,lelg    ,jt      ,maxg    ,1       , &
                lengath ,rgas    ,grav    ,cpres   ,rl      , &
                msg     ,capelmt ,nstep   ,ptendqg ,ptendtg , &
                negcapeg)
!
! limit cloud base mass flux to theoretical upper bound.
!
   do i=1,lengath
      mumax(i) = 0
   end do
   do k=msg + 2,pver
      do i=1,lengath
         mumax(i) = max(mumax(i), mu(i,k)/dp(i,k))
      end do
   end do
   do i=1,lengath
      if (mumax(i) > 0.) then
         mb(i) = min(mb(i),0.5/(delt*mumax(i)))
      else
         mb(i) = 0.
      endif
   end do
   do k=msg+1,pver
      do i=1,lengath
         mu(i,k) = mu(i,k)*mb(i)
         md(i,k) = md(i,k)*mb(i)
         mc(i,k) = mc(i,k)*mb(i)
         du(i,k) = du(i,k)*mb(i)
         eu(i,k) = eu(i,k)*mb(i)
         ed(i,k) = ed(i,k)*mb(i)
         cmeg(i,k) = cmeg(i,k)*mb(i)
!++bee
         cmfdqrg(i,k) = cmfdqrg(i,k)*mb(i)
!--bee
         cug(i,k) = cug(i,k)*mb(i)
         evpg(i,k) = evpg(i,k)*mb(i)
         pflxg(i,k+1) = pflxg(i,k+1)*mb(i)*100./grav
         mu2(i,k) = mu2(i,k)*mb(i)
         md2(i,k) = md2(i,k)*mb(i)
         du2(i,k) = du2(i,k)*mb(i)
         eu2(i,k) = eu2(i,k)*mb(i)
         ed2(i,k) = ed2(i,k)*mb(i)
      end do
   end do
   do i = 1,lengath
!
! totpcp from cldprp has the dimension of kg/kg, here it is
! converted to kg/(m^2*s), the precipitation rate
!
      totpcp(i) = totpcp(i)*mb(i)*100./grav
      totevp(i) = totevp(i)*mb(i)*100./grav
   end do
!
! compute temperature and moisture changes due to convection.
!
   call q1q2_pjr(lchnk   , &
                 dqdt    ,dsdt    ,qg      ,qs      ,qu      , &
                 su      ,du      ,qhat    ,shat    ,dp      , &
                 mu      ,md      ,sd      ,qd      ,qlg     , &
                 dsubcld ,jt      ,maxg    ,1       ,lengath , &
                 cpres   ,rl      ,msg     ,nstep   ,          &
                 dlg     ,evpg    ,cug     )
!
! compute momentum changes due to convection, if desired (i.e
! if logical switch set).
!
!
!      if(momentm)
!       then
!        call moment(dudt,dvdt,du,alpha,dp,ed,eu,mc,md,mu,
!     1             pg,qd,qu,qhat,sd,su,shat,ud,vd,tg,ug,vg,zg,zfg,
!     2             dsubcld,maxg,jd,jt,rl,
!     3             msg,2.*delt,grav,cpres,rgas,pver,1,lengath,nlond,lat(1))
!      endif
!
! move the convective transport outside of conv_cam.
!
! gather back temperature and mixing ratio.
!
   do k = msg + 1,pver
      do i = 1,lengath
         psdiss = psdiss + (dudt(i,k)*u(ideep(i),k)+dvdt(i,k)*v(ideep(i),k))* &
                  dpp(ideep(i),k)/grav
!
! q is updated to compute net precip, and then reset to old value.
! the last line is overwritten. so the input basic variables, i.e.
! q, t, u and v are updated to include convective increments.
! (5/2/95)
!
         q(ideep(i),k) = q(ideep(i),k) + 2.*delt*dqdt(i,k)
!**BAB don't update t (real state variable)
!!$         t(ideep(i),k) = t(ideep(i),k) + 2.*delt*dsdt(i,k)
         u(ideep(i),k) = u(ideep(i),k) + 2.*delt*dudt(i,k)
         v(ideep(i),k) = v(ideep(i),k) + 2.*delt*dvdt(i,k)
         cme(ideep(i),k) = cmeg(i,k)
!++bee
         cmfdqr(ideep(i),k) = cmfdqrg(i,k)
!--bee
         zdu(ideep(i),k) = du2(i,k)
         mcon(ideep(i),k) = mc(i,k)
!**BAB not this tendency
!!$         qtg(ideep(i),k) = dqdt(i,k)
!**BAB report back heating in J/kg/s
!!$         ttg(ideep(i),k) = dsdt(i,k)
         ttg(ideep(i),k) = dsdt(i,k)*cpres
!            utg(ideep(i),k) = dudt(i,k)
!            vtg(ideep(i),k) = dvdt(i,k)
         dlf(ideep(i),k) = dlg(i,k)
         pflx(ideep(i),k) = pflxg(i,k)
         ql(ideep(i),k) = qlg(i,k)
      end do
   end do
!
   do i = 1,lengath
      jctop(ideep(i)) = jt(i)
!++bee
      jcbot(ideep(i)) = maxg(i)
!--bee
      pflx(ideep(i),pverp) = pflxg(i,pverp)
      psevap = psevap + totevp(i)
      psrain = psrain + totpcp(i)
   end do
!
! convert back to specific humidity from mixing ratio.
! take into account any moisture added to ensure positiveness
! of specific humidity at start of routine.
!
   do k = msg + 1,pver
      do i = 1,ncol
         q(i,k) = q(i,k) - max((qmin-qh(i,k,1)),0._r8)
      end do
   end do
   do k = pver,msg + 1,-1
      do i = 1,ncol
         sumq(i) = sumq(i) - dpp(i,k)* (q(i,k)-qh(i,k,1))
         qtg(i,k) = (q(i,k)-qh(i,k,1)) / (2.*delt)
!
! account for the detraining cloud water in the precip
!
         sumq(i) = sumq(i) - dpp(i,k)*dlf(i,k)*2*delt
         pcpck(i,k) = max(0._r8,sumq(i))
      end do
   end do
!
! obtain final precipitation rate.
!
   do i = 1,ncol
!         llo1 = ts(i) .ge. tfreez
!
! here pcpr and pcps are in units of kg/m^2, ie. precip per
! time step
!
!         pcpr(i) = cvmgt(pcpdh(i)*max(sumq(i),0._r8),0.,llo1)
!         pcps(i) = cvmgt(0.,pcpdh(i)*max(sumq(i),0._r8),llo1)
      precip = pcpdh(i)*max(sumq(i),0._r8)
      if (ts(i) >= tfreez) then
         pcpr(i) = precip
         pcps(i) = 0.
      else
         pcpr(i) = 0.
         pcps(i) = precip
      end if
   end do

!
! accumulate precipitation, the 1000. is the density of water, so
! paprc and paprs are now in units of meters.
!
   do i = 1,ncol
      paprc(i) = paprc(i) + (pcpr(i)+pcps(i))/1000.
      paprs(i) = paprs(i) + pcps(i)/1000.
   end do
!
! convert precipitation to m/s, ie, precip rate.
!
   do i = 1,ncol
      pcpr(i) = pcpr(i)/ (2.*delt)/1000.
      pcps(i) = pcps(i)/ (2.*delt)/1000.
      pcpc(i) = pcpr(i) + pcps(i)
      psheat = psheat + (pcps(i)+pcpr(i))*rl
   end do
   do k = msg + 1,pver
      do i = 1,ncol
         pcpck(i,k) = pcpdh(i)*pcpck(i,k)/ (2.*delt)
      end do
   end do
!
! calculate conservation of quantities.
!

!       if(lat(1).eq.jlatpr)then
!        do l=msg+1,pver
!        do i=1,lengath
!          sumdq(i) = sumdq(i) + 2.*delt*(rl/cpres)*dpp(ideep(i),l)*
!     1                          dqdt(i,l)
!          sumdt(i) = sumdt(i) + 2.*delt*dpp(ideep(i),l)*dsdt(i,l)
!        end do
!        end do
!
!        write(6,*)'sumdq,sumdt,sumde in convection subroutine########'
!        do i=1,lengath
!          sumde(i) = sumdt(i) + sumdq(i)
!          write(6, 901) sumdq(i), sumdt(i),sumde(i), i, ideep(i)
!        end do
!c
!        write(6,*)'sumdq,sumdt,sumde ... all points'
!      do i=1,ncol
!         sumdq(i) = 0.0
!         sumdt(i) = 0.0
!      end do
!c
!      do l=msg+1,pver
!      do i=1,ncol
!        deltaq(i,l) = qh(i,l,1) - deltaq(i,l)
!        deltat(i,l) = t (i,l) - deltat(i,l)
!        sumdq(i) = sumdq(i) + (rl/cpres)*dpp(i,l)*deltaq(i,l)
!        sumdt(i) = sumdt(i) + dpp(i,l)*deltat(i,l)
!      end do
!      end do
!      do i=1,ncol
!        sumde(i) = sumdt(i) + sumdq(i)
!      end do
!        write(6, 902) (i,sumdq(i),sumdt(i),sumde(i),i=1,ncol)
!
!  901   format(1x,3e20.12, i10, i10)
!  902   format(1x,i10, 3e20.12)
!
!      endif

   return
end subroutine zm_convr

!===============================================================================
  subroutine zm_conv_evap(state, ptend, pflx, precc, cldo, deltat, evappct)

!-----------------------------------------------------------------------
! compute tendencies due to evaporation of rain from ZM scheme
!-----------------------------------------------------------------------

    use dycore,         only: dycore_is
    use physics_types,  only: physics_state, physics_ptend
    use wv_saturation,  only: aqsat

#include <guang.h>
!------------------------------Arguments--------------------------------
    real(r8), intent(in   ) :: cldo(pcols,pver)   ! old cloud fraction
    real(r8), intent(in   ) :: deltat             ! time step
    type(physics_state), intent(in)    :: state   ! Physics state variables
    type(physics_ptend), intent(inout) :: ptend   ! indivdual parameterization tendencies

    real(r8), intent(inout) :: precc(pcols)       ! Convective-scale preciptn rate
    real(r8), intent(inout) :: pflx (pcols,pverp) ! Conv rain flux thru out btm of lev
    real(r8), intent(out)   :: evappct(pcols)     ! Convective-scale preciptn rate
!
!---------------------------Local storage-------------------------------

    real(r8) :: dpovrg                 ! Pressure depth over gravity
    real(r8) :: envevap                ! Evaporation rate from precitation
    real(r8) :: est (pcols,pver)       ! Saturation vapor pressure
    real(r8) :: ke                     ! Tunable evaporation efficiency
    real(r8) :: pflxtmp(pcols,pver)    ! Conv rain flux thru out btm of lev (temp)
    real(r8) :: qsat(pcols,pver)       ! Specific humidity at Saturation
    real(r8) :: rh                     ! Relative humidity
    real(r8) :: sumflx(pcols)          ! flux integral

    integer :: i,k                     ! longitude,level indices

!-----------------------------------------------------------------------
!
! Evaporate some of the precip directly into the environment (Sundqvist)
! pflx = kg/m^2/s
!
    call aqsat (state%t    ,state%pmid  ,est    ,qsat    ,pcols   , &
         state%ncol ,pver  ,1       ,pver    )

    sumflx(:) = 0.

    if ( dycore_is ('LR') ) then
!====== finite-volume dynamics =============================
       ke = 7.5E-6     ! larger value to compensate for cloud fraction effect
       sumflx = 0.
       pflx(:,1) = 0.

       do k=1,pver
          do i=1,state%ncol
             if(precc(i) > 0.001*sumflx(i) .and. state%t(i,k) > tfreez) then
                dpovrg  =  state%pdel(i,k)*rgrav
                pflxtmp(i,k) = max(pflx(i,k)-sumflx(i),0._r8)
                envevap = max(ke*(1.-state%q(i,k,1)/qsat(i,k))*sqrt(pflxtmp(i,k)),0._r8)

! Sundqvist's Cloud fraction effect

                envevap = envevap * (1. - cldo(i,k))
                envevap = max(min(envevap, 0.75*(qsat(i,k)-state%q(i,k,1))/deltat),0._r8)
                envevap = min(pflxtmp(i,k)/dpovrg,envevap)

! The factor 0.999 is to prevent rounding-level negative precc

                envevap = min(envevap,0.999*(1000.*precc(i)-sumflx(i))/dpovrg)
                envevap = min(envevap,max(0._r8,0.999*(pflx(i,k)-sumflx(i))/dpovrg))
                pflx(i,k) = pflx(i,k) - envevap*dpovrg - sumflx(i)
                sumflx(i) = sumflx(i) + envevap*dpovrg
                ptend%q(i,k,1) =  envevap
                ptend%s(i,k)   = -envevap*rl
             else
                ptend%q(i,k,1) =  0.
                ptend%s(i,k)   =  0.
              endif
          end do
       end do
    else
!====== EUL/SLD dynamics =============================
        ke = 7.5E-6
!       ke = 9.0E-6 ! modified by LI Lijuan
       do k=1,pver
          do i=1,state%ncol
             dpovrg  =  state%pdel(i,k)*rgrav
             rh      = state%q(i,k,1)/qsat(i,k)
             pflxtmp(i,k)= max(pflx(i,k)-sumflx(i),0._r8)

! Determine evap amount based on rh can't be negative

             envevap   = max(ke*(1.0 - rh)*sqrt(pflxtmp(i,k)),0._r8)

! Don't let envevap supersaturate layer (approx)

             envevap   = max(min(envevap, (qsat(i,k)-state%q(i,k,1))/deltat),0._r8)

! Can't evap more than what we have at this level

             envevap   = min(pflxtmp(i,k)/dpovrg,envevap)

! Can't evap more than what we are told is in bottom level
! precc and pflx(,pver) are supposed to be identical but they differ by
! roundoff, precc is used for this check since that is what is written out

             envevap  = min((precc(i)*1.E3-sumflx(i))/dpovrg,envevap)
             pflx(i,k) = pflx(i,k) - envevap*dpovrg - sumflx(i)
             sumflx(i) = sumflx(i)+envevap*dpovrg
             ptend%q(i,k,1) =  envevap
             ptend%s(i,k)   = -envevap*rl
          end do
       end do
    endif
!
! adjust precc by the amount of precip evaporated
!
    evappct(:)=0.
    do i = 1, state%ncol
       precc(i)=precc(i)-(sumflx(i)/1.E3)
       if (precc(i) > 0) evappct(i)=(sumflx(i)/1.E3)/precc(i)
    end do

  end subroutine zm_conv_evap

    subroutine buoyan(lchnk,  ncol,   &
                      q,      t,      p,    z,    pf,    &
                      tp,     qstp,   tl,   rl,   cape,  &
                      pblt,   lcl,    lel,  lon,  mx,    &
                      rd,     grav,   cp,   msg,  nstep, &
                      tpert,  negcape)
        !-----------------------------------------------------------------------
        !
        ! Purpose:
        ! <Say what the routine does>
        !
        ! Method:
        ! <Describe the algorithm(s) used in the routine.>
        ! <Also include any applicable external references.>
        !
        ! Author:
        ! This is contributed code not fully standardized by the CCM core group.
        ! The documentation has been enhanced to the degree that we are able.
        ! Reviewed:          P. Rasch, April 1996
        !
        !-----------------------------------------------------------------------
        implicit none
#include <guang.h>

        ! input arguments
        integer, intent(in) :: lchnk              ! chunk identifier
        integer, intent(in) :: ncol               ! number of atmospheric columns

        real(r8), intent(in) :: q(pcols,pver)     ! specific humidity
        real(r8), intent(in) :: t(pcols,pver)     ! temperature
        real(r8), intent(in) :: p(pcols,pver)     ! pressure
        real(r8), intent(in) :: z(pcols,pver)     ! height
        real(r8), intent(in) :: pf(pcols,pver+1)  ! pressure at interfaces
        real(r8), intent(in) :: pblt(pcols)       ! index of pbl depth
        real(r8), intent(in) :: tpert(pcols)      ! perturbation temperature by PBL processes

        real(r8), intent(in) :: rl
        real(r8), intent(in) :: rd
        real(r8), intent(in) :: cp
        real(r8), intent(in) :: grav

        integer, intent(in) :: msg
        integer, intent(in) :: nstep

        ! output arguments
        real(r8), intent(out) :: tp(pcols,pver)   ! parcel temperature
        real(r8), intent(out) :: qstp(pcols,pver) ! saturation mixing ratio of parcel
        real(r8), intent(out) :: tl(pcols)        ! parcel temperature at lcl
        real(r8), intent(out) :: cape(pcols)      ! convective available potential energy
        real(r8), intent(out) :: negcape(pcols)   ! negative convective available potential energy
        integer lcl(pcols)                        ! lifting  condensation level
        integer lel(pcols)                        !
        integer, intent(out) :: lon(pcols)        ! level of onset of deep convection
        integer, intent(out) :: mx(pcols)         ! level of max moist static energy
!
!--------------------------Local Variables------------------------------
!
        real(r8) capeten(pcols,5)                 ! provisional value of cape
        real(r8) tv(pcols,pver)                   !
        real(r8) tpv(pcols,pver)                  !
        real(r8) buoy(pcols,pver)
        real(r8) qs(pcols,pver)                   ! saturation specific humidity

        real(r8) a1(pcols)
        real(r8) a2(pcols)
        real(r8) estp(pcols)
        real(r8) pl(pcols)
        real(r8) plexp(pcols)
        real(r8) hmax(pcols)
        !real(r8) hmax(pcols)      ! DONG Li: Is it the same with the above one?
        real(r8) hmin(pcols)
        !real(r8) hmn(pcols)       ! DONG Li: Who make this?
        real(r8) hmn(pcols,pver)   ! moist static energy of environment
        real(r8) y(pcols)

        logical plge600(pcols)
        integer knt(pcols)
        integer lelten(pcols,5)
        integer khmin(pcols)       ! level of min moist static energy

        real(r8) e

        integer i
        integer k
        integer n
        integer idx0

        real(r8) esttmp
#ifdef PERGRO
        real(r8) rhd
#endif

        !-----------------------------------------------------------------------
        ! initialize the variables
        do n = 1, 5
            do i = 1, ncol
                lelten(i,n) = pver
                capeten(i,n) = 0.
            end do
        end do
        do i = 1, ncol
            lon(i)      = pver
            knt(i)      = 0
            lel(i)      = pver
            mx(i)       = lon(i)
            khmin(i)    = pver
            cape(i)     = 0.
            negcape(i)  = 0.
            hmax(i)     = 0.
            hmin(i)     = 1.e10
        end do
        do i = 1, pcols
            do k = pver, 1, -1
                buoy(i,k) = 0.
            end do
        end do

        tp(:ncol,:)   = t(:ncol,:) ! DONG Li: The difference between "ncol" and "pcols"
        qstp(:ncol,:) = q(:ncol,:)

        !-----------------------------------------------------------------------
        ! DONG Li: What we do here?
        ! set "launching" level(mx) to be at maximum moist static energy.
        ! search for this level (stop at planetary boundary layer top).

        ! find level of max static energy below h_min level
#ifdef PERGRO
        do k = pver, msg+1, -1
            do i = 1, ncol
                hmn(i,k) = cp*t(i,k)+grav*z(i,k)+rl*q(i,k)
                ! reset max moist static energy level
                ! when relative difference exceeds 1.e-4
                ! DONG Li: Why "pblt" is real? And what is "msg"?
                rhd = (hmn(i,k)-hmax(i))/(hmn(i,k)+hmax(i))
                if (k >= nint(pblt(i)) .and. k <= lon(i) .and. rhd > -1.e-4) then
                    hmax(i) = hmn(i,k)
                    mx(i) = k
                end if
            end do
        end do
#else
        do k = pver, msg+1, -1
            do i = 1, ncol
                hmn(i,k) = cp*t(i,k)+grav*z(i,k)+rl*q(i,k)
                !hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
                if (k >= nint(pblt(i)) .and. k <= lon(i) .and. hmn(i,k) > hmax(i)) then
                    hmax(i) = hmn(i,k)
                    mx(i) = k
                end if
            end do
        end do
#endif

        ! find min hmn level above 700 mb
        do i = 1, ncol
            do k = msg+1, pver
                if (p(i,k) <= 700. .and. hmin(i) >= hmn(i,k) ) then
                    hmin(i)  = hmn(i,k)
                    khmin(i) = k
                end if
            end do
        end do

        do i = 1, ncol
            do k = pver, msg+1, -1
                ! DONG Li: calculate saturated moisture vapour pressure?
                esttmp = c1*exp((c2*(t(i,k)-tfreez))/((t(i,k)-tfreez)+c3))
                ! DONG Li: What is "qs"? Saturated moisture mixing ratio?
                if (p(i,k)-esttmp > 0.) then
                    ! DONG Li: unsaturated
                    qs(i,k) = eps1*esttmp/(p(i,k)-esttmp)
                else
                    qs(i,k) = 1.0
                end if
                ! rh threshold does not apply to midlevel convection
                if (k > khmin(i) .and. k < nint(pblt(i)) .and. &
                    q(i,k)/qs(i,k) >= 0.6 .and. hmn(i,k) > hmax(i)) then
                    hmax(i) = hmn(i,k)
                    mx(i) = k
                end if
            end do
        end do

        do i = 1, ncol
            lcl(i) = mx(i)
            e = p(i,mx(i))*q(i,mx(i))/(eps1+q(i,mx(i)))
            tl(i) = 2840./(3.5*log(t(i,mx(i)))-log(e)-4.805)+55.
            if (tl(i) < t(i,mx(i))) then
                plexp(i) = (1./(0.2854*(1.-0.28*q(i,mx(i)))))
                pl(i) = p(i,mx(i))*(tl(i)/t(i,mx(i)))**plexp(i)
            else
                tl(i) = t(i,mx(i))
                pl(i) = p(i,mx(i))
            end if
        end do

        ! calculate lifting condensation level (lcl).
        do k = pver, msg+2, -1
            do i = 1, ncol
                if (k <= mx(i) .and. (p(i,k) > pl(i) .and. p(i,k-1) <= pl(i))) then
                    lcl(i) = k-1
                end if
            end do
        end do

        ! if lcl is above the nominal level of non-divergence (600 mbs),
        ! no deep convection is permitted (ensuing calculations
        ! skipped and cape retains initialized value of zero).
        do i = 1, ncol
            plge600(i) = pl(i).ge.600.
        end do

        ! initialize parcel properties in sub-cloud layer below lcl.
        do k = pver, msg+1, -1
            do i = 1, ncol
                if (k > lcl(i) .and. k <= mx(i) .and. plge600(i)) then
                    tv(i,k) = t(i,k)*(1.+1.608*q(i,k))/(1.+q(i,k))
                    qstp(i,k) = q(i,mx(i))
                    tp(i,k) = t(i,mx(i))*(p(i,k)/p(i,mx(i)))**(0.2854*(1.-0.28*q(i,mx(i))))
!
! buoyancy is increased by 0.5 k as in tiedtke
!
!-jjh          tpv (i,k)=tp(i,k)*(1.+1.608*q(i,mx(i)))/
!-jjh     1                     (1.+q(i,mx(i)))
                    tpv(i,k) = (tp(i,k)+tpert(i))*(1.+0.608*qstp(i,k))
                    buoy(i,k) = tpv(i,k)-tv(i,k)+0.5
                end if
            end do
        end do
!
! define parcel properties at lcl (i.e. level immediately above pl).
!
        do k = pver, msg+1, -1
            do i = 1, ncol
                if (k == lcl(i) .and. plge600(i)) then
                    tv(i,k) = t(i,k)* (1.+1.608*q(i,k))/ (1.+q(i,k))
                    qstp(i,k) = q(i,mx(i))
                    tp(i,k) = tl(i)*(p(i,k)/pl(i))**(0.2854*(1.-0.28*qstp(i,k)))
!              estp(i)  =exp(a-b/tp(i,k))
! use of different formulas for est has about 1 g/kg difference
! in qs at t= 300k, and 0.02 g/kg at t=263k, with the formula
! above giving larger qs.
!
                    estp(i) = c1*exp((c2*(tp(i,k)-tfreez))/((tp(i,k)-tfreez)+c3))
                    qstp(i,k) = eps1*estp(i)/(p(i,k)-estp(i))
                    a1(i) = cp/rl+qstp(i,k)*(1.+ qstp(i,k)/eps1)*rl*eps1/ &
                        (rd*tp(i,k)**2)
                    a2(i) = .5*(qstp(i,k)*(1.+2./eps1*qstp(i,k))* &
                        (1.+qstp(i,k)/eps1)*eps1**2*rl*rl/ &
                        (rd**2*tp(i,k)**4)-qstp(i,k)* &
                        (1.+qstp(i,k)/eps1)*2.*eps1*rl/ &
                        (rd*tp(i,k)**3))
                    a1(i) = 1./a1(i)
                    a2(i) = -a2(i)*a1(i)**3
                    y(i) = q(i,mx(i))-qstp(i,k)
                    tp(i,k) = tp(i,k)+a1(i)*y(i)+a2(i)*y(i)**2
!          estp(i)  =exp(a-b/tp(i,k))
                    estp(i) = c1*exp((c2*(tp(i,k)-tfreez))/((tp(i,k)-tfreez)+c3))
                    qstp(i,k) = eps1*estp(i)/(p(i,k)-estp(i))
!
! buoyancy is increased by 0.5 k in cape calculation.
! dec. 9, 1994
!-jjh          tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/(1.+q(i,mx(i)))
!
                    tpv(i,k) = (tp(i,k)+tpert(i))*(1.+0.608*qstp(i,k))
                    buoy(i,k) = tpv(i,k)-tv(i,k)+0.5
                end if
            end do
        end do
!
! main buoyancy calculation.
!
   do k = pver - 1,msg + 1,-1
      do i=1,ncol
         if (k < lcl(i) .and. plge600(i)) then
            tv(i,k) = t(i,k)* (1.+1.608*q(i,k))/ (1.+q(i,k))
            qstp(i,k) = qstp(i,k+1)
            tp(i,k) = tp(i,k+1)* (p(i,k)/p(i,k+1))**(0.2854* (1.-0.28*qstp(i,k)))
!          estp(i) = exp(a-b/tp(i,k))
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/((tp(i,k)-tfreez)+c3))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
            a1(i) = cp/rl + qstp(i,k)* (1.+qstp(i,k)/eps1)*rl*eps1/ (rd*tp(i,k)**2)
            a2(i) = .5* (qstp(i,k)* (1.+2./eps1*qstp(i,k))* &
                    (1.+qstp(i,k)/eps1)*eps1**2*rl*rl/ &
                    (rd**2*tp(i,k)**4)-qstp(i,k)* &
                    (1.+qstp(i,k)/eps1)*2.*eps1*rl/ &
                    (rd*tp(i,k)**3))
            a1(i) = 1./a1(i)
            a2(i) = -a2(i)*a1(i)**3
            y(i) = qstp(i,k+1) - qstp(i,k)
            tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
!          estp(i)  =exp(a-b/tp(i,k))
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/ ((tp(i,k)-tfreez)+c3))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
!-jjh          tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/
!jt            (1.+q(i,mx(i)))
            tpv(i,k) = (tp(i,k)+tpert(i))* (1.+0.608*qstp(i,k))
            buoy(i,k) = tpv(i,k) - tv(i,k) + 0.5
         end if
      end do
   end do

        ! DONG Li: Between LCL and 600mb, divide the column into AT MOST 5 parts
        do k = msg+2, pver
            do i = 1, ncol
                if (k < lcl(i) .and. plge600(i)) then
                    if (buoy(i,k+1) > 0. .and. buoy(i,k) <= 0.) then
                        knt(i) = min(5, knt(i)+1)
                        lelten(i,knt(i)) = k
                    end if
                end if
            end do
        end do

        ! calculate negative CAPE
        do i = 1, pcols
            negcape(i) = 0.
            idx0 = pver
            do k = pver-1, 10, -1 ! DONG Li: Where does "10" come from? Is it OK?
                if (plge600(i) .and. &
                    (buoy(i,k+1) <= 0. .and. &
                     buoy(i,k)   >  0. .and. &
                     buoy(i,k-1) >  0. .and. &
                     buoy(i,k-2) >  0.)) then
                    idx0 = k+1
                    goto 5
                end if
            end do
5           continue
            do k = pver-1, idx0, -1
                if (plge600(i)) then
                    negcape(i) = negcape(i)+rd*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
                end if
            end do
        end do

        ! calculate CAPE
        do n = 1, 5
            do k = msg+1, pver
                do i = 1, ncol
                    if (plge600(i) .and. &
                        k <= mx(i) .and. &
                        k > lelten(i,n)) then
                        capeten(i,n) = capeten(i,n)+rd*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
                    end if
                end do
            end do
        end do

        ! find maximum cape from all possible tentative capes from one sounding,
        ! and use it as the final cape, april 26, 1995
        do n = 1, 5
            do i = 1, ncol
                if (capeten(i,n) > cape(i)) then
                    cape(i) = capeten(i,n)
                    lel(i) = lelten(i,n)
                end if
            end do
        end do

        ! put lower bound on cape for diagnostic purposes.
        do i = 1, ncol
            cape(i) = max(cape(i), 0._r8)
        end do

        return
    end subroutine buoyan

subroutine cldprp(lchnk   , &
                  q       ,t       ,u       ,v       ,p       , &
                  z       ,s       ,mu      ,eu      ,du      , &
                  md      ,ed      ,sd      ,qd      ,mc      , &
                  qu      ,su      ,zf      ,qst     ,hmn     , &
                  hsat    ,shat    ,ql      ,totpcp  ,totevp  , &
                  cmeg    ,jb      ,lel     ,jt      ,jlcl    , &
                  mx      ,j0      ,jd      ,rl      ,il2g    , &
                  rd      ,grav    ,cp      ,msg     ,nstep   , &
                           pflx    ,evp     ,cu      ,mu2     , &
                  eu2     ,du2     ,md2     ,ed2     ,cmfdqr  , &
                  limcnv  )
!-----------------------------------------------------------------------
!
! Purpose:
! <Say what the routine does>
!
! Method:
! may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
!             original version cldprop.
!
! Author: P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! there are debug statements left strewn about and code segments disabled
! these are to facilitate future development. We expect to release a
! cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
!
!**** PLEASE NOTE ****
!
! we are aware of a specific problem in this code
! (identified by the string ---> PROBLEM ONE)
! during the calculation of the updraft cloud properties,
! rather than adding a perturbation to the updraft temperature of
! half a degree, (there was an inadvertant addition of cp*0.5) degrees
! or about 500 degrees. (This problem was in the code prior to its
! contribution to the NCAR effort)

! Fortunately, the erroneous values
! are overwritten later in the code. The problem is quite subtle.
! The erroneous values would persist between cloud base and the lifting
! condensation level. The addition of the very high perturbation to the updraft
! temperature causes the saturation mixing ratio to be set to zero,
! and later the lcl to be set to one level above cloud base.
! There are therefore no levels between cloud base and the lcl. Therefore
! all erroneous values are overwritten.

! The only manifestation we are aware of with respect to this problem
! is that the lifting condensation level is constrained to be one level above
! cloud base.

! We discovered the problem after too much had been invested in
! very long integrations (in terms of computer time)
! to allow for a modification and model retuning. It is our expectation that
! this problem will be fixed in the next release of the model.
!
!-----------------------------------------------------------------------

   implicit none

#include <guang.h>
!------------------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                  ! chunk identifier

   real(r8), intent(in) :: q(pcols,pver)         ! spec. humidity of env
   real(r8), intent(in) :: t(pcols,pver)         ! temp of env
   real(r8), intent(in) :: p(pcols,pver)         ! pressure of env
   real(r8), intent(in) :: z(pcols,pver)         ! height of env
   real(r8), intent(in) :: s(pcols,pver)         ! normalized dry static energy of env
   real(r8), intent(in) :: zf(pcols,pverp)       ! height of interfaces
   real(r8), intent(in) :: u(pcols,pver)         ! zonal velocity of env
   real(r8), intent(in) :: v(pcols,pver)         ! merid. velocity of env

   integer, intent(in) :: jb(pcols)              ! updraft base level
   integer, intent(in) :: lel(pcols)             ! updraft launch level
   integer, intent(out) :: jt(pcols)              ! updraft plume top
   integer, intent(out) :: jlcl(pcols)            ! updraft lifting cond level
   integer, intent(in) :: mx(pcols)              ! updraft base level (same is jb)
   integer, intent(out) :: j0(pcols)              ! level where updraft begins detraining
   integer, intent(out) :: jd(pcols)              ! level of downdraft
   integer, intent(in) :: limcnv                 ! convection limiting level
   integer, intent(in) :: il2g                   !CORE GROUP REMOVE
   integer, intent(in) :: msg                    ! missing moisture vals (always 0)
   integer, intent(in) :: nstep                  ! time step index
   real(r8), intent(in) :: rl                    ! latent heat of vap
   real(r8), intent(in) :: shat(pcols,pver)      ! interface values of dry stat energy
!
! output
!
   real(r8), intent(out) :: cmfdqr(pcols,pver)   ! rate of production of precip at that layer
   real(r8), intent(out) :: du(pcols,pver)       ! detrainement rate of updraft
   real(r8), intent(out) :: ed(pcols,pver)       ! entrainment rate of downdraft
   real(r8), intent(out) :: eu(pcols,pver)       ! entrainment rate of updraft
   real(r8), intent(out) :: hmn(pcols,pver)      ! moist stat energy of env
   real(r8), intent(out) :: hsat(pcols,pver)     ! sat moist stat energy of env
   real(r8), intent(out) :: mc(pcols,pver)       ! net mass flux
   real(r8), intent(out) :: md(pcols,pver)       ! downdraft mass flux
   real(r8), intent(out) :: mu(pcols,pver)       ! updraft mass flux
   real(r8), intent(out) :: pflx(pcols,pverp)    ! precipitation flux thru layer
   real(r8), intent(out) :: qd(pcols,pver)       ! spec humidity of downdraft
   real(r8), intent(out) :: ql(pcols,pver)       ! liq water of updraft
   real(r8), intent(out) :: qst(pcols,pver)      ! saturation spec humidity of env.
   real(r8), intent(out) :: qu(pcols,pver)       ! spec hum of updraft
   real(r8), intent(out) :: sd(pcols,pver)       ! normalized dry stat energy of downdraft
   real(r8), intent(out) :: su(pcols,pver)       ! normalized dry stat energy of updraft
!
!     these version of the mass fluxes conserve mass (used in tracer transport)
!
   real(r8), intent(out) :: mu2(pcols,pver)      ! updraft mass flux
   real(r8), intent(out) :: eu2(pcols,pver)      ! updraft entrainment
   real(r8), intent(out) :: du2(pcols,pver)      ! updraft detrainment
   real(r8), intent(out) :: md2(pcols,pver)      ! downdraft mass flux
   real(r8), intent(out) :: ed2(pcols,pver)      ! downdraft entrainment


   real(r8) rd                   ! gas constant for dry air
   real(r8) grav                 ! gravity
   real(r8) cp                   ! heat capacity of dry air

!
! Local workspace
!
   real(r8) gamma(pcols,pver)
   real(r8) dz(pcols,pver)
   real(r8) iprm(pcols,pver)
   real(r8) hu(pcols,pver)
   real(r8) hd(pcols,pver)
   real(r8) eps(pcols,pver)
   real(r8) f(pcols,pver)
   real(r8) k1(pcols,pver)
   real(r8) i2(pcols,pver)
   real(r8) ihat(pcols,pver)
   real(r8) i3(pcols,pver)
   real(r8) idag(pcols,pver)
   real(r8) i4(pcols,pver)
   real(r8) qsthat(pcols,pver)
   real(r8) hsthat(pcols,pver)
   real(r8) gamhat(pcols,pver)
   real(r8) cu(pcols,pver)
   real(r8) evp(pcols,pver)
   real(r8) cmeg(pcols,pver)
   real(r8) qds(pcols,pver)
   real(r8) hmin(pcols)
   real(r8) expdif(pcols)
   real(r8) expnum(pcols)
   real(r8) ftemp(pcols)
   real(r8) eps0(pcols)
   real(r8) rmue(pcols)
   real(r8) zuef(pcols)
   real(r8) zdef(pcols)
   real(r8) epsm(pcols)
   real(r8) ratmjb(pcols)
   real(r8) est(pcols)
   real(r8) totpcp(pcols)
   real(r8) totevp(pcols)
   real(r8) alfa(pcols)
   real(r8) beta
   real(r8) c0
   real(r8) ql1
   real(r8) weight
   real(r8) tu
   real(r8) estu
   real(r8) qstu

   real(r8) small
   real(r8) mdt
!      real(r8) cu2

   integer khighest
   integer klowest
   integer kount
   integer i,k

   logical doit(pcols)
   logical done(pcols)
!
!------------------------------------------------------------------------------
!
   do i = 1,il2g
      ftemp(i) = 0.
      expnum(i) = 0.
      expdif(i) = 0.
   end do
!
!jr Change from msg+1 to 1 to prevent blowup
!
   do k = 1,pver
      do i = 1,il2g
         dz(i,k) = zf(i,k) - zf(i,k+1)
      end do
   end do

!
! initialize many output and work variables to zero
!
   pflx(:il2g,1) = 0

   do k = msg + 1,pver
      do i = 1,il2g
         k1(i,k) = 0.
         i2(i,k) = 0.
         i3(i,k) = 0.
         i4(i,k) = 0.
         mu(i,k) = 0.
         f(i,k) = 0.
         eps(i,k) = 0.
         eu(i,k) = 0.
         du(i,k) = 0.
         ql(i,k) = 0.
         cu(i,k) = 0.
         evp(i,k) = 0.
         cmeg(i,k) = 0.
         qds(i,k) = q(i,k)
         md(i,k) = 0.
         ed(i,k) = 0.
         sd(i,k) = s(i,k)
         qd(i,k) = q(i,k)
         mc(i,k) = 0.
         qu(i,k) = q(i,k)
         su(i,k) = s(i,k)
!        est(i)=exp(a-b/t(i,k))
         est(i) = c1*exp((c2* (t(i,k)-tfreez))/((t(i,k)-tfreez)+c3))
!++bee
         if ( p(i,k)-est(i) > 0. ) then
            qst(i,k) = eps1*est(i)/ (p(i,k)-est(i))
         else
            qst(i,k) = 1.0
         end if
!--bee
         gamma(i,k) = qst(i,k)*(1. + qst(i,k)/eps1)*eps1*rl/(rd*t(i,k)**2)*rl/cp
         hmn(i,k) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
         hsat(i,k) = cp*t(i,k) + grav*z(i,k) + rl*qst(i,k)
         hu(i,k) = hmn(i,k)
         hd(i,k) = hmn(i,k)
         mu2(i,k) = 0.
         eu2(i,k) = 0.
         du2(i,k) = 0.
         md2(i,k) = 0.
         ed2(i,k) = 0.
         pflx(i,k) = 0.
         cmfdqr(i,k) = 0.
      end do
   end do
!
!jr Set to zero things which make this routine blow up
!
   do k=1,msg
      do i=1,il2g
         cmfdqr(i,k) = 0.
         mu2(i,k) = 0.
         eu2(i,k) = 0.
         du2(i,k) = 0.
         md2(i,k) = 0.
         ed2(i,k) = 0.
      end do
   end do
!
! interpolate the layer values of qst, hsat and gamma to
! layer interfaces
!
   do i = 1,il2g
      hsthat(i,msg+1) = hsat(i,msg+1)
      qsthat(i,msg+1) = qst(i,msg+1)
      gamhat(i,msg+1) = gamma(i,msg+1)
      totpcp(i) = 0.
      totevp(i) = 0.
   end do
   do k = msg + 2,pver
      do i = 1,il2g
         if (abs(qst(i,k-1)-qst(i,k)) > 1.E-6) then
            qsthat(i,k) = log(qst(i,k-1)/qst(i,k))*qst(i,k-1)*qst(i,k)/ (qst(i,k-1)-qst(i,k))
         else
            qsthat(i,k) = qst(i,k)
         end if
         hsthat(i,k) = cp*shat(i,k) + rl*qsthat(i,k)
         if (abs(gamma(i,k-1)-gamma(i,k)) > 1.E-6) then
            gamhat(i,k) = log(gamma(i,k-1)/gamma(i,k))*gamma(i,k-1)*gamma(i,k)/ &
                                (gamma(i,k-1)-gamma(i,k))
         else
            gamhat(i,k) = gamma(i,k)
         end if
      end do
   end do
!
! initialize cloud top to highest plume top.
!jr changed hard-wired 4 to limcnv+1 (not to exceed pver)
!
   do i = 1,il2g
      jt(i) = max(lel(i),limcnv+1)
      jt(i) = min(jt(i),pver)
      jd(i) = pver
      jlcl(i) = lel(i)
      hmin(i) = 1.E6
   end do
!
! find the level of minimum hsat, where detrainment starts
!
   do k = msg + 1,pver
      do i = 1,il2g
         if (hsat(i,k) <= hmin(i) .and. k >= jt(i) .and. k <= jb(i)) then
            hmin(i) = hsat(i,k)
            j0(i) = k
         end if
      end do
   end do
   do i = 1,il2g
      j0(i) = min(j0(i),jb(i)-2)
      j0(i) = max(j0(i),jt(i)+2)
!
! Fix from Guang Zhang to address out of bounds array reference
!
      j0(i) = min(j0(i),pver)
   end do
!
! Initialize certain arrays inside cloud
!
   do k = msg + 1,pver
      do i = 1,il2g
         if (k >= jt(i) .and. k <= jb(i)) then
            hu(i,k) = hmn(i,mx(i)) + cp*0.5
            su(i,k) = s(i,mx(i)) + 0.5
         end if
      end do
   end do
!
! *********************************************************
! compute taylor series for approximate eps(z) below
! *********************************************************
!
   do k = pver - 1,msg + 1,-1
      do i = 1,il2g
         if (k < jb(i) .and. k >= jt(i)) then
            k1(i,k) = k1(i,k+1) + (hmn(i,mx(i))-hmn(i,k))*dz(i,k)
            ihat(i,k) = 0.5* (k1(i,k+1)+k1(i,k))
            i2(i,k) = i2(i,k+1) + ihat(i,k)*dz(i,k)
            idag(i,k) = 0.5* (i2(i,k+1)+i2(i,k))
            i3(i,k) = i3(i,k+1) + idag(i,k)*dz(i,k)
            iprm(i,k) = 0.5* (i3(i,k+1)+i3(i,k))
            i4(i,k) = i4(i,k+1) + iprm(i,k)*dz(i,k)
         end if
      end do
   end do
!
! re-initialize hmin array for ensuing calculation.
!
   do i = 1,il2g
      hmin(i) = 1.E6
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if (k >= j0(i) .and. k <= jb(i) .and. hmn(i,k) <= hmin(i)) then
            hmin(i) = hmn(i,k)
            expdif(i) = hmn(i,mx(i)) - hmin(i)
         end if
      end do
   end do
!
! *********************************************************
! compute approximate eps(z) using above taylor series
! *********************************************************
!
   do k = msg + 2,pver
      do i = 1,il2g
         expnum(i) = 0.
         ftemp(i) = 0.
         if (k < jt(i) .or. k >= jb(i)) then
            k1(i,k) = 0.
            expnum(i) = 0.
         else
            expnum(i) = hmn(i,mx(i)) - (hsat(i,k-1)*(zf(i,k)-z(i,k)) + &
                        hsat(i,k)* (z(i,k-1)-zf(i,k)))/(z(i,k-1)-z(i,k))
         end if
         if ((expdif(i) > 100. .and. expnum(i) > 0.) .and. &
         k1(i,k) > expnum(i)*dz(i,k)) then
            ftemp(i) = expnum(i)/k1(i,k)
            f(i,k) = ftemp(i) + i2(i,k)/k1(i,k)*ftemp(i)**2 + &
                     (2.*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2* &
                     ftemp(i)**3 + (-5.*k1(i,k)*i2(i,k)*i3(i,k)+ &
                     5.*i2(i,k)**3+k1(i,k)**2*i4(i,k))/ &
                     k1(i,k)**3*ftemp(i)**4
            f(i,k) = max(f(i,k),0._r8)
            f(i,k) = min(f(i,k),0.002_r8)
         end if
      end do
   end do
   do i = 1,il2g
      if (j0(i) < jb(i)) then
         if (f(i,j0(i)) < 1.E-6 .and. f(i,j0(i)+1) > f(i,j0(i))) j0(i) = j0(i) + 1
      end if
   end do
   do k = msg + 2,pver
      do i = 1,il2g
         if (k >= jt(i) .and. k <= j0(i)) then
            f(i,k) = max(f(i,k),f(i,k-1))
         end if
      end do
   end do
   do i = 1,il2g
      eps0(i) = f(i,j0(i))
      eps(i,jb(i)) = eps0(i)
   end do
!
! This is set to match the Rasch and Kristjansson paper
!
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k >= j0(i) .and. k <= jb(i)) then
            eps(i,k) = f(i,j0(i))
         end if
      end do
   end do
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k < j0(i) .and. k >= jt(i)) eps(i,k) = f(i,k)
      end do
   end do
!
! specify the updraft mass flux mu, entrainment eu, detrainment du
! and moist static energy hu.
! here and below mu, eu,du, md and ed are all normalized by mb
!
   do i = 1,il2g
      if (eps0(i) > 0.) then
!          mu(i,jb(i)) = 1.
!          eu(i,jb(i)) = eps0(i)/2.
!**pjr NOTE TO CORE GROUP mu and mu2 (and eu and eu2 should have the same vals now)
         mu2(i,jb(i)) = 1.
         eu2(i,jb(i)) = mu2(i,jb(i))/dz(i,jb(i))
         mu(i,jb(i)) = mu2(i,jb(i))
         eu(i,jb(i)) = eu2(i,jb(i))
      end if
   end do
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (eps0(i) > 0. .and. (k >= jt(i) .and. k < jb(i))) then
            zuef(i) = zf(i,k) - zf(i,jb(i))
            rmue(i) = (1./eps0(i))* (exp(eps(i,k+1)*zuef(i))-1.)/zuef(i)
            mu(i,k) = (1./eps0(i))* (exp(eps(i,k)*zuef(i))-1.)/zuef(i)
            eu(i,k) = (rmue(i)-mu(i,k+1))/dz(i,k)
            du(i,k) = (rmue(i)-mu(i,k))/dz(i,k)
            mu2(i,k) = mu(i,k)
            eu2(i,k) = eu(i,k)
            du2(i,k) = du(i,k)
         end if
      end do
   end do
!
   khighest = pverp
   klowest = 1
   do i=1,il2g
      khighest = min(khighest,lel(i))
      klowest = max(klowest,jb(i))
   end do
   do k = klowest-1,khighest,-1
!cdir$ ivdep
      do i = 1,il2g
         if (k <= jb(i)-1 .and. k >= lel(i) .and. eps0(i) > 0.) then
            if (mu(i,k) < 0.01) then
               hu(i,k) = hu(i,jb(i))
               mu(i,k) = 0.
               mu2(i,k) = mu(i,k)
               eu2(i,k) = 0.
               du2(i,k) = mu2(i,k+1)/dz(i,k)
               eu(i,k) = eu2(i,k)
               du(i,k) = du2(i,k)
            else
               hu(i,k) = mu(i,k+1)/mu(i,k)*hu(i,k+1) + &
                         dz(i,k)/mu(i,k)* (eu(i,k)*hmn(i,k)- du(i,k)*hsat(i,k))
            end if
         end if
      end do
   end do
!
! reset cloud top index beginning from two layers above the
! cloud base (i.e. if cloud is only one layer thick, top is not reset
!
   do i=1,il2g
      doit(i) = .true.
   end do
   do k=klowest-2,khighest-1,-1
      do i=1,il2g
         if (doit(i) .and. k <= jb(i)-2 .and. k >= lel(i)-1) then
       if (hu(i,k) <= hsthat(i,k) .and. hu(i,k+1) > hsthat(i,k+1) &
           .and. mu(i,k) >= 0.02) then
               if (hu(i,k)-hsthat(i,k) < -2000.) then
                  jt(i) = k + 1
                  doit(i) = .false.
               else
                  jt(i) = k
                  if (eps0(i) <= 0.) doit(i) = .false.
               end if
            else if (hu(i,k) > hu(i,jb(i)) .or. mu(i,k) < 0.02) then
               jt(i) = k + 1
               doit(i) = .false.
            end if
         end if
      end do
   end do
   do k = pver,msg + 1,-1
!cdir$ ivdep
      do i = 1,il2g
         if (k >= lel(i) .and. k <= jt(i) .and. eps0(i) > 0.) then
            mu(i,k) = 0.
            eu(i,k) = 0.
            du(i,k) = 0.
            mu2(i,k) = 0.
            eu2(i,k) = 0.
            du2(i,k) = 0.
            hu(i,k) = hu(i,jb(i))
         end if
         if (k == jt(i) .and. eps0(i) > 0.) then
            du(i,k) = mu(i,k+1)/dz(i,k)
            du2(i,k) = mu2(i,k+1)/dz(i,k)
            eu2(i,k) = 0.
            mu2(i,k) = 0.
            eu(i,k) = 0.
            mu(i,k) = 0.
         end if
      end do
   end do
!
! specify downdraft properties (no downdrafts if jd.ge.jb).
! scale down downward mass flux profile so that net flux
! (up-down) at cloud base in not negative.
!
   do i = 1,il2g
!
! in normal downdraft strength run alfa=0.2.  In test4 alfa=0.1
!
      alfa(i) = 0.1
      jt(i) = min(jt(i),jb(i)-1)
      jd(i) = max(j0(i),jt(i)+1)
      jd(i) = min(jd(i),jb(i))
      hd(i,jd(i)) = hmn(i,jd(i)-1)
      if (jd(i) < jb(i) .and. eps0(i) > 0.) then
         epsm(i) = eps0(i)
!          alfa(i)=2.*epsm(i)*( zf(i,jd(i))-zf(i,jb(i)) )/
!     1         (  exp(2.*epsm(i)*( zf(i,jd(i))-
!               zf(i,jb(i)) ))-1.  )
         md(i,jd(i)) = -alfa(i)*epsm(i)/eps0(i)
         md2(i,jd(i)) = md(i,jd(i))
      end if
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if ((k > jd(i) .and. k <= jb(i)) .and. eps0(i) > 0.) then
            zdef(i) = zf(i,jd(i)) - zf(i,k)
            md(i,k) = -alfa(i)/ (2.*eps0(i))*(exp(2.*epsm(i)*zdef(i))-1.)/zdef(i)
            md2(i,k) = md(i,k)
         end if
      end do
   end do
   do k = msg + 1,pver
!CDIR$ ivdep
      do i = 1,il2g
         if ((k >= jt(i) .and. k <= jb(i)) .and. eps0(i) > 0. .and. jd(i) < jb(i)) then
            ratmjb(i) = min(abs(mu2(i,jb(i))/md2(i,jb(i))),1._r8)
            md2(i,k) = md2(i,k)*ratmjb(i)
!            ratmjb(i) = min(abs(mu(i,jb(i))/md(i,jb(i))),1._r8)
!            md(i,k) = md(i,k)*ratmjb(i)
            md(i,k) = md2(i,k)
         end if
      end do
   end do

   small = 1.e-20
   do k = msg + 1,pver
      do i = 1,il2g
         if ((k >= jt(i) .and. k <= pver) .and. eps0(i) > 0.) then
            ed2(i,k-1) = (md2(i,k-1)-md2(i,k))/dz(i,k-1)
            ed(i,k-1) = ed2(i,k-1)
            mdt = min(md2(i,k),-small)
            hd(i,k) = (md(i,k-1)*hd(i,k-1) - dz(i,k-1)*ed(i,k-1)*hmn(i,k-1))/mdt
         end if
      end do
   end do
!
! calculate updraft and downdraft properties.
!
   do k = msg + 2,pver
      do i = 1,il2g
         if ((k >= jd(i) .and. k <= jb(i)) .and. eps0(i) > 0. .and. jd(i) < jb(i)) then
!         sd(i,k) = shat(i,k)
!    1             +              (hd(i,k)-hsthat(i,k))/
!    2               (cp    *(1.+gamhat(i,k)))
            qds(i,k) = qsthat(i,k) + gamhat(i,k)*(hd(i,k)-hsthat(i,k))/ &
               (rl*(1. + gamhat(i,k)))
         end if
      end do
   end do
!
   do i = 1,il2g
      done(i) = .false.
   end do
   kount = 0
   do k = pver,msg + 2,-1
      do i = 1,il2g
         if (( .not. done(i) .and. k > jt(i) .and. k < jb(i)) .and. eps0(i) > 0.) then
            su(i,k) = mu(i,k+1)/mu(i,k)*su(i,k+1) + &
                      dz(i,k)/mu(i,k)* (eu(i,k)-du(i,k))*s(i,k)
            qu(i,k) = mu(i,k+1)/mu(i,k)*qu(i,k+1) + dz(i,k)/mu(i,k)* (eu(i,k)*q(i,k)- &
                            du(i,k)*qst(i,k))
            tu = su(i,k) - grav/cp*zf(i,k)
            estu = c1*exp((c2* (tu-tfreez))/ ((tu-tfreez)+c3))
            qstu = eps1*estu/ ((p(i,k)+p(i,k-1))/2.-estu)
            if (qu(i,k) >= qstu) then
               jlcl(i) = k
               kount = kount + 1
               done(i) = .true.
            end if
         end if
      end do
      if (kount >= il2g) goto 690
   end do
690 continue
   do k = msg + 2,pver
      do i = 1,il2g
         if (k == jb(i) .and. eps0(i) > 0.) then
            qu(i,k) = q(i,mx(i))
            su(i,k) = (hu(i,k)-rl*qu(i,k))/cp
         end if
         if ((k > jt(i) .and. k <= jlcl(i)) .and. eps0(i) > 0.) then
            su(i,k) = shat(i,k) + (hu(i,k)-hsthat(i,k))/(cp* (1.+gamhat(i,k)))
            qu(i,k) = qsthat(i,k) + gamhat(i,k)*(hu(i,k)-hsthat(i,k))/ &
                     (rl* (1.+gamhat(i,k)))
         end if
      end do
   end do
!
   do k = pver,msg + 2,-1
      do i = 1,il2g
         if (k >= jt(i) .and. k < jb(i) .and. eps0(i) > 0.) then
            cu(i,k) = ((mu(i,k)*su(i,k)-mu(i,k+1)*su(i,k+1))/ &
                      dz(i,k)- (eu(i,k)-du(i,k))*s(i,k))/(rl/cp)
            if (k == jt(i)) cu(i,k) = 0.
!               cu(i,k) = max(0._r8,cu(i,k))
!               cu2     = max(0._r8,
!               cu2     = max(-1.e99_r8,
!     $                   +(eu(i,k)*q(i,k) - du(i,k)*qst(i,k))
!     $                   -(mu(i,k)*qu(i,k)-mu(i,k+1)*qu(i,k+1))/dz(i,k)
!     $                   )
!
!               if (abs(cu(i,k)-cu2)/(abs(cu(i,k))+abs(cu2)+1.e-50)
!     $              .gt.0.0000001) then
!                  write (6,*) ' inconsistent condensation rates ',
!     $                 i, k, lat(i),
!     $                 cu(i,k), cu2, jt(i), jb(i), jlcl(i), lel(i)
!     $                 ,mu(i,k)
!               endif
            cu(i,k) = max(0._r8,cu(i,k))
         end if
      end do
   end do
!
   beta = 0.
!   c0 = 2.E-3
   c0=3.E-4
   do k = pver,msg + 2,-1
      do i = 1,il2g
         cmfdqr(i,k) = 0.
! this modification is for test3 run, modified on 6/20/1995
!cc        if(t(i,jt(i) ).gt.tfreez)    c0=0.
!cc        if(t(i,jt(i) ).le.tfreez   )    c0=2.e-3
         if (k >= jt(i) .and. k < jb(i) .and. eps0(i) > 0. .and. mu(i,k) >= 0.0) then
            if (mu(i,k) > 0.) then
               ql1 = 1./mu(i,k)* (mu(i,k+1)*ql(i,k+1)- &
                     dz(i,k)*du(i,k)*ql(i,k+1)+dz(i,k)*cu(i,k))
               ql(i,k) = ql1/ (1.+dz(i,k)*c0)
            else
               ql(i,k) = 0.
            end if
            totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)* &
                        (beta*ql(i,k) + (1. - beta)*ql(i,k+1)))
            cmfdqr(i,k) = c0*mu(i,k)*ql(i,k)
         end if
      end do
   end do
!
   do i = 1,il2g
      qd(i,jd(i)) = qds(i,jd(i))
      sd(i,jd(i)) = (hd(i,jd(i)) - rl*qd(i,jd(i)))/cp
   end do
!
   do k = msg + 2,pver
      do i = 1,il2g
         if (k >= jd(i) .and. k < jb(i) .and. eps0(i) > 0.) then
            qd(i,k+1) = qds(i,k+1)
            evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k)-md(i,k+1)*qd(i,k+1))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._r8)
            mdt = min(md(i,k+1),-small)
            sd(i,k+1) = ((rl/cp*evp(i,k)-ed(i,k)*s(i,k))*dz(i,k) + md(i,k)*sd(i,k))/mdt
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
         end if
      end do
   end do
   do i = 1,il2g
!*guang         totevp(i) = totevp(i) + md(i,jd(i))*q(i,jd(i)-1) -
      totevp(i) = totevp(i) + md(i,jd(i))*qd(i,jd(i)) - md(i,jb(i))*qd(i,jb(i))
   end do
   if (.true.) then
      do i = 1,il2g
         k = jb(i)
         if (eps0(i) > 0.) then
            evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._r8)
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
         end if
      end do
   endif

   do i = 1,il2g
      totpcp(i) = max(totpcp(i),0._r8)
      totevp(i) = max(totevp(i),0._r8)
   end do
!
   weight = 1.0
   do k = msg + 2,pver
      do i = 1,il2g
         if (totevp(i) > 0. .and. totpcp(i) > 0.) then
            md2(i,k) = md2(i,k)*min(1._r8,weight*totpcp(i)/(totevp(i)+weight*totpcp(i)))
            ed2(i,k) = ed2(i,k)*min(1._r8,weight*totpcp(i)/(totevp(i)+weight*totpcp(i)))
            evp(i,k) = evp(i,k)*min(1._r8,weight*totpcp(i)/ (totevp(i)+weight*totpcp(i)))
         else
            md2(i,k) = 0.
            ed2(i,k) = 0.
            evp(i,k) = 0.
         end if
         md(i,k) = md2(i,k)
         ed(i,k) = ed2(i,k)
         cmeg(i,k) = cu(i,k) - evp(i,k)
!           cmeg is the cloud water condensed - rain water evaporated
!           cmfdqr  is the cloud water converted to rain - (rain evaporated)
         cmfdqr(i,k) = cmfdqr(i,k)-evp(i,k)
      end do
   end do
   do k = 2,pverp
      do i = 1,il2g
         pflx(i,k) = pflx(i,k-1) + cmfdqr(i,k-1)*dz(i,k-1)
      end do
   end do
   do i = 1,il2g
      if (totevp(i) > 0. .and. totpcp(i) > 0.) then
         totevp(i) = totevp(i)*min(1._r8,weight*totpcp(i)/(totevp(i) + weight*totpcp(i)))
      else
         totevp(i) = 0.
      end if
   end do
!
   do k = msg + 1,pver
      do i = 1,il2g
         mc(i,k) = mu(i,k) + md(i,k)
      end do
   end do
!
   return
end subroutine cldprp

end module zm_conv
