#include <misc.h>
#include <params.h>
#define PCWDETRAIN

!-----------------------------------------------------------------------
!
! Purpose:
! Tendency physics BEFORE coupling to land, sea, and ice models.
!
! Method:
! Call physics subroutines and compute the following:
!     o cloud calculations (cloud fraction, emissivity, etc.)
!     o radiation calculations
! Pass surface fields for separate surface flux calculations
! Dump appropriate fields to history file.
!
! Author: CCM1, CMS Contact: J. Truesdale
!
!-----------------------------------------------------------------------

subroutine tphysbc(ztodt,   pblht,   tpert,   ts,      &
                   qpert,                              &
                   precl,   precc,   precsl,  precsc,  &
                   asdir,   asdif,   aldir,   aldif,   &
                   snowh,                              &
                   qrs,     qrl,                       & ! radiation heating rate
                   flwds,                              &
                   fsns,    fsnt,    flns,    flnt,    & !
                   lwup,    srfrad,                    & ! surface radiations
                   sols,    soll,    solsd,   solld,   & ! solar radiations
                   cldo,    cldn,                      &
                   tcwato,  tcwatn,                    &
                   qcwato,  qcwatn,                    &
                   lcwato,  lcwatn,                    &
                   state,   tend,                      & ! state and tendencies
                   icefrac, landfrac, ocnfrac,         & ! fractions
                   tin,                                & ! DONG Li: clean this
                   cflx,    prcsnw,   state0 ,         & ! added by WAN Hui
                   sst, pbuf)                            ! added by SHI Xiangjun

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,     only: get_rlat_all_p, get_rlon_all_p
   use cldwat,        only: pcond
   use geopotential,  only: geopotential_t
   use physics_types, only: physics_state, physics_tend, physics_ptend, physics_update, physics_ptend_init
   use diagnostics,   only: diag_dynvar
   use history,       only: outfld
   use physconst,     only: gravit, latvap, cpair, tmelt, cappa, zvir, rair, rga
   use radheat,       only: radheat_net
   use constituents,  only: pcnst, pnats, ppcnst, qmin,cnst_get_ind
   use tracers,       only: dcconnam, ixcldw
   use zm_conv,       only: zm_conv_evap, zm_convr
   use zm_conv_3,     only: zm_convr_3,momtran_3
   use time_manager,  only: is_first_step, get_nstep, get_curr_calday
!! use moistconvection, only: cmfmca                           !!(wh)
   use moistconvection, only: cmfmca,limcnv,convection_scheme  !!(wh,2004.12.17)
   use m_cucall, only: cucall                                  !!(wh)

   use cloudsimulator,  only: doisccp, cloudsimulator_run      !!(wh 2005.01.28,following cam3.0)

   use MG,             only:stratiform_tend                 !! sxj
   use phys_buffer,    only: pbuf_fld,pbuf_size_max,pbuf_old_tim_idx, &
                           pbuf_get_fld_idx              !!

   USE pmgrid,                  ONLY: masterproc,iam     !sxj-2008-11-08
#ifdef SPMD                                              !sxj
   USE mpishorthand,      only:mpicom
#endif                                                   !sxj 2008-11-09

   implicit none

#include <comctl.h>
#include <RK_or_MG.h>

   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in) :: cflx(pcols)                    ! surface moisture !!(wh according to Ping Liu 2003)
   real(r8), intent(in) :: ts(pcols)                      ! surface temperature
   real(r8), intent(in) :: tcwato(pcols,pver)             !cloud water old temperature
   real(r8), intent(in) :: qcwato(pcols,pver)             ! cloud water old q
   real(r8), intent(in) :: lcwato(pcols,pver)             ! cloud liquid water old q
   real(r8), intent(in) :: icefrac(pcols)                 ! sea ice fraction (fraction)
   real(r8), intent(in) :: landfrac(pcols)                ! land fraction (fraction)
   real(r8), intent(in) :: ocnfrac(pcols)                 ! ocean fraction (fraction)
   real(r8), intent(inout) :: pblht(pcols)                ! Planetary boundary layer height
   real(r8), intent(inout) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(inout) :: qpert(pcols,ppcnst)         ! Thermal humidity & constituent excess
   real(r8), intent(in) :: asdir(pcols)                  ! Albedo: shortwave, direct
   real(r8), intent(in) :: asdif(pcols)                  ! Albedo: shortwave, diffuse
   real(r8), intent(in) :: aldir(pcols)                  ! Albedo: longwave, direct
   real(r8), intent(in) :: aldif(pcols)                  ! Albedo: longwave, diffuse
   real(r8), intent(in) :: snowh(pcols)                  ! Snow depth (liquid water equivalent)
   real(r8), intent(inout) :: qrs(pcols,pver)            ! Shortwave heating rate
   real(r8), intent(inout) :: qrl(pcols,pver)            ! Longwave  heating rate
   real(r8), intent(inout) :: flwds(pcols)               ! Surface longwave down flux
   real(r8), intent(inout) :: fsns(pcols)                   ! Surface solar absorbed flux
   real(r8), intent(inout) :: fsnt(pcols)                   ! Net column abs solar flux at model top
   real(r8), intent(inout) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
   real(r8), intent(inout) :: flnt(pcols)                   ! Net outgoing lw flux at model top
   real(r8), intent(in) :: lwup(pcols)                    ! Surface longwave up flux
   real(r8), intent(out) :: srfrad(pcols)                 ! Net surface radiative flux (watts/m**2)
   real(r8), intent(inout) :: sols(pcols)                   ! Direct beam solar rad. onto srf (sw)
   real(r8), intent(inout) :: soll(pcols)                   ! Direct beam solar rad. onto srf (lw)
   real(r8), intent(inout) :: solsd(pcols)                  ! Diffuse solar radiation onto srf (sw)
   real(r8), intent(inout) :: solld(pcols)                  ! Diffuse solar radiation onto srf (lw)
   real(r8), intent(out) :: precl(pcols)                  ! Large-scale precipitation rate
   real(r8), intent(out) :: precc(pcols)                  ! Convective-scale preciptn rate
   real(r8), intent(out) :: precsl(pcols)                 ! L.S. snowfall rate
   real(r8), intent(out) :: precsc(pcols)                 ! C.S. snowfall rate
   real(r8), intent(inout) :: cldo(pcols,pver)            !old cloud fraction
   real(r8), intent(out) :: cldn(pcols,pver)              !new cloud fraction
   real(r8), intent(out) :: tcwatn(pcols,pver)            !cloud water new temperature
   real(r8), intent(out) :: qcwatn(pcols,pver)            ! cloud water new q
   real(r8), intent(out) :: lcwatn(pcols,pver)            ! cloud liq. water new q
   real(r8), intent(out) :: tin(pcols,pver)               ! input T, to compute FV output T
   real(r8), intent(out) :: prcsnw(pcols)                 ! snowfall rate (precsl + precsc)
   real(r8), intent(in) :: sst(pcols)

!! type(physics_state), intent(inout) :: state
   type(physics_state), intent(inout) :: state,state0     !!(wh)
   type(physics_tend ), intent(inout) :: tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max)     :: pbuf  !!-sxj---
!
!---------------------------Local workspace-----------------------------
!
   real(r8) :: preclint(pcols)               ! Large-scale precipitation rate (>0.1mm/hr)
   real(r8) :: preccint(pcols)               ! Convective-scale preciptn rate (>0.1mm/hr)
   real(r8) :: preclfrq(pcols)               ! Large-scale precipitation frequency (>0.1mm/hr)
   real(r8) :: preccfrq(pcols)               ! Convective-scale preciptn frequency (>0.1mm/hr)
   real(r8) :: rhdfda(pcols,pver)            ! dRh/dcloud, old
   real(r8) :: rhu00 (pcols,pver)            ! Rh threshold for cloud, old

   type(physics_ptend)   :: ptend                  ! indivdual parameterization tendencies

   integer :: nstep                          ! current timestep number

   real(r8) :: calday                        ! current calendar day
   real(r8) :: clat(pcols)                   ! current latitudes(radians)
   real(r8) :: clon(pcols)                   ! current longitudes(radians)

   real(r8) :: zdu(pcols,pver)               ! detraining mass flux from deep convection
   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) :: cmfdqr(pcols,pver)            ! dq/dt due to moist convective rainout
   real(r8) :: cmfmc(pcols,pverp)            ! Convective mass flux--m sub c
   real(r8) :: cmfsl(pcols,pver)             ! Moist convection lw stat energy flux
   real(r8) :: cmflq(pcols,pver)             ! Moist convection total water flux
   real(r8) :: dtcond(pcols,pver)            ! dT/dt due to moist processes
   real(r8) :: dqcond(pcols,pver,ppcnst)     ! dq/dt due to moist processes

   real(r8) cldst(pcols,pver)
   real(r8) cltot(pcols)                      ! Diagnostic total cloud cover
   real(r8) cllow(pcols)                      !       "     low  cloud cover
   real(r8) clmed(pcols)                      !       "     mid  cloud cover
   real(r8) clhgh(pcols)                      !       "     hgh  cloud cover
   real(r8) cmfcme(pcols,pver)                ! cmf condensation - evaporation
   real(r8) cmfdqr2(pcols,pver)               ! dq/dt due to moist convective rainout
   real(r8) cmfmc2(pcols,pver)                ! Moist convection cloud mass flux
   real(r8) cmfsl2(pcols,pver)                ! Moist convection lw stat energy flux
   real(r8) cmflq2(pcols,pver)                ! Moist convection total water flux
   real(r8) cnt(pcols)                        ! Top level of convective activity
   real(r8) cnb(pcols)                        ! Lowest level of convective activity
   real(r8) cnt2(pcols)                       ! Top level of convective activity
   real(r8) cnb2(pcols)                       ! Bottom level of convective activity
   real(r8) concld(pcols,pver)
   real(r8) coszrs(pcols)                     ! Cosine solar zenith angle
   real(r8) dlf(pcols,pver)                   ! Detraining cld H20 from convection
   real(r8) fwaut(pcols,pver)
   real(r8) fsaut(pcols,pver)
   real(r8) fracw(pcols,pver)
   real(r8) fsacw(pcols,pver)
   real(r8) fsaci(pcols,pver)
   real(r8) nevapr(pcols,pver)                ! local evaporation of precipitation
   real(r8) prain(pcols,pver)                 ! local formation of precipitation
   real(r8) pflx(pcols,pverp)                 ! Conv rain flux thru out btm of lev
   real(r8) precc2(pcols)                     ! Convective-scale preciptn rate
   real(r8) preclp(pcols)                     ! sfc flux of precip from pcond
   real(r8) prect(pcols)                      ! total (conv+large scale) precip rate
   real(r8) qc(pcols,pver)                    ! dq/dt due to rainout terms
   real(r8) qc2(pcols,pver)                   ! dq/dt due to rainout terms
   real(r8) qme(pcols,pver)                   ! local condensation of cloud water
   real(r8) ptendlaq(pcols,pver)
   real(r8) ptendlat(pcols,pver)              !ljli20090830
   real(r8) qpert2(pcols,ppcnst)              ! Perturbation q
   ! DONG Li: Is "rtdt" needed?
   real(r8) rtdt                              ! 1./ztodt
   real(r8) tpert2(pcols)                     ! Perturbation T
   real(r8) tvm(pcols,pver)                   ! Virtual temperature
   real(r8) pmxrgn(pcols,pverp)               ! Maximum values of pressure for each
!                                             !    maximally overlapped region.
!                                             !    0->pmxrgn(i,1) is range of pressure for
!                                             !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!                                             !    2nd region, etc
   integer lchnk                              ! chunk identifier
   integer ncol                               ! number of atmospheric columns

   integer nmxrgn(pcols)                      ! Number of maximally overlapped regions
   integer  i,k,m                             ! Longitude, level, constituent indices

!  real(r8) engt                              ! Thermal   energy integral
!  real(r8) engk                              ! Kinetic   energy integral
!  real(r8) engp                              ! Potential energy integral
   real(r8) clwp(pcols,pver)                  ! Presribed cloud liq. h2o path
   real(r8) rel(pcols,pver)                   ! Liquid cloud particle effective radius
   real(r8) rei(pcols,pver)                   ! Ice effective drop size (microns)
   real(r8) fice(pcols,pver)                  ! Fractional ice content within cloud
   real(r8) fice_MG(pcols,pver)               !sxj        ! Fractional ice content within cloud
   real(r8) effcld(pcols,pver)                ! Effective cloud=cld*emis
   real(r8) emis(pcols,pver)                  ! Cloud longwave emissivity
   real(r8) clc(pcols)                        ! Total convective cloud (cloud scheme)
   real(r8) qtend(pcols,pver)                 ! moisture tendencies
   real(r8) ttend(pcols,pver)                 ! temp tendencies
   real(r8) lctend(pcols,pver)                ! cloud liquid water tendencies
   real(r8) rmelt(pcols,pver)                 ! heating rate due to phase change of precip
   real(r8) zero(pcols,pverp)                 ! a dummy array
   real(r8) clwp2(pcols,pver)                 ! in-cloud cloud water path
   real(r8) cliqwp(pcols,pver)                ! in-cloud cloud ice water path   !!(wh 2005.01.28)
   real(r8) cicewp(pcols,pver)                ! in-cloud cloud liquid water path!!(wh following cam3.0)
   real(r8) gclwp2(pcols,pver)                ! grid-box cloud water path
   real(r8) tgcwp(pcols)                      ! Vertically integrated cloud water path
   real(r8) tgiwp(pcols)                      ! Vertically integrated ice water path
   real(r8) tglwp(pcols)                      ! Vertically integrated liquid water path
   real(r8) tpw(pcols)                        ! Total precipitable water
   real(r8) hl (pcols)                        ! Liquid water scale height
!
   real(r8) dellow(pcols)                     ! delta p for bottom three levels of model
   real(r8) tavg(pcols)                       ! mass weighted average temperature for
!
! Used for OUTFLD only
!
   real(r8) icimr(pcols,pver)                 ! in cloud ice mixing ratio
   real(r8) icwmr(pcols,pver)                 ! in cloud water mixing ratio
   real(r8) cldv(pcols,pver)                  ! cloud volume (fraction) occupied by rain or cloud water
   ! DONG Li: "rain" is not being used
   real(r8) rain(pcols,pver)                  ! total precip mixing ratio
   real(r8) icwmr1(pcols,pver)                ! in cloud water mixing ration for zhang scheme
   real(r8) icwmr2(pcols,pver)                ! in cloud water mixing ration for hack scheme
   real(r8) fracis(pcols,pver,ppcnst)         ! fraction of transported species that are insoluble
   real(r8) evappct(pcols)                    ! Convective-scale preciptn rate
   real(r8) timestep(pcols)
!
!     Variables for doing deep convective transport outside of zm_convr
!
   real(r8) mu2(pcols,pver)
   real(r8) eu2(pcols,pver)
   real(r8) du2(pcols,pver)
   real(r8) md2(pcols,pver)
   real(r8) ed2(pcols,pver)
   real(r8) dp(pcols,pver)
   real(r8) dsubcld(pcols)
!  for Zhang _ Neale
   real(r8) cape(pcols)   ! w  convective available potential energy.                   !sxj^M
   real(r8) rliq(pcols)   ! reserved liquid (not yet in cldliq) for energy integrals    !sxj^M
   real(r8) winds(pcols, pver, 2)
   real(r8) wind_tends(pcols, pver, 2)
   real(r8) pguall(pcols, pver, 2)
   real(r8) pgdall(pcols, pver, 2)
   real(r8) icwu(pcols,pver, 2)
   real(r8) icwd(pcols,pver, 2)
   real(r8) seten(pcols, pver)    
   logical l_windt(2)

   integer jt(pcols)
   integer maxg(pcols)
   integer ideep(pcols)
   integer lengath

! For Tiedtke scheme              !!
   integer itype(pcols)           !!
   real(r8) xtec(pcols,pver)      !!
   logical loland(pcols)          !!(wh 2004.12.17,according to Ping Liu 2003)
   real(r8) topmaxm(pcols)        !!
   integer index(pcols)           !!

! stratiform precipitation variables   !! sxj---
   real(r8) :: prec_str(pcols)     ! sfc flux of precip from stratiform (m/s)
   real(r8) :: snow_str(pcols)     ! sfc flux of snow from stratiform   (m/s)
   real(r8) :: prec_pcw(pcols)     ! total precip from prognostic cloud scheme
   real(r8) :: snow_pcw(pcols)     ! snow from prognostic cloud scheme
   real(r8) :: prec_sed(pcols)     ! total precip from cloud sedimentation
   real(r8) :: snow_sed(pcols)     ! snow from cloud ice sedimentation
   integer ifld,itim
   integer ix_cldice, ix_cldliq         !sxj
   logical aer_indirect                 !sxj
   logical detrain_to_cld               !sxj

!##############################--------------------------
    detrain_to_cld=.true.
!##############################-------------------------
!
!-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   nstep = get_nstep()
   calday = get_curr_calday()
!
! Output NSTEP for debugging
!
   timestep(:ncol) = nstep
   call outfld ('NSTEP   ',timestep, pcols, lchnk)

!--------------------------------------------------------------------
        !CALL mpibarrier (mpicom)  !debug
        !if (masterproc) then
   !write(6,*) "pbuf test--tphysbc.F90"
   !write(6,*) "pbuf(2)%name:",pbuf(2)%name
   !write(6,*) "pbuf(2)%scope:",pbuf(2)%scope
   !write(6,*) "pbuf(2)%fdim:",pbuf(2)%fdim
   !write(6,*) "pbuf(2)%mdim:",pbuf(2)%mdim
   !write(6,*) "pbuf(2)%ldim:",pbuf(2)%ldim
   !write(6,*) pbuf(12)%fld_ptr(1,1,1,1,1)
!endif
!CALL mpibarrier (mpicom)
!call endrun               !debug
!---------------------------------------------
!
!*** BAB's FV kludge
!
   tin(:ncol,:pver) = state%t(:ncol,:pver)
!
! Convert mixing ratio of non-water tracers to mass fraction of total
! atmospheric mass (Overwrite non-water portions of q3m1).
!
   if (ppcnst > 1) then
      call mr2mf (lchnk, ncol, state%q)
   end if
!
! Set reciprocal of layer thickness
! Set physics tendencies to 0
!
   state%rpdel(:ncol,:pver) = 1./state%pdel(:ncol,:pver)
   tend %dTdt(:ncol,:pver)  = 0.
   tend %dudt(:ncol,:pver)  = 0.
   tend %dvdt(:ncol,:pver)  = 0.
!
! Compute initial geopotential heights and dry static energies
!
   call geopotential_t (state%lnpint, state%lnpmid  , state%pint,  &
                        state%pmid  , state%pdel    , state%rpdel, &
                        state%t     , state%q(1,1,1), rair , gravit , zvir   , &
                        state%zi    , state%zm      , ncol      )

   state%s(:ncol,:pver) = cpair*state%t(:ncol,:pver) + gravit*state%zm(:ncol,:pver)
   call physics_ptend_init (ptend) ! Initialize parameterization tendency structure
!
! Make sure that input tracers are all positive (probably unnecessary)
!
   call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
              ppcnst,qmin  ,state%q )
!
! Setup q and t accumulation fields
!
   dqcond(:ncol,:,:) = state%q(:ncol,:,:)
   dtcond(:ncol,:)   = state%s(:ncol,:)
!
! Zero out precip and convective fields before accumulating terms
!
   precl (:ncol)   = 0.
   preclp(:ncol)   = 0.
   precc (:ncol)   = 0.
   precsl(:ncol)   = 0.
   precsc(:ncol)   = 0.
   qc    (:ncol,:) = 0.
   cmfdqr(:ncol,:) = 0.
   cmfmc (:ncol,:) = 0.
   cmfsl (:ncol,:) = 0.
   cmflq (:ncol,:) = 0.
   dqcond(ncol+1:pcols,:,:) = 0.
   dtcond(ncol+1:pcols,:)   = 0.

   fracis (:ncol,:,1:ppcnst) = 1.
!
!===================================================
! Dry adjustment
!===================================================

! Copy state info for input to dadadj
! This is a kludge, so that dadadj does not have to be correctly reformulated in dry static energy

   ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
   ptend%q(:ncol,:pver,1) = state%q(:ncol,:pver,1)

   call t_startf ('dadadj')
   call dadadj (lchnk, ncol, state%pmid,  state%pint,  state%pdel,  &
                ptend%s, ptend%q(1,1,1))
   ptend%name  = 'dadadj'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.
   ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/ztodt * cpair
   ptend%q(:ncol,:,1) = (ptend%q(:ncol,:,1) - state%q(:ncol,:,1))/ztodt
   call t_stopf ('dadadj')
   call physics_update (state, tend, ptend, ztodt)
!
!===================================================
! Moist convection
!===================================================
!
! Since the PBL doesn't pass constituent perturbations, they
! are zeroed here for input to the moist convection routine
!
   qpert(:ncol,2:ppcnst) = 0.0
!
! JR Set some arrays to zero @ nstep=0. Otherwise random junk off the heap
! or stack will be used in zm_convr
!
   if (is_first_step()) then
      pblht(:ncol)  = 0.
      tpert(:ncol)  = 0.
      cldo (:ncol,:)= 0.
   end if

   if (convection_scheme == 'Zhang-Hack') then
!      write(6,*) 'Zhang-Hack moistconvection' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Begin with Zhang-McFarlane (1996) convection parameterization
!
!gz070721
   ptendlaq(:ncol,:pver)=(state%q(:ncol,:pver,1)-state%q1(:ncol,:pver))/ztodt
   ptendlat(:ncol,:pver)=(state%t(:ncol,:pver)-state%t1(:ncol,:pver))/ztodt
!endgz070721
   call t_startf ('zm_convr')
   call zm_convr( lchnk,    ncol, &
                  state%t,   state%q,    precc,   cnt,     cnb,      &
                  pblht,   state%zm, state%phis,    state%zi,   ptend%q(:,:,1),     &
                  ptend%s, state%pmid,   state%pint,  state%pdel,  ts,       &
                  .5*ztodt,cmfmc,    cmfcme,  nstep,             &
                  tpert,   dlf,      pflx,    zdu,     cmfdqr,   &
                  mu2,      md2,     du2,     eu2,     ed2,      &
                  dp,       dsubcld, jt,      maxg,    ideep,    &
                  lengath, icwmr1,ptendlaq,ptendlat )
   ptend%name  = 'zm_convr'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.

   ftem(:ncol,:pver) = ptend%s(:ncol,:pver)/cpair
   call outfld('ZMDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('ZMDQ    ',ptend%q(1,1,1) ,pcols   ,lchnk   )
   call t_stopf('zm_convr')

   call physics_update(state, tend, ptend, ztodt)
!
! Evaporate some of the precip directly into the environment (Sundqvist)
!
   call zm_conv_evap(state, ptend, pflx, precc, cldo, ztodt, evappct)
   ptend%name  = 'zm_conv_evap'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.
   call outfld('EVAPPCT ',evappct,pcols,state%lchnk)
   call physics_update(state, tend, ptend, ztodt)
!ljli20090830
   state%t1(:ncol,:)=state%t(:ncol,:)
   state%q1(:ncol,:)=state%q(:ncol,:,1)
!ljli20090830
!
! Transport cloud water only
!
   ptend%name = 'convtran1'
   do m=2,ppcnst
      if (m == ixcldw) ptend%lq(m) = .true.
   end do
   call t_startf ('convtran1')
   call convtran (lchnk,                                        &
                  ptend%lq,state%q, ppcnst,  mu2,     md2,   &
                  du2,     eu2,     ed2,     dp,      dsubcld,  &
                  jt,      maxg,    ideep,   1,       lengath,  &
                  nstep,   fracis,  ptend%q   )
   call t_stopf ('convtran1')
   call physics_update (state, tend, ptend, ztodt)
!
! Convert mass flux from reported mb/s to kg/m^2/s
!
   cmfmc(:ncol,:pver) = cmfmc(:ncol,:pver) * 100./gravit
!
! Add production of rain by zm_convr to qc.  Added 1 to k-index of pflx
! at instruction of PJR
!
   do k=2,pver
      do i=1,ncol
         qc(i,k) = qc(i,k) + (pflx(i,k+1) - pflx(i,k))*gravit/state%pdel(i,k)
      end do
   end do
!
! Call Hack (1994) convection scheme to deal with shallow/mid-level convection
! Begin by zeroing local copies of mass flux, energy fluxes, etc.
!
   cmfmc2 (:ncol,:pver) = 0.
   cmfdqr2(:ncol,:pver) = 0.
   cmfsl2 (:ncol,:pver) = 0.
   cmflq2 (:ncol,:pver) = 0.
   qc2    (:ncol,:pver) = 0.
!
! At PJR's instruction, deleted kludge to get past a once in a lifetime
! problem in cmfmca's transport of liq water due to reliance on m=2 being
! hard-wired to cloud water--JR.  Put back in after run bombed.
!
   where (abs(state%q(:ncol,:pver,ixcldw)) < 1.e-36)
      state%q(:ncol,:pver,ixcldw) = 0.
   end where

   call t_startf('cmfmca')
   tpert2(:ncol  ) =0.
   qpert2(:ncol,:) = qpert(:ncol,:)  ! BAB Why is this not zero, if tpert2=0???
   call cmfmca (lchnk,   ncol, &
                nstep,   ztodt,   state%pmid,  state%pdel,   &
                state%rpdel,   state%zm,      tpert2,  qpert2,  state%phis,     &
                pblht,   state%t,   state%q,   ptend%s,   ptend%q,      &
                cmfmc2,  cmfdqr2, cmfsl2,  cmflq2,  precc2,   &
                qc2,     cnt2,    cnb2,    icwmr2   )
   ptend%name  = 'cmfmca'
   ptend%ls    = .TRUE.
   ptend%lq(:) = .TRUE.

   ftem(:ncol,:pver) = ptend%s(:ncol,:pver)/cpair
   call outfld('CMFDT   ',ftem          ,pcols   ,lchnk   )
   call outfld('CMFDQ   ',ptend%q(1,1,1),pcols   ,lchnk   )
   call t_stopf('cmfmca')
   call physics_update (state, tend, ptend, ztodt)
!
! Merge shallow/mid-level output with prior results from Zhang-McFarlane
!
   do i=1,ncol
      precc(i) = precc(i) + precc2(i)
      if (cnt2(i) < cnt(i)) cnt(i) = cnt2(i)
      if (cnb2(i) > cnb(i)) cnb(i) = cnb2(i)
   end do
!
   cmfmc(:ncol,:pver)  = cmfmc(:ncol,:pver)  + cmfmc2(:ncol,:pver)
   cmfdqr(:ncol,:pver) = cmfdqr(:ncol,:pver) + cmfdqr2(:ncol,:pver)
   cmfsl(:ncol,:pver)  = cmfsl(:ncol,:pver)  + cmfsl2(:ncol,:pver)
   cmflq(:ncol,:pver)  = cmflq(:ncol,:pver)  + cmflq2(:ncol,:pver)
   qc(:ncol,:pver)     = qc(:ncol,:pver)     + qc2(:ncol,:pver)


   if (RK_or_MG=='RK') then  !!sxj--MG scheme include pcwdetrain-2008-11-20
   if (detrain_to_cld) then !sxj---
     ! write(*,*) "detrain_to_cld"
     do k = 1,pver
      do i = 1,ncol
         ptend%q(i,k,ixcldw) = dlf(i,k)
      end do
     end do
     ptend%name  = 'pcwdetrain'
     ptend%lq(ixcldw) = .TRUE.
     call physics_update(state, tend, ptend, ztodt)
   else
#ifndef PCWDETRAIN
!
! put the detraining cloud water into precip to conserve
! mass
!
   do k = 1,pver
      do i = 1,ncol
         precc(i) = precc(i) + dlf(i,k)*state%pdel(i,k)/(gravit*1000.)
      end do
   end do
#else
!
! put the detraining cloud water into the cloud and environment in
! proportion to the cloud fraction
!
   do k = 1,pver
      do i = 1,ncol
         ptend%q(i,k,1)      = dlf(i,k)*(1.-cldo(i,k))
         ptend%s(i,k)        =-dlf(i,k)*(1.-cldo(i,k))*latvap
         ptend%q(i,k,ixcldw) = dlf(i,k)*cldo(i,k)
      end do
   end do
   ptend%name  = 'pcwdetrain'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.
   ptend%lq(ixcldw) = .TRUE.
   call physics_update(state, tend, ptend, ztodt)
#endif
   endif ! detrain_to_cld
   endif ! RK_or_MG=='RK'

   end if    !( convection_scheme = 'Zhang-Hack' ) !!(wh 2004.12)

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
    if (convection_scheme == 'Zhang_Neale') then   !20101116ljliaccording to SXJ
    call t_startf ('zm_convr')
       call zm_convr_3(   lchnk   ,ncol    , &
                          state%t       ,state%q     ,precc    ,cnt   ,cnb   , &
                          pblht    ,state%zm      ,state%phis    ,state%zi      ,ptend%q(:,:,1)    , &
                          ptend%s    ,state%pmid     ,state%pint    ,state%pdel     , &
                          .5_r8*ztodt    ,cmfmc  ,cmfcme    , cape,      &
                          tpert   ,dlf     ,pflx    ,zdu     ,cmfdqr   , &
                          mu2   ,md2    ,du2    ,eu2      ,ed2        , &
                          dp      ,dsubcld     ,jt   ,maxg   ,ideep   , &
                          lengath ,icwmr1 ,rliq    )
       call outfld('CAPE', cape, pcols, lchnk)
       ptend%name  = 'zm_convr'
       ptend%ls    = .TRUE.
       ptend%lq(1) = .TRUE.
       ftem(:ncol,:pver) = ptend%s(:ncol,:pver)/cpair
       call outfld('ZMDT    ',ftem           ,pcols   ,lchnk   )
       call outfld('ZMDQ    ',ptend%q(1,1,1) ,pcols   ,lchnk   )
       call t_stopf('zm_convr')
       call physics_update(state, tend, ptend, ztodt)
! Evaporate some of the precip directly into the environment (Sundqvist)^
       call zm_conv_evap(state, ptend, pflx, precc, cldo, ztodt, evappct)
       ptend%name  = 'zm_conv_evap'
       ptend%ls    = .TRUE.
       ptend%lq(1) = .TRUE.
       call outfld('EVAPPCT ',evappct,pcols,state%lchnk)
       call physics_update(state, tend, ptend, ztodt)
       winds(:ncol,:pver,1) = state%u(:ncol,:pver)
       winds(:ncol,:pver,2) = state%v(:ncol,:pver)  
       l_windt(1) = .true.
       l_windt(2) = .true.
       call t_startf ('momtran')
       call momtran_3(lchnk, ncol,                                        &
                      l_windt,winds, 2,  mu2, md2,   &
                      du2, eu2, ed2, dp, dsubcld,  &
                      jt ,maxg , ideep , 1, lengath,  &
                      nstep,  wind_tends, pguall, pgdall, icwu, icwd, ztodt, seten )  
       call t_stopf ('momtran')
       ptend%name  = 'momtran'
       ptend%lu = .TRUE.
       ptend%lv = .TRUE.
       ptend%ls = .TRUE.
       ptend%u(:ncol,:pver) = wind_tends(:ncol,:pver,1)
       ptend%v(:ncol,:pver) = wind_tends(:ncol,:pver,2)
       ptend%s(:ncol,:pver) = seten(:ncol,:pver)
       call physics_update(state, tend, ptend, ztodt)
       call outfld('ZMMTU', wind_tends(1,1,1), pcols, lchnk)
       call outfld('ZMMTV', wind_tends(1,1,2), pcols, lchnk)
! Output apparent force from  pressure gradient
       call outfld('ZMUPGU', pguall(1,1,1), pcols, lchnk)
       call outfld('ZMUPGD', pgdall(1,1,1), pcols, lchnk)
       call outfld('ZMVPGU', pguall(1,1,2), pcols, lchnk)
       call outfld('ZMVPGD', pgdall(1,1,2), pcols, lchnk)
! Output in-cloud winds^M
       call outfld('ZMICUU', icwu(1,1,1), pcols, lchnk)
       call outfld('ZMICUD', icwd(1,1,1), pcols, lchnk)
       call outfld('ZMICVU', icwu(1,1,2), pcols, lchnk)
       call outfld('ZMICVD', icwd(1,1,2), pcols, lchnk)
! Transport cloud water only
       ptend%name = 'convtran1'
       do m=2,ppcnst
           if (m == ixcldw) ptend%lq(m) = .true.
       end do
       call t_startf ('convtran1')
       call convtran (lchnk,                                        &
                      ptend%lq,state%q, ppcnst,  mu2,     md2,   &
                      du2,     eu2,     ed2,     dp,      dsubcld,  &
                      jt,      maxg,    ideep,   1,       lengath,  &
                      nstep,   fracis,  ptend%q   )
       call t_stopf ('convtran1')
       call physics_update (state, tend, ptend, ztodt)
! Convert mass flux from reported mb/s to kg/m^2/s
        cmfmc(:ncol,:pver) = cmfmc(:ncol,:pver) * 100./gravit
! Add production of rain by zm_convr to qc.  Added 1 to k-index of pflx
! at instruction of PJR
!
        do k=2,pver
          do i=1,ncol
               qc(i,k) = qc(i,k) + (pflx(i,k+1) - pflx(i,k))*gravit/state%pdel(i,k)
          end do
        end do
!
! Call Hack (1994) convection scheme to deal with shallow/mid-level convection
! Begin by zeroing local copies of mass flux, energy fluxes, etc.
!
        cmfmc2 (:ncol,:pver) = 0.
        cmfdqr2(:ncol,:pver) = 0.
        cmfsl2 (:ncol,:pver) = 0.
        cmflq2 (:ncol,:pver) = 0.
        qc2    (:ncol,:pver) = 0.
!
! At PJR's instruction, deleted kludge to get past a once in a lifetime
! problem in cmfmca's transport of liq water due to reliance on m=2 being
! hard-wired to cloud water--JR.  Put back in after run bombed.
!
       where (abs(state%q(:ncol,:pver,ixcldw)) < 1.e-36)
         state%q(:ncol,:pver,ixcldw) = 0.
       end where
    
       call t_startf('cmfmca')
       tpert2(:ncol  ) =0.
       qpert2(:ncol,:) = qpert(:ncol,:)  ! BAB Why is this not zero, if tpert2=0???
       call cmfmca (lchnk,   ncol, &
                    nstep,   ztodt,   state%pmid,  state%pdel,   &
                    state%rpdel,   state%zm,      tpert2,  qpert2,  state%phis,     &
                    pblht,   state%t,   state%q,   ptend%s,   ptend%q,      &
                    cmfmc2,  cmfdqr2, cmfsl2,  cmflq2,  precc2,   &
                    qc2,     cnt2,    cnb2,    icwmr2   )
       ptend%name  = 'cmfmca'
       ptend%ls    = .TRUE.
       ptend%lq(:) = .TRUE.
     

      ftem(:ncol,:pver) = ptend%s(:ncol,:pver)/cpair
      call outfld('CMFDT   ',ftem          ,pcols   ,lchnk   )
      call outfld('CMFDQ   ',ptend%q(1,1,1),pcols   ,lchnk   )
      call t_stopf('cmfmca')
      call physics_update (state, tend, ptend, ztodt)
!
! Merge shallow/mid-level output with prior results from Zhang-McFarlane
!
      do i=1,ncol
         precc(i) = precc(i) + precc2(i)
         if (cnt2(i) < cnt(i)) cnt(i) = cnt2(i)
         if (cnb2(i) > cnb(i)) cnb(i) = cnb2(i)
      end do
!
      cmfmc(:ncol,:pver)  = cmfmc(:ncol,:pver)  + cmfmc2(:ncol,:pver)
      cmfdqr(:ncol,:pver) = cmfdqr(:ncol,:pver) + cmfdqr2(:ncol,:pver)
      cmfsl(:ncol,:pver)  = cmfsl(:ncol,:pver)  + cmfsl2(:ncol,:pver)
      cmflq(:ncol,:pver)  = cmflq(:ncol,:pver)  + cmflq2(:ncol,:pver)
      qc(:ncol,:pver)     = qc(:ncol,:pver)     + qc2(:ncol,:pver)
   
   
      if (RK_or_MG=='RK') then  !!sxj--MG scheme include pcwdetrain-2008-11-20
      if (detrain_to_cld) then !sxj-
! write(*,*) "detrain_to_cld"
        do k = 1,pver
         do i = 1,ncol
            ptend%q(i,k,ixcldw) = dlf(i,k)
         end do
        end do
             ptend%name  = 'pcwdetrain'
             ptend%lq(ixcldw) = .TRUE.
        call physics_update(state, tend, ptend, ztodt)
        else
#ifndef PCWDETRAIN
!
! put the detraining cloud water into precip to conserve
! mass
!
                  do k = 1,pver
                      do i = 1,ncol
                         precc(i) = precc(i) + dlf(i,k)*state%pdel(i,k)/(gravit*1000.)
                      end do
                  end do
#else
!
! put the detraining cloud water into the cloud and environment in
! proportion to the cloud fraction
!
              do k = 1,pver
                   do i = 1,ncol
                     ptend%q(i,k,1)      = dlf(i,k)*(1.-cldo(i,k))
                     ptend%s(i,k)        =-dlf(i,k)*(1.-cldo(i,k))*latvap
                     ptend%q(i,k,ixcldw) = dlf(i,k)*cldo(i,k)
                   end do
               end do
               ptend%name  = 'pcwdetrain'
               ptend%ls    = .TRUE.
               ptend%lq(1) = .TRUE.
               ptend%lq(ixcldw) = .TRUE.
               call physics_update(state, tend, ptend, ztodt)
#endif
               endif ! detrain_to_cld
               endif ! RK_or_MG=='RK'
           
          end if    ! end (Zhang_Neale)
!------------------------------------------------------------------------------------------
!---------------------------------------END------------------------------------------------
   if ( convection_scheme == 'Tiedtke' ) then     !!(wh 2004.12)
!      write(6,*) 'Tiedtke moistconvection'
!
! Begin the Tiedtke scheme from ECHAM4.61, Ping Liu, Apr. 2003
!
   state0%rpdel(:ncol,:pver) = 1./state0%pdel(:ncol,:pver)
!
! Compute initial geopotential heights
!
   call geopotential_t (state0%lnpint, state0%lnpmid  , state0%pint,  &
                        state0%pmid  , state0%pdel    , state0%rpdel, &
                        state0%t     , state0%q(1,1,1), rair, gravit, &
                        zvir, state0%zi    , state0%zm      , ncol    )

   itype(:ncol)=0 !integer itype(ncol)
   loland(:ncol)= nint(landfrac(:ncol)) == 1 !logical loland(pcols)
   topmaxm(:ncol)=99999. !top convective cloud max

   xtec(:ncol,:pver)=0.
   precc(:ncol)=0.
   mu2(:ncol,:pver)=0.
   md2(:ncol,:pver)=0.
   du2(:ncol,:pver)=0.
   eu2(:ncol,:pver)=0.
   ed2(:ncol,:pver)=0.
   zdu(:ncol,:pver)=0.

   ptend%q(:ncol,:pver,1)=(state%q(:ncol,:pver,1)-state0%q(:ncol,:pver,1))/ztodt
   do k=1,pver
     do i=1,ncol
      state0%zm(i,k)=state0%zm(i,k)*gravit+state0%phis(i)
     enddo
   enddo


   call t_startf ('cucall')
   CALL cucall(ncol,pver,pver+1,pver-1, &
     state%t,state0%q(1,1,1),state%u,state%v,state%q(1,1,ixcldw),& ! IN
             !here only q is from stat0, which can eliminate negative q
     ptend%s,ptend%q(1,1,1),ptend%u,ptend%v,& ! OUT
     ptend%q(1,1,ixcldw),    &! cloud water tendency  OUT
     state%omega, & !IN
     cflx,             &!moisture flux at the surface, IN
     xtec,             &!temp array (nglpx,nlev), INOUT, actu. for cond.
                        !tendency of detrained cloud water, Kg/Kg/s
     state%pmid,       &!full-level pressures, Pa IN
     state%pint,       &!half-level pressures, Pa IN
     state0%zm,        &!geopotential height, m^2/s^2 IN
     precc,            &!accum. conv. precipitation, kg/m^2/s INOUT
     itype,            &!temp array (nglpx), INOUT, convection type
     loland,           &!land logical, IN
     cnt,              &!Maximum convective cloud tops level number, INOUT
     topmaxm,          &!99999.0, IN
     cnb,              &!cloud base level
     cmfmc,            &!net mass flux, mu2+md2
     zdu,ztodt,        &
     mu2,              & !up mass flux   -----------|
     md2,              & !down mass flux            |
     du2,              & !detrainment in updraft    |--NT: flowing conversions
     eu2,              & !entrainment in updraft    |
     ed2               & !entrainment in downdraft -|
        )
   call t_stopf('cucall')


   precc(:ncol)=precc(:ncol)/ztodt
   dp(:ncol,:pver)=state0%pdel(:ncol,:pver)*0.01
   mu2(:ncol,:pver)=mu2(:ncol,:pver)* gravit/100.
   md2(:ncol,:pver)=md2(:ncol,:pver)* gravit/100.
   du2(:ncol,:pver)=du2(:ncol,:pver)* gravit/100.
   eu2(:ncol,:pver)=eu2(:ncol,:pver)* gravit/100.
   ed2(:ncol,:pver)=ed2(:ncol,:pver)* gravit/100.
   zdu(:ncol,:pver)=du2(:ncol,:pver)/dp(:ncol,:pver)

   jt(:ncol)=pver
   maxg(:ncol)=1
   dsubcld(:ncol)=0.

   lengath = 0
   do i=1,ncol
     if(itype(i)==1)then
        lengath=lengath+1
        index(lengath)=i

        jt(lengath)=nint(cnt(i))
        maxg(lengath)=nint(cnb(i))

        do k=1,pver
        dp(lengath,k)=dp(i,k)
        mu2(lengath,k)=mu2(i,k)
        md2(lengath,k)=md2(i,k)
        du2(lengath,k)=du2(i,k)/dp(i,k)
        eu2(lengath,k)=eu2(i,k)/dp(i,k)
        ed2(lengath,k)=ed2(i,k)/dp(i,k)


        if(i .ne. lengath)then
        mu2(i,k)=0.
        md2(i,k)=0.
        du2(i,k)=0.
        eu2(i,k)=0.
        ed2(i,k)=0.
        endif
        enddo

     else
        do k=1,pver
        mu2(i,k)=0.
        md2(i,k)=0.
        du2(i,k)=0.
        eu2(i,k)=0.
        ed2(i,k)=0.
        enddo
     endif
   enddo
   do i=1,lengath
     k=index(i)
     ideep(i)=k
   enddo
   do i=1,lengath
     do k=limcnv,pver
        if(k>=maxg(i))then
          dsubcld(i)=dsubcld(i)+dp(i,k)
        endif
     enddo
   enddo
!
! End of Tiedtke scheme
!
   ptend%name  = 'cucall'
   ptend%ls    = .TRUE.
   ptend%lu    = .TRUE.
   ptend%lv    = .TRUE.
   ptend%lq(1) = .TRUE.

   ptend%s(:ncol,:pver)=ptend%s(:ncol,:pver)*cpair
   ptend%q(:ncol,:pver,1)=ptend%q(:ncol,:pver,1)- &
            (state%q(:ncol,:pver,1)-state0%q(:ncol,:pver,1))/ztodt
   ptend%q(:ncol,:pver,ixcldw)=0.
   ftem(:ncol,:pver) = ptend%s(:ncol,:pver)/cpair
   call outfld('ZMDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('ZMDQ    ',ptend%q(1,1,1) ,pcols   ,lchnk   )

   call physics_update(state, tend, ptend, ztodt)


   if (RK_or_MG=='RK') then  !!sxj-- MG scheme include pcwdetrain-2008-11-20
! MG scheme include pcwdetrain  all detrain cloud water conver to stratiform cloud water/ice
   if (detrain_to_cld) then
     do k = 1,pver
      do i = 1,ncol
         ptend%q(i,k,ixcldw) = xtec(i,k)*cldo(i,k)
      end do
     end do
     ptend%name  = 'pcwdetrain'
     ptend%lq(ixcldw) = .TRUE.
     call physics_update(state, tend, ptend, ztodt)
   else
#ifndef PCWDETRAIN
!
! put the detraining cloud water into precip to conserve
! mass
!
   do k = 1,pver
      do i = 1,ncol
         precc(i) = precc(i) + xtec(i,k)*state%pdel(i,k)/(gravit*1000.)
      end do
   end do
#else
!
! put the detraining cloud water into the cloud and environment in
! proportion to the cloud fraction
!
   do k = 1,pver
      do i = 1,ncol
         ptend%q(i,k,1)      = xtec(i,k)*(1.-cldo(i,k))         !!(wh)
         ptend%s(i,k)        =-xtec(i,k)*(1.-cldo(i,k))*latvap  !!(wh)
         ptend%q(i,k,ixcldw) = xtec(i,k)*cldo(i,k)              !!(wh)
      end do
   end do
   ptend%name  = 'pcwdetrain'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.
   ptend%lq(ixcldw) = .TRUE.
   call physics_update(state, tend, ptend, ztodt)
#endif
   endif ! detrain_to_cld
   endif ! RK_or_MG=='RK'

!
! Transport cloud water only
!
   ptend%name = 'convtran1'
   do m=2,ppcnst
      if (m == ixcldw) ptend%lq(m) = .true.
   end do
   call t_startf ('convtran1')
   call convtran (lchnk,                                        &
                  ptend%lq,state%q, ppcnst,  mu2,     md2,   &
                  du2,     eu2,     ed2,     dp,      dsubcld,  &
                  jt,      maxg,    ideep,   1,       lengath,  &
                  nstep,   fracis,  ptend%q   )
   call t_stopf ('convtran1')
   call physics_update (state, tend, ptend, ztodt)


   end if    !( convection_scheme = 'Tiedtke')     !!(wh 2004.12)



997 continue
!---------------------------------------------------------------------------------
!-----------------------------------------------------------XXX
! ----stratiform cloud computer-----------------------------XXX
!---------stratiform cloud----------------------------------XXX
!-----------------------------------------------------------XXX

   if (RK_or_MG=='RK') then
!-------RK scheme------------sxj-2008-11-09----------------

!     if (masterproc)  write(6,*) "RK"
!
! cloud fraction after transport and convection,
! derive the relationship between rh and cld from
! the employed cloud scheme
!
   call t_startf('cldnrh')
   call cldnrh(lchnk,   ncol,                                &
               state%pmid,    state%t,   state%q(1,1,1),   state%omega, &
               cnt,     cnb,     cldn,    clc,     state%pdel,   &
               cmfmc,   landfrac,snowh,   concld,  cldst,    &
               ts,      state%pint(1,pverp),       zdu,  ocnfrac, &
               rhdfda,   rhu00 )
   call t_stopf('cldnrh')
!
! calculate the tendencies for moisture, temperature and cloud fraction
!
   rtdt = 1./ztodt
   qtend(:ncol,:pver) = (state%q(:ncol,:pver,1)       - qcwato(:ncol,:pver))*rtdt
   ttend(:ncol,:pver) = (state%t(:ncol,:pver)         - tcwato(:ncol,:pver))*rtdt
   lctend(:ncol,:pver) = (state%q(:ncol,:pver,ixcldw) - lcwato(:ncol,:pver))*rtdt
!
! strat condensation via prognostic cloud water
! calculate tendencies
!
   call t_startf('pcond')
   zero(:ncol,:pverp) = 0.
   call pcond (lchnk,   ncol, &
               state%t,   ttend,   state%q(1,1,1), qtend,       state%omega,     &
               state%q(1,1,ixcldw),state%pmid,     state%pdel,  cldn,     &
               qme,     nevapr,    prain,          rmelt,    &
               ztodt,   zero,      fwaut,          fsaut,       fracw,    &
               fsacw,   fsaci,     lctend,         rhdfda,      rhu00, icefrac)
   call t_stopf('pcond')
!
   call outfld('FWAUT',fwaut, pcols,lchnk)
   call outfld('FSAUT',fsaut, pcols,lchnk)
   call outfld('FRACW',fracw, pcols,lchnk)
   call outfld('FSACW',fsacw, pcols,lchnk)
   call outfld('FSACI',fsaci, pcols,lchnk)
!
! make it interactive
!
   do k = 1,pver
      do i = 1,ncol
         ptend%s(i,k)        = (qme(i,k) - nevapr(i,k))*latvap + rmelt(i,k)
         ptend%q(i,k,1)      =-(qme(i,k) - nevapr(i,k))
         ptend%q(i,k,ixcldw) = (qme(i,k) - prain(i,k))
         preclp(i) = preclp(i) + (prain(i,k)-nevapr(i,k))*state%pdel(i,k)/gravit
      end do
   end do
   ptend%name  = 'pcond'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.
   ptend%lq(ixcldw) = .TRUE.
   call physics_update (state, tend, ptend, ztodt)
!
! save off q and t after cloud water
!
   do k=1,pver
      qcwatn(:ncol,k) = state%q(:ncol,k,1)
      tcwatn(:ncol,k) = state%t(:ncol,k)
      lcwatn(:ncol,k) = state%q(:ncol,k,ixcldw)
   end do

    else if (RK_or_MG == 'MG') then

        ! calculate the tendencies for moisture, temperature and cloud fraction
        rtdt = 1./ztodt
        qtend(:ncol,:pver)  = (state%q(:ncol,:pver,1)     -qcwato(:ncol,:pver))*rtdt
        ttend(:ncol,:pver)  = (state%t(:ncol,:pver)       -tcwato(:ncol,:pver))*rtdt
        lctend(:ncol,:pver) = (state%q(:ncol,:pver,ixcldw)-lcwato(:ncol,:pver))*rtdt

        call stratiform_tend(state, ptend, ztodt, icefrac, landfrac, ocnfrac, &
                             snowh, dlf, cmfmc, cmfmc2, ts, sst, zdu, &
                             prec_str, snow_str, prec_sed, snow_sed, &
                             prec_pcw, snow_pcw, qtend, ttend, lctend, &
                             fice_MG, pbuf)

        call physics_update(state, tend, ptend, ztodt)
! prec_str(pcols)  ! [Total] sfc flux of precip from stratiform (m/s)
! preclp is in kg/m2/s
    preclp(:ncol)=prec_str(:ncol)*1000.
    itim = pbuf_old_tim_idx()  ! now itim =1
    if (itim>1) then
       write(*,*) " itim should be 1  in gamil now",itim
       call endrun
    endif
    ifld = pbuf_get_fld_idx('CLD')
    cldn(1:pcols,1:pver)= pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)


   do k=1,pver
      qcwatn(:ncol,k) = state%q(:ncol,k,1)
      tcwatn(:ncol,k) = state%t(:ncol,k)
      lcwatn(:ncol,k) = state%q(:ncol,k,ixcldw)
   end do
   !write(*,*) state%q(1,2,ixcldw),state%q(1,2,-1)  -
   !stop 22
   else
#ifdef SPMD
      CALL mpibarrier (mpicom)
      if (masterproc) write(6,*) "RK_or_MG must be RK or MG --tphysbc.F90"
#endif
      call endrun
   endif

!-----------------stratiform done  --sxj-------------@@@@@@@@@@@@@@@@@@@@@
!-----------------------------------------------------@@@@@@@@@@@@@@@@@@@@@@

   if (RK_or_MG=='MG')  fice=fice_MG !----------------sxj-------


!
!     Convective transport of all trace species except water vapor and
!     cloud water done here because we need to do the scavenging first
!     to determine the interstitial fraction.
!
   ptend%name  = 'convtran2'
   do m=2,ppcnst
      if (m /= ixcldw) ptend%lq(m) = .true.
   end do
   call t_startf ('convtran2')
   call convtran (lchnk,                                           &
                  ptend%lq,state%q, ppcnst,     mu2,     md2,      &
                  du2,     eu2,     ed2,        dp,      dsubcld,  &
                  jt,      maxg,    ideep,      1,       lengath,  &
                  nstep,   fracis,  ptend%q)
   call t_stopf ('convtran2')

   call physics_update (state, tend, ptend, ztodt)

   if (RK_or_MG=='RK') then
      call outfld('CMFDQR  ',cmfdqr, pcols,lchnk)
      call outfld('CLDST   ',cldst,  pcols,lchnk)
      call outfld('CNVCLD  ',clc,    pcols,lchnk)  !!!! sxj
      call outfld('CONCLD  ',concld, pcols,lchnk)
      call outfld('CME     ',qme,    pcols,lchnk)
      call outfld('PRAIN   ',prain,  pcols,lchnk)
      call outfld('EVAPR   ',nevapr, pcols,lchnk)
      call outfld('DQP     ',qc,      pcols,   lchnk  )
   endif
!
! and add the precip of pcond and cond together. Note
! since preclp is in kg/m2/s,  We need to renormalize to m/s
! like other precips
!
   precl(:ncol) = precl(:ncol) + preclp(:ncol)/1000.
!
! Compute rate of temperature change due to moist processes
!
   dtcond(:ncol,:) = (state%s(:ncol,:) - dtcond(:ncol,:))*rtdt
   call outfld('DTCOND  ',dtcond / cpair  ,pcols   ,lchnk   )
!
! Compute rate of constituent change due to moist processes
!
   dqcond(:ncol,:,:) = (state%q(:ncol,:,:) - dqcond(:ncol,:,:))*rtdt
   do m=1,ppcnst
      call outfld(dcconnam(m),dqcond(1,1,m),pcols   ,lchnk )
   end do
!
!===================================================
! Moist physical parameteriztions complete:
! send dynamical variables, and derived variables to history file
!===================================================
!
   call diag_dynvar (lchnk, ncol, state)
!
!===================================================
! Radiation computations
!===================================================
!
! Cosine solar zenith angle for current time step
!
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)
   call zenith (calday, clat, clon, coszrs, ncol)
!
! Compute liquid water paths
!
   tgcwp(:ncol) = 0.
   do k=1,pver
      do i = 1,ncol
         if (RK_or_MG=='RK') then
           gclwp2(i,k) = state%q(i,k,ixcldw)*state%pdel(i,k)/gravit*1000.0 ! Grid box liquid water path.
         endif
         if (RK_or_MG=='MG') then  ! cldice
           call cnst_get_ind('CLDICE',ix_cldice)
           call cnst_get_ind('CLDLIQ',ix_cldliq)
           !gclwp2(i,k) = state%q(i,k,ix_cldliq)*state%pdel(i,k)/gravit*1000.0
           gclwp2(i,k) = (state%q(i,k,ix_cldliq)+state%q(i,k,ix_cldice))*state%pdel(i,k)/gravit*1000.0
         endif
        ! write(*,*) ix_cldice,ix_cldliq

!!(2004)
!!(ljli0726)         gclwp2(i,k) = gclwp2(i,k)*1.67
          gclwp2(i,k) = gclwp2(i,k)   !!(ljli0821)
!!(ljli1215)                     gclwp2(i,k) = gclwp2(i,k)*1.25   !!(ljli1215)
!!(2004)
         clwp2(i,k) = gclwp2(i,k) / max(0.01_r8,cldn(i,k))               ! In-cloud liquid water path.
         tgcwp(i) = tgcwp(i) + gclwp2(i,k)
      end do
   end do

   call outfld('GCLDLWP  ',gclwp2,  pcols,lchnk)
   call outfld('ICLDLWP  ',clwp2,   pcols,lchnk)
   call outfld('TGCLDCWP ',tgcwp,  pcols,lchnk)

   if (dosw .or. dolw) then
!
! Compute cloud properties for input to radiation
!
      call t_startf('cldint')
      call virtem (ncol, pcols, pver, state%t, state%q(1,1,1), zvir, tvm)

      call cldint (lchnk, ncol, state%pmid, state%t, state%q(1,1,1), &
                   state%pint, state%lnpint, state%lnpmid, tvm, state%zi,  &
                   cldn, clwp, emis, effcld, landfrac,                    &
                   rel, rei, fice, state%pdel, tpw,                  &
                   hl, state%ps, nmxrgn, pmxrgn, clwp2)
!   real(r8), intent(in)  :: clwp2(pcols,pver)    ! prognostic cloud liquid water path
!   real(r8), intent(out) :: clwp(pcols,pver)     ! evaluate cloud liquid water path
!   real(r8), intent(out) :: emis(pcols,pver)     ! cloud emissivity
!   real(r8), intent(out) :: effcld(pcols,pver)   ! effective cloud=cld*emis
!   real(r8), intent(out) :: rel(pcols,pver)      ! effective drop radius (microns)
!   real(r8), intent(out) :: rei(pcols,pver)      ! ice effective drop size (microns)
!   real(r8), intent(out) :: fice(pcols,pver)     ! fractional ice content within cloud
!   real(r8), intent(out) :: tpw(pcols)           ! total precipitable water (in mm)  from q
!   real(r8), intent(out) :: hl(pcols)            ! liquid water scale height
!   real(r8), intent(out) :: pmxrgn(pcols,pverp) ! Maximum values of pressure for each
!
!###############################################################
   aer_indirect=.false.
!    aer_indirect=.true.
!#################################################################
    if ((RK_or_MG=='MG').and.(aer_indirect)) then
      ifld = pbuf_get_fld_idx('REL')
      !ifld = pbuf_get_fld_idx('REL_Fn')
      rel(1:ncol,1:pver)=pbuf(ifld)%fld_ptr(1,1:ncol,1:pver,lchnk, 1)
      ifld = pbuf_get_fld_idx('REI')
      rei(1:ncol,1:pver)=pbuf(ifld)%fld_ptr(1,1:ncol,1:pver,lchnk, 1)
      !ifld = pbuf_get_fld_idx('REL_FN')
      ! => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
      do k=1,pver
      do i=1,ncol
         !if (rel(i,k)<2.)   write(*,*) rel(i,k)
         !if (rel(i,k)>50.)  write(*,*) rel(i,k)
         !if (rei(i,k)<10.)  write(*,*) rei(i,k)
         !if (rei(i,k)>400.) write(*,*) rei(i,k)
         if (rel(i,k)<10.)  rel(i,k)=10.
         if (rel(i,k)>20.) rel(i,k)=20.
         if (rei(i,k)<10.) rei(i,k)=10.
         if (rei(i,k)>30.) rei(i,k)=30.
      enddo
      enddo
      call cldems(lchnk   ,ncol    ,clwp2   ,fice    ,rei     ,emis    )
       !intent(in) :: clwp2(pcols,pver)       ! cloud liquid water path (g/m**2)
       !intent(in) :: rei(pcols,pver)         ! ice effective drop size (microns)
      do k=1,pver
      do i=1,ncol
         effcld(i,k) = cldn(i,k)*emis(i,k)
      end do
      end do
    endif

      call t_stopf('cldint')
!
! Dump cloud field information to history tape buffer (diagnostics)
!
      call outfld('CLOUD   ',cldn,  pcols,lchnk)
      call outfld('EFFCLD  ',effcld, pcols,lchnk)
      call outfld('LWSH    ',hl,     pcols,lchnk)
!
! Compute in cloud ice mixing ratio and in cloud liquid mixing ratio
!
      if (RK_or_MG=='RK') then
      do k=1,pver
         do i = 1,ncol
            icimr(i,k) = state%q(i,k,ixcldw)*fice(i,k) / max(0.01_r8,cldn(i,k))
            icwmr(i,k) = state%q(i,k,ixcldw)*(1.-fice(i,k)) / max(0.01_r8,cldn(i,k))
         end do
      end do
      call outfld('ICIMR ',icimr,  pcols,lchnk)
      call outfld('ICWMR ',icwmr,  pcols,lchnk)
      call outfld('FICE  ',fice,   pcols,lchnk)
      endif
!
! Special diagnostic cloud water fields:
!
      call outfld('SETLWP  ',clwp,   pcols,lchnk)
!
! Output pure ice and pure liquid water paths
!
      tgiwp(:ncol) = 0.
! tgiwp(pcols)                      ! Vertically integrated ice water path
! tglwp(pcols)                      ! Vertically integrated liquid water path
      do k=1,pver
         tgiwp(:ncol) = tgiwp(:ncol) + gclwp2(:ncol,k)*fice(:ncol,k)
      end do

      tglwp(:ncol) = tgcwp(:ncol) - tgiwp(:ncol)

      call outfld ('TGCLDLWP',tglwp,  pcols,lchnk)
      call outfld ('TGCLDIWP',tgiwp,  pcols,lchnk)
!
! Complete radiation calculations
!
      call t_startf ('radctl')
!qrs(pcols,pver)            ! Shortwave heating rate
!qrl(pcols,pver)            ! Longwave  heating rate
!flwds(pcols)               ! Surface longwave down flux
!fsns(pcols)                   ! Surface solar absorbed flux
!fsnt(pcols)                   ! Net column abs solar flux at model top
!flns(pcols)                   ! Srf longwave cooling (up-down) flux
!flnt(pcols)                   ! Net outgoing lw flux at model top

      call radctl (lchnk, ncol, lwup, emis, state%pmid,             &
                   state%pint, state%lnpmid, state%lnpint, state%t, state%q,   &
                   cldn, clwp2, coszrs, asdir, asdif,               &
                   aldir, aldif, pmxrgn, nmxrgn, fsns, fsnt    ,flns    ,flnt    , &
                   qrs, qrl, flwds, rel, rei,                       &
                   fice, sols, soll, solsd, solld,                  &
                   landfrac, state%zm)
      call t_stopf ('radctl')
!intent(in) :: lwup(pcols)          ! Longwave up flux at surface
!intent(in) :: emis(pcols,pver)     ! Cloud emissivity
!intent(in) :: pmid(pcols,pver)     ! Model level pressures
!intent(in) :: pint(pcols,pverp)    ! Model interface pressures
!intent(in) :: pmln(pcols,pver)     ! Natural log of pmid
!intent(in) :: rel(pcols,pver)      ! liquid effective drop size (microns)
!intent(in) :: rei(pcols,pver)      ! ice effective drop size (microns)
!intent(in) :: fice(pcols,pver)     ! fractional ice content within cloud
!intent(in) :: piln(pcols,pverp)    ! Natural log of pint
!intent(in) :: t(pcols,pver)        ! Model level temperatures
!intent(in) :: qm1(pcols,pver,ppcnst) ! Specific humidity and tracers
!intent(in) :: cld(pcols,pver)      ! Fractional cloud cover
!intent(in) :: clwp(pcols,pver)     ! Cloud liquid water path
!intent(in) :: coszrs(pcols)        ! Cosine solar zenith angle
!intent(in) :: asdir(pcols)         ! albedo shortwave direct
!intent(in) :: asdif(pcols)         ! albedo shortwave diffuse
!intent(in) :: aldir(pcols)         ! albedo longwave direct
!intent(in) :: aldif(pcols)         ! albedo longwave diffuse
!intent(in) :: landfrac(pcols)      ! land fraction
!intent(in) :: zm(pcols,pver)       ! Height of midpoints (above surface)
!intent(inout) :: pmxrgn(pcols,pverp) ! Maximum values of pmid for each
!intent(inout) :: nmxrgn(pcols)     ! Number of maximally overlapped regions
!intent(out) :: fsns(pcols)          ! Surface absorbed solar flux
!intent(out) :: fsnt(pcols)          ! Net column abs solar flux at model top
!intent(out) :: flns(pcols)          ! Srf longwave cooling (up-down) flux
!intent(out) :: flnt(pcols)          ! Net outgoing lw flux at model top
!intent(out) :: sols(pcols)          ! Downward solar rad onto surface (sw direct)
!intent(out) :: soll(pcols)          ! Downward solar rad onto surface (lw direct)
!intent(out) :: solsd(pcols)         ! Downward solar rad onto surface (sw diffuse)
!intent(out) :: solld(pcols)         ! Downward solar rad onto surface (lw diffuse)
!intent(out) :: qrs(pcols,pver)      ! Solar heating rate
!intent(out) :: qrl(pcols,pver)      ! Longwave cooling rate
!intent(out) :: flwds(pcols)         ! Surface down longwave flux
!
! Cloud cover diagnostics
! radctl can change pmxrgn and nmxrgn so cldsav needs to follow
! radctl.
!
      call cldsav (lchnk, ncol, cldn, state%pmid, cltot, &
                   cllow, clmed, clhgh, nmxrgn, pmxrgn)
!
! Dump cloud field information to history tape buffer (diagnostics)
!
      call outfld('CLDTOT  ',cltot  ,pcols,lchnk)
      call outfld('CLDLOW  ',cllow  ,pcols,lchnk)
      call outfld('CLDMED  ',clmed  ,pcols,lchnk)
      call outfld('CLDHGH  ',clhgh  ,pcols,lchnk)
!                                                                          !!(wh)
! ISCCP Cloud Simulator                                                    !!(wh)
!                                                                          !!(wh)
      if (doisccp) then                                                    !!(wh)
                                                                           !!(wh)
         do k =1 , pver                                                    !!(wh)
            do i = 1, ncol                                                 !!(wh)
               cliqwp(i,k) = clwp2(i,k)*(1.0-fice(i,k))                    !!(wh)
               cicewp(i,k) = clwp2(i,k)*fice(i,k)                          !!(wh)
            end do                                                         !!(wh)
         end do                                                            !!(wh)
                                                                           !!(wh)
         call cloudsimulator_run(state, ts, concld, cldn,cliqwp, &         !!(wh 2005.01.28)
                                 cicewp, rel, rei, emis, coszrs  )         !!(wh following cam3.0)
      end if                                                               !!(wh)

   end if
!
! Compute net flux (for use in SLD energy fixer; not used in other dyn cores)
! Since fsns, fsnt, flns, and flnt are in the buffer, array values will be carried across
! timesteps when the radiation code is not invoked.
!
   do i=1,ncol
      tend%flx_net(i) = fsnt(i) - fsns(i) - flnt(i) + flns(i)
   end do
!
! Compute net radiative heating
!
   call radheat_net (state, ptend, qrl, qrs)
!
! Add radiation tendencies to cummulative model tendencies and update profiles
!
   call physics_update(state, tend, ptend, ztodt)
!
! Compute net surface radiative flux for use by surface temperature code.
! Note that units have already been converted to mks in RADCTL.  Since
! fsns and flwds are in the buffer, array values will be carried across
! timesteps when the radiation code is not invoked.
!
   srfrad(:ncol) = fsns(:ncol) + flwds(:ncol)
   call outfld('SRFRAD  ',srfrad,pcols,lchnk)
!
! determine whether precipitation, prec, is frozen (snow) or not
! by taking the mass-weighted temperature of the bottom most layers.
!
!  determine mass weighted average temperature in the bottom three
!  levels of the model (approximately the lowest 900 meters)
!
   dellow(:ncol) = 0.0
   tavg (:ncol) = 0.0
!
   do k=pver-2, pver
      dellow(:ncol) = dellow(:ncol) + state%pdel(:ncol,k)
      tavg(:ncol) = tavg(:ncol) + state%t(:ncol,k )*state%pdel(:ncol,k)
   end do
!
   tavg (:ncol) = tavg(:ncol)/dellow(:ncol)
!
   where (tavg(:ncol) > (tmelt-2.0) )
      precsc(:ncol) = 0.
      precsl(:ncol) = 0.
   elsewhere
      precsc(:ncol) = precc(:ncol)
      precsl(:ncol) = precl(:ncol)
   end where
!   if (RK_or_MG=='MG') then  !!! sxj-------
!      precsl(1:ncol) = snow_str(1:ncol)  !!!!!!!  lnd cannot deal with both rain and snow
!   endif
  prcsnw(:ncol) = precsc(:ncol) + precsl(:ncol)   ! total snowfall rate: needed by slab ocean model
!
! Precipitation intensity and frequency (where precip exceeds given threshold)
!
  where( precc(:ncol) >= precc_thresh )
     preccfrq(:ncol) = 1.0
     preccint(:ncol) = precc(:ncol)
  elsewhere
     preccint(:ncol) = 0.0
     preccfrq(:ncol) = 0.0
  end where
  where( precl(:ncol) >= precl_thresh )
     preclfrq(:ncol) = 1.0
     preclint(:ncol) = precl(:ncol)
  elsewhere
     preclint(:ncol) = 0.0
     preclfrq(:ncol) = 0.0
  end where
  preclint(:ncol) = preclint(:ncol)*(3600.*1000.) ! convert from m/sec to mm/hr
  preccint(:ncol) = preccint(:ncol)*(3600.*1000.) ! convert from m/sec to mm/hr
!
! Save atmospheric fields to force surface models
!
   call srfxfer (lchnk, ncol, state%ps, state%u(1,pver), state%v(1,pver),    &
                 state%t(1,pver), state%q(1,pver,1), state%exner(1,pver), state%zm(1,pver), &
                    state%pmid,      &
                 state%rpdel(1,pver))

!---------------------------------------------------------------------------------------
! Save history variables. These should move to the appropriate parameterization interface
!---------------------------------------------------------------------------------------

   call outfld('CMFMC   ',cmfmc   ,pcols   ,lchnk   )
   call outfld('CMFSL   ',cmfsl   ,pcols   ,lchnk   )
   call outfld('CMFLQ   ',cmflq   ,pcols   ,lchnk   )
   call outfld('PRECL   ',precl   ,pcols   ,lchnk       )
   call outfld('PRECC   ',precc   ,pcols   ,lchnk       )
   call outfld('PRECLINT',preclint,pcols   ,lchnk       )
   call outfld('PRECCINT',preccint,pcols   ,lchnk       )
   call outfld('PRECLFRQ',preclfrq,pcols   ,lchnk       )
   call outfld('PRECCFRQ',preccfrq,pcols   ,lchnk       )
   call outfld('PRECSL  ',precsl  ,pcols   ,lchnk       )
   call outfld('PRECSC  ',precsc  ,pcols   ,lchnk       )

! prect total
   prect(:ncol) = precc(:ncol) + precl(:ncol)
   call outfld('PRECT   ',prect   ,pcols   ,lchnk       )
   call outfld('PRECTMX ',prect   ,pcols   ,lchnk       )

#ifdef COUP_CSM
   call outfld('PRECLav ',precl   ,pcols   ,lchnk   )
   call outfld('PRECCav ',precc   ,pcols   ,lchnk   )
#endif
!
! Compute heating rate for dtheta/dt
!
   do k=1,pver
      do i=1,ncol
         ftem(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5/state%pmid(i,k))**cappa
      end do
   end do
   call outfld('HR      ',ftem    ,pcols   ,lchnk   )
!
! Convert mass fractions of non-water tracers back to mixing ratios.
! (Overwrite non-water portions of q3m1).
!
   if (ppcnst > 1) then
      call mf2mr (lchnk, ncol, state%q)
   end if

   return
 end subroutine tphysbc

