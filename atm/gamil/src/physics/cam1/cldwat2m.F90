#include<misc.h>   !    Use  SPMD

!---------------------------------------------------------------------------------
! Purpose:
!
!   CAM Interface for microphysics
!
! Author:
!
!   Andrew Gettelman and Hugh Morrison, December 2005
!
!   Description in Morrison and Gettelman, 2007. J. Climate, (submitted)  
!   for questions contact Hugh Morrison 
!   e-mail: morrison@ucar.edu
!   phone: 303-497-8916
!
!   SHI Xiangjun
!     bring this code to GAMIL 2008/11/12
!     conceal aerosol 2008/11/12
!     add aerosol 2009/03/08
!---------------------------------------------------------------------------------

module cldwat2m

    ! TAG-3
    ! debug modules
    use pmgrid,         only: masterproc, iam
#ifdef spmd           
    use mpishorthand,   only:mpicom
#endif

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use ppgrid,         only: pcols, pver
    use physconst,      only: gravit, rair, tmelt, cpair, rh2o, r_universal, mwh2o, rhoh2o
    use physconst,      only: latvap, latice
    use aerosol_index,  only: naer_all 
    !++ag v1.68  new error function
    use error_function, only: erf, erfc

    ! this needed for 'subroutine findsp1'
    use wv_saturation, only: estblf, hlatv, tmin, hlatf, rgasv, pcf, &
        cp, epsqs, ttrice,vqsatd

    implicit none

    private

    save

    public inimmc, mmicro_pcond, gamma, polysvp, derf  

    ! constants remaped
    real(r8) g             ! gravity
    real(r8) mw            ! molecular weight of water
    real(r8) r             ! dry air gas constant
    real(r8) rv            ! water vapor gas contstant
    real(r8) rr            ! universal gas constant
    real(r8) cpp           ! specific heat of dry air
    real(r8) rhow          ! density of liquid water
    real(r8) xxlv          ! latent heat of vaporization
    real(r8) xlf           ! latent heat of freezing
    real(r8) xxls          ! latent heat of sublimation

    ! from 'microconstants'
    real(r8) rhosn         ! bulk density snow
    real(r8) rhoi          ! bulk density ice

    real(r8) ac, bc, as, bs, ai, bi, ar, br  ! fall speed parameters 
    real(r8) ci,di         ! ice mass-diameter relation parameters
    real(r8) cs,ds         ! snow mass-diameter relation parameters
    real(r8) cr,dr         ! drop mass-diameter relation parameters
    real(r8) f1s,f2s       ! ventilation param for snow
    real(r8) Eii           ! collection efficiency aggregation of ice
    real(r8) Ecc           ! collection efficiency
    real(r8) Ecr           ! collection efficiency cloud droplets/rain
    real(r8) f1r,f2r       ! ventilation param for rain
    real(r8) DCS           ! autoconversion size threshold
    real(r8) qsmall        ! min mixing ratio 
    real(r8) bimm,aimm     ! immersion freezing
    real(r8) rhosu         ! typical 850mn air density
    real(r8) mi0           ! new crystal mass
    real(r8) rin           ! radius of contact nuclei
    real(r8) qcvar         ! 1/relative variance of sub-grid qc
    real(r8) pi

    ! maximum, minimum temperature for ice, liq
    real(r8) tmax_fice
    real(r8) tmin_fice

    ! parameters for snow/rain fraction for convective clouds
    real(r8), parameter :: tmax_fsnow = tmelt       ! max temperature for transition to convective snow
    real(r8), parameter :: tmin_fsnow = tmelt-5._r8 ! min temperature for transition to convective snow

    ! from 'microconstants, aerosol properties
    ! hm, note, block of code below no longer needed and can be removed
    real(r8) osm, vi, epsm, rhoa, map, ma, bact
    real(r8) rm1,sig1,nanew1,f11,f21  !mode 1
    real(r8) rm2,sig2,nanew2,f12,f22  !mode 2
    !..........................................................................

    !needed for findsp
    real(r8), private:: t0       ! Freezing temperature

    ! activate parameters

    integer, private:: psat
    parameter (psat=6) ! number of supersaturations to calc ccn concentration
    real(r8), private:: aten
    real(r8), private:: alogsig(naer_all) ! natl log of geometric standard dev of aerosol
    real(r8), private:: exp45logsig(naer_all)
    real(r8), private:: argfactor(naer_all)
    real(r8), private:: amcube(naer_all) ! cube of dry mode radius (m)
    real(r8), private:: smcrit(naer_all) ! critical supersatuation for activation
    real(r8), private:: lnsm(naer_all) ! ln(smcrit)
    real(r8), private:: alogten,alog2,alog3,alogaten
    real(r8), private, parameter :: supersat(psat)= & ! supersaturation (%) to determine ccn concentration
        (/0.02,0.05,0.1,0.2,0.5,1.0/)
    real(r8), private:: ccnfact(psat,naer_all)

    real(r8), private:: f1(naer_all),f2(naer_all) ! abdul-razzak functions of width
    real(r8), private:: third, sixth,zero
    real(r8), private:: sq2, sqpi

contains

    !===============================================================================

    subroutine inimmc

        !----------------------------------------------------------------------- 
        ! 
        ! Purpose: 
        ! initialize constants for the morrison microphysics
        ! called from stratiform.F90
        ! 
        ! Author: Andrew Gettelman Dec 2005
        ! 
        !-----------------------------------------------------------------------

        use pmgrid, only: plev, plevp
        use aerosol_index, only: naer_all
        use prescribed_aerosols, only: dryrad_aer, density_aer, hygro_aer, dispersion_aer, num_to_mass_aer

        integer k

        integer l,m
        real(r8) surften       ! surface tension of water w/respect to air (N/m)
        real(r8) super(psat)   ! supersaturation
        real(r8) arg

        ! hm modify to use my error function

	real(r8) derf

        !declarations for morrison codes (transforms variable names)

        g= gravit                  !gravity
        mw = mwh2o / 1000._r8      !molecular weight of water
        r= rair       	      !Dry air Gas constant: note units(phys_constants are in J/K/kmol)
        rv= rh2o                   !water vapor gas contstant
        rr = r_universal           !universal gas constant
        cpp = cpair                 !specific heat of dry air
        rhow = rhoh2o              !density of liquid water

        !NOTE:
        ! latent heats should probably be fixed with temperature 
        ! for energy conservation with the rest of the model
        ! (this looks like a +/- 3 or 4% effect, but will mess up energy balance)

        xxlv = latvap         ! latent heat vaporization
        xlf = latice          ! latent heat freezing
        xxls = xxlv + xlf     ! latent heat of sublimation

        ! from microconstants

        ! parameters below from Reisner et al. (1998)
        ! density parameters (kg/m3)

        rhosn = 100._r8    ! bulk density snow
        rhoi = 500._r8     ! bulk density ice
        rhow = 1000._r8    ! bulk density liquid


        ! fall speed parameters, V = aD^b
        ! V is in m/s

        ! droplets
	ac = 3.e7_r8
	bc = 2._r8

        ! snow
	as = 11.72_r8
	bs = 0.41_r8

        ! cloud ice
	ai = 700._r8
	bi = 1._r8

        ! rain
	ar = 841.99667_r8
	br = 0.8_r8

        ! particle mass-diameter relationship
        ! currently we assume spherical particles for cloud ice/snow
        ! m = cD^d

        pi= 3.1415927_r8

        ! cloud ice mass-diameter relationship

	ci = rhoi*pi/6._r8
	di = 3._r8

        ! snow mass-diameter relationship

	cs = rhosn*pi/6._r8
	ds = 3._r8

        ! drop mass-diameter relationship

	cr = rhow*pi/6._r8
	dr = 3._r8

        ! ventilation parameters for snow
        ! hall and prupacher

	f1s = 0.86_r8
	f2s = 0.28_r8

        ! collection efficiency, aggregation of cloud ice and snow

        !ljli20090830        Eii = 0.1_r8
        Eii = 0.2_r8   !ljli20090830

        ! collection efficiency, accretion of cloud water by rain

	Ecr = 1.0_r8

        ! ventilation constants for rain

	f1r = 0.78_r8
	f2r = 0.32_r8

        ! autoconversion size threshold for cloud ice to snow (m)

	Dcs = 200.e-6_r8

        ! smallest mixing ratio considered in microphysics

	qsmall = 1.e-18_r8  

        ! immersion freezing parameters, bigg 1953

	bimm = 100._r8
	aimm = 0.66_r8

        ! typical air density at 850 mb

	rhosu = 85000._r8/(rair * tmelt)

        ! mass of new crystal due to aerosol freezing and growth (kg)

	mi0 = 4._r8/3._r8*pi*rhoi*(10.e-6_r8)*(10.e-6_r8)*(10.e-6_r8)

        ! radius of contact nuclei aerosol (m)

        rin = 0.1e-6_r8

        ! 1 / relative variance of sub-grid cloud water distribution
        ! see morrison and gettelman, 2007, J. Climate for details

	qcvar = 1._r8

        ! max/min temperature parameters for ice deposition/bergeron process
        tmax_fice = 268.15_r8
        tmin_fice = 238.15_r8

        ! freezing temperature
	t0=273.15_r8

        !-----------------------------------------------------------------------
        ! aerosol properties needed by the scheme are specified below.
        ! for now assume 2 mode lognormal aerosol, with soluble portion
        ! consisting of ammonium sulfate
!!!!!!!!!!!!!!following block of code no longer used!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        mw = 0.018_r8           ! molecular weight water
        osm = 1._r8             ! osmotic coefficient
        vi = 3._r8              ! number of ions produced durin dissociation in h20
        epsm = 0.5_r8           ! soluble fraction
        rhoa = 1777._r8         ! bulk density of aerosol
        map = 0.132_r8          ! molecular weight of aerosol
        ma = 0.0284_r8          ! molecular weight of 'air'
        rr = 8.3187_r8          ! universal gas constant
        bact = vi*osm*epsm*mw*rhoa/(map*rhow)  ! 'big B' in Abdul-Razzak et al. (1998)

        ! mode 1

        rm1 = 0.05e-6_r8        ! geometric mean radius (m)
        sig1 = 1.5_r8           ! standard deviation
        nanew1 = 500.e6_r8      ! total concentration (m-3)
        f11 = 0.5_r8*exp(2.5_r8*(log(sig1))**2)  ! fi factor in Abdul-Razak and Ghan (2000)
        f21 = 1._r8+0.25_r8*log(sig1)  ! gi factor in Abdul-Razak and Ghan (2000)

        ! mode 2
        ! definitions as in mode 1

        rm2 = 1.e-6_r8
        sig2 = 2._r8
        nanew2 = 0.1e6_r8
        f12 = 0.5_r8*exp(2.5_r8*(log(sig2))**2)
        f22 = 1._r8+0.25_r8*log(sig2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! set parameters for droplet activation, following abdul-razzak and ghan 2000, JGR

        !      mathematical constants

        zero=0._r8
        third=1./3._r8
        sixth=1./6._r8
        sq2=sqrt(2._r8)
        pi=4._r8*atan(1.0_r8)
        sqpi=sqrt(pi)


        surften=0.076_r8
        aten=2.*mwh2o*surften/(r_universal*t0*rhoh2o)
        alogaten=log(aten)
        alog2=log(2._r8)
        alog3=log(3._r8)
        super(:)=0.01*supersat(:)

        do m=1,naer_all
            !         use only if width of size distribution is prescribed
            alogsig(m)=log(dispersion_aer(m))
            exp45logsig(m)=exp(4.5*alogsig(m)*alogsig(m))
            argfactor(m)=2./(3.*sqrt(2.)*alogsig(m))
            f1(m)=0.5*exp(2.5*alogsig(m)*alogsig(m))
            f2(m)=1.+0.25*alogsig(m)
            !         use only if mode radius of size distribution is prescribed
            amcube(m)=3./(4.*pi*exp45logsig(m)*density_aer(m)*num_to_mass_aer(m))
            !         use only if only one component per mode
            if(hygro_aer(m).gt.1.e-10) then
                smcrit(m)=2.*aten*sqrt(aten/(27.*hygro_aer(m)*amcube(m)))
            else
                smcrit(m)=100.
            endif
            lnsm(m)=log(smcrit(m))

            do l=1,psat  ! psat=6
                arg=argfactor(m)*log(smcrit(m)/super(l))
                if(arg<2)then
                    if(arg<-2)then
                        ccnfact(l,m)=1.e-6
                    else
                        !  ++ag revise error function...
                        ccnfact(l,m)=1.e-6_r8*0.5_r8*erfc(arg)
                        !  ++old error function pre V1.68
                        ! 		   ccnfact(l,m)=1.e-6*0.5*(1.-derf(arg))
                        !  original ghan code
                        ccnfact(l,m)=1.e-6*0.5*ERFC(arg)
                    endif
                else
                    ccnfact(l,m) = 0.
                endif
            enddo

        end do

        ! CALL mpibarrier (mpicom)         !! for debug-test   
        ! write(6,*) iam,"inimmc_MG done"  !!  sxj-
        ! call endrun
        !write(*,*) alogsig
        !write(*,*) dispersion_aer
        !write(*,*) "cldwat2.F90 inimmc done ln363"
        !call endrun
        return
    end subroutine inimmc


    subroutine mmicro_pcond (lchnk,ncol,deltatin,tn,ttend,qn,qtend,cwtend,qc,qi, &
                                !nc,ni,p,pdel,cldn,aer_mmr,rhdfda,rhu00,fice,tlat,qvlat, &  !! 
                                !nc,ni,p,pdel,cldn,         rhdfda,rhu00,fice,tlat,qvlat, &  !! sxj--20081113  aerosol not read
	nc,ni,p,pdel,cldn,aer_mmr,rhdfda,rhu00,fice,tlat,qvlat, &    !! sxj-20090311 take aerol into accout
        qctend,qitend,nctend,nitend,effc,effc_fn,effi,prect,preci,kkvh,wsub, &
        ! hm added variables for trop_mozart
        nevapr,evapsnow,prain,prodsnow,cmeout,landfrac,snowh)
        !Author: Hugh Morrison Dec 2005, NCAR
        ! e-mail: morrison@ucar.edu
        ! phone: 303-497-8916
        use wv_saturation,       only: vqsatd
        use aerosol_index,       only: naer_all, idxSUL,n_aer=>naer   !! naer_all=13 naer=11 !sxj
        use prescribed_aerosols, only: dryrad_aer, density_aer, hygro_aer, dispersion_aer, num_to_mass_aer
        use history,         only: outfld

        ! variable declarations

        implicit none

        ! input variables
        integer, intent(in) :: lchnk
        integer, intent(in) :: ncol

        real(r8), intent(in) :: tn(pcols,pver)       ! input temperature (K)
        real(r8), intent(in) :: ttend(pcols,pver)    ! non-microphysical temperature tendency (K/s)
        real(r8), intent(in) :: qn(pcols,pver)       ! input h20 vapor mixing ratio (kg/kg)
        real(r8), intent(in) :: qtend(pcols,pver)    ! non-microphysical qv tendency (1/s)
        real(r8), intent(in) :: cwtend(pcols,pver)   ! non-microphysical tendency of cloud water (1/s)	
        ! note: all input cloud variables are grid-averaged
        real(r8), intent(inout) :: qc(pcols,pver)    ! cloud water mixing ratio (kg/kg)
        real(r8), intent(inout) :: qi(pcols,pver)    ! cloud ice mixing ratio (kg/kg)
        real(r8), intent(inout) :: nc(pcols,pver)    ! cloud water number conc (1/kg)
        real(r8), intent(inout) :: ni(pcols,pver)    ! cloud ice number conc (1/kg)
        !note sxj  qc,qi,nc,ni   vary.  
        real(r8), intent(in) :: p(pcols,pver)        ! air pressure (pa)
        real(r8), intent(in) :: pdel(pcols,pver)     ! pressure difference across level (pa)
        real(r8), intent(in) :: cldn(pcols,pver)     ! cloud fraction
        real(r8), intent(in) :: rhdfda(pcols,pver)   ! dA/dRH
        real(r8), intent(in) :: rhu00(pcols,pver)    ! threshold rh for cloud
        !real(r8), intent(in) :: fice(pcols,pver)     ! fraction of ice/liquid
        real(r8), intent(inout) :: fice(pcols,pver)     ! fraction of ice/liquid
        real(r8), intent(in) :: deltatin             ! time step (s)
        real(r8), intent(in) :: aer_mmr(pcols,pver,naer_all) ! aerosol mass mixing ratio
        real(r8), intent(in) :: kkvh(pcols,pver+1)           ! vertical eddy diff coef (m2 s-1)
        real(r8), intent(in)  :: landfrac(pcols)      !20101114ljli  land fraction (fraction)^M
        real(r8), intent(in) :: snowh(pcols)

        ! output variables

        real(r8), intent(out) :: tlat(pcols,pver)    ! latent heating rate     (K/s)
        real(r8), intent(out) :: qvlat(pcols,pver)   ! microphysical tendency qv (1/s)
        real(r8), intent(out) :: qctend(pcols,pver)  ! microphysical tendency qc (1/s) 
        real(r8), intent(out) :: qitend(pcols,pver)  ! microphysical tendency qi (1/s)
        real(r8), intent(out) :: nctend(pcols,pver)  ! microphysical tendency nc (1/(kg*s))
        real(r8), intent(out) :: nitend(pcols,pver)  ! microphysical tendency ni (1/(kg*s))
        real(r8), intent(out) :: effc(pcols,pver)    ! droplet effective radius (micron)
        real(r8), intent(out) :: effc_fn(pcols,pver) ! droplet effective radius, assuming nc = 1.e8 kg-1
        real(r8), intent(out) :: effi(pcols,pver)    ! cloud ice effective radius (micron)
        real(r8), intent(out) :: prect(pcols)        ! surface precip rate (m/s)
        real(r8), intent(out) :: preci(pcols)        ! cloud ice/snow precip rate (m/s)
        real(r8), intent(out) :: wsub(pcols,pver)    ! diagnosed sub-grid vertical velocity st. dev. (m/s)

        ! hm added variables for trop-mozart

        real(r8), intent(out) :: nevapr(pcols,pver)  ! evaporation rate of rain
        real(r8), intent(out) :: evapsnow(pcols,pver)! sublimation rate of snow
        real(r8), intent(out) :: prain(pcols,pver)   ! production of rain
        real(r8), intent(out) :: prodsnow(pcols,pver)! production of snow
        real(r8), intent(out) :: cmeout(pcols,pver)  ! evap/sub of cloud

        ! hm add 3/19 new vars for sub-step solution
        ! note, '1' variables below are the same as above, but needed for calculating
        ! average rates over time step, since sub-stepping is used here

        real(r8) :: t1(pcols,pver)
        real(r8) :: q1(pcols,pver)
        real(r8) :: qc1(pcols,pver)
        real(r8) :: qi1(pcols,pver)
        real(r8) :: nc1(pcols,pver)
        real(r8) :: ni1(pcols,pver)
        real(r8) :: tlat1(pcols,pver)
        real(r8) :: qvlat1(pcols,pver)
        real(r8) :: qctend1(pcols,pver)
        real(r8) :: qitend1(pcols,pver)
        real(r8) :: nctend1(pcols,pver)
        real(r8) :: nitend1(pcols,pver)
        real(r8) :: prect1(pcols)
        real(r8) :: preci1(pcols)
        real(r8) :: deltat     ! sub-time step (s)

        ! trop-mozart variables for sub-stepping
        real(r8) :: nevapr1(pcols,pver)
        real(r8) :: evapsnow1(pcols,pver)
        real(r8) :: prain1(pcols,pver)
        real(r8) :: prodsnow1(pcols,pver)
        real(r8) :: cmeout1(pcols,pver)

        ! local workspace
        ! all units mks unless otherwise stated
        real(r8) :: omsm    ! number near unity for round-off issues
        real(r8) :: dto2    ! dt/2 (s)
        real(r8) :: mincld  ! minimum allowed cloud fraction
        real(r8) :: q(pcols,pver) ! water vapor mixing ratio (kg/kg)
        real(r8) :: t(pcols,pver) ! temperature (K)
        real(r8) :: rho(pcols,pver) ! air density (kg m-3)
        real(r8)  :: aer_act(pcols,pver) ! aerosol activation number (#/cm3)
        real(r8) :: dv(pcols,pver)  ! diffusivity of water vapor in air
        real(r8) :: mu(pcols,pver)  ! viscocity of air
        real(r8) :: sc(pcols,pver)  ! schmidt number
        real(r8) :: kap(pcols,pver) ! thermal conductivity of air
        real(r8) :: dap(pcols,pver) ! effecvtive diffusivity of contact ice nuclei
        real(r8) :: cldmax(pcols,pver) ! precip fraction assuming maximum overlap
        real(r8) :: cldm(pcols,pver)   ! cloud fraction
        real(r8) :: icwc(pcols)    ! in cloud water content (liquid+ice)
        real(r8) :: calpha(pcols)  ! parameter for cond/evap (Zhang et al. 2003)
        real(r8) :: cbeta(pcols) ! parameter for cond/evap (Zhang et al. 2003)
        real(r8) :: cbetah(pcols) ! parameter for cond/evap (Zhang et al. 2003)
        real(r8) :: cgamma(pcols) ! parameter for cond/evap (Zhang et al. 2003)
        real(r8) :: cgamah(pcols) ! parameter for cond/evap (Zhang et al. 2003)
        real(r8) :: rcgama(pcols) ! parameter for cond/evap (Zhang et al. 2003)
        real(r8) :: cmec1(pcols) ! parameter for cond/evap (Zhang et al. 2003)
        real(r8) :: cmec2(pcols) ! parameter for cond/evap (Zhang et al. 2003)
        real(r8) :: cmec3(pcols) ! parameter for cond/evap (Zhang et al. 2003)
        real(r8) :: cmec4(pcols) ! parameter for cond/evap (Zhang et al. 2003)
        real(r8) :: qtmp ! dummy qv 
        real(r8) :: dum  ! temporary dummy variable
        real(r8) :: cme(pcols,pver)  ! total (liquid+ice) cond/evap rate of cloud
        real(r8) :: cmei(pcols,pver) ! dep/sublimation rate of cloud ice
        real(r8) :: cmel(pcols,pver) ! cond/evap rate of cloud liquid
        real(r8) :: cwml(pcols,pver) ! cloud water mixing ratio
        real(r8) :: cwmi(pcols,pver) ! cloud ice mixing ratio
        real(r8) :: nnuccd(pver)   ! ice nucleation rate from deposition/cond.-freezing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! variables below are no longer used
        real(r8) :: sigvl
        real(r8) :: aact
        real(r8) :: alpha
        real(r8) :: gamm
        real(r8) :: gg
        real(r8) :: psi
        real(r8) :: eta1
        real(r8) :: eta2
        real(r8) :: sm1
        real(r8) :: sm2
        real(r8) :: dum1
        real(r8) :: dum2
        real(r8) :: smax
        real(r8) :: uu1
        real(r8) :: uu2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(r8) :: npccn(pver)     ! droplet activation rate
        real(r8) :: qcic(pcols,pver) ! in-cloud cloud liquid mixing ratio
        real(r8) :: qiic(pcols,pver) ! in-cloud cloud ice mixing ratio
        real(r8) :: qniic(pcols,pver) ! in-precip snow mixing ratio
        real(r8) :: qric(pcols,pver) ! in-precip rain mixing ratio
        real(r8) :: ncic(pcols,pver) ! in-cloud droplet number conc
        real(r8) :: niic(pcols,pver) ! in-cloud cloud ice number conc
        real(r8) :: nsic(pcols,pver) ! in-precip snow number conc
        real(r8) :: nric(pcols,pver) ! in-precip rain number conc
        real(r8) :: lami(pver) ! slope of cloud ice size distr
        real(r8) :: n0i(pver) ! intercept of cloud ice size distr
        real(r8) :: lamc(pver) ! slope of cloud liquid size distr
        real(r8) :: n0c(pver) ! intercept of cloud liquid size distr
        real(r8) :: lams(pver) ! slope of snow size distr
        real(r8) :: n0s(pver) ! intercept of snow size distr
        real(r8) :: lamr(pver) ! slope of rain size distr
        real(r8) :: n0r(pver) ! intercept of rain size distr
        real(r8) :: cdist1(pver) ! size distr parameter to calculate droplet freezing
        real(r8) :: pgam(pver) ! spectral width parameter of droplet size distr
        real(r8) :: lammax  ! maximum allowed slope of size distr
        real(r8) :: lammin  ! minimum allowed slope of size distr
        real(r8) :: nacnt   ! number conc of contact ice nuclei
        real(r8) :: mnuccc(pver) ! mixing ratio tendency due to freezing of cloud water
        real(r8) :: nnuccc(pver) ! number conc tendency due to freezing of cloud water
        real(r8) :: prc(pver) ! qc tendency due to autoconversion of cloud droplets
        real(r8) :: nprc(pver) ! number conc tendency due to autoconversion of cloud droplets
        real(r8) :: nprc1(pver) ! qr tendency due to autoconversion of cloud droplets
        real(r8) :: nsagg(pver) ! ns tendency due to self-aggregation of snow
        real(r8) :: dc0  ! mean size droplet size distr
        real(r8) :: ds0  ! mean size snow size distr (area weighted)
        real(r8) :: eci  ! collection efficiency for riming of snow by droplets
        real(r8) :: psacws(pver) ! mixing rat tendency due to collection of droplets by snow
        real(r8) :: npsacws(pver) ! number conc tendency due to collection of droplets by snow
        real(r8) :: uni ! number-weighted cloud ice fallspeed
        real(r8) :: umi ! mass-weighted cloud ice fallspeed
        real(r8) :: uns(pver) ! number-weighted snow fallspeed
        real(r8) :: ums(pver) ! mass-weighted snow fallspeed
        real(r8) :: unr(pver) ! number-weighted rain fallspeed
        real(r8) :: umr(pver) ! mass-weighted rain fallspeed
        real(r8) :: unc ! number-weighted cloud droplet fallspeed
        real(r8) :: umc ! mass-weighted cloud droplet fallspeed
        real(r8) :: pracs(pver) ! mixing rat tendency due to collection of rain	by snow
        real(r8) :: npracs(pver) ! number conc tendency due to collection of rain by snow
        real(r8) :: mnuccr(pver) ! mixing rat tendency due to freezing of rain
        real(r8) :: nnuccr(pver) ! number conc tendency due to freezing of rain
        real(r8) :: pra(pver) ! mixing rat tendnency due to accretion of droplets by rain
        real(r8) :: npra(pver) ! nc tendnency due to accretion of droplets by rain
        real(r8) :: nragg(pver) ! nr tendency due to self-collection of rain
        real(r8) :: prci(pver) ! mixing rat tendency due to autoconversion of cloud ice to snow
        real(r8) :: nprci(pver) ! number conc tendency due to autoconversion of cloud ice to snow
        real(r8) :: prai(pver) ! mixing rat tendency due to accretion of cloud ice by snow
        real(r8) :: nprai(pver) ! number conc tendency due to accretion of cloud ice by snow
        real(r8) :: qvs ! liquid saturation vapor mixing ratio
        real(r8) :: qvi ! ice saturation vapor mixing ratio
        real(r8) :: dqsdt ! change of sat vapor mixing ratio with temperature
        real(r8) :: dqsidt ! change of ice sat vapor mixing ratio with temperature
        real(r8) :: ab ! correction factor for rain evap to account for latent heat
        real(r8) :: qclr ! water vapor mixing ratio in clear air
        real(r8) :: abi ! correction factor for snow sublimation to account for latent heat
        real(r8) :: epss ! 1/ sat relaxation timescale for snow
        real(r8) :: epsr ! 1/ sat relaxation timescale for rain
        real(r8) :: pre(pver) ! rain mixing rat tendency due to evaporation
        real(r8) :: prds(pver) ! snow mixing rat tendency due to sublimation
        real(r8) :: qce ! dummy qc for conservation check
        real(r8) :: qie ! dummy qi for conservation check
        real(r8) :: nce ! dummy nc for conservation check
        real(r8) :: nie ! dummy ni for conservation check
        real(r8) :: ratio ! parameter for conservation check
        real(r8) :: dumc(pcols,pver) ! dummy in-cloud qc
        real(r8) :: dumnc(pcols,pver) ! dummy in-cloud nc
        real(r8) :: dumi(pcols,pver) ! dummy in-cloud qi
        real(r8) :: dumni(pcols,pver) ! dummy in-cloud ni
        real(r8) :: dums(pcols,pver) ! dummy in-cloud snow mixing rat
        real(r8) :: dumns(pcols,pver) ! dummy in-cloud snow number conc
        real(r8) :: dumr(pcols,pver) ! dummy in-cloud rain mixing rat
        real(r8) :: dumnr(pcols,pver) ! dummy in-cloud rain number conc
        ! below are parameters for cloud water and cloud ice sedimentation calculations
        real(r8) :: fr(pver)
        real(r8) :: fnr(pver)
        real(r8) :: fc(pver)
        real(r8) :: fnc(pver)
        real(r8) :: fi(pver)
        real(r8) :: fni(pver)
        real(r8) :: fs(pver)
        real(r8) :: fns(pver)
        real(r8) :: faloutr(pver)
        real(r8) :: faloutnr(pver)
        real(r8) :: faloutc(pver)
        real(r8) :: faloutnc(pver)
        real(r8) :: falouti(pver)
        real(r8) :: faloutni(pver)
        real(r8) :: falouts(pver)
        real(r8) :: faloutns(pver)
        real(r8) :: faltndr
        real(r8) :: faltndnr
        real(r8) :: faltndc
        real(r8) :: faltndnc
        real(r8) :: faltndi
        real(r8) :: faltndni
        real(r8) :: faltnds
        real(r8) :: faltndns
        real(r8) :: faltndqie
        real(r8) :: faltndqce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(r8) :: relhum(pcols) ! relative humidity
        real(r8) :: csigma(pcols) ! parameter for cond/evap of cloud water/ice
        real(r8) :: rgvm ! max fallspeed for all species
        real(r8) :: arn(pcols,pver) ! air density corrected rain fallspeed parameter
        real(r8) :: asn(pcols,pver) ! air density corrected snow fallspeed parameter
        real(r8) :: acn(pcols,pver) ! air density corrected cloud droplet fallspeed parameter
        real(r8) :: ain(pcols,pver) ! air density corrected cloud ice fallspeed parameter
        real(r8) :: nsubi(pver) ! evaporation of cloud ice number
        real(r8) :: nsubc(pver) ! evaporation of droplet number
        real(r8) :: nsubs(pver) ! evaporation of snow number
        real(r8) :: nsubr(pver) ! evaporation of rain number
        real(r8) :: mtime ! factor to account for droplet activation timescale
        real(r8) :: dz(pcols,pver) ! height difference across model vertical level

        !++ag 1.90 new fice variable
        real(r8) :: nfice(pcols,pver)

        ! returns from function/subroutine calls
        real(r8) :: tsp(pcols,pver)      ! saturation temp (K)
        real(r8) :: qsp(pcols,pver)      ! saturation mixing ratio (kg/kg)
        real(r8) :: qs(pcols)            ! liquid-ice weighted sat mixing rat (kg/kg)
        real(r8) :: es(pcols)            ! liquid-ice weighted sat vapor press (pa)
        real(r8) :: esl(pcols,pver)      ! liquid sat vapor pressure (pa)
        real(r8) :: esi(pcols,pver)      ! ice sat vapor pressure (pa)
        real(r8) :: gammas(pcols)        ! parameter for cond/evap of cloud water

        ! sum of source/sink terms for diagnostic precip

        real(r8) :: qnitend(pcols,pver) ! snow mixing ratio source/sink term
        real(r8) :: nstend(pcols,pver)  ! snow number concentration source/sink term
        real(r8) :: qrtend(pcols,pver) ! rain mixing ratio source/sink term
        real(r8) :: nrtend(pcols,pver)  ! rain number concentration source/sink term
        real(r8) :: qrtot ! vertically-integrated rain mixing rat source/sink term
        real(r8) :: nrtot ! vertically-integrated rain number conc source/sink term
        real(r8) :: qstot ! vertically-integrated snow mixing rat source/sink term
        real(r8) :: nstot ! vertically-integrated snow number conc source/sink term

        ! new terms for Bergeron process

        real(r8) :: dumnnuc ! provisional ice nucleation rate (for calculating bergeron)
        real(r8) :: ninew  ! provisional cloud ice number conc (for calculating bergeron)
        real(r8) :: qinew ! provisional cloud ice mixing ratio (for calculating bergeron)
        real(r8) :: qvl  ! liquid sat mixing ratio   
        real(r8) :: epsi ! 1/ sat relaxation timecale for cloud ice
        real(r8) :: prd ! provisional deposition rate of cloud ice at water sat 
        real(r8) :: berg(pcols,pver) ! mixing rat tendency due to bergeron process for cloud ice
        real(r8) :: bergs(pver) ! mixing rat tendency due to bergeron process for snow

        ! new aerosol variables
        real(r8) naermod(naer_all) ! aerosol number concentration (/m3)
	real(r8) naer2(pcols,pver,naer_all)   ! new aerosol number concentration (/m3)

        real(r8) maerosol(1,naer_all)   ! aerosol mass conc (kg/m3)
	real(r8) naer(pcols)
        real(r8) ccn(pcols,pver,psat)        ! number conc of aerosols activated at supersat
        character*8, parameter :: ccn_name(psat)=(/'CCN1','CCN2','CCN3','CCN4','CCN5','CCN6'/)

        ! diagnostic rain/snow for output to history
        ! values are in-precip (local) !!!!

        real(r8) :: qrout(pcols,pver) ! rain mixing ratio (kg/kg)
        real(r8) :: qsout(pcols,pver) ! snow mixing ratio (kg/kg)
        real(r8) :: nrout(pcols,pver) ! rain number concentration (1/m3)
        real(r8) :: nsout(pcols,pver) ! snow number concentration (1/m3)

        !++ag averageed rain/snow for history
        real(r8) :: qrout2(pcols,pver)
        real(r8) :: qsout2(pcols,pver)
        real(r8) :: nrout2(pcols,pver)
        real(r8) :: nsout2(pcols,pver)
        real(r8) :: freqs(pcols,pver)
        real(r8) :: freqr(pcols,pver)
        real(r8) :: dumfice
        real(r8) :: drout2(pcols,pver) ! mean rain particle diameter (m)
        real(r8) :: dsout2(pcols,pver) ! mean snow particle diameter (m)

        ! hm add 3/19/07, ice nucleation, droplet activation
        real(r8) :: dum2i(pcols,pver) ! number conc of ice nuclei available (1/kg)
        real(r8) :: dum2l(pcols,pver) ! number conc of CCN (1/kg)
        real(r8) :: ncmax
        real(r8) :: nimax

        ! loop array variables
        integer i,k,nstep,n, l
        integer ii,kk, m  , ntype(naer_all)
        integer ftrue  ! integer used to determine cloud depth

        ! hm add 3/19/07, new loop variables for sub-step solution
        integer iter,it,ltrue(pcols)

        ! variabels to check for RH after rain evap

        real(r8) :: esn
        real(r8) :: qsn
        real(r8) :: ttmp
        !CCCCCCCCCCCCC
        !    write(*,*) "aer_mmr(pcols,pver,naer_all) "      
        !     write(*,*) aer_mmr


        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !++ag assign variable deltat for sub-stepping...
        deltat=deltatin 	!deltatin             ! time step (s)
        ! parameters for scheme
        omsm=0.99999_r8
        dto2=0.5_r8*deltat
        mincld=0.0001_r8
        ! initialize multi-level fields
        do i=1,ncol
            do k=1,pver
                q(i,k)=qn(i,k) ! input h20 vapor mixing ratio
                t(i,k)=tn(i,k) ! input temperature (K)
	    end do
        end do
        ! calculat sub-grid vertical velocity wsub ####
        do i=1,ncol
	    ftrue=0  !ftrue -- integer used to determine cloud depth
            do k=1,pver
                ! get sub-grid vertical velocity from diff coef.
                ! following morrison et al. 2005, JAS
                ! assume mixing length of 30 m       
                dum=(kkvh(i,k)+kkvh(i,k+1))/2._r8/30._r8
                ! kkvh(pcols,pver+1)   ! vertical eddy diff coef (m2 s-1) 
                ! use maximum sub-grid vertical vel of 10 m/s
                dum=min(dum,10._r8)
                ! set wsub to value at current vertical level
                dum=max(dum,0.10_r8)  !!! sxj-add MG2008 a minimun sub-grid vertical velocity of 10cm/s is assumped.
                wsub(i,k)=dum
                if ((wsub(i,k)<0.10).or.(wsub(i,k)>10.)) then 
                    write(*,*) "wsub:",wsub(i,k)     !!!!!sxj
                    call endrun
                endif
            end do
        end do

        ! initialize time-varying parameters
        do k=1,pver
            do i=1,ncol
                rho(i,k)=p(i,k)/(r*t(i,k))
                dv(i,k) = 8.794E-5_r8*t(i,k)**1.81_r8/p(i,k) 
                ! dv(pcols,pver) - diffusivity of water vapor in air
                mu(i,k) = 1.496E-6_r8*t(i,k)**1.5_r8/ &
                    (t(i,k)+120._r8)
                ! mu(pcols,pver) - viscosity of air	
                sc(i,k) = mu(i,k)/(rho(i,k)*dv(i,k))
                ! sc(pcols,pver) - schmidt number
                kap(i,k) = 1.414e3_r8*1.496e-6_r8*t(i,k)** &
                    1.5_r8/(t(i,k)+120._r8) 
                ! kap(pcols,pver)- thermal conductivity of air

                ! air density adjustment for fallspeed parameters
                ! hm added 11/18/06, add air density correction factor to the
                ! power of 0.54 following Heymsfield and Bansemer 2006

                ! ac,bc,as,bs,ai,bi,ar,br-fall speed parameters 
                !arn(pcols,pver) ! air density corrected rain fallspeed parameter
                !asn(pcols,pver) ! air density corrected snow fallspeed parameter
                !acn(pcols,pver) ! air density corrected cloud droplet fallspeed parameter
                !ain(pcols,pver) ! air density corrected cloud ice fallspeed parameter
                ! rhosu-typical 850mn air density
                arn(i,k)=ar*(rhosu/rho(i,k))**0.54  
                asn(i,k)=as*(rhosu/rho(i,k))**0.54
                acn(i,k)=ac*(rhosu/rho(i,k))**0.54
                ain(i,k)=ai*(rhosu/rho(i,k))**0.54
                ! typical air density at 850 mb  -- rhosu = 85000._r8/(rair * tmelt) 

                ! mean free path
                dum = 7.37_r8*t(i,k)/(288._r8*10._r8*p(i,k))/100._r8
                ! effective diffusivity based on Brownian collection

                dap(i,k) = 4._r8*pi*1.38e-23_r8*t(i,k)*(1._r8+dum/rin)/ &
                    (6._r8*pi*rin*mu(i,k))
                ! dap---effecvtive diffusivity of contact ice nuclei
                ! pi ---pi=4._r8*atan(1.0_r8)
                ! rin---radius of contact nuclei aerosol (m):::: rin = 0.1e-6_r8

                ! get dz from dp and hydrostatic approx
                ! keep dz positive (define as layer k-1 - layer k)
                dz(i,k)= pdel(i,k)/(rho(i,k)*g)
            end do
        end do

        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! calculate condensation based on cloud fraction, t,q tendencies
        ! this is uses same logic as old pcond (zhang et al. 2003)
        do i=1,ncol
            do k=1,pver
                ! store original variables for sub-stepping
                t1(i,k) = t(i,k)
                q1(i,k) = q(i,k)
                qc1(i,k) = qc(i,k)
                qi1(i,k) = qi(i,k)
                nc1(i,k) = nc(i,k)
                ni1(i,k) = ni(i,k)
                ! initialize tendencies to zero
                tlat1(i,k)=0._r8
                qvlat1(i,k)=0._r8
                ! tlat(pcols,pver)  -- latent heating rate     (K/s)
                ! qvlat(pcols,pver) -- microphysical tendency qv (1/s)
                qctend1(i,k)=0._r8
                qitend1(i,k)=0._r8
                nctend1(i,k)=0._r8
                nitend1(i,k)=0._r8
                ! initialize precip output
                qrout(i,k)=0._r8
                qsout(i,k)=0._r8
                nrout(i,k)=0._r8
                nsout(i,k)=0._r8
                ! hm add, initialize variables for trop_mozart
                nevapr(i,k) = 0._r8
                evapsnow(i,k) = 0._r8
                prain(i,k) = 0._r8
                prodsnow(i,k) = 0._r8
                cmeout(i,k) = 0._r8
            end do
        end do
        ! initialize precip fraction and output tendencies
        do k=1,pver
            do i=1,ncol
                cldmax(i,k)=mincld
            end do
        end do
        ! initialize avg precip rate
        do i=1,ncol
	    prect1(i)=0._r8
	    preci1(i)=0._r8
            !prect(pcols) -- surface precip rate (m/s)
            !preci(pcols) -- cloud ice/snow precip rate (m/s)
        end do



        call findsp1 (lchnk,ncol,qn,tn,p,tsp,qsp)
        ! tsp(pcols,pver) -- saturation temp (K)
        ! qsp(pcols,pver) -- saturation mixing ratio (kg/kg) 
        do k=1,pver
            ! find wet bulk temperature and saturation value for provisional t and q without
            ! condensation
            call vqsatd(t(1,k),p(1,k),es,qs,gammas,ncol)
            ! es(len) -- saturation vapor pressure
            ! qs(len) -- saturation specific humidity
            ! gam(len) -- (l/cp)*(d(qs)/dt)
            do i=1,ncol
                ! esl(pcols,pver) - liquid sat vapor pressure (pa)
                ! esi(pcols,pver) - ice sat vapor pressure (pa)
                esl(i,k)=polysvp(t(i,k),0)
                esi(i,k)=polysvp(t(i,k),1)
                relhum(i)=q(i,k)/qs(i)
                ! get cloud fraction, check for minimum
                cldm(i,k)=max(cldn(i,k),mincld)
                ! calculate provisional ice nucleation!!!!!!!!!@@@@
                ! this is done prior to condensation/deposition calculation
                ! in order to calculate bergeron process
                if (t(i,k).lt.268.15_r8) then
                    ! fletcher curve
                    !             dum= 0.01_r8*exp(0.6_r8*(273.15_r8-t(i,k)))
                    ! cooper curve (factor of 1000 is to convert from L-1 to m-3)
                    dum=0.005_r8*exp(0.304_r8*(273.15_r8-t(i,k)))*1000._r8
                    ! put limit on number of nucleated crystals, set to number at T=-30 C
                    ! fletcher curve
                    !           dum=min(dum,6.56e5_r8)/rho(i,k)  ! convert from m-3 to kg-1
                    ! cooper (limit to value at -35 C)
                    dum=min(dum,208.9e3_r8)/rho(i,k)  ! convert from m-3 to kg-1
                    dumnnuc=(dum-ni(i,k)/cldm(i,k))/deltat*cldm(i,k)
                    dumnnuc=max(dumnnuc,0._r8)
                    ! get provisional ni and qi after nucleation in order to calculate
                    ! Bergeron process below
                    ninew=ni(i,k)+dumnnuc*deltat
                    qinew=qi(i,k)+dumnnuc*deltat*mi0 
                    !  ninew  - provisional cloud ice number conc (for calculating bergeron)
                    !  qinew -- provisional cloud ice mixing ratio (for calculating bergeron)
                    !  (inout) :: ni(pcols,pver)    ! cloud ice number conc (1/kg)
                    !  mass of new crystal due to aerosol freezing and growth (kg) 
                    ! ::::mi0 = 4._r8/3._r8*pi*rhoi*(10.e-6_r8)*(10.e-6_r8)*(10.e-6_r8)
                else
                    ninew=ni(i,k)
                    qinew=qi(i,k)
                end if
                !----------@@@@@@@@@@@@@
                ! define coefficients for cond/evap of cloud water/ice calculations
                ! cond/evap closure follows that of Zhang et al. 2003, JGR
                ! provsional in-cloud water content for calculating cme
                icwc(i)=(qc(i,k)+qi(i,k))/cldm(i,k)
                ! xxlv = latvap    ----   latent heat vaporization
                calpha(i) = 1.0_r8/qs(i)
                cbeta(i)  = q(i,k)/qs(i)**2*gammas(i)*cpp/xxlv
                ! ! gammas -- (l/cp)*(d(qs)/dt)
                cbetah(i) = 1.0_r8/qs(i)*gammas(i)*cpp/xxlv
                cgamma(i) = calpha(i)+xxlv*cbeta(i)/cpp
                cgamah(i) = calpha(i)+xxlv*cbetah(i)/cpp
                rcgama(i) = cgamma(i)/cgamah(i)
                if (relhum(i).gt. rhu00(i,k)) then
                    if (relhum(i).ge.0.999_r8.or.cldm(i,k).ge.0.999_r8) then
                        cme(i,k)=(calpha(i)*qtend(i,k)- &
                            cbetah(i)*ttend(i,k))/cgamah(i)
                        !if ((i==13).and.(k==24))    write(*,*) "1"
                    else
                        !if ((i==13).and.(k==24))    write(*,*) "2"
                        csigma(i)=1.0_r8/(rhdfda(i,k)+cgamma(i)*icwc(i))
                        cmec1(i) =(1.0_r8-cldm(i,k))*csigma(i)*rhdfda(i,k)
                        cmec2(i) =cldm(i,k)*calpha(i)/cgamah(i)+ &
                            (1.0_r8-rcgama(i)*cldm(i,k))*csigma(i)*calpha(i)*icwc(i)
                        cmec3(i) =cldm(i,k)*cbetah(i)/cgamah(i)+ &
                            (cbeta(i)-rcgama(i)*cldm(i,k)*cbetah(i))*csigma(i)*icwc(i)
                        cmec4(i) =csigma(i)*cgamma(i)*icwc(i)
                        cme(i,k) =-cmec1(i)*cwtend(i,k)+cmec2(i)*qtend(i,k) &
                            -cmec3(i)*ttend(i,k) !!! sxj---????????no evaporation term		
                    end if
                    ! need to evaporate
                else if (qc(i,k)+qi(i,k).gt.0.0_r8) then
                    !if ((i==13).and.(k==24))    write(*,*) "3"
                    cme(i,k)=-min(max(0._r8,qsp(i,k)-q(i,k)),qc(i,k)+qi(i,k))/deltat
                else
                    !if ((i==13).and.(k==24))    write(*,*) "4"
                    cme(i,k)=0.0_r8
	            cmel(i,k)=0._r8
	            cmei(i,k)=0._r8
                end if

                !if ((i==13).and.(k==24)) then
                !   write(*,*) "cmel(i,k)-927-cmei(i,k)"
                !   write(*,*)  cmel(i,k)  ,cmei(i,k) ,cme(i,k)
                !endif 
                ! make sure excess supersaturation is removed
                qtmp=q(i,k)-cme(i,k)*deltat
                if (qtmp .gt. qsp(i,k)) then
                    cme(i,k)=cme(i,k)+(qtmp-qsp(i,k))/deltat
                end if
                ! ensure that condensation does not reduce grid mean rh below rhu00
                if (cme(i,k).gt.0.0_r8.and.relhum(i).gt.rhu00(i,k)) &
                    cme(i,k)=min(cme(i,k),(q(i,k)-qs(i)*rhu00(i,k))/deltat)
                ! for condensation
                ! for T < tmin_fice, assume all new condensate is ice
                ! for T > tmin_fice, put all new condensate into liquid temporarily,
                ! then calculate transfer to ice through Bergeron process
                ! note: this approach assumes that ice/liquid is mixed throughout the cloudy
                ! portion of the grid cell


                ! make sure to initialize bergeron process to zero   
                ! separate cme   and calculate berg ( vaper deposition expense previous cloud water)--sxj--
                berg(i,k)=0._r8
                if (cme(i,k).ge.0._r8)   then
                    if (t(i,k).ge.tmin_fice) then
                        if (t(i,k).lt.273.15)    then
                            ! calculate Bergeron process
                            qvi=0.622_r8*esi(i,k)/(p(i,k)-esi(i,k))
                            qvl=0.622_r8*esl(i,k)/(p(i,k)-esl(i,k))
                            dqsidt =  xxls*qvi/(rv*t(i,k)**2)
                            abi = 1._r8+dqsidt*xxls/cpp
                            ! get ice size distribution parameters
                            ! get in-cloud qi and ni after nucleation
                            qiic(i,k)=qinew/cldm(i,k)
                            niic(i,k)=ninew/cldm(i,k)
                            if (qiic(i,k).ge.qsmall) then
                                ! smallest mixing ratio considered in microphysics  qsmall = 1.e-18_r8 
                                ! ! m = cD^d  cloud ice mass-diameter relationship --ci = rhoi*pi/6._r8 ----di = 3._r8
                                lami(k) = (gamma(1._r8+di)*ci* &
                                    niic(i,k)/qiic(i,k))**(1._r8/di)
                                n0i(k) = niic(i,k)*lami(k)
                                ! check for slope (mean cloud ice diameter 10e-6 ~ 400e-6)
                                lammax = 1._r8/10.e-6_r8
                                lammin = 1._r8/(2._r8*dcs)
                                ! autoconversion size threshold for cloud ice to snow (m)	Dcs = 200.e-6_r8
                                ! adjust vars
                                if (lami(k).lt.lammin) then
                                    lami(k) = lammin
                                    n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*gamma(1._r8+di))
                                else if (lami(k).gt.lammax) then
                                    lami(k) = lammax
                                    n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*gamma(1._r8+di))
                                end if
                                epsi = 2._r8*pi*n0i(k)*rho(i,k)*Dv(i,k) &   !!!!? ---paper wrong-
                                    /(lami(k)*lami(k))                   !!!!?--sxj----
                                prd = epsi*(qvl-qvi)/abi                    !!!!? --- code right
                            else
                                prd = 0._r8
                            end if
                            ! multiply by cloud fraction
                            prd = prd*cldm(i,k)
                            ! get fraction of new condensate converted to ice
                            dum=0._r8
                            if (cme(i,k).ne.0._r8) then
                                dum=prd/cme(i,k)
                            endif
                            dum=max(dum,0._r8)
                            dum=min(dum,1._r8)
                            cmei(i,k)=dum*cme(i,k)
                            cmel(i,k)=(1._r8-dum)*cme(i,k)
                            ! transfer of existing cloud liquid to ice
                            if (prd.gt.cme(i,k)) then
                                berg(i,k)=min(prd-cme(i,k),qc(i,k)/deltat*omsm)
                            end if
                        else ! t>273.15
                            cmel(i,k)=cme(i,k)
                            cmei(i,k)=0._r8
                        end if ! t>273.15
                    else ! t,tmin_fice
                        cmel(i,k)=0._r8
                        cmei(i,k)=cme(i,k)
                    end if ! t,tmin_fice
                    ! for evaporation, evaporate liquid first, then cloud ice if 
                    ! ice is present 
	        else  ! (cme(i,k)<0._r8) 
                    if (t(i,k) > 273.15_r8) then
                        cmei(i,k)=0.0_r8
                        cmel(i,k)=cme(i,k)
                    else
                        cmel(i,k)=max(cme(i,k),-qc(i,k)/deltat)
                        cmei(i,k)=min(cme(i,k)-cmel(i,k),0._r8)
                    end if
                end if ! (cme(i,k)<0._r8) 
                ! evaporation should not exceed available water
                if (cmel(i,k).lt.-qc(i,k)/deltat) &
                    cmel(i,k)=-qc(i,k)/deltat
                ! sublimation should not exceed available ice
                if (cmei(i,k).lt.-qi(i,k)/deltat) &
                    cmei(i,k)=-qi(i,k)/deltat
                ! limit cmel,cmei due for roundoff error
                cmel(i,k)=cmel(i,k)*omsm
                cmei(i,k)=cmei(i,k)*omsm


                if (qi(i,k).ge.qsmall.or.cmei(i,k)+berg(i,k).ge.qsmall) then
                    ! cooper curve (factor of 1000 is to convert from L-1 to m-3)
                    !  dum2i(pcols,pver) ------------number conc of ice nuclei available (1/kg)
                    dum2i(i,k)=0.005_r8*exp(0.304_r8*(273.15_r8-t(i,k)))*1000._r8
                    ! put limit on number of nucleated crystals, set to number at T=-30 C
                    ! cooper (limit to value at -35 C)
                    dum2i(i,k)=min(dum2i(i,k),208.9e3_r8)/rho(i,k)  ! convert from m-3 to kg-1
                else
                    dum2i(i,k)=0._r8
                end if

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!
                !-------------------------aerosol activation
                !-------------------------------------------------
                if (qc(i,k).ge.qsmall.or.cmel(i,k).ge.qsmall) then
                    ! assign minimum value for sub-grid vertical velocity, 0.1 m/s
                    wsub(i,k)=max(0.10_r8,wsub(i,k))
                    !    0.1<wsub<10
                    ! more general treatment with hooks into simulated aerosol
                    !do m=1,naer_all !sxj--------
                    do m=1,n_aer
                        ntype(m)=1
                        maerosol(1,m)=aer_mmr(i,k,m)*rho(i,k)
                        ! set number nucleated for sulfate based on Lohmann et al 2000 (JGR) eq 2
                        !           Na=340.*(massSO4)^0.58  where Na=cm-3 and massSO4=ug/m3
                        ! convert units to Na [m-3] and SO4 [kgm-3]
                        !   Na(m-3)= 1.e6 cm3 m-3 Na(cm-3)=340. *(massSO4[kg/m3]*1.e9ug/kg)^0.58
                        !   or Na(m-3)= 1.e6* 340.*(1.e9ug/kg)^0.58 * (massSO4[kg/m3])^0.58 
                        if(m .eq. idxSUL) then
                            !	   naermod(m)= 5.64259e13_r8 * maerosol(1,m)**0.58   
                            naer2(i,k,m)= 5.64259e13_r8 * maerosol(1,m)**0.58   
                            !   naer2(pcols,pver,naer_all)  ---  new aerosol number concentration (/m3)
                        else
                            !	   naermod(m)=maerosol(1,m)*num_to_mass_aer(m)
                            naer2(i,k,m)=maerosol(1,m)*num_to_mass_aer(m)
                            !  num_to_mass_aer(naer_all) ! ratio of number concentration to mass concentration (#/kg)
                        endif
                    enddo
                    ! get droplet activation rate- -! out  --dum2 
                    call activate(wsub(i,k),t(i,k),rho(i,k), &
                                ! naer2(i,k,:), 1,ntype, naer_all,naer_all, maerosol,        &  !! sxj---
                        naer2(i,k,1:n_aer), 1,ntype(1:n_aer), n_aer,n_aer, maerosol(1,1:n_aer),  &  !! 
                                !dispersion_aer,hygro_aer, density_aer, dum2)
                        dispersion_aer(1:n_aer),hygro_aer(1:n_aer), density_aer(1:n_aer), dum2)

                    if ((landfrac(i)>=0.1).and. (snowh(i) <= 0.000001_r8)) then 
                        dum2l(i,k) = max(dum2,40e6/rho(i,k))
                    else
                        dum2l(i,k) = max(dum2,10e6/rho(i,k))
                    end if  !add by Lijuan Li according to SXJ
                    ! dum2l(i,k) = nc(i,k)      !!sxj add20081113  transient.. on activation !!
                    !write(*,*) dum2
                    !write(*,*) "cldwat2m.F90-line1047  after activation"
                    !call endrun
                else
                    dum2l(i,k) = 0._r8
                end if
                aer_act(i,k)=dum2l(i,k)*rho(i,k)*1.e-6  !sxj--
            enddo ! i loop
        enddo ! k loop
        ! dum2l(pcols,pver) ! number conc of CCN (1/kg)
        ! aer_act (#/cm3)
        call outfld('aer_act ', aer_act , pcols, lchnk   )  !!--- sxj


        !---------
        do i=1,ncol
	    ltrue(i)=0 
	    do k=1,pver
                ! hm add 3/19/07 skip microphysical calculations if no cloud water
                if (qc(i,k).ge.qsmall.or.qi(i,k).ge.qsmall.or. &
                    cmel(i,k).ge.qsmall.or.cmei(i,k).ge.qsmall) ltrue(i)=1
	    end do
        end do

        iter = 2
        ! get sub-step time step
        deltat=deltat/real(iter)
!!!! skip calculations if no cloud water
        do i=1,ncol
	    if (ltrue(i).eq.0) then
                do k=1,pver
                    tlat(i,k)=0._r8
                    qvlat(i,k)=0._r8
                    qctend(i,k)=0._r8
                    qitend(i,k)=0._r8
                    qnitend(i,k)=0._r8
                    qrtend(i,k)=0._r8
                    nctend(i,k)=0._r8
                    nitend(i,k)=0._r8
                    nrtend(i,k)=0._r8
                    nstend(i,k)=0._r8
                    prect(i)=0._r8
                    preci(i)=0._r8
                    qniic(i,k)=0._r8
                    qric(i,k)=0._r8
                    nsic(i,k)=0._r8
                    nric(i,k)=0._r8
                end do
                goto 300
	    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@@
            !----begin sub-step------------------------------
            do it=1,iter
                do k=1,pver
                    ! initialize sub-step microphysical tendencies
                    tlat(i,k)=0._r8
                    qvlat(i,k)=0._r8
                    qctend(i,k)=0._r8
                    qitend(i,k)=0._r8
                    qnitend(i,k)=0._r8
                    qrtend(i,k)=0._r8
                    nctend(i,k)=0._r8
                    nitend(i,k)=0._r8
                    nrtend(i,k)=0._r8
                    nstend(i,k)=0._r8
                    ! initialize diagnostic precipitation to zero
                    qniic(i,k)=0._r8
                    qric(i,k)=0._r8
                    nsic(i,k)=0._r8
                    nric(i,k)=0._r8
                    ! hm add, initialize variables for trop_mozart
                    nevapr1(i,k) = 0._r8
                    evapsnow1(i,k) = 0._r8
                    prain1(i,k) = 0._r8
                    prodsnow1(i,k) = 0._r8
                    cmeout1(i,k) = 0._r8
                end do
                ! begin new i,k loop, calculate new cldmax after adjustment to cldm above
                ! initialize vertically-integrated rain and snow tendencies
                qrtot = 0._r8
                nrtot = 0._r8
                qstot = 0._r8
                nstot = 0._r8
                ! initialize precip at surface
                prect(i)=0._r8
                preci(i)=0._r8
                do k=1,pver
                    ! set cwml and cwmi to current qc and qi
                    cwml(i,k)=qc(i,k)
                    cwmi(i,k)=qi(i,k)
                    ! initialize precip fallspeeds to zero
                    ums(k)=0._r8 
                    uns(k)=0._r8 
                    umr(k)=0._r8 
                    unr(k)=0._r8
                    ! calculate precip fraction based on maximum overlap assumption
                    if (k.eq.1) then
                        cldmax(i,k)=cldm(i,k)
                    else
                        ! hm add sep 6, 2006, if rain or snow mix ratio is smaller than
                        ! threshold, then set cldmax to cloud fraction at current level
                        if (qric(i,k-1).ge.qsmall.or.qniic(i,k-1).ge.qsmall) then
                            !qric(pcols,pver) ! in-precip rain mixing ratio
                            !qniic(pcols,pver) ! in-precip snow mixing ratio
                            cldmax(i,k)=max(cldmax(i,k-1),cldm(i,k))
                        else
                            cldmax(i,k)=cldm(i,k)
                        end if
                    end if
                    ! decrease in number concentration due to sublimation/evap
                    ! divide by cloud fraction to get in-cloud decrease
                    ! don't reduce Nc due to bergeron process
                    !nsubi(pver) ! evaporation of cloud ice number
                    !nsubc(pver) ! evaporation of droplet number
                    if (cmei(i,k) < 0._r8) then
                        nsubi(k)=cmei(i,k)/qi(i,k)*ni(i,k)/cldm(i,k)
                    else
                        nsubi(k)=0._r8
                    end if
                    if (cmel(i,k) < 0._r8) then
                        nsubc(k)=cmel(i,k)/qc(i,k)*nc(i,k)/cldm(i,k)
                    else
                        nsubc(k)=0._r8
                    end if
                    !----------
                    ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
                    ! for microphysical process calculations
                    ! units are kg/kg for mixing ratio, 1/kg for number conc
                    ! limit in-cloud values to 0.005 kg/kg  !????????--sxj-- 

                    qcic(i,k)=min(cwml(i,k)/cldm(i,k),5.e-3_r8)
                    qiic(i,k)=min(cwmi(i,k)/cldm(i,k),5.e-3_r8)
                    ncic(i,k)=max(nc(i,k)/cldm(i,k),0._r8)
                    niic(i,k)=max(ni(i,k)/cldm(i,k),0._r8)
                    if (qc(i,k)+(cmel(i,k)-berg(i,k))*deltat.lt.qsmall) then
                        qcic(i,k)=0._r8
                        ncic(i,k)=0._r8
                        if (qc(i,k)+(cmel(i,k)-berg(i,k))*deltat.lt.0._r8) then
                            if (cmel(i,k).lt.0._r8) then
                                dum=-cmel(i,k)+berg(i,k)
                                cmel(i,k)=cmel(i,k)*qc(i,k)/deltat/dum*omsm
                                berg(i,k)=berg(i,k)*qc(i,k)/deltat/dum*omsm
                            else
                                berg(i,k)=qc(i,k)/deltat*omsm
                            end if
                        end if
                    end if
                    if (qi(i,k)+(cmei(i,k)+berg(i,k))*deltat.lt.qsmall) then
                        qiic(i,k)=0._r8
                        niic(i,k)=0._r8
                        if (qi(i,k)+(cmei(i,k)+berg(i,k))*deltat.lt.0._r8) then
                            cmei(i,k)=(-qi(i,k)/deltat-berg(i,k))*omsm
                        end if
                    end if
                    ! add to cme output
                    cmeout1(i,k) = cmeout1(i,k)+cmel(i,k)+cmei(i,k)     !!!

                    !----------
                    ! ice nucleation (based on cooper curve as function of T)
                    ! calculate ice nucleation if cloud ice is present and T < -5 C
                    ! dum2i(pcols,pver) ------------number conc of ice nuclei available (1/kg)
                    ! dum2i(i,k)=0.005_r8*exp(0.304_r8*(273.15_r8-t(i,k)))*1000._r8
                    if (qiic(i,k).ge.qsmall.and.t(i,k).lt.268.15_r8) then
                        nnuccd(k)=(dum2i(i,k)-ni(i,k)/cldm(i,k))/deltat*cldm(i,k)
                        nnuccd(k)=max(nnuccd(k),0._r8)
                        nimax = dum2i(i,k)*cldm(i,k)
                    else
                        nnuccd(k)=0._r8
                        nimax = 0._r8
                    end if
                    !----------
                    ! droplet activation
                    ! calculate potential for droplet activation if cloud water is present
                    ! formulation from Abdul-Razzak and Ghan (2000) and Abdul-Razzak et al. (1998), AR98
                    if (qcic(i,k).ge.qsmall) then
                        ! assume aerosols already activated are equal number of existing droplets for simplicity
                        ! multiply by cloud fraction to obtain grid-average tendency
                        npccn(k) = (dum2l(i,k)-nc(i,k)/cldm(i,k))/deltat*cldm(i,k)
                        ! make sure number activated > 0
                        npccn(k) = max(0._r8,npccn(k))
                        ncmax = dum2l(i,k)*cldm(i,k)
                    else
                        npccn(k)=0._r8
                        ncmax = 0._r8
                    end if
                    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! get size distribution parameters based on in-cloud cloud water/ice 
                    ! these calculations also ensure consistency between number and mixing ratio
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    !----------
                    ! cloud ice
                    if (qiic(i,k).ge.qsmall) then
                        ! add upper limit to in-cloud number concentration to prevent numerical error
                        niic(i,k)=min(niic(i,k),qiic(i,k)*1.e20_r8)  !  ####?? 
                        lami(k)  =(gamma(1._r8+di)*ci* &
                            niic(i,k)/qiic(i,k))**(1._r8/di)
                        n0i(k)   = niic(i,k)*lami(k)
                        ! check for slope
                        lammax = 1._r8/10.e-6_r8     ! 
                        lammin = 1._r8/(2._r8*dcs)   ! 
                        ! adjust vars
                        if (lami(k).lt.lammin) then
                            lami(k) = lammin
                            n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*gamma(1._r8+di))
                            niic(i,k) = n0i(k)/lami(k)
                            !di = 3._r8
                        else if (lami(k).gt.lammax) then
                            lami(k) = lammax
                            n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*gamma(1._r8+di))
                            niic(i,k) = n0i(k)/lami(k)
                        end if
                    else
                        lami(k) = 0._r8               
                        n0i(k) = 0._r8                
                    end if
                    !----------
                    !cloud droplet
                    if (qcic(i,k).ge.qsmall) then
                        ! add upper limit to in-cloud number concentration to prevent numerical error
                        ncic(i,k)=min(ncic(i,k),qcic(i,k)*1.e20_r8)
                        ! get pgam from fit to observations of martin et al. 1994
                        pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8/rho(i,k))+0.2714_r8
                        pgam(k)=1._r8/(pgam(k)**2)-1._r8
                        pgam(k)=max(pgam(k),2._r8)
                        pgam(k)=min(pgam(k),15._r8)
                        ! calculate lamc
                        lamc(k) = (pi/6._r8*rhow*ncic(i,k)*gamma(pgam(k)+4._r8)/ &
                            (qcic(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
                        ! lammin, 40 micron diameter max mean size
                        lammin = (pgam(k)+1._r8)/50.e-6_r8
                        lammax = (pgam(k)+1._r8)/2.e-6_r8
                        if (lamc(k).lt.lammin) then
                            lamc(k) = lammin
                            ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)* &
                                gamma(pgam(k)+1._r8)/ &
                                (pi*rhow*gamma(pgam(k)+4._r8))
                        else if (lamc(k).gt.lammax) then
                            lamc(k) = lammax
                            ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)* &
                                gamma(pgam(k)+1._r8)/ &
                                (pi*rhow*gamma(pgam(k)+4._r8))
                        end if
                        ! parameter to calculate droplet freezing
                        cdist1(k) = ncic(i,k)/gamma(pgam(k)+1._r8)   
                        ! cdist1(pver) -- size distr parameter to calculate droplet freezing
                    else
                        lamc(k) = 0._r8
                        cdist1(k) = 0._r8
                    end if
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! begin micropysical process calculations 
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    !--------------
                    ! autoconversion of cloud liquid water to rain+++++
                    ! formula from Khrouditnov and Kogan (2000)
                    ! minimum qc of 1 x 10^-8 prevents floating point error
                    if (qcic(i,k).ge.1.e-8_r8) then
                        ! nprc is increase in rain number conc due to autoconversion
                        ! nprc1 is decrease in cloud droplet conc due to autoconversion
                        ! assume exponential sub-grid distribution of qc, resulting in additional
                        ! factor related to qcvar below
                        ! nprc(pver) ------ number conc tendency due to autoconversion of cloud droplets:: rain 
                        ! nprc1(pver) ----- qr tendency due to autoconversion of cloud droplets         :: droplet 
                        ! qcvar---1
                        prc(k) = gamma(qcvar+2.47_r8)/(gamma(qcvar)*qcvar**2.47_r8)* &
                            1350._r8*qcic(i,k)**2.47_r8* &
                            (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
                        nprc(k) = prc(k)/(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
                        nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))
                    else
                        prc(k)=0._r8
                        nprc(k)=0._r8
                        nprc1(k)=0._r8
                    end if
                    ! add autoconversion to precip from above to get provisional rain mixing ratio
                    ! and number concentration (qric and nric) 
                    ! 0.45 m/s is fallspeed of new rain drop (80 micron diameter)
                    dum=0.45_r8
                    dum1=0.45_r8
                    ! hm modify 6/12
                    if (k.eq.1) then
                        qric(i,k)=prc(k)*cldm(i,k)*dz(i,k)/cldmax(i,k)/dum
                        nric(i,k)=nprc(k)*cldm(i,k)*dz(i,k)/cldmax(i,k)/dum
                    else
                        if (qric(i,k-1).ge.qsmall) then
                            dum=umr(k-1)
                            dum1=unr(k-1)
                            !				  if (dum<0) then          !--------
                            !				  write(*,*) "cldwat2.f90 line1350"
                            !				  call endrun
                            !				  endif
                        end if
                        ! precip fallspeeds --umr(k)-unr(k)
                        ! hm add 4/17/06, no autoconversion of rain number if rain/snow falling from above
                        ! this assumes that new drizzle drops formed by autoconversion are rapidly collected
                        ! by the existing rain/snow particles from above
                        if (qric(i,k-1).ge.1.e-9_r8.or.qniic(i,k-1).ge.1.e-9_r8) then
                            nprc(k)=0._r8
                        end if
                        !   K>1
                        ! pra(pver) ! mixing rat tendnency due to accretion of droplets by rain  first step =0
                        ! prc is increase in rain  due to autoconversion
                        ! pre(pver) ! rain mixing rat tendency due to evaporation   first step =0
                        ! pracs(pver) ! mixing rat tendency due to collection of rain	by snow  first step =0
                        ! mnuccr(pver) ! mixing rat tendency due to freezing of rain  first step =0
                        qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*cldmax(i,k-1)+ &
                            (rho(i,k)*dz(i,k)* &
                            ((pra(k-1)+prc(k))*cldm(i,k)+(pre(k-1)-pracs(k-1)-mnuccr(k-1))*cldmax(i,k))))&
                            /(dum*rho(i,k)*cldmax(i,k))
                        nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*cldmax(i,k-1)+ &
                            (rho(i,k)*dz(i,k)* &
                            (nprc(k)*cldm(i,k)+(nsubr(k-1)-npracs(k-1)-nnuccr(k-1)+nragg(k-1))*cldmax(i,k))))&
                            /(dum1*rho(i,k)*cldmax(i,k))
                    end if

                    !-------------------------
                    ! Autoconversion of cloud ice to snow
                    ! similar to Ferrier (1994)
                    if (t(i,k).le.273.15_r8.and.qiic(i,k).ge.qsmall) then
                        ! note: assumes autoconversion timescale of 180 sec
                        ! ! autoconversion size threshold for cloud ice to snow (m)	Dcs = 200.e-6_r8
                        nprci(k) = n0i(k)/(lami(k)*180._r8)*exp(-lami(k)*dcs)
                        prci(k) = pi*rhoi*n0i(k)/(6._r8*180._r8)* &
                            (dcs**3/lami(k)+3._r8*dcs**2/lami(k)**2+ &
                            6._r8*dcs/lami(k)**3+6._r8/lami(k)**4)*exp(-lami(k)*dcs)          
                    else
                        prci(k)=0._r8
                        nprci(k)=0._r8
                    end if
                    ! add autoconversion to flux from level above to get provisional snow mixing ratio
                    ! and number concentration (qniic and nsic)
                    dum=(asn(i,k)*dcs**bs)
                    dum1=(asn(i,k)*dcs**bs)
                    if (k.eq.1) then
                        qniic(i,k)=0.5_r8*prci(k)*cldm(i,k)*dz(i,k)/cldmax(i,k)/dum
                        nsic(i,k)=0.5_r8*nprci(k)*cldm(i,k)*dz(i,k)/cldmax(i,k)/dum
                    else
                        if (qniic(i,k-1).ge.qsmall) then
                            dum=ums(k-1)
                            dum1=uns(k-1)
                        end if
                        qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*cldmax(i,k-1)+ &
                            (rho(i,k)*dz(i,k)*(prci(k)*cldm(i,k)+(prai(k-1)+ &
                            psacws(k-1)+prci(k-1)+bergs(k-1))*cldmax(i,k))))&
                            /(dum*rho(i,k)*cldmax(i,k))
                        nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*cldmax(i,k-1)+ &
                            (rho(i,k)*dz(i,k)*(nprci(k)*cldm(i,k)+(nsubs(k-1)+ &
                            nsagg(k-1)+nnuccr(k-1))*cldmax(i,k))))&
                            /(dum1*rho(i,k)*cldmax(i,k))
                    end if
                    ! if precip mix ratio is zero so should number concentration
                    if (qniic(i,k).lt.qsmall) then
                        qniic(i,k)=0._r8
                        nsic(i,k)=0._r8
                    end if
                    if (qric(i,k).lt.qsmall) then
                        qric(i,k)=0._r8
                        nric(i,k)=0._r8
                    end if
                    ! make sure number concentration is a positive number to avoid 
                    ! taking root of negative later
                    nric(i,k)=max(nric(i,k),0._r8)
                    nsic(i,k)=max(nsic(i,k),0._r8)
                    !.-------------------
                    ! get size distribution parameters for precip
                    !-----
                    ! rain
                    if (qric(i,k).ge.qsmall) then
                        lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
                        n0r(k) = nric(i,k)*lamr(k)
                        ! check for slope
                        lammax = 1._r8/20.e-6_r8
                        lammin = 1._r8/500.e-6_r8
                        ! adjust vars
                        if (lamr(k).lt.lammin) then
                            lamr(k) = lammin
                            n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                            nric(i,k) = n0r(k)/lamr(k)
                        else if (lamr(k).gt.lammax) then
                            lamr(k) = lammax
                            n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                            nric(i,k) = n0r(k)/lamr(k)
                        end if
                        ! provisional rain number and mass weighted mean fallspeed (m/s)
                        unr(k) = min(arn(i,k)*gamma(1._r8+br)/lamr(k)**br,9.1_r8)
                        umr(k) = min(arn(i,k)*gamma(4._r8+br)/(6._r8*lamr(k)**br),9.1_r8)
                    else
                        lamr(k) = 0._r8
                        n0r(k) = 0._r8
                        umr(k) = 0._r8
                        unr(k) = 0._r8
                    end if
                    !-----
                    ! snow
                    ! snow mass-diameter relationship ::cs = rhosn*pi/6._r8 ::ds = 3._r8
                    if (qniic(i,k).ge.qsmall) then
                        lams(k) = (gamma(1._r8+ds)*cs*nsic(i,k)/ &
                            qniic(i,k))**(1._r8/ds)
                        n0s(k) = nsic(i,k)*lams(k)
                        ! check for slope
                        lammax = 1._r8/10.e-6_r8
                        lammin = 1._r8/2000.e-6_r8
                        ! adjust vars
                        if (lams(k).lt.lammin) then
                            lams(k) = lammin
                            n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*gamma(1._r8+ds))
                            nsic(i,k) = n0s(k)/lams(k)
                        else if (lams(k).gt.lammax) then
                            lams(k) = lammax
                            n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*gamma(1._r8+ds))
                            nsic(i,k) = n0s(k)/lams(k)
                        end if
                        ! provisional snow number and mass weighted mean fallspeed (m/s)
                        ums(k) = min(asn(i,k)*gamma(4._r8+bs)/(6._r8*lams(k)**bs),1.2_r8)
                        uns(k) = min(asn(i,k)*gamma(1._r8+bs)/lams(k)**bs,1.2_r8)
                    else
                        lams(k) = 0._r8
                        n0s(k) = 0._r8
                        ums(k) = 0._r8
                        uns(k) = 0._r8
                    end if
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! heterogeneous freezing of cloud water
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    if (qcic(i,k).ge.qsmall .and. t(i,k).lt.269.15_r8) then
                        !  immersion freezing (Bigg, 1953)
                        !  mnuccc(pver) ! mixing ratio tendency due to freezing of cloud water
                        !  nnuccc(pver) ! number conc tendency due to freezing of cloud water
                        !  cdist1(k) = ncic(i,k)/gamma(pgam(k)+1._r8) 
                        !  qcvar=1
                        ! immersion freezing parameters, bigg 1953  bimm = 100._r8  aimm = 0.66_r8
                        mnuccc(k) = gamma(qcvar+2._r8)/(gamma(qcvar)*qcvar**2)* &
                            pi*pi/36._r8*rhow* &
                            cdist1(k)*gamma(7._r8+pgam(k))* &
                            bimm*exp(aimm*(273.15_r8-t(i,k)))/ &
                            lamc(k)**3/lamc(k)**3
                        nnuccc(k) = gamma(qcvar+1._r8)/(gamma(qcvar)*qcvar)* &
                            pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
                            *bimm* &
                            exp(aimm*(273.15_r8-t(i,k)))/lamc(k)**3
                        ! hm add 11/17/06
                        ! make sure number of droplets frozen does not exceed available ice nuclei concentration
                        ! this prevents 'runaway' droplet freezing
                        ! nnuccd--- activation
                        if (nnuccc(k).gt.nnuccd(k)/cldm(i,k)) then
                            dum=(nnuccd(k)/cldm(i,k))/nnuccc(k)
                            ! scale mixing ratio of droplet freezing with limit
                            mnuccc(k)=mnuccc(k)*dum
                            nnuccc(k)=nnuccd(k)/cldm(i,k)
                        end if
                    else
                        mnuccc(k)=0._r8
                        nnuccc(k)=0._r8
                    end if
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! snow self-aggregation from passarelli, 1978, used by reisner, 1998
                    ! this is hard-wired for bs = 0.4 for now
                    ! ignore self-collection of cloud ice
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    if (qniic(i,k).ge.qsmall .and. t(i,k).le.273.15_r8) then
                        nsagg(k) = -1108._r8*asn(i,k)*Eii* &
                            pi**((1._r8-bs)/3._r8)*rhosn**((-2._r8-bs)/3._r8)*rho(i,k)** &
                            ((2._r8+bs)/3._r8)*qniic(i,k)**((2._r8+bs)/3._r8)* &
                            (nsic(i,k)*rho(i,k))**((4._r8-bs)/3._r8)/ &
                            (4._r8*720._r8*rho(i,k))
                    else
                        nsagg(k)=0._r8
                    end if
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! accretion of cloud droplets onto snow/graupel
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! here use continuous collection equation with simple gravitational collection kernel
                    ! ignore collisions between droplets/cloud ice
                    ! since minimum size ice particle for accretion is 50 - 150 micron
                    ! ignore collision of snow with droplets above freezing
                    !   if (qniic(i,k).ge.qsmall .and. t(i,k).le.273.15_r8 .and. &
                    !      qcic(i,k).ge.qsmall) then  !! so odd -sxj
                    if ( qniic(i,k)==qsmall) then
                        write(*,*) "qniic(i,k)==qsmall"
                        write(*,*) "cldwat2m.F90"  !!!
                    endif
                    if (qniic(i,k).gt.qsmall .and. t(i,k).le.273.15_r8 .and. &
                        qcic(i,k).ge.qsmall) then  !! so odd-sxj  @XXXXXXXXX   
                        ! put in size dependent collection efficiency mean diameter of snow is area-weighted, since
                        ! accretion is function of crystal geometric area
                        ! collection efficiency is from stoke's law (Thompson et al. 2004)
                        dc0 = (pgam(k)+1._r8)/lamc(k)
                        ds0 = 1._r8/lams(k)
                        dum = dc0*dc0*uns(k)*rhow/(9._r8*mu(i,k)*ds0)
                        eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))
                        eci = max(eci,0._r8)
                        eci = min(eci,1._r8)
                        ! no impact of sub-grid distribution of qc since psacws
                        ! is linear in qc
                        psacws(k) = pi/4._r8*asn(i,k)*qcic(i,k)*rho(i,k)* &
                            n0s(k)*Eci*gamma(bs+3._r8)/ &
                            lams(k)**(bs+3._r8)	
                        npsacws(k) = pi/4._r8*asn(i,k)*ncic(i,k)*rho(i,k)* &
                            n0s(k)*Eci*gamma(bs+3._r8)/ &
                            lams(k)**(bs+3._r8)
                    else
                        psacws(k)=0._r8
                        npsacws(k)=0._r8
                    end if
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! accretion of rain water by snow
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! formula from ikawa and saito, 1991, used by reisner et al., 1998
                    if (qric(i,k).ge.1.e-8_r8 .and. qniic(i,k).ge.1.e-8_r8 .and. & 
                        t(i,k).le.273.15_r8) then
                        pracs(k) = pi*pi*ecr*(((1.2_r8*umr(k)-0.95_r8*ums(k))**2+ &
                            0.08_r8*ums(k)*umr(k))**0.5_r8*rhow*rho(i,k)* &
                            n0r(k)*n0s(k)* &
                            (5._r8/(lamr(k)**6*lams(k))+ &
                            2._r8/(lamr(k)**5*lams(k)**2)+ &
                            0.5_r8/(lamr(k)**4*lams(k)**3)))
                        npracs(k) = pi/2._r8*rho(i,k)*ecr*(1.7_r8*(unr(k)-uns(k))**2+ &
                            0.3_r8*unr(k)*uns(k))**0.5_r8*n0r(k)*n0s(k)* &
                            (1._r8/(lamr(k)**3*lams(k))+ &
                            1._r8/(lamr(k)**2*lams(k)**2)+ &
                            1._r8/(lamr(k)*lams(k)**3))
                    else
                        pracs(k)=0._r8
                        npracs(k)=0._r8
                    end if
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! heterogeneous freezing of rain drops
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! follows from Bigg (1953)
                    if (t(i,k).lt.269.15_r8 .and. qric(i,k).ge.qsmall) then
                        mnuccr(k) = 20._r8*pi*pi*rhow*nric(i,k)*bimm* &
                            exp(aimm*(273.15_r8-t(i,k)))/lamr(k)**3 &
                            /lamr(k)**3
                        nnuccr(k) = pi*nric(i,k)*bimm* &
                            exp(aimm*(273.15_r8-t(i,k)))/lamr(k)**3
                    else
                        mnuccr(k)=0._r8
                        nnuccr(k)=0._r8
                    end if
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! accretion of cloud liquid water by rain
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! formula from Khrouditnov and Kogan (2000)
                    ! gravitational collection kernel, droplet fall speed neglected
                    if (qric(i,k).ge.qsmall .and. qcic(i,k).ge.qsmall) then
                        ! include sub-grid distribution of cloud water
                        pra(k) = gamma(qcvar+1.15_r8)/(gamma(qcvar)*qcvar**1.15_r8)* &
                            67._r8*(qcic(i,k)*qric(i,k))**1.15_r8
                        npra(k) = pra(k)/(qcic(i,k)/ncic(i,k))
                    else
                        pra(k)=0._r8
                        npra(k)=0._r8
                    end if
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! Self-collection of rain drops
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! from Beheng(1994)
                    if (qric(i,k)==qsmall) write(*,*) "qric(i,k)==qsmall,cldwat2.f90"
                    if (qric(i,k).gt.qsmall) then   !! change ge to gt  sxj!!!!!!!!!
                        nragg(k) = -8._r8*nric(i,k)*qric(i,k)*rho(i,k)
                    else
                        nragg(k)=0._r8
                    end if
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! Accretion of cloud ice by snow
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! For this calculation, it is assumed that the Vs >> Vi
                    ! and Ds >> Di for continuous collection
                    if (qniic(i,k).ge.qsmall.and.qiic(i,k).ge.qsmall &
                        .and.t(i,k).le.273.15_r8) then
                        prai(k) = pi/4._r8*asn(i,k)*qiic(i,k)*rho(i,k)* &
                            n0s(k)*Eii*gamma(bs+3._r8)/ &
                            lams(k)**(bs+3._r8)	
                        nprai(k) = pi/4._r8*asn(i,k)*niic(i,k)* &
                            rho(i,k)*n0s(k)*Eii*gamma(bs+3._r8)/ &
                            lams(k)**(bs+3._r8)
                    else
                        prai(k)=0._r8
                        nprai(k)=0._r8
                    end if
                    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! calculate evaporation/sublimation of rain and snow
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
                    ! in-cloud condensation/deposition of rain and snow is neglected
                    ! except for transfer of cloud water to snow through bergeron process
                    ! initialize evap/sub tendncies  ???? why only mixing ratio 
                    ! pre(pver) ! rain mixing rat tendency due to evaporation          
                    pre(k)=0._r8
                    prds(k)=0._r8
                    ! evaporation of rain
                    ! only calculate if there is some precip fraction > cloud fraction
                    if (qcic(i,k)+qiic(i,k).lt.1.e-6_r8.or.cldmax(i,k).gt.cldm(i,k)) then
                        ! set temporary cloud fraction to zero if cloud water + ice is very small
                        ! this will ensure that evaporation/sublimation of precip occurs over
                        ! entire grid cell, since min cloud fraction is specified otherwise
                        if (qcic(i,k)+qiic(i,k).lt.1.e-6_r8) then
                            dum=0._r8
                        else
                            dum=cldm(i,k)
                        end if
                        esn=estblf(t(i,k))
                        qsn=min(epsqs*esn/(p(i,k)-(1._r8-epsqs)*esn),1._r8)
                        ! recalculate saturation vapor pressure for liquid and ice
                        esl(i,k)=polysvp(t(i,k),0)
                        esi(i,k)=polysvp(t(i,k),1)
                        ! calculate q for out-of-cloud region
                        qclr=(q(i,k)-dum*qsn)/(1._r8-dum)
                        if (qric(i,k).ge.qsmall) then
                            qvs=0.622_r8*esl(i,k)/(p(i,k)-esl(i,k))
                            dqsdt = xxlv*qvs/(rv*t(i,k)**2)
                            ! dqsdt--change of sat vapor mixing ratio with temperature
                            ab = 1._r8+dqsdt*xxlv/cpp
                            epsr = 2._r8*pi*n0r(k)*rho(i,k)*Dv(i,k)* &
                                (f1r/(lamr(k)*lamr(k))+ &
                                f2r*(arn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                                sc(i,k)**(1._r8/3._r8)*gamma(5._r8/2._r8+br/2._r8)/ &
                                (lamr(k)**(5._r8/2._r8+br/2._r8)))
                            pre(k) = epsr*(qclr-qvs)/ab
                            ! only evaporate in out-of-cloud region
                            ! and distribute across cldmax
                            pre(k)=min(pre(k)*(cldmax(i,k)-dum),0._r8)
                            pre(k)=pre(k)/cldmax(i,k)
                        end if
                        ! sublimation of snow
                        if (qniic(i,k).ge.qsmall) then
                            qvi=0.622_r8*esi(i,k)/(p(i,k)-esi(i,k))
                            dqsidt =  xxls*qvi/(rv*t(i,k)**2)
                            abi = 1._r8+dqsidt*xxls/cpp
                            epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
                                (f1s/(lams(k)*lams(k))+ &
                                f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                                sc(i,k)**(1._r8/3._r8)*gamma(5._r8/2._r8+bs/2._r8)/ &
                                (lams(k)**(5._r8/2._r8+bs/2._r8)))
                            prds(k) = epss*(qclr-qvi)/abi	
                            ! only sublimate in out-of-cloud region and distribute over cldmax!!!!!!!
                            prds(k)=min(prds(k)*(cldmax(i,k)-dum),0._r8)
                            prds(k)=prds(k)/cldmax(i,k)
                        end if
                        ! hm add 2/2/07, make sure RH not pushed above 100%
                        ! get updated RH at end of time step based on cloud water/ice condensation/evap
                        qtmp=q(i,k)-(cmel(i,k)+cmei(i,k)+(pre(k)+prds(k))*cldmax(i,k))*deltat   
                        ttmp=t(i,k)+((cmel(i,k)+pre(k)*cldmax(i,k))*xxlv+ &
                            (cmei(i,k)+prds(k))*cldmax(i,k)*xxls)*deltat/cpp
                        esn=estblf(ttmp)
                        qsn=min(epsqs*esn/(p(i,k)-(1._r8-epsqs)*esn),1._r8)
                        ! modify precip evaporation rate if q > qsat
                        if (qtmp.gt.qsn) then !!!!----
                            if (pre(k)+prds(k).lt.-1.e-20) then       !!!! ???pre prds minus denote evaporation
                                dum1=pre(k)/(pre(k)+prds(k))
                                ! recalculate q and t after cloud water cond but without precip evap
                                qtmp=q(i,k)-(cmel(i,k)+cmei(i,k))*deltat
                                ttmp=t(i,k)+(cmel(i,k)*xxlv+ &
                                    cmei(i,k)*xxls)*deltat/cpp
                                esn=estblf(ttmp)
                                qsn=min(epsqs*esn/(p(i,k)-(1._r8-epsqs)*esn),1._r8)
                                dum=(qtmp-qsn)/(1._r8+(xxlv*dum1+xxls*(1._r8-dum1))**2*qsn &
                                    /(cpp*rv*ttmp**2))
                                dum=min(dum,0._r8)
                                ! modify rates if needed, divide by cldmax to get local (in-precip) value
                                pre(k)=dum*dum1/deltat/cldmax(i,k)
                                prds(k)=dum*(1._r8-dum1)/deltat/cldmax(i,k)
                            end if
                        end if
                    end if

                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! bergeron process - evaporation of droplets and deposition onto snow
                    ! bergeron process for snow is neglected for now.............
                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    !	if (qniic(i,k).ge.qsmall.and.qcic(i,k).ge.qsmall.and.t(i,k).lt.273.15_r8) then
                    !        qvs=0.622_r8*esl(i,k)/(p(i,k)-esl(i,k))
                    !        qvi=0.622_r8*esi(i,k)/(p(i,k)-esi(i,k))
                    !        dqsidt =  xxls*qvi/(rv*t(i,k)**2)
                    !        abi = 1._r8+dqsidt*xxls/cpp
                    !        epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
                    !                   (f1s/(lams(k)*lams(k))+ &
                    !                    f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                    !                    sc(i,k)**(1._r8/3._r8)*gamma(5._r8/2._r8+bs/2._r8)/ &
                    !               (lams(k)**(5._r8/2._r8+bs/2._r8)))
                    !	bergs(k)=epss*(qvs-qvi)/abi
                    !++ag
                    bergs(k)=0._r8
                    !--ag
                    !	else
                    !	bergs(k)=0._r8
                    !	end if

                    ! sensitivity - no rain/snow evaporation
                    !	pre(k)=0._r8
                    !	prds(k)=0._r8

                    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    ! conservation to ensure no negative values of cloud water/precipitation
                    ! in case microphysical process rates are large

                    ! make sure and use end-of-time step values for cloud water, ice, due
                    ! condensation/deposition

                    ! note: for check on conservation, processes are multiplied by omsm
                    ! to prevent problems due to round off error

                    ! since activation/nucleation processes are fast, need to take into account
                    ! factor mtime = mixing timescale in cloud / model time step
                    ! mixing time can be interpreted as cloud depth divided by sub-grid vertical velocity
                    ! for now mixing timescale is assumed to be 20 min
                    ! could possibly be estimated better from model variables

                    mtime=deltat/1200._r8     !!! time step 1200!!!!!!!!!!!!
                    qce=(qc(i,k)+(cmel(i,k)-berg(i,k))*deltat)
                    nce=(nc(i,k)+npccn(k)*deltat*mtime)/cldm(i,k)
                    qie=(qi(i,k)+(cmei(i,k)+berg(i,k))*deltat)
                    nie=(ni(i,k)+nnuccd(k)*deltat*mtime)/cldm(i,k)
                    ! conservation of qc
                    dum = (prc(k)+pra(k)+mnuccc(k)+ &
                        psacws(k)+bergs(k))*cldm(i,k)*deltat
                    if (dum.gt.qce) then
                        ratio = qce/deltat/cldm(i,k)/(prc(k)+pra(k)+ &
                            mnuccc(k)+psacws(k)+bergs(k))*omsm
                        prc(k) = prc(k)*ratio
                        pra(k) = pra(k)*ratio
                        mnuccc(k) = mnuccc(k)*ratio
                        psacws(k) = psacws(k)*ratio	
                        bergs(k) = bergs(k)*ratio	
                    end if
                    ! conservation of nc
                    dum = (nprc1(k)+npra(k)+nnuccc(k)+ &
                        npsacws(k)-nsubc(k))*deltat
                    if (dum.gt.nce) then
                        ratio = nce/deltat/(nprc1(k)+npra(k)+nnuccc(k)+&
                            npsacws(k)-nsubc(k))*omsm
                        nprc1(k) = nprc1(k)*ratio
                        npra(k) = npra(k)*ratio
                        nnuccc(k) = nnuccc(k)*ratio
                        npsacws(k) = npsacws(k)*ratio	
                        nsubc(k)=nsubc(k)*ratio
                    end if
                    ! conservation of qi
                    dum = (-mnuccc(k)+prci(k)+ &
                        prai(k))*cldm(i,k)*deltat
                    if ((dum.gt.qie).and.(abs(prci(k)+prai(k))>1.e-150)) then  !ljli according to Liu Li
                        ratio = (qie/deltat/cldm(i,k)+mnuccc(k))/(prci(k)+prai(k))*omsm
                        prci(k) = prci(k)*ratio
                        prai(k) = prai(k)*ratio
                    end if
                    ! conservation of ni
                    dum = (nprci(k)+ &
                        nprai(k)-nsubi(k))*deltat
                    if (dum.gt.nie) then
                        ratio = (nie/deltat)/ &
                            (nprci(k)+nprai(k)-nsubi(k))*omsm
                        nprci(k) = nprci(k)*ratio
                        nprai(k) = nprai(k)*ratio
                        nsubi(k) = nsubi(k)*ratio
                    end if
                    ! for preciptiation conservation, use logic that vertical integral 
                    ! of tendency from current level to top of model (i.e., qrtot) cannot be negative
                    ! conservation of rain mixing rat
                    !   pre(pver) ------ rain mixing rat tendency due to evaporation
                    !   mnuccr(pver) - mixing rat tendency due to freezing of rain
                    !   pracs(pver)-- mixing rat tendency due to collection of rain	by snow
                    if (((prc(k)+pra(k))*cldm(i,k)+(-mnuccr(k)+pre(k)-pracs(k))*&
                        cldmax(i,k))*dz(i,k)*rho(i,k)+qrtot.lt.0._r8) then
                        if (-pre(k)+pracs(k)+mnuccr(k).ge.qsmall) then	     
                            ratio = (qrtot/(dz(i,k)*rho(i,k))+(prc(k)+pra(k))*cldm(i,k))/&
                                ((-pre(k)+pracs(k)+mnuccr(k))*cldmax(i,k))*omsm 
                            pre(k) = pre(k)*ratio
                            pracs(k) = pracs(k)*ratio
                            mnuccr(k) = mnuccr(k)*ratio
                            !		      else
                            !			    write(*,*) "clwat2.f90-1823"
                            !				write(*,*) "-pre(k)+pracs(k)+mnuccr(k).ge.qsmall"
                        end if
                    end if
                    ! conservation of nr
                    ! for now neglect evaporation of nr
                    ! calculate evaporation of nr
                    !       if (pre(k).lt.0._r8.and.qric(i,k).ge.qsmall) then
                    !       nsubr(k)=pre(k)/qric(i,k)*nric(i,k)        
                    !       else
                    nsubr(k)=0._r8  !
                    !       end if
                    if ((nprc(k)*cldm(i,k)+(-nnuccr(k)+nsubr(k)-npracs(k)&
                        +nragg(k))*cldmax(i,k))*dz(i,k)*rho(i,k)+nrtot.lt.0._r8) then
                        if (-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k).ge.qsmall) then
                            ratio = (nrtot/(dz(i,k)*rho(i,k))+nprc(k)*cldm(i,k))/&
                                ((-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k))*cldmax(i,k))*omsm
                            nsubr(k) = nsubr(k)*ratio
                            npracs(k) = npracs(k)*ratio
                            nnuccr(k) = nnuccr(k)*ratio
                            nragg(k) = nragg(k)*ratio
                            !			  else
                            !			    write(*,*) "-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k).<.qsmall" !sxj
                        end if
                    end if
                    ! conservation of snow mix ratio
                    if (((bergs(k)+psacws(k)+prai(k)+prci(k))*cldm(i,k)+(pracs(k)+&
                        mnuccr(k)+prds(k))*cldmax(i,k))*dz(i,k)*rho(i,k)+qstot.lt.0._r8) then
                        if (-prds(k).ge.qsmall) then
                            ratio = (qstot/(dz(i,k)*rho(i,k))+(bergs(k)+psacws(k)+prai(k)+prci(k)) &
                                *cldm(i,k)+(pracs(k)+mnuccr(k))*cldmax(i,k))/(-prds(k)*cldmax(i,k))*omsm
                            prds(k) = prds(k)*ratio
                            !		   else
                            !			    write(*,*) "-prds(k).<.qsmall"  !sxj
                        end if
                    end if
                    ! conservation of ns
                    ! calculate loss of number due to sublimation
                    ! for now neglect sublimation of ns
                    !       if (prds(k).lt.0._r8.and.qniic(i,k).ge.qsmall) then
                    !       nsubs(k)=prds(k)/qniic(i,k)*nsic(i,k)     !!!! 
                    !       else
                    nsubs(k)=0._r8
                    !       end if
                    if ((nprci(k)*cldm(i,k)+(nnuccr(k)+nsubs(k)+nsagg(k))*cldmax(i,k))*&
                        dz(i,k)*rho(i,k)+nstot.lt.0._r8) then
                        if (-nsubs(k)-nsagg(k).ge.qsmall) then
                            ratio = (nstot/(dz(i,k)*rho(i,k))+nprci(k)*cldm(i,k)+&
                                nnuccr(k)*cldmax(i,k))/((-nsubs(k)-nsagg(k))*cldmax(i,k))*omsm
                            nsubs(k) = nsubs(k)*ratio
                            nsagg(k) = nsagg(k)*ratio
                            !			  else
                            !			    write(*,*) "clwat2.f90-1877"
                            !				write(*,*) "-nsubs(k)-nsagg(k).<.qsmall"
                        end if
                    end if
                    ! get tendencies due to microphysical conversion processes
                    ! note: tendencies are multiplied by appropaiate cloud/precip 
                    ! fraction to get grid-scale values
                    ! note: cmei,cmel are already grid-average values
                    qvlat(i,k) = qvlat(i,k)- &
                        (pre(k)+prds(k))*cldmax(i,k)-cmel(i,k)-cmei(i,k)
                    tlat(i,k) = tlat(i,k)+((pre(k)*cldmax(i,k)+cmel(i,k)) &
                        *xxlv+(prds(k)*cldmax(i,k)+cmei(i,k))*xxls+ &
                        ((bergs(k)+psacws(k)+mnuccc(k))*cldm(i,k)+(mnuccr(k)+ &
                        pracs(k))*cldmax(i,k)+berg(i,k))*xlf)	   
                    qctend(i,k) = qctend(i,k)+ &
                        (-pra(k)-prc(k)-mnuccc(k)- &
                        psacws(k)-bergs(k))*cldm(i,k)+cmel(i,k)-berg(i,k)
                    qitend(i,k) = qitend(i,k)+ &
                        (mnuccc(k)-prci(k)- &
                        prai(k))*cldm(i,k)+cmei(i,k)+berg(i,k)
                    qrtend(i,k) = qrtend(i,k)+ &
                        (pra(k)+prc(k))*cldm(i,k)+(pre(k)-pracs(k)- &
                        mnuccr(k))*cldmax(i,k)
                    qnitend(i,k) = qnitend(i,k)+ &
                        (prai(k)+psacws(k)+prci(k)+bergs(k))*cldm(i,k)+(prds(k)+ &
                        pracs(k)+mnuccr(k))*cldmax(i,k)

                    ! assign variables for trop_mozart, these are grid-average
                    ! evaporation/sublimation is stored here as positive term
                    evapsnow1(i,k) = evapsnow1(i,k)-prds(k)*cldmax(i,k)
                    nevapr1(i,k) = nevapr1(i,k)-pre(k)*cldmax(i,k)
                    !++ag 2.3 change to make sure prain is positive: do not remove snow from
                    ! prain used for wet deposition
                    !	prain1(i,k) = prain1(i,k)+(pra(k)+prc(k))*cldm(i,k)+(-pracs(k)- &
                    !                  mnuccr(k))*cldmax(i,k)
                    prain1(i,k) = prain1(i,k)+(pra(k)+prc(k))*cldm(i,k)
                    !--ag
                    prodsnow1(i,k) = prodsnow1(i,k)+(prai(k)+psacws(k)+prci(k)+bergs(k))*cldm(i,k)+(&
                        pracs(k)+mnuccr(k))*cldmax(i,k)
                    ! multiply activation/nucleation by mtime to account for fast timescale
                    nctend(i,k) = nctend(i,k)+ npccn(k)*mtime+&
                        (-nnuccc(k)-npsacws(k)+nsubc(k) &
                        -npra(k)-nprc1(k))*cldm(i,k)      
                    nitend(i,k) = nitend(i,k)+ nnuccd(k)*mtime+&
                        (nsubi(k)-nprci(k)- &
                        nprai(k))*cldm(i,k)
                    nstend(i,k) = nstend(i,k)+(nsubs(k)+ &
                        nsagg(k)+nnuccr(k))*cldmax(i,k)+nprci(k)*cldm(i,k)
                    nrtend(i,k) = nrtend(i,k)+ &
                        nprc(k)*cldm(i,k)+(nsubr(k)-npracs(k)-nnuccr(k) &
                        +nragg(k))*cldmax(i,k)
                    ! make sure that nc and ni at advanced time step do not exceed
                    ! maximum (existing N + source terms*dt), which is possible due to
                    ! fast nucleation timescale

                    ! ncmax/nimax ~ activation process
                    if (nctend(i,k).gt.0._r8.and.nc(i,k)+nctend(i,k)*deltat.gt.ncmax) then
                        nctend(i,k)=max(0._r8,(ncmax-nc(i,k))/deltat)
                    end if
                    if (nitend(i,k).gt.0._r8.and.ni(i,k)+nitend(i,k)*deltat.gt.nimax) then
                        nitend(i,k)=max(0._r8,(nimax-ni(i,k))/deltat)
                    end if
                    ! get final values for precipitation q and N, based on
                    ! flux of precip from above, source/sink term, and terminal fallspeed
                    ! see eq. 15-16 in Morrison and Gettelman, 2007, J. Climate
                    ! rain

                    if (qric(i,k).ge.qsmall) then
                        if (k.eq.1) then
                            qric(i,k)=qrtend(i,k)*dz(i,k)/cldmax(i,k)/umr(k)
                            nric(i,k)=nrtend(i,k)*dz(i,k)/cldmax(i,k)/unr(k)
                        else
                            qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*cldmax(i,k-1)+ &
                                (rho(i,k)*dz(i,k)*qrtend(i,k)))/(umr(k)*rho(i,k)*cldmax(i,k))
                            nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*cldmax(i,k-1)+ &
                                (rho(i,k)*dz(i,k)*nrtend(i,k)))/(unr(k)*rho(i,k)*cldmax(i,k))
                        end if
                    else
                        qric(i,k)=0._r8
                        nric(i,k)=0._r8
                    end if
                    if (qric(i,k)<0.) then            !!sxj 
                        !		       write(*,*) "cldwat2.f90 1961","qric<0",qric(i,k)
                        qric(i,k)=0._r8
                        nric(i,k)=0._r8		   
                    endif
                    ! snow
                    if (qniic(i,k).ge.qsmall) then
                        if (k.eq.1) then
                            qniic(i,k)=qnitend(i,k)*dz(i,k)/cldmax(i,k)/ums(k)	
                            nsic(i,k)=nstend(i,k)*dz(i,k)/cldmax(i,k)/uns(k)
                        else
                            qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*cldmax(i,k-1)+ &
                                (rho(i,k)*dz(i,k)*qnitend(i,k)))/(ums(k)*rho(i,k)*cldmax(i,k))
                            nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*cldmax(i,k-1)+ &
                                (rho(i,k)*dz(i,k)*nstend(i,k)))/(uns(k)*rho(i,k)*cldmax(i,k))
                        end if
                    else
                        qniic(i,k)=0._r8
                        nsic(i,k)=0._r8
                    end if
                    if (qniic(i,k)<0) then   ! sxj----
                        !		  !    write(*,*) "cldwat2.f90 1980","qniic <0",qniic(i,k)
                        qniic(i,k)=0._r8
                        nsic(i,k)=0._r8			   
                    endif

                    ! calculate precipitation flux at surface
                    ! divide by density of water to get units of m/s
                    prect(i) = prect(i)+(qrtend(i,k)*dz(i,k)*rho(i,k)+&
                        qnitend(i,k)*dz(i,k)*rho(i,k))/rhow  ! ##### rain and snow   rhow 
                    preci(i) = preci(i)+qnitend(i,k)*dz(i,k)*rho(i,k)/rhow
                    ! vertically-integrated precip source/sink terms (note: grid-averaged)
                    qrtot = max(qrtot+qrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
                    qstot = max(qstot+qnitend(i,k)*dz(i,k)*rho(i,k),0._r8)
                    nrtot = max(nrtot+nrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
                    nstot = max(nstot+nstend(i,k)*dz(i,k)*rho(i,k),0._r8)
                    ! calculate melting and freezing of precip
                    ! melt snow at +2 C

                    if (t(i,k)+tlat(i,k)/cp*deltat > 275.15_r8) then
                        if (qstot > 0._r8) then
                            ! make sure melting snow doesn't reduce temperature below threshold
                            dum = -xlf/cp*qstot/(dz(i,k)*rho(i,k))
                            if (t(i,k)+tlat(i,k)/cp*deltat+dum.lt.275.15_r8) then
                                dum = (t(i,k)+tlat(i,k)/cp*deltat-275.15_r8)*cp/xlf
                                dum = dum/(xlf/cp*qstot/(dz(i,k)*rho(i,k)))
                                dum = max(0._r8,dum)
                                dum = min(1._r8,dum)
                            else
                                dum = 1._r8
                            end if
                            qric(i,k)=qric(i,k)+dum*qniic(i,k)
                            nric(i,k)=nric(i,k)+dum*nsic(i,k)
                            qniic(i,k)=(1._r8-dum)*qniic(i,k)
                            nsic(i,k)=(1._r8-dum)*nsic(i,k)
                            tlat(i,k)=tlat(i,k)-xlf*dum*qstot/(dz(i,k)*rho(i,k))
                            qrtot=qrtot+dum*qstot   ! 
                            nrtot=nrtot+dum*nstot
                            qstot=(1._r8-dum)*qstot
                            nstot=(1._r8-dum)*nstot
                            preci(i)=(1._r8-dum)*preci(i)
                        end if
                    end if


                    ! freeze rain homogeneously at -40 C
                    if (t(i,k)+tlat(i,k)/cp*deltat < 233.15_r8) then
                        if (qrtot > 0._r8) then
                            ! make sure freezing rain doesn't increase temperature above threshold
                            dum = xlf/cp*qrtot/(dz(i,k)*rho(i,k))
                            if (t(i,k)+tlat(i,k)/cp*deltat+dum.gt.233.15_r8) then
                                dum = -(t(i,k)+tlat(i,k)/cp*deltat-233.15_r8)*cp/xlf
                                dum = dum/(xlf/cp*qrtot/(dz(i,k)*rho(i,k)))
                                dum = max(0._r8,dum)
                                dum = min(1._r8,dum)
                            else
                                dum = 1._r8
                            end if
                            qniic(i,k)=qniic(i,k)+dum*qric(i,k)
                            nsic(i,k)=nsic(i,k)+dum*nric(i,k)
                            qric(i,k)=(1._r8-dum)*qric(i,k)
                            nric(i,k)=(1._r8-dum)*nric(i,k)
                            tlat(i,k)=tlat(i,k)+xlf*dum*qrtot/(dz(i,k)*rho(i,k))
                            qstot=qstot+dum*qrtot
                            qrtot=(1._r8-dum)*qrtot
                            nstot=nstot+dum*nrtot
                            nrtot=(1._r8-dum)*nrtot
                            preci(i)=preci(i)+dum*(prect(i)-preci(i))
                        end if
                    end if
                    ! if rain/snow mix ratio is zero so should number concentration
                    if (qniic(i,k).lt.qsmall) then
                        qniic(i,k)=0._r8
                        nsic(i,k)=0._r8
                    end if
                    if (qric(i,k).lt.qsmall) then
                        qric(i,k)=0._r8
                        nric(i,k)=0._r8
                    end if
                    ! make sure number concentration is a positive number to avoid 
                    ! taking root of negative
                    nric(i,k)=max(nric(i,k),0._r8)
                    nsic(i,k)=max(nsic(i,k),0._r8)
                    !.......................................................................
                    ! get size distribution parameters for fallspeed calculations
                    !......................................................................
                    ! rain
                    if (qric(i,k).ge.qsmall) then
                        lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
                        n0r(k) = nric(i,k)*lamr(k)
                        ! check for slope
                        ! hm 4/5/07, change lammax and lammin for rain and snow
                        lammax = 1._r8/20.e-6_r8
                        lammin = 1._r8/500.e-6_r8
                        ! adjust vars
                        if (lamr(k).lt.lammin) then
                            lamr(k) = lammin
                            n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                            nric(i,k) = n0r(k)/lamr(k)
                        else if (lamr(k).gt.lammax) then
                            lamr(k) = lammax
                            n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                            nric(i,k) = n0r(k)/lamr(k)
                        end if
                        ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)
                        unr(k) = min(arn(i,k)*gamma(1._r8+br)/lamr(k)**br,9.1_r8)
                        umr(k) = min(arn(i,k)*gamma(4._r8+br)/(6._r8*lamr(k)**br),9.1_r8)
                    else
                        lamr(k) = 0._r8
                        n0r(k) = 0._r8
                        umr(k)=0._r8
                        unr(k)=0._r8
                    end if
                    ! snow
                    if (qniic(i,k).ge.qsmall) then
                        lams(k) = (gamma(1._r8+ds)*cs*nsic(i,k)/ &
                            qniic(i,k))**(1._r8/ds)
                        n0s(k) = nsic(i,k)*lams(k)
                        ! check for slope
                        lammax = 1._r8/10.e-6_r8
                        lammin = 1._r8/2000.e-6_r8
                        ! adjust vars
                        if (lams(k).lt.lammin) then
                            lams(k) = lammin
                            n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*gamma(1._r8+ds))
                            nsic(i,k) = n0s(k)/lams(k)
                        else if (lams(k).gt.lammax) then
                            lams(k) = lammax
                            n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*gamma(1._r8+ds))
                            nsic(i,k) = n0s(k)/lams(k)
                        end if
                        ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)
                        ums(k) = min(asn(i,k)*gamma(4._r8+bs)/(6._r8*lams(k)**bs),1.2_r8)
                        uns(k) = min(asn(i,k)*gamma(1._r8+bs)/lams(k)**bs,1.2_r8)
                    else
                        lams(k) = 0._r8
                        n0s(k) = 0._r8
                        ums(k) = 0._r8
                        uns(k) = 0._r8
                    end if
                    !------------
                    ! sum over sub-step for average process rates
                    ! convert rain/snow q and N for output to history, note, 
                    ! output is for in-precip (local) values
                    qrout(i,k)=qrout(i,k)+qric(i,k)
                    qsout(i,k)=qsout(i,k)+qniic(i,k)
                    nrout(i,k)=nrout(i,k)+nric(i,k)*rho(i,k)
                    nsout(i,k)=nsout(i,k)+nsic(i,k)*rho(i,k)
                    tlat1(i,k)=tlat1(i,k)+tlat(i,k)
                    qvlat1(i,k)=qvlat1(i,k)+qvlat(i,k)
                    qctend1(i,k)=qctend1(i,k)+qctend(i,k)
                    qitend1(i,k)=qitend1(i,k)+qitend(i,k)
                    nctend1(i,k)=nctend1(i,k)+nctend(i,k)
                    nitend1(i,k)=nitend1(i,k)+nitend(i,k)
                    t(i,k)=t(i,k)+tlat(i,k)*deltat/cpp
                    q(i,k)=q(i,k)+qvlat(i,k)*deltat
                    qc(i,k)=qc(i,k)+qctend(i,k)*deltat
                    qi(i,k)=qi(i,k)+qitend(i,k)*deltat
                    nc(i,k)=nc(i,k)+nctend(i,k)*deltat
                    ni(i,k)=ni(i,k)+nitend(i,k)*deltat
                end do ! k loop

                prect1(i)=prect1(i)+prect(i)
                preci1(i)=preci1(i)+preci(i)
            end do    ! it loop, sub-step

300         continue  ! continue if no cloud water
        end do    ! i loop

        !#ifdef SPMD
        !	    CALL mpibarrier (mpicom)                                       !! for debug-test   
        !       write(6,*) iam,"mmicro_pcond  MG done **********sxj*********"  !!  sxj-
        !#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        ! convert dt from sub-step back to full time step
	deltat=deltat*real(iter)

        !c.............................................................................

        do i=1,ncol

            ! skip all calculations if no cloud water
            if (ltrue(i).eq.0) then

                do k=1,pver
                    ! assign default values for effective radius
                    effc(i,k)=10._r8
                    effi(i,k)=25._r8
                    effc_fn(i,k)=10._r8
                    !effc(pcols,pver)    ! droplet effective radius (micron)
                    !effc_fn(pcols,pver) ! droplet effective radius, assuming nc = 1.e8 kg-1
                end do
                goto 500
            end if

            ! initialize nstep for sedimentation sub-steps
            nstep = 1

            ! divide precip rate by number of sub-steps to get average over time step

            prect(i)=prect1(i)/real(iter)
            preci(i)=preci1(i)/real(iter)

            do k=1,pver
                ! assign variables back to start-of-timestep values before updating after sub-steps 

                t(i,k)=t1(i,k)  !
                q(i,k)=q1(i,k)
                qc(i,k)=qc1(i,k)
                qi(i,k)=qi1(i,k)
                nc(i,k)=nc1(i,k)
                ni(i,k)=ni1(i,k)

                ! divide microphysical tendencies by number of sub-steps to get average over time step

                tlat(i,k)=tlat1(i,k)/real(iter)
                qvlat(i,k)=qvlat1(i,k)/real(iter)
                qctend(i,k)=qctend1(i,k)/real(iter)
                qitend(i,k)=qitend1(i,k)/real(iter)
                nctend(i,k)=nctend1(i,k)/real(iter)
                nitend(i,k)=nitend1(i,k)/real(iter)

                ! divide output precip q and N by number of sub-steps to get average over time step

                qrout(i,k)=qrout(i,k)/real(iter)
                qsout(i,k)=qsout(i,k)/real(iter)
                nrout(i,k)=nrout(i,k)/real(iter)
                nsout(i,k)=nsout(i,k)/real(iter)

                ! divide trop_mozart variables by number of sub-steps to get average over time step 

                nevapr(i,k) = nevapr1(i,k)/real(iter)
                evapsnow(i,k) = evapsnow1(i,k)/real(iter)
                prain(i,k) = prain1(i,k)/real(iter)
                prodsnow(i,k) = prodsnow1(i,k)/real(iter)
                cmeout(i,k) = cmeout1(i,k)/real(iter)

                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! calculate sedimentation for cloud water and ice
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

                ! update in-cloud cloud mixing ratio and number concentration 
                ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
                ! note: these are in-cloud values***, hence we divide by cloud fraction

                dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)/cldm(i,k)
                dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)/cldm(i,k)
                dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat)/cldm(i,k),0._r8)
                dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat)/cldm(i,k),0._r8)

                ! obtain new slope parameter to avoid possible singularity

                if (dumi(i,k).ge.qsmall) then
                    ! add upper limit to in-cloud number concentration to prevent numerical error
                    dumni(i,k)=min(dumni(i,k),dumi(i,k)*1.e20_r8)

                    lami(k) = (gamma(1._r8+di)*ci* &
                        dumni(i,k)/dumi(i,k))**(1._r8/di)
                    lammax = 1._r8/10.e-6_r8
                    lammin = 1._r8/(2._r8*dcs)
                    lami(k)=max(lami(k),lammin)
                    lami(k)=min(lami(k),lammax)
                else
                    lami(k)=0._r8
                end if

                if (dumc(i,k).ge.qsmall) then
                    ! add upper limit to in-cloud number concentration to prevent numerical error
                    dumnc(i,k)=min(dumnc(i,k),dumc(i,k)*1.e20_r8)

                    pgam(k)=0.0005714_r8*(dumnc(i,k)/1.e6_r8/rho(i,k))+0.2714_r8
                    pgam(k)=1._r8/(pgam(k)**2)-1._r8
                    pgam(k)=max(pgam(k),2._r8)
                    pgam(k)=min(pgam(k),15._r8)

                    lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
                        (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
                    lammin = (pgam(k)+1._r8)/50.e-6_r8
                    lammax = (pgam(k)+1._r8)/2.e-6_r8
                    lamc(k)=max(lamc(k),lammin)
                    lamc(k)=min(lamc(k),lammax)
                else
                    lamc(k)=0._r8
                end if

                ! calculate number and mass weighted fall velocity for droplets
                ! include effects of sub-grid distribution of cloud water


                !if (dumc(i,k).ge.qsmall) then
                if (dumc(i,k).gt.qsmall) then           !!!sxj---------
                    unc = gamma(qcvar+bc/3._r8)/(gamma(qcvar)*qcvar**(bc/3._r8))* &
                        acn(i,k)*gamma(1._r8+bc+pgam(k))/ &
                        (lamc(k)**bc*gamma(pgam(k)+1._r8))
                    umc = gamma(qcvar+bc/3._r8)/(gamma(qcvar)*qcvar**(bc/3._r8))* &
                        acn(i,k)*gamma(4._r8+bc+pgam(k))/ &
                        (lamc(k)**bc*gamma(pgam(k)+4._r8))
                else
                    umc = 0._r8
                    unc = 0._r8
                end if

                ! calculate number and mass weighted fall velocity for cloud ice

                !if (dumi(i,k).ge.qsmall) then       !!!!sxj----------
                if (dumi(i,k).gt.qsmall) then
                    uni =  ain(i,k)*gamma(1._r8+bi)/lami(k)**bi
                    umi = ain(i,k)*gamma(4._r8+bi)/(6._r8*lami(k)**bi)
                    uni=min(uni,1.2_r8)
                    umi=min(umi,1.2_r8)
                else
                    umi = 0._r8
                    uni = 0._r8
                end if

                fi(k) = g*rho(i,k)*umi
                fni(k) = g*rho(i,k)*uni
                fc(k) = g*rho(i,k)*umc
                fnc(k) = g*rho(i,k)*unc

                ! calculate number of split time steps to ensure courant stability criteria
                ! for sedimentation calculations

                rgvm = max(fi(k),fc(k),fni(k),fnc(k))
                nstep = max(int(rgvm*deltat/pdel(i,k)+1._r8),nstep)

                ! redefine dummy variables - sedimentation is calculated over grid-scale
                ! quantities to ensure conservation

                dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
                dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
                dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat),0._r8)
                dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat),0._r8)

                if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
                if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8

            end do       !!! vertical loop

!!!!!!!!!!
            do n = 1,nstep  !! loop over sub-time step to ensure stability		
                do k = 1,pver
                    falouti(k) = fi(k)*dumi(i,k)
                    faloutni(k) = fni(k)*dumni(i,k)
                    faloutc(k) = fc(k)*dumc(i,k)
                    faloutnc(k) = fnc(k)*dumnc(i,k)
                end do

                ! top of model

                k = 1
                faltndi = falouti(k)/pdel(i,k)
                faltndni = faloutni(k)/pdel(i,k)
                faltndc = faloutc(k)/pdel(i,k)
                faltndnc = faloutnc(k)/pdel(i,k)

                ! add fallout terms to microphysical tendencies	

                qitend(i,k) = qitend(i,k)-faltndi/nstep
                nitend(i,k) = nitend(i,k)-faltndni/nstep
                qctend(i,k) = qctend(i,k)-faltndc/nstep
                nctend(i,k) = nctend(i,k)-faltndnc/nstep

                dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
                dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep
                dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
                dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep

                do k = 2,pver

                    ! for cloud liquid and ice, if cloud fraction increases with height
                    ! then add flux from above to both vapor and cloud water of current level
                    ! this means that flux entering clear portion of cell from above evaporates
                    ! instantly

                    dum=cldm(i,k)/cldm(i,k-1)
                    dum=min(dum,1._r8)

                    faltndqie=(falouti(k)-falouti(k-1))/pdel(i,k) ! 
                    faltndi=(falouti(k)-dum*falouti(k-1))/pdel(i,k)
                    faltndni=(faloutni(k)-dum*faloutni(k-1))/pdel(i,k)
                    faltndqce=(faloutc(k)-faloutc(k-1))/pdel(i,k)
                    faltndc=(faloutc(k)-dum*faloutc(k-1))/pdel(i,k)
                    faltndnc=(faloutnc(k)-dum*faloutnc(k-1))/pdel(i,k)

                    ! add fallout terms to eulerian tendencies

                    qitend(i,k) = qitend(i,k)-faltndi/nstep
                    nitend(i,k) = nitend(i,k)-faltndni/nstep
                    qctend(i,k) = qctend(i,k)-faltndc/nstep
                    nctend(i,k) = nctend(i,k)-faltndnc/nstep

                    ! add terms to to evap/sub of cloud water

                    qvlat(i,k)=qvlat(i,k)-(faltndqie-faltndi)/nstep
                    qvlat(i,k)=qvlat(i,k)-(faltndqce-faltndc)/nstep
                    tlat(i,k)=tlat(i,k)+(faltndqie-faltndi)*xxls/nstep
                    tlat(i,k)=tlat(i,k)+(faltndqce-faltndc)*xxlv/nstep

                    dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
                    dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep
                    dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
                    dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep

                    Fni(K)=MAX(Fni(K)/pdel(i,K),Fni(K-1)/pdel(i,K-1))*pdel(i,K)
                    FI(K)=MAX(FI(K)/pdel(i,K),FI(K-1)/pdel(i,K-1))*pdel(i,K)
                    fnc(k)=max(fnc(k)/pdel(i,k),fnc(k-1)/pdel(i,k-1))*pdel(i,k)
                    Fc(K)=MAX(Fc(K)/pdel(i,K),Fc(K-1)/pdel(i,K-1))*pdel(i,K)

                end do   !! k loop

                ! units below are m/s
                ! cloud water/ice sedimentation flux at surface 
                ! is added to precip flux at surface to get total precip (cloud + precip water)
                ! rate

                prect(i) = prect(i)+(faloutc(pver)+falouti(pver)) &
                    /g/nstep/1000._r8	
                preci(i) = preci(i)+(falouti(pver)) &
                    /g/nstep/1000._r8

                if (preci(i)<0.) then 					    
                    ! write(*,*) "preci(i)",preci(i)
                    preci(i)=0.
                endif
                if (prect(i)<0.) then 					    
                    ! write(*,*) "preci(i)",preci(i)
                    prect(i)=0.
                endif

            end do   !! nstep loop

            ! end sedimentation
            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



            ! get new update for variables that includes sedimentation tendency
            ! note : here dum variables are grid-average, NOT in-cloud

            do k=1,pver	

                dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)
                dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)
                dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)
                dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)

                if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
                if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8

                ! calculate instantaneous processes (melting, homogeneous freezing)

                if (t(i,k)+tlat(i,k)/cp*deltat > 273.15_r8) then
                    if (dumi(i,k) > 0._r8) then

                        ! limit so that melting does not push temperature below freezing
                        dum = -dumi(i,k)*xlf/cp
                        if (t(i,k)+tlat(i,k)/cp*deltat+dum.lt.273.15_r8) then
                            dum = (t(i,k)+tlat(i,k)/cp*deltat-273.15_r8)*cp/xlf
                            dum = dum/dumi(i,k)*xlf/cp 
                            dum = max(0._r8,dum)
                            dum = min(1._r8,dum)
                        else
                            dum = 1._r8
                        end if

                        qctend(i,k)=qctend(i,k)+dum*dumi(i,k)/deltat

                        ! hm add, 9/15/06, assume melting ice produces droplet
                        ! mean volume radius of 8 micron

                        nctend(i,k)=nctend(i,k)+3._r8*dum*dumi(i,k)/deltat/ &
                            (4._r8*pi*5.12e-16_r8*rhow)

                        qitend(i,k)=((1._r8-dum)*dumi(i,k)-qi(i,k))/deltat
                        nitend(i,k)=((1._r8-dum)*dumni(i,k)-ni(i,k))/deltat
                        tlat(i,k)=tlat(i,k)-xlf*dum*dumi(i,k)/deltat
                    end if
                end if

                ! homogeneously freeze droplets at -40 C

                if (t(i,k)+tlat(i,k)/cp*deltat < 233.15_r8) then
                    if (dumc(i,k) > 0._r8) then

                        ! limit so that freezing does not push temperature above threshold
                        dum = dumc(i,k)*xlf/cp
                        if (t(i,k)+tlat(i,k)/cp*deltat+dum.gt.233.15_r8) then
                            dum = -(t(i,k)+tlat(i,k)/cp*deltat-233.15_r8)*cp/xlf
                            dum = dum/dumc(i,k)*xlf/cp
                            dum = max(0._r8,dum)
                            dum = min(1._r8,dum)
                        else
                            dum = 1._r8
                        end if

                        qitend(i,k)=qitend(i,k)+dum*dumc(i,k)/deltat
                        ! hm add 11/18/06
                        ! assume 25 micron mean volume radius of homogeneously frozen droplets
                        ! consistent with size of detrained ice in stratiform.F90
                        nitend(i,k)=nitend(i,k)+dum*3._r8*dumc(i,k)/(4._r8*3.14_r8*1.563e-14_r8* &
                            500._r8)/deltat
                        qctend(i,k)=((1._r8-dum)*dumc(i,k)-qc(i,k))/deltat
                        nctend(i,k)=((1._r8-dum)*dumnc(i,k)-nc(i,k))/deltat
                        tlat(i,k)=tlat(i,k)+xlf*dum*dumc(i,k)/deltat
                    end if
                end if

                ! hm add 2/2/07.......................
                ! remove any excess over-saturation, which is possible due to non-linearity when adding 
                ! together all microphysical processes
                ! follow code similar to old CAM scheme

                qtmp=q(i,k)+qvlat(i,k)*deltat
                ttmp=t(i,k)+tlat(i,k)/cpp*deltat
                esn = estblf(ttmp)
                qsn = min(epsqs*esn/(p(i,k)-(1._r8-epsqs)*esn),1._r8)

                if (qtmp > qsn .and. qsn > 0) then
                    ! expression below is approximate since there may be ice deposition
                    dum = (qtmp-qsn)/(1._r8+xxlv**2*qsn/(cpp*rv*ttmp**2))/deltat
                    ! add to output cme
                    cmeout(i,k) = cmeout(i,k)+dum
                    ! now add to tendencies, partition between liquid and ice based on temperature
                    if (ttmp > 268.15_r8) then
                        dum1=0.0_r8
                        ! now add to tendencies, partition between liquid and ice based on te
                    else if (ttmp < 238.15_r8) then
                        dum1=1.0_r8
                    else
                        dum1=(268.15_r8-ttmp)/30._r8
                    end if

                    dum = (qtmp-qsn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
                        *qsn/(cpp*rv*ttmp**2))/deltat
                    qctend(i,k)=qctend(i,k)+dum*(1._r8-dum1)
                    qitend(i,k)=qitend(i,k)+dum*dum1
                    qvlat(i,k)=qvlat(i,k)-dum
                    tlat(i,k)=tlat(i,k)+dum*(1._r8-dum1)*xxlv+dum*dum1*xxls
                end if

                !...............................................................................
                ! calculate effective radius for pass to radiation code
                ! if no cloud water, default value is 10 micron for droplets,
                ! 25 micron for cloud ice

                ! update cloud variables after instantaneous processes to get effective radius
                ! variables are in-cloud to calculate size dist parameters

                dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)/cldm(i,k)
                dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)/cldm(i,k)
                dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)/cldm(i,k)
                dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)/cldm(i,k)

                ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1

                dumc(i,k)=min(dumc(i,k),5.e-3_r8)
                dumi(i,k)=min(dumi(i,k),5.e-3_r8)

                !...................
                ! cloud ice effective radius

                if (dumi(i,k).ge.qsmall) then
                    ! add upper limit to in-cloud number concentration to prevent numerical error
                    dumni(i,k)=min(dumni(i,k),dumi(i,k)*1.e20_r8)
                    lami(k) = (gamma(1._r8+di)*ci* &
                        dumni(i,k)/dumi(i,k))**(1._r8/di)
                    lammax = 1._r8/10.e-6_r8
                    lammin = 1._r8/(2._r8*dcs)

                    if (lami(k).lt.lammin) then
                        lami(k) = lammin
                        n0i(k) = lami(k)**(di+1._r8)*dumi(i,k)/(ci*gamma(1._r8+di))
                        niic(i,k) = n0i(k)/lami(k)
                        ! adjust number conc if needed to keep mean size in reasonable range
                        nitend(i,k)=(niic(i,k)*cldm(i,k)-ni(i,k))/deltat
                    else if (lami(k).gt.lammax) then
                        lami(k) = lammax
                        n0i(k) = lami(k)**(di+1._r8)*dumi(i,k)/(ci*gamma(1._r8+di))
                        niic(i,k) = n0i(k)/lami(k)
                        ! adjust number conc if needed to keep mean size in reasonable range
                        nitend(i,k)=(niic(i,k)*cldm(i,k)-ni(i,k))/deltat
                    end if
                    effi(i,k) = 1.5_r8/lami(k)*1.e6_r8
                else
                    effi(i,k) = 25._r8
                end if

                !...................
                ! cloud droplet effective radius

                if (dumc(i,k).ge.qsmall) then
                    ! add upper limit to in-cloud number concentration to prevent numerical error
                    dumnc(i,k)=min(dumnc(i,k),dumc(i,k)*1.e20_r8)
                    pgam(k)=0.0005714_r8*(dumnc(i,k)/1.e6_r8/rho(i,k))+0.2714_r8
                    pgam(k)=1._r8/(pgam(k)**2)-1._r8
                    pgam(k)=max(pgam(k),2._r8)
                    pgam(k)=min(pgam(k),15._r8)

                    lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
                        (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
                    lammin = (pgam(k)+1._r8)/50.e-6_r8
                    lammax = (pgam(k)+1._r8)/2.e-6_r8
                    if (lamc(k).lt.lammin) then
                        lamc(k) = lammin
                        ncic(i,k) = 6._r8*lamc(k)**3*dumc(i,k)* &
                            gamma(pgam(k)+1._r8)/ &
                            (pi*rhow*gamma(pgam(k)+4._r8))
                        ! adjust number conc if needed to keep mean size in reasonable range
                        nctend(i,k)=(ncic(i,k)*cldm(i,k)-nc(i,k))/deltat

                    else if (lamc(k).gt.lammax) then
                        lamc(k) = lammax
                        ncic(i,k) = 6._r8*lamc(k)**3*dumc(i,k)* &
                            gamma(pgam(k)+1._r8)/ &
                            (pi*rhow*gamma(pgam(k)+4._r8))
                        ! adjust number conc if needed to keep mean size in reasonable range
                        nctend(i,k)=(ncic(i,k)*cldm(i,k)-nc(i,k))/deltat
                    end if
                    effc(i,k) = gamma(qcvar+1._r8/3._r8)/(gamma(qcvar)*qcvar**(1._r8/3._r8))* &  ! Sub-grid variability is considered for cloud water
                        gamma(pgam(k)+4._r8)/ &
                        gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8
                else
                    effc(i,k) = 10._r8
                end if


!!! recalculate effective radius for constant number, in order to separate
                ! first and second indirect effects
                ! assume constant number of 10^8 kg-1

                dumnc(i,k)=1.e8

                if (dumc(i,k).ge.qsmall) then
                    pgam(k)=0.0005714_r8*(dumnc(i,k)/1.e6_r8/rho(i,k))+0.2714_r8
                    pgam(k)=1._r8/(pgam(k)**2)-1._r8
                    pgam(k)=max(pgam(k),2._r8)
                    pgam(k)=min(pgam(k),15._r8)

                    lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
                        (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
                    lammin = (pgam(k)+1._r8)/50.e-6_r8
                    lammax = (pgam(k)+1._r8)/2.e-6_r8
                    if (lamc(k).lt.lammin) then
                        lamc(k) = lammin
                    else if (lamc(k).gt.lammax) then
                        lamc(k) = lammax
                    end if
                    effc_fn(i,k) = gamma(qcvar+1._r8/3._r8)/ &
                        (gamma(qcvar)*qcvar**(1._r8/3._r8))* &
                        gamma(pgam(k)+4._r8)/ &
                        gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8
                else
                    effc_fn(i,k) = 10._r8
                end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!

            end do ! vertical k loop



500         continue

            do k=1,pver
                ! if updated q (after microphysics) is zero, then ensure updated n is also zero

                if (qc(i,k)+qctend(i,k)*deltat.lt.qsmall) nctend(i,k)=-nc(i,k)/deltat
                if (qi(i,k)+qitend(i,k)*deltat.lt.qsmall) nitend(i,k)=-ni(i,k)/deltat
            end do

	end do ! i loop

        !        ccn concentration as diagnostic
        ! ccn(pcols,pver,psat) ---number conc of aerosols activated at supersat

        do k=1,pver
            ccn(:ncol,k,:) = 0.
            do m=1,naer_all
                !	       naer(:ncol)=aer_mmr(:ncol,k,m)*num_to_mass_aer(m)*rho(:ncol,k)
                do l=1,psat ! psat=6
                    !                  ccn(:ncol,k,l)=ccn(:ncol,k,l)+naer(:ncol)*ccnfact(l,m)
                    ccn(:ncol,k,l)=ccn(:ncol,k,l)+naer2(:ncol,k,m)*ccnfact(l,m)
                enddo
            enddo
        enddo

        do l=1,psat
            call outfld(ccn_name(l), ccn(1:pcols,1:pver,l)    , pcols, lchnk   )
        enddo

        ! hm add rain/snow mixing ratio and number concentration as diagnostic

        call outfld('QRAIN',qrout,   pcols, lchnk) !
        call outfld('QSNOW',qsout,   pcols, lchnk)
        call outfld('NRAIN',nrout,   pcols, lchnk)
        call outfld('NSNOW',nsout,   pcols, lchnk)

        !++ag (1.90) averaging for snow and rain number
        !added 2.02: rain and snow mean diameter.

        qrout2(:,:)=0._r8
        qsout2(:,:)=0._r8
        nrout2(:,:)=0._r8
        nsout2(:,:)=0._r8	
        drout2(:,:)=0._r8
        dsout2(:,:)=0._r8	
        freqs(:,:)=0._r8
        freqr(:,:)=0._r8
        do i = 1,ncol
            do k=1,pver
                if (qrout(i,k).gt.1.e-7_r8.and.nrout(i,k).gt.0._r8) then
                    qrout2(i,k)=qrout(i,k)
                    nrout2(i,k)=nrout(i,k)
                    drout2(i,k)=(pi * rhow * nrout(i,k)/qrout(i,k))**(-1._r8/3._r8)
                    freqr(i,k)=1._r8
                endif
                if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
                    qsout2(i,k)=qsout(i,k)
                    nsout2(i,k)=nsout(i,k)
                    dsout2(i,k)=(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)
                    freqs(i,k)=1._r8
                endif
            end do
        end do

        !ag add averaged out fields.
        call outfld('AQRAIN',qrout2,    pcols,lchnk)
        call outfld('AQSNOW',qsout2,    pcols,lchnk)
        call outfld('ANRAIN',nrout2,    pcols,lchnk)
        call outfld('ANSNOW',nsout2,    pcols,lchnk)
        call outfld('ADRAIN',drout2,    pcols,lchnk)
        call outfld('ADSNOW',dsout2,    pcols,lchnk)
        call outfld('FREQR',freqr,    pcols,lchnk)
        call outfld('FREQS',freqs,    pcols,lchnk)

        !++ag redefine fice here....
	nfice(:,:)=0._r8
	do k=1,pver
            do i=1,ncol
        	dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
        	dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
                dumfice=qsout(i,k) + qrout(i,k) + dumc(i,k) + dumi(i,k)  
		if (dumfice.gt.0._r8) then
		    nfice(i,k)=(qsout(i,k) + dumi(i,k))/dumfice
                    fice(i,k)=nfice(i,k) !sxj-------
	        endif
            enddo
	enddo


        call outfld('FICE',nfice,   pcols, lchnk)


        return

    endsubroutine mmicro_pcond


    !##############################################################################

    subroutine findsp1 (lchnk, ncol, q, t, p, tsp, qsp)
        !----------------------------------------------------------------------- 
        ! 
        ! Purpose: 
        !     find the wet bulb temperature for a given t and q
        !     in a longitude height section
        !     wet bulb temp is the temperature and spec humidity that is 
        !     just saturated and has the same enthalpy
        !     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
        !     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
        !
        ! Method: 
        ! a Newton method is used
        ! first guess uses an algorithm provided by John Petch from the UKMO
        ! we exclude points where the physical situation is unrealistic
        ! e.g. where the temperature is outside the range of validity for the
        !      saturation vapor pressure, or where the water vapor pressure
        !      exceeds the ambient pressure, or the saturation specific humidity is 
        !      unrealistic
        ! 
        ! Author: P. Rasch
        ! 
        !-----------------------------------------------------------------------
        !
        !     input arguments
        !
        integer, intent(in) :: lchnk                 ! chunk identifier
        integer, intent(in) :: ncol                  ! number of atmospheric columns

        real(r8), intent(in) :: q(pcols,pver)        ! water vapor (kg/kg)
        real(r8), intent(in) :: t(pcols,pver)        ! temperature (K)
        real(r8), intent(in) :: p(pcols,pver)        ! pressure    (Pa)
        !
        ! output arguments
        !
        real(r8), intent(out) :: tsp(pcols,pver)      ! saturation temp (K)
        real(r8), intent(out) :: qsp(pcols,pver)      ! saturation mixing ratio (kg/kg)
        !
        ! local variables
        !
        integer i                 ! work variable
        integer k                 ! work variable
        logical lflg              ! work variable
        integer iter              ! work variable
        integer l                 ! work variable
        logical :: error_found

        real(r8) omeps                ! 1 minus epsilon
        real(r8) trinv                ! work variable
        real(r8) es                   ! sat. vapor pressure
        real(r8) desdt                ! change in sat vap pressure wrt temperature
        !     real(r8) desdp                ! change in sat vap pressure wrt pressure
        real(r8) dqsdt                ! change in sat spec. hum. wrt temperature
        real(r8) dgdt                 ! work variable
        real(r8) g                    ! work variable
        real(r8) weight(pcols)        ! work variable
        real(r8) hlatsb               ! (sublimation)
        real(r8) hlatvp               ! (vaporization)
        real(r8) hltalt(pcols,pver)   ! lat. heat. of vap.
        real(r8) tterm                ! work var.
        real(r8) qs                   ! spec. hum. of water vapor
        real(r8) tc                   ! crit temp of transition to ice

        ! work variables
        real(r8) t1, q1, dt, dq
        real(r8) dtm, dqm
        real(r8) qvd, a1, tmp
        real(r8) rair
        real(r8) r1b, c1, c2, c3
        real(r8) denom
        real(r8) dttol
        real(r8) dqtol
        integer doit(pcols) 
        real(r8) enin(pcols), enout(pcols)
        real(r8) tlim(pcols)

        omeps = 1.0_r8 - epsqs
        trinv = 1.0_r8/ttrice
        a1 = 7.5_r8*log(10._r8)
        rair =  287.04_r8
        c3 = rair*a1/cp
        dtm = 0._r8    ! needed for iter=0 blowup with f90 -ei
        dqm = 0._r8    ! needed for iter=0 blowup with f90 -ei
        dttol = 1.e-4_r8 ! the relative temp error tolerance required to quit the iteration
        dqtol = 1.e-4_r8 ! the relative moisture error tolerance required to quit the iteration
        !  tmin = 173.16 ! the coldest temperature we can deal with
        !
        ! max number of times to iterate the calculation
        iter = 8
        !
        do k = 1,pver

            !
            ! first guess on the wet bulb temperature
            !
            do i = 1,ncol

#ifdef DEBUG
                if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
                    write (6,*) ' '
                    write (6,*) ' level, t, q, p', k, t(i,k), q(i,k), p(i,k)
                endif
#endif

                ! limit the temperature range to that relevant to the sat vap pres tables
#if ( ! defined WACCM_MOZART )
                tlim(i) = min(max(t(i,k),173._r8),373._r8)
#else
                tlim(i) = min(max(t(i,k),128._r8),373._r8)
#endif
                es = estblf(tlim(i))
                denom = p(i,k) - omeps*es
                qs = epsqs*es/denom
                doit(i) = 0
                enout(i) = 1._r8
                ! make sure a meaningful calculation is possible
                if (p(i,k) > 5._r8*es .and. qs > 0._r8 .and. qs < 0.5_r8) then
                    !
                    ! Saturation specific humidity
                    !
                    qs = min(epsqs*es/denom,1._r8)
                    !
                    ! "generalized" analytic expression for t derivative of es
                    !  accurate to within 1 percent for 173.16 < t < 373.16
                    !
                    ! Weighting of hlat accounts for transition from water to ice
                    ! polynomial expression approximates difference between es over
                    ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
                    ! -40): required for accurate estimate of es derivative in transition
                    ! range from ice to water also accounting for change of hlatv with t
                    ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
                    !
                    tc     = tlim(i) - t0
                    lflg   = (tc >= -ttrice .and. tc < 0.0_r8)
                    weight(i) = min(-tc*trinv,1.0_r8)
                    hlatsb = hlatv + weight(i)*hlatf
                    hlatvp = hlatv - 2369.0_r8*tc
                    if (tlim(i) < t0) then
                        hltalt(i,k) = hlatsb
                    else
                        hltalt(i,k) = hlatvp
                    end if
                    enin(i) = cp*tlim(i) + hltalt(i,k)*q(i,k)

                    ! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
                    tmp =  q(i,k) - qs
                    c1 = hltalt(i,k)*c3
                    c2 = (tlim(i) + 36._r8)**2
                    r1b    = c2/(c2 + c1*qs)
                    qvd   = r1b*tmp
                    tsp(i,k) = tlim(i) + ((hltalt(i,k)/cp)*qvd)
#ifdef DEBUG
                    if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
                        write (6,*) ' relative humidity ', q(i,k)/qs
                        write (6,*) ' first guess ', tsp(i,k)
                    endif
#endif
                    es = estblf(tsp(i,k))
                    qsp(i,k) = min(epsqs*es/(p(i,k) - omeps*es),1._r8)
                else
                    doit(i) = 1
                    tsp(i,k) = tlim(i)
                    qsp(i,k) = q(i,k)
                    enin(i) = 1._r8
                endif
            end do   ! end do i
            !
            ! now iterate on first guess
            !
            do l = 1, iter
                dtm = 0
                dqm = 0
                do i = 1,ncol
                    if (doit(i) == 0) then
                        es = estblf(tsp(i,k))
                        !
                        ! Saturation specific humidity
                        !
                        qs = min(epsqs*es/(p(i,k) - omeps*es),1._r8)
                        !
                        ! "generalized" analytic expression for t derivative of es
                        ! accurate to within 1 percent for 173.16 < t < 373.16
                        !
                        ! Weighting of hlat accounts for transition from water to ice
                        ! polynomial expression approximates difference between es over
                        ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
                        ! -40): required for accurate estimate of es derivative in transition
                        ! range from ice to water also accounting for change of hlatv with t
                        ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
                        !
                        tc     = tsp(i,k) - t0
                        lflg   = (tc >= -ttrice .and. tc < 0.0_r8)
                        weight(i) = min(-tc*trinv,1.0_r8)
                        hlatsb = hlatv + weight(i)*hlatf
                        hlatvp = hlatv - 2369.0_r8*tc
                        if (tsp(i,k) < t0) then
                            hltalt(i,k) = hlatsb
                        else
                            hltalt(i,k) = hlatvp
                        end if
                        if (lflg) then
                            tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3)+tc*(pcf(4) + tc*pcf(5))))
                        else
                            tterm = 0.0_r8
                        end if
                        desdt = hltalt(i,k)*es/(rgasv*tsp(i,k)*tsp(i,k)) + tterm*trinv
                        dqsdt = (epsqs + omeps*qs)/(p(i,k) - omeps*es)*desdt
                        !              g = cp*(tlim(i)-tsp(i,k)) + hltalt(i,k)*q(i,k)- hltalt(i,k)*qsp(i,k)
                        g = enin(i) - (cp*tsp(i,k) + hltalt(i,k)*qsp(i,k))
                        dgdt = -(cp + hltalt(i,k)*dqsdt)
                        t1 = tsp(i,k) - g/dgdt
                        dt = abs(t1 - tsp(i,k))/t1
                        tsp(i,k) = max(t1,tmin)
                        es = estblf(tsp(i,k))
                        q1 = min(epsqs*es/(p(i,k) - omeps*es),1._r8)
                        dq = abs(q1 - qsp(i,k))/max(q1,1.e-12_r8)
                        qsp(i,k) = q1
#ifdef DEBUG
                        if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
                            write (6,*) ' rel chg lev, iter, t, q ', k, l, dt, dq, g
                        endif
#endif
                        dtm = max(dtm,dt)
                        dqm = max(dqm,dq)
                        ! if converged at this point, exclude it from more iterations
                        if (dt < dttol .and. dq < dqtol) then
                            doit(i) = 2
                        endif
                        enout(i) = cp*tsp(i,k) + hltalt(i,k)*qsp(i,k)
                        ! bail out if we are too near the end of temp range
#if ( ! defined WACCM_MOZART )
                        if (tsp(i,k) < 174.16_r8) then
#else
                            if (tsp(i,k) < 130.16_r8) then
#endif
                                doit(i) = 4
                            endif
                        else
                        endif
                    end do              ! do i = 1,ncol

                    if (dtm < dttol .and. dqm < dqtol) then
                        go to 10
                    endif

                end do                 ! do l = 1,iter
10              continue

                error_found = .false.
                if (dtm > dttol .or. dqm > dqtol) then
                    do i = 1,ncol
                        if (doit(i) == 0) error_found = .true.
                    end do
                    if (error_found) then
                        do i = 1,ncol
                            if (doit(i) == 0) then
                                write (6,*) ' findsp not converging at point i, k ', i, k
                                write (6,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
                                write (6,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
                                call endrun ('FINDSP')
                            endif
                        end do
                    endif
                endif
                do i = 1,ncol
                    if (doit(i) == 2 .and. abs((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) then
                        error_found = .true.
                    endif
                end do
                if (error_found) then
                    do i = 1,ncol
                        if (doit(i) == 2 .and. abs((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) then
                            write (6,*) ' the enthalpy is not conserved for point ', &
                                i, k, enin(i), enout(i)
                            write (6,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
                            write (6,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
                            call endrun ('FINDSP')
                        endif
                    end do
                endif

            end do                    ! level loop (k=1,pver)

            return
        end subroutine findsp1

        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! error function in single precision
        !
        !    Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
        !    You may use, copy, modify this code for any purpose and 
        !    without fee. You may distribute this ORIGINAL package.

        function derf(x)
            implicit real (a - h, o - z)
            real(r8) a,b,x
            dimension a(0 : 64), b(0 : 64)
            integer i,k
            data (a(i), i = 0, 12) / & 
                0.00000000005958930743d0, -0.00000000113739022964d0, & 
                0.00000001466005199839d0, -0.00000016350354461960d0, &
                0.00000164610044809620d0, -0.00001492559551950604d0, &
                0.00012055331122299265d0, -0.00085483269811296660d0, &
                0.00522397762482322257d0, -0.02686617064507733420d0, &
                0.11283791670954881569d0, -0.37612638903183748117d0, &
                1.12837916709551257377d0 / 
            data (a(i), i = 13, 25) / &
                0.00000000002372510631d0, -0.00000000045493253732d0, &
                0.00000000590362766598d0, -0.00000006642090827576d0, &
                0.00000067595634268133d0, -0.00000621188515924000d0, &
                0.00005103883009709690d0, -0.00037015410692956173d0, &
                0.00233307631218880978d0, -0.01254988477182192210d0, &
                0.05657061146827041994d0, -0.21379664776456006580d0, &
                0.84270079294971486929d0 / 
            data (a(i), i = 26, 38) / &
                0.00000000000949905026d0, -0.00000000018310229805d0, &
                0.00000000239463074000d0, -0.00000002721444369609d0, &
                0.00000028045522331686d0, -0.00000261830022482897d0, &
                0.00002195455056768781d0, -0.00016358986921372656d0, &
                0.00107052153564110318d0, -0.00608284718113590151d0, &
                0.02986978465246258244d0, -0.13055593046562267625d0, &
                0.67493323603965504676d0 / 
            data (a(i), i = 39, 51) / &
                0.00000000000382722073d0, -0.00000000007421598602d0, &
                0.00000000097930574080d0, -0.00000001126008898854d0, &
                0.00000011775134830784d0, -0.00000111992758382650d0, &
                0.00000962023443095201d0, -0.00007404402135070773d0, &
                0.00050689993654144881d0, -0.00307553051439272889d0, &
                0.01668977892553165586d0, -0.08548534594781312114d0, &
                0.56909076642393639985d0 / 
            data (a(i), i = 52, 64) / &
                0.00000000000155296588d0, -0.00000000003032205868d0, &
                0.00000000040424830707d0, -0.00000000471135111493d0, &
                0.00000005011915876293d0, -0.00000048722516178974d0, &
                0.00000430683284629395d0, -0.00003445026145385764d0, &
                0.00024879276133931664d0, -0.00162940941748079288d0, &
                0.00988786373932350462d0, -0.05962426839442303805d0, &
                0.49766113250947636708d0 / 
            data (b(i), i = 0, 12) / &
                -0.00000000029734388465d0, 0.00000000269776334046d0, &
                -0.00000000640788827665d0, -0.00000001667820132100d0, &
                -0.00000021854388148686d0, 0.00000266246030457984d0, &
                0.00001612722157047886d0, -0.00025616361025506629d0, &
                0.00015380842432375365d0, 0.00815533022524927908d0, &
                -0.01402283663896319337d0, -0.19746892495383021487d0,& 
                0.71511720328842845913d0 / 
            data (b(i), i = 13, 25) / &
                -0.00000000001951073787d0, -0.00000000032302692214d0, &
                0.00000000522461866919d0, 0.00000000342940918551d0, &
                -0.00000035772874310272d0, 0.00000019999935792654d0, &
                0.00002687044575042908d0, -0.00011843240273775776d0, &
                -0.00080991728956032271d0, 0.00661062970502241174d0, &
                0.00909530922354827295d0, -0.20160072778491013140d0, &
                0.51169696718727644908d0 / 
            data (b(i), i = 26, 38) / &
                0.00000000003147682272d0, -0.00000000048465972408d0, &
                0.00000000063675740242d0, 0.00000003377623323271d0, &
                -0.00000015451139637086d0, -0.00000203340624738438d0,& 
                0.00001947204525295057d0, 0.00002854147231653228d0, &
                -0.00101565063152200272d0, 0.00271187003520095655d0, &
                0.02328095035422810727d0, -0.16725021123116877197d0, &
                0.32490054966649436974d0 / 
            data (b(i), i = 39, 51) / &
                0.00000000002319363370d0, -0.00000000006303206648d0, &
                -0.00000000264888267434d0, 0.00000002050708040581d0, &
                0.00000011371857327578d0, -0.00000211211337219663d0, &
                0.00000368797328322935d0, 0.00009823686253424796d0, &
                -0.00065860243990455368d0, -0.00075285814895230877d0,& 
                0.02585434424202960464d0, -0.11637092784486193258d0, &
                0.18267336775296612024d0 / 
            data (b(i), i = 52, 64) / &
                -0.00000000000367789363d0, 0.00000000020876046746d0, &
                -0.00000000193319027226d0, -0.00000000435953392472d0, &
                0.00000018006992266137d0, -0.00000078441223763969d0, &
                -0.00000675407647949153d0, 0.00008428418334440096d0, &
                -0.00017604388937031815d0, -0.00239729611435071610d0, &
                0.02064129023876022970d0, -0.06905562880005864105d0, &
                0.09084526782065478489d0 / 
            w = abs(x)
            if (w .lt. 2.2d0) then
                t = w * w
                k = int(t)
                t = t - k
                k = k * 13
                y = ((((((((((((a(k) * t + a(k + 1)) * t + &
                    a(k + 2)) * t + a(k + 3)) * t + a(k + 4)) * t + &
                    a(k + 5)) * t + a(k + 6)) * t + a(k + 7)) * t + &
                    a(k + 8)) * t + a(k + 9)) * t + a(k + 10)) * t + &
                    a(k + 11)) * t + a(k + 12)) * w
            else if (w .lt. 6.9d0) then
                k = int(w)
                t = w - k
                k = 13 * (k - 2)
                y = (((((((((((b(k) * t + b(k + 1)) * t + &
                    b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + &
                    b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + &
                    b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + &
                    b(k + 11)) * t + b(k + 12)
                y = y * y
                y = y * y
                y = y * y
                y = 1 - y * y
            else
                y = 1
            end if
            if (x .lt. 0) y = -y
            derf = y
        end function derf
        !
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        function polysvp (T,type)

            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            !  Compute saturation vapor pressure by using 
            ! function from Goff and Gatch (1946)

            !  Polysvp returned in units of pa.
            !  T is input in units of K.
            !  type refers to saturation with respect to liquid (0) or ice (1)

            implicit none

            real(r8) dum

            real(r8) T,polysvp

            integer type

            ! ice

            if (type.eq.1) then

                ! Goff Gatch equation (good down to -100 C)

                polysvp = 10._r8**(-9.09718_r8*(273.16_r8/t-1._r8)-3.56654_r8* &
                    log10(273.16_r8/t)+0.876793_r8*(1._r8-t/273.16_r8)+ &
                    log10(6.1071_r8))*100._r8

            end if

            ! Goff Gatch equation, uncertain below -70 C

            if (type.eq.0) then
                polysvp = 10._r8**(-7.90298_r8*(373.16_r8/t-1._r8)+ &
                    5.02808_r8*log10(373.16_r8/t)- &
                    1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t/373.16_r8))-1._r8)+ &
                    8.1328e-3_r8*(10._r8**(-3.49149_r8*(373.16_r8/t-1._r8))-1._r8)+ &
                    log10(1013.246_r8))*100._r8                      
            end if


        end function polysvp

        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        FUNCTION GAMMA(X)

            !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            !D    DOUBLE PRECISION FUNCTION DGAMMA(X)
            !----------------------------------------------------------------------
            !
            ! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
            !   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
            !   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
            !   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
            !   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
            !   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
            !   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
            !   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
            !   MACHINE-DEPENDENT CONSTANTS.
            !
            !
            !*******************************************************************
            !*******************************************************************
            !
            ! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
            !
            ! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
            ! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
            ! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
            !          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
            !                  GAMMA(XBIG) = BETA**MAXEXP
            ! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
            !          APPROXIMATELY BETA**MAXEXP
            ! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
            !          1.0+EPS .GT. 1.0
            ! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
            !          1/XMININ IS MACHINE REPRESENTABLE
            !
            !     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
            !
            !                            BETA       MAXEXP        XBIG
            !
            ! CRAY-1         (S.P.)        2         8191        966.961
            ! CYBER 180/855
            !   UNDER NOS    (S.P.)        2         1070        177.803
            ! IEEE (IBM/XT,
            !   SUN, ETC.)   (S.P.)        2          128        35.040
            ! IEEE (IBM/XT,
            !   SUN, ETC.)   (D.P.)        2         1024        171.624
            ! IBM 3033       (D.P.)       16           63        57.574
            ! VAX D-FORMAT   (D.P.)        2          127        34.844
            ! VAX G-FORMAT   (D.P.)        2         1023        171.489
            !
            !                            XINF         EPS        XMININ
            !
            ! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
            ! CYBER 180/855
            !   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
            ! IEEE (IBM/XT,
            !   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
            ! IEEE (IBM/XT,
            !   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
            ! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
            ! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
            ! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
            !
            !*******************************************************************
            !*******************************************************************
            !
            ! ERROR RETURNS
            !
            !  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
            !     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
            !     TO BE FREE OF UNDERFLOW AND OVERFLOW.
            !
            !
            !  INTRINSIC FUNCTIONS REQUIRED ARE:
            !
            !     INT, DBLE, EXP, LOG, REAL, SIN
            !
            !
            ! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
            !              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
            !              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
            !              (ED.), SPRINGER VERLAG, BERLIN, 1976.
            !
            !              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
            !              SONS, NEW YORK, 1968.
            !
            !  LATEST MODIFICATION: OCTOBER 12, 1989
            !
            !  AUTHORS: W. J. CODY AND L. STOLTZ
            !           APPLIED MATHEMATICS DIVISION
            !           ARGONNE NATIONAL LABORATORY
            !           ARGONNE, IL 60439
            !
            !----------------------------------------------------------------------
            INTEGER I,N
            LOGICAL PARITY

            real(r8) gamma
            REAL(r8) &
                !D    DOUBLE PRECISION
                C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE, &
                TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
            DIMENSION C(7),P(8),Q(8)
            !----------------------------------------------------------------------
            !  MATHEMATICAL CONSTANTS
            !----------------------------------------------------------------------
            DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0_r8,0.5E0_r8,12.0E0_r8,2.0E0_r8,0.0E0_r8/, &
                SQRTPI/0.9189385332046727417803297E0_r8/, &
                PI/3.1415926535897932384626434E0_r8/
            !D    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
            !D   1     SQRTPI/0.9189385332046727417803297D0/,
            !D   2     PI/3.1415926535897932384626434D0/
            !----------------------------------------------------------------------
            !  MACHINE DEPENDENT PARAMETERS
            !----------------------------------------------------------------------
            DATA XBIG,XMININ,EPS/35.040E0_r8,1.18E-38_r8,1.19E-7_r8/, &
                XINF/3.4E38_r8/
            !D    DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
            !D   1     XINF/1.79D308/
            !----------------------------------------------------------------------
            !  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
            !     APPROXIMATION OVER (1,2).
            !----------------------------------------------------------------------
            DATA P/-1.71618513886549492533811E+0_r8,2.47656508055759199108314E+1_r8,&
                -3.79804256470945635097577E+2_r8,6.29331155312818442661052E+2_r8,&
                8.66966202790413211295064E+2_r8,-3.14512729688483675254357E+4_r8,&
                -3.61444134186911729807069E+4_r8,6.64561438202405440627855E+4_r8/
            DATA Q/-3.08402300119738975254353E+1_r8,3.15350626979604161529144E+2_r8,&
                -1.01515636749021914166146E+3_r8,-3.10777167157231109440444E+3_r8,&
                2.25381184209801510330112E+4_r8,4.75584627752788110767815E+3_r8,&
                -1.34659959864969306392456E+5_r8,-1.15132259675553483497211E+5_r8/
            !D    DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
            !D   1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
            !D   2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
            !D   3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
            !D    DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
            !D   1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
            !D   2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
            !D   3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
            !----------------------------------------------------------------------
            !  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
            !----------------------------------------------------------------------
            DATA C/-1.910444077728E-03_r8,8.4171387781295E-04_r8, &
                -5.952379913043012E-04_r8,7.93650793500350248E-04_r8,&
                -2.777777777777681622553E-03_r8,8.333333333333333331554247E-02_r8,&
                5.7083835261E-03_r8/
            !D    DATA C/-1.910444077728D-03,8.4171387781295D-04,
            !D   1     -5.952379913043012D-04,7.93650793500350248D-04,
            !D   2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
            !D   3      5.7083835261D-03/
            !----------------------------------------------------------------------
            !  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
            !----------------------------------------------------------------------
            CONV(I) = REAL(I,r8)
            !D    CONV(I) = DBLE(I)
            PARITY=.FALSE.
            FACT=ONE
            N=0
            Y=X
            IF(Y.LE.ZERO)THEN
                !----------------------------------------------------------------------
                !  ARGUMENT IS NEGATIVE
                !----------------------------------------------------------------------
                Y=-X
                Y1=AINT(Y)
                RES=Y-Y1
                IF(RES.NE.ZERO)THEN
                    IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
                    FACT=-PI/SIN(PI*RES)
                    Y=Y+ONE
                ELSE
                    RES=XINF
                    GOTO 900
                ENDIF
            ENDIF
            !----------------------------------------------------------------------
            !  ARGUMENT IS POSITIVE
            !----------------------------------------------------------------------
            IF(Y.LT.EPS)THEN
                !----------------------------------------------------------------------
                !  ARGUMENT .LT. EPS
                !----------------------------------------------------------------------
                IF(Y.GE.XMININ)THEN
                    RES=ONE/Y
                ELSE
                    RES=XINF
                    GOTO 900
                ENDIF
            ELSEIF(Y.LT.TWELVE)THEN
                Y1=Y
                IF(Y.LT.ONE)THEN
                    !----------------------------------------------------------------------
                    !  0.0 .LT. ARGUMENT .LT. 1.0
                    !----------------------------------------------------------------------
                    Z=Y
                    Y=Y+ONE
                ELSE
                    !----------------------------------------------------------------------
                    !  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
                    !----------------------------------------------------------------------
                    N=INT(Y)-1
                    Y=Y-CONV(N)
                    Z=Y-ONE
                ENDIF
                !----------------------------------------------------------------------
                !  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
                !----------------------------------------------------------------------
                XNUM=ZERO
                XDEN=ONE
                DO 260 I=1,8
                    XNUM=(XNUM+P(I))*Z
                    XDEN=XDEN*Z+Q(I)
260                 CONTINUE
                    RES=XNUM/XDEN+ONE
                    IF(Y1.LT.Y)THEN
                        !----------------------------------------------------------------------
                        !  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
                        !----------------------------------------------------------------------
                        RES=RES/Y1
                    ELSEIF(Y1.GT.Y)THEN
                        !----------------------------------------------------------------------
                        !  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
                        !----------------------------------------------------------------------
                        DO 290 I=1,N
                            RES=RES*Y
                            Y=Y+ONE
290                         CONTINUE
                        ENDIF
                    ELSE
                        !----------------------------------------------------------------------
                        !  EVALUATE FOR ARGUMENT .GE. 12.0,
                        !----------------------------------------------------------------------
                        IF(Y.LE.XBIG)THEN
                            YSQ=Y*Y
                            SUM=C(7)
                            DO 350 I=1,6
                                SUM=SUM/YSQ+C(I)
350                             CONTINUE
                                SUM=SUM/Y-Y+SQRTPI
                                SUM=SUM+(Y-HALF)*LOG(Y)
                                RES=EXP(SUM)
                            ELSE
                                RES=XINF
                                GOTO 900
                            ENDIF
                        ENDIF
                        !----------------------------------------------------------------------
                        !  FINAL ADJUSTMENTS AND RETURN
                        !----------------------------------------------------------------------
                        IF(PARITY)RES=-RES
                        IF(FACT.NE.ONE)RES=FACT/RES
900                     GAMMA=RES
                        !D900 DGAMMA = RES
                        RETURN
                        ! ---------- LAST LINE OF GAMMA ----------
                    END function gamma


                    subroutine activate(wbar, tair, rhoair,  &
                        na, ptype, ntype, pmode, nmode, ma, sigman, hygro, rhodry, &
                        nact)
                        !      calculates number, surface, and mass fraction of aerosols activated as CCN
                        !      calculates flux of cloud droplets, surface area, and aerosol mass into cloud
                        !      assumes an internal mixture within each of up to pmode multiple aerosol modes
                        !      a gaussiam spectrum of updrafts can be treated.

                        !      mks units

                        !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
                        !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

                        use physconst, only: rair, epsilo, cpair, rh2o, latvap, gravit,   &
                            rhoh2o, mwh2o, r_universal
                        use wv_saturation, only: estblf, epsqs

                        implicit none


                        !      input

                        integer pmode,ptype ! dimension of modes, types in modes
                        real(r8) wbar          ! grid cell mean vertical velocity (m/s)
                        real(r8) tair          ! air temperature (K)
                        real(r8) rhoair        ! air density (kg/m3)
                        real(r8) na(pmode)           ! aerosol number concentration (/m3)
                        integer ntype(pmode)      ! number of aerosol types
                        integer nmode      ! number of aerosol modes
                        real(r8) ma(ptype,pmode)     ! aerosol mass concentration (kg/m3)
                        real(r8) sigman(pmode)  ! geometric standard deviation of aerosol size distribution
                        real(r8) hygro(pmode)  ! hygroscopicity of aerosol mode
                        real(r8) rhodry(ptype,pmode) ! density of aerosol material

                        !      output
                        real(r8) nact      ! out number fraction of aerosols activated

                        !      local

                        !      external erf,erfc
                        !      real(r8) erf,erfc
#if (defined AIX)
#define ERF erf
#define ERFC erfc
#else
#define ERF derf
#define ERFC derfc
                        real(r8) derf,derfc
#endif

                        integer, parameter:: nx=200
                        integer, parameter:: maxmodes=naer_all !!sxj
                        integer iquasisect_option, isectional
                        real(r8) integ,integf
                        real(r8) surften       ! surface tension of water w/respect to air (N/m)
                        data surften/0.076/
                        save surften
                        real(r8) p0     ! reference pressure (Pa)
                        data p0/1013.25e2/
                        save p0
                        real(r8) xmin(maxmodes),xmax(maxmodes) ! ln(r) at section interfaces
                        real(r8) surfmin(maxmodes),surfmax(maxmodes) ! surface area at interfaces
                        real(r8) volmin(maxmodes),volmax(maxmodes) ! volume at interfaces
                        real(r8) surfc(maxmodes) ! surface concentration (m2/m3)
                        real(r8) volc(maxmodes) ! total aerosol volume  concentration (m3/m3)
                        real(r8) tmass ! total aerosol mass concentration (g/cm3)
                        real(r8) sign(maxmodes)    ! geometric standard deviation of size distribution
                        real(r8) rm ! number mode radius of aerosol at max supersat (cm)
                        real(r8) pres ! pressure (Pa)
                        real(r8) path ! mean free path (m)
                        real(r8) diff ! diffusivity (m2/s)
                        real(r8) conduct ! thermal conductivity (Joule/m/sec/deg)
                        real(r8) diff0,conduct0
                        real(r8) qs ! water vapor saturation mixing ratio
                        real(r8) dqsdt ! change in qs with temperature
                        real(r8) dqsdp ! change in qs with pressure
                        real(r8) gloc ! thermodynamic function (m2/s)
                        real(r8) zeta, eta(maxmodes)
                        real(r8) smc(maxmodes)
                        real(r8) lnsmax ! ln(smax)
                        real(r8) alpha
                        real(r8) gammaloc
                        real(r8) beta
                        real(r8) sqrtg
                        real(r8) alogam
                        real(r8) rlo,rhi,xint1,xint2,xint3,xint4
                        real(r8) w,wnuc,wb
                        real(r8) dfmin,dfmax,fnew,fold,fnmin,fnbar,fsbar,fmbar
                        real(r8) alw,sqrtalw
                        real(r8) smax
                        real(r8) x,arg
                        real(r8) xmincoeff,xcut,volcut,surfcut
                        real(r8) z,z1,z2,wf1,wf2,zf1,zf2,gf1,gf2,gf
                        real(r8) etafactor1,etafactor2(maxmodes),etafactor2max
                        real(r8) es
                        integer m,n
                        !++ag
                        real(r8) amcubeloc(maxmodes)
                        real(r8) lnsmloc(maxmodes)
                        !--ag
                        !      numerical integration parameters
                        real(r8) eps,fmax,sds
                        data eps/0.3/,fmax/0.99/,sds/3./
                        save eps,fmax,sds

                        if(maxmodes<pmode)then
                            write(6,*)'maxmodes,pmode in activate =',maxmodes,pmode
                            call endrun('activate')
                        endif
                        nact=0._r8
                        if(nmode.eq.1.and.na(1).lt.1.e-20)return
                        if(wbar.le.0.)return
                        pres=rair*rhoair*tair
                        diff0=0.211e-4*(p0/pres)*(tair/t0)**1.94
                        conduct0=(5.69+0.017*(tair-t0))*4.186e2*1.e-5 ! convert to J/m/s/deg

                        es = estblf(tair)
                        qs = epsilo*es/(pres-(1.0_r8 - epsqs)*es)
                        dqsdt=latvap/(rh2o*tair*tair)*qs
                        alpha=gravit*(latvap/(cpair*rh2o*tair*tair)-1./(rair*tair))
                        ! alpha gammaloc::  Abdul-Razzak & Ghan 1998 eqn 11 . 12
                        gammaloc=(1+latvap/cpair*dqsdt)/(rhoair*qs)
                        !     growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
                        !     should depend on mean radius of mode to account for gas kinetic effects
                        gloc=1./(rhoh2o/(diff0*rhoair*qs)                                    &
                            +latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair)-1.))
                        sqrtg=sqrt(gloc)
                        ! gloc ::growth coefficent
                        beta=4.*pi*rhoh2o*gloc*gammaloc
                        etafactor2max=1.e10/(alpha*wbar)**1.5 ! this should make eta big if na is very small.
                        do m=1,nmode
                            !         internal mixture of aerosols
                            !++ag: turn on variable size
                            volc(m)=0
                            do n=1,ntype(m)
                                volc(m)=volc(m)+ma(n,m)/(rhodry(n,m)) ! only if variable size dist
                            enddo

                            if(volc(m).gt.1.e-39.and.na(m).gt.1.e-39)then
                                !ag         if(na(m).gt.1.e-39)then    !only if fixed size dist
                                etafactor2(m)=1./(na(m)*beta*sqrtg)  !fixed or variable size dist
                                !            number mode radius (m)
                                amcubeloc(m)=(3.*volc(m)/(4.*pi*exp45logsig(m)*na(m)))  ! only if variable size dist
                                smc(m)=smcrit(m) ! only for prescribed size dist

                                if(hygro(m).gt.1.e-10)then   ! loop only if variable size dist
                                    smc(m)=2.*aten*sqrt(aten/(27.*hygro(m)*amcubeloc(m))) 
                                else
                                    smc(m)=100.
                                endif
                            else
                                smc(m)=1.
                                etafactor2(m)=etafactor2max ! this should make eta big if na is very small.
                            endif
                            lnsmloc(m)=log(smc(m)) ! only if variable size dist
                        enddo
                        !--ag
                        !         single  updraft
                        wnuc=wbar
                        !        write(6,*)'uniform updraft =',wnuc
                        w=wbar
                        alw=alpha*wnuc
                        sqrtalw=sqrt(alw)
                        zeta=2.*sqrtalw*aten/(3.*sqrtg)
                        etafactor1=2.*alw*sqrtalw

                        do m=1,nmode
                            eta(m)=etafactor1*etafactor2(m)
                        enddo
                        call maxsat(zeta,eta,nmode,smc,smax)
                        lnsmax=log(smax)
                        xmincoeff=alogaten-2.*third*(lnsmax-alog2)-alog3

                        nact=0.
                        do m=1,nmode
                            x=2*(lnsmloc(m)-lnsmax)/(3*sq2*alogsig(m))
                            !original ghan code
                            !               nact=nact+0.5*(1.-ERF(x))*na(m)
                            !++ag replace sg erf with hm derf pre 1.68
                            !	        nact=nact+0.5*(1.-derf(x))*na(m)
                            !++ag 1.68 new error function
                            nact=nact+0.5*(1.-erf(x))*na(m)
                        enddo
                        nact=nact/rhoair ! convert from #/m3 to #/kg
                        !write(*,*) nact,rhoair,x,lnsmax
                        !write(*,*) lnsmloc
                        !write(*,*) alogsig
                        !write(*,*) volc
                        !write(6,*) "cldwat2.F90 activation-ln2994"
                        !call endrun 
                        return
                    end subroutine activate


                    subroutine maxsat(zeta,eta,nmode,smc,smax)
                        !      calculates maximum supersaturation for multiple
                        !      competing aerosol modes.
                        !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
                        !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.
                        implicit none
                        integer pmode
                        parameter (pmode=naer_all)
                        integer nmode ! number of modes
                        real(r8) smc(pmode) ! critical supersaturation for number mode radius
                        real(r8) zeta, eta(pmode)
                        real(r8) smax ! maximum supersaturation    ! intent out
                        integer m  ! mode index
                        real(r8) sum, g1, g2

                        do m=1,nmode
                            if(zeta.gt.1.e5*eta(m).or.smc(m)*smc(m).gt.1.e5*eta(m))then
                                !            weak forcing. essentially none activated
                                smax=1.e-20
                            else
                                !            significant activation of this mode. calc activation all modes.
                                go to 1
                            endif
                        enddo

                        return

1                       continue

                        sum=0
                        do m=1,nmode
                            if(eta(m).gt.1.e-20)then
                                g1=sqrt(zeta/eta(m))
                                g1=g1*g1*g1
                                g2=smc(m)/sqrt(eta(m)+3*zeta)
                                g2=sqrt(g2)
                                g2=g2*g2*g2
                                sum=sum+(f1(m)*g1+f2(m)*g2)/(smc(m)*smc(m))
                            else
                                sum=1.e20
                            endif
                        enddo

                        smax=1./sqrt(sum)

                        return

                    end subroutine maxsat


                end module cldwat2m
