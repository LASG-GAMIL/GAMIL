#include <misc.h>
#include <params.h>

subroutine radcswmx(lchnk   ,ncol    ,                            &
                    pint    ,pmid    ,h2ommr  ,rh      ,o3mmr   , &
                    aermmr  ,cld     ,clwp    ,rel     ,rei     , &
                    fice    ,eccf    ,coszrs  ,scon    ,solin   , &
                    asdir   ,asdif   ,aldir   ,aldif   ,nmxrgn  , &
                    pmxrgn  ,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                    fsntoac ,fsutoac ,fsnirtoa,fsnrtoac,fsnrtoaq, &
                    fsns    , fsnsc   ,fsdsc   ,fsusc   ,fsds   , &
                    sols    , soll    ,solsd   ,solld   ,fsus    )
!-----------------------------------------------------------------------
!
! Purpose:
! Solar radiation code
!
! Method:
! Basic method is Delta-Eddington as described in:
!
! Briegleb, Bruce P., 1992: Delta-Eddington
! Appoximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
!
! Four changes to the basic method described above are:
! (1) addition of sulfate aerosols (Kiehl and Briegleb, 1993)
! (2) the distinction between liquid and ice particle clouds
! (Kiehl et al, 1996);
! (3) provision for calculating TOA fluxes with spectral response to
! match Nimbus-7 visible/near-IR radiometers (Collins, 1998);
! (4) max-random overlap (Collins and Truesdale, 2000)
!
! The treatment of maximum-random overlap is described in the
! comment block "INDEX CALCULATIONS FOR MAX OVERLAP".
!
! Divides solar spectrum into 19 intervals from 0.2-5.0 micro-meters.
! solar flux fractions specified for each interval. allows for
! seasonally and diurnally varying solar input.  Includes molecular,
! cloud, aerosol, and surface scattering, along with h2o,o3,co2,o2,cloud,
! and surface absorption. Computes delta-eddington reflections and
! transmissions assuming homogeneously mixed layers. Adds the layers
! assuming scattering between layers to be isotropic, and distinguishes
! direct solar beam from scattered radiation.
!
! Longitude loops are broken into 1 or 2 sections, so that only daylight
! (i.e. coszrs > 0) computations are done.
!
! Note that an extra layer above the model top layer is added.
!
! cgs units are used.
!
! Special diagnostic calculation of the clear sky surface and total column
! absorbed flux is also done for cloud forcing diagnostics.
!
!-----------------------------------------------------------------------
    use shr_kind_mod, only: r8 => shr_kind_r8
    use ppgrid
    use history,      only: outfld ! added by DONG Li for outputting

   implicit none

   integer nspint            ! Num of spctrl intervals across solar spectrum

   parameter ( nspint = 19 )
!-----------------------Constants for new band (640-700 nm)-------------
   real(r8) v_raytau_35
   real(r8) v_raytau_64
   real(r8) v_abo3_35
   real(r8) v_abo3_64
   real(r8) v_ksa_35
   real(r8) v_ksa_64
   real(r8) v_gsa_35
   real(r8) v_gsa_64
   parameter( &
        v_raytau_35 = 0.155208, &
        v_raytau_64 = 0.0392, &
        v_abo3_35 = 2.4058030e+01, &
        v_abo3_64 = 2.210e+01, &
        v_ksa_35 = 5.64884, &
        v_ksa_64 = 3.6771, &
        v_gsa_35 = .699326, &
        v_gsa_64 = .663642 &
        )


!-------------Parameters for accelerating max-random solution-------------
!
! The solution time scales like prod(j:1->N) (1 + n_j) where
! N   = number of max-overlap regions (nmxrgn)
! n_j = number of unique cloud amounts in region j
!
! Therefore the solution cost can be reduced by decreasing n_j.
! cldmin reduces n_j by treating cloud amounts < cldmin as clear sky.
! cldeps reduces n_j by treating cloud amounts identical to log(1/cldeps)
! decimal places as identical
!
! areamin reduces the cost by dropping configurations that occupy
! a surface area < areamin of the model grid box.  The surface area
! for a configuration C(j,k_j), where j is the region number and k_j is the
! index for a unique cloud amount (in descending order from biggest to
! smallest clouds) in region j, is
!
! A = prod(j:1->N) [C(j,k_j) - C(j,k_j+1)]
!
! where C(j,0) = 1.0 and C(j,n_j+1) = 0.0.
!
! nconfgmax reduces the cost and improves load balancing by setting an upper
! bound on the number of cloud configurations in the solution.  If the number
! of configurations exceeds nconfgmax, the nconfgmax configurations with the
! largest area are retained, and the fluxes are normalized by the total area
! of these nconfgmax configurations.  For the current max/random overlap
! assumption (see subroutine cldovrlap), 30 levels, and cloud-amount
! parameterization, the mean and RMS number of configurations are
! both roughly 5.  nconfgmax has been set to the mean+2*RMS number, or 15.
!
! Minimum cloud amount (as a fraction of the grid-box area) to
! distinguish from clear sky
!
   real(r8) cldmin
   parameter (cldmin = 1.0e-80_r8)
!
! Minimimum horizontal area (as a fraction of the grid-box area) to retain
! for a unique cloud configuration in the max-random solution
!
   real(r8) areamin
   parameter (areamin = 0.01_r8)
!
! Decimal precision of cloud amount (0 -> preserve full resolution;
! 10^-n -> preserve n digits of cloud amount)
!
   real(r8) cldeps
   parameter (cldeps = 0.0_r8)
!
! Maximum number of configurations to include in solution
!
   integer nconfgmax
   parameter (nconfgmax = 15)
!------------------------------Commons----------------------------------
#include <crdcon.h>
!
! Input arguments
!
   integer, intent(in) :: lchnk             ! chunk identifier
   integer, intent(in) :: ncol              ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver) ! Level pressure
   real(r8), intent(in) :: pint(pcols,pverp) ! Interface pressure
   real(r8), intent(in) :: h2ommr(pcols,pver) ! Specific humidity (h2o mass mix ratio)
   real(r8), intent(in) :: o3mmr(pcols,pver) ! Ozone mass mixing ratio
   real(r8), intent(in) :: aermmr(pcols,pver) ! Aerosol mass mixing ratio
   real(r8), intent(in) :: rh(pcols,pver)   ! Relative humidity (fraction)
!
   real(r8), intent(in) :: cld(pcols,pver)  ! Fractional cloud cover
   real(r8), intent(in) :: clwp(pcols,pver) ! Layer liquid water path
   real(r8), intent(in) :: rel(pcols,pver)  ! Liquid effective drop size (microns)
   real(r8), intent(in) :: rei(pcols,pver)  ! Ice effective drop size (microns)
   real(r8), intent(in) :: fice(pcols,pver) ! Fractional ice content within cloud
!
   real(r8), intent(in) :: eccf             ! Eccentricity factor (1./earth-sun dist^2)
   real(r8), intent(in) :: coszrs(pcols)    ! Cosine solar zenith angle
   real(r8), intent(in) :: asdir(pcols)     ! 0.2-0.7 micro-meter srfc alb: direct rad
   real(r8), intent(in) :: aldir(pcols)     ! 0.7-5.0 micro-meter srfc alb: direct rad
   real(r8), intent(in) :: asdif(pcols)     ! 0.2-0.7 micro-meter srfc alb: diffuse rad
   real(r8), intent(in) :: aldif(pcols)     ! 0.7-5.0 micro-meter srfc alb: diffuse rad

   real(r8), intent(in) :: scon             ! solar constant
!
! IN/OUT arguments
!
   real(r8), intent(inout) :: pmxrgn(pcols,pverp) ! Maximum values of pressure for each
!                                                 !    maximally overlapped region.
!                                                 !    0->pmxrgn(i,1) is range of pressure for
!                                                 !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!                                                 !    2nd region, etc
   integer, intent(inout) ::  nmxrgn(pcols)    ! Number of maximally overlapped regions
!
! Output arguments
!

   real(r8), intent(out) :: solin(pcols)     ! Incident solar flux
   real(r8), intent(out) :: qrs(pcols,pver)  ! Solar heating rate
   real(r8), intent(out) :: fsns(pcols)      ! Surface absorbed solar flux
   real(r8), intent(out) :: fsnt(pcols)      ! Total column absorbed solar flux
   real(r8), intent(out) :: fsntoa(pcols)    ! Net solar flux at TOA
   real(r8), intent(out) :: fsds(pcols)      ! Flux shortwave downwelling surface
   real(r8), intent(out) :: fsus(pcols)      ! short wave upward flux at surface, added by DONG Li
!
   real(r8), intent(out) :: fsnsc(pcols)     ! Clear sky surface absorbed solar flux
   real(r8), intent(out) :: fsdsc(pcols)     ! Clear sky surface downwelling solar flux
   real(r8), intent(out) :: fsusc(pcols)     ! Clear sky surface upward solar flux
   real(r8), intent(out) :: fsntc(pcols)     ! Clear sky total column absorbed solar flx
   real(r8), intent(out) :: fsntoac(pcols)   ! Clear sky net solar flx at TOA
   real(r8), intent(out) :: fsutoac(pcols)   ! Clear sky upward solar flx at TOA
   real(r8), intent(out) :: sols(pcols)      ! Direct solar rad on surface (< 0.7)
   real(r8), intent(out) :: soll(pcols)      ! Direct solar rad on surface (>= 0.7)
   real(r8), intent(out) :: solsd(pcols)     ! Diffuse solar rad on surface (< 0.7)
   real(r8), intent(out) :: solld(pcols)     ! Diffuse solar rad on surface (>= 0.7)
   real(r8), intent(out) :: fsnirtoa(pcols)  ! Near-IR flux absorbed at toa
   real(r8), intent(out) :: fsnrtoac(pcols)  ! Clear sky near-IR flux absorbed at toa
   real(r8), intent(out) :: fsnrtoaq(pcols)  ! Net near-IR flux at toa >= 0.7 microns
!
!---------------------------Local variables-----------------------------
!
! Max/random overlap variables
!
   real(r8) asort(pverp)     ! 1 - cloud amounts to be sorted for max ovrlp.
   real(r8) atmp             ! Temporary storage for sort when nxs = 2
   real(r8) cld0             ! 1 - (cld amt) used to make wstr, cstr, nstr
   real(r8) totwgt           ! Total of xwgts = total fractional area of
!   grid-box covered by cloud configurations
!   included in solution to fluxes

   real(r8) wgtv(nconfgmax)  ! Weights for fluxes
!   1st index is configuration number
   real(r8) wstr(pverp,pverp) ! area weighting factors for streams
!   1st index is for stream #,
!   2nd index is for region #

   real(r8) xexpt            ! solar direct beam trans. for layer above
   real(r8) xrdnd            ! diffuse reflectivity for layer above
   real(r8) xrupd            ! diffuse reflectivity for layer below
   real(r8) xrups            ! direct-beam reflectivity for layer below
   real(r8) xtdnt            ! total trans for layers above

   real(r8) xwgt             ! product of cloud amounts

   real(r8) yexpt            ! solar direct beam trans. for layer above
   real(r8) yrdnd            ! diffuse reflectivity for layer above
   real(r8) yrupd            ! diffuse reflectivity for layer below
   real(r8) ytdnd            ! dif-beam transmission for layers above
   real(r8) ytupd            ! dif-beam transmission for layers below

   real(r8) zexpt            ! solar direct beam trans. for layer above
   real(r8) zrdnd            ! diffuse reflectivity for layer above
   real(r8) zrupd            ! diffuse reflectivity for layer below
   real(r8) zrups            ! direct-beam reflectivity for layer below
   real(r8) ztdnt            ! total trans for layers above

   logical new_term          ! Flag for configurations to include in fluxes
   logical region_found      ! flag for identifying regions

   integer ccon(0:pverp,nconfgmax)
! flags for presence of clouds
!   1st index is for level # (including
!    layer above top of model and at surface)
!   2nd index is for configuration #
   integer cstr(0:pverp,pverp)
! flags for presence of clouds
!   1st index is for level # (including
!    layer above top of model and at surface)
!   2nd index is for stream #
   integer icond(0:pverp,nconfgmax)
! Indices for copying rad. properties from
!     one identical downward cld config.
!     to another in adding method (step 2)
!   1st index is for interface # (including
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer iconu(0:pverp,nconfgmax)
! Indices for copying rad. properties from
!     one identical upward configuration
!     to another in adding method (step 2)
!   1st index is for interface # (including
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer iconfig           ! Counter for random-ovrlap configurations
   integer irgn              ! Index for max-overlap regions
   integer is0               ! Lower end of stream index range
   integer is1               ! Upper end of stream index range
   integer isn               ! Stream index
   integer istr(pverp+1)     ! index for stream #s during flux calculation
   integer istrtd(0:pverp,0:nconfgmax+1)
! indices into icond
!   1st index is for interface # (including
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer istrtu(0:pverp,0:nconfgmax+1)
! indices into iconu
!   1st index is for interface # (including
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer j                 ! Configuration index
   integer k1                ! Level index
   integer k2                ! Level index
   integer ksort(pverp)      ! Level indices of cloud amounts to be sorted
   integer ktmp              ! Temporary storage for sort when nxs = 2
   integer kx1(0:pverp)      ! Level index for top of max-overlap region
   integer kx2(0:pverp)      ! Level index for bottom of max-overlap region
   integer l                 ! Index
   integer l0                ! Index
   integer mrgn              ! Counter for nrgn
   integer mstr              ! Counter for nstr
   integer n0                ! Number of configurations with ccon(k,:)==0
   integer n1                ! Number of configurations with ccon(k,:)==1
   integer nconfig           ! Number of random-ovrlap configurations
   integer nconfigm          ! Value of config before testing for areamin,
!    nconfgmax
   integer npasses           ! number of passes over the indexing loop
   integer nrgn              ! Number of max overlap regions at current
!    longitude
   integer nstr(pverp)       ! Number of unique cloud configurations
!   ("streams") in a max-overlapped region
!   1st index is for region #
   integer nuniq             ! # of unique cloud configurations
   integer nuniqd(0:pverp)   ! # of unique cloud configurations: TOA
!   to level k
   integer nuniqu(0:pverp)   ! # of unique cloud configurations: surface
!   to level k
   integer nxs               ! Number of cloudy layers between k1 and k2
   integer ptr0(nconfgmax)   ! Indices of configurations with ccon(k,:)==0
   integer ptr1(nconfgmax)   ! Indices of configurations with ccon(k,:)==1
   integer ptrc(nconfgmax)   ! Pointer for configurations sorted by wgtv
   integer findvalue         ! Function for finding kth smallest element
!   in a vector
   external findvalue

!
! Other
!
   integer ns                ! Spectral loop index
   integer i                 ! Longitude loop index
   integer k                 ! Level loop index
   integer km1               ! k - 1
   integer kp1               ! k + 1
   integer n                 ! Loop index for daylight
   integer ndayc             ! Number of daylight columns
   integer idayc(pcols)      ! Daytime column indices
   integer indxsl            ! Index for cloud particle properties
!
! A. Slingo's data for cloud particle radiative properties (from 'A GCM
! Parameterization for the Shortwave Properties of Water Clouds' JAS
! vol. 46 may 1989 pp 1419-1427)
!
   real(r8) abarl(4)         ! A coefficient for extinction optical depth
   real(r8) bbarl(4)         ! B coefficient for extinction optical depth
   real(r8) cbarl(4)         ! C coefficient for single scat albedo
   real(r8) dbarl(4)         ! D coefficient for single  scat albedo
   real(r8) ebarl(4)         ! E coefficient for asymmetry parameter
   real(r8) fbarl(4)         ! F coefficient for asymmetry parameter

   save abarl, bbarl, cbarl, dbarl, ebarl, fbarl

   data abarl/ 2.817e-02, 2.682e-02,2.264e-02,1.281e-02/
   data bbarl/ 1.305    , 1.346    ,1.454    ,1.641    /
   data cbarl/-5.62e-08 ,-6.94e-06 ,4.64e-04 ,0.201    /
   data dbarl/ 1.63e-07 , 2.35e-05 ,1.24e-03 ,7.56e-03 /
   data ebarl/ 0.829    , 0.794    ,0.754    ,0.826    /
   data fbarl/ 2.482e-03, 4.226e-03,6.560e-03,4.353e-03/

   real(r8) abarli           ! A coefficient for current spectral band
   real(r8) bbarli           ! B coefficient for current spectral band
   real(r8) cbarli           ! C coefficient for current spectral band
   real(r8) dbarli           ! D coefficient for current spectral band
   real(r8) ebarli           ! E coefficient for current spectral band
   real(r8) fbarli           ! F coefficient for current spectral band
!
! Caution... A. Slingo recommends no less than 4.0 micro-meters nor
! greater than 20 micro-meters
!
! ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
!
   real(r8) abari(4)         ! a coefficient for extinction optical depth
   real(r8) bbari(4)         ! b coefficient for extinction optical depth
   real(r8) cbari(4)         ! c coefficient for single scat albedo
   real(r8) dbari(4)         ! d coefficient for single scat albedo
   real(r8) ebari(4)         ! e coefficient for asymmetry parameter
   real(r8) fbari(4)         ! f coefficient for asymmetry parameter

   save abari, bbari, cbari, dbari, ebari, fbari

   data abari/ 3.448e-03, 3.448e-03,3.448e-03,3.448e-03/
   data bbari/ 2.431    , 2.431    ,2.431    ,2.431    /
   data cbari/ 1.00e-05 , 1.10e-04 ,1.861e-02,.46658   /
   data dbari/ 0.0      , 1.405e-05,8.328e-04,2.05e-05 /
   data ebari/ 0.7661   , 0.7730   ,0.794    ,0.9595   /
   data fbari/ 5.851e-04, 5.665e-04,7.267e-04,1.076e-04/

   real(r8) abarii           ! A coefficient for current spectral band
   real(r8) bbarii           ! B coefficient for current spectral band
   real(r8) cbarii           ! C coefficient for current spectral band
   real(r8) dbarii           ! D coefficient for current spectral band
   real(r8) ebarii           ! E coefficient for current spectral band
   real(r8) fbarii           ! F coefficient for current spectral band
!
   real(r8) delta            ! Pressure (in atm) for stratos. h2o limit
   real(r8) o2mmr            ! O2 mass mixing ratio:

   save delta, o2mmr

   data delta /  1.70e-3 /
   data o2mmr / .23143 /

   real(r8) albdir(pcols,nspint) ! Current spc intrvl srf alb to direct rad
   real(r8) albdif(pcols,nspint) ! Current spc intrvl srf alb to diffuse rad
!
! Next series depends on spectral interval
!
   real(r8) frcsol(nspint)   ! Fraction of solar flux in spectral interval
   real(r8) wavmin(nspint)   ! Min wavelength (micro-meters) of interval
   real(r8) wavmax(nspint)   ! Max wavelength (micro-meters) of interval
   real(r8) raytau(nspint)   ! Rayleigh scattering optical depth
   real(r8) abh2o(nspint)    ! Absorption coefficiant for h2o (cm2/g)
   real(r8) abo3 (nspint)    ! Absorption coefficiant for o3  (cm2/g)
   real(r8) abco2(nspint)    ! Absorption coefficiant for co2 (cm2/g)
   real(r8) abo2 (nspint)    ! Absorption coefficiant for o2  (cm2/g)
   real(r8) ph2o(nspint)     ! Weight of h2o in spectral interval
   real(r8) pco2(nspint)     ! Weight of co2 in spectral interval
   real(r8) po2 (nspint)     ! Weight of o2  in spectral interval
   real(r8) nirwgt(nspint)   ! Spectral Weights to simulate Nimbus-7 filter
   real(r8) wgtint           ! Weight for specific spectral interval

   save frcsol ,wavmin ,wavmax ,raytau ,abh2o ,abo3 , &
        abco2  ,abo2   ,ph2o   ,pco2   ,po2   ,nirwgt

   data frcsol / .001488, .001389, .001290, .001686, .002877, &
                 .003869, .026336, .360739, .065392, .526861, &
                 .526861, .526861, .526861, .526861, .526861, &
                 .526861, .006239, .001834, .001834/
!
! weight for 0.64 - 0.7 microns  appropriate to clear skies over oceans
!
   data nirwgt /  0.0,   0.0,   0.0,      0.0,   0.0, &
                  0.0,   0.0,   0.0, 0.320518,   1.0,  1.0, &
                  1.0,   1.0,   1.0,      1.0,   1.0, &
                  1.0,   1.0,   1.0 /

   data wavmin / .200,  .245,  .265,  .275,  .285, &
                 .295,  .305,  .350,  .640,  .700,  .701, &
                 .701,  .701,  .701,  .702,  .702, &
                 2.630, 4.160, 4.160/

   data wavmax / .245,  .265,  .275,  .285,  .295, &
                 .305,  .350,  .640,  .700, 5.000, 5.000, &
                 5.000, 5.000, 5.000, 5.000, 5.000, &
                 2.860, 4.550, 4.550/

   data raytau / 4.020, 2.180, 1.700, 1.450, 1.250, &
                 1.085, 0.730, v_raytau_35, v_raytau_64, 0.020, &
                 .0001, .0001, .0001, .0001, .0001, .0001, &
                 .0001, .0001, .0001/
!
! Absorption coefficients
!
   data abh2o /    .000,     .000,    .000,    .000,    .000, &
                   .000,     .000,    .000,    .000,    .002, &
                   .035,     .377,   1.950,   9.400,  44.600, &
                   190.000,     .000,    .000,    .000/

   data abo3  /5.370e+04, 13.080e+04,  9.292e+04, 4.530e+04, 1.616e+04, &
               4.441e+03,  1.775e+02, v_abo3_35, v_abo3_64,      .000, &
               .000,   .000    ,   .000   ,   .000   ,      .000, &
               .000,   .000    ,   .000   ,   .000    /

   data abco2  /   .000,     .000,    .000,    .000,    .000, &
                   .000,     .000,    .000,    .000,    .000, &
                   .000,     .000,    .000,    .000,    .000, &
                   .000,     .094,    .196,   1.963/

   data abo2  /    .000,     .000,    .000,    .000,    .000, &
                   .000,     .000,    .000,1.11e-05,6.69e-05, &
                   .000,     .000,    .000,    .000,    .000, &
                   .000,     .000,    .000,    .000/
!
! Spectral interval weights
!
   data ph2o  /    .000,     .000,    .000,    .000,    .000, &
        .000,     .000,    .000,    .000,    .505,     &
        .210,     .120,    .070,    .048,    .029,     &
        .018,     .000,    .000,    .000/

   data pco2  /    .000,     .000,    .000,    .000,    .000, &
        .000,     .000,    .000,    .000,    .000,     &
        .000,     .000,    .000,    .000,    .000,     &
        .000,    1.000,    .640,    .360/

   data po2   /    .000,     .000,    .000,    .000,    .000, &
        .000,     .000,    .000,   1.000,   1.000,     &
        .000,     .000,    .000,    .000,    .000,     &
        .000,     .000,    .000,    .000/
!
! Diagnostic and accumulation arrays; note that sfltot, fswup, and
! fswdn are not used in the computation,but are retained for future use.
!
   real(r8) solflx           ! Solar flux in current interval
   real(r8) sfltot           ! Spectrally summed total solar flux
   real(r8) totfld(0:pver)   ! Spectrally summed flux divergence
   real(r8) fswup(0:pverp)   ! Spectrally summed up flux
   real(r8) fswdn(0:pverp)   ! Spectrally summed down flux
!
! Cloud radiative property arrays
!
   real(r8) tauxcl(pcols,0:pver) ! water cloud extinction optical depth
   real(r8) tauxci(pcols,0:pver) ! ice cloud extinction optical depth
   real(r8) wcl(pcols,0:pver) ! liquid cloud single scattering albedo
   real(r8) gcl(pcols,0:pver) ! liquid cloud asymmetry parameter
   real(r8) fcl(pcols,0:pver) ! liquid cloud forward scattered fraction
   real(r8) wci(pcols,0:pver) ! ice cloud single scattering albedo
   real(r8) gci(pcols,0:pver) ! ice cloud asymmetry parameter
   real(r8) fci(pcols,0:pver) ! ice cloud forward scattered fraction
!
! Aerosol radiative property arrays
!
   real(r8) tauxar(pcols,0:pver) ! aerosol extinction optical depth
   real(r8) wa(pcols,0:pver) ! aerosol single scattering albedo
   real(r8) ga(pcols,0:pver) ! aerosol assymetry parameter
   real(r8) fa(pcols,0:pver) ! aerosol forward scattered fraction
!
! Sulphate aerosol properties taken from:
!
! Kiehl, J.T., B.P.Briegleb, 1993. The Relative Roles of Sulfate Aerosols
! and Greenhouse Gases in Climate Forcing. Science, Vol. 260, pp. 311-314.
!
   real(r8) ksa(nspint)      ! aerosol spectral mass abs. coeff(m2/g)
   real(r8) wsa(nspint)      ! aerosol spectral single scat. albedo
   real(r8) gsa(nspint)      ! aerosol spectral asymmetry parameter
!
   data ksa /11.1163, 10.5472, 10.2468, 10.0392,  9.8292, &
        9.6199,  9.0407,v_ksa_35,v_ksa_64,  1.9169,   &
        0.3780,  0.3780,  0.3780,  0.3780,  0.5704,   &
        0.5704,  0.5704,  0.5704,  0.5704 /

   data wsa / .999999, .999999, .999999, .999999, .999999, &
        .999999, .999999, .999999, .999999, .999991,  &
        .989772, .989772, .989772, .989772, .847061,  &
        .847061, .847061, .847061, .847061 /

   data gsa / .719161, .719012, .718453, .717820, .716997, &
        .715974, .712743,v_gsa_35,v_gsa_64, .618115,  &
        .485286, .485286, .485286, .485286, .295557,  &
        .295557, .295557, .295557, .295557 /

!
! Other variables and arrays needed for aerosol:
!
   real(r8) rhfac            ! multiplication factor for kaer
   real(r8) rhpc             ! level relative humidity in %

   real(r8) a0               ! constant in rh mult factor
   real(r8) a1               ! constant in rh mult factor
   real(r8) a2               ! constant in rh mult factor
   real(r8) a3               ! constant in rh mult factor

   save a0,a1,a2,a3

   data a0 / -9.2906106183    /
   data a1 /  0.52570211505   /
   data a2 / -0.0089285760691 /
   data a3 /  5.0877212432e-05/

!
! Various arrays and other constants:
!
   real(r8) pflx(pcols,0:pverp) ! Interface press, including extra layer
   real(r8) zenfac(pcols)    ! Square root of cos solar zenith angle
   real(r8) sqrco2           ! Square root of the co2 mass mixg ratio
   real(r8) tmp1             ! Temporary constant array
   real(r8) tmp2             ! Temporary constant array
   real(r8) pdel             ! Pressure difference across layer
   real(r8) path             ! Mass path of layer
   real(r8) ptop             ! Lower interface pressure of extra layer
   real(r8) ptho2            ! Used to compute mass path of o2
   real(r8) ptho3            ! Used to compute mass path of o3
   real(r8) pthco2           ! Used to compute mass path of co2
   real(r8) pthh2o           ! Used to compute mass path of h2o
   real(r8) h2ostr           ! Inverse sq. root h2o mass mixing ratio
   real(r8) wavmid(nspint)   ! Spectral interval middle wavelength
   real(r8) trayoslp         ! Rayleigh optical depth/standard pressure
   real(r8) tmp1l            ! Temporary constant array
   real(r8) tmp2l            ! Temporary constant array
   real(r8) tmp3l            ! Temporary constant array
   real(r8) tmp1i            ! Temporary constant array
   real(r8) tmp2i            ! Temporary constant array
   real(r8) tmp3i            ! Temporary constant array
   real(r8) rdenom           ! Multiple scattering term
   real(r8) rdirexp          ! layer direct ref times exp transmission
   real(r8) tdnmexp          ! total transmission - exp transmission
   real(r8) psf(nspint)      ! Frac of solar flux in spect interval
!
! Layer absorber amounts; note that 0 refers to the extra layer added
! above the top model layer
!
   real(r8) uh2o(pcols,0:pver) ! Layer absorber amount of h2o
   real(r8) uo3(pcols,0:pver) ! Layer absorber amount of  o3
   real(r8) uco2(pcols,0:pver) ! Layer absorber amount of co2
   real(r8) uo2(pcols,0:pver) ! Layer absorber amount of  o2
   real(r8) uaer(pcols,0:pver) ! Layer aerosol amount
!
! Total column absorber amounts:
!
   real(r8) uth2o(pcols)     ! Total column  absorber amount of  h2o
   real(r8) uto3(pcols)      ! Total column  absorber amount of  o3
   real(r8) utco2(pcols)     ! Total column  absorber amount of  co2
   real(r8) uto2(pcols)      ! Total column  absorber amount of  o2
!
! These arrays are defined for pver model layers; 0 refers to the extra
! layer on top:
!
   real(r8) rdir(nspint,pcols,0:pver) ! Layer reflectivity to direct rad
   real(r8) rdif(nspint,pcols,0:pver) ! Layer reflectivity to diffuse rad
   real(r8) tdir(nspint,pcols,0:pver) ! Layer transmission to direct rad
   real(r8) tdif(nspint,pcols,0:pver) ! Layer transmission to diffuse rad
   real(r8) explay(nspint,pcols,0:pver) ! Solar beam exp trans. for layer

   real(r8) rdirc(nspint,pcols,0:pver) ! Clear Layer reflec. to direct rad
   real(r8) rdifc(nspint,pcols,0:pver) ! Clear Layer reflec. to diffuse rad
   real(r8) tdirc(nspint,pcols,0:pver) ! Clear Layer trans. to direct rad
   real(r8) tdifc(nspint,pcols,0:pver) ! Clear Layer trans. to diffuse rad
   real(r8) explayc(nspint,pcols,0:pver) ! Solar beam exp trans. clear layer

   real(r8) flxdiv           ! Flux divergence for layer
!
!
! Radiative Properties:
!
! There are 1 classes of properties:
! (1. All-sky bulk properties
! (2. Clear-sky properties
!
! The first set of properties are generated during step 2 of the solution.
!
! These arrays are defined at model interfaces; in 1st index (for level #),
! 0 is the top of the extra layer above the model top, and
! pverp is the earth surface.  2nd index is for cloud configuration
! defined over a whole column.
!
   real(r8) exptdn(0:pverp,nconfgmax) ! Sol. beam trans from layers above
   real(r8) rdndif(0:pverp,nconfgmax) ! Ref to dif rad for layers above
   real(r8) rupdif(0:pverp,nconfgmax) ! Ref to dif rad for layers below
   real(r8) rupdir(0:pverp,nconfgmax) ! Ref to dir rad for layers below
   real(r8) tdntot(0:pverp,nconfgmax) ! Total trans for layers above
!
! Bulk properties used during the clear-sky calculation.
!
   real(r8) exptdnc(0:pverp) ! clr: Sol. beam trans from layers above
   real(r8) rdndifc(0:pverp) ! clr: Ref to dif rad for layers above
   real(r8) rupdifc(0:pverp) ! clr: Ref to dif rad for layers below
   real(r8) rupdirc(0:pverp) ! clr: Ref to dir rad for layers below
   real(r8) tdntotc(0:pverp) ! clr: Total trans for layers above

   real(r8) fluxup(0:pverp)  ! Up   flux at model interface
   real(r8) fluxdn(0:pverp)  ! Down flux at model interface
   real(r8) wexptdn          ! Direct solar beam trans. to surface

    real(r8) :: swuptoa(pcols) = 0.0_r8  ! short wave up at top of atmosphere
    real(r8) :: swdntoa(pcols) = 0.0_r8  ! short wave down at top of atmosphere

!-----------------------------------------------------------------------
! START OF CALCULATION
!-----------------------------------------------------------------------
!
   do i=1, ncol
!
! Initialize output fields:
!
      fsds(i)     = 0.0_r8

      fsnirtoa(i) = 0.0_r8
      fsnrtoac(i) = 0.0_r8
      fsnrtoaq(i) = 0.0_r8

      fsns(i)     = 0.0_r8
      fsnsc(i)    = 0.0_r8
      fsdsc(i)    = 0.0_r8
      fsusc(i)    = 0.0_r8

      fsnt(i)     = 0.0_r8
      fsntc(i)    = 0.0_r8
      fsntoa(i)   = 0.0_r8
      fsntoac(i)  = 0.0_r8
      fsutoac(i)  = 0.0_r8

      solin(i)    = 0.0_r8

      sols(i)     = 0.0_r8
      soll(i)     = 0.0_r8
      solsd(i)    = 0.0_r8
      solld(i)    = 0.0_r8

      do k=1, pver
         qrs(i,k) = 0.0_r8
      end do
   end do
!
! Compute starting, ending daytime loop indices:
!  *** Note this logic assumes day and night points are contiguous so
!  *** will not work in general with chunked data structure.
!
   ndayc = 0
   do i=1,ncol
      if (coszrs(i) > 0.0_r8) then
         ndayc = ndayc + 1
         idayc(ndayc) = i
      end if
   end do
!
! If night everywhere, return:
!
   if (ndayc == 0) return
!
! Perform other initializations
!
   tmp1   = 0.5_r8/(gravit*sslp)
   tmp2   = delta/gravit
   sqrco2 = sqrt(co2mmr)

   do n=1,ndayc
      i=idayc(n)
!
! Define solar incident radiation and interface pressures:
!
         solin(i)  = scon*eccf*coszrs(i)
         pflx(i,0) = 0._r8
         do k=1,pverp
            pflx(i,k) = pint(i,k)
         end do
!
! Compute optical paths:
!
         ptop      = pflx(i,1)
         ptho2     = o2mmr * ptop / gravit
         ptho3     = o3mmr(i,1) * ptop / gravit
         pthco2    = sqrco2 * (ptop / gravit)
         h2ostr    = sqrt( 1._r8 / h2ommr(i,1) )
         zenfac(i) = sqrt(coszrs(i))
         pthh2o    = ptop**2*tmp1 + (ptop*rga)* &
                    (h2ostr*zenfac(i)*delta)
         uh2o(i,0) = h2ommr(i,1)*pthh2o
         uco2(i,0) = zenfac(i)*pthco2
         uo2 (i,0) = zenfac(i)*ptho2
         uo3 (i,0) = ptho3
         uaer(i,0) = 0.0_r8

         do k=1,pver
            pdel      = pflx(i,k+1) - pflx(i,k)
            path      = pdel / gravit
            ptho2     = o2mmr * path
            ptho3     = o3mmr(i,k) * path
            pthco2    = sqrco2 * path
            h2ostr    = sqrt(1.0_r8/h2ommr(i,k))
            pthh2o    = (pflx(i,k+1)**2 - pflx(i,k)**2)*tmp1 + pdel*h2ostr*zenfac(i)*tmp2
            uh2o(i,k) = h2ommr(i,k)*pthh2o
            uco2(i,k) = zenfac(i)*pthco2
            uo2 (i,k) = zenfac(i)*ptho2
            uo3 (i,k) = ptho3
!
! Adjust aerosol amount by relative humidity factor:

            if (rh(i,k) > .90) then
               rhfac  = 2.8
            else if (rh(i,k) .lt. .60 ) then
               rhfac  = 1.0
            else
               rhpc   = 100. * rh(i,k)
               rhfac  = (a0 + a1*rhpc + a2*rhpc**2 + a3*rhpc**3)
            endif
            uaer(i,k) = aermmr(i,k)*rhfac*path
         end do
!
! Compute column absorber amounts for the clear sky computation:
!
         uth2o(i) = 0.0_r8
         uto3(i)  = 0.0_r8
         utco2(i) = 0.0_r8
         uto2(i)  = 0.0_r8

         do k=1,pver
            uth2o(i) = uth2o(i) + uh2o(i,k)
            uto3(i)  = uto3(i)  + uo3(i,k)
            utco2(i) = utco2(i) + uco2(i,k)
            uto2(i)  = uto2(i)  + uo2(i,k)
         end do
!
! Set cloud properties for top (0) layer; so long as tauxcl is zero,
! there is no cloud above top of model; the other cloud properties
! are arbitrary:
!
         tauxcl(i,0)  = 0._r8
         wcl(i,0)     = 0.999999_r8
         gcl(i,0)     = 0.85_r8
         fcl(i,0)     = 0.725_r8
         tauxci(i,0)  = 0._r8
         wci(i,0)     = 0.999999_r8
         gci(i,0)     = 0.85_r8
         fci(i,0)     = 0.725_r8
!
! Aerosol
!
         tauxar(i,0)  = 0._r8
         wa(i,0)      = 0.925_r8
         ga(i,0)      = 0.850_r8
         fa(i,0)      = 0.7225_r8
!
! End  do n=1,ndayc
!
   end do
!
! Begin spectral loop
!
   do ns=1,nspint
!
! Set index for cloud particle properties based on the wavelength,
! according to A. Slingo (1989) equations 1-3:
! Use index 1 (0.25 to 0.69 micrometers) for visible
! Use index 2 (0.69 - 1.19 micrometers) for near-infrared
! Use index 3 (1.19 to 2.38 micrometers) for near-infrared
! Use index 4 (2.38 to 4.00 micrometers) for near-infrared
!
! Note that the minimum wavelength is encoded (with .001, .002, .003)
! in order to specify the index appropriate for the near-infrared
! cloud absorption properties
!
      if(wavmax(ns) <= 0.7_r8) then
         indxsl = 1
      else if(wavmin(ns) == 0.700_r8) then
         indxsl = 2
      else if(wavmin(ns) == 0.701_r8) then
         indxsl = 3
      else if(wavmin(ns) == 0.702_r8 .or. wavmin(ns) > 2.38_r8) then
         indxsl = 4
      end if
!
! Set cloud extinction optical depth, single scatter albedo,
! asymmetry parameter, and forward scattered fraction:
!
      abarli = abarl(indxsl)
      bbarli = bbarl(indxsl)
      cbarli = cbarl(indxsl)
      dbarli = dbarl(indxsl)
      ebarli = ebarl(indxsl)
      fbarli = fbarl(indxsl)
!
      abarii = abari(indxsl)
      bbarii = bbari(indxsl)
      cbarii = cbari(indxsl)
      dbarii = dbari(indxsl)
      ebarii = ebari(indxsl)
      fbarii = fbari(indxsl)
!
! adjustfraction within spectral interval to allow for the possibility of
! sub-divisions within a particular interval:
!
      psf(ns) = 1.0_r8
      if(ph2o(ns)/=0._r8) psf(ns) = psf(ns)*ph2o(ns)
      if(pco2(ns)/=0._r8) psf(ns) = psf(ns)*pco2(ns)
      if(po2 (ns)/=0._r8) psf(ns) = psf(ns)*po2 (ns)

      do n=1,ndayc
         i=idayc(n)
            do k=1,pver
!
! liquid
!
               tmp1l = abarli + bbarli/rel(i,k)
               tmp2l = 1._r8 - cbarli - dbarli*rel(i,k)
               tmp3l = fbarli*rel(i,k)
!
! ice
!
               tmp1i = abarii + bbarii/rei(i,k)
               tmp2i = 1._r8 - cbarii - dbarii*rei(i,k)
               tmp3i = fbarii*rei(i,k)

               if (cld(i,k) >= cldmin .and. cld(i,k) >= cldeps) then
                  tauxcl(i,k) = clwp(i,k)*tmp1l*(1._r8-fice(i,k))
                  tauxci(i,k) = clwp(i,k)*tmp1i*fice(i,k)
               else
                  tauxcl(i,k) = 0.0
                  tauxci(i,k) = 0.0
               endif
!
! Do not let single scatter albedo be 1.  Delta-eddington solution
! for non-conservative case has different analytic form from solution
! for conservative case, and raddedmx is written for non-conservative case.
!
               wcl(i,k) = min(tmp2l,.999999_r8)
               gcl(i,k) = ebarli + tmp3l
               fcl(i,k) = gcl(i,k)*gcl(i,k)
!
               wci(i,k) = min(tmp2i,.999999_r8)
               gci(i,k) = ebarii + tmp3i
               fci(i,k) = gci(i,k)*gci(i,k)
!
! Set aerosol properties
! Conversion factor to adjust aerosol extinction (m2/g)
!
               tauxar(i,k) = 1.e4 * ksa(ns) * uaer(i,k)
!
               wa(i,k)     = wsa(ns)
               ga(i,k)     = gsa(ns)
               fa(i,k)     = gsa(ns)*gsa(ns)
!
! End do k=1,pver
!
            end do
!
! End do n=1,ndayc
!
      end do

!
! Set reflectivities for surface based on mid-point wavelength
!
      wavmid(ns) = 0.5_r8*(wavmin(ns) + wavmax(ns))
!
! Wavelength less  than 0.7 micro-meter
!
      if (wavmid(ns) < 0.7_r8 ) then
         do n=1,ndayc
            i=idayc(n)
               albdir(i,ns) = asdir(i)
               albdif(i,ns) = asdif(i)
         end do
!
! Wavelength greater than 0.7 micro-meter
!
      else
         do n=1,ndayc
            i=idayc(n)
               albdir(i,ns) = aldir(i)
               albdif(i,ns) = aldif(i)
         end do
      end if
      trayoslp = raytau(ns)/sslp
!
! Layer input properties now completely specified; compute the
! delta-Eddington solution reflectivities and transmissivities
! for each layer
!
      call raddedmx(coszrs   ,ndayc    ,idayc   , &
              abh2o(ns),abo3(ns) ,abco2(ns),abo2(ns) , &
              uh2o     ,uo3      ,uco2     ,uo2      , &
              trayoslp ,pflx     ,ns       , &
              tauxcl   ,wcl      ,gcl      ,fcl      , &
              tauxci   ,wci      ,gci      ,fci      , &
              tauxar   ,wa       ,ga       ,fa       , &
              rdir     ,rdif     ,tdir     ,tdif     ,explay  , &
              rdirc    ,rdifc    ,tdirc    ,tdifc    ,explayc )
!
! End spectral loop
!
   end do
!
!----------------------------------------------------------------------
!
! Solution for max/random cloud overlap.
!
! Steps:
! (1. delta-Eddington solution for each layer (called above)
!
! (2. The adding method is used to
! compute the reflectivity and transmissivity to direct and diffuse
! radiation from the top and bottom of the atmosphere for each
! cloud configuration.  This calculation is based upon the
! max-random overlap assumption.
!
! (3. to solve for the fluxes, combine the
! bulk properties of the atmosphere above/below the region.
!
! Index calculations for steps 2-3 are performed outside spectral
! loop to avoid redundant calculations.  Index calculations (with
! application of areamin & nconfgmax conditions) are performed
! first to identify the minimum subset of terms for the configurations
! satisfying the areamin & nconfgmax conditions. This minimum set is
! used to identify the corresponding minimum subset of terms in
! steps 2 and 3.
!

   do n=1,ndayc
      i=idayc(n)

!----------------------------------------------------------------------
! INDEX CALCULATIONS FOR MAX OVERLAP
!
! The column is divided into sets of adjacent layers, called regions,
! in which the clouds are maximally overlapped.  The clouds are
! randomly overlapped between different regions.  The number of
! regions in a column is set by nmxrgn, and the range of pressures
! included in each region is set by pmxrgn.
!
! The following calculations determine the number of unique cloud
! configurations (assuming maximum overlap), called "streams",
! within each region. Each stream consists of a vector of binary
! clouds (either 0 or 100% cloud cover).  Over the depth of the region,
! each stream requires a separate calculation of radiative properties. These
! properties are generated using the adding method from
! the radiative properties for each layer calculated by raddedmx.
!
! The upward and downward-propagating streams are treated
! separately.
!
! We will refer to a particular configuration of binary clouds
! within a single max-overlapped region as a "stream".  We will
! refer to a particular arrangement of binary clouds over the entire column
! as a "configuration".
!
! This section of the code generates the following information:
! (1. nrgn    : the true number of max-overlap regions (need not = nmxrgn)
! (2. nstr    : the number of streams in a region (>=1)
! (3. cstr    : flags for presence of clouds at each layer in each stream
! (4. wstr    : the fractional horizontal area of a grid box covered
! by each stream
! (5. kx1,2   : level indices for top/bottom of each region
!
! The max-overlap calculation proceeds in 3 stages:
! (1. compute layer radiative properties in raddedmx.
! (2. combine these properties between layers
! (3. combine properties to compute fluxes at each interface.
!
! Most of the indexing information calculated here is used in steps 2-3
! after the call to raddedmx.
!
! Initialize indices for layers to be max-overlapped
!
! Loop to handle fix in totwgt=0. For original overlap config
! from npasses = 0.
!
         npasses = 0
         do
            do irgn = 0, nmxrgn(i)
               kx2(irgn) = 0
            end do
            mrgn = 0
!
! Outermost loop over regions (sets of adjacent layers) to be max overlapped
!
            do irgn = 1, nmxrgn(i)
!
! Calculate min/max layer indices inside region.
!
               region_found = .false.
               if (kx2(irgn-1) < pver) then
                  k1 = kx2(irgn-1)+1
                  kx1(irgn) = k1
                  kx2(irgn) = k1-1
                  do k2 = pver, k1, -1
                     if (pmid(i,k2) <= pmxrgn(i,irgn)) then
                        kx2(irgn) = k2
                        mrgn = mrgn+1
                        region_found = .true.
                        exit
                     end if
                  end do
               else
                  exit
               endif

               if (region_found) then
!
! Sort cloud areas and corresponding level indices.
!
                  nxs = 0
                  if (cldeps > 0) then
                     do k = k1,k2
                        if (cld(i,k) >= cldmin .and. cld(i,k) >= cldeps) then
                           nxs = nxs+1
                           ksort(nxs) = k
!
! We need indices for clouds in order of largest to smallest, so
! sort 1-cld in ascending order
!
                           asort(nxs) = 1.0_r8-(floor(cld(i,k)/cldeps)*cldeps)
                        end if
                     end do
                  else
                     do k = k1,k2
                        if (cld(i,k) >= cldmin) then
                           nxs = nxs+1
                           ksort(nxs) = k
!
! We need indices for clouds in order of largest to smallest, so
! sort 1-cld in ascending order
!
                           asort(nxs) = 1.0_r8-cld(i,k)
                        end if
                     end do
                  endif
!
! If nxs eq 1, no need to sort.
! If nxs eq 2, sort by swapping if necessary
! If nxs ge 3, sort using local sort routine
!
                  if (nxs == 2) then
                     if (asort(2) < asort(1)) then
                        ktmp = ksort(1)
                        ksort(1) = ksort(2)
                        ksort(2) = ktmp

                        atmp = asort(1)
                        asort(1) = asort(2)
                        asort(2) = atmp
                     endif
                  else if (nxs >= 3) then
                     call sortarray(nxs,asort,ksort)
                  endif
!
! Construct wstr, cstr, nstr for this region
!
                  cstr(k1:k2,1:nxs+1) = 0
                  mstr = 1
                  cld0 = 0.0_r8
                  do l = 1, nxs
                     if (asort(l) /= cld0) then
                        wstr(mstr,mrgn) = asort(l) - cld0
                        cld0 = asort(l)
                        mstr = mstr + 1
                     endif
                     cstr(ksort(l),mstr:nxs+1) = 1
                  end do
                  nstr(mrgn) = mstr
                  wstr(mstr,mrgn) = 1.0_r8 - cld0
!
! End test of region_found = true
!
               endif
!
! End loop over regions irgn for max-overlap
!
            end do
            nrgn = mrgn
!
! Finish construction of cstr for additional top layer
!
            cstr(0,1:nstr(1)) = 0
!
! INDEX COMPUTATIONS FOR STEP 2-3
! This section of the code generates the following information:
! (1. totwgt     step 3     total frac. area of configurations satisfying
! areamin & nconfgmax criteria
! (2. wgtv       step 3     frac. area of configurations
! (3. ccon       step 2     binary flag for clouds in each configuration
! (4. nconfig    steps 2-3  number of configurations
! (5. nuniqu/d   step 2     Number of unique cloud configurations for
! up/downwelling rad. between surface/TOA
! and level k
! (6. istrtu/d   step 2     Indices into iconu/d
! (7. iconu/d    step 2     Cloud configurations which are identical
! for up/downwelling rad. between surface/TOA
! and level k
!
! Number of configurations (all permutations of streams in each region)
!
            nconfigm = product(nstr(1: nrgn))
!
! Construction of totwgt, wgtv, ccon, nconfig
!
            istr(1: nrgn) = 1
            nconfig = 0
            totwgt = 0.0_r8
            new_term = .true.
            do iconfig = 1, nconfigm
               xwgt = 1.0_r8
               do mrgn = 1,  nrgn
                  xwgt = xwgt * wstr(istr(mrgn),mrgn)
               end do
               if (xwgt >= areamin) then
                  nconfig = nconfig + 1
                  if (nconfig <= nconfgmax) then
                     j = nconfig
                     ptrc(nconfig) = nconfig
                  else
                     nconfig = nconfgmax
                     if (new_term) then
                        j = findvalue(1,nconfig,wgtv,ptrc)
                     endif
                     if (wgtv(j) < xwgt) then
                        totwgt = totwgt - wgtv(j)
                        new_term = .true.
                     else
                        new_term = .false.
                     endif
                  endif
                  if (new_term) then
                     wgtv(j) = xwgt
                     totwgt = totwgt + xwgt
                     do mrgn = 1, nrgn
                        ccon(kx1(mrgn):kx2(mrgn),j) = cstr(kx1(mrgn):kx2(mrgn),istr(mrgn))
                     end do
                  endif
               endif

               mrgn =  nrgn
               istr(mrgn) = istr(mrgn) + 1
               do while (istr(mrgn) > nstr(mrgn) .and. mrgn > 1)
                  istr(mrgn) = 1
                  mrgn = mrgn - 1
                  istr(mrgn) = istr(mrgn) + 1
               end do
!
! End do iconfig = 1, nconfigm
!
            end do
!
! If totwgt = 0 implement maximum overlap and make another pass
! if totwgt = 0 on this second pass then terminate.
!
            if (totwgt > 0.) then
               exit
            else
               npasses = npasses + 1
               if (npasses >= 2 ) then
                  write(6,*)'RADCSWMX: Maximum overlap of column ','failed'
                  call endrun
               endif
               nmxrgn(i)=1
               pmxrgn(i,1)=1.0e30
            end if
!
! End npasses = 0, do
!
         end do
!
!
! Finish construction of ccon
!
         ccon(0,:) = 0
         ccon(pverp,:) = 0
!
! Construction of nuniqu/d, istrtu/d, iconu/d using binary tree
!
         nuniqd(0) = 1
         nuniqu(pverp) = 1

         istrtd(0,1) = 1
         istrtu(pverp,1) = 1

         do j = 1, nconfig
            icond(0,j)=j
            iconu(pverp,j)=j
         end do

         istrtd(0,2) = nconfig+1
         istrtu(pverp,2) = nconfig+1

         do k = 1, pverp
            km1 = k-1
            nuniq = 0
            istrtd(k,1) = 1
            do l0 = 1, nuniqd(km1)
               is0 = istrtd(km1,l0)
               is1 = istrtd(km1,l0+1)-1
               n0 = 0
               n1 = 0
               do isn = is0, is1
                  j = icond(km1,isn)
                  if (ccon(k,j) == 0) then
                     n0 = n0 + 1
                     ptr0(n0) = j
                  endif
                  if (ccon(k,j) == 1) then
                     n1 = n1 + 1
                     ptr1(n1) = j
                  endif
               end do
               if (n0 > 0) then
                  nuniq = nuniq + 1
                  istrtd(k,nuniq+1) = istrtd(k,nuniq)+n0
                  icond(k,istrtd(k,nuniq):istrtd(k,nuniq+1)-1) =  ptr0(1:n0)
               endif
               if (n1 > 0) then
                  nuniq = nuniq + 1
                  istrtd(k,nuniq+1) = istrtd(k,nuniq)+n1
                  icond(k,istrtd(k,nuniq):istrtd(k,nuniq+1)-1) =  ptr1(1:n1)
               endif
            end do
            nuniqd(k) = nuniq
         end do

         do k = pver, 0, -1
            kp1 = k+1
            nuniq = 0
            istrtu(k,1) = 1
            do l0 = 1, nuniqu(kp1)
               is0 = istrtu(kp1,l0)
               is1 = istrtu(kp1,l0+1)-1
               n0 = 0
               n1 = 0
               do isn = is0, is1
                  j = iconu(kp1,isn)
                  if (ccon(k,j) == 0) then
                     n0 = n0 + 1
                     ptr0(n0) = j
                  endif
                  if (ccon(k,j) == 1) then
                     n1 = n1 + 1
                     ptr1(n1) = j
                  endif
               end do
               if (n0 > 0) then
                  nuniq = nuniq + 1
                  istrtu(k,nuniq+1) = istrtu(k,nuniq)+n0
                  iconu(k,istrtu(k,nuniq):istrtu(k,nuniq+1)-1) =  ptr0(1:n0)
               endif
               if (n1 > 0) then
                  nuniq = nuniq + 1
                  istrtu(k,nuniq+1) = istrtu(k,nuniq)+n1
                  iconu(k,istrtu(k,nuniq):istrtu(k,nuniq+1)-1) = ptr1(1:n1)
               endif
            end do
            nuniqu(k) = nuniq
         end do
!
!----------------------------------------------------------------------
! End of index calculations
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! Start of flux calculations
!----------------------------------------------------------------------
!
! Initialize spectrally integrated totals:
!
         do k=0,pver
            totfld(k) = 0.0_r8
            fswup (k) = 0.0_r8
            fswdn (k) = 0.0_r8
         end do

         sfltot        = 0.0_r8
         fswup (pverp) = 0.0_r8
         fswdn (pverp) = 0.0_r8
!
! Start spectral interval
!
         do ns = 1,nspint
            wgtint = nirwgt(ns)
!----------------------------------------------------------------------
! STEP 2
!
!
! Apply adding method to solve for radiative properties
!
! First initialize the bulk properties at TOA
!
            rdndif(0,1:nconfig) = 0.0_r8
            exptdn(0,1:nconfig) = 1.0_r8
            tdntot(0,1:nconfig) = 1.0_r8
!
! Solve for properties involving downward propagation of radiation.
! The bulk properties are:
!
! (1. exptdn   Sol. beam dwn. trans from layers above
! (2. rdndif   Ref to dif rad for layers above
! (3. tdntot   Total trans for layers above
!
            do k = 1, pverp
               km1 = k - 1
               do l0 = 1, nuniqd(km1)
                  is0 = istrtd(km1,l0)
                  is1 = istrtd(km1,l0+1)-1

                  j = icond(km1,is0)

                  xexpt   = exptdn(km1,j)
                  xrdnd   = rdndif(km1,j)
                  tdnmexp = tdntot(km1,j) - xexpt

                  if (ccon(km1,j) == 1) then
!
! If cloud in layer, use cloudy layer radiative properties
!
                     ytdnd = tdif(ns,i,km1)
                     yrdnd = rdif(ns,i,km1)

                     rdenom  = 1._r8/(1._r8-yrdnd*xrdnd)
                     rdirexp = rdir(ns,i,km1)*xexpt

                     zexpt = xexpt * explay(ns,i,km1)
                     zrdnd = yrdnd + xrdnd*(ytdnd**2)*rdenom
                     ztdnt = xexpt*tdir(ns,i,km1) + ytdnd*(tdnmexp + xrdnd*rdirexp)*rdenom
                  else
!
! If clear layer, use clear-sky layer radiative properties
!
                     ytdnd = tdifc(ns,i,km1)
                     yrdnd = rdifc(ns,i,km1)

                     rdenom  = 1._r8/(1._r8-yrdnd*xrdnd)
                     rdirexp = rdirc(ns,i,km1)*xexpt

                     zexpt = xexpt * explayc(ns,i,km1)
                     zrdnd = yrdnd + xrdnd*(ytdnd**2)*rdenom
                     ztdnt = xexpt*tdirc(ns,i,km1) + ytdnd* &
                                            (tdnmexp + xrdnd*rdirexp)*rdenom
                  endif

!
! If 2 or more configurations share identical properties at a given level k,
! the properties (at level k) are computed once and copied to
! all the configurations for efficiency.
!
                  do isn = is0, is1
                     j = icond(km1,isn)
                     exptdn(k,j) = zexpt
                     rdndif(k,j) = zrdnd
                     tdntot(k,j) = ztdnt
                  end do
!
! end do l0 = 1, nuniqd(k)
!
               end do
!
! end do k = 1, pverp
!
            end do
!
! Solve for properties involving upward propagation of radiation.
! The bulk properties are:
!
! (1. rupdif   Ref to dif rad for layers below
! (2. rupdir   Ref to dir rad for layers below
!
! Specify surface boundary conditions (surface albedos)
!
            rupdir(pverp,1:nconfig) = albdir(i,ns)
            rupdif(pverp,1:nconfig) = albdif(i,ns)

            do k = pver, 0, -1
               do l0 = 1, nuniqu(k)
                  is0 = istrtu(k,l0)
                  is1 = istrtu(k,l0+1)-1

                  j = iconu(k,is0)

                  xrupd = rupdif(k+1,j)
                  xrups = rupdir(k+1,j)

                  if (ccon(k,j) == 1) then
!
! If cloud in layer, use cloudy layer radiative properties
!
                     yexpt = explay(ns,i,k)
                     yrupd = rdif(ns,i,k)
                     ytupd = tdif(ns,i,k)

                     rdenom  = 1._r8/( 1._r8 - yrupd*xrupd)
                     tdnmexp = (tdir(ns,i,k)-yexpt)
                     rdirexp = xrups*yexpt

                     zrupd = yrupd + xrupd*(ytupd**2)*rdenom
                     zrups = rdir(ns,i,k) + ytupd*(rdirexp + xrupd*tdnmexp)*rdenom
                  else
!
! If clear layer, use clear-sky layer radiative properties
!
                     yexpt = explayc(ns,i,k)
                     yrupd = rdifc(ns,i,k)
                     ytupd = tdifc(ns,i,k)

                     rdenom  = 1._r8/( 1._r8 - yrupd*xrupd)
                     tdnmexp = (tdirc(ns,i,k)-yexpt)
                     rdirexp = xrups*yexpt

                     zrupd = yrupd + xrupd*(ytupd**2)*rdenom
                     zrups = rdirc(ns,i,k) + ytupd*(rdirexp + xrupd*tdnmexp)*rdenom
                  endif

!
! If 2 or more configurations share identical properties at a given level k,
! the properties (at level k) are computed once and copied to
! all the configurations for efficiency.
!
                  do isn = is0, is1
                     j = iconu(k,isn)
                     rupdif(k,j) = zrupd
                     rupdir(k,j) = zrups
                  end do
!
! end do l0 = 1, nuniqu(k)
!
               end do
!
! end do k = pver,0,-1
!
            end do
!
!----------------------------------------------------------------------
!
! STEP 3
!
! Compute up and down fluxes for each interface k.  This requires
! adding up the contributions from all possible permutations
! of streams in all max-overlap regions, weighted by the
! product of the fractional areas of the streams in each region
! (the random overlap assumption).  The adding principle has been
! used in step 2 to combine the bulk radiative properties
! above and below the interface.
!
            do k = 0,pverp
!
! Initialize the fluxes
!
               fluxup(k)=0.0_r8
               fluxdn(k)=0.0_r8

               do iconfig = 1, nconfig
                  xwgt = wgtv(iconfig)
                  xexpt = exptdn(k,iconfig)
                  xtdnt = tdntot(k,iconfig)
                  xrdnd = rdndif(k,iconfig)
                  xrupd = rupdif(k,iconfig)
                  xrups = rupdir(k,iconfig)
!
! Flux computation
!
                  rdenom = 1._r8/(1._r8 - xrdnd * xrupd)

                  fluxup(k) = fluxup(k) + xwgt *  &
                              ((xexpt * xrups + (xtdnt - xexpt) * xrupd) * rdenom)
                  fluxdn(k) = fluxdn(k) + xwgt *  &
                              (xexpt + (xtdnt - xexpt + xexpt * xrups * xrdnd) * rdenom)
!
! End do iconfig = 1, nconfig
!
               end do
!
! Normalize by total area covered by cloud configurations included
! in solution
!
               fluxup(k)=fluxup(k) / totwgt
               fluxdn(k)=fluxdn(k) / totwgt
!
! End do k = 0,pverp
!
            end do
!
! Initialize the direct-beam flux at surface
!
            wexptdn = 0.0_r8

            do iconfig = 1, nconfig
               wexptdn =  wexptdn + wgtv(iconfig) * exptdn(pverp,iconfig)
            end do

            wexptdn = wexptdn / totwgt
!
! Monochromatic computation completed; accumulate in totals
!
            solflx   = solin(i)*frcsol(ns)*psf(ns)
            fsnt(i)  = fsnt(i) + solflx*(fluxdn(1) - fluxup(1))
            fsntoa(i)= fsntoa(i) + solflx*(fluxdn(0) - fluxup(0))
            fsns(i)  = fsns(i) + solflx*(fluxdn(pverp)-fluxup(pverp))
            swuptoa(i) = swuptoa(i)+solflx*fluxup(0) ! added by DONG Li
            swdntoa(i) = swdntoa(i)+solflx*fluxdn(0) !
            sfltot   = sfltot + solflx
            fswup(0) = fswup(0) + solflx*fluxup(0)
            fswdn(0) = fswdn(0) + solflx*fluxdn(0)
!
! Down spectral fluxes need to be in mks; thus the .001 conversion factors
!
            if (wavmid(ns) < 0.7_r8) then
               sols(i)  = sols(i) + wexptdn*solflx*0.001_r8
               solsd(i) = solsd(i)+(fluxdn(pverp)-wexptdn)*solflx*0.001_r8
            else
               soll(i)  = soll(i) + wexptdn*solflx*0.001_r8
               solld(i) = solld(i)+(fluxdn(pverp)-wexptdn)*solflx*0.001_r8
               fsnrtoaq(i) = fsnrtoaq(i) + solflx*(fluxdn(0) - fluxup(0))
            end if
            fsnirtoa(i) = fsnirtoa(i) + wgtint*solflx*(fluxdn(0) - fluxup(0))

            do k=0,pver
!
! Compute flux divergence in each layer using the interface up and down
! fluxes:
!
               kp1 = k+1
               flxdiv = (fluxdn(k  ) - fluxdn(kp1)) + (fluxup(kp1) - fluxup(k  ))
               totfld(k)  = totfld(k)  + solflx*flxdiv
               fswdn(kp1) = fswdn(kp1) + solflx*fluxdn(kp1)
               fswup(kp1) = fswup(kp1) + solflx*fluxup(kp1)
            end do
!
! Perform clear-sky calculation
!
            exptdnc(0) =   1.0_r8
            rdndifc(0) =   0.0_r8
            tdntotc(0) =   1.0_r8
            rupdirc(pverp) = albdir(i,ns)
            rupdifc(pverp) = albdif(i,ns)

            do k = 1, pverp
               km1 = k - 1
               xexpt = exptdnc(km1)
               xrdnd = rdndifc(km1)
               yrdnd = rdifc(ns,i,km1)
               ytdnd = tdifc(ns,i,km1)

               exptdnc(k) = xexpt*explayc(ns,i,km1)

               rdenom  = 1._r8/(1._r8 - yrdnd*xrdnd)
               rdirexp = rdirc(ns,i,km1)*xexpt
               tdnmexp = tdntotc(km1) - xexpt

               tdntotc(k) = xexpt*tdirc(ns,i,km1) + ytdnd*(tdnmexp + xrdnd*rdirexp)* &
                                rdenom
               rdndifc(k) = yrdnd + xrdnd*(ytdnd**2)*rdenom
            end do

            do k=pver,0,-1
               xrupd = rupdifc(k+1)
               yexpt = explayc(ns,i,k)
               yrupd = rdifc(ns,i,k)
               ytupd = tdifc(ns,i,k)

               rdenom = 1._r8/( 1._r8 - yrupd*xrupd)

               rupdirc(k) = rdirc(ns,i,k) + ytupd*(rupdirc(k+1)*yexpt + &
                            xrupd*(tdirc(ns,i,k)-yexpt))*rdenom
               rupdifc(k) = yrupd + xrupd*ytupd**2*rdenom
            end do

            do k=0,1
               rdenom    = 1._r8/(1._r8 - rdndifc(k)*rupdifc(k))
               fluxup(k) = (exptdnc(k)*rupdirc(k) + (tdntotc(k)-exptdnc(k))*rupdifc(k))* &
                           rdenom
               fluxdn(k) = exptdnc(k) + &
                           (tdntotc(k) - exptdnc(k) + exptdnc(k)*rupdirc(k)*rdndifc(k))* &
                           rdenom
            end do
            k = pverp
            rdenom      = 1._r8/(1._r8 - rdndifc(k)*rupdifc(k))
            fluxup(k)   = (exptdnc(k)*rupdirc(k) + (tdntotc(k)-exptdnc(k))*rupdifc(k))* &
                           rdenom
            fluxdn(k)   = exptdnc(k) + (tdntotc(k) - exptdnc(k) + &
                          exptdnc(k)*rupdirc(k)*rdndifc(k))*rdenom

            fsntc(i)    = fsntc(i)+solflx*(fluxdn(1)-fluxup(1))
            fsntoac(i)  = fsntoac(i)+solflx*(fluxdn(0)-fluxup(0))
            fsutoac(i)  = fsutoac(i)+solflx*fluxup(0)
            fsnsc(i)    = fsnsc(i)+solflx*(fluxdn(pverp)-fluxup(pverp))
            fsdsc(i)    = fsdsc(i)+solflx*(fluxdn(pverp))
            fsusc(i)    = fsusc(i)+solflx*(fluxup(pverp))
            fsnrtoac(i) = fsnrtoac(i)+wgtint*solflx*(fluxdn(0)-fluxup(0))
!
! End of clear sky calculation
!

!
! End of spectral interval loop
!
         end do
!
! Compute solar heating rate (J/kg/s)
!
         do k=1,pver
            qrs(i,k) = -1.E-4*gravit*totfld(k)/(pint(i,k) - pint(i,k+1))
         end do
!
! Set the downwelling flux at the surface
!
        fsds(i) = fswdn(pverp)
        fsus(i) = fswup(pverp)
!
! End do n=1,ndayc
!
    end do

    call outfld('FSUTOA  ', swuptoa*1.0e-3,       pcols, lchnk)
    call outfld('FSDTOA  ', swdntoa*1.0e-3,       pcols, lchnk)

   return
end subroutine radcswmx
