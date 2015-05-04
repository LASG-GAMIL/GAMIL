#include <misc.h>
#include <params.h>

module dynamics_vars
!BOP
!
! !MODULE: dynamics_vars --- Lin-Rood specific variables and methods
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid

! !PUBLIC MEMBER FUNCTIONS:
   public dynamics_init, dynamics_clean

! !PUBLIC DATA MEMBERS:

!--------------------------------------------------------
! Local variables specific to the Lin-Rood dynamical core
!--------------------------------------------------------

! Temporary variables to be removed at a later time
   real(r8), allocatable, save :: xyt(:,:,:)
   real(r8), allocatable, save :: yzt(:,:,:)
   real(r8), allocatable, save :: q3t(:,:,:,:)

! Variables set by dynamics_init
   real(r8) dtime        ! large time step
   integer  iord         ! order of LR scheme in X
   integer  jord         ! order of LR scheme in Y

! General variables set by setrig
   real (r8) pi
   real(r8) dl           ! Radians per (finite-volume) longitude
   real(r8) dp           ! Radians per (finite-volume) latitude

! Geometric arrays set by setrig
   real(r8) cosp(plat)   ! Cos of latitude angle - volume mean
   real(r8) sinp(plat)   ! Sin of latitude angle - volume mean
   real(r8) cose(plat)   ! Cos of finite-volume edges
   real(r8) sine(plat)   ! Sin of finite-volume edges
   real(r8) gw(plat)     ! Delta sine

! Variables set by set_eta
   real(r8) ptop         ! pressure at top of atmosphere
   real(r8) pint         ! pressure at sig-p interface
   real(r8) ak(plev+1)   ! A's of the ETA-coordinate
   real(r8) bk(plev+1)   ! B's of the ETA-coordinate
   integer  ks           ! Total number of pure-P layers

! Variables set by rayf_init
   real(r8) rfac(plev)   ! Rayleigh friction factor

! Needed in Held-Suarez (hswf) set by hswf_init
   real (r8) sinp2(plat)
   real (r8) cosp2(plat)
   real (r8) cosp4(plat)
   real (r8) rf(plev)

! Scalars set in dynpkg_init

   integer icd, jcd

   integer ng_c       ! ghost zone needed by the c-gird dynamics
   integer ng_d       ! ghost zone needed by the d-gird dynamics
   integer ng_s       ! for certain arrays, max(ng_c+1,ng_d)

   real (r8) acap     ! scaled polar cap area
   real (r8) rcap     ! inverse of scaled polar cap area

! Arrays initialized by dynpkg_init

   real(r8) coslon(plon) ! Cos of longitudes - volume center
   real(r8) sinlon(plon) ! Sin of longitudes - volume center

   real(r8) cosl5(plon)  ! Cos of longitudes - volume center
   real(r8) sinl5(plon)  ! Sin of longitudes - volume center

   real(r8) acosp(plat)

! Scalars initialized by d_split
   integer   ns       ! total number of splits for Lagrangian dynamics

! Constants  -- these might be user-configurable at a later date
   logical rayf         ! Rayleigh friction flag (off by default)
   parameter (rayf = .false.)         ! off
   logical dcaf         ! Dry convection flag (use only in ideal_phys case)
   parameter (dcaf = .false.)         ! off

!
! !DESCRIPTION:
!
!      This module provides variables which are specific to the Lin-Rood
!      dynamical core.  Most of them were previously SAVE variables in 
!      different routines and were set with an "if (first)" statement.
!
!      \begin{tabular}{|l|l|} \hline \hline
!        lr\_init    &  Initialize the Lin-Rood variables  \\ \hline
!        lr\_clean   &  Deallocate all internal data structures \\ \hline 
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   01.06.06   Sawyer     Consolidated from various code snippets
!   01.07.12   Sawyer     Removed CCM common blocks comtim.h and commap.h
!
! !BUGS:
!   o rayf_init and hswf_init should only be performed if they are used....
!   o Where does the value of ns0 come from??
!   o Should gw be here??  It's available in commap.h (physics variable)
!     and set in stepon.  Needed in gmean.
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: dynamics_init --- initialize the lin-rood dynamical core
!
! !INTERFACE: 
   subroutine dynamics_init( dtime_in, iord_in, jord_in, nsplit_in )

! !USES:
      use constituents, only : pcnst, pnats  ! TEMPORARY!!!!
      implicit none

! !INPUT PARAMETERS:
      real (r8), intent(in) :: dtime_in   !  Large time step
      integer, intent(in)   :: iord_in    !  Order of LR scheme in X
      integer, intent(in)   :: jord_in    !  Order of LR scheme in Y
      integer, intent(in)   :: nsplit_in  !  Small time steps in dtime

! !DESCRIPTION:
!
!   Initialize Lin-Rood specific variables
!
! !REVISION HISTORY:
!
!   01.06.06   Sawyer     Create
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
      dtime = dtime_in
      iord  = iord_in
      jord  = jord_in

      call setrig
      call set_eta
      if ( rayf ) call rayf_init
      call hswf_init
      call dynpkg_init
      call d_split(nsplit_in)
#if defined( SPMD )
! Temporary data structures
      allocate(   xyt(beglonxy:endlonxy,beglatxy:endlatxy,plev) )
      allocate(   yzt(plon,beglat:endlat,beglev:endlev) )
      allocate(   q3t(plon,beglat:endlat,beglev:endlev,pcnst+pnats) )
#endif
      return
!EOC
   end subroutine dynamics_init
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: dynamics_clean -- clean up Lin-Rood-specific variables
!
! !INTERFACE: 
   subroutine dynamics_clean

      implicit none

! !DESCRIPTION:
!
! Clean up (deallocate) Lin-Rood-specific variables
!
! !REVISION HISTORY:
!
!   01.06.06   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! Temporary data structures
   deallocate( q3t )
   deallocate( yzt )
   deallocate( xyt )
   return
!EOC
   end subroutine dynamics_clean
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  setrig --- Specify the grid attributes
!
! !INTERFACE:
      subroutine setrig

! !USES:
      implicit none

!
! !DESCRIPTION:
!
!   Specify the grid attributes, such as the spacing between
!   grid points in latitude and longitude, the sines and cosines of
!   latitude at cell midpoints and edges.
!
! !REVISION HISTORY: 
!   ??.??.??    Lin?       Creation
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer j
      real*8 pi, ph5      ! This is to ensure 64-bit for any choice of r8

      pi  = 4.0 * atan(1.)
      dl  = (pi+pi)/plon
      dp  = pi/(plat-1)

      do j=2,plat
         ph5  = -0.5d0*pi + ((j-1)-0.5d0)*(pi/(plat-1))
         sine(j) = sin(ph5)
      enddo

      do j=2,plat-1
         gw(j) = sine(j+1) - sine(j)
      end do
      gw( 1)   =  1. + sine(2)
      gw(plat) =  1. - sine(plat)

      cosp( 1) =  0.
      cosp(plat) =  0.

      do j=2,plat-1
         cosp(j) = (sine(j+1)-sine(j)) / dp
      enddo

! Define cosine at edges..

      do j=2,plat
         cose(j) = 0.5 * (cosp(j-1) + cosp(j))
      enddo
         cose(1) = cose(2)

         sinp( 1) = -1.
         sinp(plat) =  1.

      do j=2,plat-1
         sinp(j) = 0.5 * (sine(j) + sine(j+1))
      enddo

      return
!EOC
      end subroutine setrig
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  set_eta --- Define vertical coordinate
!
! !INTERFACE:
      subroutine set_eta

! !USES:
      implicit none

!
! !DESCRIPTION:
!
!   Specify the vertical coordinate system.  Currently this is a 
!   dual pressure - sigma coordinate system (???), which transitions at
!   level ks, but it could be just about anything reasonable.
!
!   Choices for vertical resolutions are as follows:
!   \begin{tabular}{l}
!     NCAR: 18, 26, and 30 \\
!     NASA DAO: smoothed version of CAM's 30-level, 32, 55, 64, and 96 \\
!     New 32-layer setup with top at 0.4 mb for high horizontal
!     resolution runs. pint = 176.93 mb \\
!     Revised 55-level eta with pint at 176.93 mb  SJL: 2000-03-20
!   \end{tabular}
! 
! !REVISION HISTORY: 
!   98.01.15    Lin        Creation
!   ongoing     Lin        Fine tuning
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! NCAR specific
      real(r8) a18(19),b18(19)              ! CCM3
      real(r8) a26(27),b26(27)              ! CAM
      real(r8) a30(31),b30(31)              ! CCM3

! NASA only
      real(r8) a30m(31),b30m(31)            ! smoothed CAM 30-L
      real(r8) a32(33),b32(33)
      real(r8) a32old(33),b32old(33)
      real(r8) a55(56),b55(56)
      real(r8) a55old(56),b55old(56)
      real(r8) a64(65),b64(65)
      real(r8) a96(97),b96(97)

      integer k

! *** NCAR settings ***

      data a18 /291.70,  792.92,  2155.39,  4918.34,  8314.25, &
               7993.08, 7577.38,  7057.52,  6429.63,  5698.38, &
               4879.13, 3998.95,  3096.31,  2219.02,  1420.39, &
               754.13,  268.38,   0.0000,   0.0000 /

      data b18 /0.0000,    0.0000,    0.0000,   0.0000,   0.0000, &
                0.0380541, 0.0873088, 0.1489307, 0.2232996,       &
                0.3099406, 0.4070096, 0.5112977, 0.6182465,       &
                0.7221927, 0.8168173, 0.8957590, 0.9533137,       &
                0.9851122, 1.0  /
     
      data a26 /219.4067,  489.5209,   988.2418,   1805.201,      &
                2983.724,  4462.334,   6160.587,   7851.243,      &
                7731.271,  7590.131,   7424.086,   7228.744,      &
                6998.933,  6728.574,   6410.509,   6036.322,      &
                5596.111,  5078.225,   4468.96,    3752.191,      &
                2908.949,  2084.739,   1334.443,   708.499,       &
                252.136,   0.,         0. /

      data b26 /0.,         0.,         0.,         0.,           &
                0.,         0.,         0.,         0.,           &
                0.01505309, 0.03276228, 0.05359622, 0.07810627,   &
                0.1069411,  0.14086370, 0.180772,   0.227722,     &
                0.2829562,  0.3479364,  0.4243822,  0.5143168,    &
                0.6201202,  0.7235355,  0.8176768,  0.8962153,    &
                0.9534761,  0.9851122,  1.        /

      data a30 /225.523952394724, 503.169186413288, 1015.79474285245,  &
               1855.53170740604, 3066.91229343414,  4586.74766123295,  &
               6332.34828710556, 8070.14182209969,  9494.10423636436,  &
              11169.321089983,  13140.1270627975,  15458.6806893349,   &
              18186.3352656364, 17459.799349308,   16605.0657629967,   &
              15599.5160341263, 14416.541159153,   13024.8308181763,   &
              11387.5567913055,  9461.38575673103,  7534.44507718086,  &
               5765.89405536652, 4273.46378564835,  3164.26791250706,  &
               2522.12174236774, 1919.67375576496,  1361.80268600583,  &
                853.108894079924, 397.881818935275,    0.,             &
                  0.  /

      data b30 /0.,                 0.,                                 &
                0.,                 0.,                0.,              &
                0.,                 0.,                0.,              &
                0.,                 0.,                0.,              &
                0.,                 0.,                0.03935482725501,&
                0.085653759539127,  0.140122056007385, 0.20420117676258,&
                0.279586911201477,  0.368274360895157, 0.47261056303978,&
                0.576988518238068,  0.672786951065063, 0.75362843275070,&
                0.813710987567902,  0.848494648933411, 0.88112789392471,&
                0.911346435546875,  0.938901245594025, 0.96355980634689,&
                0.985112190246582,  1.   /

! *** NASA DAO settings ***

! Smoothed CAM's 30-Level setup
      data a30m / 300.00000,     725.00000,    1500.00000,     &
             2600.00000,    3800.00000,    5050.00000,         &
             6350.00000,    7750.00000,    9300.00000,         &
            11100.00000,   13140.00000,   15458.00000,         &
            18186.33580,   20676.23761,   22275.23783,         &
            23025.65071,   22947.33569,   22038.21991,         &
            20274.24578,   17684.31619,   14540.98138,         &
            11389.69990,    8795.97971,    6962.67963,         &
             5554.86684,    4376.83633,    3305.84967,         &
             2322.63910,    1437.78398,     660.76994,         &
                0.00000 /

      data b30m / 0.00000,       0.00000,       0.00000,       &
                  0.00000,       0.00000,       0.00000,       &
                  0.00000,       0.00000,       0.00000,       &
                  0.00000,       0.00000,       0.00000,       &
                  0.00000,       0.00719,       0.02895,       &
                  0.06586,       0.11889,       0.18945,       &
                  0.27941,       0.38816,       0.50692,       &
                  0.61910,       0.70840,       0.77037,       &
                  0.81745,       0.85656,       0.89191,       &
                  0.92421,       0.95316,       0.97850,       &
                  1.00000 /

      data a32/40.00000,     100.00000,     200.00000,         &
            370.00000,     630.00000,    1000.00000,           &
           1510.00000,    2160.00000,    2900.00000,           &
           3680.00000,    4535.00000,    5505.00000,           &
           6607.26750,    7851.22980,    9236.56610,           &
          10866.34270,   12783.70000,   15039.30000,           &
          17693.00000,   20119.20876,   21686.49129,           &
          22436.28749,   22388.46844,   21541.75227,           &
          19873.78342,   17340.31831,   13874.44006,           &
          10167.16551,    6609.84274,    3546.59643,           &
           1270.49390,       0.00000,       0.00000   /

      data b32/0.00000,       0.00000,       0.00000,          &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00696,       0.02801,            &
             0.06372,       0.11503,       0.18330,            &
             0.27033,       0.37844,       0.51046,            &
             0.64271,       0.76492,       0.86783,            &
             0.94329,       0.98511,       1.00000   /

      data a32old /300.0000,  454.1491,  652.5746,  891.9637, 1159.7102, &
             1492.8248, 1902.5026, 2400.4835, 2998.6740, 3708.6584,      &
             4541.1041, 5505.0739, 6607.2675, 7851.2298, 9236.5661,      &
            10866.3427, 12420.400, 13576.500, 14365.400, 14807.800,      &
             14915.500, 14691.400, 14129.400, 13214.800, 11923.200,      &
             10220.700,  8062.000,  5849.500,  3777.000,  2017.200,      &
               720.600,     0.000,     0.000 /

      data b32old /0.00, 0.0000000, 0.0000000, 0.0000000, 0.0000000,     &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,     &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,     &
              0.0000000, 0.003633 , 0.014628 , 0.033276 , 0.060071 ,     &
              0.095722 , 0.141171 , 0.197623 , 0.266571 , 0.349839 ,     &
              0.449632 , 0.568589 , 0.685887 , 0.793252 , 0.883128 ,     &
              0.948792 , 0.9851119, 1.0000000 /

      data a55/ 1.00000,       2.00000,       3.27000,              &
              4.75850,       6.60000,       8.93450,                &
             11.97030,      15.94950,      21.13490,                &
             27.85260,      36.50410,      47.58060,                &
             61.67790,      79.51340,     101.94420,                &
            130.05080,     165.07920,     208.49720,                &
            262.02120,     327.64330,     407.65670,                &
            504.68050,     621.68000,     761.98390,                &
            929.29430,    1127.68880,    1364.33920,                &
           1645.70720,    1979.15540,    2373.03610,                &
           2836.78160,    3380.99550,    4017.54170,                &
           4764.39320,    5638.79380,    6660.33770,                &
           7851.22980,    9236.56610,   10866.34270,                &
          12783.70000,   15039.30000,   17693.00000,                &
          20119.20876,   21686.49129,   22436.28749,                &
          22388.46844,   21541.75227,   19873.78342,                &
          17340.31831,   13874.44006,   10167.16551,                &
           6609.84274,    3546.59643,    1270.49390,                &
              0.00000,       0.00000   /

      data b55 / 0.00000,       0.00000,       0.00000,         &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00696,       0.02801,       0.06372,           &
               0.11503,       0.18330,       0.27033,           &
               0.37844,       0.51046,       0.64271,           &
               0.76492,       0.86783,       0.94329,           &
               0.98511,       1.00000  /

      data a55old /1.0000,    2.0000,    3.2700,    4.7585,     6.6000, &
              8.9345,   11.9703,   15.9495,   21.1349,    27.8526,      &
             36.5041,   47.5806,   61.6779,   79.5134,   101.9442,      &
            130.0508,  165.0792,  208.4972,  262.0212,   327.6433,      &
            407.6567,  504.6805,  621.6800,  761.9839,   929.2943,      &
           1127.6888, 1364.3392, 1645.7072, 1979.1554,  2373.0361,      &
           2836.7816, 3380.9955, 4017.5417, 4764.3932,  5638.7938,      &
           6660.3377, 7851.2298, 9236.5661,10866.3427, 12420.400 ,      &
          13576.500 , 14365.400, 14807.800, 14915.500 , 14691.400,      &
          14129.400 , 13214.800, 11923.200, 10220.700 ,  8062.000,      &
           5849.500 ,  3777.000,  2017.200,   720.600,      0.000,      &
              0.000 /

      data b55old /   0.0000000, 0.0000000, 0.0000000,                  &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.003633 , 0.014628 , 0.033276 , 0.060071 ,    &
              0.095722 , 0.141171 , 0.197623 , 0.266571 , 0.349839 ,    &
              0.449632 , 0.568589 , 0.685887 , 0.793252 , 0.883128 ,    &
              0.948792 , 0.9851119, 1.0000000 /

      data a64/1.00000,       1.80   ,       2.80086,     &
             3.93309,       5.20139,       6.77626,       &
             8.69654,      10.99483,      13.81736,       &
            17.26058,      21.43286,      26.45448,       &
            32.45730,      39.58402,      47.98678,       &
            57.82525,      69.26401,      82.46925,       &
            97.60468,     114.82686,     135.08787,       &
           158.92390,     186.96575,     219.95555,       &
           258.76633,     304.42522,     358.14053,       &
           421.33383,     495.67748,     583.13893,       &
           686.03282,     807.08215,     949.49044,       &
          1117.02644,    1314.12387,    1545.99882,       &
          1818.78771,    2139.70974,    2517.25793,       &
          2961.42386,    3483.96212,    4098.70138,       &
          4821.91034,    5672.72831,    6673.67169,       &
          7851.22983,    9236.56613,   10866.34270,       &
         12783.69059,   15039.35130,   17693.01955,       &
         20814.92310,   23609.16397,   25271.17281,       &
         25844.93368,   25345.63084,   23760.05052,       &
         21046.23129,   17132.35351,   12832.00555,       &
          8646.27815,    5012.23907,    2299.34286,       &
           781.15294,       0.00000  /

      data b64/0.00000,       0.00000,       0.00000,     &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00879,       0.03537,       &
             0.08047,       0.14526,       0.23147,       &
             0.34138,       0.47789,       0.61606,       &
             0.74456,       0.85318,       0.93300,       &
             0.97730,       1.00000  /

      data a96/  1.00000,       2.32782,       3.34990,   &
               4.49484,       5.62336,       6.93048,     &
               8.41428,      10.06365,      11.97630,     &
              14.18138,      16.70870,      19.58824,     &
              22.84950,      26.52080,      30.62845,     &
              35.19588,      40.24273,      45.78375,     &
              51.82793,      58.43583,      65.62319,     &
              73.40038,      81.77154,      90.73373,     &
             100.27628,     110.82243,     122.47773,     &
             135.35883,     149.59464,     165.32764,     &
             182.71530,     201.93164,     223.16899,     &
             246.63988,     272.57922,     301.24661,     &
             332.92902,     367.94348,     406.64044,     &
             449.40720,     496.67181,     548.90723,     &
             606.63629,     670.43683,     740.94727,     &
             818.87329,     904.99493,    1000.17395,     &
            1105.36304,    1221.61499,    1350.09326,     &
            1492.08362,    1649.00745,    1822.43469,     &
            2014.10168,    2225.92627,    2460.02905,     &
            2718.75195,    3004.68530,    3320.69092,     &
            3669.93066,    4055.90015,    4482.46240,     &
            4953.88672,    5474.89111,    6050.68994,     &
            6687.04492,    7390.32715,    8167.57373,     &
            9026.56445,    9975.89648,   11025.06934,     &
           12184.58398,   13466.04785,   14882.28320,     &
           16447.46289,   18177.25781,   20088.97461,     &
           21886.89453,   23274.16602,   24264.66602,     &
           24868.31641,   25091.15430,   24935.41016,     &
           24399.52148,   23478.13281,   22162.01758,     &
           20438.00586,   18288.83984,   15693.01172,     &
           12624.54199,    9584.35352,    6736.55713,     &
            4231.34326,    2199.57910,     747.11890,     &
              0.00000 /

      data b96/0.00000,       0.00000,       0.00000,     &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00315,       0.01263,       0.02853,      &
              0.05101,       0.08030,       0.11669,      &
              0.16055,       0.21231,       0.27249,      &
              0.34169,       0.42062,       0.51005,      &
              0.61088,       0.70748,       0.79593,      &
              0.87253,       0.93400,       0.97764,      &
              1.00000 /

      select case (plev)

! *** Original CCM3 18-Level setup ***
        case (18)
          ks = 4
          do k=1,plev+1
            ak(k) = a18(k)
            bk(k) = b18(k)
          enddo

        case (26)
! CAM 26-Level setup ***
          ks = 7
          do k=1,plev+1
            ak(k) = a26(k)
            bk(k) = b26(k)
          enddo

        case (30)
! CAM 30-Level setup ***
          ks = 12
          do k=1,plev+1
            ak(k) = a30(k)
            bk(k) = b30(k)
          enddo

! *** Revised 32-L setup with ptop at 0.4 mb ***
        case (32)
          ks = 18
          do k=1,plev+1
            ak(k) = a32(k)
            bk(k) = b32(k)
          enddo

! *** Revised 55-L setup with ptop at 0.01 mb ***
        case (55)
          ks = 41
          do k=1,plev+1
            ak(k) = a55(k)
            bk(k) = b55(k)
          enddo

! *** Others ***
        case (64)
          ks = 51
          do k=1,plev+1
            ak(k) = a64(k)
            bk(k) = b64(k)
          enddo

        case (96)
          ks = 77
          do k=1,plev+1
            ak(k) = a96(k)
            bk(k) = b96(k)
          enddo

      end select

          ptop = ak(1)
          pint = ak(ks+1) 

      return
!EOC
      end subroutine set_eta
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  rayf_init --- Initialization for Rayleigh friction
!
! !INTERFACE:
subroutine rayf_init

! !USES:
   implicit none

!------------------------------Commons----------------------------------


! !DESCRIPTION:
! 
!   {\bf Purpose:} Initialization of the Rayleigh friction
! 
! !REVISION HISTORY: 
!   00.01.10    Grant        Creation using code from SJ Lin
!   01.03.26    Sawyer       Added ProTeX documentation
!   01.06.06    Sawyer       Modified for dynamics_vars
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
   integer k, tdt

   real(r8) c1
   real(r8) pc
   real(r8) press(plev)

   tdt = int(dtime)   ! dtime is a variable internal to this module

   write(6,*) 'Time step (seconds) for Rayleigh friction =',tdt
   write(6,*) 'Level, pressure, rfac:'


! e-folding time

   if (ak(1) .le. 50.) then
      c1 = 1. / (5.*24*3600.)
   else
      c1 = 1. / (14.*24.*3600.)
   endif

   pc = max(10., ak(1))

   do k = 1, ks
      press(k) = 0.5*(ak(k) + ak(k+1))
      rfac(k) = tdt*c1*(1.+tanh(1.5*log10(pc/press(k))))
      write(6,*) k, press(k), rfac(k)
   enddo

   return
!EOC
end subroutine rayf_init
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  hswf_init --- Initialization for Held-Suarez
!
! !INTERFACE:
subroutine hswf_init

! !USES:
   implicit none

! !DESCRIPTION:
! 
!   {\bf Purpose:} Initialization of the Held-Suarez Forcing
! 
! !REVISION HISTORY: 
!   01.06.06    Sawyer       Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer j, k
      real (r8) c1, pc, tmp
      real (r8)   pdt                       ! Time-step in seconds

      pdt = int(dtime)   ! dtime is a variable internal to this module
      do j=2,plat-1
        sinp2(j) = sinp(j)**2
        cosp2(j) = cosp(j)**2
      enddo
 
      sinp2(1) = ( 0.5*(-1.+sine(2)) )**2
      sinp2(plat) = sinp2(1)
      cosp2(1) = ( 0.5*cose(2) ) **2
      cosp2(plat) = cosp2(1)

      do j=1,plat
        cosp4(j) = cosp2(j)**2
      enddo

      if ( rayf ) then
        c1 = 1. / (12.*3600)
        pc = 1.
        do k=1,ks           ! ks is a dynamics_vars variable set by set_eta
          tmp = 0.5*(ak(k) + ak(k+1))
          rf(k) = c1*(1.+tanh(1.5*log10(pc/tmp)))
          rf(k) = 1./(1.+pdt*rf(k))
        enddo
      endif
      return
!EOC
end subroutine hswf_init
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  dynpkg_init --- Initialization for dynamics package
!
! !INTERFACE:
subroutine dynpkg_init

! !USES:
   implicit none

! !DESCRIPTION:
! 
!   {\bf Purpose:} Initialization of the Rayleigh friction
! 
! !REVISION HISTORY: 
!   00.01.10    Grant        Creation using code from SJ Lin
!   01.03.26    Sawyer       Added ProTeX documentation
!   01.06.06    Sawyer       Modified for dynamics_vars
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, imh
      real (r8) zam5, zamda

      if( iord <= 2 ) then
         icd =  1
      else
         icd = -2
      endif
 
      if( jord <= 2 ) then
         jcd =  1
      else
         jcd =  -2
      endif

#if defined( SPMD )
!
! Calculate the ghost region sizes for the SPMD version (tricky stuff)
!
      ng_c = min(abs(jcd ), 2)
      ng_d = min(abs(jord), 3)    ! SJL: number of max ghost latitudes
      ng_d = max(ng_d, 2)
      ng_s = max( ng_c+1, ng_d )
#else
      ng_c = 0
      ng_d = 0                   ! No ghosting necessary for pure SMP runs
      ng_s = 0
#endif

!
! Pole cap area and inverse
      acap = plon*(1.+sine(2)) / dp
      rcap = 1.d0 / acap
 
      imh = plon/2
      if(plon .ne. 2*imh) then
         write(6,*) 'plon must be an even integer'
         stop
      endif
 
! Define logitude at the center of the volume
! i=1, Zamda = -pi
 
      do i=1,imh
         zam5          = ((i-1)-0.5d0) * dl
         cosl5(i)      =  cos(zam5)
         cosl5(i+imh)  = -cosl5(i)
         sinl5(i)      =  sin(zam5)
         sinl5(i+imh)  = -sinl5(i)
         zamda         = (i-1)*dl
         coslon(i)     =  cos(zamda)
         coslon(i+imh) = -coslon(i)
         sinlon(i)     =  sin(zamda)
         sinlon(i+imh) = -sinlon(i)
      enddo

      do j=2,plat-1
         acosp(j) = 1.d0 / cosp(j)
      enddo
      acosp( 1) = rcap * plon
      acosp(plat) = rcap * plon
      return
!EOC
end subroutine dynpkg_init
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  d_split --- find proper value for nsplit if not specified
!
! !INTERFACE:
      subroutine d_split(nsplit_in)
!
! !USES:
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in)   :: nsplit_in  !  Small time steps in dtime

! !DESCRIPTION:
!
!    If nsplit=0 (module variable) then determine a good value 
!    for ns (used in dynpkg) based on resolution and the large-time-step 
!    (pdt). The user may have to set this manually if instability occurs.
!
! !REVISION HISTORY:
!   00.10.19   Lin     Creation
!   01.03.26   Sawyer  ProTeX documentation
!   01.06.10   Sawyer  Modified for dynamics_init framework
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      real (r8)   pdt                       ! Time-step in seconds
                                            ! Negative dt (backward in time
                                            ! integration) is allowed
      real (r8)   dim
      real (r8)   dim0                      ! base dimension
      real (r8)   dt0                       ! base time step
      real (r8)   ns0                       ! base nsplit for base dimension

      parameter ( dim0 = 180.  )
      parameter ( dt0  = 1800. )
      parameter ( ns0  = 4.    )

      if ( nsplit_in == 0 ) then
          pdt = int(dtime)   ! dtime is a variable internal to this module
          dim    = max ( plon, 2*(plat-1) )
          ns = int ( ns0*abs(pdt)*dim/(dt0*dim0) + 0.75 )
          ns = max ( 1, ns )   ! for cases in which dt or dim is too small
      else
          ns = nsplit_in
      endif
      write(6,*) 'Lagrangian time splits (NSPLIT) =', ns

      return
!EOC
      end subroutine d_split
!---------------------------------------------------------------------

end module dynamics_vars

