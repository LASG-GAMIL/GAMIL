#include <misc.h>
#include <params.h>
module ramp_so4_mod
!----------------------------------------------------------------------- 
! 
! Purpose: 
!
! Module to handle the SO4 ramping.
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   implicit none

   save

! Time 

   integer, private, parameter:: ntim = 130    ! Number of SO4 sample times stored

! Ouput type declaration

   character*64, private, parameter:: ramp_type =  &
   'RAMP_SO4: using ramp_so4 data'
   logical, private:: ramp_write = .true.      ! Flag to write out information on ramping

! Input data values

!
! yearly data values
!
   integer, private, parameter:: yrdata(ntim) =  (/ &
           1870  ,1871  ,1872  ,1873  ,1874  , &
           1875  ,1876  ,1877  ,1878  ,1879  , &
           1880  ,1881  ,1882  ,1883  ,1884  , &
           1885  ,1886  ,1887  ,1888  ,1889  , &
           1890  ,1891  ,1892  ,1893  ,1894  , &
           1895  ,1896  ,1897  ,1898  ,1899  , &
           1900  ,1901  ,1902  ,1903  ,1904  , &
           1905  ,1906  ,1907  ,1908  ,1909  , &
           1910  ,1911  ,1912  ,1913  ,1914  , &
           1915  ,1916  ,1917  ,1918  ,1919  , &
           1920  ,1921  ,1922  ,1923  ,1924  , &
           1925  ,1926  ,1927  ,1928  ,1929  , &
           1930  ,1931  ,1932  ,1933  ,1934  , &
           1935  ,1936  ,1937  ,1938  ,1939  , &
           1940  ,1941  ,1942  ,1943  ,1944  , &
           1945  ,1946  ,1947  ,1948  ,1949  , &
           1950  ,1951  ,1952  ,1953  ,1954  , &
           1955  ,1956  ,1957  ,1958  ,1959  , &
           1960  ,1961  ,1962  ,1963  ,1964  , &
           1965  ,1966  ,1967  ,1968  ,1969  , &
           1970  ,1971  ,1972  ,1973  ,1974  , &
           1975  ,1976  ,1977  ,1978  ,1979  , &
           1980  ,1981  ,1982  ,1983  ,1984  , &
           1985  ,1986  ,1987  ,1988  ,1989  , &
           1990  ,1991  ,1992  ,1993  ,1994  , &
           1995  ,1996  ,1997  ,1998  ,1999  /)
!
! input Global sulfer emissisions (Tg S/yr)
!
   real(r8), private, parameter ::    semis(ntim) = (/ &
           2.18     ,2.36     ,2.54     ,2.72     ,2.9      , &
           3.08     ,3.26     ,3.44     ,3.62     ,3.8      , &
           3.98     ,4.222    ,4.464    ,4.706    ,4.948    , &
           5.19     ,5.432    ,5.674    ,5.916    ,6.158    , &
           6.4      ,6.761    ,7.122    ,7.483    ,7.844    , &
           8.205    ,8.566    ,8.927    ,9.288    ,9.649    , &
           10.01    ,10.545   ,11.08    ,11.615   ,12.15    , &
           12.685   ,13.22    ,13.755   ,14.29    ,14.825   , &
           15.36    ,15.601   ,15.842   ,16.083   ,16.324   , &
           16.565   ,16.806   ,17.047   ,17.288   ,17.529   , &
           17.77    ,18.018   ,18.266   ,18.514   ,18.762   , &
           19.01    ,19.258   ,19.506   ,19.754   ,20.002   , &
           20.25    ,20.525   ,20.8     ,21.075   ,21.35    , &
           21.625   ,21.9     ,22.175   ,22.45    ,22.725   , &
           23.      ,23.528   ,24.056   ,24.584   ,25.112   , &
           25.64    ,26.168   ,26.696   ,27.224   ,27.752   , &
           28.28    ,29.794   ,31.308   ,32.822   ,34.336   , &
           35.85    ,37.364   ,38.878   ,40.392   ,41.906   , &
           43.42    ,45.384   ,47.348   ,49.312   ,51.276   , &
           53.24    ,55.204   ,57.168   ,59.132   ,61.096   , &
           63.06    ,63.821   ,64.582   ,65.343   ,66.104   , &
           66.865   ,67.626   ,68.387   ,69.148   ,69.909   , &
           70.67    ,71.1242  ,71.5785  ,72.0327  ,72.4869  , &
           72.9412  ,73.1671  ,73.393   ,73.6189  ,73.8448  , &
           74.0707  ,74.0707  ,74.8608  ,75.2558  ,75.6509  , &
           76.0459  ,76.2434  ,76.4410  ,76.6385  ,76.8360  /)

contains

subroutine rampnl_so4( year )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize the ramp options that are controlled by namelist input.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: masterproc
   implicit none

#include <ramp.h>

! Input args.
   integer, intent(in) :: year ! Ramped gases fixed at this year
!-----------------------------------------------------------------------
   rampYear_so4 = year
   fixYear_so4 = .false.
   if ( year > 0 ) then
      fixYear_so4 = .true.
      if (masterproc) &
         write(6,*) 'RAMP_SO4: Ramped gases being fixed at year ',rampYear_so4
   end if
   return
end subroutine rampnl_so4

!##############################################################################

subroutine ramp_so4
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Computes scale factor for ramping sulfate mass mixing ratios
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: <Who is primarily responsible for the code> 
! 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use so4bnd
   use time_manager, only: get_curr_date, get_curr_calday
   use pmgrid,       only: masterproc

   implicit none

#include <ramp.h>
!---------------------------Local variables-----------------------------

   real(r8) semis_ref                ! reference value for sulfer emissions (1985)
   parameter (semis_ref = 65.0)  ! Hardwired as per discussion with Byron (1May98)

   integer yrmodel           ! model year
   integer nyrm              ! year index
   integer nyrp              ! year index
   integer :: yr, mon, day   ! components of a date
   integer :: ncdate         ! current date in integer format [yyyymmdd]
   integer :: ncsec          ! current time of day [seconds]

   real(r8) :: calday            ! current calendar day
   real(r8) doymodel             ! model day of year
   real(r8) doydatam             ! day of year for input data yrdata(nyrm)
   real(r8) doydatap             ! day or year for input data yrdata(nyrp)
   real(r8) deltat               ! delta time
   real(r8) fact1, fact2         ! time interpolation factors
!
! ---------------------------------------------------------------------
!
   calday = get_curr_calday()
   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day

   if (ramp_write) then
      if (masterproc) &
         write(6,*) ramp_type
      ramp_write = .false.
   endif
!
! determine index into input data
!
   if ( fixYear_so4 ) then
      yrmodel  = rampYear_so4
   else
      yrmodel  = ncdate/10000
   end if

   nyrm       = yrmodel - yrdata(1) + 1
   nyrp       = nyrm + 1
!
! if current date is before 1870, quit
!
   if (nyrm < 1) then
      write(6,*)'RAMP_SO4: data time index is out of bounds'
      write(6,*)'nyrm = ',nyrm,' nyrp= ',nyrp, ' ncdate= ', ncdate
      call endrun
   endif
!
! if current date later than ntim, quit
! if want to use just use ntim values - uncomment the following lines
! below and comment the call to endrun and previous write
!
   if (nyrp > ntim) then
      write(6,*)'RAMP_SO4: error - current date is past the end of ', &
                ' valid sulfate scale factor data'
      call endrun
!         write(6,*)'RAMP_SO4: using sulfate scale factor for ',yrdata(ntim)
!         call setso4ramp( semis(ntim)/semis_ref )
!         return
   endif
!
! determine time interpolation factors, check sanity
! of interpolation factors to within 32-bit roundoff
! assume that day of year is 1 for all input data
!
   doymodel = yrmodel*365.    + calday
   doydatam = yrdata(nyrm)*365. + 1.
   doydatap = yrdata(nyrp)*365. + 1.
   deltat   = doydatap - doydatam
   fact1    = (doydatap - doymodel)/deltat
   fact2    = (doymodel - doydatam)/deltat

   if (abs(fact1+fact2-1.) > 1.e-6 .or. &
       fact1 > 1.000001 .or. &
       fact1 < -1.e-6 .or. &
       fact2 > 1.000001 .or. &
       fact2 < -1.e-6) then
      write(6,*)'RAMP_SO4: Bad fact1 and/or fact2=',fact1,fact2
      call endrun
   end if
!
! do time interpolation:
!
   call setso4ramp((semis(nyrm)*fact1 + semis(nyrp)*fact2)/semis_ref)

   return
end subroutine ramp_so4

end module ramp_so4_mod
