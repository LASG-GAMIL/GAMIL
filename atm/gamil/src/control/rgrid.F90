#include <misc.h>

module rgrid

    use pmgrid, only: plat, platd
    use pspect, only: pmmax, pmax
    use infnan, only: bigint

    implicit none

    integer :: nlon(plat)        = bigint ! num longitudes per latitude
    integer :: nlonex(platd)     = bigint ! num longitudes per lat (extended grid)
    integer :: beglatpair(pmmax) = bigint
    integer :: nmmax(plat/2)     = bigint
    integer :: wnummax(plat)     = bigint ! cutoff Fourier wavenumber
#ifdef PVP
    integer :: nmreduced(pmax,plat/2) = bigint
#endif

end module rgrid
