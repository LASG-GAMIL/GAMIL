#include <misc.h>
module fv_prints
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: fv_prints --- print maxima and minima of dycore varibles
!
! !PUBLIC MEMBER FUNCTIONS:
      PUBLIC     fv_out
!
! !DESCRIPTION:
!
!   This module provides basic utilities to evaluate the dynamics state
!
! !REVISION HISTORY:
!   00.08.01   Lin     Creation
!   01.01.05   Boville Modifications
!   01.03.26   Sawyer  Added ProTex documentation
!
!EOP
!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!BOP
! !IROUTINE: fv_out --- Write out maxima and minima of dynamics state
!
! !INTERFACE: 
  subroutine  fv_out( im,     jm,    km,     jfirst,    jlast,             &
                      ng,     kfirst,klast,  pk,        pt,                &
                      ptop,   ps,    q3,     nc,        nq,                &
                      delp,   pe, surf_state, phys_state, ncdate,          &
                      ncsec, full_phys  )

! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use dynamics_vars,  only: gw, cosp
    use ppgrid,         only: begchunk, endchunk, pcols, pver
    use phys_grid,      only: gather_chunk_to_field, get_ncols_p
    use physics_types,  only: physics_state
    use comsrf,         only: surface_state
    use pmgrid,         only: iam, plat, plon, strip3zaty, strip2d,  &
                              myid_y, myid_z
#if defined( SPMD )
    use spmd_dyn,       only: comm_y, comm_z
    use parutilitiesmodule, only : pargatherreal
#endif
    implicit none

! !INPUT PARAMETERS:
    integer im                          ! Total longitudes
    integer jm                          ! Total latitudes
    integer km                          ! Total levels
    integer jfirst                      ! First latitude on this PE
    integer jlast                       ! Last latitude on this PE
    integer ng                          ! Latitude ghost width (D-grid)
    integer kfirst                      ! First level on this PE
    integer klast                       ! Last level on this PE
    integer nc, nq                      ! No. of non-advected and adv. tracers
    integer ncdate                      ! Date
    integer ncsec                       ! Time

    real(r8) ptop                       ! Pressure at top
    real(r8) ps(im,jfirst:jlast)        ! Surface pressure
    real(r8) pk(im,jfirst:jlast,kfirst:klast+1)   ! Pe**kappa
    real(r8) pt(im,jfirst-ng:jlast+ng,kfirst:klast)     ! Potential temp.
    real(r8) delp(im,jfirst:jlast,kfirst:klast)   ! Layer thickness (pint(k+1) - pint(k))
    real(r8)   q3(im,jfirst-ng:jlast+ng,kfirst:klast,nc)! Tracers
    real(r8)   pe(im,kfirst:klast+1,jfirst:jlast) ! Edge pressure

    type(surface_state), intent(in), dimension(begchunk:endchunk) :: surf_state

    type(physics_state), intent(in), dimension(begchunk:endchunk) :: phys_state
    logical full_phys                   ! Full physics on?

!
! !DESCRIPTION:
!
!   Determine maxima and minima of dynamics state and write them out
!
! !REVISION HISTORY:
!   00.08.01   Lin     Creation
!   01.01.05   Boville Modifications
!   01.03.26   Sawyer  Added ProTex documentation
!   01.06.27   Mirin   Converted to 2D yz decomposition
!   01.12.18   Mirin   Calculate average height (htsum) metric
!   02.02.13   Eaton   Pass precc and precl via surface_state type
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    integer i, j, k, ic, nj, lchnk, nck, ncol
    real(r8), dimension(begchunk:endchunk)    :: pmax, tmax, umax, vmax, wmax
    real(r8), dimension(begchunk:endchunk)    :: pmin, tmin, umin, vmin, wmin
    real(r8), dimension(pcols,begchunk:endchunk) :: precc    ! convective precip rate
    real(r8), dimension(pcols,begchunk:endchunk) :: precl    ! large-scale precip rate
    real(r8), dimension(begchunk:endchunk)    :: preccmax, preclmax
    real(r8), dimension(begchunk:endchunk)    :: preccmin, preclmin
    real(r8), dimension(jfirst:jlast,nc) :: qmax
    real(r8), dimension(jfirst:jlast,nc) :: qmin
    real(r8) fac, precmax, precmin
    real(r8) pcon, pls
    real(r8) qtmp(im,kfirst:klast,jfirst:jlast)
    real(r8) p1, p2,dtmp, apcon, htsum
    real(r8), dimension(plon,plat) :: dfield
    real(r8), allocatable :: htgeo(:,:,:),htgeoz(:,:,:)
    real(r8), allocatable :: htg(:,:),htgy(:,:)
    

    integer n, nhmsf
! statement function for hour minutes seconds of day
    nhmsf(n)  = n/3600*10000 + mod(n,3600 )/ 60*100 + mod(n, 60)

    if (iam .eq. 0) then
       write(6,*) ' '
       write(6,*) nhmsf(ncsec), ncdate
    endif
!
! Check total air and dry air mass.

    call dryairm( im,     jm,    km,    jfirst, jlast,     &
                  ng,     kfirst,klast,                    &
                  .true., ptop,  ps,    q3,     nc,        &
                  nq,     delp,  pe,    .true.)

!$omp parallel do private(lchnk, ncol)
    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       pmax(lchnk) = maxval(phys_state(lchnk)%ps(1:ncol))
       pmin(lchnk) = minval(phys_state(lchnk)%ps(1:ncol))
       tmax(lchnk) = maxval(phys_state(lchnk)%t(1:ncol,1:pver))
       tmin(lchnk) = minval(phys_state(lchnk)%t(1:ncol,1:pver))
       umax(lchnk) = maxval(phys_state(lchnk)%u(1:ncol,1:pver))
       umin(lchnk) = minval(phys_state(lchnk)%u(1:ncol,1:pver))
       vmax(lchnk) = maxval(phys_state(lchnk)%v(1:ncol,1:pver))
       vmin(lchnk) = minval(phys_state(lchnk)%v(1:ncol,1:pver))
       wmax(lchnk) = maxval(phys_state(lchnk)%omega(1:ncol,1:pver))
       wmin(lchnk) = minval(phys_state(lchnk)%omega(1:ncol,1:pver))
    end do

    nck = endchunk - begchunk + 1
    call pmaxmin2('PS',         pmin, pmax, nck, 0.01)
    call pmaxmin2('U ',         umin, umax, nck, 1.)
    call pmaxmin2('V ',         vmin, vmax, nck, 1.)
    call pmaxmin2('T ',         tmin, tmax, nck, 1.)
    call pmaxmin2('W (mb/day)', wmin, wmax, nck, 864.)

    nj = jlast - jfirst + 1
    do ic=1,nc
!$omp parallel do private(i, j, k)
       do j=jfirst,jlast
          do k=kfirst,klast
             do i=1,im
                qtmp(i,k,j) = q3(i,j,k,ic)
             enddo
          enddo
       enddo
    call pmaxmin('Q3', qtmp, p1, p2, im*(klast-kfirst+1), nj, 1.)
    end do

    allocate (htgeoz(im,jfirst:jlast,km))
    allocate (htgeo(im,jfirst:jlast,kfirst:klast))
    allocate (htgy(im,jm))
    allocate (htg(im,jfirst:jlast))
    apcon = 1./9.80616
!$omp parallel do private(i, j, k)
    do k=kfirst,klast
      do j=jfirst,jlast
        do i=1,im
          htgeo(i,j,k) = apcon * pt(i,j,k) * (pk(i,j,k+1)-pk(i,j,k))
        enddo
      enddo
    enddo
#if defined( SPMD )
    call pargatherreal(comm_z, 0, htgeo, strip3zaty, htgeoz) 
#else
!$omp parallel do private(i, j, k)
    do k=1,km
      do j=jfirst,jlast
        do i=1,im
          htgeoz(i,j,k) = htgeo(i,j,k)
        enddo
      enddo
    enddo
#endif
    if (myid_z .eq. 0) then
!$omp parallel do private(i, j, k)
       do j=jfirst,jlast
         do i=1,im
           htg(i,j) = 0.
         enddo
         do k=1,km
           do i=1,im
             htg(i,j) = htg(i,j) + htgeoz(i,j,k)
           enddo
         enddo
       enddo
#if defined( SPMD )
       call pargatherreal(comm_y, 0, htg, strip2d, htgy)
#else
       do j=1,jm
          do i=1,im
             htgy(i,j) = htg(i,j)
          enddo
       enddo
#endif
       if (myid_y .eq. 0) then
          htsum = 0.
          do j=1,jm
            do i=1,im
              htsum = htsum + htgy(i,j)*cosp(j)
            enddo
          enddo
          htsum = htsum / (2.*im)
          print *, 'Average Height (geopotential units) = ', htsum
       endif
    endif
    deallocate (htgeoz)
    deallocate (htgeo)
    deallocate (htgy)
    deallocate (htg)

    if ( .not. full_phys ) return

! Global means:

    fac = 86400000.                     ! convert to mm/day

!$omp parallel do private(lchnk, ncol)
    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       precc(:,lchnk) = surf_state(lchnk)%precc(:)
       precl(:,lchnk) = surf_state(lchnk)%precl(:)
       preccmax(lchnk) = maxval(precc(1:ncol,lchnk))
       preccmin(lchnk) = minval(precc(1:ncol,lchnk))
       preclmax(lchnk) = maxval(precl(1:ncol,lchnk))
       preclmin(lchnk) = minval(precl(1:ncol,lchnk))
    end do

    nck = endchunk - begchunk + 1
    call pmaxmin2('PRECC', preccmin, preccmax, nck, fac)
    call pmaxmin2('PRECL', preclmin, preclmax, nck, fac)

    call gather_chunk_to_field(1,1,1,plon,precc,dfield)
    if (iam .eq. 0) then
       pcon = 0.0
       do j=1,plat
          dtmp = dfield(1,j)
          do i=2,plon
             dtmp = dtmp + dfield(i,j)
          enddo
          pcon = pcon + dtmp*gw(j)
       enddo
       pcon = pcon / (2*plat)
    endif

    call gather_chunk_to_field(1,1,1,plon,precl,dfield)
    if (iam .eq. 0) then
       pls = 0.0
       do j=1,plat
          dtmp = dfield(1,j)
          do i=2,plon
             dtmp = dtmp + dfield(i,j)
          enddo
          pls = pls + dtmp*gw(j)
       enddo
       pls = pls / (2*plat)
    endif

    if (iam .eq. 0) then
       pcon = pcon * fac
       pls  = pls  * fac
       write(6,*) 'Total precp=',pcon+pls,' CON=', pcon,' LS=',pls
       write(6,*) ' '
    endif

!EOC
  end subroutine fv_out
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pmaxmin --- Find and print the maxima and minima of a field
!
! !INTERFACE: 
  subroutine pmaxmin( qname, a, pmin, pmax, im, jm, fac )

! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
#if defined( SPMD )
#define CPP_PRT_PREFIX  if(gid.eq.0)
    use parutilitiesmodule, only : commglobal,gid,maxop,parcollective
#else
#define CPP_PRT_PREFIX
#endif
    implicit none

! !INPUT PARAMETERS:
    character*(*)  qname             ! Name of field
    integer  im                      ! Total longitudes
    integer  jm                      ! Total latitudes
    real(r8) a(im,jm)                ! 2D field
    real(r8) fac                     ! multiplication factor

! !OUTPUT PARAMETERS:
    real(r8) pmax                    ! Field maximum
    real(r8) pmin                    ! Field minimum

! !DESCRIPTION:
!
!   Parallelized utility routine for computing/printing global 
!   max/min from input lists of max/min's (usually for each latitude).  
! 
! !REVISION HISTORY:
!   00.03.01   Lin     Creation
!   00.05.01   Mirin   Coalesce variables to minimize collective ops
!   01.08.05   Sawyer  Modified to use parcollective
!   01.03.26   Sawyer  Added ProTex documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

    integer  i, j
    real(r8) qmin(jm), qmax(jm)
    real(r8) pm(2)

!$omp  parallel do default(shared) private(i,j, pmax, pmin)

    do j=1,jm
       pmax = a(1,j)
       pmin = a(1,j)
       do i=2,im
          pmax = max(pmax, a(i,j))
          pmin = min(pmin, a(i,j))
       enddo
       qmax(j) = pmax
       qmin(j) = pmin
    enddo
!
! Now find max/min of qmax/qmin
!
    pmax = qmax(1)
    pmin = qmin(1)
    do j=2,jm
       pmax = max(pmax, qmax(j))
       pmin = min(pmin, qmin(j))
    enddo

#if defined( SPMD )
    pm(1) = pmax
    pm(2) = -pmin
    call parcollective( commglobal, maxop, 2, pm )
    pmax = pm(1)
    pmin = -pm(2)
#endif

    CPP_PRT_PREFIX write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

    return
!EOC
  end subroutine pmaxmin
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pmaxmin2 --- Find and print the maxima and minima of 1-D array
!
! !INTERFACE: 
  subroutine pmaxmin2( qname, qmin, qmax, nj, fac )

! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
#if defined( SPMD )
#define CPP_PRT_PREFIX  if(gid.eq.0)
    use parutilitiesmodule, only : commglobal,gid,maxop,parcollective
#else
#define CPP_PRT_PREFIX
#endif
    implicit none

! !INPUT PARAMETERS:
    character*(*)  qname
    integer nj
    real(r8), intent(in), dimension(nj) :: qmax, qmin      ! Fields
    real(r8) fac                     ! multiplication factor

! !DESCRIPTION:
!
!   Parallelized utility routine for computing/printing global max/min from 
!   input lists of max/min's (usually for each latitude). The primary purpose 
!   is to allow for the original array and the input max/min arrays to be 
!   distributed across nodes.
! 
! !REVISION HISTORY:
!   00.10.01   Lin     Creation from pmaxmin
!   01.03.26   Sawyer  Added ProTex documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    real(r8) pm(2)
    real(r8) pmin, pmax

    pmax = maxval(qmax)
    pmin = minval(qmin)

#if defined( SPMD )
    pm(1) = pmax
    pm(2) = -pmin
    call parcollective( commglobal, maxop, 2, pm )
    pmax = pm(1)
    pmin = -pm(2)
#endif

    CPP_PRT_PREFIX write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

    return
!EOC
  end subroutine pmaxmin2
!-----------------------------------------------------------------------

end module fv_prints

