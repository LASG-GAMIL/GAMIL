#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: dryairm --- Check dry air mass; set to a predefined value if
!                       nlres is false (initialization run)
!
! !INTERFACE:

subroutine dryairm( im,    jm,    km,   jfirst,   jlast,     &
                    ng,    kfirst,klast,                     &
                    moun,  ptop,  ps,   q,       nc,         &
                    nq,    delp,  pe,   nlres )

! !USES:
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid, only: myid_z, npr_z, strip3zaty

#if defined( SPMD )
#define CPP_PRT_PREFIX  if( gid == 0 )
  use spmd_dyn, only: comm_z
  use parutilitiesmodule, only : gid, parcollective, pargatherreal, sumop
#else
#define CPP_PRT_PREFIX
#endif

 implicit   none

 integer, intent(in):: im, jm, km     ! Dimensions
 integer, intent(in):: jfirst, jlast  ! Latitude strip
 integer, intent(in):: ng             ! Ghost latitudes
 integer, intent(in):: kfirst, klast  ! Vertical strip
 integer, intent(in):: nc             ! Total number of tracers         
 integer, intent(in):: nq             ! Number of advective tracers         
 logical, intent(in):: nlres
 logical, intent(in):: moun
 real(r8), intent(in) :: ptop

 real(r8), intent(inout) :: q(im,jfirst-ng:jlast+ng,kfirst:klast,nc) 
 real(r8), intent(inout) :: ps(im,jfirst:jlast)          ! surface pressure
 real(r8), intent(inout) :: delp(im,jfirst:jlast,kfirst:klast)     ! pressure thickness
 real(r8), intent(inout) :: pe(im,kfirst:klast+1,jfirst:jlast)     ! edge pressure

! !DESCRIPTION:
!  Perform adjustment of the total dry-air-mass while preserving total
!  tracer mass
!  Developer: S.-J. Lin, Aug 2000
!
! !REVISION HISTORY:
!   AAM   01.06.27       Assure agreement thru roundoff for 2D decomp.
!
!EOP
!---------------------------------------------------------------------
!BOC

! Use work arrays psdk/psdkg to assure identical answers through roundoff
!    for different z decompositions

      real(r8), allocatable :: psdk(:,:,:)     ! local work array
      real(r8), allocatable :: psdkg(:,:,:)    ! global work array
      real(r8)    psd(im,jfirst:jlast)    ! dry surface pressure
      real(r8)   drym                     ! global mean dry air mass in pascals

#if defined ( NAVY10 )
      parameter (drym = 98222.)           ! For US NAVY 10-min terrain
#else
      parameter (drym = 98288.)           ! For USGS terrain
#endif

      integer  i, j, k, ic
      real(r8) psm0, psm1
      real(r8) psdry
      real(r8) dpd

! Check global maximum/minimum

    call gmean ( im,   jm,    jfirst,    jlast,   ps(1,jfirst), psm0 )

    allocate (psdk(im,jfirst:jlast,kfirst:klast))
    allocate (psdkg(im,jfirst:jlast,km))

!$omp  parallel do private(i,j,k)
    do k=kfirst,klast
       do j=jfirst,jlast
          do i=1,im
             psdk(i,j,k) = 0.
          enddo
       enddo
    enddo
!$omp  parallel do private(i,j,k)
    do k=1,km
       do j=jfirst,jlast
          do i=1,im
             psdkg(i,j,k) = 0.
          enddo
       enddo
    enddo

    if (kfirst .eq. 1) then

!$omp  parallel do private(i,j)

       do j=jfirst,jlast
          do i=1,im
             psdk(i,j,1) = ptop
          enddo
       enddo

    endif

    if( nq .ne. 0 ) then

!$omp  parallel do private(i,j,k)

       do j=jfirst,jlast
          do k=kfirst,klast
             do i=1,im
                psdk(i,j,k) = psdk(i,j,k) + (1.-q(i,j,k,1))*(pe(i,k+1,j)-pe(i,k,j))
             enddo
          enddo
       enddo

    else

!$omp  parallel do private(i,j,k)

       do j=jfirst,jlast
          do k=kfirst,klast
             do i=1,im
                psdk(i,j,k) = psdk(i,j,k) +  pe(i,k+1,j) - pe(i,k,j)
             enddo
          enddo
       enddo

    endif

    if (npr_z .gt. 1) then
#if defined (SPMD)
       call pargatherreal(comm_z, 0, psdk, strip3zaty, psdkg)
#endif
    else
!$omp  parallel do private(i,j,k)
       do k=kfirst,klast
          do j=jfirst,jlast
             do i=1,im
                psdkg(i,j,k) = psdk(i,j,k)
             enddo
          enddo
       enddo
    endif

!$omp  parallel do private(i,j)
    do j=jfirst,jlast
       do i=1,im
          psd(i,j) = 0.
       enddo
    enddo

    if (myid_z .eq. 0) then

!$omp  parallel do private(i,j,k)

       do j=jfirst,jlast
          do k=1,km
             do i=1,im
                psd(i,j) = psd(i,j) + psdkg(i,j,k)
             enddo
          enddo
       enddo

    endif

#if defined (SPMD)
    if (npr_z .gt. 1) then
       call parcollective(comm_z, sumop, im, jlast-jfirst+1, psd)
    endif
#endif

    call gmean( im,  jm,  jfirst,  jlast,  psd(1,jfirst), psdry )
 
 CPP_PRT_PREFIX write(6,*) 'Total Mass=', 0.01*psm0, '(mb), Dry Mass=', 0.01*psdry, '(mb)'
 CPP_PRT_PREFIX write(6,*) 'Total Precipitable Water =', (psm0-psdry)/9.80616, '(kg/m**2)'

    deallocate (psdk)
    deallocate (psdkg)

    if( nlres ) return

    if(moun) then
       dpd = drym - psdry
    else
       dpd = 1000.*100. - psdry
    endif
 CPP_PRT_PREFIX write(6,*) 'dry mass to be added =', 0.01*dpd

    if (klast .eq. km) then

!$omp  parallel do private(i, j, ic)

       do j=jfirst,jlast

          do ic=1,nq
             do i=1,im
                q(i,j,km,ic) =  q(i,j,km,ic)*delp(i,j,km)/(delp(i,j,km)+dpd)
             enddo
          enddo

! Adjust the lowest Lagrangian layer
          do i=1,im
             delp(i,j,km) = delp(i,j,km) + dpd
             pe(i,km+1,j) = pe(i,km,j) + delp(i,j,km)
             ps(i,j) = pe(i,km+1,j)
          enddo
       enddo
    endif

    if (npr_z .gt. 1) then

       if (myid_z .ne. npr_z-1) then
!$omp  parallel do private(i,j)
          do j=jfirst,jlast
             do i=1,im
                ps(i,j) = 0.
             enddo
          enddo
       endif

#if defined (SPMD)
       call parcollective(comm_z, sumop, im, jlast-jfirst+1, ps)
#endif

    endif

    call gmean(im, jm, jfirst, jlast, ps(1,jfirst), psm1)

 CPP_PRT_PREFIX write(6,*) 'Total moist surface pressure after adjustment (mb) = ',0.01*psm1 

 return

!EOC
end subroutine dryairm
!---------------------------------------------------------------------
