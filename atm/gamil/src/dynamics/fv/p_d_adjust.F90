#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: p_d_adjust --- complete full physics update
!
! !INTERFACE: 
  subroutine p_d_adjust(pe,   pt,   q3,   q_tmp,   delp,   ps,          &
                        peln, pk,   pkz,  zvir,    cappa,  delpxy,      &
                        pexy, ptop, full_phys  )

! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid, only : iam, npr_y, npr_z, myid_z, plat, plon, plev,     &
                       beglonxy, endlonxy, beglatxy, endlatxy,          &
                       beglat, endlat, beglev, endlev, endlevp
    use dynamics_vars, only : ng_d
    use constituents, only: ppcnst, pcnst
#if defined( SPMD )
    use spmd_dyn, only: comm_z, inter_ijk, inter_ikjp
    use mod_comm, only: mp_send3d, mp_recv3d
    use parutilitiesmodule, only: parcollective2d, sumop
    use redistributemodule, only: redistributetype, redistributestart, &
                                  redistributefinish
#endif
!-----------------------------------------------------------------------
    implicit none

! !INPUT PARAMETERS:
    real(r8), intent(in) :: zvir
    real(r8), intent(in) :: ptop
    real(r8), intent(in) :: cappa
    logical,  intent(in) :: full_phys
    real(r8), intent(in) :: q_tmp(plon,beglev:endlev,beglat:endlat) ! temporary moisture storage

! !INPUT/OUTPUT PARAMETERS:
    real(r8), intent(inout) :: pt(plon,beglat-ng_d:endlat+ng_d,beglev:endlev)
    real(r8), intent(inout) :: q3(plon,beglat-ng_d:endlat+ng_d,beglev:endlev,ppcnst) ! constituents

    real(r8), intent(inout) :: pe(plon,beglev:endlev+1,beglat:endlat) ! interface pres
    real(r8), intent(inout) :: delpxy(beglonxy:endlonxy,beglatxy:endlatxy,plev)
    real(r8), intent(inout) :: pexy(beglonxy:endlonxy,plev+1,beglatxy:endlatxy)

! !OUTPUT PARAMETERS
    real(r8), intent(out) :: ps(plon,beglat:endlat)                 ! surf. press
    real(r8), intent(out) :: delp(plon,beglat:endlat,beglev:endlev)   ! pressure difference
    real(r8), intent(out) :: peln(plon,beglev:endlev+1,beglat:endlat) ! interface pres
    real(r8), intent(out) :: pkz(plon,beglat:endlat,beglev:endlev)    ! Layer-mean value of PK
    real(r8), intent(out) :: pk(plon,beglat:endlat,beglev:endlev+1)  ! PE**cappa

! !DESCRIPTION:
!
!   Complete adjustment of quantities after physics update
!
! !REVISION HISTORY:
!   00.06.01   Grant?     Creation
!   01.06.08   AAM        Created from p_d_coupling
!   02.04.24   WS         New mod_comm interface
!   02.05.01   WS         Fix of S.-J. and Phil to peln, pk update
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
    integer :: i, k, m, j ! indices
    real(r8) dp2(plon,beglev:endlev)
    real(r8) fdq(plon)
    real(r8) pekp(plon,beglev:endlevp,beglat:endlat)
    integer :: incount, outcount
!---------------------------End Local workspace-------------------------

!
! ----------------------------------------------------
! Complete update of dynamics variables
! ----------------------------------------------------
!

    if (full_phys) then

!$omp parallel do private(i, j, k, m, fdq, dp2)

       do j=beglat,endlat
!
! Update interface pressure (pe), and del-pressure modified by gain/loss in water vapor.
! Update q3, pe, delp, ps, peln, pk, and pkz
!
          do k = beglev, endlev
             do i = 1, plon
                dp2(i,k) = pe(i,k+1,j) - pe(i,k,j)
             enddo
          enddo
         
          do k = beglev, endlev

             if (j==1 .or. j==plat) then
!
! Perform average at both poles: all scalars at poles should be single-valued
!
                call xpavg( pt(1,j,k),   plon )
                do m=1, ppcnst
                   call xpavg( q3(1,j,k,m), plon )
                enddo
             endif
            
             do i = 1, plon
                fdq(i) = 1. + q3(i,j,k,1) - q_tmp(i,k,j)
             enddo

             if (j==1 .or. j==plat) then
                call xpavg( fdq,   plon ) ! This is needed for bit-for-bit restartability; but why?
             endif
            
             do i = 1, plon
                q3  (i,j,k,1) = q3(i,j,k,1) / fdq(i)
                dp2 (i,k)     = dp2(i,k) * fdq(i)
                pe  (i,k+1,j) = pe(i,k,j) + dp2(i,k)
!
! Incorrect pe in most subdomains if 2D decomposition, since z-communication
!   is required; this is taken care of below.
!
                delp(i,j,k)   = dp2(i,k)
             enddo
            
!
! Update advected tracers
!
             if ( pcnst > 1 ) then
                do m=2, pcnst            ! advective tracers
                   do i = 1,plon
                      q3(i,j,k,m) = q3(i,j,k,m) / fdq(i)
                   enddo
                enddo
             endif
            
          enddo
!
! Set surface pressure (incorrect if 2D decomposition, since in that case
!   pe is incorrect in most subdomains; this is taken care of below.)
!
          do i = 1, plon
             ps(i,j) = pe(i,endlev+1,j)
          enddo
       end do     ! beglat:endlat loop
!
! Fix pe,ps if nontrivial z decomposition
! Transpose delp, pe - change to better method (16-byte?) later on
!
#if defined( SPMD )
       if (npr_z .gt. 1) then
!
! Transform delp to delpxy
!
          call redistributestart (inter_ijk, .true., delp)
          call redistributefinish(inter_ijk, .true., delpxy)
!
! Compute pexy
!
!$omp parallel do private(i, j)
          do j = beglatxy,endlatxy
             do i = beglonxy, endlonxy
                pexy(i,1,j) = ptop
             enddo
          enddo
!$omp parallel do private(i, j, k)
          do j = beglatxy,endlatxy
             do k = 1, plev
                do i = beglonxy, endlonxy
                   pexy(i,k+1,j) = pexy(i,k,j) + delpxy(i,j,k)
                enddo
             enddo
          enddo
!
! Transform pexy to pe
!
          call redistributestart (inter_ikjp, .false., pexy)
          call redistributefinish(inter_ikjp, .false., pekp)
!$omp parallel do private(i,j,k)
          do j = beglat,endlat
             do k = beglev,endlevp
                do i = 1,plon
                   pe(i,k,j) = pekp(i,k,j)
                enddo
             enddo
          enddo
!
! Send to PE below and receive from PE above
!
          call mp_send3d( iam-npr_y, iam+npr_y, plon, plev+1, plat,           &
                          1, plon, beglev, endlev+1, beglat, endlat,          &
                          1, plon, beglev, beglev, beglat, endlat, pe )
          call mp_recv3d( iam+npr_y, plon, plev+1, plat,                      &
                          1, plon, beglev, endlev+1, beglat, endlat,          &
                          1, plon, endlev+1, endlev+1, beglat, endlat, pe )
!
! Broadcast ps
!
          if (myid_z .eq. npr_z-1) then
             do j=beglat,endlat
                do i=1,plon
                   ps(i,j) = pe(i,endlev+1,j)
                enddo
             enddo
          else
             do j=beglat,endlat
                do i=1,plon
                   ps(i,j) = 0.
                enddo
             enddo
          endif
          call parcollective2d(comm_z, sumop, plon, endlat-beglat+1, ps)
       endif
#endif

    endif

!$omp parallel do private(i, j, k)

    do j=beglat,endlat

!
! Update peln and pk
!
       do k = beglev, endlev+1
          do i = 1, plon
             peln(i,k,j) = log( pe(i,k,j) )
             pk  (i,j,k) = pe(i,k,j) ** cappa
          enddo
       enddo

!
! Update pkz
!
       do k = beglev,endlev
          do i = 1,plon
             pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(cappa*(peln(i,k+1,j)-peln(i,k,j)))
          enddo
       enddo
    enddo     ! beglat:endlat loop

    if (full_phys) then
!
! Calculate virtual potential temperature
!

!$omp parallel do private(i, j, k)

       do j=beglat,endlat
          do k = beglev,endlev
             do i = 1,plon
                pt(i,j,k) = pt(i,j,k)*(1.+zvir*q3(i,j,k,1))/pkz(i,j,k)
             enddo
          enddo
       enddo     ! beglat:endlat loop

    endif

!EOC
  end subroutine p_d_adjust
!-----------------------------------------------------------------------
