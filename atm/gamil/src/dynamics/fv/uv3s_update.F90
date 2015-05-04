#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: uv3s_update --  update u3s, v3s
!
! !INTERFACE:

   subroutine uv3s_update(dua, u3s, dva, v3s, dt5, im, jm,                 &
                          km, jfirst, jlast, ngus, ngun, ngvs, ngvn,       &
                          kfirst, klast)
! !USES:

      use shr_kind_mod, only: r8 => shr_kind_r8

#if defined( SPMD )
      use parutilitiesmodule, only : pargatherreal
      use mod_comm, only : mp_send3d, mp_recv3d
#endif
      use pmgrid, only : myid_y, npr_y, twod_decomp, myid_z, npr_z, iam,   &
                         strip3zatyt4
      use history, only: outfld

      implicit none
! !INPUT PARAMETERS:
      integer, intent(in)  :: im      ! Dimensions longitude
      integer, intent(in)  :: jm      ! Dimensions latitude  (total)
      integer, intent(in)  :: km      ! Dimensions vertical (total)
      integer, intent(in)  :: jfirst  ! latitude strip start
      integer, intent(in)  :: jlast   ! latitude strip finish
      integer, intent(in)  :: ngus    ! ghost latitudes U south
      integer, intent(in)  :: ngun    ! ghost latitudes U north
      integer, intent(in)  :: ngvs    ! ghost latitudes V south
      integer, intent(in)  :: ngvn    ! ghost latitudes V north
      integer, intent(in)  :: kfirst  ! vertical strip start
      integer, intent(in)  :: klast   ! vertical strip finish
      real(r8),intent(in)  :: dua(im,kfirst:klast,jfirst:jlast)    ! dudt on A-grid 
      real(r8),intent(in)  :: dva(im,kfirst:klast,jfirst:jlast)    ! dvdt on A-grid 
      real(r8),intent(in)  :: dt5     ! weighting factor

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: u3s(im,jfirst-ngus:jlast+ngun,kfirst:klast)  ! U-Wind on D Grid
      real(r8), intent(inout) :: v3s(im,jfirst-ngvs:jlast+ngvn,kfirst:klast)  ! V-Wind on D Grid

! !DESCRIPTION:
!
!     This routine performs the update for the N-S staggered u-wind
!       and the E-W staggered v-wind
!
! !REVISION HISTORY:
!    WS   00.12.22 : Creation from d2a3d
!    SJL  01.01.20 : modifications
!    AAM  01.06.08 : Name change; folding in of v3s update and outfld calls
!    WS   02.04.25 : New mod_comm interfaces
!    WS   02.07.04 : Fixed 2D decomposition bug dest/src for mp_send3d
!
!EOP
!-----------------------------------------------------------------------
!BOC

      integer i, j, k

#if defined( SPMD )
   real(r8) duasouth(im,kfirst:klast)
   integer  dest, src
#endif
   real(r8) u3s_tmp(im,kfirst:klast), v3s_tmp(im,kfirst:klast)

#if defined( SPMD )
!
! Transfer dua(:,jlast) to the node directly to the north
!
      dest = iam+1
      src  = iam-1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      if ( mod(iam,npr_y) == 0 ) src = -1
      call mp_send3d( dest, src, im, jm, km,                       &
                      1, im, kfirst, klast, jfirst, jlast,            &
                      1, im, kfirst, klast, jlast, jlast, dua )
      call mp_recv3d( src, im, jm, km,                              &
                      1, im, jfirst-1, jfirst-1, kfirst, klast,       &
                      1, im, jfirst-1, jfirst-1, kfirst, klast, duasouth )
#endif

!$omp parallel do private (i, j, k)

      do k = kfirst, klast

!
! Adjust D-grid winds by interpolating A-grid tendencies.
!

        do j = jfirst+1, jlast
          do i = 1, im
             u3s(i,j,k) = u3s(i,j,k) + dt5*(dua(i,k,j)+dua(i,k,j-1))
          enddo
        enddo

#if defined( SPMD )
        if ( jfirst .gt. 1 ) then
          do i = 1, im
             u3s(i,jfirst,k) = u3s(i,jfirst,k) +                         &
                         dt5*( dua(i,k,jfirst) + duasouth(i,k) )
          enddo
        endif
#endif

        do j = max(jfirst,2), min(jlast,jm-1)
           v3s(1,j,k) = v3s(1,j,k) + dt5*(dva(1,k,j)+dva(im,k,j))
           do i=2,im
              v3s(i,j,k) = v3s(i,j,k) + dt5*(dva(i,k,j)+dva(i-1,k,j))
           enddo
        enddo

      enddo

!$omp parallel do private (i, j, k, u3s_tmp, v3s_tmp)

      do j = jfirst, jlast
         do k = kfirst, klast
            do i = 1, im
               u3s_tmp(i,k) = u3s(i,j,k)
               v3s_tmp(i,k) = v3s(i,j,k)
            enddo
         enddo

         call outfld ('FU      ', dua(1,kfirst,j), im, j )
         call outfld ('FV      ', dva(1,kfirst,j), im, j )
         call outfld ('US      ', u3s_tmp, im, j )
         call outfld ('VS      ', v3s_tmp, im, j )

      enddo

      return
!EOC
      end
!-----------------------------------------------------------------------
