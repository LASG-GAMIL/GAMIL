#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: rayl_fric --- Rayleigh friction
!
! !INTERFACE: 
  subroutine rayl_fric(phys_state, phys_tend, dt, pe11k, pe11kln,     &
                       cpair, cappa, rfac, rayf )

! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid
    use ppgrid, only: pcols
    use phys_grid
    use physics_types, only: physics_state, physics_tend
!-----------------------------------------------------------------------
    implicit none

! !INPUT PARAMETERS:

    real(r8), intent(in) :: dt                             ! time step size
    real(r8), intent(in) :: cpair                          !
    real(r8), intent(in) :: cappa                          !
    real(r8), intent(in) :: rfac(plev)                     ! Rayleigh friction factor
    real(r8), intent(in) :: pe11k(plev+1)                  ! Reference pressure
    real(r8), intent(in) :: pe11kln(plev+1)                ! Reference log
    logical,  intent(in) :: rayf                           ! Rayleigh friction flag

! !INPUT/OUTPUT PARAMETERS:

    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend),  intent(inout), dimension(begchunk:endchunk) :: phys_tend

! !DESCRIPTION:
!
!   Computes Rayleigh friction
!
! !REVISION HISTORY:
!   01.06.06   Mirin     Creation (from p_d_coupling)
!   01.07.13   Mirin     Accommodate multi-2D decomposition
!   01.07.20   Mirin     Use reference values for pe, peln
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
    integer :: nlonl, i, k, lchnk     ! indices
    integer :: ncol                   ! number of columns in current chunk
    integer :: lats(pcols)            ! array of latitude indices
    integer :: lons(pcols)            ! array of longitude indices
    real(r8) rdt                      ! Inverse of time step
    real(r8) dt5                      ! 0.5 * dt
    real(r8) rcp                      ! Inverse of cpair
    real(r8) rtmp
    real(r8) durf, dvrf
!---------------------------End Local workspace-------------------------

    rdt = 1. / dt
    dt5 = 0.5*dt

! -----------------------------
! Perform Rayleigh friction
! -----------------------------
!$omp parallel do private(lchnk, i, k, ncol, rtmp, durf, dvrf, rcp)

       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)
          call get_lat_all_p(lchnk, ncol, lats)

         do k = 1, plev
            rcp = 1. / ( cpair * ( 1. - cappa*pe11k(k)*    &
                        ( pe11kln(k+1)-pe11kln(k) ) /      &
                        ( pe11k(k+1)-pe11k(k) )  ) )

            do i = 1, ncol
               if ( rayf .and.  pe11k(k) < 3000. ) then    ! only above 30 mb
                  rtmp = - rfac(k) / (1. + rfac(k))
! Implicit-in-time
                  durf = phys_state(lchnk)%u(i,k)*rtmp
                  dvrf = phys_state(lchnk)%v(i,k)*rtmp
! -----------------------------
! Update t, u, v
! -----------------------------
                  phys_state(lchnk)%t(i,k) = phys_state(lchnk)%t(i,k)       &
                            - (durf*(phys_state(lchnk)%u(i,k)+0.5*durf)     &
                            +  dvrf*(phys_state(lchnk)%v(i,k)+0.5*dvrf))*rcp
                  phys_state(lchnk)%u(i,k) = phys_state(lchnk)%u(i,k) + durf
                  phys_state(lchnk)%v(i,k) = phys_state(lchnk)%v(i,k) + dvrf

! -----------------------------
! Update (u, v) tendencies
! -----------------------------
                  phys_tend(lchnk)%dudt(i,k) = phys_tend(lchnk)%dudt(i,k) + durf*rdt
                  phys_tend(lchnk)%dvdt(i,k) = phys_tend(lchnk)%dvdt(i,k) + dvrf*rdt
               endif
            enddo
         enddo                 ! k-loop

       enddo                   ! lchnk-loop


!EOC
  end subroutine rayl_fric
!-----------------------------------------------------------------------
