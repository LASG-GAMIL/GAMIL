#include <misc.h>
#include <params.h>

subroutine phys_idealized (phys_state, phys_tend, ztodt, etamid)

!----------------------------------------------------------------------- 
! 
! Purpose: 
!   Driver for idealized physics
! 
! Method: 
! 
! Author: 
!   B. A. Boville
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid       , only: pcols, pver, pverp, begchunk, endchunk
  use phys_grid    , only: get_ncols_p
  use physics_types, only: physics_state, physics_tend
  use buffer       , only: pblht, tpert, qpert
  use comsrf       , only: srfflx_state2d
  use diagnostics  , only: diag_dynvar
  use geopotential , only: geopotential_t
  use physconst,     only: zvir, rair, gravit
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  real(r8), intent(in) :: etamid(pver)     ! vertical coords at midpoints 
  real(r8), intent(in) :: ztodt            ! time step
!
! Output arguments
!
  type(physics_tend ), intent(out  ), dimension(begchunk:endchunk) :: phys_tend
!
!---------------------------Local workspace-----------------------------
  integer  :: i,k,lchnk      ! indices
  integer  :: ncol             ! number of columns
  real(r8) rpdel(pcols,pver)   ! 1./(pintm1(k+1)-pintm1(k))
!-----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE (I, K, LCHNK, NCOL, RPDEL)
  do lchnk=begchunk, endchunk
     ncol = get_ncols_p(lchnk)
!
! Dump dynamics variables to H.T.
!
     do k=1,pver
        do i=1,ncol
           rpdel(i,k) = 1./phys_state(lchnk)%pdel(i,k)
        end do
     end do
     call geopotential_t(                                                             &
          phys_state(lchnk)%lnpint  , phys_state(lchnk)%lnpmid, phys_state(lchnk)%pint    , &
          phys_state(lchnk)%pmid    , phys_state(lchnk)%pdel  , rpdel, phys_state(lchnk)%t, &
          phys_state(lchnk)%q(:,:,1), rair , gravit, zvir , phys_state(lchnk)%zi        , &
          phys_state(lchnk)%zm      , ncol )

     call diag_dynvar (lchnk, ncol, phys_state(lchnk))
!
! Set tendencies to 0
!
     do k=1,pver
        do i=1,ncol
           phys_tend(lchnk)%dTdt(i,k) = 0.
           phys_tend(lchnk)%dudt(i,k) = 0.
           phys_tend(lchnk)%dvdt(i,k) = 0.
        end do
     end do
     do i=1,ncol
        phys_tend(lchnk)%flx_net(i) = 0.
     end do
!
     call t_startf('tphysidl')
     call tphysidl ( &
          ztodt, srfflx_state2d(lchnk)%wsx, srfflx_state2d(lchnk)%wsy, &
             etamid, phys_state(lchnk), phys_tend(lchnk) )
     call t_stopf('tphysidl')

  end do                    ! chunk loop
  return	
end subroutine phys_idealized
