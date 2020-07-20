#include <misc.h>
#include <params.h>

subroutine dynpkg(dtdy, nseq)

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use prognostics
  use commap
  use stdatm  
  use comfm1
  use comhd
  use fspan

  implicit none

#include <comhyb.h>
#include <comctl.h>

  real(r8), intent(in)  :: dtdy ! timestep size ( dyn core )
  integer , intent(in)  :: nseq

  real(r8), allocatable :: tkk(:,:,:)
  real(r8), allocatable :: ukk(:,:,:)
  real(r8), allocatable :: vkk(:,:,:)

  integer i, j, k, begj

  allocate(tkk(beglonex:endlonex,beglatexdyn:endlatexdyn,plev))
  allocate(ukk(beglonex:endlonex,beglatexdyn:endlatexdyn,plev))
  allocate(vkk(beglonex:endlonex,beglatexdyn:endlatexdyn,plev))

  begj = beglatexdyn

  call t_startf('dyfram')

  call dyfram2(nseq, dtdy, itime                                  &
              ,u, v, t, q, ws, pes, wpa, ghs, ghi, ply, tb        &
              ,su, sv, st, sq                                     &
              ,pmtop, sig, sigl, dsig                             &
              ,tbb, hbb, cbb, dcbb, psb, tsb                      &
              ,dy, wtgu(begj),wtgv(begj)                          &
              ,dx, sinu, sinv, oux, ouy, ovx, ovy, ff, cur        &
              ,mm1, mp1, mm2, mp2, mm3, mp3, mdj)

  call t_stopf('dyfram')

  call t_startf('hdifus')

  if (.not. aqua_planet)  then
!$omp parallel do collapse(3)
    do k = 1, plev
      do j = jbeg0, jend0
        do i = beglonex, endlonex
          ukk(i,j,k) = u(i,j,k)
          vkk(i,j,k) = v(i,j,k)
          tkk(i,j,k) = t(i,j,k)
          qk (i,j,k) = q(i,j,k)
        end do
      end do
    end do

    call hdifus(u, v, t, q, frdt, frds, frdu, frdv, frdp, tb, ply, dxvpn, dxvps)


!$omp parallel do collapse(3)
    do k = 1, plev
      do j = jbeg0, jend0
        do i = beglonex, endlonex
          su(i,j,k) = (u(i,j,k) - ukk(i,j,k)) / dthdfs
          sv(i,j,k) = (v(i,j,k) - vkk(i,j,k)) / dthdfs
          st(i,j,k) = (t(i,j,k) - tkk(i,j,k)) / dthdfs
          u (i,j,k) = ukk(i,j,k)
          v (i,j,k) = vkk(i,j,k)
          t (i,j,k) = tkk(i,j,k)
          q (i,j,k) = qk (i,j,k)
        end do
      end do
    end do
  else
    su = 0.0
    sv = 0.0
    st = 0.0
  end if

  deallocate(ukk)
  deallocate(vkk)
  deallocate(tkk)

  call t_stopf('hdifus')

end subroutine dynpkg
