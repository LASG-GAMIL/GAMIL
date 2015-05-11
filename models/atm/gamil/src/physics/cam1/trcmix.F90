#include <misc.h>
#include <params.h>
subroutine trcmix(lchnk   ,ncol     , &
                  pmid    ,n2o      ,ch4     ,          &
                  cfc11   , cfc12   )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Specify zonal mean mass mixing ratios of CH4, N2O, CFC11 and
! CFC12
! 
! Method: 
! Distributions assume constant mixing ratio in the troposphere
! and a decrease of mixing ratio in the stratosphere. Tropopause
! defined by ptrop. The scale height of the particular trace gas
! depends on latitude. This assumption produces a more realistic
! stratospheric distribution of the various trace gases.
! 
! Author: J. Kiehl
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,    only: get_rlat_all_p
   use physconst,    only: mwdry, mwch4, mwn2o, mwf11, mwf12
   use constituents, only: ch4vmr, n2ovmr, f11vmr, f12vmr

   implicit none

!-----------------------------Arguments---------------------------------
!
! Input
!
   integer, intent(in) :: lchnk                    ! chunk identifier
   integer, intent(in) :: ncol                     ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver)        ! model pressures
!
! Output
!
   real(r8), intent(out) :: n2o(pcols,pver)         ! nitrous oxide mass mixing ratio
   real(r8), intent(out) :: ch4(pcols,pver)         ! methane mass mixing ratio
   real(r8), intent(out) :: cfc11(pcols,pver)       ! cfc11 mass mixing ratio
   real(r8), intent(out) :: cfc12(pcols,pver)       ! cfc12 mass mixing ratio

!
!--------------------------Local Variables------------------------------

   real(r8), parameter :: rmwn2o = mwn2o/mwdry      ! ratio of molecular weight n2o   to dry air
   real(r8), parameter :: rmwch4 = mwch4/mwdry      ! ratio of molecular weight ch4   to dry air
   real(r8), parameter :: rmwf11 = mwf11/mwdry      ! ratio of molecular weight cfc11 to dry air
   real(r8), parameter :: rmwf12 = mwf12/mwdry      ! ratio of molecular weight cfc12 to dry air
!
   integer i                ! longitude loop index
   integer k                ! level index
!
   real(r8) clat(pcols)         ! latitude in radians for columns
   real(r8) coslat(pcols)       ! cosine of latitude
   real(r8) dlat                ! latitude in degrees
   real(r8) ptrop               ! pressure level of tropopause
   real(r8) pratio              ! pressure divided by ptrop
!
   real(r8) xn2o                ! pressure scale height for n2o
   real(r8) xch4                ! pressure scale height for ch4
   real(r8) xcfc11              ! pressure scale height for cfc11
   real(r8) xcfc12              ! pressure scale height for cfc12
!
   real(r8) ch40                ! tropospheric mass mixing ratio for ch4
   real(r8) n2o0                ! tropospheric mass mixing ratio for n2o
   real(r8) cfc110              ! tropospheric mass mixing ratio for cfc11
   real(r8) cfc120              ! tropospheric mass mixing ratio for cfc12
!
!-----------------------------------------------------------------------
!
! get latitudes
!
   call get_rlat_all_p(lchnk, ncol, clat)
   do i = 1, ncol
      coslat(i) = cos(clat(i))
   end do
!
! set tropospheric mass mixing ratios
!
   ch40   = rmwch4 * ch4vmr
   n2o0   = rmwn2o * n2ovmr
   cfc110 = rmwf11 * f11vmr
   cfc120 = rmwf12 * f12vmr

   do i = 1, ncol
      coslat(i) = cos(clat(i))
   end do
!
   do k = 1,pver
      do i = 1,ncol
!
!        set stratospheric scale height factor for gases
         dlat = abs(57.2958 * clat(i))
         if(dlat.le.45.0) then
            xn2o = 0.3478 + 0.00116 * dlat
            xch4 = 0.2353
            xcfc11 = 0.7273 + 0.00606 * dlat
            xcfc12 = 0.4000 + 0.00222 * dlat
         else
            xn2o = 0.4000 + 0.013333 * (dlat - 45)
            xch4 = 0.2353 + 0.0225489 * (dlat - 45)
            xcfc11 = 1.00 + 0.013333 * (dlat - 45)
            xcfc12 = 0.50 + 0.024444 * (dlat - 45)
         end if
!
!        pressure of tropopause
         ptrop = 250.0e2 - 150.0e2*coslat(i)**2.0
!
!        determine output mass mixing ratios
         if (pmid(i,k) >= ptrop) then
            ch4(i,k) = ch40
            n2o(i,k) = n2o0
            cfc11(i,k) = cfc110
            cfc12(i,k) = cfc120
         else
            pratio = pmid(i,k)/ptrop
            ch4(i,k) = ch40 * (pratio)**xch4
            n2o(i,k) = n2o0 * (pratio)**xn2o
            cfc11(i,k) = cfc110 * (pratio)**xcfc11
            cfc12(i,k) = cfc120 * (pratio)**xcfc12
         end if
      end do
   end do
!
   return
!
end subroutine trcmix

