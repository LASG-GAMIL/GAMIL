#include <misc.h>
#include <params.h>
subroutine cldint(lchnk   ,ncol    , &
                  pmid    ,t       ,q       ,pint    ,piln    , &
                  pmln    ,tvm     ,zi      ,cld     ,clwp    , &
                  emis    ,effcld  ,landfrac,rel     ,rei     , &
                  fice    ,pdel    ,tpw     ,hl      ,ps      , &
                  nmxrgn  ,pmxrgn  ,clwp2   )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface routine for cloud fraction evaluation
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Kiehl
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use physconst, only: gravit

   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver)     ! midpoint pressures
   real(r8), intent(in) :: pint(pcols,pverp)     ! midpoint pressures
   real(r8), intent(in) :: pdel(pcols,pver)     ! pressure depth of layer
   real(r8), intent(in) :: t(pcols,pver)        ! temperature
   real(r8), intent(in) :: q(pcols,pver)        ! specific humidity
   real(r8), intent(in) :: clwp2(pcols,pver)    ! prognostic cloud liquid water path
   real(r8), intent(in) :: piln(pcols,pverp)    ! log of interface pressures
   real(r8), intent(in) :: pmln(pcols,pver)     ! log of midpoint pressures
   real(r8), intent(in) :: tvm(pcols,pver)      ! virtual temperature
   real(r8), intent(in) :: zi(pcols,pverp)      ! height of interfaces (above surface)
   real(r8), intent(in) :: ps(pcols)            ! surface pressure
!
! Output arguments
!
   real(r8), intent(in) :: cld(pcols,pver)       ! cloud fraction
   real(r8), intent(out) :: clwp(pcols,pver)     ! cloud liquid water path
   real(r8), intent(out) :: emis(pcols,pver)     ! cloud emissivity
   real(r8), intent(out) :: effcld(pcols,pver)   ! effective cloud=cld*emis
   real(r8), intent(in)  :: landfrac(pcols)      ! Land fraction
   real(r8), intent(out) :: rel(pcols,pver)      ! effective drop radius (microns)
   real(r8), intent(out) :: rei(pcols,pver)      ! ice effective drop size (microns)
   real(r8), intent(out) :: fice(pcols,pver)     ! fractional ice content within cloud
   real(r8), intent(out) :: tpw(pcols)           ! total precipitable water
   real(r8), intent(out) :: hl(pcols)            ! liquid water scale height
   real(r8), intent(out) :: pmxrgn(pcols,pverp) ! Maximum values of pressure for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pressure for
!    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc
   integer nmxrgn(pcols)        ! Number of maximally overlapped regions
!
!---------------------------Local workspace-----------------------------
!
   integer i,k               ! longitude,level indices

   real(r8) rgrav                ! inverse gravitational acceleration
!
!-----------------------------------------------------------------------
!
! Cloud liquid water path
! Begin by diagnosing total preciptable water in column (in mm)
!
   do i=1,ncol
      tpw(i) = 0.0
   end do
   rgrav = 1.0/gravit
   do k=1,pver
      do i=1,ncol
         tpw(i) = tpw(i) + pdel(i,k)*q(i,k)*rgrav
      end do
   end do
!
   call cldclw(lchnk   ,ncol    ,zi      ,clwp    ,tpw     ,hl      )
!
! Cloud particle size and fraction of ice
!
   call cldefr(lchnk   ,ncol    , &
               landfrac,t       ,rel     ,rei     ,fice    , &
               ps      ,pmid    )
!
! Cloud emissivity.
!
   call cldems(lchnk   ,ncol    ,clwp2   ,fice    ,rei     ,emis    )
!
! Effective cloud cover
!
   do k=1,pver
      do i=1,ncol
         effcld(i,k) = cld(i,k)*emis(i,k)
      end do
   end do
!
! Determine parameters for maximum/random overlap
!
   call cldovrlap(lchnk   ,ncol    ,pint    ,cld     ,nmxrgn  ,pmxrgn  )
!
   return
end subroutine cldint


