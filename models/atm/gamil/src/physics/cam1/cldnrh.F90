#include <misc.h>
#include <params.h>
subroutine cldnrh(lchnk   ,ncol    , &
                  pmid    ,t       ,q       ,omga    , &
                  cnt     ,cnb     ,cldn    ,clc     ,pdel    , &
                  cmfmc   ,landfrac,snowh   ,concld  ,cldst   , &
                  ts      ,ps      ,zdu     ,ocnfrac ,          &
                  rhdfda  ,rhu00   )  
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to compute cloud, and its response to rh perturbation
! 
!
! Author: W. Lin 
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid

   implicit none
!------------------------------Parameters-------------------------------
   real(r8) pnot                  ! reference pressure
   parameter (pnot = 1.e5)

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                  ! chunk identifier
   integer, intent(in) :: ncol                   ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver)      ! midpoint pressures
   real(r8), intent(in) :: t(pcols,pver)         ! temperature
   real(r8), intent(in) :: q(pcols,pver)         ! specific humidity
   real(r8), intent(in) :: omga(pcols,pver)      ! vertical pressure velocity
   real(r8), intent(in) :: cnt(pcols)            ! top level of convection
   real(r8), intent(in) :: cnb(pcols)            ! bottom level of convection
   real(r8), intent(in) :: cmfmc(pcols,pverp)    ! convective mass flux--m sub c
   real(r8), intent(in) :: snowh(pcols)          ! snow depth (liquid water equivalent)
   real(r8), intent(in) :: pdel(pcols,pver)      ! pressure depth of layer
   real(r8), intent(in) :: landfrac(pcols)       ! Land fraction
   real(r8), intent(in) :: ocnfrac(pcols)        ! Ocean fraction
   real(r8), intent(in) :: ts(pcols)             ! surface temperature
   real(r8), intent(in) :: ps(pcols)             ! surface pressure
   real(r8), intent(in) :: zdu(pcols,pver)       ! detrainment rate from deep convection

!
! Output arguments
!
   real(r8), intent(out) :: cldn(pcols,pver)      ! cloud fraction
   real(r8), intent(out) :: clc(pcols)            ! column convective cloud amount
   real(r8), intent(out) :: cldst(pcols,pver)     ! cloud fraction
   real(r8), intent(out) :: rhdfda(pcols,pver)    ! d_RH/d_cloud_fraction    ====wlin
   real(r8), intent(out) :: rhu00(pcols,pver)     ! RH limit, U00             ====wlin
   real(r8), intent(out) :: concld(pcols,pver)    ! convective cloud cover

!
! local work space
!
   real(r8) relhum(pcols,pver)         ! RH, output to determine drh/da
   real(r8) rhu002(pcols,pver)         ! same as rhu00 but for perturbed rh 
   real(r8) cldn2(pcols,pver)          ! same as cldn but for perturbed rh  
   real(r8) concld2(pcols,pver)        ! same as concld but for perturbed rh 
   real(r8) cldst2(pcols,pver)         ! same as cldst but for perturbed rh 
   real(r8) relhum2(pcols,pver)        ! RH after  perturbation            

   integer i,k                    ! longitude, level indices

!
   call t_startf("cldfrc")
   call cldfrc(lchnk,   ncol,                                &
               pmid,      t,        q,     omga, &
               cnt,     cnb,     cldn,    clc,     pdel,   &
               cmfmc,   landfrac,snowh,   concld,  cldst,    &
               ts,      ps,      zdu,     ocnfrac, rhu00, &
	       relhum,  0  )    

! re-calculate cloud with perturbed rh             add call cldfrc  
   call cldfrc(lchnk,   ncol,                                &
               pmid,       t,      q,      omga, &
               cnt,     cnb,     cldn2,   clc,     pdel,   &
               cmfmc,   landfrac,snowh,   concld2, cldst2,   &
               ts,      ps,        zdu,   ocnfrac, rhu002, &
	       relhum2, 1  )              

   call t_stopf("cldfrc")

! cldfrc does not define layer cloud for model layer at k=1
! so set rhu00(k=1)=2.0 to not calculate cme for this layer

  rhu00(:ncol,1) = 2.0 

! Add following to estimate rhdfda                       

  do k=1,pver
     do i=1,ncol
        if(relhum(i,k) < rhu00(i,k) ) then
              rhdfda(i,k)=0.0
        else if (relhum(i,k) >= 1.0 ) then
              rhdfda(i,k)=0.0
        else
              !under certain circumstances, rh+ cause cld not to changed
              !when at an upper limit, or w/ strong subsidence
              !need to further check whether this if-block is necessary

              if((cldn2(i,k) - cldn(i,k) ) < 1.e-4 ) then
                 rhdfda(i,k) = 0.01*relhum(i,k)*1.e+4   !instead of 0.0
              else
                 rhdfda(i,k)=0.01*relhum(i,k)/(cldn2(i,k)-cldn(i,k))
              endif
        endif
     enddo
  enddo

  return
end subroutine cldnrh
