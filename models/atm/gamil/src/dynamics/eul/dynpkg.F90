#include <misc.h>
#include <params.h>

!! (2003.05.18)
!! (2003.11.13)
!! (2003.12.01)
!! (2005.01.28)
!!------------------------

!!subroutine dynpkg  ( u,v,t,q,pes,ghs,ws,wpa, su,sv,st,sq, dtdy,itime,nseq, dsghl )
  subroutine dynpkg  (  dtdy,nseq, dsghl )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Driving routines for dynamics and transport.
! 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use prognostics
   use qadv       
   use commap
   use stdatm  
   use comfm1
   use comhd

   use fspan      !!(wh 2003.11.04)


!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comhyb.h>
!----------------------------------------------------------------------
#include <comctl.h>


!------------------------------Arguments--------------------------------
!
! Input arguments
!

   real(r8),intent(in   )  :: dtdy               ! timestep size ( dyn core )
   integer, intent(in   )  :: nseq
   real(r8),intent(in   )  :: dsghl (plev)


!---------------------------Local workspace-----------------------------
   real(r8),allocatable :: tkk (:,:,:)      !  for fm2003
   real(r8),allocatable :: ukk (:,:,:)      !  for fm2003
   real(r8),allocatable :: vkk (:,:,:)      !  for fm2003

   integer  :: i,j,k
   integer  :: begj


!----------------------------------------------------------
!....  PERFORM THE DYNAMIC INTEGRATON CYCLE
!----------------------------------------------------------


      allocate (tkk (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
      allocate (ukk (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
      allocate (vkk (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
!
      begj = beglatexdyn

     call t_startf('dyfram')

!!      CALL DYFRAM2(NSEQ,DSGHL,IPq                             &
!!                ,DSNP,DSSP,DTDLN,DTDLT,GC,DTDSG               &
!!                ,PMTOP,GHS,SINU,SINV,WTGU,WTGV,SIG,SIGL,DSIG  &
!!                ,CBB,TBB,dtdy,ITIME,SU,SV,ST,SQ              &
!!                ,PeS,U,V,WS,WPA,T,Q,GHI,PLY,TB)

     

        call dyfram2( nseq,dtdy,itime                                          &
                     ,u,v,t,q,ws,pes,wpa,ghs,ghi,ply,tb                        &
                     ,su,sv,st,sq                                              &
                     ,nonos,iord,isor,ep,ipq,dsnp,dssp,dtdln,dtdlt,gc,dtdsg,dsghl  &
                     ,pmtop,sig,sigl,dsig                                      &
                     ,tbb,hbb,cbb,dcbb,psb,tsb                                 &
                     ,dy,wtgu(begj),wtgv(begj)                                 &
                     ,dx,sinu,sinv,oux,ouy,ovx,ovy,ff,cur                      &
                     ,mm1,mp1,mm2,mp2,mm3,mp3,mdj )


      call t_stopf('dyfram')
!
!----------------------------------------------------------
!....  DO FIRST HALF-STEP HORIZONTAL DIFFUSION
!----------------------------------------------------------

   call t_startf('hdifus')

   if (.not.aqua_planet)  then
   !if ((.not.aqua_planet).and.(.not.adiabatic))  then

!$OMP PARALLEL DO PRIVATE (I, J, K)
      DO K=1,plev
        DO J=jbeg0,jend0
          DO I=beglonex,endlonex
             UKK(I,J,K)=U(I,J,K)
             VKK(I,J,K)=V(I,J,K)
             tkk(I,J,K)=T(I,J,K)
             QK(I,J,K)=Q(I,J,K)
          END DO
        END DO
      END DO

      CALL HDIFUS(U,V,T,Q,FRDT,FRDS,FRDU,FRDV,FRDP,TB,PLY,DXVPN,DXVPS)


!$OMP PARALLEL DO PRIVATE (I, J, K)
      DO K=1,plev
        DO J=jbeg0,jend0
          DO I=beglonex,endlonex
             SU(I,J,K)=(U(I,J,K)-UKK(I,J,K))/DTHDFS
             SV(I,J,K)=(V(I,J,K)-VKK(I,J,K))/DTHDFS
             ST(I,J,K)=(T(I,J,K)-tkk(I,J,K))/DTHDFS
             U(I,J,K)=UKK(I,J,K)
             V(I,J,K)=VKK(I,J,K)
             T(I,J,K)=tkk(I,J,K)
             Q(I,J,K)=QK(I,J,K)
          END DO
        END DO
      END DO

   else
      su(:,:,:) = 0.0
      sv(:,:,:) = 0.0
      st(:,:,:) = 0.0

   endif

      deallocate (ukk)
      deallocate (vkk)
      deallocate (tkk)

      call t_stopf('hdifus')

   return

end subroutine dynpkg

