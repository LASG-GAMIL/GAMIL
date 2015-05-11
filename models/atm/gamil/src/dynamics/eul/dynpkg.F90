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

   real(r8),allocatable :: ply(:,:,:)      !

   real(r8),allocatable :: tb (:,:,:)      !
   real(r8),allocatable :: uk (:,:,:)      !
   real(r8),allocatable :: vk (:,:,:)      !
   real(r8),allocatable :: tk (:,:,:)      !  for fm2003
   real(r8),allocatable :: qk (:,:,:)      !

   integer  :: i,j,k
   integer  :: begj


!----------------------------------------------------------
!....  PERFORM THE DYNAMIC INTEGRATON CYCLE
!----------------------------------------------------------


      allocate (ply (plond, beglatexdyn:endlatexdyn, plevp))

      allocate (tb  (plond, beglatexdyn:endlatexdyn, plev))
      allocate (uk  (plond, beglatexdyn:endlatexdyn, plev))
      allocate (vk  (plond, beglatexdyn:endlatexdyn, plev))
      allocate (tk  (plond, beglatexdyn:endlatexdyn, plev))
      allocate (qk  (plond, beglatexdyn:endlatexdyn, plev))
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

!!   write(6,*) 'calling hdifus.....'

!! if (.not.aqua_planet)  then
   if ((.not.aqua_planet).and.(.not.adiabatic))  then

      DO K=1,plev
        DO J=beglatexdyn+numbnd,endlatexdyn-numbnd
          DO I=1,plond
             UK(I,J,K)=U(I,J,K)
             VK(I,J,K)=V(I,J,K)
             TK(I,J,K)=T(I,J,K)
             QK(I,J,K)=Q(I,J,K)
          END DO
        END DO
      END DO


      CALL HDIFUS(U,V,T,Q,FRDT,FRDS,FRDU,FRDV,FRDP,TB,PLY,DXVPN,DXVPS)

      DO K=1,plev
        DO J=beglatexdyn+numbnd,endlatexdyn-numbnd
          DO I=1,plond
             SU(I,J,K)=(U(I,J,K)-UK(I,J,K))/DTHDFS
             SV(I,J,K)=(V(I,J,K)-VK(I,J,K))/DTHDFS
             ST(I,J,K)=(T(I,J,K)-TK(I,J,K))/DTHDFS
             U(I,J,K)=UK(I,J,K)
             V(I,J,K)=VK(I,J,K)
             T(I,J,K)=TK(I,J,K)
             Q(I,J,K)=QK(I,J,K)
          END DO
        END DO
      END DO

   else
      write(*,*) 'no hdiffus'
      su(:,:,:) = 0.0
      sv(:,:,:) = 0.0
      st(:,:,:) = 0.0

   endif

      call t_stopf('hdifus')

      deallocate (ply)

      deallocate (tb)
      deallocate (uk)
      deallocate (vk)
      deallocate (tk)
      deallocate (qk)


   return

end subroutine dynpkg

