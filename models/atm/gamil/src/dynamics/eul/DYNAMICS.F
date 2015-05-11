#include <misc.h>
#include <params.h>

!! (2003.07)
!! (2003.10.23-24)
!! (2003.11.03-04)
!! (2003.10.09)
!! (wb 2004.02.15)
!------------------------

!!      SUBROUTINE DYNAMICS(US,VS,WS,W0,PS,TS,H,DP,SU,SV,ST,FAC,FBC,
!!     _                   PMTOP,HS,SIGL,DSIG,DTDY,ITIME,IHDIFUS,contn)


!!      subroutine dynamics (us,vs,ws,w0,ps,ts,h,dp,su,sv,st,fac,fbc,  &
!!                             hs,          dtdy,itime,ihdifus,contn)

        subroutine dynamics( dtdy,itime,ihdifus,
     _                       us,vs,ws,w0,ps,ts,h,hs,dp,su,sv,st,fac,fbc,
     _                       pmtop,sigl,dsig,
     _                       tbb,hbb,cbb,dcbb,psb,tsb,
     _                       dy,wtgu,wtgv,
     _                       dx,sinu,sinv,oux,ouy,ovx,ovy,ff,cur,
     _                       mm1,mp1,mm2,mp2,mm3,mp3,mdj )

!---------------------------------------------------------------------------------
!
!	    This subroutine is the main part of the model dynamical
!	framework developed by Dr. Bin Wang from LASG/IAP in September
!	of 2001 recoded in May of 2002.
!         In this framework, some new numerical methods are used,
!     including:
!         1) the explicit difference scheme with exact linear and square
!            conservations, developed by Bin Wang, Zhongzhen Ji and
!            Qingcun Zeng, used to solve the atmopsheric equations. By
!            using this scheme, the model can conserve the total mass
!            exactly all the time and keep the conservation of the total
!            available energy when ignoring the outer forcing and friction,
!            and DLT1=DLT2=0.0
!         2) the weighted even-area coordinate along the latitude,
!            developed by Bin Wang, used to reduce the instability
!            of the model at the poles. The coordinate is produced by
!            the subroutine LATMESH
!         3) the flexible leaping-grid method, developed by Bin Wang to
!            further reduce the instability at the ploes. By using the
!            weighted even-area coordinate and the flexible leap-grid method,
!            no filter or smoother is needed in the model.
!         4) the reduction of the standard atmosphere, developed by
!            Qingcun Zeng, to improve the prediction accuracy.
!         This dynamica framework is normalized by Dr. Bin Wang and Rucong Yu,
!     Parallelized by Xin Zhang, translated to F90 by Pu Ye.
!
!---------------------------------------------------------------------------------
        implicit none

#include <PARADYN>

#if (defined SPMD)
#include <commpi.h>
#endif

#include <comfm2.h>
!
!     1) INPUT CONSTANTS:	PMTOP,HS,SIGL,DSIG,DTDY,ITIME,IHDIFUS
!
        real*8   PMTOP          !  PMTOP=10hPa, the pressure of the model top
        real*8   SIGL(NL )      !  half-move vertical layers (sigma(k+1/2)) (unit: 1)
        real*8  DSIG(NL )      !  vertical stepsizes (unit: 1)

        real*8   HS(NX,NY)      !  geopotential of the ground elevation  (unit: m^2/s^2)
        real*8   DTDY           !  time stepsize      (unit: s)

        INTEGER*4  ITIME    !  integrated time
        INTEGER   IHDIFUS   !  control parameter related to the horizontal


! 2) INPUT ARRAYS:SU,SV,ST

        real*8  SU(NX,NY,NL)  !  the forcing on u
        real*8  SV(NX,NY,NL)  !  the forcing on v
        real*8  ST(NX,NY,NL)  !  the forcing on T


! 3) BOTH INPUT AND OUTPUT ARRAYS: US, VS, PS, TS

        real*8  US (NX,NY,NL )  !  zonal wind        (unit: m/s)
        real*8  VS (NX,NY,NL )  !  meridional wind   (unit: m/s)
        real*8  TS (NX,NY,NL)   !  Temperature       (unit: K  )
        real*8  PS (NX,NY   )   !  Surface pressure  (unit: hPa)

! 4) OUTPUT ARRAYS ONLY: WS, W0, H, DPS

        real*8  WS (NX,NY,NZ)  !  diagnosed vertical velocity at the present step (unit: 1/s)
        real*8  W0 (NX,NY,NL)  !  diagnosed vertical velocity at the previous step (unit: 1/s)
        real*8  H  (NX,NY,NZ)  !  diagnosed geopotential height (unit: m^2/s^2)
        real*8  DP (NX,NY   )  !  the tendency of the surface pressure (unit: hPa/s)


! 5) WORKING VARIABLES AND ARRAYS:

!   (5-1) FOR IAP TRANSFORMATION:                      ! in module comfm2

!   (5-2) FOR STANDARD ATMOSPHERE AT THE P-LAYERS
        real*8  TBB (NA)       !  temperature         (unit: K)
        real*8  HBB (NA)       !  geopotential height (unit: m^2/s^2)
        real*8  CBB (NA)       !
        real*8  DCBB(NA)       !

!   (5-3) FOR STANDARD ATMOSPHERE AT THE SIGMA-LAYERS  ! in module comfm2

!   (5-4) SURFACE VARIBALES
        real*8  PSB(NX,NY  )   !  Surface pressure of the standard atmopshere  (unit: hPa)
        real*8  TSB(NX,NY  )   !  Surface temperature of the standard atmopshere (unit: K  )

!   (5-5) FOR FELXIBLE LEAPING-POINT ZONAL DIFFERENCE  ! in module commap
        integer  MM1(NX,NY)     !    MM1=i-MDJ      , if MM1<2, MM1=NX-2+MM1
        integer  MP1(NX,NY)     !    MP1=i+MDJ      , if MP1>NX-1, MP1=MP1-NX+2
        integer  MM2(NX,NY)     !    MM2=i-(MDJ-1)/2, if MM2<2, MM2=NX-2+MM2
        integer  MP2(NX,NY)     !    MP2=i+(MDJ+1)/2, if MP2>NX-1, MP2=MP2-NX+2
        integer  MM3(NX,NY)     !    MM3=i-(MDJ+1)/2, if MM3<2, MM3=NX-2+MM3
        integer  MP3(NX,NY)     !    MP3=i+(MDJ-1)/2, if MP3>NX-1, MP3=MP3-NX+2
        integer  MDJ(NY)        !    leaping span of the difference

!   (5-6) FOR EVEN-AREA MERIDIONAL PARTITION
        real*8  WTGU(NY)       !    area weighting at the normal grid    (unit: 1)
        real*8  WTGV(NY)       !    area weighting at the half-move grid (unit: 1)
        real*8  DY             !    meridional stepsize in computing mesh

!   (5-7) TENDENCIES                                   ! in module comfm2

!   (5-8) OTHERS
        real*8  DX             !  input constant (unit: s)
        real*8  SINU(NY)       !  input constant, sin(theta) at INTEGER grid j
        real*8  SINV(NY)       !  input constant, sin(theta) at half grid j+1/2
        real*8  OUX (NY)       !  input constant, OUX=1/(RAD*SINU*DX*MDJ)
        real*8 OVX (NY)       !  input constant, OUX=1/(RAD*SINV*DX*MDJ)
        real*8  OUY (NY)       !  input constant, OUY=1/(RAD*SINU*DY*WTGU)
        real*8  OVY (NY)       !  input constant, OUY=1/(RAD*SINV*DY*WTGU)
        real*8  FF  (NY)       !  input constant
        real*8  CUR (NY)       !  input constant

        real*8  WK0,WK5,WK6,BYY1,BYY2,BYY3,DT2,TE,TM,TMS
        real*8  UUK(NX,NY),HHK(NX,NY),DUS(NX,NY),DPS2(NX,NY)
        real*8  INNER
        real*8  FAC(nx,ny,nz),FBC(nx,ny,nl)
        integer  I,J,K,KPP,KWB
        integer  begj,endj

        EXTERNAL INNER

!--------------------------------------------------------------------------
! AS FOLLOWS, THE IMPROVED LEAP-FROG AND THE REGENERATION OF
! VARIABLES WILL BE FINISHED IN THE SAME CYCLES.
!--------------------------------------------------------------------------

#if (defined SPMD)
        begj = 2
        endj = ny-1
#else
        begj = 1
        endj = ny
#endif

        DT2=0.5D0*DTDY
!
	IF (ITIME.EQ.0) THEN
!
!	FOR THE NORMAL RUN, DLT1 AND DTL2 MUST BE 1.0D0
!
 	DLT1=1.0d0
 	DLT2=1.0d0
!
!       DLT1=0.OD0 AND DLT2=0.0D0 ONLY FOR TESTING
!       THE CONSERVATION OF AVAILABLE ENERGY
!
!	DLT1=0.0d0
!	DLT2=0.0d0
!
        DO J=begj,endj
	   PS(1,J)=PS(NX-1,J)
	   PS(NX,J)=PS(2,J)
	   TMS=0.0
           DO I=1,nx
	     P (I,J)=SQRT(PS(I,J))
	     C0(I,J)=SQRT(RD*TSB(I,J)/PSB(I,J))
             DO K=1,nl
               PLY(I,J,K) = PS(I,J)*SIGL(K) + PMTOP
               WK5=PLY(I,J,K)/DPALIB
               KPP=INT(WK5)
               WK6=WK5-DFLOAT(KPP)
               CB(I,J,K)=(1.D0-WK6)*CBB(KPP)+WK6*CBB(KPP+1)
               DCB(I,J,K)=(1.D0-WK6)*DCBB(KPP)+WK6*DCBB(KPP+1)
               TB(I,J,K)=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
               CB0(I,J,K)=CB(I,J,K)*P(I,J)/PLY(I,J,K)
               WK5=CB0(I,J,K)*P(I,J)*DSIG(K)
               IF (TMS.LT.WK5) TMS=WK5
             ENDDO
           ENDDO
           KPP=DTDY*TMS*OUX(J)*2.0+0.0001
           NIGW(J)=KPP+1
!
         ENDDO

!
        CALL IAPTRSF(US,VS,WS,PS,TS,UU,VV,P,TT,TB,CB,1)
!
        CALL DIFPS0(UU,P,DPS0,DSIG,OUX,MP1,MP2,MM1,MM2)
        CALL DIFPS1(UU,VV,P,PS,WS,DPS1,DPS0
     _             ,DSIG,DY,OUX,OUY,SINV,MP1,MP2,MM1,MM2,WTGV)
!
        CALL DIAG(UU,VV,P,PS,PLY,TT,US,VS,TS,H
     _           ,HPS,PMTOP,PSB,TSB,TB,CB,FAC,DSIG)
!
	DO K=1,nl
	DO J=begj,endj
	DO I=1,nx
	   W0(I,J,K)=WS(I,J,K)
	END DO
	END DO
	END DO
!
	END IF


!
        DO J=begj,endj
        DO I=1,NX
           PSK(I,J)=PS(I,J)
           DPS(I,J)=DPS0(I,J)+DPS1(I,J)
        ENDDO
        ENDDO
!
        DO K=1,NL
        DO J=begj,endj
        DO I=1,NX
           UK (I,J,K)=UU(I,J,K)
           VK (I,J,K)=VV(I,J,K)
           TTK(I,J,K)=TT(I,J,K)
        ENDDO
        ENDDO
        ENDDO
!
#if (defined SPMD)
        DO J=jbeg,jend
#else
        DO J=2,NY-1
#endif
           DO I=2,NX-1
              HHK(I,J)=C0(I,J)*(PS(I,J)+PMTOP-PSB(I,J))
              UUK(I,J)=0.0d0
              DO K=1,NL
                 UUK(I,J)=UUK(I,J)+UU(I,J,K)*DSIG(K)
              ENDDO
           ENDDO
        ENDDO
!
        DO KWB=1,2
!
!-------------------------------------------------------------------------------
! CALCULATING DU/DT, DV/DT AND DTT/DT.
!  FORMING THE FACTORS OUT OF ALL THE CYCLES AND CLEANING SOME SPECIAL ARRAIES.
!-------------------------------------------------------------------------------
        CALL DIFUVT(UU,US,VV,VS,WS,P,PS,PLY,DPS,TT,H,HPS,CB,DCB
     _           ,SIGL,DSIG,DY,OUX,OVX,OUY,OVY,SINV,FF
     _           ,CUR,DLT1,DLT2,MP1,MP2,MP3,MM1,MM2,MM3,WTGV
     _           ,DU,DV,DTT,SU,SV,ST,FBC)
!
        CALL SEMIU(UUK,HHK,P,DU,DPS1,DUS,DPS0,DTDY,OUX,DSIG,C0)

!
        DO J=begj,endj
        DO I=1,NX
           DPS2(I,J)=DPS1(I,J)+DPS0(I,J)
           PS (I,J)=PSK (I,J)+DT2*DPS2(I,J)
        ENDDO
        ENDDO
!
        !write(*,*) "dynamics_sxj_2"
        DO K=1,nl
        DO J=begj,endj
        DO I=1,nx
           DU(I,J,K)=DU (I,J,K)+    DUS(I,J  )
           TT(I,J,K)=TTK(I,J,K)+DT2*DTT(I,J,K)
           UU(I,J,K)=UK (I,J,K)+DT2*DU (I,J,K)
           VV(I,J,K)=VK (I,J,K)+DT2*DV (I,J,K)
!
           PLY(I,J,K) = PS(I,J)*SIGL(K) + PMTOP
           WK5=PLY(I,J,K)/DPALIB
           KPP=INT(WK5)
           WK6=WK5-DFLOAT(KPP)
           CB(I,J,K)=(1.D0-WK6)*CBB(KPP)+WK6*CBB(KPP+1)
           DCB(I,J,K)=(1.D0-WK6)*DCBB(KPP)+WK6*DCBB(KPP+1)
	     CB0(I,J,K)=CB(I,J,K)*P(I,J)/PLY(I,J,K)
!          DU(I,J,K)=DUS(I,J)
!          DV(I,J,K)=0.0
!          DTT(I,J,K)=0.0
        ENDDO
        ENDDO
        ENDDO
!
!	BYY2=INNER(DU,DV,DTT,DPS0,0,UU,VV,TT,PS,1
!    _             ,DSIG,TSB,PSB,PMTOP,SINU,SINV,WTGU,WTGV)
! 	print *,BYY2
! 	STOP
       ! write(*,*) "dynamics_sxj_3"
!
        IF (KWB.LT.2) THEN
!
        DO J=begj,endj
        DO I=1,NX
           P  (I,J)=SQRT(PS(I,J))
        ENDDO
        ENDDO
!
        ENDIF
!
        CALL DIFPS1(UU,VV,P,PS,WS,DPS1,DPS0
     _             ,DSIG,DY,OUX,OUY,SINV,MP1,MP2,MM1,MM2,WTGV)
!
        CALL DIAG(UU,VV,P,PS,PLY,TT,US,VS,TS,H
     _           ,HPS,PMTOP,PSB,TSB,TB,CB,FAC,DSIG)
!
        DO J=begj,endj
        DO I=1,NX
           DPS (I,J)=DPS0(I,J)+DPS1(I,J)
        ENDDO
        ENDDO
!
        ENDDO

!     CALCULATING DU/DT, DV/DT AND DTT/DT.
!
!     FORMING THE FACTORS OUT OF ALL THE CYCLES AND CLEANING SOME
!     SPECIAL ARRAIES.
!
!
      CALL DIFUVT(UU,US,VV,VS,WS,P,PS,PLY,DPS,TT,H,HPS,CB,DCB
     _           ,SIGL,DSIG,DY,OUX,OVX,OUY,OVY,SINV,FF
     _           ,CUR,DLT1,DLT2,MP1,MP2,MP3,MM1,MM2,MM3,WTGV
     _           ,DU1,DV1,DTT1,SU,SV,ST,FBC)
!
        DO K=1,nl
#if (defined SPMD)
        DO J=jbeg,jend
#else
        DO J=2,ny-1
#endif
          DO I=1,nx
             DU1(I,J,K)=DU1(I,J,K)+DUS(I,J)
          ENDDO
        ENDDO
        END DO
!
        DO J=begj,endj
        DO I=1,nx
           DPS1(I,J)=DPS (I,J)
           DPS (I,J)=DPS2(I,J)
        ENDDO
        ENDDO
!
!     To deduct the inner gravity waves from the tendences of the zonal wind
!     and the temperature during the long-time-stepsize integrations
!
       CALL MINUS_INGW(UU,P,TT,CB0,DSIG,OUX,DU1,DTT1,NIGW)
!
        BYY1=INNER(DU1,DV1,DTT1,DPS1,0,DU ,DV ,DTT ,DPS ,0
     _            ,DSIG,TSB,PSB,PMTOP,SINU,SINV,WTGU,WTGV)
!
        DO J=begj,endj
        DO I=1,nx
           PS(I,J)=PSK(I,J)+DT2*DPS1(I,J)
           P(I,J)=SQRT(PS(I,J))
        ENDDO
        ENDDO
!
!       To add the inner gravity waves to the tendences of the zonal wind
!       and the temperature by the method of short-time-stepsize integrations
!
        CALL PLUS_INGW(UK,P,TTK,CB0,DSIG,OUX,DTDY,SINU,WTGU,DU1,DTT1,BYY1,NIGW)
!
	BYY3=INNER(DU1,DV1,DTT1,DPS1,0,DU1,DV1,DTT1,DPS1,0
     _            ,DSIG,TSB,PSB,PMTOP,SINU,SINV,WTGU,WTGV)

        DT2=DTDY*BYY1/BYY3

        DO J=begj,endj
        DO I=1,nx
           PS(I,J)=PSK(I,J)+DT2*DPS1(I,J)
           P(I,J)=SQRT(PS(I,J))
        ENDDO
!
	TMS=0.0
        DO K=1,nl
        DO I=1,nx
           PLY(I,J,K) = PS(I,J)*SIGL(K) + PMTOP
           WK5=PLY(I,J,K)/DPALIB
           KPP=INT(WK5)
           WK6=WK5-DFLOAT(KPP)
           TB(I,J,K)=(1.D0-WK6)*TBB(KPP)+WK6*TBB(KPP+1)
           CB(I,J,K)=(1.D0-WK6)*CBB(KPP)+WK6*CBB(KPP+1)
           DCB(I,J,K)=(1.D0-WK6)*DCBB(KPP)+WK6*DCBB(KPP+1)
           CB0(I,J,K)=CB(I,J,K)*P(I,J)/PLY(I,J,K)
           WK5=CB0(I,J,K)*P(I,J)*DSIG(K)
           IF (TMS.LT.WK5) TMS=WK5
!
           TT(I,J,K)=TTK(I,J,K)+DT2*DTT1(I,J,K)
           UU(I,J,K)=UK (I,J,K)+DT2*DU1 (I,J,K)
           VV(I,J,K)=VK (I,J,K)+DT2*DV1 (I,J,K)
        ENDDO
        ENDDO
!
        KPP=DTDY*TMS*OUX(J)*2.0+0.0001
        NIGW(J)=KPP+1
!
        ENDDO
 	I=ITIME
        k=i/86400
!!**  	IF (K*86400.EQ.I) THEN
	   IF (DLT1*DLT2.EQ.0) THEN
 	      BYY2=INNER(UK,VK,TTK,PSK,1,UK,VK,TTK,PSK,1
     _                ,DSIG,TSB,PSB,PMTOP,SINU,SINV,WTGU,WTGV)
	      BYY1=INNER(UU,VV,TT,PS,1,UU,VV,TT,PS,1
     _                ,DSIG,TSB,PSB,PMTOP,SINU,SINV,WTGU,WTGV)
#if (defined SPMD)
          if (myrank.eq.0) then
#endif
            WRITE(6,*) float(ITIME)+DTDY,BYY1,DT2
!           WRITE(6,*) BYY2,BYY1,DT2
!          WRITE(112,*) ITIME,BYY1,DT2
#if (defined SPMD)
          endif
#endif
	   ELSE
            CALL TEM(UU,VV,PS,P,TT,HS,CB,TB,DX,DY,DSIG
     _              ,SINU,SINV,WTGU,WTGV,TE,TM)

#if (defined SPMD)
          if (myrank.eq.0) then
#endif
            WRITE(6,'(1x,3e25.18)') TE,TM,DT2  !---sxj--
           ! WRITE(112,'(1x,3e25.18)') TE,TM,DT2
#if (defined SPMD)
          endif
#endif
	   END IF
!!** 	END IF
!
        CALL DIFPS0(UU,P,DPS0,DSIG,OUX,MP1,MP2,MM1,MM2)
        CALL DIFPS1(UU,VV,P,PS,WS,DPS1,DPS0
     _                ,DSIG,DY,OUX,OUY,SINV,MP1,MP2,MM1,MM2,WTGV)
!
        CALL DIAG(UU,VV,P,PS,PLY,TT,US,VS,TS,H
     _           ,HPS,PMTOP,PSB,TSB,TB,CB,FAC,DSIG)
!
        DO J=begj,endj
        DO I=1,nx
 	   DP(I,J)=DPS0(I,J)+DPS1(I,J)
 	END DO
 	END DO
      RETURN
      END
