#include <misc.h>
#include <params.h>
subroutine trcab(lchnk   ,ncol    ,                            &
                 k1      ,k2      ,ucfc11  ,ucfc12  ,un2o0   , &
                 un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
                 uco221  ,uco222  ,uco223  ,bn2o0   ,bn2o1   , &
                 bch4    ,to3co2  ,pnm     ,dw      ,pnew    , &
                 s2c     ,uptype  ,dplh2o  ,abplnk1 ,tco2    , &
                 th2o    ,to3     ,abstrc  )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate absorptivity for non nearest layers for CH4, N2O, CFC11 and
! CFC12.
! 
! Method: 
! See CCM3 description for equations.
! 
! Author: J. Kiehl
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid

   implicit none

#include <crdcon.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                    ! chunk identifier
   integer, intent(in) :: ncol                     ! number of atmospheric columns
   integer, intent(in) :: k1,k2                    ! level indices
!
   real(r8), intent(in) :: to3co2(pcols)           ! pressure weighted temperature
   real(r8), intent(in) :: pnm(pcols,pverp)        ! interface pressures
   real(r8), intent(in) :: ucfc11(pcols,pverp)     ! CFC11 path length
   real(r8), intent(in) :: ucfc12(pcols,pverp)     ! CFC12 path length
   real(r8), intent(in) :: un2o0(pcols,pverp)      ! N2O path length
!
   real(r8), intent(in) :: un2o1(pcols,pverp)      ! N2O path length (hot band)
   real(r8), intent(in) :: uch4(pcols,pverp)       ! CH4 path length
   real(r8), intent(in) :: uco211(pcols,pverp)     ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco212(pcols,pverp)     ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco213(pcols,pverp)     ! CO2 9.4 micron band path length
!
   real(r8), intent(in) :: uco221(pcols,pverp)     ! CO2 10.4 micron band path length
   real(r8), intent(in) :: uco222(pcols,pverp)     ! CO2 10.4 micron band path length
   real(r8), intent(in) :: uco223(pcols,pverp)     ! CO2 10.4 micron band path length
   real(r8), intent(in) :: bn2o0(pcols,pverp)      ! pressure factor for n2o
   real(r8), intent(in) :: bn2o1(pcols,pverp)      ! pressure factor for n2o
!
   real(r8), intent(in) :: bch4(pcols,pverp)       ! pressure factor for ch4
   real(r8), intent(in) :: dw(pcols)               ! h2o path length
   real(r8), intent(in) :: pnew(pcols)             ! pressure
   real(r8), intent(in) :: s2c(pcols,pverp)        ! continuum path length
   real(r8), intent(in) :: uptype(pcols,pverp)     ! p-type h2o path length
!
   real(r8), intent(in) :: dplh2o(pcols)           ! p squared h2o path length
   real(r8), intent(in) :: abplnk1(14,pcols,pverp) ! Planck factor
   real(r8), intent(in) :: tco2(pcols)             ! co2 transmission factor
   real(r8), intent(in) :: th2o(pcols)             ! h2o transmission factor
   real(r8), intent(in) :: to3(pcols)              ! o3 transmission factor
!
!  Output Arguments
!
   real(r8), intent(out) :: abstrc(pcols)           ! total trace gas absorptivity
!
!--------------------------Local Variables------------------------------
!
   integer  i,l                     ! loop counters

   real(r8) sqti(pcols)             ! square root of mean temp
   real(r8) du1                     ! cfc11 path length
   real(r8) du2                     ! cfc12 path length
   real(r8) acfc1                   ! cfc11 absorptivity 798 cm-1
   real(r8) acfc2                   ! cfc11 absorptivity 846 cm-1
!
   real(r8) acfc3                   ! cfc11 absorptivity 933 cm-1
   real(r8) acfc4                   ! cfc11 absorptivity 1085 cm-1
   real(r8) acfc5                   ! cfc12 absorptivity 889 cm-1
   real(r8) acfc6                   ! cfc12 absorptivity 923 cm-1
   real(r8) acfc7                   ! cfc12 absorptivity 1102 cm-1
!
   real(r8) acfc8                   ! cfc12 absorptivity 1161 cm-1
   real(r8) du01                    ! n2o path length
   real(r8) dbeta01                 ! n2o pressure factor
   real(r8) dbeta11                 !         "
   real(r8) an2o1                   ! absorptivity of 1285 cm-1 n2o band
!
   real(r8) du02                    ! n2o path length
   real(r8) dbeta02                 ! n2o pressure factor
   real(r8) an2o2                   ! absorptivity of 589 cm-1 n2o band
   real(r8) du03                    ! n2o path length
   real(r8) dbeta03                 ! n2o pressure factor
!
   real(r8) an2o3                   ! absorptivity of 1168 cm-1 n2o band
   real(r8) duch4                   ! ch4 path length
   real(r8) dbetac                  ! ch4 pressure factor
   real(r8) ach4                    ! absorptivity of 1306 cm-1 ch4 band
   real(r8) du11                    ! co2 path length
!
   real(r8) du12                    !       "
   real(r8) du13                    !       "
   real(r8) dbetc1                  ! co2 pressure factor
   real(r8) dbetc2                  ! co2 pressure factor
   real(r8) aco21                   ! absorptivity of 1064 cm-1 band
!
   real(r8) du21                    ! co2 path length
   real(r8) du22                    !       "
   real(r8) du23                    !       "
   real(r8) aco22                   ! absorptivity of 961 cm-1 band
   real(r8) tt(pcols)               ! temp. factor for h2o overlap factor
!
   real(r8) psi1                    !                 "
   real(r8) phi1                    !                 "
   real(r8) p1                      ! h2o overlap factor
   real(r8) w1                      !        "
   real(r8) ds2c(pcols)             ! continuum path length
!
   real(r8) duptyp(pcols)           ! p-type path length
   real(r8) tw(pcols,6)             ! h2o transmission factor
   real(r8) g1(6)                   !         "
   real(r8) g2(6)                   !         "
   real(r8) g3(6)                   !         "
!
   real(r8) g4(6)                   !         "
   real(r8) ab(6)                   ! h2o temp. factor
   real(r8) bb(6)                   !         "
   real(r8) abp(6)                  !         "
   real(r8) bbp(6)                  !         "
!
   real(r8) tcfc3                   ! transmission for cfc11 band
   real(r8) tcfc4                   ! transmission for cfc11 band
   real(r8) tcfc6                   ! transmission for cfc12 band
   real(r8) tcfc7                   ! transmission for cfc12 band
   real(r8) tcfc8                   ! transmission for cfc12 band
!
   real(r8) tlw                     ! h2o transmission
   real(r8) tch4                    ! ch4 transmission
!
!--------------------------Data Statements------------------------------
!
   data g1 /0.0468556,0.0397454,0.0407664,0.0304380,0.0540398,0.0321962/
   data g2 /14.4832,4.30242,5.23523,3.25342,0.698935,16.5599/
   data g3 /26.1898,18.4476,15.3633,12.1927,9.14992,8.07092/
   data g4 /0.0261782,0.0369516,0.0307266,0.0243854,0.0182932,0.0161418/
   data ab /3.0857e-2,2.3524e-2,1.7310e-2,2.6661e-2,2.8074e-2,2.2915e-2/
   data bb /-1.3512e-4,-6.8320e-5,-3.2609e-5,-1.0228e-5,-9.5743e-5,-1.0304e-4/
   data abp/2.9129e-2,2.4101e-2,1.9821e-2,2.6904e-2,2.9458e-2,1.9892e-2/
   data bbp/-1.3139e-4,-5.5688e-5,-4.6380e-5,-8.0362e-5,-1.0115e-4,-8.8061e-5/
!
!--------------------------Statement Functions--------------------------
!
   real(r8) func, u, b
   func(u,b) = u/sqrt(4.0 + u*(1.0 + 1.0 / b))
!
!------------------------------------------------------------------------
!
   do i = 1,ncol
      sqti(i) = sqrt(to3co2(i))
!
! h2o transmission
!
      tt(i) = abs(to3co2(i) - 250.0)
      ds2c(i) = abs(s2c(i,k1) - s2c(i,k2))
      duptyp(i) = abs(uptype(i,k1) - uptype(i,k2))
   end do
!
   do l = 1,6
      do i = 1,ncol
         psi1 = exp(abp(l)*tt(i) + bbp(l)*tt(i)*tt(i))
         phi1 = exp(ab(l)*tt(i) + bb(l)*tt(i)*tt(i))
         p1 = pnew(i)*(psi1/phi1)/sslp
         w1 = dw(i)*phi1
         tw(i,l) = exp(-g1(l)*p1*(sqrt(1.0 + g2(l)*(w1/p1)) - 1.0) - &
                   g3(l)*ds2c(i)-g4(l)*duptyp(i))
      end do
   end do
!
   do i = 1,ncol
      du1 = abs(ucfc11(i,k1) - ucfc11(i,k2))
      du2 = abs(ucfc12(i,k1) - ucfc12(i,k2))
!
! cfc transmissions
!
      tcfc3 = exp(-175.005*du1)
      tcfc4 = exp(-1202.18*du1)
      tcfc6 = exp(-5786.73*du2)
      tcfc7 = exp(-2873.51*du2)
      tcfc8 = exp(-2085.59*du2)
!
! Absorptivity for CFC11 bands
!
      acfc1 =  50.0*(1.0 - exp(-54.09*du1))*tw(i,1)*abplnk1(7,i,k2)
      acfc2 =  60.0*(1.0 - exp(-5130.03*du1))*tw(i,2)*abplnk1(8,i,k2)
      acfc3 =  60.0*(1.0 - tcfc3)*tw(i,4)*tcfc6*abplnk1(9,i,k2)
      acfc4 = 100.0*(1.0 - tcfc4)*tw(i,5)*abplnk1(10,i,k2)
!
! Absorptivity for CFC12 bands
!
      acfc5 = 45.0*(1.0 - exp(-1272.35*du2))*tw(i,3)*abplnk1(11,i,k2)
      acfc6 = 50.0*(1.0 - tcfc6)* tw(i,4) * abplnk1(12,i,k2)
      acfc7 = 80.0*(1.0 - tcfc7)* tw(i,5) * tcfc4*abplnk1(13,i,k2)
      acfc8 = 70.0*(1.0 - tcfc8)* tw(i,6) * abplnk1(14,i,k2)
!
! Emissivity for CH4 band 1306 cm-1
!
      tlw = exp(-1.0*sqrt(dplh2o(i)))
      duch4 = abs(uch4(i,k1) - uch4(i,k2))
      dbetac = abs(bch4(i,k1) - bch4(i,k2))/duch4
      ach4 = 6.00444*sqti(i)*log(1.0 + func(duch4,dbetac))*tlw*abplnk1(3,i,k2)
      tch4 = 1.0/(1.0 + 0.02*func(duch4,dbetac))
!
! Absorptivity for N2O bands
!
      du01 = abs(un2o0(i,k1) - un2o0(i,k2))
      du11 = abs(un2o1(i,k1) - un2o1(i,k2))
      dbeta01 = abs(bn2o0(i,k1) - bn2o0(i,k2))/du01
      dbeta11 = abs(bn2o1(i,k1) - bn2o1(i,k2))/du11
!
! 1285 cm-1 band
!
      an2o1 = 2.35558*sqti(i)*log(1.0 + func(du01,dbeta01) &
              + func(du11,dbeta11))*tlw*tch4*abplnk1(4,i,k2)
      du02 = 0.100090*du01
      du12 = 0.0992746*du11
      dbeta02 = 0.964282*dbeta01
!
! 589 cm-1 band
!
      an2o2 = 2.65581*sqti(i)*log(1.0 + func(du02,dbeta02) + &
              func(du12,dbeta02))*th2o(i)*tco2(i)*abplnk1(5,i,k2)
      du03 = 0.0333767*du01
      dbeta03 = 0.982143*dbeta01
!
! 1168 cm-1 band
!
      an2o3 = 2.54034*sqti(i)*log(1.0 + func(du03,dbeta03))* &
              tw(i,6)*tcfc8*abplnk1(6,i,k2)
!
! Emissivity for 1064 cm-1 band of CO2
!
      du11 = abs(uco211(i,k1) - uco211(i,k2))
      du12 = abs(uco212(i,k1) - uco212(i,k2))
      du13 = abs(uco213(i,k1) - uco213(i,k2))
      dbetc1 = 2.97558*abs(pnm(i,k1) + pnm(i,k2))/(2.0*sslp*sqti(i))
      dbetc2 = 2.0*dbetc1
      aco21 = 3.7571*sqti(i)*log(1.0 + func(du11,dbetc1) &
              + func(du12,dbetc2) + func(du13,dbetc2)) &
              *to3(i)*tw(i,5)*tcfc4*tcfc7*abplnk1(2,i,k2)
!
! Emissivity for 961 cm-1 band
!
      du21 = abs(uco221(i,k1) - uco221(i,k2))
      du22 = abs(uco222(i,k1) - uco222(i,k2))
      du23 = abs(uco223(i,k1) - uco223(i,k2))
      aco22 = 3.8443*sqti(i)*log(1.0 + func(du21,dbetc1) &
              + func(du22,dbetc1) + func(du23,dbetc2)) &
              *tw(i,4)*tcfc3*tcfc6*abplnk1(1,i,k2)
!
! total trace gas absorptivity
!
      abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 + &
                  acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 + &
                  aco21 + aco22
   end do
!
   return
!
end subroutine trcab



