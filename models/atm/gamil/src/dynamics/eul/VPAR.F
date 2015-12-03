#include <misc.h>
#include <params.h>

      SUBROUTINE VPAR(PMTOP, P0, SIG, SIGL, DSIG)

      IMPLICIT NONE

#include <PARADYN>

      REAL*8 P0          ! MODEL BOTTOM PRESSURE,     INTENT( IN)
      REAL*8 PMTOP       ! MODEL TOP PRESSURE,        INTENT(OUT)
      REAL*8 SIG (NZ)    ! SIGMA AT MODEL INTERFACES, INTENT(OUT)
      REAL*8 SIGL(NL)    ! SIGMA AT MODEL LEVELS,     INTENT(OUT)
      REAL*8  DSIG(NL)   ! VERTICAL STEPSIZE,         INTENT(OUT)
      REAL*8  PS0        ! REFERRENCE SURFACE PRESSRE OF CAM2
      REAL*8  HYAI(NZ)   ! PARAMETER 'A' IN CAM2 PRESSURE FORMULA
      REAL*8  HYBI(NZ)   ! PARAMETER 'B' IN CAM2 PRESSURE FORMULA
      INTEGER K          ! LAYER INDEX

      DATA PS0  /1.0D3/  ! UNIT: HPA

      DATA HYAI /0.00219406700000001D0, 0.00489520900000001D0,
     _           0.009882418D0,         0.01805201D0,
     _           0.02983724D0,          0.0446233400000002D0,
     _           0.0616058700000002D0,  0.0785124300000004D0,
     _           0.0773127100000002D0,  0.0759013100000003D0,
     _           0.0742408600000002D0,  0.0722874400000002D0,
     _           0.0699893299999998D0,  0.06728574D0,
     _           0.06410509D0,          0.0603632200000002D0,
     _           0.0559611100000001D0,  0.0507822500000001D0,
     _           0.0446896000000001D0,  0.0375219099999999D0,
     _           0.0290894900000001D0,  0.02084739D0,
     _           0.01334443D0,          0.00708499000000001D0,
     _           0.00252136D0,          0.0D0,
     _           0.0D0/

      DATA HYBI /0.0D0,                 0.0D0,
     _           0.0D0,                 0.0D0,
     _           0.0D0,                 0.0D0,
     _           0.0D0,                 0.0D0,
     _           0.01505309D0,          0.03276228D0,
     _           0.05359622D0,          0.0781062700000006D0,
     _           0.1069411D0,           0.140863700000001D0,
     _           0.180772D0,            0.227722D0,
     _           0.282956200000001D0,   0.347936400000002D0,
     _           0.4243822D0,           0.514316800000003D0,
     _           0.620120200000002D0,   0.723535500000004D0,
     _           0.817676800000001D0,   0.896215300000001D0,
     _           0.953476100000003D0,   0.9851122D0,
     _           1.0D0/

      PMTOP = HYAI(1)*PS0

      SIG(1) = 0.0D0
      DO K = 1, NL
        SIG (K+1) = (HYAI(K+1)*PS0+HYBI(K+1)*P0-PMTOP)/(P0-PMTOP)
        SIGL(K  ) = 0.5D0*(SIG(K)+SIG(K+1))
        DSIG(K  ) = SIG(K+1)-SIG(K)
      END DO

      RETURN
      END
