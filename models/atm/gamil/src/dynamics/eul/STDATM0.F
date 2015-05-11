#include <misc.h>
#include <params.h>

      SUBROUTINE STDATM0(TBB,CBB,DCBB,HBB,P00SL,T00SL)
C     *****************
C     *****************
C
C     CREATE LOOK-UP TABLES FOR THE MODEL STANDARD ATMOSPHERE
C                 BY THE USE OF THE SPLINE-FITTING METHOD
C
      IMPLICIT NONE
!
#include <PARADYN>
!
      REAL*8  TBB(NA),CBB(NA),DCBB(NA),HBB(NA),P00SL,T00SL,x,y,z,q
!
      INTEGER KS,KC,NG
C
      REAL*8  PBALIB,PGALIB
      PARAMETER ( KS=13,KC=KS-1 )
      PARAMETER ( PBALIB=0.5E0,PGALIB=1060.0E0,NG=PGALIB/DPALIB )
C
C     SK    == LOG( P )         = PRESSURE [PBALIB,PEALIB] MB
C     HM0   == THE GEOPOTENTIAL HEIGHT AT 1062MB          GPM
C
      REAL*8    SK(KS), TK(KS), CK(KC,3)
C
CWB   Add the data above 10hPa From CIRA climate data
C
C  Pressure (P) SK (Log(P))  TK (T)    CK(N,1)    CK(N,2)    CK(N,3)
C
C   2.3999412E-05  -10.637480   380.7196  -726.4936   313.7478  -34.77775
C   7.5467549E-04  -7.1892230   180.2578   196.7008  -46.02016   2.807630
C   0.7959216E+00  -0.2282546   266.5730  -35.85713   12.61131  -1.790231

CWB
      DATA SK /-0.1063748E+02 ,-0.7189223E+01 ,-0.2282546E+00
CWB
     &        , 0.2302585E+01 , 0.4554327E+01 , 0.5163915E+01
     &        , 0.5664595E+01 , 0.5807347E+01 , 0.6057099E+01
     &        , 0.6265922E+01 , 0.6498247E+01 , 0.6685708E+01
     &        , 0.6967910E+01 /
CWB
      DATA TK / 0.3807196E+03 , 0.1802579E+03 , 0.2665730E+03
CWB
     &        , 0.2275813E+03 , 0.2156167E+03 , 0.2180869E+03
     &        , 0.2288967E+03 , 0.2339356E+03 , 0.2442222E+03
     &        , 0.2540081E+03 , 0.2656011E+03 , 0.2752952E+03
     &        , 0.2900925E+03 /
CWB CK(N,1)
      DATA CK /-726.49360E+00 , 196.70080E+00 ,-35.857130E+00
CWB CK(N,1)
     &        ,-6.4228150E+00 ,-0.8858135E+00 , 11.831290E+00
     &        , 32.385420E+00 , 37.749270E+00 , 44.579530E+00
     &        , 48.629640E+00 , 51.063030E+00 , 52.167430E+00
CWB CK(N,2)
     &        , 313.74780E+00 ,-46.020160E+00 ,  12.61131E+00
CWB CK(N,2)
     &        ,-0.9810501E+00 ,  3.440529E+00 , 17.422030E+00
     &        , 23.633940E+00 , 13.947370E+00 , 13.402850E+00
     &        ,  5.924198E+00 ,  4.503868E+00 ,  1.418971E+00
CWB CK(N,3)
     &        , -34.77775E+00 ,  2.807630E+00 ,-1.7902310E+00
CWB CK(N,3)
     &        ,  0.654469E+00 ,  7.644918E+00 ,  4.133282E+00
     &        ,-22.634740E+00 ,-0.7294737E+00 ,-11.829330E+00
     &        ,-1.9718339E+00 ,-5.5413730E+00 ,-1.6525690E+00/
      REAL*8  ZERO ,HALF ,ONE  ,TWO  ,THREE,FOUR ,SIX
      DATA     ZERO ,HALF ,ONE  ,TWO  ,THREE,FOUR ,SIX  /
     &         0.0E0,0.5E0,1.0E0,2.0E0,3.0E0,4.0E0,6.0E0/
      REAL*8   THIRD,FOURTH,SIXTH,HM0,P,S,DD,D,TW,CW,DW,HW
     &        ,C02,HMW,CKIII,DS,HICAOS,TICAOS,P0J,P0I,DP0,FI0,DP1
      INTEGER K,N,IS,NP,I,II,J0B,I0B,I00,N0
      THIRD      = ONE  / THREE
      FOURTH     = ONE  / FOUR
      SIXTH      = ONE  / SIX
      CALL SAICAO
C
      HM0        = GRAVIT * (-399.0E0)
      P          = PBALIB-DPALIB
      N0         = PBALIB/DPALIB+0.001
      DO 1000 K  = N0,NA
      P          = P    + DPALIB
C
      IF( P.LT.PBALIB .OR. P.GT.PEALIB ) THEN
        PRINT*,'P IN STD.ATM IS OUT OF [10MB,1160MB] ==> STOP !'
        STOP
      ELSE IF( P.LE.PGALIB ) THEN
        S        = LOG( P )
        DD       = S    - SK(1)
        DO 100 N = 1 ,KC
        IS       = N
        D        = DD
        IF( D.EQ.ZERO ) THEN
          TW     = TK(N)
          CW     = CK(N,1)
          DW     = TWO  * CK(N,2)
          HW     = ZERO
          GOTO 105
        ELSE
          NP     = N    + 1
          DD     = S    - SK(NP)
          IF( DD.EQ.ZERO ) THEN
            TW   = TK(NP)
            CW   = CK(NP,1)
            DW   = TWO  * CK(NP,2)
            HW   = ZERO
            IF( N.EQ.KC ) THEN
              IS = KC
            ELSE
              IS = NP
            ENDIF
            GOTO 105
          ELSE IF( DD.LT.ZERO ) THEN
            TW   = ((CK(N,3)*D+CK(N,2))*D+CK(N,1))*D + TK(N)
            CW   = (THREE*CK(N,3)*D + TWO*CK(N,2))*D + CK(N,1)
            DW   =    SIX*CK(N,3)*D + TWO*CK(N,2)
C           INTEGRATE FROM S TO SK(IS=N)
            HW   = ( ( (FOURTH*CK(N,3)*D + THIRD*CK(N,2) )*D
     &           +        HALF*CK(N,1) )*D + TK(N) )*D
            GOTO 105
          ENDIF
        ENDIF
100     CONTINUE
C
105     TBB (K)  = TW
        C02      = RD * (CAPA*TW - CW)
        CBB (K)  = SQRT( C02 )
        DCBB(K)  = RD * (CAPA*CW - DW) / (C02+C02)
C       INTEGRATE FROM SK(IS) TO S(KC)
        HMW      = ZERO
        DO 300 I = IS,KC
        II       = I  + 1
        IF( I.EQ.KC ) THEN
          CKIII  = CK(I ,2) + CK(I,2)
        ELSE
          CKIII  = CK(II,2) + CK(I,2)
        ENDIF
        DS       = SK(II) - SK(I)
        HMW      = HMW + HALF*DS*( TK(II)+TK(I)
     &           -                 CKIII*DS*DS*SIXTH )
300     CONTINUE
        HBB (K)  = RD * ( HMW - HW ) + HM0
      ELSE
        CALL ICAOSA( P,HICAOS,TICAOS )
        HBB(K)   = HICAOS*GRAVIT
        TBB(K)   = TICAOS
        CBB(K)   = CBB(NG)
        DCBB(K)  = 0.0E0
      ENDIF
1000  CONTINUE
      DO 1001 K  = 1 ,N0-1
      HBB(K)     = HBB(N0)
      TBB(K)     = TBB(N0)
      CBB(K)     = CBB(N0)
      DCBB(K)    = DCBB(N0)
1001  CONTINUE
C     CALCULATE PRESSURE P00 & TEMPERATURE T00 AT THE SEA LEVEL
      P0J        = 1010.0E0
      P0I        = 1020.0E0
      J0B        = P0J/DPALIB + 1.E-4
      I0B        = P0I/DPALIB + 1.E-4
      P00SL      =(P0I*HBB(J0B)-P0J*HBB(I0B))/(HBB(J0B)-HBB(I0B))
      DP0        = P00SL  / DPALIB
      I00        = DP0    + 1.E-4
      FI0        = I00
      DP1        = DP0    - FI0
      T00SL      = (ONE - DP1)*TBB(I00) + DP1*TBB(I00+1)
      RETURN
      END
