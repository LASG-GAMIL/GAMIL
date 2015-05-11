      SUBROUTINE ICAOSA(P,HICAOS,TICAOS)
C     **********************************
C     **********************************
C
C     SET THE International Civil Aviation Organization
C             (ICAO) STANDARD ATMOSPHERE PROFILES
C
      IMPLICIT NONE
      REAL*8  HTROP0,HSTRA0,HINVE0,TTROP0,FINVE0
      COMMON /ICAOCM/ HTROP0,HSTRA0,HINVE0,TTROP0,FINVE0
      REAL*8  P,HICAOS,TICAOS,HWK
C
      IF( P.GE.226.32E0 ) THEN
        HICAOS = 44330.77E0 - HTROP0*(P**0.19026E0)
        TICAOS =   288.15E0 - TTROP0*HICAOS
      ELSE IF( P.GE.54.7487E0 ) THEN
        HICAOS =     HSTRA0 - 6341.55E0*LOG(P)
        TICAOS =   216.65E0
      ELSE IF( P.GE.8.68014E0 ) THEN
        HICAOS = HINVE0/(P**0.029271E0) - 1.9665E5
        HWK    =     FINVE0 + 4.615736E-6*HICAOS
        TICAOS =   216.65E0 * HWK
      ELSE
        print '(A22)','No SA data available !'
      ENDIF
      RETURN
C     ++++++++++++
C     ::::::::::::
      ENTRY SAICAO
C     ::::::::::::
C     ++++++++++++
      HTROP0   = 44330.77E0 / (1013.25E0**0.19026E0)
      HSTRA0   = 6341.55E0*LOG(226.32E0) + 11000.E0
      HINVE0   = 2.1665E+5 * (54.7487E0**0.029271E0)
      TTROP0   = 288.150E0 * 2.25577E-5
      FINVE0   = 1.00E0    - (4.615736E-6*20000.0E0)
      RETURN
      END
