#include "params.h"

      SUBROUTINE LATMESH(DY, YTHU, YTHV, WTGU, WTGV)

      USE COMMAP, ONLY: LATMESH_B

      IMPLICIT NONE

      INTEGER, PARAMETER :: NY = PLAT
      REAL(8), PARAMETER :: PI = ATAN(1.0D0) * 4.0D0

      REAL(8), INTENT(OUT) :: DY
      REAL(8), INTENT(OUT) :: YTHU(NY)
      REAL(8), INTENT(OUT) :: YTHV(NY)
      REAL(8), INTENT(OUT) :: WTGU(NY)
      REAL(8), INTENT(OUT) :: WTGV(NY)

      REAL(8) B
      REAL(8) AS
      REAL(8) A  ! A  = B/(0.5*PI) A CONTROL PARAMETER RELATED TO THE AREA-SIZE COMPUTING
      REAL(8) A2 ! A2 = A*2
      REAL(8) DA ! DA = 1/A
      REAL(8) S  ! S  = AS(PI/2)   THE TOTAL AREA-SIZE OF THE WHOLE REGION
      REAL(8) S1 ! S1 = AS(-PI/3)  THE AREA-SIZE AT THE POINT THETA=-PI/3
      REAL(8) S2 ! S2 = AS(PI/3)   THE AREA-SIZE AT THE POINT THETA=PI/3
      REAL(8) U1 ! U1 = 1+2B/3     NEEDED WHEN CALCULATE YTHU, YTHV
      REAL(8) U2 ! U2 = (1-B/3)^2  NEEDED WHEN CALCULATE YTHU, YTHV
      INTEGER I, J
      INTEGER M1 ! M1 = NY-1

      M1 = NY - 1
      B  = LATMESH_B
      A  = B * 2.0D0 / PI
      A2 = A * 2.0D0
      DA = 1.0D0 / A
      S  = PI * (1.0D0 - B / 18.0D0)
      S1 = PI * (1.0D0 - B / 6.0D0) / 6.0D0
      S2 = PI * (5.0D0 - B / 6.0D0) / 6.0D0
      U1 = 1.0D0 + 2.0D0 * B / 3.0D0
      U2 = (1.0D0 - B / 3.0D0) * (1.0D0 - B / 3.0D0)
      DY = S / DFLOAT(M1)
      DO J = 0, M1
        AS = DY * J
        IF (AS <= S1) THEN
          YTHU(J+1) = (SQRT(AS * A2 + U2) - U1) * DA + PI * 0.5D0
          IF (YTHU(J+1) < 0.0) YTHU(J+1) = 0.0D0
          WTGU(J+1) = 1.0D0 - A * (ABS(YTHU(J+1) - PI * 0.5D0) - PI / 3.0D0)
        ELSE IF (AS <= S2) THEN
          YTHU(J+1) = AS - S1 - PI / 3.0D0 + PI * 0.5D0
          WTGU(J+1) = 1.0D0
        ELSE
          YTHU(J+1) = (U1 - SQRT(1.0D0 - MIN(1.0D0, (AS - S2) * A2))) * DA + PI * 0.5D0
          WTGU(J+1) = 1.0D0 - A * (ABS(YTHU(J+1) - PI * 0.5D0) - PI / 3.0D0)
        END IF
        AS = DY * (J + 0.5D0)
        IF (AS <= S1) THEN
          YTHV(J+1) = (SQRT(AS * A2 + U2) - U1) * DA + PI * 0.5D0
          WTGV(J+1) = 1.0D0 - A * (ABS(YTHV(J+1) - PI * 0.5D0) - PI / 3.0D0)
        ELSE IF (AS <= S2) THEN
          YTHV(J+1) = AS - S1 - PI / 3.0D0 + PI * 0.5D0
          WTGV(J+1) = 1.0D0
        ELSE IF (J < M1) THEN
          YTHV(J+1) = (U1 - SQRT(1.0D0 - MIN(1.0D0, (AS - S2) * A2))) * DA + PI * 0.5D0
          WTGV(J+1) = 1.0D0 - A * (ABS(YTHV(J+1) - PI * 0.5D0) - PI / 3.0D0)
        END IF
      END DO

      YTHV(NY) = PI
      WTGV(NY) = 1.0D0 - A * (PI * 0.5D0 - PI / 3.0D0)

      END
