module dynconst
!
! Constants used in dynamics
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_const_mod, only: shr_const_omega, shr_const_rearth
   implicit none

   private

   public dynconsti

   real(r8), public, parameter :: omega = shr_const_omega     ! Angular velocity of Earth's rotation
   real(r8), public, parameter :: rearth = shr_const_rearth   ! Radius of the earth
   real(r8), public, parameter :: ra = 1.0 / rearth           ! Reciprocal of earth radius
   real(r8), public, save :: ez                    ! Coriolis expansion coeff -> omega/sqrt(0.375)

contains

subroutine dynconsti
   ez = omega / sqrt(0.375)
end subroutine dynconsti


end module dynconst
