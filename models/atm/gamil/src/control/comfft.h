!----------------------------------------------------------------------- 
! 
! Purpose: Permanent set-up space for fft
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
      integer pcray          ! length of vector register (words)
      parameter (pcray=64)
!
      common /comfft/ trig(3*plon/2+1,plat)
      common /comfft/ ifax(19,plat)
!
      real(r8) trig          ! trigonometric funct values used by fft
      integer ifax           ! fft factorization of plon/2
!
 
