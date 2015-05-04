subroutine sortarray(n, ain, indxa) 
!-----------------------------------------------
!
! Purpose:
!       Sort an array
! Alogrithm:
!       Based on Shell's sorting method.
!
! Author: T. Craig
!-----------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   implicit none
!
!  Arguments
!
   integer , intent(in) :: n             ! total number of elements
   integer , intent(inout) :: indxa(n)   ! array of integers
   real(r8), intent(inout) :: ain(n)     ! array to sort
!
!  local variables
!
   integer :: i, j                ! Loop indices
   integer :: ni                  ! Starting increment
   integer :: itmp                ! Temporary index
   real(r8):: atmp                ! Temporary value to swap
 
   ni = 1 
   do while(.TRUE.) 
      ni = 3*ni + 1 
      if (ni <= n) cycle  
      exit  
   end do 
 
   do while(.TRUE.) 
      ni = ni/3 
      do i = ni + 1, n 
         atmp = ain(i) 
         itmp = indxa(i) 
         j = i 
         do while(.TRUE.) 
            if (ain(j-ni) <= atmp) exit  
            ain(j) = ain(j-ni) 
            indxa(j) = indxa(j-ni) 
            j = j - ni 
            if (j > ni) cycle  
            exit  
         end do 
         ain(j) = atmp 
         indxa(j) = itmp 
      end do 
      if (ni > 1) cycle  
      exit  
   end do 
   return  
 
end subroutine sortarray
