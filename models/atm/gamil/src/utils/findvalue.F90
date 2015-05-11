function findvalue(ix,n,ain,indxa)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Subroutine for finding ix-th smallest value in the array
! The elements are rearranged so that the ix-th smallest
! element is in the ix place and all smaller elements are
! moved to the elements up to ix (with random order).
!
! Algorithm: Based on the quicksort algorithm.
!
! Author:       T. Craig
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   implicit none
!
! arguments
!
   integer, intent(in) :: ix                ! element to search for
   integer, intent(in) :: n                 ! total number of elements
   integer, intent(inout):: indxa(n)        ! array of integers
   real(r8), intent(in) :: ain(n)           ! array to search
!
   integer findvalue                        ! return value
!
! local variables
!
   integer i,j
   integer il,im,ir

   integer ia
   integer itmp
!
!---------------------------Routine-----------------------------
!
   il=1
   ir=n
   do
      if (ir-il <= 1) then
         if (ir-il == 1) then
            if (ain(indxa(ir)) < ain(indxa(il))) then
               itmp=indxa(il)
               indxa(il)=indxa(ir)
               indxa(ir)=itmp
            endif
         endif
         findvalue=indxa(ix)
         return
      else
         im=(il+ir)/2
         itmp=indxa(im)
         indxa(im)=indxa(il+1)
         indxa(il+1)=itmp
         if (ain(indxa(il+1)) > ain(indxa(ir))) then
            itmp=indxa(il+1)
            indxa(il+1)=indxa(ir)
            indxa(ir)=itmp
         endif
         if (ain(indxa(il)) > ain(indxa(ir))) then
            itmp=indxa(il)
            indxa(il)=indxa(ir)
            indxa(ir)=itmp
         endif
         if (ain(indxa(il+1)) > ain(indxa(il))) then
            itmp=indxa(il+1)
            indxa(il+1)=indxa(il)
            indxa(il)=itmp
         endif
         i=il+1
         j=ir
         ia=indxa(il)
         do
            do
               i=i+1
               if (ain(indxa(i)) >= ain(ia)) exit
            end do
            do
               j=j-1
               if (ain(indxa(j)) <= ain(ia)) exit
            end do
            if (j < i) exit
            itmp=indxa(i)
            indxa(i)=indxa(j)
            indxa(j)=itmp
         end do
         indxa(il)=indxa(j)
         indxa(j)=ia
         if (j >= ix)ir=j-1
         if (j <= ix)il=i
      endif
   end do
end function findvalue




