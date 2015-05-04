#include <misc.h>
#include <params.h>

subroutine iostop (iostat  ,nunit   ,nrec    ,clabel)
!----------------------------------------------------------------------- 
! 
! Purpose: Explain the CRAY FORTRAN I/O error, then call endrun
! 
! Method: Print input diagnostic message.  If Cray, also call "explain function to
!         provide further diagnosis
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: iostat           ! I/O error number from err=iostat option
   integer, intent(in) :: nrec             ! Number of current record (ignored if <=0)
   integer, intent(in) :: nunit            ! I/O Unit number
   character (len=*), intent(in) :: clabel ! Users written diagnostic
!
!---------------------------Local variables-----------------------------
!
   integer i              ! index
#if ( defined CRAY )
   integer iret           ! Return code for the ishell call
   character*16 iolabel   ! String to store error explanation
!
!------------------------------Externals--------------------------------
!
   integer ishell
   external ishell,endrun
#endif
!
!-----------------------------------------------------------------------
!
   if (iostat /= 0) then
      write (6,*) 'IOSTOP:',('*',i=1,30),'  I/O ERROR  ',('*',i=1,29)
      write (6,*) '       ',clabel 
      if (nrec.ge.1) then
         write (6,*) 'I/O Unit = ',nunit,'   Record number = ',nrec,'  Error number = ',iostat
      else
         write (6,*) 'I/O Unit = ',nunit,'   Error number = ',iostat
      end if
      if (iostat.gt.0) then
#if ( defined CRAY )
         write(iolabel(1:16),'(a12,i4)') 'explain lib ',iostat
         write (6,*) iolabel(1:16)
         iret = ishell(iolabel(1:16))
#endif  
      else
         write (6,*) 'End Of File (EOF) was encountered.'
      end if
      call endrun
   end if
!
   return
end subroutine iostop

