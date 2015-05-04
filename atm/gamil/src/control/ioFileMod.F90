#include <misc.h>
module ioFileMod
!---------------------------------------------------------------------
!
! Purpose:
!
!	Input/Output file manipulations. Mind file on Mass Store, or local
!	disk etc.
!
! Author: Mariana Vertenstein
!
!---------------------------------------------------------------------
 
   use shr_kind_mod, only: r8 => shr_kind_r8
   implicit none

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

   private

   public getfil      ! Get file from archive
   public putfil      ! Put file to archive
   public opnfil      ! Open file

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

   private shell_cmd  ! Use system command
 
!=======================================================================
   contains
!=======================================================================
 
   subroutine getfil (fulpath, locfn, iflag)
 
! --------------------------------------------------------------------
! obtain local copy of file
! o first check current working directory
! o next check full pathname[fulpath] on disk
! o finally check full pathname[fulpath] on mass store
! --------------------------------------------------------------------
 
! ------------------------ arguments -----------------------------------
   character(len=*), intent(in)  :: fulpath !MSS or permanent disk full pathname
   character(len=*), intent(out) :: locfn   !output local file name
   integer, optional, intent(in) :: iflag   !0=>abort if file not found 1=>do not abort
! --------------------------------------------------------------------
 
! ------------------------ local variables ---------------------------
   integer i               !loop index
   integer klen            !length of fulpath character string
   integer ierr            !error status
   logical lexist          !true if local file exists
   character(len=512) text !mswrite command
! --------------------------------------------------------------------
 
 
! get local file name from full name: start at end. look for first "/"
 
   klen = len_trim(fulpath)
   do i = klen, 1, -1
      if (fulpath(i:i).eq.'/') go to 100
   end do
   i = 0
  100 locfn = fulpath(i+1:klen)
   if (len_trim(locfn) == 0) then
      write(6,*)'(GETFIL): local filename has zero length'
      call endrun
   else
      write(6,*)'(GETFIL): attempting to find local file ',          &
     &     trim(locfn)
   endif
 
! first check if file is in current working directory.
 
   inquire (file=locfn,exist=lexist)
   if (lexist) then
      write (6,*) '(GETFIL): using ',trim(locfn),                    &
     &     ' in current working directory'
      RETURN
   endif
 
! second check for full pathname on disk
 
   inquire(file=fulpath,exist=lexist)
   if (lexist) then
      locfn = trim(fulpath)
      write(6,*)'(GETFIL): using ',trim(fulpath)
      return
   endif
 
! finally check on mass store
 
   text='msread '//trim(locfn)//' '//trim(fulpath)
   call shell_cmd(text, ierr)
   if (ierr==0) then
      write(6,*)'(GETFIL): File ',trim(locfn),' read from MSS'
   else  ! all tries to get file have been unsuccessful
      write(6,*)'(GETFIL): failed cmd=',trim(text)
      if (present(iflag) .and. iflag==0) then
         call endrun
      else
         RETURN
      endif
   end if
 
   return
   end subroutine getfil
 
!=======================================================================
 
   subroutine putfil(locfn   ,mssfpn  ,pass    , &
     &               irt     ,lremov  )
 
!-----------------------------------------------------------------------
! Dispose model output file to Mass Store
!-----------------------------------------------------------------------
 
!------------------------------Arguments--------------------------------
   integer, intent(in) :: irt              ! Mass Store retention time
   character(len=*), intent(in) :: locfn   ! Local filename
   character(len=*), intent(in) :: mssfpn  ! Mass Store full pathname
   character(len=*), intent(in) :: pass    ! write password
   logical, intent(in) :: lremov           ! true=>remove local file
!-----------------------------------------------------------------------
 
!---------------------------Local workspace-----------------------------
   character(len=512) cmd     ! Command string
   character(len=512) cmdtem  ! Temporary for command string
   character(len=  4) crt     ! Retention time as characters
   character(len= 16) wpass   ! Write password
   integer ier                ! error number
!-----------------------------------------------------------------------
 
! Dispose to Mass Store only if nonzero retention period.
 
   if (irt==0) return
 
 
! Non-NCAR users can change the "cmd" below from mswrite to access the
! appropriate archival command for their system.
 
   wpass = ' '
   if (pass(1:1) /= ' ') wpass = ' -w ' // trim(pass)
   write (crt,'(i4)') irt
!! (wh 03.04)
!   write(cmd,'(100a)') 'mswrite ',                                &
!  &     ' -t ',crt,trim(wpass),' ',trim(locfn),' ',trim(mssfpn)
   write(cmd,'(100a)') 'cp ',trim(locfn),' ',trim(mssfpn)
!!
 
! Put mswrite command in background for asynchronous behavior.
 
   if (lremov) then
      cmdtem = '('//trim(cmd)//'; /bin/rm '//trim(locfn)//' )&'
   else
      cmdtem = '('//trim(cmd)//' )&'
   end if
   write(6,*)'PUTFIL: Issuing shell cmd:',trim(cmdtem)
   call shell_cmd(cmdtem, ier)
   if (ier /= 0) then
      write(6,*)'PUTFIL: Error from shell cmd'
!!      call endrun
   end if
 
   return
   end subroutine putfil
 
!=======================================================================
 
   subroutine opnfil (locfn, iun, form, status)
 
!-----------------------------------------------------------------------
! open file locfn in unformatted or formatted form on unit iun
!-----------------------------------------------------------------------
 
! ------------------------ input variables ---------------------------
   character(len=*), intent(in):: locfn  !file name
   integer, intent(in):: iun             !fortran unit number
   character(len=1), intent(in):: form   !file format: u = unformatted. f = formatted
   character(len=*), optional, intent(in):: status !file status
! --------------------------------------------------------------------
 
! ------------------------ local variables ---------------------------
   integer ioe             !error return from fortran open
   character(len=11) ft    !format type: formatted. unformatted
   character(len=11) st    !file status: old or unknown
! --------------------------------------------------------------------
 
   if (len_trim(locfn) == 0) then
      write(6,*)'(OPNFIL): local filename has zero length'
      call endrun
   endif
   if (form=='u' .or. form=='U') then
      ft = 'unformatted'
   else
      ft = 'formatted  '
   end if
   if ( present(status) ) then
      st = status
   else
      st = "unknown"
   end if
   open (unit=iun,file=locfn,status=st, form=ft,iostat=ioe)
   if (ioe /= 0) then
      write(6,*)'(OPNFIL): failed to open file ',trim(locfn),        &
     &     ' on unit ',iun,' ierr=',ioe
      call endrun
   else
      write(6,*)'(OPNFIL): Successfully opened file ',trim(locfn),   &
     &     ' on unit= ',iun
   end if
 
   return
   end subroutine opnfil
 
!=======================================================================
 
   subroutine shell_cmd(text, ier)
 
! ------------------------ arguments -----------------------------------
   character(len=*), intent(in) :: text
   integer         , intent(out):: ier
! ----------------------------------------------------------------------
#if ( defined CRAY )
   integer, external :: ishell ! System routine, execute shell command
#elif (defined IRIX64)
   integer, external :: system_cmd  ! Wrapper to "C" system command
#elif (!defined AIX)
   integer, external :: system ! System routine, execute shell command
#endif
! ----------------------------------------------------------------------
 
#if ( defined CRAY )
   ier = ishell(trim(text))
#elif ( defined AIX )
   call system(trim(text), ier)
#elif ( defined IRIX64 )
   ier = system_cmd(trim(text))
#else
   ier = system(trim(text))
#endif

 
   return
   end subroutine shell_cmd
 
!=======================================================================
 
end module ioFileMod
