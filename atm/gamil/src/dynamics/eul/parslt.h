!
! $Id: parslt.h,v 1.3 2000/03/02 01:44:09 rosinski Exp $
! $Author: rosinski $
!
!
! Parameters common to many SLT routines
!
      integer ppdy              ! length of interpolation grid stencil
      logical plimdr            ! flag to limit derivatives
!
      parameter(ppdy = 4,  plimdr = .true.)
!
 
