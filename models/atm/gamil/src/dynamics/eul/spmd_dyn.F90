#include <misc.h>
#include <params.h>

module spmd_dyn

!! (wanhui 2003.10.18)
!! (zx,wh  2003.10.21)
!! (wanhui 2003.10.21)
!! (wanhui 2003.10.24)
!----------------------------------------------------------------------- 
! 
! Purpose: SPMD implementation of CAM.  Currently used for both dynamics
!          and physics, but ultimately the physics part should be broken off.
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------

#if (defined SPMD)

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plat, masterproc, iam, beglatex, endlatex, numbnd, numlats, numlatsex, &
!!                         beglat, endlat, begirow, endirow, plev
!!                         beglat, endlat,                   plev, plond      !!(2003.10.21)
                           beglat, endlat, beglatexdyn, endlatexdyn, plev, plond !!(2003.10.24)
   use constituents, only: pcnst
   use mpishorthand, only: mpir8, mpicom
   use infnan, only: inf
   implicit none

   private
!! public spmdinit_dyn, compute_gsfactors,                        pair, ceil2
   public spmdinit_dyn, compute_gsfactors, compute_gsfactors_dyn, pair, ceil2
   save
   integer, public :: npes                 ! Total number of MPI tasks
   integer, public :: cut(2,0:plat-1)      ! partition for MPI tasks
   integer, public :: cutex(2,0:plat-1)    ! extended partition 
   integer, public :: proc(plat)           ! MPI task id associated with a given lat.
   integer, public :: neighs               ! number of south neighbors to comm guardcells
   integer, public :: neighn               ! number of north neighbors to comm guardcells
   integer, public :: npessp               ! number of MPI tasks in spectral space
   integer, public :: maxm                 ! max number of Fourier wavenumbers per MPI task
   integer, public :: numm(0:plat-1)       ! number of Fourier wavenumbers owned per task
   integer, public :: bsiz                 ! buffer size
   integer, public :: maxlats              ! max number of lats on any MPI task
!  integer, public, allocatable :: nlat_p(:)    ! number of latitudes per MPI task
   integer, public :: nlat_p(0:1000)    ! number of latitudes per MPI task

   real(r8), public, allocatable :: buf1(:),buf2(:) ! buffers for packing MPI msgs

CONTAINS

!========================================================================

   subroutine spmdinit_dyn ()
!----------------------------------------------------------------------- 
! 
! Purpose: Distribute latitudes among available processors
! 
! Method: Distribution is S->N for processors 0->npes
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!     use pspect, only: maxlats
!-----------------------------------------------------------------------
!
! Local workspace
!
      integer procid    ! processor id
      integer procids   ! processor id SH
      integer procidn   ! processor id NH
      integer lat       ! latitude index
      integer iend      ! ending latitude band of work for a given proc
      integer workleft  ! amount of work still to be parcelled out
      integer actual    ! actual amount of work parcelled out
      integer ideal     ! ideal amt of work to parcel out
      integer pesleft   ! number of procs still to be given work
      integer isum      ! running total of work parcelled out
      integer smostlat  ! southern-most latitude index
      integer nmostlat  ! northern-most latitude index
      integer m2,m3,m5  ! 2, 3, 5 prime factors for problem decomposition
!
!-----------------------------------------------------------------------
!
! Allocate memory for number of lats per proc
!
!     allocate (nlat_p (0:npes-1))
      nlat_p(0:npes-1) = 0
!
! Make sure number of PEs and number of latitudes are kosher
!
      call factor (plat, m2, m3, m5)

      if (m2 < 1) then
         write(6,*) 'FACTOR: Problem size is not divisible by 2'
         call endrun
      end if

      if (masterproc) then
         write (6,*) 'Problem factors: 2**',m2,' * 3**',m3,' * 5**',m5
      end if
      call factor (npes, m2, m3, m5)
      
      if (mod(npes,2) /= 0) then
         write(6,*)'SPMDINIT_DYN: nprocs(',npes,') must be a multiple of 2'
         call endrun
      end if

      workleft = plat/2
      pesleft = npes/2
      iend = 0
   
      maxlats = 0
      do procids=0,npes/2-1
         procidn = npes - procids - 1
         if (workleft > 0) then
            ideal = workleft/pesleft
            cut(1,procids) = iend + 1
            lat = cut(1,procids)
            actual = 1
10          if (lat+1 <= plat/2) then
               if (actual+1 <= ideal .or. pesleft == 1) then
                  lat = lat + 1
                  actual = actual + 1
                  goto 10
               end if
            end if
            cut(2,procids) = lat
!
! Assign mirror latitudes
!
            cut(1,procidn) = plat - cut(2,procids) + 1
            cut(2,procidn) = plat - cut(1,procids) + 1
         else
            write(6,*)'SPMDINIT_DYN: Ran out of work to parcel to processors'
            call endrun
         end if
         
         nlat_p(procids) = actual
         nlat_p(procidn) = actual
         maxlats = max (maxlats, actual)
         
         if (iam == procids .or. iam == procidn) then
            beglat = cut(1,iam)
            endlat = cut(2,iam)
            numlats = actual
!!          begirow = cut(1,procids)
!!          endirow = cut(2,procids)
         end if
!
! Prepare for next iteration
!
         iend = lat
         workleft = workleft - actual
         pesleft = pesleft - 1
      end do

      if (workleft /= 0) then
         write(6,*)'SPMDINIT_DYN: Workleft not zero.  Value is ',workleft
         call endrun
      end if
   
      do procid=0,npes-1
         if (masterproc) then
            write(6,*)'procid ',procid,' assigned ', &
                      cut(2,procid)-cut(1,procid)+1,' latitude values from', &
                      cut(1,procid),' through ',cut(2,procid)
         end if
!
! Determine which processor is responsible for the defined latitudes
!
         do lat=cut(1,procid),cut(2,procid)
            proc(lat) = procid
         end do
!
! The extended regions are simply "numbnd" wider at each
! side. The extended region do not go beyond 1 and plat, though
!
         cutex(1,procid) = cut(1,procid) - numbnd
         cutex(2,procid) = cut(2,procid) + numbnd
         if (iam == procid) then
!!          beglatex = cutex(1,procid) + numbnd
!!          endlatex = cutex(2,procid) + numbnd
            beglatex = cutex(1,procid)         !!(wh 2003.10.24)
            endlatex = cutex(2,procid)         !!
            beglatexdyn = plat + 1 - endlatex  !!
            endlatexdyn = plat + 1 - beglatex  !!
            numlatsex = endlatex - beglatex + 1
         end if
      end do
!
! Number of neighbor processors needed for boundary communication.  North
! first.
!
      isum = 0
      neighn = 0
      
      do procid=iam+1,npes-1
         nmostlat = cut(2,procid)
         isum = isum + cut(2,procid) - cut(1,procid) + 1
         neighn = neighn + 1
         if (isum >= numbnd) goto 20
      end do
      
20    if (iam /= npes-1 .and. isum < numbnd .and. nmostlat /= plat) then
         write (6,*) 'SPMDINIT_DYN: Something wrong in computation of northern neighbors'
         call endrun
      end if
      
      isum = 0
      neighs = 0
      
      do procid=iam-1,0,-1
         smostlat = cut(1,procid)
         isum = isum + cut(2,procid) - cut(1,procid) + 1
         neighs = neighs + 1
         if (isum >= numbnd) goto 30
      end do

30    if (iam /= 0 .and. isum < numbnd .and. smostlat /= 1) then
         write(6,*)'Something wrong in computation of southern neighbors'
         call endrun
      end if

      if (masterproc) then
         write(6,*)'-----------------------------------------'
         write(6,*)'Number of lats passed north & south = ',numbnd
         write(6,*)'Node  Partition  Extended Partition'
         write(6,*)'-----------------------------------------'
         do procid=0,npes-1
            write(6,200) procid,cut(1,procid),cut(2,procid) ,cutex(1,procid), cutex(2,procid)
200         format(i3,4x,i3,'-',i3,7x,i3,'-',i3)
         end do
      end if
!      write(6,*)'iam=',iam,'Number of south neighbors needed for bndry exchange = ',neighs
!      write(6,*)'iam=',iam,'Number of north neighbors needed for bndry exchange = ',neighn

!!    call decomp_wavenumbers ()
!!    call spmdbuf ()

    write(6,11) '*** p',iam,'beglat:endlat=',beglat,':',endlat,'numlats=',numlats
    write(6,11) '*** p',iam,'beglatex:endlatex=',beglatex,':',endlatex,'numlatsex=',numlatsex
    write(6,11) '*** p',iam,'beglatexdyn:endlatexdyn=',beglatexdyn,':',endlatexdyn
    write(6,11) '*** p',iam,'numbnd=',numbnd
11  format(1x,a6,i2,a25,i3,a2,i3,a12,i3)
!!    stop

      return
   end subroutine spmdinit_dyn

!========================================================================

   subroutine factor (nitems, m2, m3, m5)
!----------------------------------------------------------------------- 
! 
! Purpose: Factor a given number into powers of 2,3,5
! 
! Method: Brute force application of "mod" function
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nitems      ! Number to be factored into powers of 2,3,5
      integer, intent(out) :: m2,m3,m5   ! Powers of 2, 3, and 5 respectively
!
! Local workspace
!
      integer num                        ! current number to be factored
!
!-----------------------------------------------------------------------
!
      num = nitems
      m2 = 0
      m3 = 0
      m5 = 0
      
2     if (mod(num,2) == 0) then
         m2 = m2 + 1
         num = num/2
         goto 2
      end if
      
3     if (mod(num,3) == 0) then
         m3 = m3 + 1
         num = num/3
         goto 3
      end if
      
5     if (mod(num,5) == 0) then
         m5 = m5 + 1
         num = num/5
         goto 5
      end if
      
      if (num /= 1) then
         write(6,*) 'FACTOR: ',nitems,' has a prime factor other than 2, 3, or 5.  Aborting...'
         call endrun
      end if
      
      return
   end subroutine factor

!========================================================================

   subroutine decomp_wavenumbers
!----------------------------------------------------------------------- 
! 
! Purpose: partition the spectral work among the given number of processors
! 
! Method: Make the labor division as equal as possible given loop lengths
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      use pspect, only: pmmax
      use comspe, only: nlen, begm, endm, nstart
!
! Local workspace
!
      integer endfourier  ! ending fourier wavenumber         
      integer procid      ! processor id
      integer m           ! fourier wavenumber index
      integer workleft    ! amt of work still to be parceled out
      integer actual      ! actual amt of work given to a particular proc
      integer ideal       ! ideal amt of work given to a particular proc
      integer pesleft     ! number of pes still to be given work
      integer test        ! test value to compare to ideal amt of work

!-----------------------------------------------------------------------

      workleft = nstart(pmmax) + nlen(pmmax)
      pesleft = min(pmmax,npes)
      endfourier = 0
      npessp = 0
      maxm = 0

      do procid = 0,npes-1
         if (workleft > 0) then
            npessp = npessp + 1
            ideal = workleft / pesleft
            begm(procid) = endfourier + 1
            m = begm(procid)
            actual = nlen(m)
            
1           if (m+1 <= pmmax) then
               test = actual + nlen(m+1)
               if (test <= ideal) then
                  m = m + 1
                  actual = test
                  goto 1
               else if (test > ideal) then
                  if (test-ideal < ideal-actual) then
                     m = m + 1
                     actual = test
                  end if
               end if
            end if
            
            endm(procid) = m
            endfourier = m
            workleft = workleft - actual
            pesleft = pesleft - 1

            if (masterproc) then
               write(6,*)'procid ',procid,' assigned ', endm(procid)-begm(procid)+1, &
                         ' m values from ', begm(procid),' through ',endm(procid)
            end if
         else
            begm(procid) = 0
            endm(procid) = -1
         end if
         numm(procid) = endm(procid) - begm(procid) + 1
         if (numm(procid) > maxm) maxm = numm(procid)
      end do

      if (workleft/=0) then
         write(6,*)'MCUTS: Workleft not zero.  Value is ',workleft
         call endrun
      end if
   
      return
   end subroutine decomp_wavenumbers

!========================================================================

   integer function pair(np,p,k)

      integer np,p,k,q
      q = ieor(p,k)
      if(q.gt.np-1) then
         pair = -1
      else
         pair = q
      endif
      return

   end function pair

!========================================================================

  integer function ceil2(n)
     integer n,p
     p=1
     do while(p.lt.n)
        p=p*2
     enddo
     ceil2=p
     return
  end function ceil2
  
!========================================================================

  subroutine spmdbuf
!----------------------------------------------------------------------- 
! 
! Purpose: allocate spmd pack buffers used in pairwise all-all exchanges
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
     use comspe, only: begm, endm, nlen

     integer maxcount(4),m
     integer length,i
!
! realloc4 max: 16 2 plev*numm*numlats (e.g. grt1)
!               2  2     *numm*numlats (grlps1, grlps2)
!               2             *numlats (begirow, endirow)
!
     maxcount(1) = maxlats*(2*maxm*(plev*16 + 2) + 2)

!
! realloc6 max: 3 plev*(nlen*numm)  (e.g. vz)
!               1     *(nlen*numm)  (alps)
!               2                   (length, mstrt)
!
     length = 0
     do i=0,npessp-1
        do m=begm(i),endm(i)
           length = length + 2*nlen(m)
        end do
     end do

     maxcount(2) = length*(1 + 3*plev) + 2      
!
! realloc5 max: 3                 (len,beglat,numlats)
!               1 numlats         (tmass)
!               5 numlats  *pcnst (e.g. hw1lat)
!               2 4*numlats*pcnst (e.g.hw2al)
!
     maxcount(3) = 3 + maxlats*(1 + (5 + 2*4)*pcnst)
!
! realloc7 max: 2                  (beglat, numlats)
!               3 plev *numlats    (e.g. vmax2d)
!               5      *numlats    (e.g. psurf)
!
     maxcount(4) = maxlats*(3*plev + 5) + 2
     m = maxval(maxcount)
     call mpipack_size (m, mpir8, mpicom, bsiz)
     write(6,*) 'SPMDBUF: Allocating SPMD buffers of size ',bsiz
     allocate(buf1(bsiz/8+1))
     allocate(buf2(bsiz/8+1))
     return
  end subroutine spmdbuf

  subroutine compute_gsfactors (numperlat, numtot, numperproc, displs)
!----------------------------------------------------------------------- 
! 
! Purpose: Compute arguments for gatherv, scatterv
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
     integer, intent(in) :: numperlat    ! number of elements per latitude
!
! Output arguments
!
     integer, intent(out) :: numtot               ! total number of elements (to send or recv)
     integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
     integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements
!
! Local variables
!
     integer :: p                    ! index
   
     numtot = numperlat*numlats
   
     do p=0,npes-1
        numperproc(p) = numperlat*nlat_p(p)
     end do
     
     displs(0) = 0
     do p=1,npes-1
        displs(p) = displs(p-1) + numperproc(p-1)
     end do
     
  end subroutine compute_gsfactors


!############################(wanhui)##################################

!!subroutine compute_gsfactors (numperlat, numtot, numperproc, displs)
  subroutine compute_gsfactors_dyn (numhor, numver, numtot, numperproc, displs)
!----------------------------------------------------------------------- 
! 
! Purpose: Compute arguments for gatherv, scatterv
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
!!   integer, intent(in) :: numperlat    ! number of elements per latitude
     integer, intent(in) :: numhor
     integer, intent(in) :: numver
!
! Output arguments
!
     integer, intent(out) :: numtot               ! total number of elements (to send or recv)
     integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
     integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements
!
! Local variables
!
     integer :: p                    ! index
   
!!   numtot = numperlat*numlats
     numtot = numhor*numver
   
     do p=0,npes-1
!!      numperproc(p) = numperlat*nlat_p(p)
        numperproc(p) = plond*( nlat_p(p)+numlatsex-numlats )*numver
     end do
     
     displs(0) = 0
     do p=1,npes-1
        displs(p) = displs(p-1) + numperproc(p-1)
     end do
     
!!end subroutine compute_gsfactors
  end subroutine compute_gsfactors_dyn

!############################(wanhui)##################################

#endif

end module spmd_dyn
