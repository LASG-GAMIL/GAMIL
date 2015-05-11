#include <misc.h>
#include <params.h>

subroutine bndexch

!----------------------------------------------------------------------- 
! 
! Purpose: Pack and Exchange initial prognostic information among all the 
!          processors
!
! Method: 
! 
! Author: 
! 
!-----------------------------------------------------------------------
! $Id: bndexch.F90,v 1.5 2001/10/05 16:39:41 boville Exp $
! $Author: boville $
!
!----------------------------Parameters---------------------------------

#ifdef SPMD
  use pmgrid, only: iam
  use spmd_dyn, only: npes, cut, cutex, neighs, neighn

  implicit none
!
! Local workspace
!
  integer ns, nn
  integer inreg( 2 )
  integer outreg( 2 )
  integer others,othern   ! Other node
!
! Return if number of processors is less than 2
!
  if (npes .lt. 2) return
!
! For each partition (south and north) communicate boundaries
! on each side of partition among however many neighbors necessary
!
! send south, receive north
!
  ns = 1
  nn = 1
  do while (ns .le. neighs .or. nn .le. neighn)
     if (ns .le. neighs) then
        others = iam - ns
!
! Intersection of my cuts and neighbor processor's extended
! cuts tells if this node needs to send data to neighbor 
!
        call intersct(cut(1,iam),cutex(1,others),outreg)  
        ns = ns + 1
     else
        others = -1
        outreg(1) = 0
        outreg(2) = 0
     end if

     if (nn .le. neighn) then
        othern = iam + nn
!
! Intersection of neighbor cuts and this node's extended
! cut tells if this node recieves data from neighbor 
!
        call intersct(cut(1,othern),cutex(1,iam),inreg)  
        nn = nn + 1
     else
        othern = -1
        inreg(1) = 0
        inreg(2) = 0
     end if

     call bndexch_mpi(others,outreg,othern,inreg)
  end do

!
! send north, receive south
!
  ns = 1
  nn = 1
  do while (ns .le. neighs .or. nn .le. neighn)
     if (nn .le. neighn) then
        othern = iam + nn
!
! Intersection of my cuts and neighbor processor's extended
! cuts tells if this node needs to send data to neighbor 
!
        call intersct(cut(1,iam),cutex(1,othern),outreg)  
        nn = nn + 1
     else
        othern = -1
        outreg(1) = 0
        outreg(2) = 0
     end if

     if (ns .le. neighs) then
        others = iam - ns
!
! Intersection of neighbor cuts and this node's extended
! cut tells if this node recieves data from neighbor 
!
        call intersct(cut(1,others),cutex(1,iam),inreg)  
        ns = ns + 1
     else
        others = -1
        inreg(1) = 0
        inreg(2) = 0
     end if

     call bndexch_mpi(othern,outreg,others,inreg)
  end do
#endif
  return
end subroutine bndexch

#ifdef SPMD
subroutine bndexch_mpi(othero,outreg,otheri,inreg)
!-----------------------------------------------------------------------
! Send initial prognostic information to my peer process
!-----------------------------------------------------------------------
  use pmgrid, only: plat, plndlv, j1
  use constituents, only: pcnst
  use prognostics,  only: u3, v3, qminus, n3m1
  use mpishorthand

  implicit none

  integer, parameter :: msgtype = 6000
  integer, parameter :: j1m = j1 - 1
  integer, parameter :: siz = (2 + pcnst)*plndlv

  integer othero,outreg(2),otheri,inreg(2)
  integer num
  integer msg

  integer reqs(3*(plat+1))
  integer stats(MPI_STATUS_SIZE, 3*(plat+1))

  integer reqr(3*(plat+1))
  integer statr(MPI_STATUS_SIZE, 3*(plat+1))

  integer i,j
  integer reqs_i,reqr_i

  reqs_i = 0
  if (othero .ne. -1) then
     do i = outreg(1), outreg(2)
        j = 3*(i-outreg(1))

        msg = msgtype + j
        reqs_i = reqs_i + 1
        call mpiisend (u3(1,1,j1m+i,n3m1),plndlv,mpir8, othero,msg,mpicom,reqs(reqs_i))

        msg = msgtype + j + 1
        reqs_i = reqs_i + 1
        call mpiisend (v3(1,1,j1m+i,n3m1),plndlv,mpir8, othero,msg,mpicom,reqs(reqs_i))

        msg = msgtype + j + 2
        reqs_i = reqs_i + 1
        num = pcnst*plndlv
        call mpiisend (qminus(1,1,1,j1m+i),num,mpir8, othero,msg,mpicom,reqs(reqs_i))
             
     end do
  end if

  reqr_i = 0
  if (otheri .ne. -1) then
     do i = inreg(1), inreg(2)
        j = 3*(i-inreg(1))
        msg = msgtype + j
        reqr_i = reqr_i + 1
        call mpiirecv (u3(1,1,j1m+i,n3m1),plndlv,mpir8, otheri,msg,mpicom,reqr(reqr_i))

        msg = msgtype + j + 1
        reqr_i = reqr_i + 1
        call mpiirecv (v3(1,1,j1m+i,n3m1),plndlv,mpir8, otheri,msg,mpicom,reqr(reqr_i))

        msg = msgtype + j + 2
        reqr_i = reqr_i + 1
        num = pcnst*plndlv
        call mpiirecv (qminus(1,1,1,j1m+i),num,mpir8, otheri,msg,mpicom,reqr(reqr_i))
             
     end do
  end if

  if (reqs_i .ne. 0) then
     call mpiwaitall(reqs_i,reqs,stats)
  end if

  if (reqr_i .ne. 0) then
     call mpiwaitall(reqr_i,reqr,statr)
  end if

  return
end subroutine bndexch_mpi

subroutine intersct (regiona, regionb, regionc)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! Given two regions (a,b) output the intersection (common latitudes)  
! of these two sets.  The routine is used in bndexch to determine which
! latitudes need to be communicated to neighboring processors.  Typically
! this routine is invoked as the intersection of the set of resident 
! latitudes on processor A with the set of extended latitudes (needed for 
! the SLT) of processor B.  Any common latitudes will need to be 
! communicated to B to complete SLT processing. 
! 
! Author: 
! Original version:  CCM2
! Standardized:      J. Rosinski, Oct 1995
!                    J. Truesdale, Feb. 1996
!-----------------------------------------------------------------------
!
! $Id: bndexch.F90,v 1.5 2001/10/05 16:39:41 boville Exp $
! $Author: boville $
!
!----------------------------Commons------------------------------------
  implicit none
!
!---------------------------Local workspace-----------------------------
!
  integer regiona( 2 ),regionb( 2 ),regionc( 2 )
!
!-----------------------------------------------------------------------
!
  regionc( 1 ) = max( regiona( 1 ), regionb( 1 ) )
  regionc( 2 ) = min( regiona( 2 ), regionb( 2 ) )

  return
end subroutine intersct
#endif
