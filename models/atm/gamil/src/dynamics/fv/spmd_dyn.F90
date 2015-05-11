#include <misc.h>
#include <params.h>

module spmd_dyn
!BOP
!
! !MODULE: Subroutines to initialize SPMD implementation of CAM
!

#if (defined SPMD)

!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plat, plon, masterproc, iam, numbnd, &
                     numlats, beglat, endlat, &
                     plev, beglev, endlev, endlevp1, &
                     endlevp, myid_y, myid_z, npr_y, npr_z, plevp, &
                     myidxy_x, myidxy_y, nprxy_x, nprxy_y, &
                     beglonxy, endlonxy, beglatxy, endlatxy, twod_decomp
   use constituents, only: ppcnst
   use mpishorthand, only: mpir8, mpicom, mpiint
   use decompmodule, only: decomptype, decompcreate
   use redistributemodule, only: redistributetype, redistributecreate
   use infnan, only: inf

   implicit none

! !PUBLIC MEMBER FUNCTIONS:

   public spmdinit_dyn, decomp_wavenumbers
   public compute_gsfactors, spmdbuf

! !PUBLIC DATA MEMBERS:

   integer :: npes           ! Total number of MPI tasks
   integer, allocatable :: cut(:,:)   ! partition for MPI tasks
   integer proc(plat)        ! processor id associated with a given lat.
   integer, allocatable :: nlat_p(:)  ! number of latitudes per subdomain

   integer comm_y            ! communicator in latitude
   integer comm_z            ! communicator in vertical
   integer commxy_x          ! communicator in longitude (xy second. decomp.)
   integer commxy_y          ! communicator in latitude (xy second. decomp.)
   integer, allocatable :: lonrangexy(:,:)   ! global xy-longitude subdomain index
   integer, allocatable :: latrangexy(:,:)   ! global xy-latitude subdomain index

   type (redistributetype) :: inter_ijk, inter_ikj, inter_ijkp,     &
                              inter_ikjp, inter_q3

!
! !DESCRIPTION: 
!   {\bf Purpose:} Subroutines to initialize SPMD implementation of CAM
!
! !REVISION HISTORY:
!   ??.??.??  CCM Core Group     Creation
!   00.09.30  Sawyer             Alterations for LR SPMD mode
!   01.05.09  Mirin              2-D yz decomposition
!   01.06.27  Mirin              Secondary 2-D xy decomposition
!   01.12.20  Sawyer             Changed index order of Q3 decomposition
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: spmdinit_dyn --- SPMD initialization for dynamics
!
! !INTERFACE:

   subroutine spmdinit_dyn ()

! !USES:
      use parutilitiesmodule, only : parinit, parsplit
      use decompmodule, only : decompcreate
      use pmgrid, only: strip2d, strip3dxyz, strip3dxzy,                     &
                        strip3dxyzp, strip3dxzyp, strip3zaty,                &
                        strip3yatz, strip3yatzp,                             &
                        strip3kxyz, strip3kxzy, strip3kxyzp,                 &
                        strip3kxzyp, strip3zatypt, strip3zatyj2,             &
                        strip3zatypj1, strip3zatypj2, strip3zatyt4,          &
                        strip3dq3, strip3kq3, strip3dq3old

! !DESCRIPTION:
!
!   SPMD initialization routine: get number of cpus, processes, tids, etc
!
! !REVISION HISTORY:
!   ??.??.??  CCM Core Group     Creation
!   00.09.30  Sawyer             Added LR-specific initialization
!   01.03.26  Sawyer             Added ProTeX documentation
!   01.06.27  Mirin              Secondary 2-D xy decomposition
!   01.10.16  Sawyer             Added Y at each Z decompositions
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
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
      integer xdist(1)  ! number of lons per subdomain
      integer, allocatable :: ydist(:) ! number of lats per subdomain
      integer, allocatable :: zdist(:) ! number of levels per subdomain
      integer, allocatable :: zdistq(:) ! number of levels per subdomain for Q3
      integer ier       ! error flag
      integer rank_y, size_y   !  rank and size wrt y-communicator
      integer rank_z, size_z   !  rank and size wrt z-communicator
      integer rankxy_x, sizexy_x   !  rank and size wrt xy x-communicator
      integer rankxy_y, sizexy_y   !  rank and size wrt xy y-communicator
      integer zdist1(1) ! used for misc. decomposition definitions
      integer, allocatable :: ydistq(:) ! number of tracer/lats per subdomain
      integer, allocatable :: xdistxy(:) ! number of xy-longs per subdomain
      integer, allocatable :: ydistxy(:) ! number of xy-lats per subdomain
      integer, allocatable :: ydistqxy(:) ! number of xy tracer/lats per subdomain
      integer zdistxy(1)  ! number of xy-verts per subdomain
      integer j, k, vert, lonn
      integer ydistk(1)

! Namelist for 2D decomposition
#if defined( TWOD_YZ )
      integer npr_yz(4)
      namelist /mprun2d/ npr_yz
#endif

! Default 2D decomposition
      beglev = 1
      endlev = plev
      endlevp1 = plev + 1
      endlevp = plev + 1
      myid_y = iam
      myid_z = 0
      npr_y = npes
      npr_z = 1
      nprxy_x = 1
      nprxy_y = npes
      myidxy_x = 0
      myidxy_y = iam

! Read namelist for actual 2D decomposition
! Namelist is read on master only; information is broadcast to other tasks.
! Define task indexing for 2D decomposition
#if defined( TWOD_YZ )
      if (masterproc) then
        npr_yz(1) = npr_y
        npr_yz(2) = npr_z
        npr_yz(3) = nprxy_x
        npr_yz(4) = nprxy_y
        read (5,mprun2d)
        write (6,*) '2-D y-z decomposition for Lin-Rood dycore'
        write (6,*) 'npr_y = ', npr_yz(1), '  npr_z = ', npr_yz(2)
        write (6,*) 'nprxy_x= ', npr_yz(3), '  nprxy_y = ', npr_yz(4)
      endif
      call mpi_bcast(npr_yz, 4, mpiint, 0, mpicom, ier)
      npr_y = npr_yz(1)
      npr_z = npr_yz(2)
      nprxy_x = npr_yz(3)
      nprxy_y = npr_yz(4)
      if (masterproc) then
        if (npr_y*npr_z /= npes .or. nprxy_x*nprxy_y /= npes) then
           write (6,*) 'SPMDINIT_DYN : incorrect domain decomposition - aborting'
           call endrun
        endif
      end if
      myid_z = iam/npr_y
      myid_y = iam - myid_z*npr_y
      myidxy_y = iam/nprxy_x
      myidxy_x = iam - myidxy_y*nprxy_x
#endif

!
! Addition for LR dynamical core to initialize PILGRIM library
!
      call parinit(mpicom)
!
! Form separate communicators
!
      call parsplit(mpicom, myid_z, iam, comm_y, rank_y, size_y)
      call parsplit(mpicom, myid_y, iam, comm_z, rank_z, size_z)
      call parsplit(mpicom, myidxy_y, iam, commxy_x, rankxy_x, sizexy_x)
      call parsplit(mpicom, myidxy_x, iam, commxy_y, rankxy_y, sizexy_y)

!
!-----------------------------------------------------------------------
!
! Compute y decomposition
!
      allocate (ydist  (npr_y))
      allocate (ydistq (npr_y))
      allocate (nlat_p (0:npes-1))
      allocate (cut    (2,0:npes-1))

      ydist(:) = 0
      ydistq(:) = 0
      nlat_p(0:npes-1) = 0

      lat = plat / npr_y
      workleft = plat - lat * npr_y
      if ( lat < 3 ) then
         write(6,*)'SPMDINIT_DYN: less than 3 latitudes per subdomain'
         call endrun
      endif
!
! Be careful:  ydist is 1-based.  NCARs arrays, e.g., cut,  are 0-based
!
      do procid=1,npr_y
         ydist(procid) = lat
      enddo

      if ( workleft /= 0 ) then
         procids = (npr_y+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = npr_y
            ydist(procids) = ydist(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               ydist(procidn) = ydist(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(ydist) /= plat ) then
         write(6,*)'SPMDINIT_DYN:', ydist,' does not add up to ', plat
         call endrun
      endif

      if (workleft/=0) then
         write(6,*)'SPMDINIT_DYN: Workleft(y) not zero.  Value is ',workleft
         call endrun
      end if

! Set the NCAR data structures

      lat  = 0
      do procid=0,npr_y-1
         cut(1,procid) = lat+1
         lat = lat + ydist(procid+1)
         cut(2,procid) = lat
         nlat_p(procid) = ydist(procid+1)

         if (masterproc) then
            write(6,*) 'nlat_p(',procid,') = ', nlat_p(procid)
         end if

         if (myid_y == procid) then
            beglat  = cut(1,myid_y)
            endlat  = cut(2,myid_y)
            numlats = ydist(procid+1)
         end if
      enddo

      do k = 1, npr_z-1
         do j = 0, npr_y-1
            procid = j + k*npr_y
            cut(1,procid) = cut(1,j)
            cut(2,procid) = cut(2,j)
            nlat_p(procid) = nlat_p(j)
         enddo
      enddo
!
! Compute z decomposition
!
      allocate (zdist (npr_z))
      allocate (zdistq(npr_z))

      zdist(:) = 0

      vert = plev / npr_z
      workleft = plev - vert * npr_z
      if ( vert < 3 ) then
         write(6,*)'SPMDINIT_DYN: less than 3 verticals per subdomain'
         call endrun
      endif

      do procid=1,npr_z
         zdist(procid) = vert
      enddo

      if ( workleft /= 0 ) then
         procids = (npr_z+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = npr_z
            zdist(procids) = zdist(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               zdist(procidn) = zdist(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(zdist) /= plev ) then
         write(6,*)'SPMDINIT_DYN:', zdist,' does not add up to ', plev
         call endrun
      endif

      if (workleft/=0) then
         write(6,*)'SPMDINIT_DYN: Workleft(z) not zero.  Value is ',workleft
         call endrun
      end if

! Compute local limits

      beglev = 1
      endlev = zdist(1)
      do procid = 1, myid_z
         beglev = endlev + 1
         endlev = beglev + zdist(procid+1) - 1
      enddo
      endlevp1 = endlev + 1
      endlevp = endlev
      if (myid_z == npr_z-1) endlevp = endlev + 1

!
! Compute x secondary decomposition
!
      allocate (xdistxy (nprxy_x))

      xdistxy(:) = 0

      lonn = plon / nprxy_x
      workleft = plon - lonn * nprxy_x
      if ( lonn < 3 ) then
         write(6,*)'SPMDINIT_DYN: less than 3 xy-longitudes per subdomain'
         call endrun
      endif

      do procid=1,nprxy_x
         xdistxy(procid) = lonn
      enddo

      if ( workleft /= 0 ) then
         procids = (nprxy_x+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = nprxy_x
            xdistxy(procids) = xdistxy(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               xdistxy(procidn) = xdistxy(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(xdistxy) /= plon ) then
         write(6,*)'SPMDINIT_DYN:', xdistxy,' does not add up to ', plon
         call endrun
      endif

      if (workleft/=0) then
         write(6,*)'SPMDINIT_DYN: Workleft(xy-x) not zero.  Value is ',workleft
         call endrun
      end if

! Compute local limits

      beglonxy = 1
      endlonxy = xdistxy(1)
      do procid = 1, myidxy_x
         beglonxy = endlonxy + 1
         endlonxy = beglonxy + xdistxy(procid+1) - 1
      enddo

! Compute global table

      allocate (lonrangexy(2,nprxy_x))
      lonrangexy(1,1) = 1
      lonrangexy(2,1) = xdistxy(1)
      do procid = 2, nprxy_x
         lonrangexy(1,procid) = lonrangexy(2,procid-1) + 1
         lonrangexy(2,procid) = lonrangexy(1,procid) + xdistxy(procid) - 1
      enddo
!
! Compute y secondary decomposition
!
      allocate (ydistxy (nprxy_y))

      ydistxy(:) = 0

      lat = plat / nprxy_y
      workleft = plat - lat * nprxy_y
      if ( lat < 3 ) then
         write(6,*)'SPMDINIT_DYN: less than 3 xy-latitudes per subdomain'
         call endrun
      endif

      do procid=1,nprxy_y
         ydistxy(procid) = lat
      enddo

      if ( workleft /= 0 ) then
         procids = (nprxy_y+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = nprxy_y
            ydistxy(procids) = ydistxy(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               ydistxy(procidn) = ydistxy(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(ydistxy) /= plat ) then
         write(6,*)'SPMDINIT_DYN:', ydistxy,' does not add up to ', plat
         call endrun
      endif

      if (workleft/=0) then
         write(6,*)'SPMDINIT_DYN: Workleft(xy-y) not zero.  Value is ',workleft
         call endrun
      end if

! Compute local limits

      beglatxy = 1
      endlatxy = ydistxy(1)
      do procid = 1, myidxy_y
         beglatxy = endlatxy + 1
         endlatxy = beglatxy + ydistxy(procid+1) - 1
      enddo

! Compute global table

      allocate (latrangexy(2,nprxy_y))
      latrangexy(1,1) = 1
      latrangexy(2,1) = ydistxy(1)
      do procid = 2, nprxy_y
         latrangexy(1,procid) = latrangexy(2,procid-1) + 1
         latrangexy(2,procid) = latrangexy(1,procid) + ydistxy(procid) - 1
      enddo
!
! WS: create decompositions for NCAR data structures
!
      xdist(1) = plon
!
! Create PILGRIM decompositions (see decompmodule)
!
      call decompcreate( 1, npr_y, xdist, ydist, strip2d )
      call decompcreate( 1, npr_y, npr_z, xdist, ydist, zdist, strip3dxyz )
      call decompcreate( "xzy", 1, npr_z, npr_y, xdist, zdist, ydist, strip3dxzy )

! In q3 the tracer number and latitude are folded together

      ydistq(:) =  ppcnst * ydist(:)
      call decompcreate( "xzy", 1, npr_z, npr_y, xdist, zdist, ydistq, strip3dq3old )

      ydistq(:) =  plon * ydist(:)
      zdist1(1) =  ppcnst
      call decompcreate( npr_y, npr_z, 1, ydistq, zdist, zdist1, strip3dq3 )

! For y communication within z subdomain (klast version)
      zdist1(1) = endlev-beglev+1
      call decompcreate( 1, npr_y, 1, xdist, ydist, zdist1, strip3yatz )

! For z communication within y subdomain

      ydistk(1) = endlat-beglat+1
      call decompcreate( 1, 1, npr_z, xdist, ydistk, zdist, strip3zaty )
      ydistk(1) = endlat-beglat+3
      call decompcreate( 1, 1, npr_z, xdist, ydistk, zdist, strip3zatyj2 )

! For uv3s_update gathering

      ydistk(1) = 4*(endlat-beglat+1)
      call decompcreate( "xzy", 1, npr_z, 1, xdist, zdist, ydistk, strip3zatyt4 )

! Arrays dimensioned plev+1

      zdist(npr_z) = zdist(npr_z) + 1
      call decompcreate( 1, npr_y, npr_z, xdist, ydist, zdist, strip3dxyzp )
      call decompcreate( "xzy", 1, npr_z, npr_y, xdist, zdist, ydist, strip3dxzyp )

! Arrays dimensioned plev+1, within y subdomain

      ydistk(1) = endlat-beglat+1
      call decompcreate( "xzy", 1, npr_z, 1, xdist, zdist, ydistk, strip3zatypt )
      ydistk(1) = endlat-beglat+2
      call decompcreate( 1, 1, npr_z, xdist, ydistk, zdist, strip3zatypj1 )
      ydistk(1) = endlat-beglat+3
      call decompcreate( 1, 1, npr_z, xdist, ydistk, zdist, strip3zatypj2 )

! For y communication within z subdomain (klast+1 version)
      zdist1(1) = endlev-beglev+2
      call decompcreate( 1, npr_y, 1, xdist, ydist, zdist1, strip3yatzp )

! Secondary xy decomposition
!
      if (twod_decomp == 1) then
        zdistxy(1) = plev
        call decompcreate( nprxy_x, nprxy_y, 1, xdistxy, ydistxy, zdistxy, strip3kxyz )
        call decompcreate( "xzy", nprxy_x, 1, nprxy_y, xdistxy, zdistxy, ydistxy, strip3kxzy )

        zdist1(1)  = ppcnst * zdistxy(1)
        call decompcreate( "xyz", nprxy_x, nprxy_y, 1, xdistxy, ydistxy, zdist1, strip3kq3 )

        zdistxy(1) = zdistxy(1) + 1
        call decompcreate( nprxy_x, nprxy_y, 1, xdistxy, ydistxy, zdistxy, strip3kxyzp )
        call decompcreate( "xzy", nprxy_x, 1, nprxy_y, xdistxy, zdistxy, ydistxy, strip3kxzyp )

! Initialize transposes
!
        call redistributecreate(strip3dxyz, strip3kxyz, inter_ijk)
        call redistributecreate(strip3dxzy, strip3kxzy, inter_ikj)
        call redistributecreate(strip3dxyzp, strip3kxyzp, inter_ijkp)
        call redistributecreate(strip3dxzyp, strip3kxzyp, inter_ikjp)
        call redistributecreate(strip3dq3, strip3kq3, inter_q3)
      endif

!
! Do generic NCAR decomposition
!
      do procid=0,npes-1
         if (iam == 0) then
            write(6,*)'procid ',procid,' assigned ', &
                 cut(2,procid)-cut(1,procid)+1,' latitude values from', &
                 cut(1,procid),' through ',cut(2,procid)
         endif
!
! Determine which processor is responsible for the defined latitudes
!
         do lat=cut(1,procid),cut(2,procid)
            proc(lat) = procid
         end do
      end do
!
! Number of neighbor processors needed for boundary communication.  North
! first.
!
      isum = 0
      do procid=myid_y+1,npr_y-1
         nmostlat = cut(2,procid)
         isum = isum + cut(2,procid) - cut(1,procid) + 1
         if (isum >= numbnd) goto 20
      end do
20    if (myid_y /= npr_y-1 .and. isum < numbnd .and. nmostlat /= plat)then
         write (6,*) 'SPMDINIT_DYN: Something wrong in computation of northern neighbors'
         call endrun
      end if

      isum = 0
      do procid=myid_y-1,0,-1
         smostlat = cut(1,procid)
         isum = isum + cut(2,procid) - cut(1,procid) + 1
         if (isum >= numbnd) goto 30
      end do
30    if (myid_y /= 0 .and. isum < numbnd .and. smostlat /= 1) then
         write(6,*)'Something wrong in computation of southern neighbors'
         call endrun
      end if

!      write(6,*)'-----------------------------------------'
!      write(6,*)'Number of lats passed north & south = ',numbnd
!      write(6,*)'Node  Partition'
!      write(6,*)'-----------------------------------------'
!      do procid=0,npes-1
!         write(6,200) procid,cut(1,procid),cut(2,procid)
!      end do
!      write(6,*)'iam=',iam,'Number of south neighbors needed for bndry exchange = ',neighs
!      write(6,*)'iam=',iam,'Number of north neighbors needed for bndry exchange = ',neighn

      deallocate (ydist)
      deallocate (zdist)

      return
!
! Formats
!
200   format(i3,4x,i3,'-',i3,7x,i3,'-',i3)

!EOC
   end subroutine spmdinit_dyn

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
      implicit none
      
      write(6,*)'decomp_wavenumbers() should never be called in LR dynamics'
      call endrun

   end subroutine decomp_wavenumbers

   subroutine spmdbuf
!----------------------------------------------------------------------- 
! 
! Purpose: placeholder for buffer allocation routine 
! 
! Method: Make the labor division as equal as possible given loop lengths
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      implicit none
      
      write(6,*)'spmdbuf() should never be called in LR dynamics'
      call endrun

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
      integer, intent(in) :: numperlat            ! number of elements per latitude
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

#endif

end module spmd_dyn

