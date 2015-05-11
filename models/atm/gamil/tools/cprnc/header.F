module header
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! Netcdf id
!
  integer, dimension(2) :: ncid = -1
!
! Dimensions
!
  integer, dimension(2) :: londimid = -1
  integer, dimension(2) :: latdimid = -1
  integer, dimension(2) :: levdimid = -1
  integer, dimension(2) :: ilevdimid = -1
  integer, dimension(2) :: unlimdimid = -1
!
! Dimension variables
!
  integer, dimension(2) :: lonid = -1
  integer, dimension(2) :: latid = -1
  integer, dimension(2) :: levid = -1
  integer, dimension(2) :: ilevid = -1
  integer, dimension(2) :: timeid = -1

  real(r8), allocatable :: lon(:,:), lat(:,:), lev(:,:), ilev(:,:)
  real(r8), allocatable :: time(:), time2(:)
!
! Character variables
!
  character*16, dimension(2) :: case
  character*80, dimension(2) :: title

  integer, dimension(2) :: date_writtenid
  integer, dimension(2) :: time_writtenid

  character(len=8), allocatable :: date_written(:,:)
  character(len=8), allocatable :: time_written(:,:)
!
! Integer variables
!
  integer, dimension(2) :: ntrmid = -1, ntrm
  integer, dimension(2) :: ntrnid = -1, ntrn
  integer, dimension(2) :: ntrkid = -1, ntrk
  integer, dimension(2) :: ndbaseid = -1, ndbase
  integer, dimension(2) :: nsbaseid = -1, nsbase
  integer, dimension(2) :: nbdateid = -1, nbdate
  integer, dimension(2) :: nbsecid = -1, nbsec
  integer, dimension(2) :: mdtid = -1, mdt
  integer, dimension(2) :: ntime

  integer, dimension(2) :: nstephid = -1
  integer, dimension(2) :: ncdateid = -1
  integer, dimension(2) :: ncsecid = -1

  integer, allocatable :: nsteph(:,:)
  integer, allocatable :: ncdate(:,:)
  integer, allocatable :: ncsec(:,:)
!
! Real variables
!
  integer, dimension(2) :: psid = -1
  integer, dimension(2) :: p0id = -1
  real(r8), dimension(2) :: p0

  integer, dimension(2) :: gwid = -1
  integer, dimension(2) :: hyaiid = -1
  integer, dimension(2) :: hybiid = -1
  integer, dimension(2) :: hyamid = -1
  integer, dimension(2) :: hybmid = -1

  real(r8), allocatable :: gw(:,:), hyai(:,:), hybi(:,:), hyam(:,:), hybm(:,:)

end module header
