module fspan

  use pmgrid, only:  plond,beglatexdyn, endlatexdyn
  use mpi_gamil
  use infnan

  implicit none

  integer,allocatable :: mdj (:)
  integer,allocatable :: mm1 (:,:)
  integer,allocatable :: mp1 (:,:)
  integer,allocatable :: mm2 (:,:)
  integer,allocatable :: mp2 (:,:)
  integer,allocatable :: mm3 (:,:)
  integer,allocatable :: mp3 (:,:)

contains

  subroutine initialize_fspan
      
    integer i, j

    allocate(mdj(beglatexdyn:endlatexdyn))
    allocate(mm1(beglonex:endlonex,beglatexdyn:endlatexdyn))
    allocate(mp1(beglonex:endlonex,beglatexdyn:endlatexdyn))
    allocate(mm2(beglonex:endlonex,beglatexdyn:endlatexdyn))
    allocate(mp2(beglonex:endlonex,beglatexdyn:endlatexdyn))
    allocate(mm3(beglonex:endlonex,beglatexdyn:endlatexdyn))
    allocate(mp3(beglonex:endlonex,beglatexdyn:endlatexdyn))

    mdj(:)   = bigint
    mm1(:,:) = bigint
    mp1(:,:) = bigint
    mm2(:,:) = bigint
    mp2(:,:) = bigint
    mm3(:,:) = bigint
    mp3(:,:) = bigint

    do j = jbeg0, jend0
      mdj(j) = 1
    end do

    do j = jbeg0, jend0
      do i = beglonex, endlonex
	      mm1(i,j) = i - mdj(j)
	      mp1(i,j) = i + mdj(j)
	      mm2(i,j) = i - (mdj(j) - 1) / 2
	      mp2(i,j) = i + (mdj(j) + 1) / 2
	      mm3(i,j) = i - (mdj(j) + 1) / 2
	      mp3(i,j) = i + (mdj(j) - 1) / 2
      end do
    end do

  end subroutine initialize_fspan

end module fspan
