!! (wanhui 2003.10.24)


module fspan

!!  use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
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

CONTAINS

   subroutine initialize_fspan
      
      allocate (mdj (beglatexdyn:endlatexdyn))
      allocate (mm1 (beglonex:endlonex,beglatexdyn:endlatexdyn))
      allocate (mp1 (beglonex:endlonex,beglatexdyn:endlatexdyn))
      allocate (mm2 (beglonex:endlonex,beglatexdyn:endlatexdyn))
      allocate (mp2 (beglonex:endlonex,beglatexdyn:endlatexdyn))
      allocate (mm3 (beglonex:endlonex,beglatexdyn:endlatexdyn))
      allocate (mp3 (beglonex:endlonex,beglatexdyn:endlatexdyn))


      mdj (:)   = bigint 
      mm1 (:,:) = bigint 
      mp1 (:,:) = bigint 
      mm2 (:,:) = bigint 
      mp2 (:,:) = bigint 
      mm3 (:,:) = bigint 
      mp3 (:,:) = bigint 

      return
   end subroutine initialize_fspan

end module fspan
