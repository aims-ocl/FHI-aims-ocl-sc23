!****s* FHI-aims/RRS-PBC/get_rrs_pbc_cell_vector()
!  NAME
!    get_rrs_pbc_cell_vector() 
!  SYNOPSIS

    subroutine get_rrs_pbc_cell_vector(i_center,i_equal,n_vec)

!  PURPOSE
!  This routine get the cell index of the equal atoms
!
!  USES

      use dimensions
      use geometry
      use runtime_choices
      use mpi_tasks
      implicit none

!  ARGUMENTS

!  INPUTS
!    none
!  OUTPUT
!    none
!  COPYRIGHT
!   
!   
!   
!  HISTORY
!   
!  SOURCE
!

      ! imported variables
      integer                :: i_center, i_equal
      real*8, dimension(3)   :: n_vec


      ! local variables

      n_vec(:) = coords(:,rrs_pbc_center_atom(1,i_center)) - &
                 coords(:,rrs_pbc_equal_atom(1,i_equal,i_center))

      end subroutine
