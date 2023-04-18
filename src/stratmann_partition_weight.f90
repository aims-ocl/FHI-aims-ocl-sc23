!****f* FHI-aims/stratmann_partition_weight
!  NAME
!   stratmann_partition_weight
!  SYNOPSIS

real*8 function stratmann_partition_weight(i_center,dist_tab,a,atom_atom_tab, &
     n_atom_list,n_compute_atoms,i_atom_list,atom_list)
!  PURPOSE
!  calculate the partition weight of a single integration point using the 
!  stratmann partition tab. Core routine for the Stratmann formalism
!
!  USES
!
  use dimensions
  use runtime_choices
  use grids
  use pbc_lists
  use spline
  use free_atoms
  use constants
!  ARGUMENTS 

  implicit none
  ! input
  integer                         :: n_atom_list, n_compute_atoms
  integer                         :: i_center, i_atom_list
  integer, dimension(n_atom_list) :: atom_list
  real*8 , dimension(n_atom_list) :: dist_tab
  real*8                          :: a
  real*8, dimension(n_atom_list, n_atom_list) :: atom_atom_tab
!  INPUTS
!    o   i_center -- atomic center index for dist_tab
!    o   dist_tab -- distances from grid point to all centers, reduced to the interesting subset available here
!    o   a -- magic partitioning parameter for weight function
!    o   atom_atom_tab -- interatomic distances for ALL atoms relevant in this system, MUST USE atom_list AS INDEX ARRAY!!!
!    o   n_atom_list --  index size for all of the lists
!    o   n_compute_atoms -- number of atoms actually used in the computation, < n_atom_list
!    o   i_atom_list -- atomic center index on atom_atom_tab
!    o   atom_list -- index list which maps all i_compute_atoms to their proper entries in atom_atom_tab
!  OUTPUT
!    o   stratmann_partition_weight -- the weight of a point if it belonged to atom i_center
!  SOURCE
  ! local variables
  real*8                          :: g_stratmann, part_weight, mu
  integer                         :: i_compute_atom
  part_weight = 1d0
  do i_compute_atom = 1, i_center-1 
     mu =  (dist_tab(i_center)- &
          dist_tab(i_compute_atom)) &
          /atom_atom_tab(i_atom_list,atom_list(i_compute_atom))
     part_weight = part_weight*g_stratmann(mu,a)
  enddo
  do i_compute_atom = i_center+1, n_compute_atoms 
     mu =  (dist_tab(i_center)- &
          dist_tab(i_compute_atom)) &
          /atom_atom_tab(i_atom_list,atom_list(i_compute_atom))
     part_weight = part_weight*g_stratmann(mu,a)
  end do
  stratmann_partition_weight = part_weight
end function stratmann_partition_weight
!******

!****f* FHI-aims/g_stratmann
!  NAME
!   g_stratmann
!  SYNOPSIS

real*8 function g_stratmann(x,a)

  !  PURPOSE
  !  evaluate the "g-function" in the original stratmann paper
  !
  !  USES
  !
  implicit none
  !  ARGUMENTS 
  real*8 :: x, a
  !  INPUTS
  !  o x -- normalized atom-point distance
  !  o a -- stratmann scaling parameter
  !  OUTPUT
  !  o g_stratmann -- stratmann g_function
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE
  real*8 :: buf, arg, leftright
  leftright =(sign(1d0,x+a)-sign(1d0,a-x))/2d0
  arg = x/a
  buf = arg*arg
  g_stratmann = (1d0-(1d0-dabs(leftright))*arg*(35d0+buf*(-35d0+buf*(21d0-5d0*buf)))/16d0-leftright)/2d0
end function g_stratmann
!******

