!****f* FHI-aims/stratmann_weight_restricted
!  NAME
!   stratmann_weight_restricted
!  SYNOPSIS

real*8 function stratmann_weight_restricted( &
                                            i_center,       &
                                            dist_tab,       &
                                            a,              &
                                            atom_atom_tab,  &
                                            n_atom_list,    &
                                            n_compute_atoms,&
                                            center_index,   &
                                            i_atom_list,    &
                                            atom_list,      &
                                            write_out)
!  PURPOSE
!  calculate the partition weight of a single integration point using the 
!  stratmann partition tab. Core routine for the Stratmann formalism
!
!  The "restriction" modifies the original stratmann weight. Regardless of 
!  the structure, atoms will now no longer modify the partition weight if
!  their distance to the current integration point is larger than the
!  "outer radius" of the atom (here, the radius of the free-atom density
!  after application of the cutoff potential). Thus, atoms whose basis functions, 
!  densities etc do not contribute at this point will also not pick up an
!  integrand here.
!
!  We note that this restriction might break things for systems whose radial 
!  functions are more extended than the original free atom. In principle, the
!  proper outer radius would have to be the maximum outermost radius of each
!  basis function. 
!
!  This should still be implemented. - VB, 18.8.2011
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
  integer, dimension(n_compute_atoms) :: center_index

  ! VB: Variable write_out does not do anything. Kept here only in order to be able
  !     to access the routine in a simple way for test output writing in the future.
  !     Will not harm production.
  integer :: write_out

  real*8 :: g_actual
  real*8 :: g_modified
  real*8 :: edge

  real*8, external :: smooth_partition_edge
  character*50 :: file_name
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
!    o   stratmann_weight_restricted -- the weight of a point if it belonged to atom i_center
!  SOURCE
  ! local variables
  real*8                          :: g_stratmann_local, part_weight, mu
  integer                         :: i_compute_atom

  part_weight = 1d0

  do i_compute_atom = 1, i_center-1 
     mu =  (dist_tab(i_center)- &
          dist_tab(i_compute_atom)) &
          /atom_atom_tab(i_atom_list,atom_list(i_compute_atom))
     g_actual = g_stratmann_local(mu,a)
     edge = smooth_partition_edge( dist_tab(i_compute_atom), species_center (center_index(i_compute_atom)) )
     g_modified = (1.d0-edge) + edge*g_actual
     part_weight = part_weight*g_modified

  enddo

  do i_compute_atom = i_center+1, n_compute_atoms 
     mu =  (dist_tab(i_center)- &
          dist_tab(i_compute_atom)) &
          /atom_atom_tab(i_atom_list,atom_list(i_compute_atom))
     g_actual = g_stratmann_local(mu,a)
     edge = smooth_partition_edge( dist_tab(i_compute_atom), species_center (center_index(i_compute_atom)) )
     g_modified = (1.d0-edge) + edge*g_actual
     part_weight = part_weight*g_modified

  end do

  stratmann_weight_restricted = part_weight
end function stratmann_weight_restricted
!******

!****f* FHI-aims/g_stratmann
!  NAME
!   g_stratmann
!  SYNOPSIS

real*8 function g_stratmann_local(x,a)

  !  PURPOSE
  !  evaluate the "g-function" in the original stratmann paper
  !
  !  Labelled "local" here only to avoid confusion with the version included in
  !  "stratmann_partition_weight.f90", the original version of the code.
  !  g_stratmann itself is identical to g_stratmann_local .
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
  g_stratmann_local = (1d0-(1d0-dabs(leftright))*arg*(35d0+buf*(-35d0+buf*(21d0-5d0*buf)))/16d0-leftright)/2d0
end function g_stratmann_local
!******

!****f* FHI-aims/stratmann_weight_restricted_LN
!  NAME
!   stratmann_weight_restricted_LN
!  SYNOPSIS

real*8 function stratmann_weight_restricted_LN( &
       i_center,       &
       dist_tab,       &
       a,              &
       n_atom_atom_tab,&
       atom_atom_dist_list, &
       atom_idx_A,     &
       atom_idx_B,     &
       n_atom_list,    &
       n_compute_atoms,&
       center_index,   &
       i_atom_list,    &
       atom_list,      &
       write_out)

!  PURPOSE
!  calculate the partition weight of a single integration point using the 
!  stratmann partition tab. Core routine for the Stratmann formalism
!
!  The "restriction" modifies the original stratmann weight. Regardless of 
!  the structure, atoms will now no longer modify the partition weight if
!  their distance to the current integration point is larger than the
!  "outer radius" of the atom (here, the radius of the free-atom density
!  after application of the cutoff potential). Thus, atoms whose basis functions, 
!  densities etc do not contribute at this point will also not pick up an
!  integrand here.
!
!  We note that this restriction might break things for systems whose radial 
!  functions are more extended than the original free atom. In principle, the
!  proper outer radius would have to be the maximum outermost radius of each
!  basis function. 
!
!  Using the memory saving version of the atom_atom_tab, namely atom_atom_dist_list
!  This should still be implemented. - VB, 18.8.2011
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
  use localorb_io, only: localorb_info

!  ARGUMENTS 

  implicit none
  ! input
  integer                         :: n_atom_list, n_compute_atoms
  integer                         :: i_center, i_atom_list
  integer, dimension(n_atom_list) :: atom_list
  real*8 , dimension(n_atom_list) :: dist_tab
  real*8                          :: a
  integer, dimension(n_compute_atoms) :: center_index
  integer, intent(in)                :: n_atom_atom_tab
  real*8, dimension(n_atom_atom_tab) :: atom_atom_dist_list
  integer, dimension(n_atom_atom_tab):: atom_idx_A
  integer, dimension(n_atom_list+1):: atom_idx_B


  integer :: write_out

  real*8 :: g_actual
  real*8 :: g_modified
  real*8 :: edge

  real*8, external :: smooth_partition_edge

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
!    o   stratmann_weight_restricted_LN -- the weight of a point if it belonged to atom i_center
!  SOURCE
  ! local variables
  real*8                          :: g_stratmann_local, part_weight, mu
  integer                         :: i_compute_atom
  integer                         :: i_atom_idx_idx
  real*8, dimension(3) :: delta
  real*8 :: dist
  character*180 :: info_str

  part_weight = 1d0
  dist = 1d0

  i_atom_idx_idx = atom_idx_B(i_atom_list)

  outer: do i_compute_atom =  1, i_center-1
    do while (atom_idx_A(i_atom_idx_idx) < atom_list(i_compute_atom))
        i_atom_idx_idx = i_atom_idx_idx + 1
        if (i_atom_idx_idx >= atom_idx_B(i_atom_list+1)) then
            exit outer
        end if
     end do

     if (atom_idx_A(i_atom_idx_idx) .eq. atom_list(i_compute_atom)) then
        dist = atom_atom_dist_list(i_atom_idx_idx)
     else
        ! missing pair in the list of interatomic distances - fallback!
        delta(1) = coords_center(1,center_index(i_center))-coords_center(1,center_index(i_compute_atom))
        delta(2) = coords_center(2,center_index(i_center))-coords_center(2,center_index(i_compute_atom))
        delta(3) = coords_center(3,center_index(i_center))-coords_center(3,center_index(i_compute_atom))

        dist = sqrt(delta(1)*delta(1)+delta(2)*delta(2)+delta(3)*delta(3) )
!    LN: In the case that the chosen cut-off radius for the atom_atom_tab is to small, the needed interatomic
!    distance is calculated to give the right  partition weight and a warning is written in the output file.
         write(info_str,'(2X,A,1X,I5,1X,A,1X,I5,1X,A,1X,F15.8,1X,A)') &
            "* Warning! The entry for atom", atom_idx_A(i_atom_idx_idx) , &
            "and atom", atom_list(i_compute_atom) , &
            "is missing in atom_atom_tab. It will be calculated on the fly. Distance:", &
            dist*bohr, "AA."
        call localorb_info(info_str)
     end if

     mu =  ( dist_tab(i_center)- &
             dist_tab(i_compute_atom)) &
             / dist
     g_actual = g_stratmann_local(mu,a)
     edge = smooth_partition_edge( dist_tab(i_compute_atom), species_center (center_index(i_compute_atom)) )
     g_modified = (1.d0-edge) + edge*g_actual
     part_weight = part_weight*g_modified
  end do outer

  outer2: do i_compute_atom =  i_center+1, n_compute_atoms
    do while (atom_idx_A(i_atom_idx_idx) < atom_list(i_compute_atom))
        i_atom_idx_idx = i_atom_idx_idx + 1
        if (i_atom_idx_idx >= atom_idx_B(i_atom_list+1)) then
            exit outer2
        end if
     end do

     if (atom_idx_A(i_atom_idx_idx) .eq. atom_list(i_compute_atom)) then
        dist = atom_atom_dist_list(i_atom_idx_idx)
     else
        ! missing pair in the list of interatomic distances - fallback!
        delta(1) = coords_center(1,center_index(i_center))-coords_center(1,center_index(i_compute_atom))
        delta(2) = coords_center(2,center_index(i_center))-coords_center(2,center_index(i_compute_atom))
        delta(3) = coords_center(3,center_index(i_center))-coords_center(3,center_index(i_compute_atom))

        dist = sqrt(delta(1)*delta(1)+delta(2)*delta(2)+delta(3)*delta(3) )
!    LN: In the case that the chosen cut-off radius for the atom_atom_tab is to small, the needed interatomic
!    distance is calculated to give the right  partition weight and a warning is written in the output file.
         write(info_str,'(2X,A,1X,I5,1X,A,1X,I5,1X,A,1X,F15.8,1XA)') &
            "* Warning! The entry for atom", atom_idx_A(i_atom_idx_idx) , &
            "and atom", atom_list(i_compute_atom) , &
            "is missing in atom_atom_tab. It will be calculated on the fly. Distance:", &
            dist*bohr, "AA."
        call localorb_info(info_str)
     end if

     mu =  ( dist_tab(i_center)- &
             dist_tab(i_compute_atom)) &
             / dist
     g_actual = g_stratmann_local(mu,a)
     edge = smooth_partition_edge( dist_tab(i_compute_atom), species_center (center_index(i_compute_atom)) )
     g_modified = (1.d0-edge) + edge*g_actual
     part_weight = part_weight*g_modified
  end do outer2


  stratmann_weight_restricted_LN = part_weight
end function stratmann_weight_restricted_LN

!******
