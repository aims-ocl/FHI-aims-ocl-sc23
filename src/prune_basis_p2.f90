!****s* FHI-aims/prune_basis_p2
!  NAME
!   prune_basis_p2
!  SYNOPSIS

subroutine prune_basis_p2 &
     ( dist_tab_sq, n_compute, i_basis,  n_basis_list, n_atom_list,  inv_list )


!  PURPOSE
!  Reduces full basis to relevant functions only, for a given set of integration points
!
!  USES

  use dimensions
  use basis
  use grids
  use geometry
  use pbc_lists
  implicit none

!  ARGUMENTS

  integer :: n_basis_list, n_atom_list
  integer :: inv_list(n_centers)
  real*8   ::dist_tab_sq(n_atom_list)
  integer :: n_compute
  integer :: i_basis(n_basis_list)

!  INPUTS
!    o n_basis_list -- number of basis functions (including periodic mirror images)
!    o n_atom_list -- number of atoms (including periodic mirror images)
!    o inv_list -- citing: atom center -> relative position in distance list
!    o dist_tab_sq -- (distance to atoms)**2
!   
!  OUTPUT
!    o n_compute_a == n_compute_c -- number of relevant basis functions
!    o i_basis -- list of relevant basis functions
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
!



  !  local variables

  !     counters

  integer :: i_basis_1
  integer :: i_compute_center
  integer :: i_compute_1

  logical :: flag = .false.



  ! tabulate total wave function value for each basis function

  if (n_compute.eq.0) then
     ! this is the first point of the batch - simply list all non-zero wave functions

     i_compute_center = 0

     do i_basis_1 = 1, n_basis_list, 1


        if (dist_tab_sq(inv_list(Cbasis_to_center(i_basis_1))) .le. &
             outer_radius_sq(basis_fn(Cbasis_to_basis(i_basis_1)))) then

           !  nonzero basis function - use it ... 
           i_compute_center = i_compute_center+1
           i_basis(i_compute_center) = i_basis_1

        end if

     enddo

     n_compute = i_compute_center

  else
     ! this is not the first integration point of the batch - check whether non-zero
     ! wave functions are already there, else add them to the list

     i_compute_center = 0

     do i_basis_1 = 1, n_basis_list, 1

        if (i_basis(i_compute_center+1).eq.i_basis_1) then
           ! basis function already in list

           i_compute_center = i_compute_center+1


        else if (dist_tab_sq(inv_list(Cbasis_to_center(i_basis_1))) .le. &
             outer_radius_sq(basis_fn(Cbasis_to_basis(i_basis_1)))) then
           ! basis function not in list - add it to the list of nonzero functions

           i_compute_center = i_compute_center+1

           do i_compute_1 = n_compute, i_compute_center, -1
              i_basis(i_compute_1+1) = i_basis(i_compute_1)
           enddo

           i_basis(i_compute_center) = i_basis_1

           n_compute = n_compute+1

        end if

     enddo

  end if

end subroutine prune_basis_p2
!---------------------------------------------------------------------
!******
