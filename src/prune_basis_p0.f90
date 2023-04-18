!****s* FHI-aims/prune_basis_p0
!  NAME
!   prune_basis_p0
!  SYNOPSIS

subroutine prune_basis_p0 &
     ( dist_tab_sq, n_compute_a, n_compute_c, i_basis, &
     n_basis_list, n_atom_list, inv_list )

!  PURPOSE
!     reduces full basis to relevant functions only, for a given set of integration points
!
!  USES

  use basis, only: basis_fn, outer_radius_sq
  use dimensions, only: n_centers, n_basis
  use pbc_lists, only: cbasis_to_center, cbasis_to_basis
  implicit none

!  ARGUMENTS


  !  imported variables

  !     input
  integer :: n_basis_list, n_atom_list
  integer :: inv_list(n_centers)
  real*8   ::dist_tab_sq(n_atom_list)

  !     output

  integer :: n_compute_c
  integer :: n_compute_a
  integer :: i_basis(n_basis_list)


!  INPUTS
!    o n_basis_list -- number of basis functions (includid periodic mirror images)
!    o n_atom_list -- number of atoms (includid periodic mirror images)
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
  integer :: i_compute_atom
  integer :: i_compute_1

  logical :: flag = .false.

  ! begin work

  ! tabulate total wave function value for each basis function

  if (n_compute_c.eq.0) then
     ! this is the first point of the batch - simply list all non-zero wave functions

     i_compute_center = 0
     i_compute_atom = 0

     do i_basis_1 = 1, n_basis_list, 1


        if (dist_tab_sq(inv_list(Cbasis_to_center(i_basis_1))) .le. &
             outer_radius_sq(basis_fn(Cbasis_to_basis(i_basis_1)))) then

           ! nonzero basis function - use it ... 
           i_compute_center = i_compute_center+1
           i_basis(i_compute_center) = i_basis_1


           if(i_basis_1 <= n_basis)then
              i_compute_atom = i_compute_atom + 1
           end if
        end if

     enddo

     n_compute_c = i_compute_center
     n_compute_a = i_compute_atom

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

           do i_compute_1 = n_compute_c, i_compute_center, -1
              i_basis(i_compute_1+1) = i_basis(i_compute_1)
           enddo

           i_basis(i_compute_center) = i_basis_1

           n_compute_c = n_compute_c+1
           if(i_basis_1 <= n_basis)then
              n_compute_a = n_compute_a + 1
           end if

        end if

     enddo

  end if



  n_compute_a = n_compute_c

end subroutine prune_basis_p0
!---------------------------------------------------------------------
!******
subroutine prune_basis_p0X &
     ( dist_tab_sq, n_compute_a, n_compute_c, i_basis, &
     n_basis_list, n_atom_list, inv_list, newcopy )

!  PURPOSE
!     reduces full basis to relevant functions only, for a given set of integration points
!
!  USES

  use basis, only: basis_fn, outer_radius_sq, outradsq_save, invlist_save
  use dimensions, only: n_centers
  use pbc_lists, only: cbasis_to_center, cbasis_to_basis
  implicit none

!  ARGUMENTS


  !  imported variables

  !     input
  integer :: n_basis_list, n_atom_list
  integer :: inv_list(n_centers)
  real*8   ::dist_tab_sq(n_atom_list)
!NEC_CB
  integer :: newcopy


  !     output

  integer :: n_compute_c
  integer :: n_compute_a
  integer :: i_basis(n_basis_list)


!  INPUTS
!    o n_basis_list -- number of basis functions (includid periodic mirror images)
!    o n_atom_list -- number of atoms (includid periodic mirror images)
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
  integer :: i_compute_atom
  integer :: i_compute_1

  logical :: flag = .false.

  ! begin work

  if (newcopy.eq.1) then
    if (allocated(outradsq_save)) then
      deallocate(outradsq_save)
      deallocate(invlist_save)
    end if
    allocate(outradsq_save(n_basis_list))
    allocate(invlist_save(n_basis_list))
    do i_basis_1 = 1, n_basis_list, 1
      outradsq_save(i_basis_1)=outer_radius_sq(basis_fn(Cbasis_to_basis(i_basis_1)))
      invlist_save(i_basis_1)=inv_list(Cbasis_to_center(i_basis_1))
    end do
  end if

  ! tabulate total wave function value for each basis function

  if (n_compute_c.eq.0) then
     ! this is the first point of the batch - simply list all non-zero wave functions

     i_compute_center = 0
     i_compute_atom = 0

     do i_basis_1 = 1, n_basis_list, 1


        if (dist_tab_sq(invlist_save(i_basis_1)) .le. &
             outradsq_save(i_basis_1)) then
!NEC_CB        if (dist_tab_sq(inv_list(Cbasis_to_center(i_basis_1))) .le. &
!NEC_CB             outer_radius_sq(basis_fn(Cbasis_to_basis(i_basis_1)))) then

           ! nonzero basis function - use it ... 
           i_compute_center = i_compute_center+1
           i_basis(i_compute_center) = i_basis_1


!NEC_CB Not used - overwritten later
!NEC_CB           if(i_basis_1 <= n_basis)then
!NEC_CB              i_compute_atom = i_compute_atom + 1
!NEC_CB           end if
        end if

     enddo
     n_compute_c = i_compute_center
     n_compute_a = i_compute_atom

  else
     ! this is not the first integration point of the batch - check whether non-zero
     ! wave functions are already there, else add them to the list

     i_compute_center = 0

     do i_basis_1 = 1, n_basis_list, 1

        if (i_basis(i_compute_center+1).eq.i_basis_1) then
           ! basis function already in list

           i_compute_center = i_compute_center+1

        else if (dist_tab_sq(invlist_save(i_basis_1)) .le. &
             outradsq_save(i_basis_1)) then
!NEC_CB        else if (dist_tab_sq(inv_list(Cbasis_to_center(i_basis_1))) .le. &
!NEC_CB             outer_radius_sq(basis_fn(Cbasis_to_basis(i_basis_1)))) then
           ! basis function not in list - add it to the list of nonzero functions

           i_compute_center = i_compute_center+1

           do i_compute_1 = n_compute_c, i_compute_center, -1
              i_basis(i_compute_1+1) = i_basis(i_compute_1)
           enddo
           i_basis(i_compute_center) = i_basis_1

           n_compute_c = n_compute_c+1
!NEC_CB Not used - overwritten later
!NEC_CB           if(i_basis_1 <= n_basis)then
!NEC_CB              n_compute_a = n_compute_a + 1
!NEC_CB           end if

        end if

     enddo

  end if



  n_compute_a = n_compute_c

end subroutine prune_basis_p0X
!---------------------------------------------------------------------
!******
