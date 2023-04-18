!****s* FHI-aims/prune_force_atoms_v1
!  NAME
!   prune_force_atoms_v1
!  SYNOPSIS

subroutine prune_force_atoms_v1(n_compute, i_basis, global_atom, &
     basis_offset, n_compute_force_atoms)

!  PURPOSE
! The subroutine determines atoms for which force contributions are to be evaluated in this shell.
! It also determines offset of local basis index for all relevant basis functions u(r)_alpha/r*Y_lm(theta,phi) 
! with respect to current atom.
!
!  USES

  use basis
  use dimensions
  implicit none

!  ARGUMENTS


  integer, intent(in) :: n_compute
  integer, intent(in) :: i_basis(n_basis)

  integer, intent(out) :: global_atom(n_atoms)
  integer, intent(out) :: basis_offset(n_atoms+1)
  integer, intent(out) :: n_compute_force_atoms


!  INPUTS
! o n_compute -- number of relevant basis functions
! o i_basis -- list of relevant basis functions
!  
!  OUTPUT
! o global_atom -- list of relevant atoms
! o basis_offset -- starting point of local spline
! o n_compute_force_atoms -- number of relevant atoms
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


	




  ! local variables
  integer :: i_basis_1
  integer :: i_atom

  ! counter
  integer :: i_compute
  integer :: i_offset
  integer :: i_compute_force_atom

  global_atom(1) = basis_atom(i_basis(1))
  n_compute_force_atoms = 1
  basis_offset(1) = 1
  i_offset = 1
  do i_compute = 2, n_compute, 1
     i_atom = basis_atom(i_basis(i_compute))
     ! found next atom ?
     if (i_atom .gt. global_atom(n_compute_force_atoms)) then
        n_compute_force_atoms = n_compute_force_atoms + 1
        global_atom(n_compute_force_atoms) = i_atom
        basis_offset(n_compute_force_atoms) = basis_offset(n_compute_force_atoms - 1) + i_offset
        i_offset = 1
     else
        i_offset = i_offset + 1
     end if
  end do
  basis_offset(n_compute_force_atoms + 1) = n_compute + 1 

end subroutine prune_force_atoms_v1
!******	
