!****s* FHI-aims/prune_density_matrix
!  NAME
!   prune_density_matrix
!  SYNOPSIS

subroutine prune_density_matrix(density_matrix, density_matrix_con, &
     n_compute, i_basis)

!  PURPOSE
!   The subroutine saves the density matrix components corresponding to the 
!   non-zero basis functions to the density_matrix_con. This subroutine works
!   for non-packed matrix format.
!
!  USES

  use dimensions
  implicit none

!  ARGUMENTS

  real*8 :: density_matrix(n_centers_basis_T, n_centers_basis_T)
  real*8 :: density_matrix_con(n_compute, n_compute)
  integer:: n_compute
  integer:: i_basis(n_compute)

!  INPUTS
!  o density_matrix -- total density matrix
!  o n_compute -- number of non-zero basis functions in this grid batch
!  o i_basis -- the list of the non-zero basis functions.
!
!  OUTPUT
!  o density_matrix_con -- the values of the density matrix corresponding to
!                          non-zero basis functions.
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


  integer :: i_compute_1, i_compute_2

 !do i_compute = 1, n_compute, 1

 !   density_matrix_con(1:i_compute, i_compute) = &
 !     density_matrix(i_basis(1:i_compute),i_basis(i_compute))

 !end do

  do i_compute_2 = 1, n_compute, 1
     do i_compute_1 = 1, n_compute, 1
        density_matrix_con(i_compute_1, i_compute_2) = &
          density_matrix(i_basis(i_compute_1),i_basis(i_compute_2))
     enddo
  end do

end subroutine prune_density_matrix
!******
