!****s* FHI-aims/add_zora_matrix_p1
!  NAME
!    add_zora_matrix_p1
!  SYNOPSIS
 
subroutine add_zora_matrix_p1 &
     ( zora_vector1, zora_vector2, n_basis_list, n_rel_points, &
     n_compute, hamiltonian_shell )

!  PURPOSE
!   The subroutine calculates the zora part of the hamiltonian evaluation
!   the vector vector product for several grid points at once. The results
!   are added to hamiltonian_shell
!
!  USES

      use dimensions
      use basis
      implicit none

!  ARGUMENTS

      integer :: n_basis_list
      integer :: n_rel_points
      real*8  :: zora_vector1(n_basis_list, 3, n_rel_points)
      real*8  :: zora_vector2(n_basis_list, 3, n_rel_points)
      integer :: n_compute
      real*8  :: hamiltonian_shell( n_compute,n_compute)

!  INPUTS
! o  n_basis_list -- total number of basis functions including periodic mirror images
! o  n_rel_points -- number of grid points relevant to relativistic calculations
! o  zora_vector1 -- first zora vector
! o  zora_vector2 -- second zora vector
! o  n_compute -- number of non-zero basis functions in current grid batch
!
!  OUTPUT
! o  hamiltonian_shell -- the results of the multiplications are added here.
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




!     Now Integrate:

      ! Notice: This will multiply a tightly packed zora vector in all three dimensions at once.

         call dgemm('N','T', n_compute, n_compute, 3*n_rel_points, &
              1.0d0, zora_vector1, n_basis_list,  &
              zora_vector2, n_basis_list, & 
              1.0d0, hamiltonian_shell, n_compute) 
         
       end subroutine add_zora_matrix_p1
!---------------------------------------------------------------------
!******
