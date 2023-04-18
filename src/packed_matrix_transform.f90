!****s* FHI-aims/packed_matrix_transform
!  NAME
!   packed_matrix_transform
!  SYNOPSIS


subroutine packed_matrix_transform &
     ( hamiltonian, n_basis, ovlp_transform, n_nonsingular, &
     trafo_hamiltonian &
     )

!  PURPOSE
!  Subroutine packed_matrix_transform is a feeble attempt at creating
!  a faster transform than a straighforward multiplication of all
!  components, using BLAS. The problem is that I do not know a fast
!  way of handling this transformation.
!
!  The goal is to perform a matrix transformation
!
!    B = T'*A*T
!
!  where B and A are both real, symmetric, and packed matrices, of order
!  n x n and m x m, respectively, while T is a real but general n x m matrix
!
!  USES

  implicit none

!  ARGUMENTS

  integer n_basis
  integer n_nonsingular
  real*8, dimension(n_basis*(n_basis+1)/2)             :: hamiltonian
  real*8, dimension(n_basis, n_nonsingular)            :: ovlp_transform
  real*8, dimension(n_nonsingular*(n_nonsingular+1)/2) :: trafo_hamiltonian

! INPUTS
! o n_basis -- number of basis functions
! o hamiltonian -- Hamiltonian matrix
! o n_nonsingular -- number of non-sigular basis functions
! o ovlp_transform -- ????????
!
! OUTPUT
! o trafo_hamiltonian -- ???????????
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


	




  !  local variables

  real*8, dimension(n_basis) :: column
  real*8 :: one = 1.0d0
  real*8 :: zero = 0.0d0

  !     counters

  integer i_trafo_1, i_trafo_2
  integer i_basis_1, i_basis_2
  integer i_index

  !     functions

  real*8 ddot

  !  begin work

  i_index = 0
  do i_trafo_2 = 1, n_nonsingular, 1
     column = 0.d0
     call dspmv('U',n_basis,one,hamiltonian, &
          ovlp_transform(1,i_trafo_2),1,zero, &
          column,1)

     do i_trafo_1 = 1, i_trafo_2,1
        i_index = i_index+1
        trafo_hamiltonian(i_index) = &
             ddot(n_basis,ovlp_transform(1,i_trafo_1),1, &
             column,1)
     enddo
  enddo

  !  end work

end subroutine packed_matrix_transform
!---------------------------------------------------------------------
!******	
