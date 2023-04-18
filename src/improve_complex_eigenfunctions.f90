!****s* FHI-aims/improve_complex_eigenfunctions
!  NAME
!   improve_complex_eigenfunctions
!  SYNOPSIS

subroutine improve_complex_eigenfunctions(overlap_matrix,hamiltonian,&
   KS_eigenvalue,KS_eigenvector,i_k_point)

!  PURPOSE
!    The subroutine calculates the complex eigenvectors and eigenvalues.
!    First the singularity of overlap matrix is tested and the eigenvectors
!    of overlapmatrix which belong to too small eigenvalues are removed.
!
!  USES

   use dimensions
   use lapack_wrapper
   use runtime_choices

   implicit none

!  ARGUMENTS

   complex*16 :: overlap_matrix(n_basis*(n_basis+1)/2)
   complex*16 :: hamiltonian(n_basis*(n_basis+1)/2)

   integer :: i_k_point

   real*8 :: KS_eigenvalue(n_states)
   complex*16 :: KS_eigenvector(n_basis,n_states)

!  INPUTS
!  o overlap_matrix -- overlap matrix
!  o hamiltonian -- Hamiltonian matrix
!  o i_k_point -- k-point index
!
!  OUTPUT
!  o KS_eigenvalue -- eigenvalues
!  o KS_eigenvector -- eigenvectors
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

   complex*16, dimension(:), allocatable, save :: ovlp_work
   complex*16, dimension(:,:), allocatable, save :: ovlp_work_full
   complex*16, dimension(:,:), allocatable, save :: hamiltonian_work

   integer, save :: n_nonsingular
   real*8, dimension(:), allocatable, save :: ovlp_eigenvalues
   complex*16, dimension(:,:), allocatable, save :: ovlp_transform
   complex*16, dimension(:), allocatable, save :: trafo_hamiltonian
   real*8 :: ev_sqrt
   complex*16, dimension(:,:), allocatable, save :: aux_eigenvector

   complex*16, parameter :: one = (1.d0,0.d0)
   complex*16, parameter :: zero = (0.d0,0.d0)

   integer :: i_basis_1
   integer :: i_basis_2
   integer :: i_index
   integer :: dim_work
   integer :: i_trafo_1

!  LAPACK version

!  first, check whether the basis is singular
   if (flag_KS_k_points(i_k_point) == BASIS_SINGULAR_NOT_TESTED) then

!    Eigenvalue decomposition of the overlap matrix
     if (.not. allocated(ovlp_eigenvalues)) then
       allocate(ovlp_eigenvalues(n_basis))
     end if
     if (.not. allocated(ovlp_transform)) then
       allocate(ovlp_transform(n_basis,n_basis))
     end if
     if (.not. allocated(ovlp_work)) then
       allocate(ovlp_work(n_basis*(n_basis+1)/2))
     end if

     ovlp_work = overlap_matrix

     call diagonalize_overlap_complex(n_basis,ovlp_work,safe_minimum,&
          basis_threshold,n_nonsingular,ovlp_eigenvalues,ovlp_transform,&
          i_k_point)

     n_states_k(i_k_point) = n_nonsingular

     if (n_nonsingular == n_basis) then
!    Basis is not singular - proceed with Cholesky decomposition instead

       flag_KS_k_points(i_k_point) = BASIS_NON_SINGULAR

       if (allocated(ovlp_eigenvalues)) then
          deallocate(ovlp_eigenvalues)
       end if
       if (allocated(ovlp_transform)) then
          deallocate(ovlp_transform)
       end if

     else
!    Basis is singular - use simultaneous diagonalisation of overlap, Hamiltonian

       flag_KS_k_points(i_k_point) = BASIS_SINGULAR

       do i_trafo_1 = 1,n_nonsingular
         ev_sqrt = sqrt(ovlp_eigenvalues(i_trafo_1))

         do i_basis_1 = 1,n_basis
           ovlp_transform(i_basis_1,i_trafo_1)&
              = ovlp_transform(i_basis_1,i_trafo_1)/ev_sqrt
         end do
       end do

     end if

   end if

   if (flag_KS_k_points(i_k_point) == BASIS_NON_SINGULAR) then

     n_states_k(i_k_point) = n_states

     if (.not. allocated(ovlp_work_full)) then
       allocate(ovlp_work_full(n_basis,n_basis))
     end if
     if (.not. allocated(hamiltonian_work)) then
       allocate(hamiltonian_work(n_basis,n_basis))
     end if

     i_index = 0

     do i_basis_2 = 1,n_basis
       do i_basis_1 = 1,i_basis_2

         i_index = i_index + 1
         ovlp_work_full(i_basis_1,i_basis_2) = overlap_matrix(i_index)
         hamiltonian_work(i_basis_1,i_basis_2) = hamiltonian(i_index)

       end do
     end do

     dim_work = 8*n_basis+10000

     call complex_lapack_solver(n_basis,n_states,ovlp_work_full,&
          hamiltonian_work,safe_minimum,dim_work,KS_eigenvalue,KS_eigenvector)

     deallocate(ovlp_work_full)
     deallocate(hamiltonian_work)

   else if (flag_KS_k_points(i_k_point) == BASIS_SINGULAR) then

     ! Eigenvalue decomposition of the overlap matrix
     if (.not. allocated(ovlp_eigenvalues)) then
       allocate(ovlp_eigenvalues(n_basis))
     end if
     if (.not. allocated(ovlp_transform)) then
       allocate(ovlp_transform(n_basis,n_basis))
     end if
     if (.not. allocated(ovlp_work)) then
       allocate(ovlp_work(n_basis*(n_basis+1)/2))
     end if

     ovlp_work = overlap_matrix

     call diagonalize_overlap_complex(n_basis,ovlp_work,safe_minimum,&
          basis_threshold,n_nonsingular,ovlp_eigenvalues,ovlp_transform,&
          i_k_point)

     n_states_k(i_k_point) = n_nonsingular

     do i_trafo_1 = 1,n_nonsingular
       ev_sqrt = sqrt(ovlp_eigenvalues(i_trafo_1))

       do i_basis_1 = 1,n_basis
         ovlp_transform(i_basis_1,i_trafo_1)&
            = ovlp_transform(i_basis_1,i_trafo_1)/ev_sqrt
       end do
     end do

!    Solve generalized eigenvalue problem by simultaneous diagonalisation,
!    using LAPACK subroutines. This approach stabilizes the overlap
!    matrix against singularities before solving the generalized EVP,
!    but I did not ensure that the implementation is maximally efficient.
!    Must check one day if there is time.

!    Transform Hamiltonian

!    This transformation is probably not ideal, but it works for now, and does not cost
!    too much time. Nevertheless, it's the biggest time-consumer
!    for the simultaenous diagonalization.

     if (.not.allocated(trafo_hamiltonian)) then
       allocate(trafo_hamiltonian(n_nonsingular*(n_nonsingular+1)/2))
     end if

     call packed_matrix_transform_complex(hamiltonian,n_basis,ovlp_transform,&
          n_nonsingular,trafo_hamiltonian)

!    now do a standard diagonalization of the transformed Hamiltonian

!  Note:  With minor effort, we should be able to avoid the auxiliary array
!         aux_eigenvector, and only use KS_eigenvector instead. This would,
!         once again, save memory.
     if (.not. allocated(aux_eigenvector)) then
       allocate(aux_eigenvector(n_nonsingular,n_states))
     end if

     call diagonalize_hamiltonian_complex(n_nonsingular,n_states,&
          trafo_hamiltonian,safe_minimum,KS_eigenvalue,aux_eigenvector)

!    finally, transform auxiliary eigenvectors back to underlying basis functions

     KS_eigenvector = zero

     call zgemm("N","N",n_basis,n_states,n_nonsingular,one,ovlp_transform,&
          n_basis,aux_eigenvector,n_nonsingular,zero,KS_eigenvector,n_basis)

     ! must deallocate since we don't know if the next call will have the same
     ! value of n_nonsigular (maybe from a different k-point)
     deallocate(trafo_hamiltonian)
     deallocate(aux_eigenvector)

   else

     call localorb_info(" * Unknown update method for KS eigenstates.")
     stop

   end if

   if (allocated(ovlp_work)) then
     deallocate(ovlp_work)
   end if

end subroutine improve_complex_eigenfunctions
!----------------------------------------------------------------------
!******
