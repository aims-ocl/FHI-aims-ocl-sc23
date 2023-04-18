!****s* FHI-aims/evaluate_first_order_DM_polar_reduce_memory
!  NAME
!    evaluate_first_order_DM_polar_reduce_memory
!  SYNOPSIS

subroutine evaluate_first_order_DM_polar_reduce_memory(   &
           first_order_H, &
           KS_eigenvalue, occ_numbers,  &
           first_order_density_matrix, n_matrix_size)
               

!  PURPOSE
!    calculate the first-order DM for polarizability 

!  shanghui,2013.12.12
!  USES

  use dimensions
  use runtime_choices

  use scalapack_wrapper, only : construct_first_order_ham_polar_reduce_memory_scalapack, & 
                                evaluate_first_order_U_polar_reduce_memory_scalapack, &
                                construct_first_order_dm_polar_reduce_memory_scalapack, &
                                get_first_order_dm_polar_reduce_memory_scalapack, &
                                get_set_full_local_matrix_scalapack_cpscf, &
                                print_ham_cpscf, &
                                first_order_ham_polar_reduce_memory_scalapack, &
                                first_order_U_polar_reduce_memory_scalapack, &
                                evaluate_first_order_U_polar_reduce_memory_scalapack_cpscf
  use synchronize_mpi_basic, only: sync_vector

  ! wyj
  use load_balancing
  use localorb_io

!  ARGUMENTS

  implicit none
  

  real*8, dimension(n_matrix_size, n_spin),intent(INOUT) :: first_order_H
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers
  real*8, dimension(n_matrix_size, n_spin),intent(INOUT) :: first_order_density_matrix
  ! wyj
  integer :: n_matrix_size

  !  INPUTS
  !   o occ_numbers -- occupations of eigenstates
  !   o KS_eigenvalue -- Kohn-Sham eigenvalues

  !
  !  OUTPUT
  !  first_order_density_matrix

  integer :: i_spin
  character*200 :: info_str

  if (use_local_index .and. use_load_balancing) then
      write(info_str,'(2X,A)') "Evaluating first-order-DM matrix with local_index and load_balancing"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)

      !---(1) first transform dense first-order H lapack to dense first-order H scalapack
      ! call construct_first_order_ham_polar_reduce_memory_scalapack(first_order_H)
       !print *, 'print_ham_cpscf'
       !call print_ham_cpscf(first_order_ham_polar_reduce_memory_scalapack)

      !---(2) evaluate first-order U scalapack
      !
      call evaluate_first_order_U_polar_reduce_memory_scalapack_cpscf(occ_numbers, KS_eigenvalue, 2)
      ! wyj: [PASS]
      !print *, 'print_U_cpscf'
      !call print_ham_cpscf(first_order_U_polar_reduce_memory_scalapack)


      !---(3)  Construct the first-order density matrix (dense, scalapack format)
      !
      call construct_first_order_dm_polar_reduce_memory_scalapack(occ_numbers) 
      ! wyj: [PASS]
      !print *, 'print_ham_cpscf_2'
      !call print_ham_cpscf(first_order_ham_polar_reduce_memory_scalapack)

      !---(4)  Put the density matrix back in natural form (n_basis dim) so that it could be used in the following
      !call get_first_order_dm_polar_reduce_memory_scalapack( first_order_density_matrix)
      do i_spin = 1, n_spin
      call get_set_full_local_matrix_scalapack_cpscf(first_order_density_matrix, 2, i_spin)
      enddo

      !---(5) Synchronize
      !call sync_vector(first_order_density_matrix,n_hamiltonian_matrix_size*n_spin)
  else
      write(info_str,'(2X,A)') "Evaluating first-order-DM matrix"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)

      !---(1) first transform dense first-order H lapack to dense first-order H scalapack
      call construct_first_order_ham_polar_reduce_memory_scalapack(first_order_H)
      ! wyj: debug [DONE]
      !print *, 'print_ham_cpscf'
      !call print_ham_cpscf(first_order_ham_polar_reduce_memory_scalapack)

      !---(2) evaluate first-order U scalapack
      !call evaluate_first_order_U_polar_reduce_memory_scalapack(occ_numbers,KS_eigenvalue)
      call evaluate_first_order_U_polar_reduce_memory_scalapack_cpscf(occ_numbers, KS_eigenvalue, 1)
      !print *, 'print_U_cpscf'
      !call print_ham_cpscf(first_order_U_polar_reduce_memory_scalapack)


      !---(3)  Construct the first-order density matrix (dense, scalapack format)
      call construct_first_order_dm_polar_reduce_memory_scalapack(occ_numbers) 
      !print *, 'print_ham_cpscf_2'
      !call print_ham_cpscf(first_order_ham_polar_reduce_memory_scalapack)

      !---(4)  Put the density matrix back in natural form (n_basis dim) so that it could be used in the following
      call get_first_order_dm_polar_reduce_memory_scalapack( first_order_density_matrix)

      !---(5) Synchronize
      call sync_vector(first_order_density_matrix,n_hamiltonian_matrix_size*n_spin)
  endif

end subroutine evaluate_first_order_DM_polar_reduce_memory
!******
