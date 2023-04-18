!****s* FHI-aims/evaluate_zero_order_EDM_reduce_memory
!  NAME
!    evaluate_zero_order_EDM_reduce_memory
!  SYNOPSIS

subroutine evaluate_zero_order_EDM_reduce_memory &
           (KS_eigenvector, KS_eigenvalue, occ_numbers, max_occ_number, energy_density_matrix)
               

!  PURPOSE
!    calculate the zero-order DM  

!  USES

  use dimensions
  use pbc_lists
  use synchronize_mpi

!  ARGUMENTS

  implicit none
  
  real*8, dimension(n_basis, n_states, n_spin,n_k_points), intent(IN) :: KS_eigenvector
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers
  integer, dimension(n_spin), intent(IN) :: max_occ_number

  real*8, dimension(n_basis, n_basis),intent(INOUT) :: energy_density_matrix

!  INPUTS
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates

!
!  OUTPUT
!  density_matrix


  integer :: i_basis, i_state  
  real*8, dimension(n_basis, n_states) :: KS_eigen_times_occnumber

  real*8 :: occu
  integer :: info,  n_occ_states

!--------------initialize----------------------       
  energy_density_matrix=0.0d0
  n_occ_states=max_occ_number(1)
 

 do i_state=1, n_states
     KS_eigen_times_occnumber(1:n_basis,i_state)=   &
     KS_eigenvector(1:n_basis,i_state,1,1)*occ_numbers(i_state,1,1)  &
     * KS_eigenvalue(i_state,1,1)  
  enddo


  CALL dgemm("N","T",n_basis,n_basis,n_basis,&
             1.0d0, KS_eigen_times_occnumber, n_basis,KS_eigenvector, n_basis,&
             0.0d0,energy_density_matrix,n_basis)


end subroutine evaluate_zero_order_EDM_reduce_memory
!******
