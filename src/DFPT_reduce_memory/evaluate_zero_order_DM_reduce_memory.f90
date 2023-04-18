!****s* FHI-aims/evaluate_zero_order_DM_reduce_memory
!  NAME
!    evaluate_zero_order_DM_reduce_memory
!  SYNOPSIS

subroutine evaluate_zero_order_DM_reduce_memory &
           (KS_eigenvector, occ_numbers, max_occ_number,density_matrix)
               

!  PURPOSE
!    calculate the zero-order DM  

!  USES

  use dimensions
  use pbc_lists
  use synchronize_mpi

!  ARGUMENTS

  implicit none
  
  real*8, dimension(n_basis, n_states, n_spin,n_k_points), intent(IN) :: KS_eigenvector
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers
  integer, dimension(n_spin), intent(IN) :: max_occ_number

  real*8, dimension(n_basis, n_basis),intent(INOUT) :: density_matrix

!  INPUTS
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates

!
!  OUTPUT
!  density_matrix


  integer :: i_basis, i_state  
  real*8, dimension(n_basis, n_states) :: KS_scaled

  real*8 :: occu
  integer :: info,  n_occ_states

!--------------initialize----------------------       
  density_matrix=0.0d0
  n_occ_states=max_occ_number(1)
 
  do i_state = 1, n_states, 1
     occu =  occ_numbers(i_state, 1, 1)
     if (occu .gt. 0.d0) then
           KS_scaled(:, i_state) = KS_eigenvector(:, i_state, 1, 1) * sqrt(occu)
     endif
  end do

  !-----here DM = Cui^(*) Cvi = Ci_aims * Ci_aims^(diagger)  the same as PM_none-----  
  call dsyrk('U', 'N', n_basis, n_occ_states, 1.d0, KS_scaled, n_basis, &
  &          0.d0, density_matrix, n_basis)

  !----zherk only give upder part, so we need add the down park
  density_matrix = density_matrix + transpose(density_matrix)
  do i_basis = 1, n_basis
     density_matrix(i_basis, i_basis) = density_matrix(i_basis, i_basis)*0.5d0
  end do

end subroutine evaluate_zero_order_DM_reduce_memory
!******
