!****f* FHI-aims/evaluate_soc_shell
!*  NAME
!*    evaluate_soc_shell
!*  SYNOPSIS
subroutine evaluate_soc_shell(n_points, partition, n_compute, potential_shell, &
                              gradient_basis_wave, soc_shell, i_pauli, i_sign)
!*  PURPOSE
!*    Calculate a selected piece of the batch matrix of the SOC operator for the
!*    current batch
!*  USES
  use dimensions, only: n_max_batch_size, n_max_compute_ham
  implicit none
!*  ARGUMENTS
  integer, intent(in) :: n_points
  real*8, intent(in) :: partition(n_max_batch_size)
  integer, intent(in) :: n_compute
  real*8, intent(in) :: potential_shell(n_max_batch_size)
  real*8, intent(in) :: &
       gradient_basis_wave(n_max_compute_ham, 3, n_max_batch_size)
  real*8, intent(out) :: soc_shell(n_compute, n_compute)
  integer, intent(in) :: i_pauli
  integer, intent(in) :: i_sign
!*  INPUTS
!*    o n_points            - Number of points in the current batch
!*    o partition           - Integration weights for current batch
!*    o n_compute           - Number of basis functions touching current batch
!*    o potential_shell     - Potential evaluated on real-space grid for current batch
!*    o gradient_basis_wave - Gradients for basis functions on current batch
!*    o i_pauli             - Selects the real-space operation whose matrix elements
!*                            we are evalulating.  They are indexed by the Pauli matrix
!*                            they couple to, hence the name.
!*    o i_sign              - Which of the two terms in the cross product we are
!*                            evaluating in the current call.  Since SOC operator has the
!*                            form of a cross product, each real-space operation has two
!*                            terms, one with a plus sign and one with a minus sign.
!*  OUTPUTS
!*    o soc_shell           - Selected piece of the batch matrix of the real-space SOC operator
!*                            for current batch
!*  AUTHORS
!*    William Huhn and Matthias Gramzow
!*  HISTORY
!*    ???            - Created
!*    September 2017 - Updated to reduce size of soc_shell and support
!*                     load balancing
!*  NOTES
!*    The implementation of second-variational SOC in FHI-aims is published in
!*      Huhn and Blum, Phys. Rev. Materials 1, 033803 (2017)
!*      https://dx.doi.org/10.1103/PhysRevMaterials.1.033803
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!*    the terms and conditions of the respective license agreement."
!*  SOURCE

  ! local variables
  integer ::  i_compute_1, i_compute_2, i_point
  real*8, dimension(n_compute, 3, n_points) :: &
       gradient_compute_a, gradient_compute_c

  do i_point = 1, n_points, 1
    gradient_compute_a(1:n_compute, 1:3, i_point) = &
         partition(i_point) * potential_shell(i_point) * &
         gradient_basis_wave(1:n_compute, 1:3, i_point)
    gradient_compute_c(1:n_compute, 1:3, i_point) = &
         gradient_basis_wave(1:n_compute, 1:3, i_point)
  end do

  ! Computes soc_shell.  soc_shell(:,:,i,j) is the matrix element that will
  ! couple to the Pauli matrix sigma_{x_{i}} with sign +1 for i_sign=1, -1 for
  ! i_sign=-1
  ! These signs will be assigned when updating the soc matrix
  if (i_pauli .eq. 3) then
    if (i_sign .eq. 1) then
      ! WPH:  soc(:,:,3,1) = <dbra/dy|partition*hartree|dket/dx> = ...
      call dgemm('N', 'T', &
                 n_compute, n_compute, n_points, &
                 1.d0, &
                 gradient_compute_a(:,2,:), n_compute, &
                 gradient_compute_c(:,1,:), n_compute, &
                 0.d0, &
                 soc_shell, n_compute)
    else if (i_sign .eq. -1) then
      ! WPH:  soc(:,:,3,2) = dbra/dx|partition*hartree|dket/dy> = ...
      call dgemm('N', 'T', &
                 n_compute, n_compute, n_points, &
                 1.d0, &
                 gradient_compute_a(:,1,:), n_compute, &
                 gradient_compute_c(:,2,:), n_compute, &
                 0.d0, &
                 soc_shell, n_compute)
    end if
  else if (i_pauli .eq. 2) then
    if (i_sign .eq. 1) then
      ! WPH:  soc(:,:,2,1) = <dbra/dx|partition*hartree|dket/dz>
      call dgemm('N', 'T', &
                 n_compute, n_compute, n_points, &
                 1.d0, &
                 gradient_compute_a(:,1,:), n_compute, &
                 gradient_compute_c(:,3,:), n_compute, &
                 0.d0, &
                 soc_shell, n_compute)
    else if (i_sign .eq. -1) then
      ! WPH:  soc(:,:,2,2) = <dbra/dz|partition*hartree|dket/dx>
      call dgemm('N', 'T', &
                 n_compute, n_compute, n_points, &
                 1.d0, &
                 gradient_compute_a(:,3,:), n_compute, &
                 gradient_compute_c(:,1,:), n_compute, &
                 0.d0, &
                 soc_shell, n_compute)
    end if
  else if (i_pauli .eq. 1) then
    if (i_sign .eq. 1) then
      ! WPH:  soc(:,:,1,1) = <dbra/dz|partition*hartree|dket/dy>
      call dgemm('N', 'T', &
                 n_compute, n_compute, n_points, &
                 1.d0, &
                 gradient_compute_a(:,3,:), n_compute, &
                 gradient_compute_c(:,2,:), n_compute, &
                 0.d0, &
                 soc_shell, n_compute)
    else if (i_sign .eq. -1) then
      ! WPH:  soc(:,:,1,2) = <dbra/dy|partition*hartree|dket/dz>
      call dgemm('N', 'T', &
                 n_compute, n_compute, n_points, &
                 1.d0, &
                 gradient_compute_a(:,2,:), n_compute, &
                 gradient_compute_c(:,3,:), n_compute, &
                 0.d0, &
                 soc_shell, n_compute)
    end if
  end if
end subroutine evaluate_soc_shell
!****
