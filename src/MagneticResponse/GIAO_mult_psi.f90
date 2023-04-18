!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Multiplies H|n> from the left by <m| and performs the
!!  cross-product between R_mn and the matrix elements of H,
!!
!!    R_mn x <m|H|n>,
!!
!!  where R_mn = R_m-R_n is the difference of the position vectors of
!!  basis functions m and n, or performs a double cross-product of the form
!!
!!  R_mn x <m|H|n> x R_mn.
!!
!!  This subroutine is called for any H that results from the use of
!!  gauge-including atomic orbitals (GIAOs). The body of the operator
!!  H is found in integrands.f90.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
pure subroutine GIAO_mult_psi(compact_indexing, do_transpose, mult_R_mn_both, &
     & i_dim, n_dims, n_compute, n_points, wave, matrix_wave, r_mn_basis, &
     & matrix_batch_aux, matrix_batch)

  use tools, only: mul
  use types, only: dp

  implicit none

  logical, intent(in) :: compact_indexing, do_transpose, mult_R_mn_both
  integer, intent(in) :: i_dim, n_dims, n_compute, n_points
  real(dp), intent(in) :: wave(n_compute,n_points)
  ! 3rd dimension depends on the perturbation
  real(dp), intent(in) :: matrix_wave(n_compute,n_points,*)
  real(dp), intent(in) :: r_mn_basis(n_compute,n_compute,3)
  ! Auxiliary work array (3rd dim depends on the perturbation)
  real(dp), intent(in out) :: matrix_batch_aux(n_compute,n_compute,*)
  real(dp), intent(out) :: matrix_batch(n_compute,n_compute)
  integer :: i_dir
  if_mult_R_mn_both: if (mult_R_mn_both) then
     ! If the response tensor is in the form
     !
     !   chi = R x H x R^T,
     !
     ! where R is a vector and H is a 3x3 tensor, the individual
     ! elements are
     !
     !   chi_11 = R_2 (H_32 R_3-H_33 R_2) - R_3 (H_22 R_3-H_23 R_2)
     !   chi_21 = R_3 (H_12 R_3-H_13 R_2) - R_1 (H_32 R_3-H_33 R_2)
     !   chi_31 = R_1 (H_22 R_3-H_23 R_2) - R_2 (H_12 R_3-H_13 R_2)
     !   chi_12 = R_2 (H_33 R_1-H_31 R_3) - R_3 (H_23 R_1-H_21 R_3)
     !   chi_22 = R_3 (H_13 R_1-H_11 R_3) - R_1 (H_33 R_1-H_31 R_3)
     !   chi_32 = R_1 (H_23 R_1-H_21 R_3) - R_2 (H_13 R_1-H_11 R_3)
     !   chi_13 = R_2 (H_31 R_2-H_32 R_1) - R_3 (H_21 R_2-H_22 R_1)
     !   chi_23 = R_3 (H_11 R_2-H_12 R_1) - R_1 (H_31 R_2-H_32 R_1)
     !   chi_33 = R_1 (H_21 R_2-H_22 R_1) - R_2 (H_11 R_2-H_12 R_1)
     !
     ! Defining
     !
     !   A = H_32 R_3-H_33 R_2,
     !   B = H_22 R_3-H_23 R_2,
     !   C = H_12 R_3-H_13 R_2
     !
     ! We can rewrite the chi[:,1] part as
     !
     !   chi_11 = R_2 A - R_3 B
     !   chi_21 = R_3 C - R_1 A
     !   chi_31 = R_1 B - R_2 C
     !
     ! This way, the intermediate results are reused as much as
     ! possible. Next, redefine
     !
     !   A = H_33 R_1-H_31 R_3
     !   B = H_23 R_1-H_21 R_3
     !   C = H_13 R_1-H_11 R_3
     !
     ! and rewrite the chi[:,2] part as
     !
     !   chi_12 = R_2 A - R_3 B
     !   chi_22 = R_3 C - R_1 A
     !   chi_32 = R_1 B - R_2 C
     !
     ! Finally,
     !
     !   A = H_31 R_2-H_32 R_1
     !   B = H_21 R_2-H_22 R_1
     !   C = H_11 R_2-H_12 R_1
     !
     ! which leads to
     !
     !   chi_13 = R_2 A - R_3 B
     !   chi_23 = R_3 C - R_1 A
     !   chi_33 = R_1 B - R_2 C
     !
     if (compact_indexing) then
        if (i_dim == 1) then
           ! The contents of matrix_batch_aux is calculated during
           ! i_dim=1 and saved for the rest of the i_dim loop in
           ! GIAO_mult_psi_and_pack.
           do i_dir = 1, 6
              ! Compact indexing:
              ! 1 - 11;  4 - 12
              ! 2 - 22;  5 - 13
              ! 3 - 33;  6 - 23
              call mul_wrapper(matrix_wave(:,:,i_dir), &
                   & matrix_batch_aux(:,:,i_dir))
              if (do_transpose) then
                 matrix_batch_aux(:,:,i_dir) = matrix_batch_aux(:,:,i_dir) + &
                      & transpose(matrix_batch_aux(:,:,i_dir))
              end if
           end do
        end if
        comp_dim_9: if (n_dims == 9) then ! If the full tensor is requested
           select case(i_dim)
           case(1)
              matrix_batch_aux(:,:,7) = & ! A = H_32 R_3-H_33 R_2
                   & matrix_batch_aux(:,:,6)*r_mn_basis(:,:,3) - &
                   & matrix_batch_aux(:,:,3)*r_mn_basis(:,:,2)
              matrix_batch_aux(:,:,8) = & ! B = H_22 R_3-H_23 R_2,
                   & matrix_batch_aux(:,:,2)*r_mn_basis(:,:,3) - &
                   & matrix_batch_aux(:,:,6)*r_mn_basis(:,:,2)
              matrix_batch_aux(:,:,9) = & ! C = H_12 R_3-H_13 R_2
                   & matrix_batch_aux(:,:,4)*r_mn_basis(:,:,3) - &
                   & matrix_batch_aux(:,:,5)*r_mn_basis(:,:,2)
              ! chi_11 = R_2 A - R_3 B
              matrix_batch = r_mn_basis(:,:,2)*matrix_batch_aux(:,:,7) - &
                   & r_mn_basis(:,:,3)*matrix_batch_aux(:,:,8)
           case(2)
              ! chi_21 = R_3 C - R_1 A
              matrix_batch = r_mn_basis(:,:,3)*matrix_batch_aux(:,:,9) - &
                   & r_mn_basis(:,:,1)*matrix_batch_aux(:,:,7)
           case(3)
              ! chi_31 = R_1 B - R_2 C
              matrix_batch = r_mn_basis(:,:,1)*matrix_batch_aux(:,:,8) - &
                   & r_mn_basis(:,:,2)*matrix_batch_aux(:,:,9)
           case(4)
              matrix_batch_aux(:,:,7) = & ! A = H_33 R_1-H_31 R_3
                   & matrix_batch_aux(:,:,3)*r_mn_basis(:,:,1) - &
                   & matrix_batch_aux(:,:,5)*r_mn_basis(:,:,3)
              matrix_batch_aux(:,:,8) = & ! B = H_23 R_1-H_21 R_3
                   & matrix_batch_aux(:,:,6)*r_mn_basis(:,:,1) - &
                   & matrix_batch_aux(:,:,4)*r_mn_basis(:,:,3)
              matrix_batch_aux(:,:,9) = & ! C = H_13 R_1-H_11 R_3
                   & matrix_batch_aux(:,:,5)*r_mn_basis(:,:,1) - &
                   & matrix_batch_aux(:,:,1)*r_mn_basis(:,:,3)
              ! chi_12 = R_2 A - R_3 B
              matrix_batch = r_mn_basis(:,:,2)*matrix_batch_aux(:,:,7) - &
                   & r_mn_basis(:,:,3)*matrix_batch_aux(:,:,8)
           case(5)
              ! chi_22 = R_3 C - R_1 A
              matrix_batch = r_mn_basis(:,:,3)*matrix_batch_aux(:,:,9) - &
                   & r_mn_basis(:,:,1)*matrix_batch_aux(:,:,7)
           case(6)
              ! chi_32 = R_1 B - R_2 C
              matrix_batch = r_mn_basis(:,:,1)*matrix_batch_aux(:,:,8) - &
                   & r_mn_basis(:,:,2)*matrix_batch_aux(:,:,9)
           case(7)
              matrix_batch_aux(:,:,7) = & ! A = H_31 R_2-H_32 R_1
                   & matrix_batch_aux(:,:,5)*r_mn_basis(:,:,2) - &
                   & matrix_batch_aux(:,:,6)*r_mn_basis(:,:,1)
              matrix_batch_aux(:,:,8) = & ! B = H_21 R_2-H_22 R_1
                   & matrix_batch_aux(:,:,4)*r_mn_basis(:,:,2) - &
                   & matrix_batch_aux(:,:,2)*r_mn_basis(:,:,1)
              matrix_batch_aux(:,:,9) = & ! C = H_11 R_2-H_12 R_1
                   & matrix_batch_aux(:,:,1)*r_mn_basis(:,:,2) - &
                   & matrix_batch_aux(:,:,4)*r_mn_basis(:,:,1)
              ! chi_13 = R_2 A - R_3 B
              matrix_batch = r_mn_basis(:,:,2)*matrix_batch_aux(:,:,7) - &
                   & r_mn_basis(:,:,3)*matrix_batch_aux(:,:,8)
           case(8)
              ! chi_23 = R_3 C - R_1 A
              matrix_batch = r_mn_basis(:,:,3)*matrix_batch_aux(:,:,9) - &
                   & r_mn_basis(:,:,1)*matrix_batch_aux(:,:,7)
           case(9)
              ! chi_33 = R_1 B - R_2 C
              matrix_batch = r_mn_basis(:,:,1)*matrix_batch_aux(:,:,8) - &
                   & r_mn_basis(:,:,2)*matrix_batch_aux(:,:,9)
           end select
        else ! comp_dim_9 (else only the diagonal elements)
           select case(i_dim)
           case(1)
              matrix_batch_aux(:,:,7) = & ! A = H_32 R_3-H_33 R_2
                   & matrix_batch_aux(:,:,6)*r_mn_basis(:,:,3) - &
                   & matrix_batch_aux(:,:,3)*r_mn_basis(:,:,2)
              matrix_batch_aux(:,:,8) = & ! B = H_22 R_3-H_23 R_2,
                   & matrix_batch_aux(:,:,2)*r_mn_basis(:,:,3) - &
                   & matrix_batch_aux(:,:,6)*r_mn_basis(:,:,2)
              ! chi_11 = R_2 A - R_3 B
              matrix_batch = r_mn_basis(:,:,2)*matrix_batch_aux(:,:,7) - &
                   & r_mn_basis(:,:,3)*matrix_batch_aux(:,:,8)
           case(2)
              matrix_batch_aux(:,:,7) = & ! A = H_33 R_1-H_31 R_3
                   & matrix_batch_aux(:,:,3)*r_mn_basis(:,:,1) - &
                   & matrix_batch_aux(:,:,5)*r_mn_basis(:,:,3)
              matrix_batch_aux(:,:,8) = & ! C = H_13 R_1-H_11 R_3
                   & matrix_batch_aux(:,:,5)*r_mn_basis(:,:,1) - &
                   & matrix_batch_aux(:,:,1)*r_mn_basis(:,:,3)
              ! chi_22 = R_3 C - R_1 A
              matrix_batch = r_mn_basis(:,:,3)*matrix_batch_aux(:,:,8) - &
                   & r_mn_basis(:,:,1)*matrix_batch_aux(:,:,7)
           case(3)
              matrix_batch_aux(:,:,7) = & ! B = H_21 R_2-H_22 R_1
                   & matrix_batch_aux(:,:,4)*r_mn_basis(:,:,2) - &
                   & matrix_batch_aux(:,:,2)*r_mn_basis(:,:,1)
              matrix_batch_aux(:,:,8) = & ! C = H_11 R_2-H_12 R_1
                   & matrix_batch_aux(:,:,1)*r_mn_basis(:,:,2) - &
                   & matrix_batch_aux(:,:,4)*r_mn_basis(:,:,1)
              ! chi_33 = R_1 B - R_2 C
              matrix_batch = r_mn_basis(:,:,1)*matrix_batch_aux(:,:,7) - &
                   & r_mn_basis(:,:,2)*matrix_batch_aux(:,:,8)
           end select
        end if comp_dim_9
     else ! compact_indexing
        if (i_dim == 1) then
           do i_dir = 1, 9
              call mul_wrapper(matrix_wave(:,:,i_dir), &
                   & matrix_batch_aux(:,:,i_dir))
              if (do_transpose) then
                 matrix_batch_aux(:,:,i_dir) = matrix_batch_aux(:,:,i_dir) + &
                      & transpose(matrix_batch_aux(:,:,i_dir))
              end if
           end do
        end if
        no_comp_dim_9: if (n_dims == 9) then ! Full tensor
           select case(i_dim)
           case(1)
              matrix_batch_aux(:,:,10) = & ! A = H_32 R_3-H_33 R_2,
                   & matrix_batch_aux(:,:,6)*r_mn_basis(:,:,3) - &
                   & matrix_batch_aux(:,:,9)*r_mn_basis(:,:,2)
              matrix_batch_aux(:,:,11) = & ! B = H_22 R_3-H_23 R_2,
                   & matrix_batch_aux(:,:,5)*r_mn_basis(:,:,3) - &
                   & matrix_batch_aux(:,:,8)*r_mn_basis(:,:,2)
              matrix_batch_aux(:,:,12) = & ! C = H_12 R_3-H_13 R_2
                   & matrix_batch_aux(:,:,4)*r_mn_basis(:,:,3) - &
                   & matrix_batch_aux(:,:,7)*r_mn_basis(:,:,2)
              ! chi_11 = R_2 A - R_3 B
              matrix_batch = r_mn_basis(:,:,2)*matrix_batch_aux(:,:,10) - &
                   & r_mn_basis(:,:,3)*matrix_batch_aux(:,:,11)
           case(2)
              ! chi_21 = R_3 C - R_1 A
              matrix_batch = r_mn_basis(:,:,3)*matrix_batch_aux(:,:,12) - &
                   & r_mn_basis(:,:,1)*matrix_batch_aux(:,:,10)
           case(3)
              ! chi_31 = R_1 B - R_2 C
              matrix_batch = r_mn_basis(:,:,1)*matrix_batch_aux(:,:,11) - &
                   & r_mn_basis(:,:,2)*matrix_batch_aux(:,:,12)
           case(4)
              matrix_batch_aux(:,:,10) = & ! A = H_33 R_1-H_31 R_3
                   & matrix_batch_aux(:,:,9)*r_mn_basis(:,:,1) - &
                   & matrix_batch_aux(:,:,3)*r_mn_basis(:,:,3)
              matrix_batch_aux(:,:,11) = & ! B = H_23 R_1-H_21 R_3
                   & matrix_batch_aux(:,:,8)*r_mn_basis(:,:,1) - &
                   & matrix_batch_aux(:,:,2)*r_mn_basis(:,:,3)
              matrix_batch_aux(:,:,12) = & ! C = H_13 R_1-H_11 R_3
                   & matrix_batch_aux(:,:,7)*r_mn_basis(:,:,1) - &
                   & matrix_batch_aux(:,:,1)*r_mn_basis(:,:,3)
              ! chi_12 = R_2 A - R_3 B
              matrix_batch = r_mn_basis(:,:,2)*matrix_batch_aux(:,:,10) - &
                   & r_mn_basis(:,:,3)*matrix_batch_aux(:,:,11)
           case(5)
              ! chi_22 = R_3 C - R_1 A
              matrix_batch = r_mn_basis(:,:,3)*matrix_batch_aux(:,:,12) - &
                   & r_mn_basis(:,:,1)*matrix_batch_aux(:,:,10)
           case(6)
              ! chi_32 = R_1 B - R_2 C
              matrix_batch = r_mn_basis(:,:,1)*matrix_batch_aux(:,:,11) - &
                   & r_mn_basis(:,:,2)*matrix_batch_aux(:,:,12)
           case(7)
              matrix_batch_aux(:,:,10) = & ! A = H_31 R_2-H_32 R_1
                   & matrix_batch_aux(:,:,3)*r_mn_basis(:,:,2) - &
                   & matrix_batch_aux(:,:,6)*r_mn_basis(:,:,1)
              matrix_batch_aux(:,:,11) = & ! B = H_21 R_2-H_22 R_1
                   & matrix_batch_aux(:,:,2)*r_mn_basis(:,:,2) - &
                   & matrix_batch_aux(:,:,5)*r_mn_basis(:,:,1)
              matrix_batch_aux(:,:,12) = & ! C = H_11 R_2-H_12 R_1
                   & matrix_batch_aux(:,:,1)*r_mn_basis(:,:,2) - &
                   & matrix_batch_aux(:,:,4)*r_mn_basis(:,:,1)
              ! chi_13 = R_2 A - R_3 B
              matrix_batch = r_mn_basis(:,:,2)*matrix_batch_aux(:,:,10) - &
                   & r_mn_basis(:,:,3)*matrix_batch_aux(:,:,11)
           case(8)
              ! chi_23 = R_3 C - R_1 A
              matrix_batch = r_mn_basis(:,:,3)*matrix_batch_aux(:,:,12) - &
                   & r_mn_basis(:,:,1)*matrix_batch_aux(:,:,10)
           case(9)
              ! chi_33 = R_1 B - R_2 C
              matrix_batch = r_mn_basis(:,:,1)*matrix_batch_aux(:,:,11) - &
                   & r_mn_basis(:,:,2)*matrix_batch_aux(:,:,12)
           end select
        else ! no_comp_dim_9 (else diagonal only)
           select case(i_dim)
           case(1)
              matrix_batch_aux(:,:,10) = & ! A = H_32 R_3-H_33 R_2,
                   & matrix_batch_aux(:,:,6)*r_mn_basis(:,:,3) - &
                   & matrix_batch_aux(:,:,9)*r_mn_basis(:,:,2)
              matrix_batch_aux(:,:,11) = & ! B = H_22 R_3-H_23 R_2,
                   & matrix_batch_aux(:,:,5)*r_mn_basis(:,:,3) - &
                   & matrix_batch_aux(:,:,8)*r_mn_basis(:,:,2)
              ! chi_11 = R_2 A - R_3 B
              matrix_batch = r_mn_basis(:,:,2)*matrix_batch_aux(:,:,10) - &
                   & r_mn_basis(:,:,3)*matrix_batch_aux(:,:,11)
           case(2)
              matrix_batch_aux(:,:,10) = & ! A = H_33 R_1-H_31 R_3
                   & matrix_batch_aux(:,:,9)*r_mn_basis(:,:,1) - &
                   & matrix_batch_aux(:,:,3)*r_mn_basis(:,:,3)
              matrix_batch_aux(:,:,11) = & ! C = H_13 R_1-H_11 R_3
                   & matrix_batch_aux(:,:,7)*r_mn_basis(:,:,1) - &
                   & matrix_batch_aux(:,:,1)*r_mn_basis(:,:,3)
              ! chi_22 = R_3 C - R_1 A
              matrix_batch = r_mn_basis(:,:,3)*matrix_batch_aux(:,:,11) - &
                   & r_mn_basis(:,:,1)*matrix_batch_aux(:,:,10)
           case(3)
              matrix_batch_aux(:,:,10) = & ! B = H_21 R_2-H_22 R_1
                   & matrix_batch_aux(:,:,2)*r_mn_basis(:,:,2) - &
                   & matrix_batch_aux(:,:,5)*r_mn_basis(:,:,1)
              matrix_batch_aux(:,:,11) = & ! C = H_11 R_2-H_12 R_1
                   & matrix_batch_aux(:,:,1)*r_mn_basis(:,:,2) - &
                   & matrix_batch_aux(:,:,4)*r_mn_basis(:,:,1)
              ! chi_33 = R_1 B - R_2 C
              matrix_batch = r_mn_basis(:,:,1)*matrix_batch_aux(:,:,10) - &
                   & r_mn_basis(:,:,2)*matrix_batch_aux(:,:,11)
           end select
        end if no_comp_dim_9
     end if
  else ! mult_R_mn_both
     if (compact_indexing) then
        select case(i_dim)
        case(1)
           ! For the compact indexing scheme, see, e.g.,
           ! f_diamagnetic_shielding_1_diag.
           call mul_wrapper(matrix_wave(:,:,1), matrix_batch)
           call mul_wrapper(matrix_wave(:,:,2), matrix_batch_aux(:,:,1))
           matrix_batch = r_mn_basis(:,:,2)*matrix_batch_aux(:,:,1) - &
                & r_mn_basis(:,:,3)*matrix_batch
        case(2)
           call mul_wrapper(matrix_wave(:,:,3), matrix_batch)
           call mul_wrapper(matrix_wave(:,:,4), matrix_batch_aux(:,:,1))
           matrix_batch = r_mn_basis(:,:,3)*matrix_batch - &
                & r_mn_basis(:,:,1)*matrix_batch_aux(:,:,1)
        case(3)
           call mul_wrapper(matrix_wave(:,:,5), matrix_batch)
           call mul_wrapper(matrix_wave(:,:,6), matrix_batch_aux(:,:,1))
           matrix_batch = r_mn_basis(:,:,1)*matrix_batch_aux(:,:,1) - &
                & r_mn_basis(:,:,2)*matrix_batch
        end select
     else
        select case(i_dim)
        case(1, 4, 7)
           do i_dir = 1, 3
              call mul_wrapper(matrix_wave(:,:,i_dir+i_dim-1), &
                   & matrix_batch_aux(:,:,i_dir))
              if (do_transpose) matrix_batch_aux(:,:,i_dir) = &
                   & matrix_batch_aux(:,:,i_dir) + &
                   & transpose(matrix_batch_aux(:,:,i_dir))
           end do
           matrix_batch = r_mn_basis(:,:,2)*matrix_batch_aux(:,:,3) - &
                & r_mn_basis(:,:,3)*matrix_batch_aux(:,:,2)
        case(2, 5, 8)
           ! The relevant contributions were already calculated in case(1)
           matrix_batch = r_mn_basis(:,:,3)*matrix_batch_aux(:,:,1) - &
                & r_mn_basis(:,:,1)*matrix_batch_aux(:,:,3)
        case(3, 6, 9)
           matrix_batch = r_mn_basis(:,:,1)*matrix_batch_aux(:,:,2) - &
                & r_mn_basis(:,:,2)*matrix_batch_aux(:,:,1)
        end select
     end if
  end if if_mult_R_mn_both
contains
  pure subroutine mul_wrapper(matrix_wave, matrix_batch_aux)
    real(dp), intent(in) :: matrix_wave(:,:)
    real(dp), intent(out) :: matrix_batch_aux(:,:)
    call mul(wave, matrix_wave, matrix_batch_aux, transb='t', &
         & M=n_compute, N=n_compute, K=n_points, do_serial=.true.)
  end subroutine mul_wrapper
end subroutine GIAO_mult_psi
