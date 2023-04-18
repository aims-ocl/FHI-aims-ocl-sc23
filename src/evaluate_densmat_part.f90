  !******
  !------------------------------------------------------------------------------
  !****s* evaluate_densmat_part
  !  NAME
  !    evaluate_densmat_part
  !  SYNOPSIS
  subroutine evaluate_densmat_part &
  ( KS_eigenvector, KS_eigenvector_complex, occ_numbers,  &
  density_matrix, density_matrix_sparse, i_spin, &
  KS_eigenvalue, energ_lower_limit, energ_upper_limit) 
    !  PURPOSE
    !    Evaluates the density matrix 
    !  USES
    use runtime_choices, only: use_symmetry_reduced_spg
    use dimensions, only: n_k_points, n_basis, n_states, n_spin, n_k_points_task,&
                          n_centers_basis_T,n_hamiltonian_matrix_size
    use sym_base, only: evaluate_densmat_sym
    use density_matrix_evaluation, only: evaluate_densmat
    implicit none

    !  ARGUMENTS

    real*8,     dimension(n_basis, n_states, n_spin,n_k_points_task):: KS_eigenvector
    complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task):: KS_eigenvector_complex
    real*8, dimension(n_states, n_spin, n_k_points), intent(in) :: occ_numbers
    ! NOTE:  For this routine to give the correct results, the occ_numbers variable passed in should have already been
    !        properly k-weighted (and are thus not the "true" occupation numbers.)
    !        In scf_solver, this is done by the kweight_occs function (and undone at the end by de_kweight_occs)
    integer :: i_spin
    real*8, dimension(n_centers_basis_T,n_centers_basis_T) :: density_matrix
    real*8 :: density_matrix_sparse(n_hamiltonian_matrix_size)
    real*8, dimension(n_states, n_spin, n_k_points) :: KS_eigenvalue
    real*8 :: energ_lower_limit, energ_upper_limit

    !  INPUTS
    !   o KS_eigenvector -- Kohn-Sham eigenvectors real format
    !   o KS_eigenvector_complex -- Kohn-Sham eigenvectors complex format
    !   o occ_numbers -- occupation of states, with k-weighting already applied
    !   o i_spin -- spin index
    !
    !  OUTPUT
    !   o density_matrix -- density matrix if non-packed matrix is in use
    !   o density_matrix_sparse -- density matrix if packed matrix is in use
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
    !

    real*8, allocatable :: zeroed_occs(:,:,:)
    integer :: i_state, i_k_point
    real*8 :: this_en

    allocate(zeroed_occs(n_states, n_spin, n_k_points))
    zeroed_occs = 0d0

    do i_k_point = 1, n_k_points
       do i_state = 1, n_states
          this_en = KS_eigenvalue(i_state, i_spin, i_k_point)
          if (this_en < energ_upper_limit .and. &
          &   this_en > energ_lower_limit) then
             zeroed_occs(i_state, i_spin, i_k_point) = &
             & occ_numbers(i_state, i_spin, i_k_point)
          end if
       end do
    end do

    if(use_symmetry_reduced_spg)then
      call evaluate_densmat_sym(KS_eigenvector, KS_eigenvector_complex, &
      & zeroed_occs, density_matrix, density_matrix_sparse, i_spin, .false.)
    else
      call evaluate_densmat(KS_eigenvector, KS_eigenvector_complex, &
      & zeroed_occs, density_matrix, density_matrix_sparse, i_spin, .false.)
    endif
    deallocate(zeroed_occs)

  end subroutine evaluate_densmat_part