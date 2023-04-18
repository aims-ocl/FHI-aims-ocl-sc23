!****s* FHI-aims/evaluate_exchange_matr_LVL_eigen
!  NAME
!    evaluate_exchange_matr_LVL_eigen
!  SYNOPSIS

subroutine evaluate_exchange_matr_LVL_eigen(KS_eigenvector, &
&                                           loc_number_of_loops, occ_numbers)

  !  PURPOSE
  !
  !    Calculate the exchange matrix using LVL, possibly making use of
  !    sparsity.  This is an edit of evaluate_exchange_matr_LVL_eigen() with
  !    an added first-order correction as specified by Dunlap.  This is for
  !    testing purposes, only.
  !
  !  USES

  use dimensions
  use prodbas
  use hartree_fock
  use runtime_choices
  use mpi_tasks
  use synchronize_mpi
  use constants
  use scalapack_utils
  use timing
  use localorb_io, only: localorb_info, OL_norm, OL_low, use_unit
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: KS_eigenvector(n_basis,n_states,n_spin,n_k_points)
  integer, intent(IN) :: loc_number_of_loops
  real*8, intent(IN) :: occ_numbers(n_states,n_spin,n_k_points)

  !  INPUTS
  !    o KS_eigenvector -- eigencoefficients
  !    o loc_number_of_loops -- current number of the self-consistent loop
  !                             use negative value for no mixing at all.
  !    o occ_numbers -- occupation numbers
  !  OUTPUTS
  !    none (fock_matr in hartree_fock.f90 is set)
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  integer, parameter :: dlen_ = 9
  integer :: n_prodbas, n_state_dim
  integer, allocatable :: basbas2col(:), basbas2row(:)
  integer, parameter :: max_n_real = 2**22 ! ~ 4e6 reals, 32 MB
  integer :: i_spin, i_block
  integer :: n_block
  integer :: i_loc_prodbas
  integer :: state_off, n_state_block, i_state, i_glob_state
  integer :: i_val
  real*8, allocatable :: KS_vec(:,:)
  real*8, allocatable :: coeff_fnKSbb(:,:,:), coeff_fnKScl(:,:,:)
  integer :: info
  integer :: n_other, ldcoeff
  integer :: desc_fnKSbb(dlen_), desc_bbbb(dlen_)
  real*8, allocatable :: glob_coul(:,:), glob_fnKScl(:,:,:), glob_fnKSbb(:,:,:)
  real*8 :: time_fnKSbb(4), time_fnKScl(4), time_final(4), time_sync(4)
  character*150 :: info_str
  logical :: correct_to_2nd
  integer :: i_basis_1, i_basis_2, i_pair
  real*8, allocatable :: o3fn_tmp(:,:)
  real*8 :: sym_err
  character(*), parameter :: func = 'evaluate_exchange_matr_LVL_eigen'

  call localorb_info('Calculating exact exchange matrix in RI-LVL', &
  &                  use_unit, '(2X,A)', OL_norm)

  time_fnKSbb = 0.d0; time_fnKScl = 0.d0; time_final = 0.d0; time_sync = 0.d0

  if (n_k_points > 1) call aims_stop('Multiple k-points not implemented', func)

  correct_to_2nd = (RI_type == RI_LVL_2nd)

  if (use_scalapack) then
     n_prodbas = n_max_loc_prodbas
  else
     n_prodbas = n_basbas    ! Account for global arrays
     if (n_tasks > 2 .and. n_basis*n_basbas > 2**20) then
        ! Large calculation (> 10^6 numbers) on truely parallel machine.
        call localorb_info('  ** Use scalapack to reduce memory footprint.')
     end if
  end if
  n_state_dim = floor(dble(max_n_real) / dble(n_basis*n_prodbas))
  n_state_dim = min(n_states, n_state_dim)
  n_state_dim = max(n_state_dim, 1)
  n_block = int(ceiling(dble(n_homo_max) / dble(n_state_dim)))

  write(info_str, "(A,I6,A)") &
  & '| Using bunches of', n_state_dim, ' states for LVL exchange.'
  call localorb_info(info_str, use_unit, '(2X,A)', OL_norm)
  write(info_str, "('Allocating about',I8,' MiB for LVL exchange.')") &
  & (8 * 2*n_basis*n_state_dim*(n_max_loc_prodbas+1)) / 2**20 + 1
  call localorb_info(info_str, use_unit, "(2X,'| ',A)", OL_low)
  call localorb_info('', use_unit, "(A)", OL_low)
  allocate(KS_vec(n_basis, n_state_dim), stat=info)
  call check_allocation(info, 'KS_vec', func)
  allocate(coeff_fnKSbb(n_basis, n_state_dim, n_loc_prodbas), stat=info)
  call check_allocation(info, 'coeff_fnKSbb', func)
  allocate(coeff_fnKScl(n_basis, n_state_dim, n_loc_prodbas),stat=info)
  call check_allocation(info, 'coeff_fnKScl', func)
  allocate(basbas2row(n_basbas), basbas2col(n_basbas), stat=info)
  call check_allocation(info, 'basbas2xxx', func)
  call get_basbas_to_rowcol(.false., basbas2row, basbas2col)

  if (correct_to_2nd) then
     allocate(o3fn_tmp(n_basis, n_basis), stat=info)
     call check_allocation(info, 'o3fn_tmp', func)  
  end if

  fock_matr = 0.d0
  do i_spin = 1, n_spin
     do i_block = 1, n_block   ! block of states
        state_off = (i_block-1) * n_state_dim
        n_state_block = min(n_state_dim, n_homo_max - state_off)

        call start_timer(time_fnKSbb)

        ! --- multiply coeff3fn (& transpose) to KS_eigenvector -> coeff2fn1KS

        do i_state = 1, n_state_block
           i_glob_state = state_off + i_state
           KS_vec(:, i_state) = KS_eigenvector(:, i_glob_state, i_spin, 1) &
           &                    * sqrt(occ_numbers(i_glob_state, i_spin, 1) &
           &                           * dble(n_spin) / 2.d0)
        end do

        call get_fnKSbb(n_state_block, n_state_dim, basbas2col, &
        &               coeff_3fn_ten, KS_vec, coeff_fnKSbb)

        call stop_timer(time_fnKSbb)
        call start_timer(time_fnKScl)

        ! --- O2fn1KS * V

        n_other = n_state_block * n_basis
        if (use_scalapack) then
           LDcoeff = n_state_dim * n_basis
           call descinit(desc_fnKSbb, n_other, n_basbas, mb_aux, nb_aux, &
           &             0, 0, my_blacs_ctxt_aux, LDcoeff, info)
           if (info /= 0) call aims_stop('descinit(desc_fnKSbb) failed', func)
           call descinit(desc_bbbb, n_basbas, n_basbas, mb_aux, nb_aux, &
           &             0, 0, my_blacs_ctxt_aux, n_basbas, info)
           if (info /= 0) call aims_stop('descinit(desc_bbbb) failed', func)

           call pdgemm('N', 'N', n_other, n_basbas, n_basbas, &
           &           1.d0, coeff_fnKSbb, 1, 1, desc_fnKSbb, &
           &                 coulomb_matr_lvl, 1, 1, desc_bbbb, &
           &           0.d0, coeff_fnKScl, 1, 1, desc_fnKSbb)
           ! ! This is slightly slower...
           ! call pdsymm('R', 'U', n_other, n_basbas, &
           ! &           1.d0, coulomb_matr_lvl, 1, 1, desc_bbbb, &
           ! &                 coeff_fnKSbb, 1, 1, desc_fnKSbb, &
           ! &           0.d0, coeff_fnKScl, 1, 1, desc_fnKSbb)

        else
           write(info_str, "(A,I8,' MB for non-scalapack LVL exchange.')") &
           & 'Additionally allocating about', &
           &  (8 * 2*n_basis*n_state_dim*(n_basbas+1)) / 2**20 + 1
           call localorb_info(info_str, use_unit, "(2X,'| ',A)", OL_low)

           allocate(glob_coul(n_basbas, n_basbas), stat=info)
           call check_allocation(info, 'glob_coul', func)
           allocate(glob_fnKSbb(n_basis, n_state_block, n_basbas), stat=info)
           call check_allocation(info, 'glob_fnKSbb', func)
           allocate(glob_fnKScl(n_basis, n_state_block, n_basbas), stat=info)
           call check_allocation(info, 'glob_fnKScl', func)

           call gather_auxmat(glob_coul, coulomb_matr_lvl, n_basbas)
           call gather_auxmat(glob_fnKSbb, coeff_fnKSbb(:,1:n_state_block,:),&
           &                  n_other)
           if (myid == 0) then
             call dgemm('N', 'N', n_other, n_basbas, n_basbas, &
             &          1.d0, glob_fnKSbb, n_other, &
             &                glob_coul, n_basbas, &
             &          0.d0, glob_fnKScl, n_other)
           end if
           call scatter_auxmat(glob_fnKScl, coeff_fnKScl(:,1:n_state_block,:),&
           &                   n_other)
           deallocate(glob_coul)
           deallocate(glob_fnKSbb)
           deallocate(glob_fnKScl)
        end if

        ! --- 2nd-order correction

        if (correct_to_2nd) then
           ! coeff_fnKScl := 2*(direct integrals) - (ordinary product)
           do i_loc_prodbas = 1, n_loc_prodbas
              i_pair = 0
              o3fn_tmp(:,:) = 0.d0
              do i_basis_1 = 1, n_basis
                 do i_basis_2 = 1, i_basis_1
                    i_pair= basis_nghr(i_basis_2, i_basis_1)
                    if(i_pair .gt.0) then
                       o3fn_tmp(i_basis_1,i_basis_2) &
                       & = ovlp_3fn(i_pair, i_loc_prodbas)
                       o3fn_tmp(i_basis_2,i_basis_1) &
                       & = ovlp_3fn(i_pair, i_loc_prodbas)
                    endif
                 enddo
              enddo
              call dgemm('N', 'N', n_basis, n_state_block, n_basis, &
              &          2.0d0, o3fn_tmp, n_basis, &
              &                 KS_vec, n_basis, &
              &         -1.d0, coeff_fnKScl(:,:,i_loc_prodbas), n_basis)
           end do
        end if

        call stop_timer(time_fnKScl)
        call start_timer(time_final)

        ! --- Contract to fock_matr

        ! Contract to fock_matrix over second and third index
        if (n_state_block == n_state_dim) then
           n_other = n_state_block * n_loc_prodbas
           call dgemm('N', 'T', n_basis, n_basis, n_other, &
           &          1.d0, coeff_fnKSbb, n_basis, &
           &                coeff_fnKScl, n_basis, &
           &          1.d0, fock_matr(:,:, i_spin), n_basis)
        else
           do i_loc_prodbas = 1, n_loc_prodbas
             call dgemm('N', 'T', n_basis, n_basis, n_state_block, &
             &          1.d0, coeff_fnKSbb(:,:, i_loc_prodbas), n_basis, &
             &                coeff_fnKScl(:,:, i_loc_prodbas), n_basis, &
             &          1.d0, fock_matr(:,:, i_spin), n_basis)
           end do
        end if
        call stop_timer(time_final)
     end do
  end do

  call start_timer(time_sync)
  call sync_vector(fock_matr, n_basis**2 * n_spin)
  call stop_timer(time_sync)

  ! Symmetrize exchange matrix if correction to 2nd order is used.
  if (correct_to_2nd) then
     sym_err = 0.d0
     do i_spin = 1, n_spin
        sym_err = max(maxval(abs(fock_matr(:,:, i_spin) - &
        &                        transpose(fock_matr(:,:, i_spin)))), sym_err)
        fock_matr(:,:, i_spin) = &
        & 0.5d0 * (fock_matr(:,:, i_spin) + &
        &          transpose(fock_matr(:,:, i_spin)))
     end do
     write(info_str, "('Maximum symmetry error in exchange matrix:',ES10.2)") &
     & sym_err
     call localorb_info(info_str, use_unit, "(2X,A)", OL_norm)
  end if

  if(use_density_matrix_hf .and. loc_number_of_loops >= 0) then
     ! As the density matrix cannot be mixed in this case, mix the fock_matrix
     ! instead.  This works because the dependence is linear.
     do i_spin = 1, n_spin
        call density_matrix_mixing(i_spin, loc_number_of_loops, &
        &                          fock_matr(:,:, i_spin))
     end do
  end if

  call output_timeheader('2X', 'Fock timings:', OL_low)
  call output_timer('Time for 2fnbb -> fnKSbb', time_fnKSbb(3:4), '2X', OL_low)
  call output_timer('Time for fnKSbb -> fnKScl', time_fnKScl(3:4), '2X', OL_low)
  call output_timer('Time for fnKScl -> fock', time_final(3:4), '2X', OL_low)
  call output_timer('Time for fock sync', time_sync(3:4), '2X', OL_low)
  call localorb_info('', use_unit, '(A)', OL_low)

  deallocate(basbas2row, basbas2col)
  deallocate(KS_vec, coeff_fnKSbb, coeff_fnKScl)

end subroutine evaluate_exchange_matr_LVL_eigen
!******
