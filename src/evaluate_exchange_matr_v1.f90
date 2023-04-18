!****s* FHI-aims/evaluate_exchange_matr_v1
!  NAME
!   evaluate_exchange_matr_v1
!  SYNOPSIS

subroutine evaluate_exchange_matr_v1 &
( KS_eigenvector, number_of_loops, occ_numbers)

  !  PURPOSE
  !  Subroutine evaluate_exchange_matrix_v1 evaluates the exchange part of
  !  the self energy. here we use a different implementation to speed up the
  !  exchange matrix evaluation when the occupied states are few.
  !
  !   Sigma^_ii(tau) = - \sum_m^occ \sum_uv O_uim V_uv O_vjm
  !                  = - \sum_m^occ \sum_u O_uim S_ujm
  !   V_uv is the bare Coulomb interaction matrix.
  !
  !  Note: this is the version 1 which sums over the occupied eigenstates when
  !  constructing the exchange matrix. No density matrix is needed.
  !
  !  JW-FIXME: We could split the occupied states into batches which are handled
  !            subsequently in order to reduce the memory footprint.
  !
  !  USES

  use dimensions
  use prodbas
  use hartree_fock
  use mpi_tasks
  use synchronize_mpi
  use runtime_choices
  use constants
  use localorb_io, only: localorb_info, use_unit, OL_norm, OL_low
  implicit none

  !  ARGUMENTS

  real*8, dimension(n_basis,n_states,n_spin,n_k_points) :: KS_eigenvector
  integer, intent(IN) :: number_of_loops
  real*8, dimension(n_states,n_spin,n_k_points) :: occ_numbers

  !  INPUTS
  !  o  number_of_loops -- integer number, the current number of the self-consistent loop
  !  o  occ_numbers -- real array,
  !       the occupation number of the electrons for each eigenstate and each spin
  !  o  KS_eigenvector -- real array,
  !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
  !  OUTPUTS
  !  none
  !  the exchange matrix (the "fock_matr" in the source code) is defined in MODULE
  !  hartree_fock
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

  integer, parameter :: max_n_real = 2**23 ! ~ 8e6 reals, 64 MB
  real*8 :: fac
  integer :: i_block, i_state, i_spin, state_off, i_glob_state
  integer :: n_prefac, n_state_block, n_state_dim, n_block
  real*8, allocatable :: O_bbbsKS(:,:,:,:)
  integer :: info, mpierr
  character*150 :: info_str
  character(*), parameter :: func = 'evaluate_exchange_matr_v1'

  integer :: i_basis_1, i_basis_2
  ! --- Broadcast KS_eigenvector

  if(use_mpi) then
     call MPI_Bcast(KS_eigenvector, n_basis*n_states*n_spin*n_k_points, &
     &              MPI_DOUBLE_PRECISION, 0, mpi_comm_global, mpierr)
  endif

  ! --- Prepare buffer (size)

  n_prefac = 8 * n_loc_prodbas * n_basis * n_spin
  n_state_dim = floor(dble(max_n_real) / dble(n_prefac))
  n_state_dim = min(n_homo_max, n_state_dim)
  n_state_dim = max(n_state_dim, 1)
  n_block = int(ceiling(dble(n_homo_max) / dble(n_state_dim)))

  write(info_str, "(A,I6,A)") &
  & '| Using bunches of', n_state_dim, ' states for eigencoefficient exchange.'
  call localorb_info(info_str, use_unit, '(2X,A)', OL_norm)
  write(info_str, "(A,I8,A)") &
  & 'Allocating about', (n_state_dim*n_prefac) / 2**20 + 1, &
  & ' MiB for eigencoefficient exchange.'
  call localorb_info(info_str, use_unit, "(2X,'| ',A)", OL_low)

  allocate(O_bbbsKS(n_loc_prodbas, n_basis, n_state_dim, n_spin), stat=info)
  call check_allocation(info, 'O_bbbsKS', func)

  ! --- Do work

  fock_matr=0.d0

  do i_block = 1, n_block   ! block of states
     state_off = (i_block-1) * n_state_dim
     n_state_block = min(n_state_dim, n_homo_max - state_off)

     call evaluate_O_bbbsKS_part(KS_eigenvector, O_bbbsKS, n_state_dim, &
     &                           n_state_block, state_off)

     do i_spin = 1, n_spin

        do i_state = 1, n_state_block
           i_glob_state = state_off + i_state

           fac = occ_numbers(i_glob_state, i_spin, 1) * dble(n_spin/2.d0)

           call dgemm('T', 'N', n_basis, n_basis, n_loc_prodbas, &
           &          fac, O_bbbsKS(:,:,i_state, i_spin), n_loc_prodbas, &
           &               O_bbbsKS(:,:,i_state, i_spin), n_loc_prodbas, &
           &          1.d0, fock_matr(:,:,i_spin), n_basis)

        end do
     end do
  end do

  ! --- Synchronization & Mixing

  ! synchronise the exchange matrix (sum over product basis functions)
  do i_spin = 1, n_spin
     call sync_matrix(fock_matr(:,:,i_spin), n_basis, n_basis)
  enddo

  if(use_density_matrix_hf .and. number_of_loops >= 0) then
     ! As the density matrix cannot be mixed in this case, mix the fock_matrix
     ! instead.  This should be OK because the dependence is linear.
     do i_spin = 1, n_spin
        call density_matrix_mixing(i_spin, number_of_loops, &
        &                          fock_matr(:,:, i_spin))
     end do
  end if
!  if(myid.eq.0) then
!       write(use_unit,*) "fock_matr:"
!       do i_basis_1 = 1, n_basis, 1
!         do i_basis_2 = 1, i_basis_1, 1
!             write(use_unit,'(2I6,f18.10)') i_basis_1, i_basis_2, &
!               fock_matr(i_basis_2, i_basis_1, i_spin)
!         enddo
!       enddo
!  endif


end subroutine evaluate_exchange_matr_v1
!---------------------------------------------------------------------
!******
