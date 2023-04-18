!****s* FHI-aims/evaluate_O_bbbsKS_part
!  NAME
!   evaluate_O_bbbsKS_part
!  SYNOPSIS

subroutine evaluate_O_bbbsKS_part(KS_eigenvector, O_bbbsKS, n_state_dim, &
&                                 n_state_block, state_off)

  !  PURPOSE
  !
  !    This routines performs the 3-center overlap integral transformation by
  !    multiplying the ovlp_3fn with the KS eigenvectors for a given subset of
  !    KS states.  The resultant O_bbbsHF matrix elements are the integrals
  !    over (orthonormalized) auxiliary (product) basis, a regular NAO basis,
  !    and a single-partilce (KS/HF) orbital.
  !
  !  USES

  use dimensions
  use prodbas
  use hartree_fock
  use mpi_tasks
  use synchronize_mpi

  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: KS_eigenvector(n_basis,n_states,n_spin,n_k_points)
  integer, intent(IN) :: n_state_dim
  real*8, intent(OUT) :: O_bbbsKS(n_loc_prodbas, n_basis, n_state_dim, n_spin)
  integer, intent(IN) :: n_state_block
  integer, intent(IN) :: state_off

  !  INPUT
  !    o KS_eigenvector -- real array,
  !                       the eigenvector of the single-particle (KS/HF)
  !                       self-consistent calculation
  !                       (already correctly broadcasted)
  !    o n_state_dim -- dimension of O_bbbsKS
  !    o n_state_block -- number of states to contract with (<= n_state_dim)
  !    o state_off -- offset
  !                   (contract with state_off+1 ... state_off+n_state_block)
  !  OUTPUT
  !    o O_bbbsKS -- ovlp_3fn, contracted once with part of KS_eigenvector.
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

  real*8, allocatable :: o3fn_tmp(:,:)    ! n_basis, n_basis
  integer :: i_loc_prodbas, i_basis_1, i_basis_2, i_pair, i_spin
  integer :: i1, iL
  integer ::  info
  character(*), parameter :: func = 'evaluate_O_bbbsKS_part'

  allocate(o3fn_tmp(n_basis,n_basis),stat=info)
  call check_allocation(info, 'o3fn_tmp', func)

  i1 = state_off+1
  iL = state_off+n_state_block

  do i_loc_prodbas = 1, n_loc_prodbas
     i_pair = 0
     o3fn_tmp(:,:) = 0.d0
     do i_basis_1 = 1, n_basis
        do i_basis_2 = 1, i_basis_1
           i_pair= basis_nghr(i_basis_2, i_basis_1)
           if(i_pair .gt.0) then
           	  o3fn_tmp(i_basis_1,i_basis_2) = ovlp_3fn(i_pair, i_loc_prodbas)
		        o3fn_tmp(i_basis_2,i_basis_1) = ovlp_3fn(i_pair, i_loc_prodbas)
           endif
        enddo
     enddo
     do i_spin = 1, n_spin
        call dgemm('N', 'N', n_basis, n_state_block, n_basis, &
        &          1.0d0, o3fn_tmp, n_basis, &
        &                 KS_eigenvector(:, i1:iL, i_spin, 1), n_basis, &
        &          0.d0, O_bbbsKS(i_loc_prodbas,:,:,i_spin), n_basis)
     enddo
  enddo
  deallocate(o3fn_tmp)

end subroutine evaluate_O_bbbsKS_part
!---------------------------------------------------------------------
!******
