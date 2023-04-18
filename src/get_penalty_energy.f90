!****s* FHI-aims/get_penalty_energy
!  NAME
!    get_penalty_energy
!  SYNOPSIS

subroutine get_penalty_energy(KS_eigenvector, KS_eigenvector_complex, &
&                             occ_numbers, penalty, penalty_energy)

  !  PURPOSE
  !    Calculate penalty energy.
  !  USES

  use mpi_tasks
  use synchronize_mpi_basic
  use dimensions
  use runtime_choices
  use pbc_lists
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: KS_eigenvector(n_basis, n_states, &
  &                                    n_spin, n_k_points_task)
  complex*16, intent(IN) :: KS_eigenvector_complex(n_basis, n_states, &
  &                                                n_spin, n_k_points_task)
  real*8, intent(IN) :: occ_numbers(n_states, n_spin, n_k_points)
  real*8, intent(IN) :: penalty
  real*8, intent(OUT) :: penalty_energy

  !  INPUTS
  !    o KS_eigenvector{_complex} -- Eigenvectors
  !    o penalty -- penalty
  !  OUTPUTS
  !    o penalty_energy -- penalty_energy
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

  integer :: i_state, i_spin
  integer :: i_k_point, i_k
  real*8 :: occ, normsq
  character(*), parameter :: func = 'get_penalty_energy'

  call check_occs(func, occ_numbers, .false.)

  penalty_energy = 0.d0
  if (penalty <= 0.d0) return
  i_k = 0
  do i_k_point = 1, n_k_points
     if (myid == modulo(i_k_point, n_tasks) .and. myid <= n_k_points) then
        i_k = i_k + 1
        do i_state = 1, n_states
           do i_spin = 1, n_spin
              occ = occ_numbers(i_state, i_spin, i_k_point) &
              &     * k_weights(i_k_point)
              if (real_eigenvectors) then
                 normsq = sum(KS_eigenvector(:, i_state, i_spin, i_k)**2)
              else
                 normsq = sum(abs(KS_eigenvector_complex(:, i_state, i_spin, i_k)**2))
              end if
              penalty_energy = penalty_energy + occ * penalty * normsq
           end do
        end do
     end if
  end do
  call sync_real_number(penalty_energy)

end subroutine get_penalty_energy
!******
