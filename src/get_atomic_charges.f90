!****s* FHI-aims/get_atomic_charges
!  NAME
!    get_atomic_charges
!  SYNOPSIS

subroutine get_atomic_charges(rho, delta_rho, partition_tab, &
&                             atomic_charges, atomic_spin_asymmetry)

  !  PURPOSE
  !
  !    Get atomic charges from the grid partitioning and additionally
  !    the spin asymmetry
  !
  !            \int d^3r |rho_1(r) - rho_2(r)|.
  !
  !    This definition has the advantage that it can be interpreted in a
  !    physical way.  If we have two radicals at different spatial places,
  !    both in different spin channels, we get a spin asymmetry of about two.
  !
  !  USES

  use grids
  use dimensions
  use synchronize_mpi_basic
  use mpi_tasks, only: aims_stop
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: rho(n_spin, n_full_points)
  real*8, intent(IN) :: delta_rho(n_full_points, n_spin)
  real*8, intent(IN) :: partition_tab(n_full_points)
  real*8, intent(OUT) :: atomic_charges(n_atoms, n_spin)
  real*8, intent(OUT) :: atomic_spin_asymmetry(n_atoms)

  !  INPUTS
  !    o rho -- Charge density on the integration grid
  !  OUTPUTS
  !    o atomic_charges -- Atomic charges as partitioned for integrations
  !    o atomic_spin_asymmetry -- \int d^3r p_i(r) |rho(:,1) - rho(:,2)|.
  !    o partition_tab -- Integration weights
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
  !    Release version, FHI-aims (2011).
  !  SOURCE

  integer :: i_my_batch, i_index
  integer :: i_full_point, i_atom
  real*8 :: weight, this_rho(n_spin)
  character(*), parameter :: func = 'get_atomic_charges'
  atomic_charges = 0.d0
  atomic_spin_asymmetry = 0.d0

  i_full_point = 0
  do i_my_batch = 1, n_my_batches
     do i_index = 1, batches(i_my_batch)%size
        i_full_point = i_full_point + 1
        i_atom = batches(i_my_batch)%points(i_index)%index_atom
        weight = partition_tab(i_full_point)
        if (n_spin == 1) then
           this_rho(1) = rho(1, i_full_point) + delta_rho(i_full_point, 1)
        else if (n_spin == 2) then
           this_rho(1) = rho(1, i_full_point) &
           & + 0.5*(delta_rho(i_full_point, 1) + delta_rho(i_full_point, 2))
           this_rho(2) = rho(2, i_full_point) &
           & + 0.5*(delta_rho(i_full_point, 1) - delta_rho(i_full_point, 2))
        else
           call aims_stop('Invalid n_spin', func)
        end if

        atomic_charges(i_atom, :) &
        & = atomic_charges(i_atom, :) + weight * this_rho
        if (n_spin == 2) then
           atomic_spin_asymmetry(i_atom) &
           & = atomic_spin_asymmetry(i_atom) &
           & + weight * abs(this_rho(2) - this_rho(1))
        end if
     end do
  end do

  call sync_vector(atomic_charges, n_atoms*n_spin)
  if (n_spin == 2) then
     call sync_vector(atomic_spin_asymmetry, n_atoms)
  end if
  
end subroutine get_atomic_charges
!******
