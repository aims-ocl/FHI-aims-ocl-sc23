!****s* FHI-aims/look_at_rho
!  NAME
!    look_at_rho
!  SYNOPSIS

subroutine look_at_rho(rho, delta_rho, partition_tab)

  !  PURPOSE
  !
  !     Call get_atomic_charges and do corresponding output.
  !
  !  USES
  use dimensions, only: n_spin, n_full_points, n_atoms
  use runtime_choices, only: out_zero_multipoles
  use localorb_io, only: use_unit, OL_norm, localorb_info
  use physics, only: n_electrons
  use species_data, only: species_z
  use geometry, only: species
  use mpi_tasks, only: check_allocation
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: rho(n_spin, n_full_points)
  real*8, intent(IN) :: delta_rho(n_spin, n_full_points)
  real*8, intent(IN) :: partition_tab(n_full_points)

  !  INPUTS
  !    o rho -- Old charge density on the integration grid
  !    o delta_rho -- Charge density update
  !    o partition_tab -- Integration weights
  !  OUTPUTS
  !    none [will write out]
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

  real*8, allocatable :: atomic_charges(:,:), atomic_spin_asymmetry(:)
  integer :: i_atom
  integer :: info
  character*150 :: info_str
  character(*), parameter :: func = 'look_at_rho'

  allocate(atomic_charges(n_atoms, n_spin), stat=info)
  call check_allocation(info, 'atomic_charges', func)
  allocate(atomic_spin_asymmetry(n_atoms), stat=info)
  call check_allocation(info, 'atomic_spin_asymmetry', func)

  call get_atomic_charges(rho, delta_rho, partition_tab, &
  &                       atomic_charges, atomic_spin_asymmetry)

  if (out_zero_multipoles) then
     call localorb_info('', use_unit, '(A)', OL_norm)
  end if

  write(info_str, &
  &     "(2X,'Integration grid: ',A,' = ',ES14.6)") &
  & 'deviation in total charge (<rho> - N_e)', &
  & sum(atomic_charges) - n_electrons
  call localorb_info(info_str, use_unit, '(A)', OL_norm)
  if (n_spin == 2) then
     write(info_str, &
     &     "(2X,'| ',A,' =',F12.4,'; ',A,' =',F12.4)") &
     & '<rho_up-rho_dn> = N_up-N_dn', &
     & sum(atomic_charges(:, 1) - atomic_charges(:, 2)), &
     & '<|rho_up-rho_dn|>', &
     & sum(atomic_spin_asymmetry)
     call localorb_info(info_str, use_unit, '(A)', OL_norm)
  end if

  if (out_zero_multipoles) then
     write(info_str, "(2X,'| ',A)") &
     & 'Partitioning above numbers into atomic contributions:'
     call localorb_info(info_str, use_unit, '(A)', OL_norm)

     if (n_spin == 1) then
        write(info_str, "(2X,'| ',A6,'  ',A12)") &
        & 'i_atom', 'part. charge'
     else
        write(info_str, "(2X,'| ',A6,'  ',A12,'  ',A12,'  ',A13)") &
        & 'i_atom', 'part. charge', 'part. magn.', 'part. <|...|>'
     end if
     call localorb_info(info_str, use_unit, '(A)', OL_norm)
     do i_atom = 1, n_atoms
        if (n_spin == 1) then
           write(info_str, &
           &     "(2X,'| ',I4,'  ',F12.6)") &
           & i_atom, &
           & sum(atomic_charges(i_atom, :)) - species_z(species(i_atom))
        else
           write(info_str, &
           &     "(2X,'| ',I4,'  ',F12.6,'  ',F12.6,'  ',F13.6)") &
           & i_atom, &
           & sum(atomic_charges(i_atom, :)) - species_z(species(i_atom)), &
           & (atomic_charges(i_atom, 1) - atomic_charges(i_atom, 2)), &
           & atomic_spin_asymmetry(i_atom)
        end if
        call localorb_info(info_str, use_unit, '(A)', OL_norm)
     end do
     call localorb_info('', use_unit, '(A)', OL_norm)
  end if

  deallocate(atomic_charges, atomic_spin_asymmetry)

end subroutine look_at_rho
!******
