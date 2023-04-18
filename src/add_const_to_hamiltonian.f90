!****s* FHI-aims/add_const_to_hamiltonian
!  NAME
!    add_const_to_hamiltonian
!  SYNOPSIS

subroutine add_const_to_hamiltonian(hamiltonian, const)

  !  PURPOSE
  !    Add a constant to the Hamiltonian diagonal
  !  USES

  use mpi_tasks
  use dimensions
  use runtime_choices
  use pbc_lists
  implicit none

  !  ARGUMENTS

  real*8, intent(INOUT) :: hamiltonian(n_hamiltonian_matrix_size, n_spin)
  real*8, intent(IN) :: const

  !  INPUTS
  !    o hamiltonian -- old Hamiltonian
  !    o const -- constant to add
  !  OUTPUTS
  !    o hamiltonian -- old Hamiltonian + const * 1
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

  integer :: i_basis, i_offset, i_index
  integer :: i_start, i_end
  character(*), parameter :: func = 'add_const_to_hamiltonian'

  select case(packed_matrix_format)
  case(PM_none)

     do i_basis = 1, n_basis
        i_offset = (i_basis-1) * i_basis / 2
        i_index = i_offset + i_basis
        hamiltonian(i_index, :) = hamiltonian(i_index, :) + const
     end do

  case(PM_index) !---------------------------------------------------

     do i_basis = 1, n_basis
        ! Even in the periodic case, only use the first cell.
        i_start =  index_hamiltonian(1,1, i_basis)
        i_end   =  index_hamiltonian(2,1, i_basis)
        do i_index = i_start, i_end
           if (column_index_hamiltonian(i_index) == i_basis) then
              hamiltonian(i_index, :) = hamiltonian(i_index, :) + const
              exit
           end if
        end do
        if (i_index > i_end) then
           call aims_stop('Diagonal element not found', func)
        end if
     end do

  case default
     call aims_stop('Invalid case', func)
  end select

end subroutine add_const_to_hamiltonian
!******
