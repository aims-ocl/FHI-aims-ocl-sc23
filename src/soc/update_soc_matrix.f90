!****f* FHI-aims/update_soc_matrix
!*  NAME
!*    update_soc_matrix
!*  SYNOPSIS
subroutine update_soc_matrix(n_compute, i_basis, soc_shell, i_pauli, i_sign, &
                             ld_soc_matrix, soc_matrix )
!*  PURPOSE
!*    Updates a component of soc_matrix with the batch matrix ("shell") computed
!*    for the current batch of integration points
!*  USES
  use dimensions, only: n_periodic
  use runtime_choices, only: packed_matrix_format, PM_index, PM_none
  use localorb_io, only: localorb_info
  use mpi_tasks, only: aims_stop
  use pbc_lists, only: n_cells, index_hamiltonian, column_index_hamiltonian, &
      Cbasis_to_basis, Cbasis_to_center, center_to_cell, position_in_hamiltonian
  use load_balancing, only: use_batch_permutation, batch_perm
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  implicit none
!*  ARGUMENTS
  integer, intent(in) :: n_compute
  integer, intent(in) :: i_basis(n_compute)
  real*8, intent(in) :: soc_shell(n_compute, n_compute)
  integer, intent(in) :: i_pauli
  integer, intent(in) :: i_sign
  integer, intent(in) :: ld_soc_matrix
  real*8, intent(inout) :: soc_matrix(ld_soc_matrix, 3)
!*  INPUTS
!*    o n_compute - the size of the compute basis; the number of Cbasis
!*                  elements that could possible contribute to the current
!*                  batch
!*    o i_basis   - list that maps from compute basis indexing to Cbasis
!*                  indexing
!*    o soc_shell - the integration of the SOC operator for the current
!*                  batch
!*    o i_pauli   - Selects the real-space operation whose matrix elements
!*                  we are evalulating.  They are indexed by the Pauli matrix
!*                  they couple to, hence the name.
!*    o i_sign    - Which of the two terms in the cross product we are
!*                  evaluating in the current call.  Since SOC operator has the
!*                  form of a cross product, each real-space operation has two
!*                  terms, one with a plus sign and one with a minus sign.
!*  OUTPUTS
!*    o soc_matrix - the real-space SOC operator matrix elements between basis
!*                   elements for a given spin combination/Pauli matrix
!*  AUTHORS
!*    William Huhn (Duke University)
!*  NOTES
!*    A rewrite of update_full_matrix_p0X()
!*
!*    The implementation of second-variational SOC in FHI-aims is published in
!*      Huhn and Blum, Phys. Rev. Materials 1, 033803 (2017)
!*      https://dx.doi.org/10.1103/PhysRevMaterials.1.033803
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is subject
!*    the terms and conditions of the respective license agreement."
!*  SOURCE
  integer :: i_compute_1, i_compute_2, i_index_real
  integer :: i_basis_1, i_basis_2, i_cell, i_cell_1, i_cell_old
  integer :: i_start, i_end, i_place
  integer, dimension(n_cells) :: offset, offset_end
  integer, dimension(n_compute) :: help_basis, help_cell

  ! Variables for load balancing
  integer :: n_bp
  integer :: i,j,i_off
  integer, allocatable :: ins_idx(:)

  ! begin work

  ! If local indexing with load balancing is in effect (use_batch_permutation
  ! > 0), the local hamiltonian is always stored in full form for the local
  ! basis functions (and ONLY the local basis functions)
  if(use_batch_permutation > 0) then
    n_bp = use_batch_permutation
    call aims_allocate(ins_idx, batch_perm(n_bp)%n_basis_local, "ins_idx")

    ! Get position of basis functions of current batch within local
    ! hamiltonian
    do i=1,n_compute
      ins_idx(i) = batch_perm(n_bp)%i_basis_glb_to_loc(i_basis(i))
    enddo

    ! Insert soc_shell of current batch
    do i=1,n_compute
      i_off = (ins_idx(i)*(ins_idx(i)-1))/2
      do j=1,i ! n_compute
        soc_matrix(ins_idx(j)+i_off,i_pauli) &
             = soc_matrix(ins_idx(j)+i_off,i_pauli) &
             + i_sign*soc_shell(j,i)
      enddo
    enddo

    ! We're done, get out of here
    call aims_deallocate(ins_idx, "ins_idx")
    return
  end if

  ! Otherwise, using one of the packing schemes for indexing the real-space
  ! Hamiltonian
  select case (packed_matrix_format)

  case(PM_index) ! Also used for local indexing without load balancing
    if (n_periodic == 0) then
      i_cell = 1
      do i_compute_1 = 1, n_compute
        i_basis_1 = i_basis(i_compute_1)
        ! The range of all possible indices corresponding to i_basis_1 as the
        ! row element
        i_start = index_hamiltonian(1, 1, i_basis_1)
        i_end   = index_hamiltonian(2, 1, i_basis_1)
        do i_compute_2 = 1, n_compute
          i_basis_2 = i_basis(i_compute_2)
          if (index_hamiltonian(1, i_cell, i_basis_1) > 0) then
            ! Go over all elements in the indexing array and see if
            ! they correspond to the current set of i_compute_1 and i_compute_2
            do i_place = i_start, i_end
              ! We've found the right location in the indexing array, update the
              ! appropriate soc_matrix element and move onto the next possible
              ! i_compute_2
              if (column_index_hamiltonian(i_place) == i_basis_2) then
                if (i_compute_2.le.i_compute_1) then
                  soc_matrix(i_place,i_pauli) = soc_matrix(i_place,i_pauli) + &
                       i_sign*soc_shell(i_compute_2,i_compute_1)
                else
                  soc_matrix(i_place,i_pauli) = soc_matrix(i_place,i_pauli) + &
                       i_sign*soc_shell(i_compute_1,i_compute_2)
                end if
                i_index_real = i_place
                exit
              ! We've gone past the set of possible indexes for i_basis_2, this
              ! particular combination of i_basis_1 and i_basis_2 does not exist
              ! in the indexing array, go onto the next possible i_basis_2
              else if (column_index_hamiltonian(i_place) > i_basis_2) then
                i_index_real = i_place
                exit
              end if
            end do
            ! Start the search for the indexing element for the next i_basis_2
            ! at the current indexing element (i_basis is ordered, a larger
            ! value for i_compute_2 cannot correspond to a smaller value for
            ! i_basis_2)
            i_start = i_index_real
          end if
        end do
      end do
    else ! n_periodic != 0
      ! Arrays to convert from i_compute to basis function/cell
      do i_compute_1 = 1, n_compute
        help_basis(i_compute_1) = Cbasis_to_basis(i_basis(i_compute_1))
        help_cell(i_compute_1)  = &
             center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))
      end do

      do i_compute_1 = 1, n_compute, 1
        i_basis_1 = help_basis(i_compute_1)
        i_cell_old = help_cell(i_compute_1)
        offset_end = -1
        offset = -1

        ! Our packed matrix has at most n_basis*n_Cbasis elements, but
        ! n_compute could have as many as n_Cbasis elements, so looping over
        ! both i_compute_1 and i_compute_2 has n_Cbasis*n_Cbasis iterations.
        ! But we don't want matrix elements between two basis elements in the
        ! supercell; we want matrix elements between one basis element in the
        ! unit cell and one in the supercell.  To get around this, we map each
        ! i_compute_2 to its position *relative* to the current i_compute_1,
        ! which by symmetry shifts i_compute_1 into the unit cell.  The usual
        ! packing scheme now applies to the shifted compute elements.  This
        ! conversion is described by the offset and offset_end lists.
        do i_cell_1 = 1, n_cells
          i_cell = position_in_hamiltonian(i_cell_old, i_cell_1)
          offset(i_cell_1)     = index_hamiltonian(1, i_cell, i_basis_1)
          offset_end(i_cell_1) = index_hamiltonian(2, i_cell, i_basis_1)
        end do

        do i_compute_2 = 1, n_compute
          i_basis_2 = help_basis(i_compute_2)
          if (i_basis_2 <= i_basis_1) then
            i_cell = help_cell(i_compute_2)
            do i_place = offset(i_cell), offset_end(i_cell), 1
              if (column_index_hamiltonian(i_place) == i_basis_2) then
                if (i_compute_2 .le. i_compute_1) then
                  soc_matrix(i_place,i_pauli) = soc_matrix(i_place,i_pauli) + &
                       i_sign*soc_shell(i_compute_2,i_compute_1)
                else
                  soc_matrix(i_place,i_pauli) = soc_matrix(i_place,i_pauli) + &
                       i_sign*soc_shell(i_compute_1,i_compute_2)
                end if
                exit
              else if (column_index_hamiltonian(i_place) > i_basis_2) then
                exit
              end if
            end do
            offset(i_cell) = i_place
          end if
        end do
      end do
    end if

  !  In the case of no packing, all basis element combinations are present in
  !  the soc matrix, of which the upper triangle is stored, so the mapping is
  !  pretty obvious.
  case(PM_none)
    do i_compute_2 = 1, n_compute, 1
      do i_compute_1 = 1, i_compute_2, 1
        i_index_real = i_basis(i_compute_1) + &
             (i_basis(i_compute_2)-1)*i_basis(i_compute_2)/2
        soc_matrix(i_index_real,i_pauli) = soc_matrix(i_index_real,i_pauli) + &
             i_sign*soc_shell(i_compute_1,i_compute_2)
      end do
    end do

  case default
    call localorb_info('Invalid packing type in update_soc_matrix')
    call aims_stop
  end select
end subroutine update_soc_matrix
!******
