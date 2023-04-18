!****s* FHI-aims/prune_density_matrix_sparse_phonon_reduce_memory
!  NAME
!    prune_density_matrix_sparse_phonon_reduce_memory
!  SYNOPSIS

! module prune_density_matrix_sparse_dielectric_inspector
!   integer, dimension(:, :), allocatable :: source_place
!   logical :: inspector_done = .false.
! end module prune_density_matrix_sparse_dielectric_inspector
subroutine prune_density_matrix_sparse_dielectric(density_matrix_sparse, & 
                                       density_matrix_con, &
                                       n_compute, i_basis_index)

!  PURPOSE
!    The subroutine change the density matrix components belongs to non-zero basis functions
!    to density_matrix_con.
!
!  USES

  use dimensions
  use pbc_lists
  !use
  implicit none

!  ARGUMENTS

  real*8 :: density_matrix_sparse(n_hamiltonian_matrix_size)
  real*8 :: density_matrix_con(n_compute, n_compute)
  integer:: n_compute
  integer:: i_basis_index(n_compute)

!  INPUTS
!    o density_matrix_sparse -- total density matrix (packed matrix format)
!    o n_compute -- number of non-zero basis function in current grid batch 
!    o i_basis -- list of the non-zero basis functions in current grid batch
!
!  OUTPUT
!   o density_matrix_con -- values of density matrix balong to non-zero basis functions.
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
! SOURCE


  integer :: i_basis,j_basis,i_compute,j_compute,i_basis_uc,j_basis_uc,i_cell,j_cell
  integer :: i_coord, i_center,i_place
  integer :: k_cell,k_atom, k_cell_new,i_center_new

  !-------initial----------------------- 
  density_matrix_con = 0.0d0
  do i_compute = 1,n_compute
     i_basis    = i_basis_index(i_compute)
     i_basis_uc = Cbasis_to_basis(i_basis)
     i_cell     = center_to_cell(Cbasis_to_center(i_basis))

     do j_compute = 1,n_compute
        j_basis    = i_basis_index(j_compute)
        j_basis_uc = Cbasis_to_basis(j_basis)

        if(j_basis_uc <= i_basis_uc) then
           j_cell = center_to_cell(Cbasis_to_center(j_basis))

           do i_place = &
                index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
                index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)

              if( column_index_hamiltonian( i_place) == j_basis_uc)then

                 density_matrix_con(i_compute,j_compute) = & 
                      density_matrix_sparse(i_place) 

                 density_matrix_con(j_compute,i_compute) = & 
                      density_matrix_sparse(i_place) 

              endif

           enddo ! i_place
        endif ! j_basis_uc <= i_basis_uc

     enddo ! j_compute 
  enddo ! i_compute 


end subroutine prune_density_matrix_sparse_dielectric

! subroutine prune_density_matrix_sparse_dielectric_inspector(n_compute, i_basis_index, source_place)

! !  PURPOSE
! !    The subroutine change the density matrix components belongs to non-zero basis functions
! !    to density_matrix_con.
! !
! !  USES

!   use dimensions
!   use pbc_lists
!   !use
!   implicit none

! !  ARGUMENTS

!   !real*8 :: density_matrix_sparse(n_hamiltonian_matrix_size)
!   !real*8 :: density_matrix_con(n_compute, n_compute)
!   integer:: n_compute
!   integer:: i_basis_index(n_compute)
!   integer,pointer :: source_place(:,:)
! !  INPUTS
! !    o density_matrix_sparse -- total density matrix (packed matrix format)
! !    o n_compute -- number of non-zero basis function in current grid batch 
! !    o i_basis -- list of the non-zero basis functions in current grid batch
! !
! !  OUTPUT
! !   o density_matrix_con -- values of density matrix balong to non-zero basis functions.
! !
! !  AUTHOR
! !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
! !  SEE ALSO
! !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
! !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
! !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
! !    Computer Physics Communications (2008), submitted.
! !  COPYRIGHT
! !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
! !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
! !   the terms and conditions of the respective license agreement."
! !  HISTORY
! !    Release version, FHI-aims (2008).
! ! SOURCE


!   integer :: i_basis,j_basis,i_compute,j_compute,i_basis_uc,j_basis_uc,i_cell,j_cell
!   integer :: i_coord, i_center,i_place, i_cnt
!   integer :: k_cell,k_atom, k_cell_new,i_center_new

!   !allocate(source_place(n_compute, n_compute))
!   do i_compute = 1, n_compute
!      do j_compute = 1, n_compute
!         source_place(i_compute, j_compute) = 0
!      end do
!   end do
!   print *, "====start of a batch===="
!   do j_compute = 1,n_compute
!      j_basis    = i_basis_index(j_compute)
!      j_basis_uc = Cbasis_to_basis(j_basis)

!      do i_compute = 1,n_compute
!         i_basis    = i_basis_index(i_compute)
!         i_basis_uc = Cbasis_to_basis(i_basis)
!         i_cell     = center_to_cell(Cbasis_to_center(i_basis))


!         if(j_basis_uc <= i_basis_uc) then
!            j_cell = center_to_cell(Cbasis_to_center(j_basis))

!            !i_cnt = 0
!            do i_place = &
!                 index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
!                 index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)
              
!               if( column_index_hamiltonian( i_place) == j_basis_uc)then
!                  if (source_place(j_compute, i_compute) == 0) source_place(j_compute, i_compute) = -i_place
!                  source_place(i_compute, j_compute) = i_place
!               endif
!            enddo ! i_place
!         endif ! j_basis_uc <= i_basis_uc
!      enddo ! j_compute 
!   enddo ! i_compute 
!   print *, "=====end of a batch====="

! end subroutine prune_density_matrix_sparse_dielectric_inspector

! subroutine prune_density_matrix_sparse_dielectric_executor(density_matrix_sparse, & 
!      density_matrix_con, &
!      n_compute, i_basis_index, source_place)
!   use dimensions
!   use pbc_lists
!   !use
!   implicit none

! !  ARGUMENTS

!   real*8 :: density_matrix_sparse(n_hamiltonian_matrix_size)
!   real*8 :: density_matrix_con(n_compute, n_compute)
!   integer:: n_compute, i_compute, j_compute, i_place
!   integer:: i_basis_index(n_compute)
!   integer, pointer :: source_place(:,:)

!   do i_compute = 1, n_compute
!      do j_compute = 1, n_compute
!         i_place = abs(source_place(i_compute, j_compute))
!         if (i_place /= 0) then
!            density_matrix_con(i_compute, j_compute) = density_matrix_sparse(i_place)
!         else
!            density_matrix_con(i_compute, j_compute) = 0.0d0
!         end if
!      end do
!   end do

! end subroutine prune_density_matrix_sparse_dielectric_executor

! subroutine prune_density_matrix_sparse_dielectric_executor_reverse(density_matrix_sparse, & 
!      density_matrix_con, &
!      n_compute, i_basis_index, source_place)
!   use dimensions
!   use pbc_lists
!   !use
!   implicit none

! !  ARGUMENTS

!   real*8 :: density_matrix_sparse(n_hamiltonian_matrix_size)
!   real*8 :: density_matrix_con(n_compute, n_compute)
!   integer:: n_compute, i_compute, j_compute, i_place
!   integer:: i_basis_index(n_compute)
!   integer, pointer :: source_place(:,:)

!   do i_compute = 1, n_compute
!      do j_compute = 1, n_compute
!         i_place = source_place(i_compute, j_compute)
!         if (i_place /= -1 .and. iand(i_place, 1) == 0) then
!            i_place = ishft(i_place, -1)
!            density_matrix_sparse(i_place) = density_matrix_sparse(i_place) + density_matrix_con(i_compute, j_compute)
!         end if
!      end do
!   end do

! end subroutine prune_density_matrix_sparse_dielectric_executor_reverse

module prune_matrix_dielectric_mod
  integer, parameter:: n_thrd = 64
  type inspector_data
     integer, allocatable, dimension(:,:) :: dense_from_sparse
     integer, allocatable, dimension(:,:) :: sparse_from_dense
     integer :: sparse_split(n_thrd + 1)
     integer :: n_sparse, n_compute
     logical :: initialized
  end type inspector_data
  type(inspector_data), allocatable :: inspectors(:)
contains
  subroutine prune_matrix_dielectric_maybe_init(n_batches)
    implicit none
    integer :: n_batches
    if (.not. allocated(inspectors)) then
       allocate(inspectors(n_batches))
       inspectors(:)%initialized = .false.
    end if
  end subroutine prune_matrix_dielectric_maybe_init

  subroutine prune_matrix_dielectric_reinit(n_batches)
    implicit none
    integer :: n_batches
    if (allocated(inspectors)) then
       deallocate(inspectors)
    end if
    allocate(inspectors(n_batches))
    inspectors(:)%initialized = .false.
  end subroutine prune_matrix_dielectric_reinit

  subroutine prune_matrix_dielectric_inspector(n_compute, i_basis_index, i_batch)

    !  PURPOSE
    !    The subroutine change the density matrix components belongs to non-zero basis functions
    !    to density_matrix_con.
    !
    !  USES

    use dimensions
    use pbc_lists
    !use
    implicit none

    !  ARGUMENTS

    !real*8 :: density_matrix_sparse(n_hamiltonian_matrix_size)
    !real*8 :: density_matrix_con(n_compute, n_compute)
    integer:: n_compute
    integer:: i_basis_index(n_compute)
    integer :: i_basis,j_basis,i_compute,j_compute,i_basis_uc,j_basis_uc,i_cell,j_cell
    integer :: i_coord, i_center,i_place, i_cnt
    integer :: k_cell,k_atom, k_cell_new,i_center_new
    integer :: n_sparse, i_batch
    integer :: work_pro_thrd, i_copy, cur_split, i_place_last
    allocate(inspectors(i_batch)%dense_from_sparse(n_compute, n_compute))
    allocate(inspectors(i_batch)%sparse_from_dense(3, n_compute * n_compute))
    !allocate(source_place(n_compute, n_compute))
    inspectors(i_batch)%dense_from_sparse(:, :) = 0
    inspectors(i_batch)%sparse_from_dense(:, :) = 0
    n_sparse = 0
    !print *, "====start of a batch===="
    do j_compute = 1,n_compute
       j_basis    = i_basis_index(j_compute)
       j_basis_uc = Cbasis_to_basis(j_basis)

       do i_compute = 1,n_compute
          i_basis    = i_basis_index(i_compute)
          i_basis_uc = Cbasis_to_basis(i_basis)
          i_cell     = center_to_cell(Cbasis_to_center(i_basis))

          if(j_basis_uc <= i_basis_uc) then
             j_cell = center_to_cell(Cbasis_to_center(j_basis))

             !i_cnt = 0
             do i_place = &
                  index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
                  index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)

                if( column_index_hamiltonian( i_place) == j_basis_uc)then
                   inspectors(i_batch)%dense_from_sparse(j_compute, i_compute) = i_place
                   inspectors(i_batch)%dense_from_sparse(i_compute, j_compute) = i_place
                   n_sparse = n_sparse + 1
                   inspectors(i_batch)%sparse_from_dense(1, n_sparse) = i_place
                   inspectors(i_batch)%sparse_from_dense(2, n_sparse) = i_compute
                   inspectors(i_batch)%sparse_from_dense(3, n_sparse) = j_compute
                endif
             enddo ! i_place
          endif ! j_basis_uc <= i_basis_uc
       enddo ! j_compute 
    enddo ! i_compute 
    !print *, "=====end of a batch====="
    call heapsort_general_integer(inspectors(i_batch)%sparse_from_dense(:, :), 3, n_sparse, 1)
    
    inspectors(i_batch)%n_sparse = n_sparse
    inspectors(i_batch)%initialized = .true.
    inspectors(i_batch)%n_compute = n_compute

    work_pro_thrd = n_sparse / n_thrd
    cur_split = 1
    inspectors(i_batch)%sparse_split(1) = 1
    do i_copy = 2, n_sparse
       if (i_copy - inspectors(i_batch)%sparse_split(cur_split) >= work_pro_thrd) then
          i_place_last = inspectors(i_batch)%sparse_from_dense(1, i_copy - 1)
          i_place = inspectors(i_batch)%sparse_from_dense(1, i_copy)
          if (i_place /= i_place_last) then
             cur_split = cur_split + 1
             inspectors(i_batch)%sparse_split(cur_split) = i_copy
          end if
       end if
    end do
    inspectors(i_batch)%sparse_split(n_thrd + 1) = n_sparse + 1
    ! do i_copy = 1, n_thrd
    !    print *, inspectors(i_batch)%sparse_split(i_copy), inspectors(i_batch)%sparse_split(i_copy + 1) - 1, &
    !         inspectors(i_batch)%sparse_split(i_copy + 1) - inspectors(i_batch)%sparse_split(i_copy)
    ! end do
    ! print *, "n_sparse=", n_sparse
    ! do i_copy = 2, n_sparse
    !    i_place_last = inspectors(i_batch)%sparse_from_dense(1, i_copy - 1)
    !    i_place = inspectors(i_batch)%sparse_from_dense(1, i_copy)
    !    if (i_place < i_place_last) then
    !       print *, "not sorted sparse_from_dense, check sorting code!"
    !    end if
    ! end do
  end subroutine prune_matrix_dielectric_inspector

  subroutine prune_matrix_dielectric_dense_from_sparse(density_matrix_sparse, & 
       density_matrix_con, &
       n_compute, i_basis_index, i_batch)
    use dimensions
    use pbc_lists
    !use
    implicit none

    !  ARGUMENTS

    real*8 :: density_matrix_sparse(n_hamiltonian_matrix_size)
    real*8 :: density_matrix_con(n_compute, n_compute)
    integer:: n_compute, i_compute, j_compute, i_place
    integer:: i_basis_index(n_compute)
    integer:: i_batch
    if (.not. inspectors(i_batch)%initialized) then
       call prune_matrix_dielectric_inspector(n_compute, i_basis_index, i_batch)
    end if
    do i_compute = 1, n_compute
       do j_compute = 1, n_compute
          i_place = inspectors(i_batch)%dense_from_sparse(i_compute, j_compute)
          if (i_place /= 0) then
             density_matrix_con(i_compute, j_compute) = density_matrix_sparse(i_place)
          else
             density_matrix_con(i_compute, j_compute) = 0.0d0
          end if
       end do
    end do
  end subroutine prune_matrix_dielectric_dense_from_sparse

  subroutine prune_matrix_dielectric_sparse_from_dense(density_matrix_sparse, & 
       density_matrix_con, &
       n_compute, i_basis_index, i_batch)
    use dimensions
    use pbc_lists
    !use
    implicit none

    !  ARGUMENTS

    real*8 :: density_matrix_sparse(n_hamiltonian_matrix_size)
    real*8 :: density_matrix_con(n_compute, n_compute)
    integer:: n_compute, i_compute, j_compute, i_place
    integer:: i_basis_index(n_compute)
    integer:: i_batch, n_sparse, i_copy
    if (.not. inspectors(i_batch)%initialized) then
       call prune_matrix_dielectric_inspector(n_compute, i_basis_index, i_batch)
    end if

    n_sparse = inspectors(i_batch)%n_sparse
    do i_copy = 1, n_sparse
       i_place = inspectors(i_batch)%sparse_from_dense(1, i_copy)
       i_compute = inspectors(i_batch)%sparse_from_dense(2, i_copy)
       j_compute = inspectors(i_batch)%sparse_from_dense(3, i_copy)
       density_matrix_sparse(i_place) = density_matrix_sparse(i_place) + density_matrix_con(i_compute, j_compute)
    end do
  end subroutine prune_matrix_dielectric_sparse_from_dense
end module prune_matrix_dielectric_mod
!******
