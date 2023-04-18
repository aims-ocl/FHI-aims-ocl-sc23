!****s* FHI-aims/RRS-PBC/get_rrs_pbc_atom2basis()
!  NAME
!    get_rrs_pbc_atom2basis() 
!  SYNOPSIS

    subroutine get_rrs_pbc_atom2basis()

!  PURPOSE
!  This routine table the basis index for RRS-PBC scheme
!
!  USES

      use localorb_io
      use dimensions
      use physics
      use geometry
      use numerical_utilities
      use basis
      use runtime_choices
      use mpi_tasks
      implicit none

!  ARGUMENTS

!  INPUTS
!    none
!  OUTPUT
!    none
!  COPYRIGHT
!   
!   
!   
!  HISTORY
!   
!  SOURCE
!

      ! imported variables

      ! local variables
      character*132 :: rrs_pbc_info
      character*132 :: func = 'get_rrs_pbc_atom2basis()'
      integer :: i, j, o
      integer :: info
      real*8, dimension(3)   :: tmp_vec
      logical :: rrs_pbc_flag_identify = .false.


      ! Now initialize the basis table
      do i=1, rrs_pbc_n_center_atom, 1
          rrs_pbc_center_atom(2,i) = 0
          rrs_pbc_center_atom(3,i) = 0
          do j=1, rrs_pbc_n_equal, 1
              rrs_pbc_equal_atom(2,j,i) = 0
          enddo
      enddo

      ! Now build the basis table
      do i=1, n_basis, 1
          ! Now table the basis index for atoms in the center unit cell
          do j=1, rrs_pbc_n_center_atom, 1
              if (basis_atom(i) .eq. rrs_pbc_center_atom(1,j)) then
                  if (rrs_pbc_center_atom(2,j) .eq. 0 .or. &
                      rrs_pbc_center_atom(2,j) > i) then
                      rrs_pbc_center_atom(2,j) = i
                  endif
                  rrs_pbc_center_atom(3,j) = rrs_pbc_center_atom(3,j) + 1
                  cycle
              endif
          enddo

          ! Now table the basis index for equal points
          rrs_pbc_flag_identify = .false.
          do j=1, rrs_pbc_n_center_atom, 1
              do o=1, rrs_pbc_n_equal, 1
                  if (rrs_pbc_equal_atom(1,o,j) .eq. basis_atom(i)) then
                      rrs_pbc_flag_identify = .true.
                      if (rrs_pbc_equal_atom(2,o,j) .eq. 0 .or. &
                          rrs_pbc_equal_atom(2,o,j) > i) then
                          rrs_pbc_equal_atom(2,o,j) = i
                          cycle
                      endif
                  endif
              enddo
              if (rrs_pbc_flag_identify) cycle
          enddo
      enddo

      rrs_pbc_n_center_basis = 0
      do i=1, rrs_pbc_n_center_atom, 1
          rrs_pbc_n_center_basis = &
              rrs_pbc_n_center_basis + rrs_pbc_center_atom(3,i)
      enddo
         
      write(rrs_pbc_info,'(4X,A,I12)') &
          '| Basis number in unit cell            :',rrs_pbc_n_center_basis
      call localorb_info(rrs_pbc_info)
      ! allocate hamiltonian of unit cell
      rrs_pbc_n_hamiltonian_k = &
          rrs_pbc_n_center_basis*(rrs_pbc_n_center_basis+1)/2
      allocate(rrs_pbc_hamiltonian_k(rrs_pbc_n_hamiltonian_k,n_spin),stat=info)
      call check_allocation(info, 'rrs_pbc_hamiltonian_k               ')
      allocate(rrs_pbc_overlap_k(rrs_pbc_n_hamiltonian_k),stat=info)
      call check_allocation(info, 'rrs_pbc_overlap_k                   ')
      allocate(rrs_pbc_KS_eigenvector(rrs_pbc_n_center_basis,&
          rrs_pbc_n_center_basis,n_spin),stat=info)
      call check_allocation(info, 'rrs_pbc_KS_eigenvector              ')
      allocate(rrs_pbc_band_vect(rrs_pbc_n_center_basis,&
          rrs_pbc_n_center_basis,n_spin),stat=info)
      call check_allocation(info, 'rrs_pbc_band_vect                   ')
      allocate(rrs_pbc_occ_num(rrs_pbc_n_center_basis,n_spin),stat=info)
      call check_allocation(info, 'rrs_pbc_occ_num                     ')

      end subroutine

