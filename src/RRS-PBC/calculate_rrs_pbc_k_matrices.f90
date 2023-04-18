!****s* FHI-aims/RRS-PBC/calculate_rrs_pbc_k_matrices()
!  NAME
!  calculate_rrs_pbc_k_matrices()
!  SYNOPSIS

    subroutine calculate_rrs_pbc_k_matrices(k_index)

!  PURPOSE
!  This routine solve the eigenvalue problem of given hamiltonian and overlap
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
      use lapack_wrapper
      use scalapack_wrapper
      implicit none

!  ARGUMENTS

!  INPUTS
!    k_index :: >= 1 stands for the index in rrs_pbc_k_point_list, which should be used to store
!                    the data for energy calculation.
!               = 0  for ploting band information
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
      integer :: k_index, iop
      complex*16, dimension(:,:), allocatable       :: ovlp_work_full
      complex*16, dimension(:,:), allocatable       :: hamiltonian_work
      real*8,dimension(:), allocatable              :: tmp_KS_eigenvalue
      real*8,dimension(:,:), allocatable              :: tmp_KS_eigenvector
      integer :: i,j,o,p,k
      integer :: i_index
      integer :: i_spin
      character*130 :: rrs_pbc_info



      if (.not.allocated(ovlp_work_full)) then
          allocate(ovlp_work_full(rrs_pbc_n_center_basis,rrs_pbc_n_center_basis))
      endif
      if (.not.allocated(hamiltonian_work)) then
          allocate(hamiltonian_work(rrs_pbc_n_center_basis,rrs_pbc_n_center_basis))
      endif
      if (.not.allocated(tmp_KS_eigenvalue)) then
          allocate(tmp_KS_eigenvalue(rrs_pbc_n_center_basis))
      endif

      if (.not.allocated(tmp_KS_eigenvector)) then
          allocate(tmp_KS_eigenvector(rrs_pbc_n_center_basis,rrs_pbc_n_center_basis))
      endif

      do i_spin = 1, n_spin, 1
          do i = 1, rrs_pbc_n_center_basis, 1
              do j = 1, i, 1
                  i_index = j + i * (i - 1) / 2
                  ovlp_work_full(i,j) = rrs_pbc_overlap_k(i_index)
                  hamiltonian_work(i,j) = rrs_pbc_hamiltonian_k(i_index,i_spin)
                  ovlp_work_full(j,i) = conjg(ovlp_work_full(i,j))
                  hamiltonian_work(j,i) = conjg(hamiltonian_work(i,j))
              enddo
          enddo


          !write(rrs_pbc_info,'(4X,A)') '| Now print overlap'
          !call localorb_info(rrs_pbc_info)
          !call output_rrs_pbc_matrix_2D(rrs_pbc_n_center_basis,&
          !    ovlp_work_full)

          !write(rrs_pbc_info,'(4X,A)') '| Now print hamiltonian'
          !call localorb_info(rrs_pbc_info)
          !call output_rrs_pbc_matrix_2D(rrs_pbc_n_center_basis,&
          !    hamiltonian_work)

          call diagnalize_rrs_pbc &
              (rrs_pbc_n_center_basis,hamiltonian_work, ovlp_work_full, &
               tmp_KS_eigenvalue,1,k)
          !call complex_lapack_solver_fast(rrs_pbc_n_center_basis,rrs_pbc_n_center_basis, &
          !    ovlp_work_full,hamiltonian_work,tmp_KS_eigenvalue,tmp_KS_eigenvector)
          if (n_spin .gt. 1) then
              write(rrs_pbc_info,'(4X,A,I5)') '| In spin state :', i_spin
              call localorb_info(rrs_pbc_info)
          endif
          if (iop > 0) then
              if (k .eq. 0) then
                  rrs_pbc_KS_eigenvalue(:,i_spin,k_index) = tmp_KS_eigenvalue(:)
                  rrs_pbc_KS_eigenvector(:,:,i_spin) = hamiltonian_work(:,:)
              else
                  rrs_pbc_KS_eigenvalue(:,i_spin,k_index) = 0.0d0
                  rrs_pbc_KS_eigenvector(:,:,i_spin) = 0.0d0
              endif
          else if (iop .eq. 0) then
              if (k .eq. 0) then
                  rrs_pbc_band_info(:,i_spin)   = tmp_KS_eigenvalue(:)
                  rrs_pbc_band_vect(:,:,i_spin) = hamiltonian_work(:,:)
              else
                  rrs_pbc_band_info(:,i_spin)   = 0.0d0
                  rrs_pbc_band_vect(:,:,i_spin) = 0.0d0
              endif
          endif
      enddo

      deallocate( ovlp_work_full)
      deallocate(hamiltonian_work)
      deallocate(tmp_KS_eigenvalue)

      end subroutine

