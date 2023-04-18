!****s* FHI-aims/RRS-PBC/prepare_rrs_pbc_energy()
!  NAME
!    run_rrs_pbc() 
!  SYNOPSIS

    subroutine prepare_rrs_pbc_energy()

!  PURPOSE
!  This routine implement the actual RRS-PBC calculation, and can only be called
!  after parse_rrs_pbc()
!
!  USES

      use localorb_io
      use constants
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
      character*132 :: func = 'parse_rrs_pbc()'
      integer :: i, j, o, p
      integer :: info
      real*8, dimension(3)   :: tmp_k_vec
      real*8, dimension(3)   :: tmp_d_vec

      do i = 1, rrs_pbc_n_k_points, 1

          call build_rrs_pbc_k_matrices(rrs_pbc_k_point_list(i,:))

          !write(rrs_pbc_info,'(4X,A)') '| Now print overlap'
          !call localorb_info(rrs_pbc_info)
          !call output_rrs_pbc_matrix(rrs_pbc_n_center_basis,&
          !    rrs_pbc_n_center_basis, rrs_pbc_n_hamiltonian_k,&
          !    rrs_pbc_overlap_k)

          !write(rrs_pbc_info,'(4X,A)') '| Now print hamiltonian'
          !call localorb_info(rrs_pbc_info)
          !call output_rrs_pbc_matrix(rrs_pbc_n_center_basis,&
          !    rrs_pbc_n_center_basis, rrs_pbc_n_hamiltonian_k,&
          !    rrs_pbc_hamiltonian_k(:,1))

          ! Print out several useful information about this k point
          write(rrs_pbc_info,'(4X,A,I5,X,A,3(F16.8),X,A)') &
              '| K-point:',i, 'at',(rrs_pbc_k_point_list(i,j),j=1,3,1), &
              '(in units of recip. lattice)'
          call localorb_info(rrs_pbc_info)
          write(rrs_pbc_info,'(4X,A,X,3(A16),1X,A)') &
              '| K-point:','VBM','CBM','Gap','(in eV)'
          call localorb_info(rrs_pbc_info)

          call calculate_rrs_pbc_k_matrices(i)

          !write(use_unit,*) '--------------------------------------------------'
          !write(use_unit,*) rrs_pbc_KS_eigenvalue(:,1,1)
          !write(use_unit,*) '--------------------------------------------------'

          if (n_spin .eq. 1) then
              write(rrs_pbc_info,'(4X,A,1X,3(F16.8))') &
                  '|         ',&
                  rrs_pbc_KS_eigenvalue(rrs_pbc_n_electron_int(1),1,i) * hartree,&
                  rrs_pbc_KS_eigenvalue(rrs_pbc_n_electron_int(1)+1,1,i) * hartree,&
                  (rrs_pbc_KS_eigenvalue(rrs_pbc_n_electron_int(1)+1,1,i) - &
                  rrs_pbc_KS_eigenvalue(rrs_pbc_n_electron_int(1),1,i)) * hartree
              call localorb_info(rrs_pbc_info)
          else
              do o = 1, n_spin, 1
                  write(rrs_pbc_info,'(4X,A,I10,A,1X,3(F16.8))') &
                      '|  spin(',o,')',&
                      rrs_pbc_KS_eigenvalue(rrs_pbc_n_electron_int(o),o,i) * hartree,&
                      rrs_pbc_KS_eigenvalue(rrs_pbc_n_electron_int(o)+1,o,i) * hartree,&
                      (rrs_pbc_KS_eigenvalue(rrs_pbc_n_electron_int(o)+1,o,i) - &
                      rrs_pbc_KS_eigenvalue(rrs_pbc_n_electron_int(o),o,i))* hartree
                  call localorb_info(rrs_pbc_info)
              enddo
          endif

          ! store the matrix of rrs_pbc_KS_eigenvector
          call output_rrs_pbc_result(i,tmp_k_vec,1)
      enddo
      ! store the matrix of rrs_pbc_KS_eigenvalue
      call output_rrs_pbc_result(0,0,2)

      end subroutine

