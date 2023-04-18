      subroutine  compute_scgw_dipole ()

      use localorb_io
      use dimensions
      use timing
      use runtime_choices
      use physics
      use species_data
      use mixing
      use grids
      use mpi_utilities
      use pbc_lists
      use hartree_fock
      use localized_basbas
!      use plot_band
      use scaled_zora_transform
      use forces_densmat
      use synchronize_mpi
      use precondition
      use restart
      use scalapack_wrapper
      use separate_core_states
      use plus_u
      use vdw_correction
      use ll_vdwdf
      use wf_extrapolation
      use KH_core_states
      use transport
      use force_occupation
      use meta_gga_postp
      use hartree_fock_p0 ! SVL
      use octree_routines
      use vdwkernel
      use load_balancing
      use hartree_potential_storage
      use numerical_stress
      use pseudodata
      use hartree_potential_recip
      use gt
 
      implicit none
      real*8 compact_density_matrix (n_hamiltonian_matrix_size)
!      real*8 g_densmat (n_basis, n_basis)
      integer i_index, i_basis, j_basis

!      i_index=0
!      compact_density_matrix (:) = 0.d0
!      do i_basis = 1, n_basis, 1
!        do j_basis =1, i_basis, 1
!          i_index = i_index+1
!          compact_density_matrix (i_index) = green_fn(i_basis, j_basis)
!        enddo
!      enddo
      call get_real_space_density &
            ( KS_eigenvector, KS_eigenvector_complex, &
            occ_numbers, &
            partition_tab, hartree_partition_tab, rho, rho_gradient, &
            l_shell_max, delta_rho_KS, delta_rho_gradient, rho_change,&
            compact_density_matrix )

      call output_dipole_moment ()

      end subroutine compute_scgw_dipole 
