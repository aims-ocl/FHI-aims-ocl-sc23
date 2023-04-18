!****h* FHI-aims/synchronize_mpi
!  NAME
!    synchronize_mpi - provides the synchronization utilities
!                      for the parallel (MPI) environment
!  SYNOPSIS
      module synchronize_mpi
!  PURPOSE
!    This module performs the high-level synchronization of MPI-tasks
!  USES
      use synchronize_mpi_basic
      use mpi_tasks
      implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
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
!
!******

      contains
!----------------------------------------------------------------------
!****s* synchronize_mpi/sync_initialize_integrals_grids
!  NAME
!    sync_initialize_integrals_grids
!  SYNOPSIS
      subroutine sync_initialize_integrals_grids(i_atom, &
           n_radial, &
           n_lebedev, n_angular, &
           r_angular, w_angular, division_boundaries, &
           n_division)
!  PURPOSE
!    Synchronize grid data in the routine initialize_integrals_p0
!  USES
      use mpi_utilities, only: radial_task_list
      use dimensions, only: n_max_angular, n_max_angular, n_max_angular_division
      implicit none
!  ARGUMENTS
      integer :: i_atom, n_radial
      integer :: n_lebedev(n_radial)
      integer :: n_angular(n_radial)
      real*8 :: r_angular(3, n_max_angular, n_radial)
      real*8 :: w_angular(n_max_angular, n_radial)
      integer :: division_boundaries(n_max_angular_division+1, &
           n_radial)
      integer :: n_division(n_radial)
!  INPUTS
!    o i_atom -- the atom that the grid data is related to
!    o n_radial -- number of radial shells in the atomistic grid
!    o n_lebedev -- index of the lebedev grid at each radial point
!    o n_angular -- number of angular points in the radial grid
!    o r_angular -- locations of the grid points in each angular shell
!    o w_angular -- integration weight of each point in atomistic grid
!    o division_boundaries -- indeces of the boundaries for the subdivision of the
!                             angular shells
!    o n_division -- number of subdivisions of each angular shell
!  OUTPUT
!    all arrays are synchorinized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables

      integer :: myid_for_bc
      integer :: i_radial
      integer :: mpierr

      if (.not.use_mpi) return

      do i_radial = 1, n_radial, 1

         myid_for_bc = radial_task_list(i_radial, i_atom)

         call MPI_Bcast(n_lebedev(i_radial), 1, &
              MPI_INTEGER, &
              myid_for_bc, mpi_comm_global, mpierr)

         call MPI_Bcast(n_angular(i_radial), 1, &
              MPI_INTEGER, &
              myid_for_bc, mpi_comm_global, mpierr)

         call MPI_Bcast(r_angular(:,:, i_radial), &
              3*n_max_angular, &
              MPI_DOUBLE_PRECISION, &
              myid_for_bc, mpi_comm_global, mpierr)

         call MPI_Bcast(w_angular(:,i_radial), &
              n_max_angular, &
              MPI_DOUBLE_PRECISION, &
              myid_for_bc, mpi_comm_global, mpierr)

         call MPI_Bcast(division_boundaries(:,i_radial), &
              n_max_angular_division+1, &
              MPI_INTEGER, &
              myid_for_bc, mpi_comm_global, mpierr)

         call MPI_Bcast(n_division(i_radial), &
              1, MPI_INTEGER, &
              myid_for_bc, mpi_comm_global, mpierr)

!     end sync over i_radial
      enddo

      end subroutine sync_initialize_integrals_grids
!******
!----------------------------------------------------------------------
!****s* synchronize_mpi/sync_initialize_integrals_grids_p0
!  NAME
!    sync_initialize_integrals_grids_p0
!  SYNOPSIS
      subroutine sync_initialize_integrals_grids_p0(i_atom, &
           n_radial, &
           n_lebedev )
!  PURPOSE
!    Synchronize grid data in the routine initialize_integrals_p0
!  USES
      use mpi_utilities, only: radial_task_list
      implicit none
!  ARGUMENTS
      integer :: i_atom, n_radial
      integer :: n_lebedev(n_radial)
!  INPUTS
!    o i_atom -- the atom that the grid data is related to
!    o n_radial -- number of radial shells in the atomistic grid
!    o n_lebedev -- index of the lebedev grid at each radial point
!  OUTPUT
!    the array n_lebedev is synchorinized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables

      integer :: myid_for_bc
      integer :: i_radial
      integer :: mpierr

      if (.not.use_mpi) return

      do i_radial = 1, n_radial, 1

         myid_for_bc = radial_task_list(i_radial, i_atom)

         call MPI_Bcast(n_lebedev(i_radial), 1, &
              MPI_INTEGER, &
              myid_for_bc, mpi_comm_global, mpierr)

!     end sync over i_radial
      enddo

      end subroutine sync_initialize_integrals_grids_p0
!******
!----------------------------------------------------------------------
!****s* synchronize_mpi/sync_initialize_integrals_matrices
!  NAME
!    sync_initialize_integrals_matrices
!  SYNOPSIS
      subroutine sync_initialize_integrals_matrices( &
           hamiltonian, overlap_matrix)
!  PURPOSE
!    Synchronize matrices for the routine initialize_integrals_p0
!  USES
      use dimensions, only: n_hamiltonian_matrix_size, n_spin
      implicit none
!  ARGUMENTS
      real*8 :: hamiltonian(n_hamiltonian_matrix_size,n_spin)
      real*8 :: overlap_matrix(n_hamiltonian_matrix_size)
!  INPUTS
!    o hamiltonian -- the Hamilton matrix
!    o overlap_matrix -- the overlap matrix
!  OUTPUT
!    the arrays are synchorinized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variable

!$$$      real*8, dimension(:,:), allocatable :: aux_matrix_mpi

      call sync_vector(hamiltonian,n_hamiltonian_matrix_size*n_spin)
      call sync_vector(overlap_matrix,n_hamiltonian_matrix_size)

!$$$      if (.not.use_mpi) return
!$$$
!$$$      if (.not.allocated(aux_matrix_mpi)) then
!$$$         allocate(aux_matrix_mpi(n_hamiltonian_matrix_size,n_spin))
!$$$      end if
!$$$
!$$$      aux_matrix_mpi = 0.0d0
!$$$      call MPI_ALLREDUCE(hamiltonian, aux_matrix_mpi,
!$$$     +     n_hamiltonian_matrix_size*n_spin, MPI_DOUBLE_PRECISION,
!$$$     +     MPI_SUM, mpi_comm_global, mpierr)
!$$$      hamiltonian = aux_matrix_mpi
!$$$
!$$$      aux_matrix_mpi(:,1) = 0.0d0
!$$$      call MPI_ALLREDUCE(overlap_matrix, aux_matrix_mpi(:,1),
!$$$     +     n_hamiltonian_matrix_size, MPI_DOUBLE_PRECISION,
!$$$     +     MPI_SUM, mpi_comm_global, mpierr)
!$$$      overlap_matrix = aux_matrix_mpi(:,1)
!$$$
!$$$      deallocate( aux_matrix_mpi )

      end subroutine sync_initialize_integrals_matrices
!******
!-------------------------------------------------------------
!****s* synchronize_mpi/sync_eigenvalues
!  NAME
!    sync_eigenvalues
!  SYNOPSIS
      subroutine sync_eigenvalues( &
         eigen_value &
           )
!  PURPOSE
!    Synchronize the eigenvalue data
!  USES
      use dimensions, only: n_states, n_spin, n_k_points
      implicit none
!  ARGUMENTS
      real*8 :: eigen_value(n_states,n_spin,n_k_points)
!  INPUTS
!    o eigen_value -- the Kohn-Sham eigenvalues
!  OUTPUT
!    the array eigen_value is synchorinized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variable

      call sync_vector(eigen_value,n_states*n_spin*n_k_points)

!$$$      real*8 :: aux_matrix_mpi(n_states,n_spin,n_k_points)
!$$$
!$$$      if (.not.use_mpi) return
!$$$
!$$$      aux_matrix_mpi = 0.0d0
!$$$      call MPI_ALLREDUCE(  eigen_value, aux_matrix_mpi,
!$$$     +     n_states*n_spin*n_k_points , MPI_DOUBLE_PRECISION,
!$$$     +     MPI_SUM, mpi_comm_global, mpierr)
!$$$      eigen_value  = aux_matrix_mpi

      end subroutine sync_eigenvalues
!******
!-------------------------------------------------------------
!****s* synchronize_mpi/sync_density
!  NAME
!    sync_density
!  SYNOPSIS
      subroutine sync_density(rho_change)
!  PURPOSE
!    Synchronize change in the electron density
!  USES
      use dimensions, only: n_spin
      implicit none
!  ARGUMENTS
      real*8, dimension(n_spin) :: rho_change
!  INPUTS
!    o rho_change -- the change in the density for both spin channels
!  OUTPUT
!    the array rho_change is synchorinized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
! SOURCE


!     local variables
      real*8 :: rho_change_mpi
      integer :: i_spin
      integer :: mpierr

      if (.not.use_mpi) return

      do i_spin = 1, n_spin, 1
        rho_change_mpi = 0.0d0
        call MPI_ALLREDUCE(rho_change(i_spin), rho_change_mpi, &
             1, MPI_DOUBLE_PRECISION, &
             MPI_SUM, mpi_comm_global, mpierr)
        rho_change(i_spin) = rho_change_mpi
      enddo

      end subroutine sync_density

!      subroutine sync_density_2(delta_rho_gradient)

!      real*8 :: delta_rho_gradient(3, n_int_points,
!     +     n_spin)

!      integer :: myid_for_bc
!      real*8, dimension(:,:,:), allocatable  :: delta_rho_gradient_mpi

!     counters
!      integer :: i_atom
!      integer :: i_radial

!      if (.not.use_mpi) return

!      if (.not.allocated(delta_rho_gradient_mpi)) then
!         allocate(delta_rho_gradient_mpi(3, n_int_points, n_spin))
!      end if

!     broadcast the result to all threads

!      delta_rho_gradient_mpi = 0.0d0
!
!      call MPI_ALLREDUCE(delta_rho_gradient,
!     +     delta_rho_gradient_mpi, 3*n_int_points*n_spin,
!     +     MPI_DOUBLE_PRECISION, MPI_SUM,
!     +     mpi_comm_global, mpierr)
!      delta_rho_gradient = delta_rho_gradient_mpi
!
!      deallocate(delta_rho_gradient_mpi)
!
!      end subroutine sync_density_2

!******
!----------------------------------------------------------
!****s* synchronize_mpi/get_core_states_density_matrix
!  NAME
!    get_core_states_density_matrix
!  SYNOPSIS
      subroutine get_core_states_density_matrix(i_atom, density_matrix, &
           size)
!  PURPOSE
!    Synchronize the density matrix for the separate core states
!  USES
      use mpi_utilities, only: task_list
      implicit none
!  ARGUMENTS
      integer:: i_atom, size
      real*8 :: density_matrix(size)
!  INPUTS
!    o i_atom -- index of the current atom
!    o size -- size of the array
!    o density_matrix -- the density matrix
!  OUTPUT
!    the array density_matrix is synchorinized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      integer :: mpierr
      if (.not.use_mpi) return

      call MPI_Bcast( density_matrix, &
           size, &
           MPI_DOUBLE_PRECISION, &
           task_list(i_atom), mpi_comm_global, mpierr)


      end subroutine get_core_states_density_matrix
!******
!----------------------------------------------------------
!****s* synchronize_mpi/sync_sum_up_whole_potential
!  NAME
!    sync_sum_up_whole_potential
!  SYNOPSIS
      subroutine sync_sum_up_whole_potential( &
           hartree_delta_energy, &
           hartree_multipole_correction, &
           hartree_multipole_error, &
           en_density_embed, &
           en_elec_delta )
!  PURPOSE
!    Synchronize the quantities in the summation of the entire potential
!  ARGUMENTS
      real*8 :: hartree_delta_energy
      real*8 :: hartree_multipole_correction
      real*8 :: hartree_multipole_error
      real*8 :: en_density_embed
      real*8 :: en_elec_delta
!  INPUTS
!    o hartree_delta_energy --
!    o hartree_multipole_correction --
!    o hartree_multipole_error --
!    o en_density_embed --
!    o en_elec_delta --
!  OUTPUT
!    the quantities are synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables
      real*8 :: hartree_delta_energy_mpi
      real*8 :: hartree_multipole_correction_mpi
      real*8 :: hartree_multipole_error_mpi
      real*8 :: en_density_embed_mpi
      real*8 :: en_elec_delta_mpi

      integer :: mpierr


!     broadcast the result to all threads

      if (.not.use_mpi) return

!      do i_atom = 1, n_atoms, 1
!
!         myid_for_bc = task_list(i_atom)
!
!         do i_radial = 1, n_radial(species(i_atom)), 1
!
!            call MPI_Bcast(potential(:,i_radial,i_atom),
!     +           n_max_angular,
!     +           MPI_DOUBLE_PRECISION,
!     +           myid_for_bc, mpi_comm_global, mpierr)
!         enddo
!      enddo

      hartree_delta_energy_mpi = 0.0d0
      en_elec_delta_mpi = 0.d0
      call MPI_ALLREDUCE(hartree_delta_energy, &
           hartree_delta_energy_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)

      call MPI_ALLREDUCE(en_elec_delta, &
           en_elec_delta_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)

      hartree_delta_energy = hartree_delta_energy_mpi
      en_elec_delta = en_elec_delta_mpi

      hartree_multipole_correction_mpi = 0.0d0
      call MPI_ALLREDUCE(hartree_multipole_correction, &
           hartree_multipole_correction_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      hartree_multipole_correction = &
           hartree_multipole_correction_mpi

      hartree_multipole_error_mpi = 0.0d0
      call MPI_ALLREDUCE(hartree_multipole_error, &
           hartree_multipole_error_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      hartree_multipole_error = &
           hartree_multipole_error_mpi

      en_density_embed_mpi = 0.0d0
      call MPI_ALLREDUCE(en_density_embed, &
           en_density_embed_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      en_density_embed = &
           en_density_embed_mpi

      end subroutine sync_sum_up_whole_potential
!******
!----------------------------------------------------------
!****s* synchronize_mpi/sync_average_potential
!  NAME
!    sync_average_potential
!  SYNOPSIS
      subroutine sync_average_potential( &
           average_delta_v_hartree_real )
!  PURPOSE
!    Synchronizes only one quantity, should be added into sync_sum_up_whole_potential if permanently useful.
!  ARGUMENTS
      real*8 :: average_delta_v_hartree_real
!  INPUTS
!    o average_delta_v_hartree_real - a real number
!  OUTPUT
!    the quantities are synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables
      real*8 :: aux_real_mpi

      integer :: mpierr


!     broadcast the result to all threads

      if (.not.use_mpi) return

      aux_real_mpi = 0.0d0
      call MPI_ALLREDUCE(average_delta_v_hartree_real, &
           aux_real_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)

      average_delta_v_hartree_real = aux_real_mpi

    end subroutine sync_average_potential
!******
!----------------------------------------------------------
!****s* synchronize_mpi/sync_integrate_hamiltonian
!  NAME
!    sync_integrate_hamiltonian
!  SYNOPSIS
      subroutine sync_integrate_hamiltonian( &
           hamiltonian )
!  PURPOSE
!    Synchronize the Hamilton matrix after integration.
!  USES
      use dimensions, only: n_hamiltonian_matrix_size, n_spin
      implicit none
!  ARGUMENTS
      real*8 :: hamiltonian(n_hamiltonian_matrix_size,n_spin)
!  INPUTS
!    o hamiltonian -- the Hamilton matrix
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables
!      real*8, dimension(:,:), allocatable :: hamiltonian_mpi

      call sync_vector(hamiltonian,n_hamiltonian_matrix_size*n_spin)

!$$$      if (.not.use_mpi) return
!$$$
!$$$      if (.not.allocated(hamiltonian_mpi)) then
!$$$         allocate(hamiltonian_mpi(n_hamiltonian_matrix_size,n_spin))
!$$$      end if
!$$$
!$$$C     here we need the extra array hamiltonian_mpi since in MPI the
!$$$C     send and recieve buffers are not allowed to be the same
!$$$C     the current solution is to compute everything in hamiltonian
!$$$C     and collect the results to hamiltonian_mpi and then copy back
!$$$C     another way to achieve this is to pick a root node that doesn't
!$$$C     compute its share of the hamiltonian but waits and in the end makes
!$$$C     one reduce followed by a broadcast
!$$$
!$$$      hamiltonian_mpi = 0.0d0
!$$$      call MPI_ALLREDUCE(hamiltonian, hamiltonian_mpi,
!$$$     +     n_hamiltonian_matrix_size*n_spin, MPI_DOUBLE_PRECISION,
!$$$     +     MPI_SUM, mpi_comm_global, mpierr)
!$$$      hamiltonian = hamiltonian_mpi
!$$$
!$$$      deallocate(hamiltonian_mpi)

      end subroutine sync_integrate_hamiltonian

!******
!------------------------------------------------
!****s* synchronize_mpi/sync_density_matrix_sparse
!  NAME
!    sync_density_matrix_sparse
!  SYNOPSIS
      subroutine sync_density_matrix_sparse( &
           hamiltonian )
!  PURPOSE
!    Synchronize a sparse density matrix.
!  USES
      use dimensions, only: n_hamiltonian_matrix_size
      implicit none
!  ARGUMENTS
      real*8 :: hamiltonian(n_hamiltonian_matrix_size)
!  INPUTS
!    o hamiltonian -- the array to be synchronized
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(hamiltonian,n_hamiltonian_matrix_size)

      end subroutine sync_density_matrix_sparse
!******
!------------------------------------------------
!****s* synchronize_mpi/sync_sparse_matrix
!  NAME
!    sync_sparse_matrix
!  SYNOPSIS
      subroutine sync_sparse_matrix( &
           sparse_matrix )
!  PURPOSE
!    Synchronize a sparse matrix.
!  USES
      use dimensions, only: n_hamiltonian_matrix_size
      implicit none
!  ARGUMENTS
      real*8 :: sparse_matrix(n_hamiltonian_matrix_size)
!  INPUTS
!    o sparse_matrix -- the array to be synchronized
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(sparse_matrix,n_hamiltonian_matrix_size)

      end subroutine sync_sparse_matrix
!******
!----------------------------------------------------
!****s* synchronize_mpi/sync_density_matrix
!  NAME
!    sync_density_matrix
!  SYNOPSIS
      subroutine sync_density_matrix( &
           density_matrix )
!  PURPOSE
!    Synchronize the density matrix.
!  USES
      use dimensions, only: n_centers_basis_T
      implicit none
!  ARGUMENTS
      real*8 :: density_matrix(n_centers_basis_T,n_centers_basis_T )
!  INPUTS
!    o density_matrix -- the array to be synchronized
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(density_matrix, &
           n_centers_basis_T*n_centers_basis_T)

      end subroutine sync_density_matrix

!******
!----------------------------------------------------------
!****s* synchronize_mpi/sync_integrate_ovlp
!  NAME
!    sync_integrate_ovlp
!  SYNOPSIS
      subroutine sync_integrate_ovlp( &
           overlap_matrix )
!  PURPOSE
!    Synchronize the overlap matrix after inegration.
!  USES
      use dimensions, only: n_hamiltonian_matrix_size
      implicit none
!  ARGUMENTS
      real*8 :: overlap_matrix(n_hamiltonian_matrix_size)
!  INPUTS
!    o overlap_matrix -- the overlap matrix to be synchronized
!  OUTPUT
!    the matrix is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(overlap_matrix,n_hamiltonian_matrix_size)

      end subroutine sync_integrate_ovlp
!******
!--------------------------------------------------------------
!****s* synchronize_mpi/sync_get_free_superpos_energy
!  NAME
!    sync_get_free_superpos_energy
!  SYNOPSIS
      subroutine sync_get_free_superpos_energy( &
           hartree_energy_free, en_xc, en_pot_xc, en_ion_ion, &
           en_elec_free, &
           en_ion_embed, en_density_embed, en_vdw, en_ll_vdw, &
           en_ll_vdw_err, en_lda_c, en_pbe_c, ionic_forces &
           )
!  PURPOSE
!    Synchronize the quantities in get_free_superpos_energy_p1
!  USES
     use dimensions, only: use_forces, use_vdw_correction, use_ll_vdwdf, &
         n_atoms, use_embedding_potential, use_embedding_pp
     implicit none
!  ARGUMENTS
      real*8 :: hartree_energy_free
      real*8 :: en_xc
      real*8 :: en_pot_xc
      real*8 :: en_ion_ion
      real*8 :: en_ion_embed
      real*8 :: en_density_embed
      real*8 :: en_vdw
      real*8 :: en_ll_vdw
      real*8 :: en_ll_vdw_err
      real*8 :: en_lda_c
      real*8 :: en_pbe_c
      real*8 :: en_elec_free
      real*8, dimension(3, n_atoms) :: ionic_forces
!  INPUTS
!    o hartree_energy_free --
!    o en_xc --
!    o en_pot_xc --
!    o en_ion_ion --
!    o en_density_embed --
!    o en_vdw --
!    o en_ll_vdw --
!    o en_ll_vdw_err --
!    o en_lda_c --
!    o en_pbe_c --
!    o en_elec_free --
!    o ionic_forces --
!  OUTPUT
!    the quantities are synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
! SOURCE


!     local variables
      real*8 :: hartree_energy_free_mpi
      real*8 :: en_xc_mpi
      real*8 :: en_pot_xc_mpi
      real*8 :: en_ion_ion_mpi
      real*8 :: en_ion_embed_mpi
      real*8 :: en_density_embed_mpi
      real*8 :: en_vdw_mpi
      real*8 :: en_ll_vdw_mpi
      real*8 :: en_ll_vdw_err_mpi
      real*8 :: en_lda_c_mpi
      real*8 :: en_pbe_c_mpi
      real*8 :: en_elec_free_mpi
      real*8, dimension(3, n_atoms) :: forces_mpi
      integer :: mpierr

!     broadcast the result to all threads

      if (.not.use_mpi) return

      hartree_energy_free_mpi = 0.0d0
      call MPI_ALLREDUCE(hartree_energy_free, &
           hartree_energy_free_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      hartree_energy_free = hartree_energy_free_mpi

      en_xc_mpi = 0.0d0
      call MPI_ALLREDUCE(en_xc, en_xc_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      en_xc = en_xc_mpi

      en_pot_xc_mpi = 0.0d0
      call MPI_ALLREDUCE(en_pot_xc, en_pot_xc_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      en_pot_xc = en_pot_xc_mpi

      en_ion_ion_mpi = 0.0d0
      call MPI_ALLREDUCE(en_ion_ion, en_ion_ion_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      en_ion_ion = en_ion_ion_mpi

      en_elec_free_mpi = 0.0d0
      call MPI_ALLREDUCE(en_elec_free, en_elec_free_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      en_elec_free = en_elec_free_mpi

      if (use_embedding_potential.or.use_embedding_pp) then
        en_ion_embed_mpi = 0.0d0
        call MPI_ALLREDUCE(en_ion_embed, en_ion_embed_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        en_ion_embed = en_ion_embed_mpi

        en_density_embed_mpi = 0.0d0
        call MPI_ALLREDUCE(en_density_embed, en_density_embed_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        en_density_embed = en_density_embed_mpi
      end if

! vdw correction
      if (use_vdw_correction) then
        en_vdw_mpi = 0.0d0
        call MPI_ALLREDUCE(en_vdw, en_vdw_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        en_vdw = en_vdw_mpi
      end if
      if (use_ll_vdwdf) then
!
        en_ll_vdw_mpi = 0.0d0
        call MPI_ALLREDUCE(en_ll_vdw, en_ll_vdw_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        en_ll_vdw = en_ll_vdw_mpi
!
        en_ll_vdw_err_mpi = 0.0d0
        call MPI_ALLREDUCE(en_ll_vdw_err, en_ll_vdw_err_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        en_ll_vdw_err = en_ll_vdw_err_mpi
!
        en_lda_c = 0.0d0
        call MPI_ALLREDUCE(en_lda_c, en_lda_c_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        en_lda_c = en_lda_c_mpi
!
        en_pbe_c = 0.0d0
        call MPI_ALLREDUCE(en_pbe_c, en_pbe_c_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        en_pbe_c = en_pbe_c_mpi
!
      end if

      if (use_forces) then
        forces_mpi = 0.0d0
        call MPI_ALLREDUCE(ionic_forces, &
           forces_mpi, 3*n_atoms, &
           MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        ionic_forces = forces_mpi
      end if

      end subroutine sync_get_free_superpos_energy
!******
!--------------------------------------------------------------
!****s* synchronize_mpi/sync_vdw_correction
!  NAME
!    sync_vdw_correction
!  SYNOPSIS
      subroutine sync_vdw_correction(en_vdw,vdw_forces)
!  PURPOSE
!    Synchronize the quantities for vdW correction.
!  USES
     use dimensions, only: use_forces, use_vdw_correction, n_atoms, &
         use_vdw_correction_hirshfeld
     implicit none
!  ARGUMENTS
      real*8 :: en_vdw
      real*8, dimension(3, n_atoms) :: vdw_forces
!  INPUTS
!    o en_vdw --
!    o vdw_forces --
!  OUTPUT
!    the quantities are synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables
      real*8 :: en_vdw_mpi
      real*8, dimension(3, n_atoms) :: forces_mpi
      integer :: mpierr

      if (.not.use_mpi) return

! vdw correction
      if (use_vdw_correction_hirshfeld) then
        en_vdw_mpi = 0.0d0
        call MPI_ALLREDUCE(en_vdw, en_vdw_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        en_vdw = en_vdw_mpi
      end if

      if (use_forces) then
        forces_mpi = 0.0d0
        call MPI_ALLREDUCE(vdw_forces, &
           forces_mpi, 3*n_atoms, &
           MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        vdw_forces = forces_mpi
      end if

      end subroutine sync_vdw_correction

      subroutine sync_mbd_finite_diff(mbd_ene_fd_components)
!  PURPOSE
!    Synchronize finite difference energy components in MBD@rsSCS.
!  USES
      use dimensions, only: n_atoms
      implicit none

      real*8, dimension(0:6*n_atoms) ::mbd_ene_fd_components
      real*8, dimension(0:6*n_atoms) ::mbd_ene_fd_components_mpi
      integer :: mpierr

      if (.not.use_mpi) return

        mbd_ene_fd_components_mpi = 0.0d0
        call MPI_ALLREDUCE(mbd_ene_fd_components, mbd_ene_fd_components_mpi, &
           1+6*n_atoms, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
       mbd_ene_fd_components = mbd_ene_fd_components_mpi


      end subroutine sync_mbd_finite_diff


!******
!--------------------------------------------------------------
!****s* synchronize_mpi/sync_ll_vdwdf
!  NAME
!    sync_ll_vdwdf
!  SYNOPSIS
      subroutine sync_ll_vdwdf &
                 (en_ll_vdw,en_ll_vdw_err,ll_vdw_forces,en_lda_c,en_pbe_c)
!  PURPOSE
!    Synchronize the quantities for LL vdWdf correction.
!  USES
     use dimensions, only: use_forces, use_ll_vdwdf, n_atoms
     implicit none
!  ARGUMENTS
      real*8 :: en_ll_vdw, en_ll_vdw_err
      real*8 :: en_lda_c, en_pbe_c
      real*8, dimension(3, n_atoms) :: ll_vdw_forces
!  INPUTS
!    o en_ll_vdw --
!    o en_ll_vdw_err --
!    o ll_vdw_forces --
!    o en_lda_c --
!    o en_pbe_c --
!  OUTPUT
!    the quantities are synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables
      real*8 :: en_ll_vdw_mpi, en_ll_vdw_err_mpi
      real*8 :: en_lda_c_mpi, en_pbe_c_mpi
      real*8, dimension(3, n_atoms) :: forces_mpi
      integer :: mpierr

      if (.not.use_mpi) return

      if (use_ll_vdwdf) then
!
        en_ll_vdw_mpi = 0.0d0
        call MPI_ALLREDUCE(en_ll_vdw, en_ll_vdw_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        en_ll_vdw = en_ll_vdw_mpi
        en_ll_vdw_err_mpi = 0.0d0
        call MPI_ALLREDUCE(en_ll_vdw_err, en_ll_vdw_err_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        en_ll_vdw_err = en_ll_vdw_err_mpi
        en_lda_c_mpi = 0.0d0
        call MPI_ALLREDUCE(en_lda_c, en_lda_c_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        en_lda_c = en_lda_c_mpi
        en_pbe_c_mpi = 0.0d0
        call MPI_ALLREDUCE(en_pbe_c, en_pbe_c_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        en_pbe_c = en_pbe_c_mpi
!
      end if

      if (use_forces) then
        forces_mpi = 0.0d0
        call MPI_ALLREDUCE(ll_vdw_forces, &
           forces_mpi, 3*n_atoms, &
           MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        ll_vdw_forces = forces_mpi
      end if

      end subroutine sync_ll_vdwdf
!******
!--------------------------------------------------------------
!****s* synchronize_mpi/sync_ext_charge_forces
!  NAME
!    sync_ext_charge_forces
!  SYNOPSIS
      subroutine sync_ext_charge_forces(ext_charge_forces)
!  PURPOSE
!    Synchronize forces on external charges (QM/MM).
!  USES
     use dimensions, only: use_forces, n_multipoles
     implicit none
!  ARGUMENTS
      real*8, dimension(3, n_multipoles) :: ext_charge_forces
!  INPUTS
!    o en_charge_forces --
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables
      real*8, dimension(3, n_multipoles) :: forces_mpi
      integer :: mpierr

      if (.not.use_mpi) return

      if (use_forces) then
        forces_mpi = 0.0d0
        call MPI_ALLREDUCE(ext_charge_forces, &
           forces_mpi, 3*n_multipoles, &
           MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        ext_charge_forces = forces_mpi
      end if

      end subroutine sync_ext_charge_forces
!******
!--------------------------------------------------------------
!****s* synchronize_mpi/sync_pp_charge_forces
!  NAME
!    sync_ext_charge_forces
!  SYNOPSIS
      subroutine sync_pp_charge_forces(pp_forces)
!  PURPOSE
!    Synchronize forces on external charges (QM/MM).
!  USES
     use dimensions, only: use_forces, n_pp_atoms
     implicit none
!  ARGUMENTS
      real*8, dimension(3, n_pp_atoms) :: pp_forces
!  INPUTS
!    o en_charge_forces --
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables
      real*8, dimension(3, n_pp_atoms) :: forces_mpi
      integer :: mpierr

      if (.not.use_mpi) return

      if (use_forces) then
        forces_mpi = 0.0d0
        call MPI_ALLREDUCE(pp_forces, &
           forces_mpi, 3*n_pp_atoms, &
           MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        pp_forces = forces_mpi
      end if

      end subroutine sync_pp_charge_forces
!******
!--------------------------------------------------------------
!****s* synchronize_mpi/sync_nlcc_forces
!  NAME
!    sync_nlcc_forces
!  SYNOPSIS
      subroutine sync_nlcc_forces(pp_forces)
!  PURPOSE
!    Synchronize forces on external charges (QM/MM).
!  USES
     use dimensions, only: use_forces, n_atoms
     implicit none
!  ARGUMENTS
      real*8, dimension(3, n_atoms) :: pp_forces
!  INPUTS
!    o en_charge_forces --
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables
      real*8, dimension(3, n_atoms) :: forces_mpi
      integer :: mpierr

      if (.not.use_mpi) return

      if (use_forces) then
        forces_mpi = 0.0d0
        call MPI_ALLREDUCE(pp_forces, &
           forces_mpi, 3*n_atoms, &
           MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        pp_forces = forces_mpi
      end if

      end subroutine sync_nlcc_forces
!******

!--------------------------------------------------------------
!****s* synchronize_mpi/sync_update_xc_potential
!  NAME
!    sync_update_xc_potential
!  SYNOPSIS
      subroutine sync_update_xc_potential( &
           en_xc, en_pot_xc &
           )
!  PURPOSE
!    Synchronize energies in the xc-potential routine
!  ARGUMENTS
      real*8 :: en_xc
      real*8 :: en_pot_xc
!  INPUTS
!    o en_xc --
!    o en_pot_xc --
!  OUTPUT
!    the quantities are synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables
      real*8 :: en_xc_mpi
      real*8 :: en_pot_xc_mpi
      integer :: mpierr

!     broadcast the result to all threads

      if (.not.use_mpi) return

      en_xc_mpi = 0.0d0
      call MPI_ALLREDUCE(en_xc, en_xc_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      en_xc = en_xc_mpi

      en_pot_xc_mpi = 0.0d0
      call MPI_ALLREDUCE(en_pot_xc, en_pot_xc_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      en_pot_xc = en_pot_xc_mpi

      end subroutine sync_update_xc_potential
!******
!--------------------------------------------------------------
!****s* synchronize_mpi/sync_evaluate_grid_quantities_rho_grad
!  NAME
!    sync_evaluate_grid_quantities_rho_grad
!  SYNOPSIS
      subroutine sync_evaluate_grid_quantities_rho_grad( &
           rho_gradient &
           )
!  PURPOSE
!    Synchronize density gradient over the grid
!  USES
      use dimensions, only: n_max_angular, n_max_radial, n_atoms
      use mpi_utilities, only: task_list
      implicit none
!  ARGUMENTS
      real*8 :: rho_gradient(3, n_max_angular, n_max_radial, n_atoms)
!  INPUTS
!    o rho_gradient -- density gradient
!  OUTPUT
!    the density gradient is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variable
      integer :: myid_for_bc
      integer :: mpierr

!     counter
      integer :: i_atom

      if (.not.use_mpi) return

      do i_atom = 1, n_atoms, 1

         myid_for_bc =  task_list(i_atom)

         call MPI_Bcast(rho_gradient(:,:,:,i_atom), &
              3*n_max_angular*n_max_radial, &
              MPI_DOUBLE_PRECISION, &
              myid_for_bc, mpi_comm_global, mpierr)

      enddo

      end subroutine sync_evaluate_grid_quantities_rho_grad
!******
!--------------------------------------------------------------
!****s* synchronize_mpi/sync_workload
!  NAME
!    sync_workload
!  SYNOPSIS
      subroutine sync_workload( &
           n_full_points, n_int_points, &
           n_full_points_total, n_int_points_total &
           )
!  PURPOSE
!    Synchronize the workload indicators.
!  ARGUMENTS
      integer :: n_full_points
      integer :: n_int_points
      integer :: n_full_points_total
      integer :: n_int_points_total
!  INPUTS
!    o n_full_points -- number of grid points in a task
!    o n_int_points -- number of active integration points in a task
!  OUTPUT
!    o n_full_points_total -- total number of grid points
!    o n_int_points_total -- total number of active integration points
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     locals
      integer :: mpierr

!     broadcast the result to all threads

      if (.not.use_mpi) then
         n_full_points_total = n_full_points
         n_int_points_total = n_int_points
         return
      end if

      call MPI_ALLREDUCE(n_full_points, &
           n_full_points_total, &
           1, MPI_INTEGER, &
           MPI_SUM, mpi_comm_global, mpierr)

      call MPI_ALLREDUCE(n_int_points, &
           n_int_points_total, &
           1, MPI_INTEGER, &
           MPI_SUM, mpi_comm_global, mpierr)

      end subroutine sync_workload
!******
!--------------------------------------------------------------
!****s* synchronize_mpi/sync_n_avg_compute
!  NAME
!    sync_n_avg_compute
!  SYNOPSIS
      subroutine sync_n_avg_compute( &
           n_avg_compute_dens )
!  PURPOSE
!    Synchronize the number of non-zero basis functions summed over all batches 
!  ARGUMENTS
      real*8 :: n_avg_compute_dens
!  INPUTS
!    o n_avg_compute_dens -- the number of non-zero basis functions over batches in
!                            a task
!  OUTPUT
!    the number of non-zero basis functions is summed up over all batches
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local
      real*8 :: n_avg_compute_dens_mpi
      integer :: mpierr

      if (.not.use_mpi) return

      n_avg_compute_dens_mpi = 0.0d0
      call MPI_ALLREDUCE(n_avg_compute_dens, &
           n_avg_compute_dens_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      n_avg_compute_dens = n_avg_compute_dens_mpi

      end subroutine sync_n_avg_compute
!******
!--------------------------------------------------------------
!****s* synchronize_mpi/sync_pulay_matrix
!  NAME
!    sync_pulay_matrix
!  SYNOPSIS
      subroutine sync_pulay_matrix( pulay_matrix )
!  PURPOSE
!    Synchronize the Pulay matrix
!  USES
      use dimensions, only: n_max_pulay
      implicit none
!  ARGUMENTS
      real*8 :: pulay_matrix(n_max_pulay, n_max_pulay)
!  INPUTS
!    o pulay_matrix -- Pulay matrix
!  OUTPUT
!    the Pulay matrix is synchronized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      call sync_vector(pulay_matrix, n_max_pulay**2)

      end subroutine sync_pulay_matrix
!******
!--------------------------------------------------------------
!****s* synchronize_mpi/sync_broyden_matrix
!  NAME
!    sync_broyden_matrix
!  SYNOPSIS
      subroutine sync_broyden_matrix( broyden_matrix )
!  PURPOSE
!    Synchronize the Broyden matrix
!  USES
      use dimensions, only: n_max_broyden
      implicit none
!  ARGUMENTS
      real*8 :: broyden_matrix(n_max_broyden, n_max_broyden)
!  INPUTS
!    o broyden_matrix -- Broyden matrix
!  OUTPUT
!    the Broyden matrix is synchronized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      call sync_vector(broyden_matrix(1,:), n_max_broyden)
      call sync_vector(broyden_matrix(2:n_max_broyden,1), n_max_broyden-1)

      end subroutine sync_broyden_matrix
!******
!--------------------------------------------------------------
!****s* synchronize_mpi/sync_pulay_matrix_first_row_and_vector
!  NAME
!    sync_pulay_matrix_first_row_and_vector
!  SYNOPSIS
      subroutine sync_pulay_matrix_first_row_and_vector( pulay_matrix, &
           pulay_vector, pulay_saved_iter )
!  PURPOSE
!    Synchronize the first row if the Pulay matrix and the Pulay vector
!  USES
      use dimensions, only: n_max_pulay
      implicit none
!  ARGUMENTS
      integer :: pulay_saved_iter
      real*8 :: pulay_matrix(n_max_pulay, n_max_pulay)
      real*8 :: pulay_vector(pulay_saved_iter)
!  INPUTS
!    o pulay_saved_iter -- number of saved iterations so far
!    o pulay_matrix -- Pulay matrix
!    o pulay_vector -- Pulay vector
!  OUTPUT
!    the first row of the Pulay matrix and the Pulay vector 
!    are synchronized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


! VB: No spin allowed here because this subroutine is called from pulay_mix.f,
!     which does not know about spin

!     local variables
      real*8, dimension(:), allocatable :: pulay_mpi
      integer :: mpierr

!      integer :: i_spin

      if (.not.use_mpi) return

      if (.not.allocated(pulay_mpi)) then
         allocate(pulay_mpi(pulay_saved_iter),stat=mpierr)
         call check_allocation(mpierr, 'pulay_mpi                     ')

      end if

!     do i_spin = 1, n_spin, 1
         pulay_mpi = 0.0d0
         call MPI_ALLREDUCE( &
              pulay_matrix(1,1:pulay_saved_iter), &
              pulay_mpi, pulay_saved_iter, &
              MPI_DOUBLE_PRECISION, &
              MPI_SUM, mpi_comm_global, mpierr)
         pulay_matrix(1,1:pulay_saved_iter) = pulay_mpi
!      enddo

      pulay_mpi = 0.0d0
      call MPI_ALLREDUCE( &
           pulay_vector(1:pulay_saved_iter), &
           pulay_mpi, pulay_saved_iter, &
           MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      pulay_vector(1:pulay_saved_iter) = pulay_mpi

      deallocate(pulay_mpi)

      end subroutine sync_pulay_matrix_first_row_and_vector
!******
!--------------------------------------------------------------
!****s* synchronize_mpi/sync_forces
!  NAME
!    sync_forces
!  SYNOPSIS
      subroutine sync_forces(pulay_forces, hellman_feynman_forces, &
           multipole_forces, gga_forces, gga_forces_on, nlcc_forces, nlcc_forces_on,&
           Gnonmf_forces,Gnonmf_forces_on)
!  PURPOSE
!    Synchronize the force components.
!  USES
      use dimensions, only: n_atoms
      implicit none
!  ARGUMENTS
      real*8, dimension(3, n_atoms) :: pulay_forces
      real*8, dimension(3, n_atoms) :: hellman_feynman_forces
      real*8, dimension(3, n_atoms) :: multipole_forces
      real*8, dimension(3, n_atoms) :: gga_forces
      real*8, dimension(3, n_atoms) :: nlcc_forces
      logical :: gga_forces_on
      logical :: nlcc_forces_on
      logical :: Gnonmf_forces_on
      real*8, dimension(3,n_atoms) :: Gnonmf_forces
!  INPUTS
!    o pulay_forces -- Pulay forces
!    o hellman_feynman_forces -- Hellman-Feynman forces
!    o multipole_forces -- multipole forces
!    o gga_forces -- gga-forces
!    o gga_forces_on -- logical indicating if the gga-forces are active
!    o nlcc_forces --  force term related to the nonlinear core correction for KB pseudocores
!    o nlcc_forces_on -- logical indicating if the nlcc-forces are active
!  OUTPUT
!    the force arrays are synchronized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     locals
      real*8, dimension(3, n_atoms) :: forces_mpi
      integer :: mpierr
      integer :: i_atom

      if (.not.use_mpi) return


      forces_mpi = 0.0d0
      call MPI_ALLREDUCE(pulay_forces, &
           forces_mpi, 3*n_atoms, &
           MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      pulay_forces = forces_mpi

      forces_mpi = 0.0d0
      call MPI_ALLREDUCE(hellman_feynman_forces, &
           forces_mpi, 3*n_atoms, &
           MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      hellman_feynman_forces = forces_mpi

      forces_mpi = 0.0d0
      call MPI_ALLREDUCE(multipole_forces, &
           forces_mpi, 3*n_atoms, &
           MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      multipole_forces = forces_mpi

      if (gga_forces_on) then
         forces_mpi = 0.0d0
         call MPI_ALLREDUCE(gga_forces, &
              forces_mpi, 3*n_atoms, &
              MPI_DOUBLE_PRECISION, &
              MPI_SUM, mpi_comm_global, mpierr)
         gga_forces = forces_mpi
      end if

      if (Gnonmf_forces_on) then
         forces_mpi = 0.0d0
         call MPI_ALLREDUCE(Gnonmf_forces, &
              forces_mpi, 3*n_atoms, &
              MPI_DOUBLE_PRECISION, &
              MPI_SUM, mpi_comm_global, mpierr)
         Gnonmf_forces = forces_mpi
      end if      
      
      if (nlcc_forces_on) then
         forces_mpi = 0.0d0
         call MPI_ALLREDUCE(nlcc_forces, &
              forces_mpi, 3*n_atoms, &
              MPI_DOUBLE_PRECISION, &
              MPI_SUM, mpi_comm_global, mpierr)
         nlcc_forces = forces_mpi
      end if


      end subroutine sync_forces
!******
!-------------------------------------------------------------
!****s* synchronize_mpi/sync_moment
!  NAME
!    sync_moment
!  SYNOPSIS
      subroutine sync_moment(moment)
!  PURPOSE
!    Synchronize the electron moment
!  ARGUMENTS
      real*8, dimension(3) :: moment
!  INPUTS
!    o moment -- the electron moment
!  OUTPUT
!    the moment array is synchronized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variable
      real*8, dimension(3) :: mpi_tmp
      integer :: mpierr

      if (.not.use_mpi) return
        mpi_tmp = 0.0d0
        call MPI_ALLREDUCE(moment,mpi_tmp, &
             3, MPI_DOUBLE_PRECISION, &
             MPI_SUM, mpi_comm_global, mpierr)
        moment=mpi_tmp

      end subroutine sync_moment
!******
!-------------------------------------------------------------
!****s* synchronize_mpi/sync_grid_partition
!  NAME
!    sync_grid_partition
!  SYNOPSIS
      subroutine sync_grid_partition(grid_partition, n_grid_batches, &
           n_points_in_grid)
!  PURPOSE
!    Synchronize the grid partitioning information.
!  ARGUMENTS
      integer :: n_points_in_grid
      integer :: n_grid_batches
      integer, dimension(n_points_in_grid) :: grid_partition
!  INPUTS
!    o n_points_in_grid -- number of points in the grid
!    o n_grid_batches -- number of batches of the grid
!    o grid_partition -- partition of the grid into batches
!  OUTPUT
!    the number of grid batches and the grid partition id broadcast to all threads
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      integer :: id_for_bc
      integer :: mpierr

      if (.not.use_mpi) return

      id_for_bc = 0

      call MPI_Bcast(n_grid_batches, 1, MPI_INTEGER, &
           id_for_bc, mpi_comm_global, mpierr)

      call MPI_Bcast(grid_partition, n_points_in_grid, MPI_INTEGER, &
           id_for_bc, mpi_comm_global, mpierr)

      end subroutine sync_grid_partition
!******
!-------------------------------------------------------------
!****s* synchronize_mpi/get_R_prec_spline_mpi
!  NAME
!    get_R_prec_spline_mpi
!  SYNOPSIS
      subroutine get_R_prec_spline_mpi( i_atom, l_dim, &
                 current_R_prec_multipole_spl, R_prec_multipole_spl, &
                 n_grid )
!  PURPOSE
!    Fetch the splined preconditioner componen 
!    for a given atom possibly from a remote storage.
      use dimensions, only: n_max_spline, n_spline_atoms, &
          use_distributed_spline_storage
      use mpi_utilities, only: spline_atom_storage,task_list
      implicit none
!  ARGUMENTS
      integer :: i_atom, l_dim, n_grid
      real*8, dimension(l_dim, n_max_spline, &
            n_grid) :: current_R_prec_multipole_spl
      real*8, dimension(l_dim, n_max_spline, n_grid, &
           n_spline_atoms) :: R_prec_multipole_spl
!  INPUTS
!    o l_dim -- the maximal multipole component the spline array is needed
!    o n_grid -- number of points in the spline grid
!    o R_prec_multipole_spl -- the spline storage array
!    o i_atom -- the atom for which the spline is needed
!  OUTPUT
!    o current_R_prec_multipole_spl -- the fetched spline array
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      integer :: mpierr, myid_for_bc

      if (.not.use_distributed_spline_storage) then
         current_R_prec_multipole_spl(:,:,:) = &
              R_prec_multipole_spl(:,:,:,spline_atom_storage(i_atom))
      else
         myid_for_bc = task_list(i_atom)
         if (myid.eq.myid_for_bc) then
            current_R_prec_multipole_spl(:,:,:) = &
              R_prec_multipole_spl(:,:,:,spline_atom_storage(i_atom))
         end if
         call MPI_Bcast(current_R_prec_multipole_spl, &
              (l_dim*n_max_spline*n_grid), &
              MPI_DOUBLE_PRECISION, &
              myid_for_bc, mpi_comm_global, mpierr)
      end if
      end subroutine get_R_prec_spline_mpi




!******
!---------------------------------------------------------------
!****s* synchronize_mpi/sync_kerker_multipole
!  NAME
!    sync_kerker_multipole
!  SYNOPSIS
      subroutine sync_kerker_multipole(R_multipole, &
                 radial_grid_count, l_dim)
!  PURPOSE
!    Synchronize the multipole arrays for the Kerker preconditioner
      use mpi_tasks, only: mpi_in_place
      use dimensions, only: n_atoms
      implicit none
!  ARGUMENTS
      integer :: radial_grid_count, l_dim
      real*8 :: R_multipole( l_dim,radial_grid_count, n_atoms )
!  INPUTS
!    o l_dim -- the maximal multipole component
!    o radial_grid_count -- number of points in the radial grid
!    o R_multipole -- the multipole storage array
!  OUTPUT
!    the array R_multipole is synchronized over all tasks.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      real*8, dimension(:,:), allocatable :: R_multipole_MPI
      integer :: i_atom, mpierr

      if (.not.use_mpi) return

      if (.not.use_mpi_in_place) then
         allocate(R_multipole_mpi( l_dim,radial_grid_count ),stat=mpierr)
         call check_allocation(mpierr, 'R_multipole_mpi               ')

         do i_atom = 1, n_atoms, 1
            R_multipole_mpi = 0.0d0
            call MPI_ALLREDUCE( &
                 R_multipole(:,:,i_atom), &
                 R_multipole_MPI, &
                 l_dim * radial_grid_count, &
                 MPI_DOUBLE_PRECISION, &
                 MPI_SUM, &
                 mpi_comm_global, &
                 mpierr )
            R_multipole(:,:,i_atom) = R_multipole_mpi
         end do
         deallocate(R_multipole_mpi)
      else
         call MPI_ALLREDUCE(MPI_IN_PLACE, &
              R_multipole(:,:,:), &
              l_dim * radial_grid_count * n_atoms, &
              MPI_DOUBLE_PRECISION, &
              MPI_SUM, &
              mpi_comm_global, &
              mpierr )
      end if

      end subroutine sync_kerker_multipole
!******


!---------------------------------------------------------------
!****s* synchronize_mpi/sync_kerker_multipole_splines
!  NAME
!    sync_kerker_multipole_splines
!  SYNOPSIS
      subroutine sync_kerker_multipole_splines(R_multipole_spl, &
                 l_dim, n_max_spline, grid_dim, n_atoms)
!  PURPOSE
!    Synchronize the multipole spline arrays for the Kerker preconditioner
      use mpi_tasks, only: mpi_in_place
      implicit none
!  ARGUMENTS
      integer :: l_dim, n_max_spline, grid_dim, n_atoms
!  INPUTS
!    o l_dim -- the maximal multipole component
!    o grid_dim -- number of points in the spline grid
!    o R_multipole_spl -- the multipole spline storage array
!    o n_max_spline -- maximal order for splines
!    o n_atoms -- number of atoms in the system
!  OUTPUT
!    the array R_multipole_spl is synchronized over all tasks.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      real*8 :: R_multipole_spl( l_dim, n_max_spline, &
                                          grid_dim, n_atoms )
      real*8, dimension(:,:,:), allocatable :: R_multipole_MPI
      integer :: i_atom, mpierr

      if (.not.use_mpi) return

      if (.not.use_mpi_in_place) then
         allocate(R_multipole_MPI(l_dim,n_max_spline,grid_dim),stat=mpierr )
         call check_allocation(mpierr, 'R_multipole_MPI               ')

         do i_atom = 1, n_atoms, 1
            R_multipole_mpi = 0.0d0
            call MPI_ALLREDUCE( &
                 R_multipole_spl(:,:,:,i_atom), &
                 R_multipole_MPI, &
                 l_dim*n_max_spline*grid_dim, &
                 MPI_DOUBLE_PRECISION, &
                 MPI_SUM, &
                 mpi_comm_global, &
                 mpierr )
            R_multipole_spl(:,:,:,i_atom) = R_multipole_mpi
         end do
         deallocate(R_multipole_mpi)
      else
         call MPI_ALLREDUCE(MPI_IN_PLACE, &
              R_multipole_spl(:,:,:,:), &
              l_dim*n_max_spline*grid_dim*n_atoms, &
              MPI_DOUBLE_PRECISION, &
              MPI_SUM, &
              mpi_comm_global, &
              mpierr )
      end if

      end subroutine sync_kerker_multipole_splines
!******
!-------------------------------------------------------------
!****s* synchronize_mpi/sync_evec_scalapack
!  NAME
!    sync_evec_scalapack
!  SYNOPSIS
      subroutine sync_evec_scalapack( KS_eigenvector, &
           n_basis, n_states, scalapack_comm )
!  PURPOSE
!    Synchronize the eigenvectors for ScaLAPACK
!  USES
      implicit none
!  ARGUMENTS
      integer :: n_basis, n_states
      real*8 :: KS_eigenvector( n_basis, n_states )
      integer :: scalapack_comm
!  INPUTS
!    o n_basis -- number of basis functions
!    o n_states -- number of states
!    o KS_eigenvector -- Kohn-Sham eigenvectors
!    o scalapack_comm -- the ScaLAPACK communicator
!  OUTPUT
!    the array KS_eigenvector is synchronized over tasks in scalapack_comm.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(KS_eigenvector,n_basis*n_states, &
           scalapack_comm)

      end subroutine sync_evec_scalapack
!******
!-------------------------------------------------------------
!****s* synchronize_mpi/sync_complex_evec_scalapack
!  NAME
!    sync_complex_evec_scalapack
!  SYNOPSIS
      subroutine sync_complex_evec_scalapack( KS_eigenvector, &
           scalapack_comm )
!  USES
      use dimensions, only: n_basis, n_states
      implicit none
!  PURPOSE
!    Synchronize the complex eigenvectors for ScaLAPACK
!  ARGUMENTS
      complex*16 :: KS_eigenvector( n_basis, n_states )
      integer :: scalapack_comm
!  INPUTS
!    o KS_eigenvector -- complex Kohn-Sham eigenvectors
!    o scalapack_comm -- the ScaLAPACK communicator
!  OUTPUT
!    the array KS_eigenvector is synchronized over tasks in scalapack_comm.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector_complex(KS_eigenvector,n_basis*n_states, &
           scalapack_comm)

      end subroutine sync_complex_evec_scalapack
!******
!------------------------------------------------------
!****s* synchronize_mpi/sync_restart
!  NAME
!    sync_restart
!  SYNOPSIS
      subroutine sync_restart(KS_eigenvector,KS_eigenvalue,occ_numbers, &
           n_basis, n_states, n_spin, start_force_occ, &
           previous_eigenvector)
!  PURPOSE
!    Synchronize restart information.
!  USES
      implicit none
!  ARGUMENTS
      integer :: n_basis, n_states, n_spin
      real*8 :: KS_eigenvector(n_basis,n_states,n_spin)
      real*8 :: KS_eigenvalue(n_states,n_spin)
      real*8 :: occ_numbers(n_states,n_spin)
      logical,optional :: start_force_occ
      real*8,optional :: previous_eigenvector(n_basis,n_states,n_spin)
!  INPUTS
!    o n_basis -- number of basis functions
!    o n_states -- number of states
!    o n_spin -- number of spin channels
!    o KS_eigenvector -- Kohn-Sham eigenvectors
!    o KS_eigenvalue -- Kohn-Sham eigenvalues
!    o occ_numbers -- occupation numbers
!  OUTPUT
!    the arrays are synchronized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      real*8 :: KS_eigenvector_MPI(n_basis,n_states,n_spin)
      real*8 :: KS_eigenvalue_MPI(n_states,n_spin)
      real*8 :: occ_numbers_MPI(n_states,n_spin)
      integer :: mpierr

      if(present(start_force_occ)) then
        call MPI_Bcast(start_force_occ, &
             1, &
             MPI_LOGICAL, &
             0 , mpi_comm_global, mpierr)
        call MPI_Bcast(previous_eigenvector, &
             n_basis*n_states*n_spin, &
             MPI_DOUBLE_PRECISION, &
             0 , mpi_comm_global, mpierr)
      end if ! synced force_occupation variables

      if (.not.use_mpi) return
!     initialize arrays to zero, just in case
      if (myid.ne.0) then
         KS_eigenvector = 0d0
         KS_eigenvalue = 0d0
         occ_numbers = 0d0
      end if
      KS_eigenvector_MPI = 0d0
      KS_eigenvalue_MPI = 0d0
      occ_numbers_MPI = 0d0
      call MPI_ALLREDUCE(KS_eigenvector,KS_eigenvector_MPI, &
          n_basis*n_states*n_spin, MPI_DOUBLE_PRECISION, MPI_SUM, &
          mpi_comm_global, mpierr)
      call MPI_ALLREDUCE(KS_eigenvalue,KS_eigenvalue_MPI, &
          n_states*n_spin, MPI_DOUBLE_PRECISION, MPI_SUM, &
          mpi_comm_global, mpierr)
      call MPI_ALLREDUCE(occ_numbers,occ_numbers_MPI, &
          n_states*n_spin, MPI_DOUBLE_PRECISION, MPI_SUM, &
          mpi_comm_global, mpierr)
      KS_eigenvector = KS_eigenvector_MPI
      KS_eigenvalue  = KS_eigenvalue_MPI
      occ_numbers    = occ_numbers_MPI
      end subroutine sync_restart
!******
!------------------------------------------------------
!****s* synchronize_mpi/sync_periodic_restart
!  NAME
!    sync_periodic_restart
!  SYNOPSIS
      subroutine sync_periodic_restart(KS_eigenvector,KS_eigenvalue,occ_numbers, &
           n_basis, n_states, n_spin, n_k_points, start_force_occ, &
           previous_eigenvector)
!  PURPOSE
!    Synchronize restart information.
      implicit none
!  ARGUMENTS
      integer :: n_basis, n_states, n_spin, n_k_points
      real*8 :: KS_eigenvector(n_basis,n_states,n_spin,n_k_points)
      real*8 :: KS_eigenvalue(n_states,n_spin, n_k_points)
      real*8 :: occ_numbers(n_states,n_spin, n_k_points)
      logical,optional :: start_force_occ
      real*8,optional :: previous_eigenvector(n_basis,n_states,n_spin,n_k_points)
!  INPUTS
!    o n_basis -- number of basis functions
!    o n_states -- number of states
!    o n_spin -- number of spin channels
!    o KS_eigenvector -- Kohn-Sham eigenvectors
!    o KS_eigenvalue -- Kohn-Sham eigenvalues
!    o occ_numbers -- occupation numbers
!  OUTPUT
!    the arrays are synchronized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      real*8 :: KS_eigenvector_MPI(n_basis,n_states,n_spin,n_k_points)
      real*8 :: KS_eigenvalue_MPI(n_states,n_spin,n_k_points)
      real*8 :: occ_numbers_MPI(n_states,n_spin,n_k_points)
      integer :: mpierr

      if(present(start_force_occ)) then
        call MPI_Bcast(start_force_occ, &
             1, &
             MPI_LOGICAL, &
             0 , mpi_comm_global, mpierr)
        call MPI_Bcast(previous_eigenvector, &
             n_basis*n_states*n_spin*n_k_points, &
             MPI_DOUBLE_PRECISION, &
             0 , mpi_comm_global, mpierr)
      end if ! synced force_occupation variables

      if (.not.use_mpi) return
!     initialize arrays to zero, just in case
      if (myid.ne.0) then
         KS_eigenvector = 0d0
         KS_eigenvalue = 0d0
         occ_numbers = 0d0
      end if
      KS_eigenvector_MPI = 0d0
      KS_eigenvalue_MPI = 0d0
      occ_numbers_MPI = 0d0
      call MPI_ALLREDUCE(KS_eigenvector,KS_eigenvector_MPI, &
          n_basis*n_states*n_spin*n_k_points, MPI_DOUBLE_PRECISION, MPI_SUM, &
          mpi_comm_global, mpierr)
      call MPI_ALLREDUCE(KS_eigenvalue,KS_eigenvalue_MPI, &
          n_states*n_spin*n_k_points, MPI_DOUBLE_PRECISION, MPI_SUM, &
          mpi_comm_global, mpierr)
      call MPI_ALLREDUCE(occ_numbers,occ_numbers_MPI, &
          n_states*n_spin*n_k_points, MPI_DOUBLE_PRECISION, MPI_SUM, &
          mpi_comm_global, mpierr)
      KS_eigenvector = KS_eigenvector_MPI
      KS_eigenvalue  = KS_eigenvalue_MPI
      occ_numbers    = occ_numbers_MPI
      end subroutine sync_periodic_restart
!******
!******
!------------------------------------------------------
!****s* synchronize_mpi/sync_eigenvector_complex
!  NAME
!    sync_eigenvector_complex
!  SYNOPSIS
      subroutine sync_eigenvector_complex(KS_eigenvector)
!  PURPOSE
!    Synchronize complex KS_eigenvectors for only one k-point.
!    This routine is only used to get one eigenvector after another
!    onto MPI task zero, for output purposes. The idea is to be
!    as memory-efficient as possible, hence only one eigenvector at a time.
      use dimensions, only: n_basis, n_states, n_spin
      implicit none
!  ARGUMENTS
      complex*16 :: KS_eigenvector(n_basis,n_states,n_spin)
!  INPUTS
!    o KS_eigenvector -- Kohn-Sham eigenvectors
!  OUTPUT
!    the array is synchronized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2010).
!  SOURCE

      complex*16 :: KS_eigenvector_MPI(n_basis,n_states,n_spin)
      integer :: mpierr

      if (.not.use_mpi) return

      KS_eigenvector_MPI = 0d0
      call MPI_ALLREDUCE(KS_eigenvector,KS_eigenvector_MPI, &
          n_basis*n_states*n_spin, MPI_DOUBLE_COMPLEX, MPI_SUM, &
          mpi_comm_global, mpierr)
      KS_eigenvector = KS_eigenvector_MPI

    end subroutine sync_eigenvector_complex
!******
!---------------------------------------------------------------
!****s* synchronize_mpi/sync_aux_selfenergy
!  NAME
!    sync_aux_selfenergy
!  SYNOPSIS
      subroutine sync_aux_selfenergy &
                 ( n_low_state,n_high_state, &
                   n_freq, aux_selfenergy &
                  )
!  PURPOSE
!    ---
!  USES
      use dimensions, only: n_states, n_spin
!  ARGUMENTS
      integer ::  n_low_state, n_high_state
      integer ::  n_freq
      real*8 :: aux_selfenergy &
                (n_low_state:n_high_state,n_states,n_freq,n_spin)
!  INPUTS
!    o n_low_state --
!    o n_high_state --
!    o n_freq --
!    o aux_selfenergy --
!  OUTPUT
!    the array is synchronized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables
      integer :: n_KS_states
      integer :: i_freq
      integer :: i_spin

      real*8, dimension(:,:), allocatable :: selfenergy_mpi
      real*8, dimension(:,:), allocatable :: selfenergy_tmp
      integer :: mpierr

      if (.not.use_mpi) return

      n_KS_states = n_high_state - n_low_state + 1

      if (.not.allocated(selfenergy_mpi)) then
         allocate(selfenergy_mpi(n_KS_states,n_states),stat=mpierr)
         call check_allocation(mpierr, 'selfenergy_mpi                ')
      end if

      if (.not.allocated(selfenergy_tmp)) then
         allocate(selfenergy_tmp(n_KS_states,n_states),stat=mpierr)
         call check_allocation(mpierr, 'selfenergy_tmp                ')
      end if

      do i_spin = 1, n_spin
        do i_freq = 1, n_freq

          selfenergy_tmp(:,:) = &
             aux_selfenergy(n_low_state:n_high_state,:,i_freq,i_spin)

          call MPI_ALLREDUCE(selfenergy_tmp, selfenergy_mpi, &
            n_KS_states*n_states, MPI_DOUBLE_PRECISION, &
            MPI_SUM, mpi_comm_global, mpierr)

            aux_selfenergy(n_low_state:n_high_state,:,i_freq,i_spin) &
              = selfenergy_mpi(:,:)
         enddo
       enddo

       deallocate(selfenergy_tmp)
       deallocate(selfenergy_mpi)

       end subroutine sync_aux_selfenergy
!******
!--------------------------------------------------------------------------
!****s* synchronize_mpi/sync_hirshfeld
!  NAME
!    sync_hirshfeld
!  SYNOPSIS
       subroutine sync_hirshfeld  &
      ( hirshfeld_charge, hirshfeld_spin_moment, &
       hirshfeld_dipole, hirshfeld_quadrupole, free_int, hirsh_int )
!  PURPOSE
!    ---
       use runtime_choices, only: spin_treatment
       use dimensions, only: n_atoms, n_atoms
       implicit none
!  ARGUMENTS
       real*8 :: hirshfeld_charge ( n_atoms )
       real*8 :: free_int ( n_atoms )
       real*8 :: hirsh_int ( n_atoms )
       real*8 :: hirshfeld_spin_moment ( n_atoms )
       real*8 :: hirshfeld_dipole ( 3, n_atoms )
       real*8 :: hirshfeld_quadrupole ( 6, n_atoms )
!  INPUTS
!    o hirshfield_charge --
!    o free_int --
!    o hirsh_int --
!    o hirshfield_spin_moment --
!    o hirshfield_dipole --
!    o hirshfield_quadrupole --
!  OUTPUT
!    the arrays are synchronized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
! SOURCE


!      local variables
       real*8 :: charge_mpi     (    n_atoms )
       real*8 :: spin_mpi     (    n_atoms )
       real*8 :: dipole_mpi     ( 3, n_atoms )
       real*8 :: quadrupole_mpi ( 6, n_atoms )
       integer :: mpierr
!      begin work

       if (.not.use_mpi) return

       charge_mpi = 0.d0
       call MPI_ALLREDUCE &
      ( hirshfeld_charge, & 
        charge_mpi, &
        n_atoms, &
        MPI_DOUBLE_PRECISION, MPI_SUM, &
        mpi_comm_global, mpierr )
       hirshfeld_charge = charge_mpi

       charge_mpi = 0.d0
       call MPI_ALLREDUCE &
      ( free_int, &
        charge_mpi, &
        n_atoms, &
        MPI_DOUBLE_PRECISION, MPI_SUM, &
        mpi_comm_global, mpierr ) 
       free_int = charge_mpi

       charge_mpi = 0.d0
       call MPI_ALLREDUCE &
      ( hirsh_int, &
        charge_mpi, &
        n_atoms, &
        MPI_DOUBLE_PRECISION, MPI_SUM, &
        mpi_comm_global, mpierr )
       hirsh_int = charge_mpi

       if (spin_treatment .eq. 1) then
        spin_mpi = 0.d0
        call MPI_ALLREDUCE &
      ( hirshfeld_spin_moment, & 
        spin_mpi, &
        n_atoms, &
        MPI_DOUBLE_PRECISION, MPI_SUM, &
        mpi_comm_global, mpierr )
        hirshfeld_spin_moment = spin_mpi
       endif

       dipole_mpi = 0.d0
       call MPI_ALLREDUCE &
      ( hirshfeld_dipole, &
        dipole_mpi, &
        3*n_atoms, &
        MPI_DOUBLE_PRECISION, MPI_SUM, &
        mpi_comm_global, mpierr )
       hirshfeld_dipole = dipole_mpi

       quadrupole_mpi = 0.d0
       call MPI_ALLREDUCE &
      ( hirshfeld_quadrupole, & 
        quadrupole_mpi, &
        6*n_atoms, &
        MPI_DOUBLE_PRECISION, MPI_SUM, &
        mpi_comm_global, mpierr )
       hirshfeld_quadrupole = quadrupole_mpi

       end subroutine sync_hirshfeld
!******
!--------------------------------------------------------------------------
!****s* synchronize_mpi/broadcast_MD_velocities
!  NAME
!    broadcast_MD_velocities
!  SYNOPSIS
      subroutine broadcast_MD_velocities(velocities, root_task)
!  PURPOSE
!    Broadcast the MD velocities from a given task to all other tasks.
!  USES
      use dimensions
      implicit none
!  ARGUMENTS
      real*8 :: velocities(3,n_atoms)
      integer :: root_task
!  INPUTS
!    o velocities -- the MD velocities
!    o root_task -- the broadcasting task
!  OUTPUT
!    the velocities are broadcast from root_task to all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      integer :: mpierr

      if (.not.use_mpi) return

      call MPI_Bcast(velocities, 3*n_atoms, &
               MPI_DOUBLE_PRECISION, &
               root_task, mpi_comm_global, mpierr)

      end subroutine

!******
!--------------------------------------------------------------------------
!****s* synchronize_mpi/broadcast_PIMD_velocities
!  NAME
!    broadcast_PIMD_velocities
!  SYNOPSIS
      subroutine broadcast_PIMD_velocities(velocities, root_task)
!  PURPOSE
!    Broadcast the PIMD velocities from a given task to all other tasks.
!  USES
      use dimensions
      implicit none
!  ARGUMENTS
      real*8 :: velocities(3,n_beads,n_atoms)
      integer :: root_task
!  INPUTS
!    o velocities -- the PIMD velocities
!    o root_task -- the broadcasting task
!  OUTPUT
!    the velocities are broadcast from root_task to all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      integer :: mpierr

      if (.not.use_mpi) return

      call MPI_Bcast(velocities, 3*n_beads*n_atoms, &
               MPI_DOUBLE_PRECISION, &
               root_task, mpi_comm_global, mpierr)

      end subroutine

      subroutine sync_tensors(temp_matrix,num_ele)
      use dimensions
      implicit none
      integer :: num_ele
      real*8,dimension(num_ele,num_ele):: temp_matrix ,temp_matrix_mpi
      integer :: mpierr
      if (.not.use_mpi) return
      call MPI_ALLREDUCE(temp_matrix,temp_matrix_mpi,&
           num_ele*num_ele,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierr)

      temp_matrix = temp_matrix_mpi
      endsubroutine

      subroutine sync_vdw_correction_in_scs(en_vdw)
      ! VVG sync vdw energy from SCS routine 
      real*8 :: en_vdw
      !local variable 
      real*8 :: en_vdw_mpi
      integer :: mpierr

      if (.not.use_mpi) return

        en_vdw_mpi = 0.0d0
        call MPI_ALLREDUCE(en_vdw, en_vdw_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
        en_vdw = en_vdw_mpi
      end subroutine sync_vdw_correction_in_scs


!----------------------------------------------------------
!****s* synchronize_mpi/sync_sum_up_whole_potential_shanghui
!  NAME
!    sync_sum_up_whole_potential_shanghui
!  SYNOPSIS
      subroutine sync_sum_up_whole_potential_shanghui( &
           hartree_multipole_error) 
!  PURPOSE
!     here shanghui only sync hartree_multipole_error
!  ARGUMENTS
      real*8 :: hartree_multipole_error
!  INPUTS
!    o hartree_multipole_error --
!  OUTPUT
!    the quantities are synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables
      real*8 :: hartree_multipole_error_mpi

      integer :: mpierr


!     broadcast the result to all threads

      if (.not.use_mpi) return

      hartree_multipole_error_mpi = 0.0d0
      call MPI_ALLREDUCE(hartree_multipole_error, &
           hartree_multipole_error_mpi, &
           1, MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      hartree_multipole_error = &
           hartree_multipole_error_mpi


      end subroutine sync_sum_up_whole_potential_shanghui


!----------------------------------------------------------
!****s* synchronize_mpi/sync_first_order_matrix
!  NAME
!    sync_first_order_matrix
!  SYNOPSIS
      subroutine  sync_first_order_matrix(first_order_M )  
!  PURPOSE
!    Synchronize the first_order_M matrix after integration.
!  USES
      use dimensions, only: n_Cbasis
      implicit none
!  ARGUMENTS
      real*8 :: first_order_M(n_Cbasis,n_Cbasis)
!  INPUTS
!    o first_order_M 
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      call sync_vector(first_order_M, n_Cbasis*n_Cbasis)

      end subroutine  sync_first_order_matrix


!----------------------------------------------------------
!****s* synchronize_mpi/sync_integrate_first_order_S
!  NAME
!    sync_integrate_first_order_S
!  SYNOPSIS
      subroutine  sync_integrate_first_order_S(first_order_S )  

!  PURPOSE
!    Synchronize the first_order_S matrix after integration.
!  USES
      use dimensions, only: n_atoms, n_Cbasis
      implicit none
!  ARGUMENTS
      real*8 :: first_order_S(3,n_atoms,n_Cbasis,n_Cbasis)
!  INPUTS
!    o first_order_S 
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(first_order_S,3*n_atoms*n_Cbasis*n_Cbasis)


      end subroutine  sync_integrate_first_order_S

!----------------------------------------------------------
!****s* synchronize_mpi/sync_integrate_second_order_S
!  NAME
!    sync_integrate_second_order_S
!  SYNOPSIS
      subroutine  sync_integrate_second_order_S(second_order_S )  

!  PURPOSE
!    Synchronize the second_order_S matrix after integration.
!  USES
      use dimensions, only: n_atoms, n_Cbasis
      implicit none
!  ARGUMENTS
      real*8 :: second_order_S(3,n_atoms,3,n_atoms,n_Cbasis,n_Cbasis)
!  INPUTS
!    o second_order_S 
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(second_order_S, & 
           3*n_atoms*3*n_atoms*n_Cbasis*n_Cbasis)


      end subroutine  sync_integrate_second_order_S

!----------------------------------------------------------
!****s* synchronize_mpi/sync_integrate_first_order_H
!  NAME
!    sync_integrate_first_order_H
!  SYNOPSIS
      subroutine  sync_integrate_first_order_H(first_order_H )  

!  PURPOSE
!    Synchronize the first_order_H matrix after integration.
!  USES
      use dimensions, only: n_atoms, n_Cbasis
      implicit none
!  ARGUMENTS
      real*8 :: first_order_H(3,n_atoms,n_Cbasis,n_Cbasis)
!  INPUTS
!    o first_order_H
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(first_order_H,3*n_atoms*n_Cbasis*n_Cbasis)


      end subroutine  sync_integrate_first_order_H

!----------------------------------------------------------
!****s* synchronize_mpi/sync_integrate_second_order_H_pulay
!  NAME
!    sync_integrate_second_order_H_pulay
!  SYNOPSIS
      subroutine  sync_integrate_second_order_H_pulay & 
                  (second_order_H_pulay )  

!  PURPOSE
!    Synchronize the second_order_H_pulay matrix after integration.
!  USES
      use dimensions, only: n_atoms, n_Cbasis
      implicit none
!  ARGUMENTS
      real*8 :: second_order_H_pulay & 
               (3,n_atoms,3,n_atoms,n_Cbasis,n_Cbasis)
!  INPUTS
!    o second_order_H_pulay
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(second_order_H_pulay, & 
           3*n_atoms*3*n_atoms*n_Cbasis*n_Cbasis)


      end subroutine  sync_integrate_second_order_H_pulay

!----------------------------------------------------------
!****s* synchronize_mpi/sync_integrate_hellman_feynman_hessian
!  NAME
!    sync_integrate_hellman_feynman_hessian
!  SYNOPSIS
      subroutine  sync_integrate_hellman_feynman_hessian & 
                  (hellman_feynman_hessian)  

!  PURPOSE
!    Synchronize the second_order_H_pulay matrix after integration.
!  USES
      use dimensions, only: n_atoms
      implicit none
!  ARGUMENTS
      real*8 :: hellman_feynman_hessian & 
               (3,n_atoms,3,n_atoms)
!  INPUTS
!    o hellman_feynman_hessian
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(hellman_feynman_hessian, & 
           3*n_atoms*3*n_atoms)


      end subroutine  sync_integrate_hellman_feynman_hessian


!----------------------------------------------------------
!****s* synchronize_mpi/sync_integrate_gradient_dipole
!  NAME
!    sync_integrate_gradient_dipole
!  SYNOPSIS
      subroutine  sync_integrate_gradient_dipole & 
                  (gradient_dipole)  

!  PURPOSE
!    Synchronize the second_order_H_pulay matrix after integration.
!  USES
      use dimensions, only: n_atoms
      implicit none
!  ARGUMENTS
      real*8 :: gradient_dipole & 
               (3*n_atoms, 3)
!  INPUTS
!    o gradient_dipole
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(gradient_dipole, & 
           3*n_atoms*3)


      end subroutine  sync_integrate_gradient_dipole


!----------------------------------------------------------
!****s* synchronize_mpi/sync_integrate_polarizability
!  NAME
!    sync_integrate_polarizability
!  SYNOPSIS
      subroutine  sync_integrate_polarizability(polarizability )  

!  PURPOSE
!    Synchronize the polarizability after integration.
!  ARGUMENTS
      real*8 :: polarizability(3,3)
!  INPUTS
!    o polarizability
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(polarizability,3*3)


      end subroutine  sync_integrate_polarizability


!----------------------------------------------------------
!****s* synchronize_mpi/sync_integrate_first_order_H_polarizability
!  NAME
!    sync_integrate_first_order_H_polarizability
!  SYNOPSIS
      subroutine  sync_integrate_first_order_H_polarizability & 
                    (first_order_H )  

!  PURPOSE
!    Synchronize the first_order_H matrix after integration.
!  USES
      use dimensions, only: n_basis, n_spin
      implicit none
!  ARGUMENTS
      real*8 :: first_order_H(3,n_basis,n_basis,n_spin)
!  INPUTS
!    o first_order_H
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(first_order_H,3*n_basis*n_basis*n_spin)


      end subroutine  sync_integrate_first_order_H_polarizability

!----------------------------------------------------------
!****s* synchronize_mpi/sync_frist_order_density_matrix
!  NAME
!    sync_frist_order_density_matrix
!  SYNOPSIS
      subroutine  sync_frist_order_density_matrix & 
                 (first_order_density_matrix )  
!  PURPOSE
!    Synchronize the first_order_density_matrix over n_k_task
!  USES
      use dimensions, only: n_atoms, n_Cbasis
      implicit none
!  ARGUMENTS
      real*8 :: first_order_density_matrix(3,n_atoms,n_Cbasis,n_Cbasis)
!  INPUTS
!    o first_order_density_matrix
!  OUTPUT
!    the array is synchronized over tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      call sync_vector(first_order_density_matrix,3*n_atoms*n_Cbasis*n_Cbasis)

      end subroutine  sync_frist_order_density_matrix

!----------------------------------------------------------
!****s* synchronize_mpi/sync_atomic_solver_outputs
!  NAME
!    sync_atomic_solver_outputs
!  SYNOPSIS
      subroutine sync_atomic_solver_outputs( n_max_ind_fns, n_max_grid, n_l_channels, n_spin, &
           eigenval, wave, wave_deriv, kinetic, &
           density, density_deriv, density_2nd_deriv, atom_pot, atom_pot_es, wave0)
!  PURPOSE
!    Synchronize the outputs from an invocation of an atomic solver by having 
!    myid.eq.0 communicate its results to all other ranks.  Only needed if
!    the relevant atomic solver implementation in aims is non-deterministic.
!  USES
      implicit none
!  ARGUMENTS
      integer, intent(in)                                                 :: n_max_ind_fns
      integer, intent(in)                                                 :: n_max_grid
      integer, intent(in)                                                 :: n_l_channels
      integer, intent(in)                                                 :: n_spin
      real*8,  intent(inout), dimension(n_max_ind_fns)                    :: eigenval
      real*8,  intent(inout), dimension(n_max_grid, n_max_ind_fns)        :: wave
      real*8,  intent(inout), dimension(n_max_grid, n_max_ind_fns)        :: wave_deriv
      real*8,  intent(inout), dimension(n_max_grid, n_max_ind_fns)        :: kinetic
      real*8,  intent(inout), dimension(n_max_grid, n_spin)               :: density
      real*8,  intent(inout), dimension(n_max_grid, n_spin)               :: density_deriv
      real*8,  intent(inout), dimension(n_max_grid, n_spin)               :: density_2nd_deriv
      real*8,  intent(inout), dimension(n_max_grid, n_l_channels*n_spin)  :: atom_pot
      real*8,  intent(inout), dimension(n_max_grid)                       :: atom_pot_es
      real*8, intent(in out), optional :: wave0(n_max_grid)
!  INPUTS
!    o n_max_ind_fns
!    o n_max_grid
!    o n_channels
!    o eigenval
!    o wave
!    o wave_deriv
!    o density
!    o density_deriv
!    o density_2nd_deriv
!    o atom_pot
!    o atom_pot_es
!    o kinetic
!    o wave0
!  OUTPUTS
!    o eigenval
!    o wave
!    o wave_deriv
!    o density
!    o density_deriv
!    o density_2nd_deriv
!    o atom_pot
!    o atom_pot_es
!    o kinetic
!    o wave0
!  NOTE
!    o As of this writing, aims supports sratom, rdirac, and atom_sphere.
!      Of these three, this subroutine is only needed for atom_sphere due to
!      its use of random_number().
!  AUTHOR
!    William Huhn, Duke University
!  HISTORY
!    Added, November 2016
!  SOURCE
        call bcast_real_vector( eigenval,          n_max_ind_fns,                  0 )
        call bcast_real_vector( wave,              n_max_grid*n_max_ind_fns,       0 )
        call bcast_real_vector( wave_deriv,        n_max_grid*n_max_ind_fns,       0 )
        call bcast_real_vector( kinetic,           n_max_grid*n_max_ind_fns,       0 )
        call bcast_real_vector( density,           n_max_grid*n_spin,              0 )
        call bcast_real_vector( density_deriv,     n_max_grid*n_spin,              0 )
        call bcast_real_vector( density_2nd_deriv, n_max_grid*n_spin,              0 )
        call bcast_real_vector( atom_pot,          n_max_grid*n_l_channels*n_spin, 0 )
        call bcast_real_vector( atom_pot_es,       n_max_grid,                     0 )
        if (present(wave0)) call bcast_real_vector(wave0, n_max_ind_fns, 0)

      end subroutine sync_atomic_solver_outputs
!------------------------------------------------------
!****s* synchronize_mpi/sync_force_occupation
!  NAME
!    sync_force_occupation
!  SYNOPSIS
      subroutine sync_force_occupation(force_occ_pr_state, &
                        force_occ_state, force_occ_spin, &
                        forced_occ_number, force_occ_min_KS_state, &
                        force_occ_max_KS_state, n_force_occ, &
                        i_k_point,force_occ_pr_state_periodic, &
                        force_occ_state_periodic)
!  PURPOSE
!    Synchronize force_occupation information.
!  USES

      use dimensions,                   only : n_periodic
      use mpi_tasks,                    only : n_tasks

      implicit none
!  ARGUMENTS
      integer, intent(in)               :: n_force_occ
      integer, intent(in) , optional    :: i_k_point
      integer           :: force_occ_pr_state(n_force_occ)
      integer, optional :: force_occ_pr_state_periodic(:,:)
      integer           :: force_occ_state(n_force_occ)
      integer, optional :: force_occ_state_periodic(:,:)
      integer           :: force_occ_spin(n_force_occ)
      real*8            :: forced_occ_number(n_force_occ)
      integer           :: force_occ_min_KS_state(n_force_occ)
      integer           :: force_occ_max_KS_state(n_force_occ)
!  INPUTS
!    o n_force_occ              -- number of force_occupations
!    o force_occ_pr_state       -- previously forced state
!    o force_occ_pr_state       -- currently forced state
!    o force_occ_spin           -- spin channel force_occupations
!    o forced_occ_number        -- the actually enforced occupation
!    o force_occ_min_KS_state   -- lower bound of MOM window
!    o force_occ_max_KS_state   -- upper bound of MOM window
!  OUTPUT
!    the arrays are synchronized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      integer :: mpierr, root_task

      if (.not. use_mpi) return

      root_task = MOD(i_k_point, n_tasks)

      if (n_periodic .eq. 0) then

        call MPI_Bcast(force_occ_pr_state, &
             n_force_occ, &
             MPI_INTEGER, &
             root_task, mpi_comm_global, mpierr)

        call MPI_Bcast(force_occ_state, &
             n_force_occ, &
             MPI_INTEGER, &
             root_task, mpi_comm_global, mpierr)
      else
        call MPI_ALLREDUCE(force_occ_pr_state_periodic, &
             force_occ_pr_state_periodic(:, i_k_point), &
             n_force_occ, MPI_INTEGER, &
             MPI_SUM, mpi_comm_global, mpierr)

        call MPI_ALLREDUCE(force_occ_state_periodic, &
             force_occ_state_periodic(:, i_k_point), &
             n_force_occ, MPI_INTEGER, &
             MPI_SUM, mpi_comm_global, mpierr)
      end if

      call MPI_Bcast(force_occ_spin, &
           n_force_occ, &
           MPI_INTEGER, &
           root_task, mpi_comm_global, mpierr)

      call MPI_Bcast(forced_occ_number, &
           n_force_occ, &
           MPI_DOUBLE_PRECISION, &
           root_task, mpi_comm_global, mpierr)

      call MPI_Bcast(force_occ_min_KS_state, &
           n_force_occ, &
           MPI_INTEGER, &
           root_task, mpi_comm_global, mpierr)

      call MPI_Bcast(force_occ_max_KS_state, &
           n_force_occ, &
           MPI_INTEGER, &
           root_task, mpi_comm_global, mpierr)

      end subroutine sync_force_occupation
!******

!------------------------------------------------------
!****s* synchronize_mpi/bcast_eigenvector_boys
!  NAME
!    sync_force_occupation
!  SYNOPSIS
      subroutine bcast_eigenvector_boys(KS_eigenvector, dim)
!  PURPOSE
!    Synchronize force_occupation information.
!  USES

      use dimensions,                   only : n_periodic
      use mpi_tasks,                    only : n_tasks

      implicit none
!  ARGUMENTS
      integer, intent(in)                       :: dim
      real*8, dimension(:,:,:,:), intent(in)    :: KS_eigenvector
!  INPUTS
!    o vec   -- KS_eigenvector to sync
!    o dim   -- dimension of KS_eigenvector
!  OUTPUT
!    the arrays are synchronized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      integer :: mpierr, root_task

      if (.not. use_mpi) return

      root_task = 0

      if (n_periodic .eq. 0) then

        call MPI_Bcast(KS_eigenvector, &
             dim, &
             MPI_DOUBLE_PRECISION, &
             root_task, mpi_comm_global, mpierr)

      else
        call aims_stop()
      end if

      end subroutine bcast_eigenvector_boys


!******
!****s* synchronize_mpi/sync_vector_trm
!  NAME
!    sync_vector_trm
!  SYNOPSIS
      subroutine sync_vector_trm(is_worker, n_worker, blk_hess, nrow_hess, &
                    ncol_hess, vector)
      implicit none
!  ARGUMENTS
      logical, intent(in) :: is_worker
      integer, intent(in) :: n_worker
      integer, intent(in) :: blk_hess
      integer, intent(in) :: nrow_hess
      integer, intent(in) :: ncol_hess
      real*8, intent(inout) :: vector(ncol_hess)
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    July 2018.
!  SOURCE
!    local variables

      integer :: i_task
      integer :: send_count
      integer :: mpierr

      integer, allocatable :: recv_count(:)
      integer, allocatable :: recv_displ(:)
      real*8, allocatable :: send_array(:)

      if (.not. use_mpi) return

      allocate(recv_count(n_tasks))
      allocate(recv_displ(n_tasks))
      allocate(send_array(nrow_hess))

      recv_count = 0
      recv_displ = 0
      send_array = 0
      send_count = 0

      if (n_worker > 1) then
         recv_count(1:n_worker-1) = blk_hess
      end if

      recv_count(n_worker) = ncol_hess-(n_worker-1)*blk_hess

      recv_displ = (n_worker-1)*blk_hess

      do i_task = 1,n_worker
         recv_displ(i_task) = (i_task-1)*blk_hess
      end do

      if (is_worker) then
         send_count = nrow_hess
         send_array = vector(recv_displ(myid+1)+1:recv_displ(myid+1)+nrow_hess)
      end if

      call MPI_Allgatherv(send_array, send_count, MPI_DOUBLE_PRECISION, &
              vector, recv_count, recv_displ, MPI_DOUBLE_PRECISION, &
              mpi_comm_global, mpierr)

      deallocate(recv_count)
      deallocate(recv_displ)
      deallocate(send_array)

      end subroutine sync_vector_trm
!******
!----------------------------------------------------------
!****s* synchronize_mpi/sync_plus_u_occupation_matrix
!  NAME
!    sync_plus_u_occupation_matrix
!  SYNOPSIS
      subroutine sync_plus_u_occupation_matrix(occ_plus_u, n_plusUorb, & 
                  dim_mi,dim_mj,n_spin)
!  PURPOSE
!    Synchronize plus_u occupation matrix in case of occupation matrix control
      implicit none
!  ARGUMENTS
      integer :: dim_mi, dim_mj, n_spin, n_plusUorb
      real*8 :: occ_plus_u(n_plusUorb,dim_mi,dim_mj,n_spin)
!  INPUTS
!    o n_plusUorb -- number of plus_u shells
!    o dim_mi -- localized subspace dimension
!    o dim_mj -- localized subspace dimension
!    o p_max  -- number of contracted basis funktion for each hubbard projector
!    o occ_plus_u -- plus_u occupation matrix for each atom
!  OUTPUT
!    the plus_u occupation matrix is synchronized over all tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2016).
!  SOURCE

      real*8 :: occ_plus_u_MPI(n_plusUorb, dim_mi, dim_mj, n_spin)
      integer :: mpierr


      if (.not.use_mpi) return
!     initialize arrays to zero, just in case
      if (myid.ne.0) then
         occ_plus_u = 0.0d0
      end if
      
      occ_plus_u_MPI = 0.0d0

      call MPI_ALLREDUCE(occ_plus_u,occ_plus_u_MPI, &
          n_plusUorb*dim_mi*dim_mj*n_spin, MPI_DOUBLE_PRECISION, MPI_SUM, &
          mpi_comm_global, mpierr)

      occ_plus_u = occ_plus_u_MPI

      end subroutine sync_plus_u_occupation_matrix
!******

      end module synchronize_mpi
