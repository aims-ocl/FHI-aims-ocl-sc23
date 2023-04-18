!****h* FHI-aims/physics
!  NAME
!    physics
!  SYNOPSIS

    module physics


!  PURPOSE
!
!     THIS MODULE SHOULD ONLY BE USED IN MAIN. ALL VARIABLES DECLARED HERE ARE SO GENERAL THAT
!     INDIVIDUAL VARIABLES SHOULD ALL BE PASSED EXPLICITLY IN SUBROUTINE CALLS SO THAT AN
!     OUTSIDE "PHYSICIST" CAN FOLLOW THEIR PATH THROUGH THE CODE!
!
!     **** IN FACT, VIOLATING THIS RULE AND TAKING VARIABLES DIRECTLY CAN LEAD TO RARE CORNER_CASE BUGS.
!          IF A FUNDAMENTAL DIMENSION VARIABLE CHANGES THROUGHOUT THE CODE FOR ANY REASON, THE ARRAY
!          IN QUESTION WILL BE RESHAPED WITH SMALLER, INTERNAL DIMENSIONS BY SOME SUBROUTINES IN THE
!          FOLLOWING, BUT WILL BE READ WITH THE ORIGINAL, LARGER DIMENSIONS BY ANYONE WHO USES PHYSICS DIRECTLY.
!
!          THIS ISSUE IS PARTICULARLY REAL FOR ARRAYS CONTAINING n_states DIMENSIONS, BUT COULD BE
!          TRUE ELSEWHERE AS WELL.
!
!     (a) Variables which are allocatable as they should be, and which I would like to keep
!     in main.f, to be passed around as explicit subroutine arguments. These variables
!     should include all central physical quantities.

!     partition_tab : Tabulated integral partition functions for each atom, radial, and
!     angular integration grid point; THIS IS ALREADY MULTIPLIED BY THE
!     INTEGRATION WEIGHTS FOR EACH POINT
!
!                -->  partition_fn * w_radial * w_angular * 4 * pi
!
!     The following three quantities encode the local Kohn-Sham potential. Because this
!     potential is only ever needed in integrals, we never calculate the full potential.
!     Reason: For gradient functionals, the XC potential requires the Hessian of the
!     density matrix, which is expensive to calculate, and expensive to store. Instead,
!     we use partial integration in all matrix elements / total energy terms, and avoid the
!     second derivative entirely.
!     The price is that we must now split the potential term into at least two separate
!     terms which are physically not so transparent:

!     * hartree_potential : only the full electrostatic part of the potential - the XC potential
!       costs no time in LDA / GGA, and is evaluated on the fly where needed.

!     * The following three are only used to refine the chemical potential when
!       using the PEXSI solver in ELSI
!     * hartree_potential_save : copy of hartree_potential
!     * dv_hartree_min: minimum value of (hartree_potential-hartree_potential_save)
!     * dv_hartree_max: maximum value of (hartree_potential-hartree_potential_save)

!     * free_hartree_superpos : Superposition of all free-atom hartree potentials,
!                 tabulated on the full integration grid.
!     * free_rho_superpos : Superposition of all free-atom densities,
!                 tabulated on the full integration grid.
!     * core_shift : shift, around each atom, between the free-atom potential and the effective
!                potential in the solid. Can use e' = [ e - v(r) + v_free(r) ], with r -> 0 ;
!                eigenvalues e' can then be directly compared to free-atom eigenvalues.

!     * overlap_matrix: Overlap matrix of current basis function set
!     * hamiltonian: Hamiltonian matrix of current basis function set
!
!         Both hamiltonian and overlap_matrix are stored in a packed format depending
!         on packed_matrix_format:
!
!         o 'none' (default for clusters):
!           Standard lapack packed format for symmetric matrices:
!                H_{i,j} = hamiltonian(i + ((j-1)*j)/2).
!           Please note that for periodic systems, i and j are the indices of all basis
!           functions within half the BvK-cell (1 <= i <= j <= n_centers_basis_I).
!
!         o 'index' (default for periodics and for
!                                use_scalapack.and.use_density_matrx):
!           Packed format with only non-zero entries actually stored.
!                H(i_basis_row, i_cell_row, i_basis_col) = hamiltonian(i_index)
!           with
!                i_index_first = index_hamiltonian(1, i_cell_row, i_basis_row)
!                i_index_last  = index_hamiltonian(2, i_cell_row, i_basis_row)
!                i_index_first <= i_index <= i_index_last
!                i_basis_col = column_index_hamiltonian(i_index).
!           In the case that there is no significant entry within the
!           corresponding row, we have i_index_first=0, i_index_last=-1.
!
!           This is also true for use_local_index.  The difference is
!           that for use_local_index, each Hamiltonian matrix element
!           is stored on exactly one node, whereas otherwise the
!           complete Hamiltonian is stored on each node.
!
!           Please note that this is effectively the compressed sparse row
!           storage scheme http://en.wikipedia.org/wiki/Sparse_matrix
!           (val=hamiltonian/overlap, cod_ind=column_index_hamiltonian,
!           row_pointer~index_hamiltonian(1,:,:)), but using explicit end
!           points in index_hamiltonian(2,:,:).
!
!             Example with cell_index(1,:) == 0,
!                          cell_index(2,:) == - cell_index(3,:):
!                       j=1   j=2
!             i=1 R=1 (  1     2  )   ! index_hamiltonian(:,1,1) = [ 1, 2]
!             i=2 R=1 (  2*    0  )   ! index_hamiltonian(:,1,2) = [ 0,-1]
!             i=1 R=2 (  0     3  )   ! index_hamiltonian(:,2,1) = [ 3, 3]
!             i=2 R=2 (  0     4  )   ! index_hamiltonian(:,2,2) = [ 4, 4]
!             i=1 R=3 (  0     0  )   ! index_hamiltonian(:,3,1) = [ 0,-1]
!             i=2 R=3 (  3*    4' )   ! index_hamiltoniah(:,3,2) = [ 5, 5]
!             hamiltonian              = [1, 2, 3, 4, 4']  # values
!             column_index_hamiltonian = [1, 2, 2, 2, 2]  # corresponding j
!             ! 2* and 3* are not stored because of symmetry [only 'U' stored]
!             ! 4  and 4' are equal because of symmetry but both stored.
!
!           Looping then looks like this:
!
!           do i_cell_row = 1, n_cells_in_hamiltonian-1   ! yes, "-1".
!             do i_basis_row = 1, n_basis
!               i_index_first = index_hamiltonian(1, i_cell_row, i_basis_row)
!               i_index_last = index_hamiltonian(2, i_cell_row, i_basis_row)
!               do i_index = i_index_first, i_index_last
!                 i_basis_col = column_index_hamiltonian(i_index)
!                 ! Use:
!                 !    hamiltonian(i_index, i_spin)
!                 !    density_matrix_sparse(i_index)
!                 !    and i_basis_row, i_cell_row, i_basis_col
!                 ! or any combination of
!                 !    (i_basis_row, i_loc_cell_row), &
!                 !    & (i_basis_col, i_loc_cell_col),
!                 ! with
!                 !    i_cell_row == &
!                 !    & position_in_hamiltonian(i_loc_cell_row, i_loc_cell_col)
!               end do
!             end do
!           end do
!
!
!           Just a few remarks:
!            * I (JW) would have named
!               + index_hamiltonian        -> row2index
!               + column_index_hamiltonian -> index2col
!              because the former specifies a range of indices for a given row
!              consisting of i_basis_row and i_cell_row.  The latter gives
!              the column (i_basis_col, as i_cell_col=0) for a given index.
!            * It is kind of unconventional in fortran that the fast index
!              corresponds to the column instead of the row.
!
!
!     * KS_eigenvector : List of Kohn-Sham eigenfunctions (single-particle states)
!     * KS_eigenvalue : List of Kohn-Sham eigenvalues (single-particle energies)
!     * occ_numbers : Occupation numbers for Kohn-Sham eigenstates
!
!          Please note that KS_eigenvector (or KS_eigenvector_complex) is distributed
!          and thus of shape
!                   (n_basis, n_states, n_spin, n_k_points_task),
!          whereas KS_eigenvalue and occ_numbers are shared/copied and have shape
!                   (n_states, n_spin, n_k_points).
!
!          The responsibility for the k-points is distributed as follows:
!
!          if (use_scalapack) then
!
!             In the ScaLapack case (n_tasks >= n_k_points), each node only
!             cares about scalapack_wrapper:my_k_point, and works on it in
!             scalapack_wrapper:eigenvec.  Obviously, n_k_points_task=1.
!             if (collect_eigenvectors) then
!                The array KS_eigenvectors{_complex}(:,:,:,1:1) contains the
!                corresponding eigencoefficients.
!             else
!                It doesn't!
!             end if
!
!          else
!
!              In this case, loops over the k-points can be done as follows:
!
!              i_k = 0
!              do i_k_point = 1, n_k_points
!                 if (myid == modulo(i_k_point, n_tasks) .and. myid <= n_k_points) then
!                    i_k = i_k + 1
!                    ... occ_numbers(:,:, i_k_point) ... KS_eigenvector(:,:,:, i_k) ...
!                 end if
!              end do
!
!              This is how it is done in get_KS_orbitals_p0(), which is /the/
!              reference by definition.  Be aware that myid starts at 0 whereas
!              i_k_point starts at 1.
!
!              Together with the condition (myid <= n_k_points) [nota bene:
!              "<="] this leads to the following distribution:
!              if (n_tasks <= n_k_points) then
!                 Each task is responsible for all k-points with
!                 myid == modulo(i_k_point, n_tasks).
!              else
!                 Only tasks with (1 <= myid <= n_k_points) are responsible for any
!                 k-point, namely for the one with (myid == i_k_point).
!              end if
!
!     * KS_eigenvalue_soc_perturbed : List of spin-orbit coupling perturbed
!                      Kohn-Sham eigenvalues (single-particle energies)
!     * KS_eigenvector_soc_perturbed : List of spin-orbit coupling perturbed
!                      Kohn-Sham eigenfunctions (single-particle states)
!     * occ_numbers_soc : SOC-perturbed occupation numbers
!     * rho : total density
!     * n_electrons: Total number of electrons in the structure
!     * delta_v_hartree_part_at_zero: delta hartree monopole potential at the nucleus ; "zero limit".

!     * rho_change: change in electron density from one self-consistency cycle to the other to check
!                 for convergence
!     * diff_forces : change in forces between different scf cycles
!     * ev_sum:     sum of eigenvalues in current sc cycle
!     * ev_sum_shifted: sum of eigenvalues in current sc cycle, shifted to free-atom-like values using
!                 core level shifts
!     * previous_ev_sum: sum of eigenvalues in preceding sc cycle
!     * previous_pot_jump: the potential jump in the preceding sc cycle
!     * en_xc: exchange correlation energy with current density
!     * en_pot_xc: mean value of exchange correlation potential for calculation of total energy
!     * en_hf   : exact exchange energy for post-processing xc evaluation
!     * en_post_xc   : xc energy for post-processing xc evaluation
!     * en_ion_ion: electrostatic interaction energy of the nuclei

!     * pot_ion_embed: potential due to presence of embedded (fixed) multipole charges
!     * en_ion_embed: electrostatic interaction energy of the nuclei with embedded charges
!     * en_vdw: empirical van der Waals correction
!     * en_ll_vdw:     LL van der Waals density functional energy
!     * en_ll_vdw_err: LL van der Waals density functional energy error
!     * en_lda_c     : total correlation energy of lda (which will be used for ll_vdwdf energy correction)
!     * en_pbe_c     : total correlation energy of pbe (which will be used for ll_vdwdf energy correction)
!     * rho_e: charge density projected on even space grids
!     * rho_e_gradient: gradient of charge density projected on even space grids
!     * rho_e_2gradient: squares of the gradient of charge density projected on even space grids
!     * hartree_energy_free: hartree energy (only of electrons!!) of superposition of free atom densities
!       hartree_energy_free is the quantity given in Eq. (61) of the
!       FHI-aims CPC paper, Blum et al., Computer Physics Communications 180 (2009) 2175-2196
!       when evaluated for the superposition density of spherical free atoms.
!     * hartree_delta_energy: Difference between Hartree energy of superpos of free atom densities
!       and full Hartree energy
!     * hartree_multipole_correction: Linear errror estimate of the Hartree energy inaccuracy due to
!                             the multipole decomposition of the Hartree potential with finite max.
!                             l_hartree(i_species)
!     * pot_jump: The potential shift in periodic calculation with dipole correction. Calculated
!                hartree_potential_recip.f90 and kept here to allow forcing convergence in
!                scf_solver.f90
!     * total_energy: total energy of system in current sc cycle
!     * chemical_potential : The electronic Fermi level, as determined by the occupation numbers
!                          and prescribed charge
!     * entropy_correction : In case of finite electronic "smearing" (Fermi, Gaussian, Methfessel-Paxton),
!                          the free energy of the system must be corrected by an "entropy" term to
!                          extrapolate down to the T=0 total energy.
!
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!     Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!
!  SOURCE

!-------------------------------------------------------------------------------


    implicit none

      real*8, dimension(:), allocatable, target :: partition_tab
      real*8, dimension(:), allocatable :: hartree_partition_tab
      real*8, dimension(:), allocatable :: hartree_partition_deriv
      real*8, dimension(:), allocatable :: weight_tab
      real*8, dimension(:), allocatable :: hartree_weight_tab

      !----------shanghui add for gradient_partition--------------
      real*8, dimension(:,:,:), allocatable :: partition_deriv_delley
      !----------shanghui end add for gradient_partition----------


      real*8, dimension(:), allocatable, target :: hartree_potential
      real*8, dimension(:), allocatable :: hartree_potential_save

!      real*8, dimension(:), allocatable :: free_xc_pot
      real*8, dimension(:),     allocatable :: free_hartree_superpos
      real*8, dimension(:),     allocatable :: free_rho_superpos
      real*8, dimension(:,:),   allocatable :: free_rho_gradient_superpos
      !SR: laplacian of nfree
!       real*8, dimension(:), allocatable :: free_rho_laplace_superpos
      real*8, dimension(:,:),   allocatable, target :: rho
      real*8, dimension(:,:,:), allocatable, target :: rho_gradient
      real*8, dimension(:,:),   allocatable :: rho_pce ! (Rundong) picture change error correction
      real*8, dimension(:,:,:), allocatable :: rho_pce_gradient
      ! the kinetic density of the system \tau_{sigma} = \sum |grad . psi|^2
      real*8, dimension(:,:), allocatable, target :: kinetic_density
      real*8, dimension(:,:), allocatable, target :: new_kinetic_density ! AJL: This is needed for delta_kinetic_density

      real*8, dimension(:),     allocatable :: core_shift
      real*8 :: av_core_shift

      real*8, dimension(:),   allocatable :: overlap_matrix
      real*8, dimension(:,:), allocatable :: hamiltonian
      real*8,     dimension(:,:,:,:), allocatable :: KS_eigenvector
      complex*16, dimension(:,:,:,:), allocatable :: KS_eigenvector_complex
      real*8,     dimension(:,:,:),   allocatable :: KS_eigenvalue
      real*8,     dimension(:,:,:,:), allocatable :: KS_eigenvector_irk ! on an irreducible set of k points
      complex*16, dimension(:,:,:,:), allocatable :: KS_eigenvector_complex_irk

      complex*16, dimension(:,:,:,:), allocatable :: KS_eigenvector_soc_perturbed
      real*8, dimension(:,:,:), allocatable :: KS_eigenvalue_soc_perturbed
      real*8, dimension(:,:,:), allocatable :: occ_numbers_soc
      ! the following flag is used to remember whether we estimate that
      ! this system has a low (or perhaps no) gap
      logical :: estimate_low_gap

      real*8 :: finite_nuclear_radius ! Used for finite nuclear models for Coulombic potentials
                                      ! Currently only used in atom_sphere

      ! RRS-PBC scheme, igor
      !  rrs_pbc_KS_eigenvector : List of Kohn-Sham eigenvalues (single-particle energies)
      !  rrs_pbc_KS_eigenvalue : List of Kohn-Sham eigenvalues (single-particle energies)
      complex*16, dimension(:,:), allocatable   :: rrs_pbc_hamiltonian_k
      complex*16, dimension(:), allocatable     :: rrs_pbc_overlap_k
      complex*16, dimension(:,:,:), allocatable :: rrs_pbc_KS_eigenvector
      real*8, dimension(:,:,:), allocatable     :: rrs_pbc_KS_eigenvalue
      real*8, dimension(:,:), allocatable       :: rrs_pbc_band_info
      complex*16, dimension(:,:,:), allocatable :: rrs_pbc_band_vect



      !atom average electrostatic potential, TZ
      !real*8,dimension(:), allocatable :: elec_atomic
      real*8,dimension(:), allocatable :: elec_free_atom
      real*8,dimension(:), allocatable :: elec_delta_atom
      !real*8,dimension(:), allocatable :: elec_delta_atom_recip
      real*8,dimension(:), allocatable :: elec_tot

      ! Beware that in the periodic case, the occ_numbers are
      ! scaled by the k_weights before the calculation of the
      ! charge density.  Use check_occs() to check correct weighting.
      real*8, dimension(:,:,:),   allocatable :: occ_numbers
      real*8, dimension(:),       allocatable :: delta_v_hartree_part_at_zero
      real*8, dimension(:,:),     allocatable :: delta_v_hartree_deriv_l0_at_zero
      real*8, dimension(:, :),    allocatable :: multipole_moments
      ! Current outer radius of multipole partitioned charge density (see
      ! species_data.f90):
      real*8, dimension(:),       allocatable :: multipole_radius_sq
      real*8, dimension(:,:),     allocatable :: outer_potential_radius
      integer,dimension(:),       allocatable :: l_hartree_max_far_distance

      real*8, dimension(:,:), allocatable :: ionic_forces
      real*8, dimension(:,:), allocatable :: hellman_feynman_forces
      real*8, dimension(:,:), allocatable :: pulay_forces
      real*8, dimension(:,:), allocatable :: multipole_forces
      real*8, dimension(:,:), allocatable :: gga_forces
      real*8, dimension(:,:), allocatable :: nlcc_forces
      real*8, dimension(:,:), allocatable :: pseudocore_forces
      real*8, dimension(:,:), allocatable :: total_forces
      real*8, dimension(:,:), allocatable :: previous_forces
      real*8, dimension(3, 3)             :: previous_stress

      ! VA: moved flag num_stress_finished from predict_new_geometry.f90  here,
      !     because need access to it in output_energy_and_forces.f90
      logical :: num_stress_finished
      real*8 :: numerical_stress_tensor(3,3)
      real*8, dimension (:,:), allocatable :: energy_gradient_on_lattice
      real*8, dimension (:,:), allocatable :: forces_lv
      real*8 :: energy_deriv_tress_tensor(3,3)
      real*8 :: numerical_pressure
      !CC: moved to own module: real*8 :: analytical_stress_tensor(3,3) = 0.d0

      ! Forces on external charges
      real*8, dimension(:,:), allocatable :: ext_charge_forces

      real*8 :: n_electrons

      real*8, dimension(:),     allocatable :: rho_change

      real*8, dimension(:),     allocatable :: rho_multipole_old
      real*8, dimension(:,:,:), allocatable :: rho_radial_multipole_old

      real*8 :: diff_forces
      real*8 :: diff_stress
      real*8 :: ev_sum
      real*8 :: pot_jump
      real*8 :: ev_sum_shifted
      real*8 :: previous_ev_sum
      real*8 :: previous_pot_jump
      real*8 :: ev_sum_change
      real*8 :: en_xc
      real*8 :: en_hf
      real*8 :: en_pot_xc
      real*8 :: en_post_xc
      real*8 :: en_ion_ion

      ! parts of the total energy needed to compute the kinetic energy
      ! in the end
      ! o en_elec_free -- free atoms energy
      !   en_elec_free is the term written in the last line of Eq. (64) of
      !   Blum et al., Computer Physics Communications 180 (2009) 2175-2196
      !   when evaluated for the superposition density of spherical free atoms.

      real*8 :: en_elec_free
      real*8 :: en_elec_delta
      real*8 :: fock_energy

      ! Energies for external embedding potential:
      ! embedding energy of nuclei, and embedding energy of charge density
      real*8 :: en_ion_embed
      real*8 :: en_density_embed

      ! in case one needs the charges later
      real*8, dimension(:), allocatable   :: hirshfeld_charge

      ! van der Waals energy correction
      real*8 :: en_vdw
      ! quantity used in the self-consistent implementation
      real*8 :: en_pot_vdw

      ! LL van der Waals density functional
      real*8 :: en_ll_vdw
      real*8 :: en_ll_vdw_err
      real*8 :: en_lda_c
      real*8 :: en_pbe_c
      real*8 :: meta_gga_total_energy
      real*8 :: vdw_total_energy
      real*8, dimension(:,:,:), allocatable :: rho_e
      real*8, dimension(:,:,:,:), allocatable :: rho_e_gradient
      real*8, dimension(:,:,:), allocatable :: rho_e_2gradient

      real*8, dimension(:), allocatable :: pot_ion_embed
      real*8 :: hartree_energy_free
      real*8 :: hartree_delta_energy
      real*8 :: hartree_multipole_correction
      real*8 :: total_energy
      real*8 :: SOC_non_sc_total_energy
      real*8 :: post_scf_total_energy
      real*8 :: previous_total_energy
      real*8 :: chemical_potential
      real*8 :: chemical_potential_soc ! For perturbative SOC
      real*8 :: vacuum_level_upper
      real*8 :: vacuum_level_lower
      real*8, dimension(:),allocatable :: chemical_potential_spin
      real*8 :: entropy_correction
      real*8, allocatable :: current_rho_multipole_spl(:,:,:)
      real*8, allocatable :: stored_rho_multipole_spl(:,:,:,:)
      real*8 :: dftpt2_total_energy_dft

      logical forces_on
      ! pulay_forces_on takes the role of gga_forces_on for the density-matrix based
      !   density update, including periodic systems
      logical pulay_forces_on
      logical gga_forces_on
      logical meta_gga_forces_on
      logical nlcc_forces_on
      logical AS_stress_on
      logical AS_pulay_stress_on

      ! External applied pressure
      real*8 :: external_pressure ! unit is eV/A**3

      ! External work due to external forces
      real*8 :: external_work ! unit is H

      logical :: grid_storage_allocated

      ! for BSSE calculations
      real*8 :: BSSE_full_energy_RPA, BSSE_full_energy_RPA_SE, BSSE_full_energy
      real*8, dimension(:), allocatable ::  BSSE_atom_RPA, BSSE_atom_RPA_SE,BSSE_per_atom
      real*8, dimension(:),allocatable :: multipole_radius_sq_public
      ! fo-dft
      real*8, dimension(:), allocatable :: local_fo_potential

      ! For PEXSI chemical potential search
      real*8 :: dv_hartree_min, dv_hartree_max

      contains

      subroutine set_physics_defaults ( )
         implicit none

         av_core_shift = 0.d0
         estimate_low_gap = .false.
         finite_nuclear_radius = 5.0d-7

         num_stress_finished = .false.
         numerical_stress_tensor = 0.0d0
         energy_deriv_tress_tensor = 0.0d0

         fock_energy = 0.0d0

         SOC_non_sc_total_energy = 0.0d0

         external_pressure = 0.0d0
         external_work = 0.0d0
         grid_storage_allocated = .false.

      end subroutine set_physics_defaults

!******
!-------------------------------------------------------------------------------
!****s* physics/allocate_physics
!  NAME
!   allocate_physics
!  SYNOPSIS
!
        subroutine allocate_physics()

!  PURPOSE
!    Allocates the variables of module physics which uses the number of atoms
!    as a dimension.
!  USES

      use dimensions
      use runtime_choices
      use mpi_tasks
      implicit none

!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



      integer:: info

      if (.not.allocated(multipole_radius_sq_public)) allocate(multipole_radius_sq_public(n_atoms))

      if (.not.allocated(outer_potential_radius))then
         allocate (outer_potential_radius(0:l_pot_max, n_atoms), stat=info)
         call check_allocation(info, 'outer_potential_radius        ')
      end if

      if (.not.allocated(rho_e))then
         allocate (rho_e(cell_edge_steps(1),cell_edge_steps(2),cell_edge_steps(3)), stat=info)
         call check_allocation(info, 'rho_e                         ')
      end if
      if (.not.allocated(rho_e_gradient))then
         allocate (rho_e_gradient(3,cell_edge_steps(1),cell_edge_steps(2),cell_edge_steps(3)), stat=info)
         call check_allocation(info, 'rho_e_gradient                ')
      end if
      if (.not.allocated(rho_e_2gradient))then
         allocate (rho_e_2gradient(cell_edge_steps(1),cell_edge_steps(2),cell_edge_steps(3)), stat=info)
         call check_allocation(info, 'rho_e_2gradient                ')
      end if

      if (.not.allocated(multipole_radius_sq))then
         allocate (multipole_radius_sq(n_atoms), stat=info)
         call check_allocation(info, 'multipole_radius_sq           ')
      end if

      if (.not.allocated(l_hartree_max_far_distance))then
         allocate (l_hartree_max_far_distance(n_atoms),stat=info)
         call check_allocation(info, 'l_hartree_max_far_distance    ')
       end if

      if (.not.allocated(delta_v_hartree_part_at_zero)) then
         allocate(delta_v_hartree_part_at_zero(n_atoms),stat=info)
         call check_allocation(info, 'delta_v_hartree_part_at_zero  ')
     end if

      if (.not.allocated(delta_v_hartree_deriv_l0_at_zero)) then
         allocate(delta_v_hartree_deriv_l0_at_zero(3,n_atoms),stat=info)
         call check_allocation(info, 'delta_v_hartree_deriv_l0_at_ze')
      end if
!  allocations for output atom average electrostatic potential
      if (flag_out_locpot_atom) then
         if (.not.allocated(elec_free_atom)) then
            allocate(elec_free_atom(n_atoms),stat=info)
            call check_allocation(info, 'elec_free_atom')
         end if
         if (.not.allocated(elec_delta_atom)) then
            allocate(elec_delta_atom(n_atoms),stat=info)
            call check_allocation(info, 'elec_delta_atom')
         end if
        !if (.not.allocated(elec_delta_atom_recip)) then
        !    allocate(elec_delta_atom_recip(n_atoms),stat=info)
        !    call check_allocation(info, 'elec_delta_atom_recip')
        ! end if

         if (.not.allocated(elec_tot)) then
            allocate(elec_tot(n_atoms),stat=info)
            call check_allocation(info, 'elec_tot')
         end if
      !   if (.not.allocated(elec_atomic)) then
      !      allocate(elec_atomic(n_atoms),stat=info)
      !      call check_allocation(info, 'elec_atomic')
      !   end if

         elec_tot = 0.0d0
      endif

!  allocations for force-calculation

      if (use_forces) then
         if (.not.allocated(ionic_forces)) then
            allocate(ionic_forces(3, n_atoms),stat=info)
            call check_allocation(info, 'ionic_forces                  ')
            ionic_forces = 0.d0
         end if
         if (.not.allocated(hellman_feynman_forces)) then
            allocate(hellman_feynman_forces(3, n_atoms),stat=info)
            call check_allocation(info, 'hellman_feynman_forces        ')
            hellman_feynman_forces = 0.d0
         end if
         if (.not.allocated(pulay_forces)) then
            allocate(pulay_forces(3, n_atoms),stat=info)
            call check_allocation(info, 'pulay_forces                  ')
            pulay_forces = 0.d0
         end if
         if (.not.allocated(multipole_forces)) then
            allocate(multipole_forces(3, n_atoms),stat=info)
            call check_allocation(info, 'multipole_forces              ')
            multipole_forces = 0.d0
         end if
         if (.not.allocated(pseudocore_forces)) then
            allocate(pseudocore_forces(3, n_atoms),stat=info)
            call check_allocation(info, 'pseudocore_forces              ')
            pseudocore_forces = 0.d0
         end if
         if (.not.allocated(total_forces)) then
            allocate(total_forces(3, n_atoms),stat=info)
            call check_allocation(info, 'total_forces                  ')
            total_forces = 0.d0
         end if
         if (.not.allocated(previous_forces)) then
            allocate(previous_forces(3, n_atoms),stat=info)
            call check_allocation(info, 'previous_forces               ')
            previous_forces = 0.d0
         end if
         previous_stress = 0.0d0
         if (.not.allocated(ext_charge_forces)) then
            if (use_qmmm) then
              allocate(ext_charge_forces(3, n_multipoles),stat=info)
              call check_allocation(info, 'ext_charge_forces             ')
              ext_charge_forces = 0.d0
            else
              allocate(ext_charge_forces(1,1)) ! dummy allocation
            end if
         end if
         if (use_gga) then
            if (.not.allocated(gga_forces)) then
               allocate(gga_forces(3, n_atoms),stat=info)
               call check_allocation(info, 'gga_forces                    ')
               gga_forces = 0.d0
            end if
         else
            if (.not.allocated(gga_forces)) then
               allocate(gga_forces(1, 1))
            end if
         end if

!         if (use_nonlinear_core) then
           if (.not.allocated(nlcc_forces)) then
              allocate(nlcc_forces(3, n_atoms),stat=info)
              call check_allocation(info, 'nlcc_forces              ')
              nlcc_forces = 0.d0
           end if
!         endif

      else
         if (.not.allocated(ionic_forces)) then
            allocate(ionic_forces(1, 1))
         end if
         if (.not.allocated(hellman_feynman_forces)) then
            allocate(hellman_feynman_forces(1, 1))
         end if
         if (.not.allocated(pulay_forces)) then
            allocate(pulay_forces(1, 1))
         end if
         if (.not.allocated(multipole_forces)) then
            allocate(multipole_forces(1, 1))
         end if
         if (.not.allocated(total_forces)) then
            allocate(total_forces(1, 1))
         end if
         if (.not.allocated(previous_forces)) then
            allocate(previous_forces(1, 1))
         end if
         if (.not.allocated(gga_forces)) then
            allocate(gga_forces(1, 1))
         end if
         if (.not.allocated(pseudocore_forces)) then
            allocate(pseudocore_forces(1, 1))
         end if
         if (.not.allocated(nlcc_forces)) then
            allocate(nlcc_forces(1, 1))
         end if
         if (.not.allocated(ext_charge_forces)) then
            allocate(ext_charge_forces(1, 1))
            ext_charge_forces = 0.d0
         end if
      end if

!  allocations for stress tensor calculation
      if (n_periodic .ne. 0) then
         if (.not.allocated(energy_gradient_on_lattice)) then
            allocate(energy_gradient_on_lattice(3,3))
            energy_gradient_on_lattice =0d0
         endif
      else
      ! Dummy variables for relaxation
         if (.not.allocated(energy_gradient_on_lattice)) then
            allocate(energy_gradient_on_lattice(1,1))
            energy_gradient_on_lattice =0d0
         endif
      end if

!  Must allocate it in ether case, because this variable (dummy or not)
!  goes into structure optimization.
!  later it is probably more convenient to directly pass the  (dummy) stress tensor
!  to the optimization routine.

      if (n_periodic .ne. 0) then
         if (.not.allocated(forces_lv)) then
            allocate(forces_lv(3,n_periodic))
            forces_lv =0d0
         endif
      else
      ! Dummy variables for relaxation
         if (.not.allocated(forces_lv)) then
            allocate(forces_lv(1,1))
            forces_lv =0d0
         endif
      end if

   ! allocations for atom_bsse
    if (calculate_atom_bsse) then
     if ( use_rpa_ene) then
      if (.not.allocated(BSSE_atom_RPA_SE)) then
         allocate(BSSE_atom_RPA_SE(n_atoms),stat=info)
         call check_allocation(info, 'BSSE_atom_RPA+SE')
      endif
      if (.not.allocated(BSSE_atom_RPA)) then
         allocate(BSSE_atom_RPA(n_atoms),stat=info)
         call check_allocation(info, 'BSSE_atom_RPA')
      endif
     elseif (use_hartree_fock) then
      if (.not.allocated(BSSE_per_atom)) then
         allocate(BSSE_per_atom(n_atoms),stat=info)
         call check_allocation(info, 'BSSE_per_atom')
      endif
     endif
    endif

   ! in order to save hirshfeld charges globally
      end subroutine allocate_physics
!******
!-------------------------------------------------------------------------------
!****s* physics/allocate_grid_storage
!  NAME
!    allocate_grid_storage
!  SYNOPSIS

      subroutine allocate_grid_storage()

!  PURPOSE
!    Allocates the grid storage variables of module physics.
!    This means variables which uses total number of gridpoint
!    as a dimension.
!  USES

      use dimensions
      use runtime_choices
      use mpi_tasks
      implicit none

!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      integer:: info


      if (grid_storage_allocated) return

!     allocate quantities for integrals in s-c loop

      if (.not.allocated(partition_tab)) then
         allocate (partition_tab(n_full_points),stat=info)
         call check_allocation(info, 'partition_tab                 ')
      end if

      if (.not.allocated(weight_tab)) then
         allocate (weight_tab(n_full_points),stat=info)
         call check_allocation(info, 'weight_tab                 ')
      end if

      if (.not.allocated(hartree_weight_tab)) then
         allocate (hartree_weight_tab(n_full_points),stat=info)
         call check_allocation(info, 'hartree_weight_tab                 ')
      end if

      if (.not.allocated(hartree_partition_tab)) then
         allocate( hartree_partition_tab(n_full_points),stat=info)
         call check_allocation(info, 'hartree_partition_tab         ')
      end if

      if (.not.allocated(rho)) then
         allocate (rho(n_spin,n_full_points),stat=info)
         call check_allocation(info, 'rho                           ')
      end if

      if (.not.allocated(rho_pce).and.(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)) then
         allocate (rho_pce(n_spin,n_full_points),stat=info)
         call check_allocation(info, 'rho_pce                       ')
      end if

      if (.not.allocated(rho_pce_gradient).and.(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)) then
         allocate (rho_pce_gradient(3,n_spin,n_full_points),stat=info)
         call check_allocation(info, 'rho_pce_gradient              ')
      end if

      if (.not.allocated(rho_multipole_old)) then
         if(flag_delta_rho_in_multipole) then
           allocate (rho_multipole_old(n_full_points),stat=info)
           call check_allocation(info, 'rho_multipole_old             ')
           rho_multipole_old = 0.d0
         else
           allocate (rho_multipole_old(1)) ! dummy allocation
         endif
      end if

      if (.not.allocated(rho_radial_multipole_old)) then
         if(flag_delta_rho_in_multipole)then
           allocate (rho_radial_multipole_old((l_pot_max+1)**2, n_max_radial+2, n_atoms),stat=info)
           call check_allocation(info, 'rho_radial_multipole_old      ')
           rho_radial_multipole_old = 0.d0
         else
           allocate (rho_radial_multipole_old(1,1,1)) ! dummy allocation
         end if
      end if


      if (.not.allocated(hartree_potential)) then
         allocate ( hartree_potential(n_full_points), stat=info)
         call check_allocation(info, 'hartree_potential             ')
      end if
      if (.not.allocated(free_hartree_superpos)) then
         allocate (free_hartree_superpos(n_full_points),stat=info)
         call check_allocation(info, 'free_hartree_superpos         ')
      end if
      if (.not.allocated(free_rho_superpos)) then
         allocate ( free_rho_superpos(n_full_points),stat=info)
         call check_allocation(info, 'free_rho_superpos             ')
      end if

      if (use_density_gradient) then
         if (.not.allocated(rho_gradient)) then
            allocate ( rho_gradient(3, n_spin, n_full_points),stat=info )
            call check_allocation(info, 'rho_gradient                  ')
         end if
      else
!     dummy allocation; this array should never be used
         allocate ( rho_gradient(1, 1, 1) )
      end if

      if (use_meta_gga .or. use_meta_gga_post .or. use_meta_gga_printout .or. use_libmbd) then
         if (.not.allocated(kinetic_density)) then
            allocate(kinetic_density(n_spin, n_full_points),stat=info)
            call check_allocation(info, 'kinetic_density               ')
         end if
         if (use_density_matrix) then
           if (.not.allocated(new_kinetic_density)) then
              allocate(new_kinetic_density(n_spin, n_full_points),stat=info)
              call check_allocation(info, 'new_kinetic_density               ')
           end if
         else
           allocate ( new_kinetic_density(1, 1) )
         end if
      else
!     dummy allocation; these arrays should never be used
         allocate ( kinetic_density(1, 1) )
         allocate ( new_kinetic_density(1, 1) )
      end if

      if ((multipole_interpolation_style.eq.1) .or. (use_density_gradient)) then
         if (.not.allocated(free_rho_gradient_superpos)) then
            allocate ( free_rho_gradient_superpos(3, n_full_points),stat=info )
            call check_allocation(info, 'free_rho_gradient_superpos    ')
         end if
!SR: laplacian of nfree
!          if (.not.allocated(free_rho_laplace_superpos)) then
!             allocate ( free_rho_laplace_superpos(n_full_points),stat=info )
!             call check_allocation(info, 'free_rho_laplace_superpos    ')
!          end if
      else
!     dummy allocation
         if (.not.allocated(free_rho_gradient_superpos)) then
            allocate ( free_rho_gradient_superpos(1, 1) )
         end if
!SR: laplacian of nfree
!          if (.not.allocated(free_rho_laplace_superpos)) then
!             allocate ( free_rho_laplace_superpos(1),stat=info )
!          end if
      end if

      if (multipole_interpolation_style.eq.1) then
         if (.not.allocated(hartree_partition_deriv)) then
            allocate ( hartree_partition_deriv(n_full_points),stat=info )
            call check_allocation(info, 'hartree_partition_deriv       ')
         end if
      else
         ! dummy allocation
         if (.not.allocated(hartree_partition_deriv)) then
            allocate ( hartree_partition_deriv(1) )
         end if
      end if


    !----------shanghui add for gradient_partition--------------
      if (use_partition_deriv) then
         if (.not.allocated(partition_deriv_delley)) then
            allocate ( partition_deriv_delley(3,n_atoms,n_full_points),stat=info )
            call check_allocation(info, 'partition_deriv_delley       ')
         end if
      else
            ! dummy allocation to allow subroutine calls without compiler warnings
            allocate ( partition_deriv_delley(1,1,1),stat=info )
      end if
    !----------shanghui end add for gradient_partition--------------



!     Allocations for embedding
!     Nadia
      if (.not.allocated(pot_ion_embed)) then
        if (use_embedding_potential) then
          allocate (pot_ion_embed(n_full_points),stat=info)
          call check_allocation(info, 'pot_ion_embed                 ')
        else
          ! dummy allocation
          allocate (pot_ion_embed(1))
        end if
      end if

      grid_storage_allocated = .true.

      end subroutine allocate_grid_storage
!******
!-------------------------------------------------------------------------------
!****s* physics/allocate_matrices
!  NAME
!    allocate_matrices
!  SYNOPSIS

      subroutine allocate_matrices()

!  PURPOSE
!    Allocates the number of the basis and the states depend variables
!    in the module physics.
!  USES

        use dimensions
        use runtime_choices
        use mpi_tasks
        use aims_memory_tracking, only : aims_allocate
        use rel_x2c_mod, only : dim_matrix_rel
        implicit none

!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      integer:: info


!     allocate overlap matrix and hamiltonian matrix

      if(.not.use_local_index .or. .not.use_scalapack) then
         if (.not.allocated(overlap_matrix)) then
            call aims_allocate( overlap_matrix, n_hamiltonian_matrix_size,   "overlap_matrix" )
         end if
         if (.not.allocated(hamiltonian)) then
            call aims_allocate( hamiltonian, n_hamiltonian_matrix_size, n_spin, "hamiltonian" )
         end if
      else
         ! Allocate dummies since otherways compiling with checking doesn't work
         if (.not.allocated(overlap_matrix)) call aims_allocate( overlap_matrix, 1, "overlap_matrix" )
         if (.not.allocated(hamiltonian))    call aims_allocate( hamiltonian, 1,1,     "hamiltonian" )
      endif


      if (.not.allocated(KS_eigenvalue)) then
        allocate( KS_eigenvalue(n_states,n_spin,n_k_points),stat=info )
        call check_allocation(info, 'KS_eigenvalue                 ')
      end if

     ! We do not allocate the SOC-perturbed eigenspectrum matrices until
     ! they are needed (during the second-variational post-processed SOC
     ! step.)  This is done to reduce memory usage during the SCF cycle.
     if (.not.calculate_perturbative_soc) then
       if (.not.allocated(KS_eigenvector_soc_perturbed))&
            call aims_allocate(KS_eigenvector_soc_perturbed,1,1,1,1,"KS_eigenvector_soc_perturbed")
       if (.not.allocated(KS_eigenvalue_soc_perturbed))&
            call aims_allocate(KS_eigenvalue_soc_perturbed,1,1,1,"KS_eigenvalue_soc_perturbed")
       if (.not.allocated(occ_numbers_soc)) then
            call aims_allocate(occ_numbers_soc,1,1,1,"occ_numbers_soc")
       end if
     end if

      if(collect_eigenvectors)then
      if(real_eigenvectors)then
         if (.not.allocated(KS_eigenvector)) then
            call aims_allocate( KS_eigenvector, n_basis,n_states,n_spin,n_k_points_task, "+KS_eigenvector")
            KS_eigenvector = 0.d0
         end if
      else
         if (.not.allocated(KS_eigenvector_complex)) then
            if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
              call aims_allocate( KS_eigenvector_complex, 2*dim_matrix_rel,n_states,n_spin,n_k_points_task, "+KS_eigenvector_complex")
            else
              call aims_allocate( KS_eigenvector_complex, n_basis,n_states,n_spin,n_k_points_task, "+KS_eigenvector_complex")
            endif
            KS_eigenvector_complex = 0.d0
         end if
      end if
      end if

      if (.not.allocated(KS_eigenvector)) call aims_allocate(KS_eigenvector, 1,1,1,1, "KS_eigenvector") ! allocate dummies
      if (.not.allocated(KS_eigenvector_complex)) call aims_allocate(KS_eigenvector_complex, 1,1,1,1, "KS_eigenvector_complex")

      if (.not.allocated(chemical_potential_spin)) then
         allocate(chemical_potential_spin(n_spin),stat=info)
         call check_allocation(info, 'chemical_potential_spin       ')
      end if

      if (.not.allocated(occ_numbers)) then
        allocate (occ_numbers(n_states,n_spin,n_k_points),stat=info)
        call check_allocation(info, 'occ_numbers                   ')
      end if

      if (.not.allocated(flag_KS_k_points)) then
        allocate (flag_KS_k_points(n_k_points),stat=info)
        call check_allocation(info, 'flag_KS_k_points              ')
      end if
      flag_KS_k_points = flag_KS

      end subroutine allocate_matrices
!******
!---------------------------------------------------------------------------
!****s* physics/reshape_matrices
!  NAME
!    reshape_matrices
!  SYNOPSIS

      subroutine reshape_matrices

!  PURPOSE
!    After the small number in the hamiltonian and overlap_matrix have removed
!    n_hamiltonian matrix size has changed. This is when packed matrix format index is in use.
!    This subroutine deallocates and allocates matrixes again to be smaller size.
!  USE

      use dimensions
      use runtime_choices
      use mpi_tasks
      use aims_memory_tracking, only : aims_allocate, aims_deallocate
      implicit none

!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      integer:: info
      real*8, allocatable,dimension(:,:):: temp_ham_ovlp

      call aims_allocate( temp_ham_ovlp, n_hamiltonian_matrix_size,n_spin, "temp_ham_ovlp" )

      temp_ham_ovlp = hamiltonian(1:n_hamiltonian_matrix_size,1:n_spin)
      call aims_deallocate( hamiltonian,                                     "hamiltonian" )
      call aims_allocate( hamiltonian, n_hamiltonian_matrix_size, n_spin,    "hamiltonian" )
      hamiltonian = temp_ham_ovlp

      temp_ham_ovlp(1:n_hamiltonian_matrix_size,1) = overlap_matrix(1:n_hamiltonian_matrix_size)
      call aims_deallocate( overlap_matrix,                               "overlap_matrix" )
      call aims_allocate( overlap_matrix, n_hamiltonian_matrix_size,      "overlap_matrix" )
      overlap_matrix = temp_ham_ovlp(1:n_hamiltonian_matrix_size,1)


      call aims_deallocate( temp_ham_ovlp,                                 "temp_ham_ovlp" )

      end subroutine reshape_matrices
!******
!---------------------------------------------------------------------------
!****s* physics/reshape_matrices_bigger
!  NAME
!    reshape_matrices_bigger
!  SYNOPSIS

      subroutine reshape_matrices_bigger

!  PURPOSE
!    The subroutine dellocates and allocates hamiltonian and overlap_matrix again.
!    This works also when n_hamiltonian_matrix_size is larger than in previous
!    allocation.
!  USES

        use dimensions
        use runtime_choices
        use mpi_tasks
        use aims_memory_tracking, only : aims_allocate, aims_deallocate
        implicit none

!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      integer:: info

      call aims_deallocate( hamiltonian, "hamiltonian" )
      call aims_allocate( hamiltonian, n_hamiltonian_matrix_size, n_spin, "hamiltonian" )

      call aims_deallocate( overlap_matrix, "overlap_matrix" )
      call aims_allocate( overlap_matrix, n_hamiltonian_matrix_size, "overlap_matrix" )

      end subroutine reshape_matrices_bigger
!---------------------------------------------------------------------------
!****s* physics/reshape_matrices_n_states
!  NAME
!    reshape_matrices_n_states
!  SYNOPSIS

      subroutine reshape_matrices_n_states(new_size)

!  PURPOSE
!    Reallocate n_states related arrays if n_states changes.
!    (1) new_size < n_states: n_states gets reduced due to ill-conditioning.
!    (2) new_size > n_states: n_states should be reset for a new geometry step.
!  USE
      use dimensions, only: n_basis,n_states,n_spin,n_k_points,n_k_points_task
      use runtime_choices, only: real_eigenvectors,collect_eigenvectors
      use aims_memory_tracking, only: aims_allocate,aims_deallocate
      implicit none
!  INPUT
      integer, intent(in) :: new_size
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Created on August 2018.
!  SOURCE
      integer :: old_size
      integer :: copy_size
      real*8, allocatable :: tmp(:,:,:)
      real*8, allocatable :: tmp_evec(:,:,:,:)
      complex*16, allocatable :: tmp_evec_cmplx(:,:,:,:)

      ! Occupation numbers
      old_size = ubound(occ_numbers,1)

      if(old_size /= new_size) then
         allocate(tmp(new_size,n_spin,n_k_points))

         copy_size = min(old_size,new_size)
         tmp = 0.d0
         tmp(1:copy_size,:,:) = occ_numbers(1:copy_size,:,:)

         deallocate(occ_numbers)
         allocate(occ_numbers(new_size,n_spin,n_k_points))

         occ_numbers = tmp

         deallocate(tmp)
      end if

      ! KS eigenvalues
      old_size = ubound(KS_eigenvalue,1)

      if(old_size /= new_size) then
         allocate(tmp(new_size,n_spin,n_k_points))

         copy_size = min(old_size,new_size)
         tmp = 0.d0
         tmp(1:copy_size,:,:) = KS_eigenvalue(1:copy_size,:,:)

         deallocate(KS_eigenvalue)
         allocate(KS_eigenvalue(new_size,n_spin,n_k_points))

         KS_eigenvalue = tmp

         deallocate(tmp)
      end if

      ! KS eigenvectors
      if(collect_eigenvectors) then
         if(real_eigenvectors) then
            old_size = ubound(KS_eigenvector,2)

            if(old_size /= new_size) then
               call aims_allocate(tmp_evec,n_basis,new_size,n_spin,&
                    n_k_points_task,"tmp_evec")

               copy_size = min(old_size,new_size)
               tmp_evec = 0.d0
               tmp_evec(:,1:copy_size,:,:) = KS_eigenvector(:,1:copy_size,:,:)

               call aims_deallocate(KS_eigenvector,"KS_eigenvector")
               call aims_allocate(KS_eigenvector,n_basis,copy_size,n_spin,&
                    n_k_points_task,"KS_eigenvector")

               KS_eigenvector = tmp_evec

               call aims_deallocate(tmp_evec,"tmp_evec")
            end if
         else
            old_size = ubound(KS_eigenvector_complex,2)

            if(old_size /= new_size) then
               call aims_allocate(tmp_evec_cmplx,n_basis,new_size,n_spin,&
                    n_k_points_task,"tmp_evec_cmplx")

               copy_size = min(old_size,new_size)
               tmp_evec_cmplx = (0.d0,0.d0)
               tmp_evec_cmplx(:,1:copy_size,:,:) &
                  = KS_eigenvector_complex(:,1:copy_size,:,:)

               call aims_deallocate(KS_eigenvector_complex,&
                    "KS_eigenvector_complex")
               call aims_allocate(KS_eigenvector_complex,n_basis,copy_size,&
                    n_spin,n_k_points_task,"KS_eigenvector_complex")

               KS_eigenvector_complex = tmp_evec_cmplx

               call aims_deallocate(tmp_evec_cmplx,"tmp_evec_cmplx")
            end if
         end if
      end if

      if(new_size > n_states) then
         n_states = new_size
      end if

      end subroutine reshape_matrices_n_states
!******
!---------------------------------------------------------------------------
!****s* physics/allocate_multipole_moments
!  NAME
!    allocate_multipole_moments
!  SYNOPSIS

      subroutine allocate_multipole_moments

!  PURPOSE
!    The subroutine allocates multipole_moments variable.
!  USES

        use dimensions
        use mpi_tasks
        implicit none

!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


       integer:: info

       if (.not.allocated(multipole_moments)) then
         allocate( multipole_moments( ( l_pot_max + 1)**2, n_atoms),stat=info)
         call check_allocation(info, 'multipole_moments             ')
       end if

      end subroutine allocate_multipole_moments
!******
!-------------------------------------------------------------------------------
!****s* physics/deallocate_physics
!  NAME
!    deallocate_physics
!  SYNOPSIS

      subroutine deallocate_physics()

!  PURPOSE
!    The subroutine deallocates variables of module physics exept the grid storagge variables and multipole_moments.
!  USES

        use runtime_choices
        use aims_memory_tracking, only: aims_deallocate
        implicit none

!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




      if (allocated(occ_numbers)) then
         deallocate( occ_numbers )
      end if
      if (allocated(chemical_potential_spin)) then
         deallocate( chemical_potential_spin )
      end if

      if(allocated(delta_v_hartree_part_at_zero))then
         deallocate(delta_v_hartree_part_at_zero)
      end if
      if (allocated(delta_v_hartree_deriv_l0_at_zero)) then
         deallocate(delta_v_hartree_deriv_l0_at_zero)
      end if
      !TZ deallocate

      if (allocated(elec_free_atom)) deallocate(elec_free_atom)
      if (allocated(elec_delta_atom)) deallocate(elec_delta_atom)
     ! if (allocated(elec_delta_atom_recip)) deallocate(elec_delta_atom_recip)
      if (allocated(elec_tot)) deallocate(elec_tot)
      !if (allocated(elec_atomic)) deallocate(elec_atomic)

      if (allocated(multipole_radius_sq_public)) deallocate(multipole_radius_sq_public)

      if (allocated(ext_charge_forces)) deallocate(ext_charge_forces)

      if (allocated(local_fo_potential)) deallocate(local_fo_potential)

      if (allocated(pseudocore_forces)) deallocate(pseudocore_forces)

      if (allocated(nlcc_forces)) deallocate(nlcc_forces)

      if (allocated(overlap_matrix)) then
         call aims_deallocate( overlap_matrix, "overlap_matrix"  )
      end if

      if (allocated(hamiltonian)) then
         call aims_deallocate( hamiltonian, "hamiltonian" )
      end if


      if (allocated(KS_eigenvalue)) then
         deallocate( KS_eigenvalue )
      end if

      if (allocated(KS_eigenvector)) then
         call aims_deallocate( KS_eigenvector, &
              "KS_eigenvector" )
      end if

      if (allocated(KS_eigenvector_complex)) then
         call aims_deallocate( KS_eigenvector_complex, &
              "KS_eigenvector_complex" )
      end if

      if (allocated(KS_eigenvector_soc_perturbed)) then
        call aims_deallocate(KS_eigenvector_soc_perturbed, &
             "KS_eigenvector_soc_perturbed")
      end if
      if (allocated(KS_eigenvalue_soc_perturbed)) then
        call aims_deallocate( KS_eigenvalue_soc_perturbed, &
             "KS_eigenvalue_soc_perturbed")
      end if
      if (allocated(occ_numbers_soc)) then
        call aims_deallocate(occ_numbers_soc, "occ_numbers_soc")
      end if

      if (allocated(rho_change)) then
        deallocate(rho_change)
      end if

      if (allocated(rho_multipole_old)) then
        deallocate(rho_multipole_old)
      end if

      if (allocated(rho_radial_multipole_old)) then
        deallocate(rho_radial_multipole_old)
      end if

      if (allocated(multipole_radius_sq))then
         deallocate (multipole_radius_sq)
      end if

      if (allocated(outer_potential_radius))then
         deallocate (outer_potential_radius)
      end if

      if (allocated(l_hartree_max_far_distance))then
         deallocate (l_hartree_max_far_distance)
      end if

      if (allocated(delta_v_hartree_part_at_zero)) then
         deallocate(delta_v_hartree_part_at_zero)
      end if

      if (allocated(ionic_forces)) then
         deallocate(ionic_forces)
      end if

      if (allocated(hellman_feynman_forces)) then
         deallocate(hellman_feynman_forces)
      end if

      if (allocated(pulay_forces)) then
         deallocate(pulay_forces)
      end if

      if (allocated(multipole_forces)) then
         deallocate(multipole_forces)
      end if

      if (allocated(total_forces)) then
         deallocate(total_forces)
      end if

      if (allocated(previous_forces)) then
         deallocate(previous_forces)
      end if

      if (allocated(gga_forces)) then
         deallocate(gga_forces)
      end if

      if (allocated(energy_gradient_on_lattice)) then
         deallocate(energy_gradient_on_lattice)
      end if

      if (allocated(forces_lv)) then
         deallocate(forces_lv)
      end if

      if (allocated(flag_KS_k_points)) then
        deallocate (flag_KS_k_points)
      end if

      if (allocated(hirshfeld_charge)) then
         deallocate(hirshfeld_charge)
      end if

      end subroutine deallocate_physics
!******
!-------------------------------------------------------------------------------
!****s* physics/deallocate_grid_storage
!  NAME
!    deallocate_grid_storage
!  SYNOPSIS

      subroutine deallocate_grid_storage()

!  PURPOSE
!    Deallocates the grid storage variables of the module physics.
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      implicit none

!     deallocate everything

      if (allocated(partition_tab)) then
         deallocate( partition_tab )
      end if

      if (allocated(weight_tab)) then
         deallocate( weight_tab )
      end if

      if (allocated(hartree_weight_tab)) then
         deallocate( hartree_weight_tab )
      end if

      if (allocated(hartree_partition_tab)) then
         deallocate( hartree_partition_tab)
      end if

      if (allocated(hartree_partition_deriv)) then
          deallocate(hartree_partition_deriv)
      end if


      !----------shanghui add for gradient_partition--------------
      if (allocated(partition_deriv_delley)) then
        deallocate(partition_deriv_delley)
      end if
      !----------shanghui end add for gradient_partition--------------



      if (allocated(free_hartree_superpos)) then
         deallocate( free_hartree_superpos)
      end if

      if (allocated(free_rho_superpos)) then
         deallocate( free_rho_superpos)
      end if

      if (allocated(free_rho_gradient_superpos)) then
         deallocate ( free_rho_gradient_superpos )
      end if

      if (allocated(hartree_potential)) then
         deallocate( hartree_potential )
      end if


      if (allocated(rho)) then
         deallocate(rho)
      end if


      if (allocated(rho_pce)) then
         deallocate(rho_pce)
      end if

      if (allocated(rho_multipole_old)) then
         deallocate(rho_multipole_old)
      end if

      if (allocated(rho_gradient)) then
         deallocate ( rho_gradient )
      end if

      if (allocated(rho_pce_gradient)) then
         deallocate(rho_pce_gradient)
      end if

      if (allocated(kinetic_density)) then
         deallocate ( kinetic_density )
      end if
      if (allocated(new_kinetic_density)) then
         deallocate ( new_kinetic_density )
      end if

      if (allocated(rho_e)) then
         deallocate(rho_e)
      end if

      if (allocated(rho_e_gradient)) then
         deallocate(rho_e_gradient)
      end if

      if (allocated(rho_e_2gradient)) then
         deallocate(rho_e_2gradient)
      end if

      if (allocated(pot_ion_embed)) then
        deallocate(pot_ion_embed)
      end if

      grid_storage_allocated = .false.

      end subroutine deallocate_grid_storage
!******
!------------------------------------------------------------------------
!****s* physics/cleanup_multipole_moments
!  NAME
!     cleanup_multipole_moments
!  SYNOPSIS

      subroutine cleanup_multipole_moments( )

!  PURPOSE
!    Deallocates multipole_moments
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    implicit none

      if(allocated(multipole_moments))then
        deallocate( multipole_moments )
      end if

      end subroutine cleanup_multipole_moments
!******
!-------------------------------------------------------------------------------
!****s* physics/deallocate_vdw_splines
!  NAME
!     deallocate_vdw_splines
!  SYNOPSIS
      subroutine deallocate_vdw_splines()

!  PURPOSE
!    Deallocates splines used by vdW-DF
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        implicit none

        if (allocated(current_rho_multipole_spl)) deallocate(current_rho_multipole_spl)
        if (allocated(stored_rho_multipole_spl)) deallocate(stored_rho_multipole_spl)

      end subroutine deallocate_vdw_splines
!******
!-------------------------------------------------------------------------------
    end module physics
