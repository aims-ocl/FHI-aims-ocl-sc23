!****s* FHI-aims/prepare_scf
!  NAME
!   prepare_scf
!  SYNOPSIS

    subroutine prepare_scf

!  PURPOSE
!  High-level wrapper around all setup of fixed data for the calculation:
!  * Basic integration grid specifications for all species
!  * Free-atom DFT solution for each species on a logarithmic grid
!  * Radial basis functions for each species, tabulated on log. grids
!  * Coefficients that prepare the computation of Ylm functions and
!    their derivatives from cartesian products x^l*y^m*z^n
!
!  USES

      use localorb_io
      use timing
      use dimensions
      use runtime_choices
      use grids
      use geometry
      use pbc_lists
      use free_atoms
      use species_data
      use cartesian_ylm
      use analytic_multipole_coefficients
      use mpi_utilities
      use constraint
      use physics, only: n_electrons
      use octree_routines
      use pseudodata
      use rel_x2c_mod
      use xc_library, only: get_xc_name
      use vdw_correction, only: mbd_calc

      implicit none

!  ARGUMENTS
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
!  SOURCE



      ! local variables

      character*100 :: info_str
      character*4 descript

      ! counters

      integer :: i_species
      integer :: i_fn

      integer :: exc_code
      character(50) :: exc_origin
      character(150) :: exc_msg

      ! begin work

      call localorb_info('',use_unit)
      call localorb_info("------------------------------------------------------------",use_unit,'(A)')
      call localorb_info("Preparing all fixed parts of the calculation.",use_unit,'(10X,A)' )
      call localorb_info("------------------------------------------------------------",use_unit,'(A)')

      call get_timestamps(time_prep, clock_time_prep)

!     Here we distribute the atoms and radial shells to different MPI-tasks for
!     the initialize_integrals routine. The argument '-1' implies that there is no
!     output in the routine.

      call get_machine_precision &
      ( safe_minimum )

      ! (Rundong) Fully-relativistic calculations only support closed shell and nutral systems for now.
      ! In this case, We only use get_free_atoms(). use_initial_rho=.F., spin_treatment=0, etc.
      ! Becareful with this.
      if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) call prepare_x2c()

!     Initialize all grids, per species:
!     * radial and angular integration grids
!     * logarithmic mesh for atomic wave functions
!     Initial integration grids are simply chosen according to the
!     angular maximum specified in control.in. Optimized integration grids
!     may be determined in initialize_integrals.

        call get_grids &
        ( out_grids )

!     Allocate all free atom data
      call allocate_free_atoms ( )
      ! Self-consistent DFT for spherical free atoms on logarithmic grids:
      ! Potential, wave functions, density.
      call get_free_atoms( )

      if (use_initial_rho) then
        ! spin-polarization or charged system, so we should not use free_rho and free_rho_deriv
        ! directly for (spin) density initialization
        ! initialize initial density
         call initialize_initial_rho()
      end if

      if (out_basis) then
        ! Test output of free atom data, to see whether we are where we think we are
        descript='free'
        call free_atoms_out &
        ( free_potential, free_rho, n_atomic, free_wave, &
          atomic_n, atomic_l, free_wave_eigenval, descript )
      end if

!     get all non-(free-atom)-derived basis functions: ionic, hydrogenic

!     Allocate all species-related basis function arrays
      call allocate_species_waves ( )

      call get_fixed_basis_fns ( )

!     add cutoff potential to free atom-like effective potentials
!     FIXME: Structurally, this belongs into get_fixed_basis_fns().
!            The routine is called separately here because the resulting pot'l can be used
!            as a basis for adaptive potential.

      call fixed_basis_potential ( )

!     FIXME: The following should go entirely into get_fixed_basis_fns - that's where it belongs.
!     FIXME: Would have to include fixed_basis_potential there, too!

!     copy free_eigenval to species_eigenval
      do i_species = 1, n_species, 1
        if (include_min_basis(i_species)) then
          do i_fn = 1, n_atomic(i_species), 1
            atomic_eigenval(i_species, i_fn) = free_wave_eigenval(i_species, i_fn)
          enddo
        end if
      enddo

      ! WPH/VB: FHI-aims by default does not use the basis functions coming out
      ! of get_free_atom as the minimal basis set (they can be constructed with different
      ! cutoff radii). FHI-aims recalculates them in the get_species_basis_fn() by
      ! integrating free_potential (which was obtained from get_free_atom) using dftseq
      ! While this works fine for LDAs and GGAs, for basis functions incorporating
      ! exact exchange, this poses difficulties, as aims assumes free_potential has
      ! no l-channel dependence.  While we should implement this in the future, for now we
      ! bypass the recalculation of the minimal basis entirely when needed and reuse the
      ! wavefunctions and derivative quantities (eigenvalues, kinetic term, and possibly wavefunction
      ! gradients) that came from the atomic solver. However, this is only justified if the
      ! original free-atom solution and the radial functions generated later use the same
      ! cutoff potential.
      if (all(include_min_basis) .and. &
           & use_atomic_solver_fns_as_minimal_basis) then
        call use_atomic_solver_basis_fns_as_species_basis_fns (  )
      else
        call get_species_basis_fns (  )
      end if

      if (out_basis) then
!       Test output of free atom data, to see whether we are where we think we are
!       Prevent this output in a very crude way if the relevant arrays are not allocated
        if (use_min_basis) then
          descript='base'
          call free_atoms_out &
          ( basis_potential, free_rho, n_atomic, atomic_wave, &
            atomic_n, atomic_l, atomic_eigenval, &
              descript )
!         alternatively, write out atomic wave functions and kinetic energy
!           by a subroutine more similar to the fixed-basis case
          do i_species = 1, n_species, 1
            call atomic_out &
            ( i_species, n_grid(i_species), n_atomic, &
              atomic_n, atomic_l, &
              r_grid(1,i_species), atomic_wave, atomic_eigenval, &
              basis_potential, descript, n_species, n_max_grid, &
              n_max_ind_fns &
            )
          enddo
        end if
      end if

!     For fully-relativistic (4c-DKS or X2C) calculations, allocate the basis index here.
      if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
        call initialize_basis_index_x2c()
       !if(flag_rel.eq.REL_x2c) call atom_x2c_pce()
        if(flag_rel.eq.REL_x2c) call atom_x2c_pce_sto()
        ! For 4c-dks/q4c, small comp. density was already saved in cubic spline
        ! in the get_free_atoms procedure.
      endif

!     Condense fixed basis to only the number of basis functions needed.
!     Repair all grids and indices for efficient integration.
!     * Copy and label atomic basis parts.
!     * Orthogonalize fixed basis parts here.
!     * store and label all relevant fixed basis parts.

!     Two versions, the first one does a potential reorganization of the
!     basis functions if the cutoff radii are chosen to vary according to the
!     basis function shape. The second, commented routine does not, and remains for testing purposes only.
!     And yes, I do think that keeping reminder notes like this in the source can be helpful.

!     Matti::
!     if auxiliary basis set is given we prepare it first and then the normal basis
!     this is because we have to fix some arrays in shrink_fixed_basis_phi_thresh
!     to constitute normal basis

!      if (use_basis_dependent_cutoff) then
         call shrink_fixed_basis_phi_thresh
!      else
!         call shrink_fixed_basis
!      end if

      if (use_ext_basis) then
      ! construct large basis set
        call shrink_fixed_large_basis
      else
      ! we prepare the matrices of auxiliary basis set from normal basis set
        call copy_fixed_basis_to_ext_arrays
      end if

      n_states_init = n_states

      ! Test onsite integration of all basis functions now, after the generation of
      ! all radial functions is complete.
      ! It must be ensured that the (sparse) radial grid in the 3D overlapping
      ! atom-centered integrations later gives the same result as the much
      ! denser, one-dimensional logarithmic grid used to generate the basis
      ! functions - or, indeed, that an analytic integration would give.
      ! This applies especially to Gaussian-type basis functions with very high
      ! exponents.
      call verify_onsite_integrals

      ! Only here do we know the extent of all basis functions around the atoms.
      ! If we use external multipoles for embedding, we must make sure that no
      ! multipoles without an integration grid are found inside the basis function
      ! radius of any atom ... or we will try to integrate the effect of a Coulomb
      ! singularity with a finite integration grid (not a good idea).
      !
      ! Integration grids on multipoles can be placed for instance by using empty sites.

      if (n_multipoles.gt.0) then
          call check_multipole_distances()
      end if

!         call prepare_embedding_grids()
      if (use_embedding_pp) then
         call get_timestamps(time_embedding, clock_time_embedding)
         call get_pseudodata()
!    for indexing reasons we prepare grids here anyway.
!     later on we still can set their partition function to zero
         if(add_embedding_grids) then
            call prepare_pseudocoregrids()
            call prepare_more_pseudocoregrids()
         endif
         call get_times(time_embedding, clock_time_embedding, &
          &              tot_time_embedding, tot_clock_time_embedding)
      end if


      ! initialize coefficients that are needed for an analytic
      ! continuation of the partitioned multipole Hartree potential
      call initialize_analytic_multipole_coefficients()

      ! Initialize coefficients for cartesian ylm's as needed for forces
      ! In principle, that treatment is redundant with the same treatment
      ! done for analytic multipole potentials, the infrastructure simply
      ! exists twice because both pieces were developed in parallel.
      call initialize_cartesian_ylm(l_wave_max)

      call allocate_constraint_projectors()

      call distribute_radial_tasks(-1)

      ! If required, prepare pieces for the partition functions used to divide
      ! the integrals and the Hartree potential into localized atom-centered pieces.

      call prepare_partition_tabs ()

      if (use_libmbd) then
          mbd_input%atom_types = species_element(species)
          if (mbd_beta > 0) then
              mbd_input%mbd_beta = mbd_beta
          else
              mbd_input%xc = get_xc_name(flag_xc)
          end if
          mbd_input%comm = mpi_comm_global
          mbd_input%calculate_forces = use_forces
          mbd_input%coords = coords
          if (mbd_input%parallel_mode == 'auto') then
              if (use_scalapack) then
                  mbd_input%parallel_mode = 'atoms'
              else if (use_mpi) then
                  mbd_input%parallel_mode = 'k_points'
              else
                  mbd_input%parallel_mode = 'none'
              end if
          end if
          if (n_periodic > 0) then
              mbd_input%lattice_vectors = lattice_vector
              if (mbd_input%k_grid(1) == -1) mbd_input%k_grid = n_k_points_xyz
          end if
          call mbd_calc%init(mbd_input)
          call mbd_calc%get_exception(exc_code, exc_origin, exc_msg)
          if (exc_code > 0) then
              call localorb_info('  Error in setting up MBD')
              call aims_stop(exc_msg, exc_origin)
          end if
      end if

      call get_times(time_prep, clock_time_prep)

      write(info_str,'(A)') ''
      call localorb_info( info_str )

      write(info_str,'(2X,A)') "Preparations completed."
      call localorb_info( info_str )
      write(info_str,'(2X,A,F10.3,A)') "max(cpu_time)          : ", time_prep, " s."
      call localorb_info( info_str )
      write(info_str,'(2X,A,F10.3,A)') "Wall clock time (cpu1) : ", clock_time_prep, " s."
      call localorb_info( info_str )

      write(info_str,'(A)') "------------------------------------------------------------"
      call localorb_info( info_str )


    end subroutine prepare_scf
!******
