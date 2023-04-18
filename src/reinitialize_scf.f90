!****s* FHI-aims/reinitialize_scf
!  NAME
!    reinitialize_scf
!  SYNOPSIS

    subroutine reinitialize_scf &
    ( converged )

!  PURPOSE
!  High-level wrapper around the steps needed to reinitialize
!  the scf cycle from given species-dependent grids and wave functions
!
!  USES

      use localorb_io
      use dimensions
      use basis !
      use timing
      use runtime_choices
      use physics
      use species_data
      use cartesian_ylm !
      use analytic_multipole_coefficients !
      use mixing
      use grids
      use mpi_utilities
      use pbc_lists
      use hartree_non_periodic_ewald
      use hartree_fock
      use plus_u
      use lapack_wrapper
      use scalapack_wrapper
      use separate_core_states !
      use KH_core_states !
      use hartree_potential_real_p0
      use wf_extrapolation
      use thermodynamic_integration
      use octree_routines, only: cube_ready
      use hartree_potential_storage, only: &
          initialize_hartree_potential_storage, cleanup_hartree_potential_storage
      use load_balancing
      use mbd_std_wrapper, only: &
          mbd_std_initialize, mbd_std_finalize, mbd_self_consistent
      use pseudodata
      use vdw_correction
      use lpb_solver_utilities, only: dielecfunc_from_density, &
          dielecfunc_gradient_from_density_gradient,alphafunc_from_density,&
          dynamic_cavity,allocate_lpbsolver_dummy,&
          allocate_lpbsolver,mpb_start_crit,start_at_scf_step,dielec_func_mpb,&
          dielec_func_gradient_mpb,&
          solve_lpbe_only,decrement_kind,&
          correct_dielectric_decrement,initialize_eps_s_func,dynamic_ions,deps_drho,&
          d2eps_drho2,dielec_func_grad2,d2eps2_drho2,deps2_drho,dielec_func2,evaluate_v_free_gradient,&
          mpb_converged, initialized_lpb_solver,mpb_solver_started, deallocate_lpbsolver, v_free_gradient_public,&
          mpb_forces_on,Gnonmf_forces_on
      use SPE_solver, only : prepare_SPE_solver
      use mpb_solver_utilities, only: deallocate_mpbsolver, allocate_mpbsolver, ln_cosh_function,&
          allocate_mpbsolver_dummy
      use elsi_wrapper, only: aims_elsi_reinit_scf
      use initialize_grid_storage, only: initialize_grid_storage_p1
      use aims_memory_tracking, only : aims_allocate, aims_deallocate

      use free_atoms, only: free_rho_spl,free_rho_spl,free_drho_dr_spl,free_pot_es_spl,free_pot_es_at_zero
      use rlsy_interface, only: rlsy_h,generate_rlsy_handle
      use generate_index_map, only: rl_generate_index_map,rl_collect_overlap,rl_collect_hamiltonian,rl_check_if_reinit_needed
      use batch_rlsy_interface, only: batch_irr, compute_batch_rlsy_interface

      implicit none

!  ARGUMENTS

      logical :: converged

!  INPUTS
!    none
!  OUTPUT
!    o converged -- initialized to .false.
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
!

      ! local variables

      integer :: task_distribution_method
      character*100 :: info_str

      real*8 :: penalty_energy
      integer :: n_full_points_old, max_full_points

      real*8, allocatable :: ham_ovlp_work(:)

      character*8  :: cdate
      character*10 :: ctime

      ! only for timing output
      character(*), parameter :: deffmt = '2X'

      ! counters

      integer :: i_spin, info

      integer :: i_point, i_batch

      ! for bsse: may need cleaning up later
      integer :: i_stat,i_sp,i_p
      logical :: info_t

      integer :: i_coords

      real*8, allocatable, dimension(:,:) :: free_rho_superpos_copy
      real*8, allocatable, dimension(:,:,:) :: free_rho_gradient_superpos_copy

      integer :: info2

      real*8 :: clock_ke_time, ke_time

      integer :: prev_n_kpoints

      integer :: rlsy_verbosity
      logical :: rlsy_need_reinit

      ! begin work

      cube_ready = .false.  !vdw cube not defined yet for this MD step. SAG

      if(output_level .eq. 'MD_light' .and. MD_stepcount .gt. 1) output_priority = 2
!     First scf iteration is started separately - if we need to adapt all grids,
!     or starting density / potential, etc etc, then this must be done here

      write(info_str,'(A)') ''
      call localorb_info ( info_str,use_unit,'(A)', OL_norm  )

      write(info_str,'(A)') "------------------------------------------------------------"
      call localorb_info( info_str,use_unit,'(A)', OL_norm )

      write(info_str,'(10X,A)') "Begin self-consistency loop: Re-initialization."
      call localorb_info( info_str,use_unit,'(A)', OL_high )

      write(info_str,'(A)') ''
      call localorb_info ( info_str,use_unit,'(A)', OL_high )

      call date_and_time(cdate, ctime)
      write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, ", Time     :  ", ctime
      call localorb_info ( info_str,use_unit,'(A)', OL_high )

      write(info_str,'(A)') "------------------------------------------------------------"
      call localorb_info( info_str,use_unit,'(A)', OL_high)
      call get_timestamps( time_sc_loop, clock_time_sc_loop)
      ! minor repetition - this counter is used for relaxation or MD, to keep track of the
      ! time used for one fully converged (single-point) scf cycle. Must be initialized here
      ! as well as in reinitialize_scf()
      call get_timestamps ( time_sc_tot, clock_time_sc_tot )

      ! This is the place where we restart an scf cycle from scratch for the first time.
      scf_restarted = .true.
      restarting_scf = .false. ! reset to default.

      call reprepare_pseudocoregrids()
      call distribute_radial_tasks(-1)

      if (use_embedding_pp) then
          call get_timestamps(time_embedding, clock_time_embedding)
          call baswave_pp_overlap3()

          call get_nonlocal_pot()
          call get_times(time_embedding, clock_time_embedding, &
          &              tot_time_embedding, tot_clock_time_embedding)
      end if

      ! AJL/ July 2016
      ! This is necessary to prevent some weird spurious numbers from kinetic_density
      ! ruining the XC energy. Seems to be a problem for stress optimisation on some machines,
      ! perhaps because the kinetic_density is being deallocated and then reallocated?
      ! Ideally this'd be tracked down to error source but this work around will do for now.
      if (use_meta_gga) first_integration = .true.

      ! steps for atom_bsse from prepare scf
      if (calculate_atom_bsse) then
         write(info_str,'(2X,A)') "Preparing basis functions for fresh BSSE geometry"
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
         call shrink_fixed_basis_phi_thresh
         call cleanup_analytic_multipole_coefficients
         call initialize_analytic_multipole_coefficients()

         call cleanup_cartesian_ylm()
         call initialize_cartesian_ylm(l_wave_max)
         call prepare_partition_tabs () ! this step is really important
      endif

      ! These are the conditions under which an Ewald summation is needed
      ! for the Hartree potential and the Ewald range separation radius
      ! might have changed
      if ( (n_periodic.gt.0) .and. unit_cell_shape_changed .and. &
           Ewald_radius_automatic ) &
         then
         ! in this case, update the real-space part of the Ewald sum
         ! to reflect the new Ewald_radius
         call cleanup_analytic_multipole_coefficients
         call initialize_analytic_multipole_coefficients
      end if

      if (use_molecular_dynamics.and.use_wf_extrapolation.and..not.wf_extra_use_densmat) then
        write(info_str,'(2X,A)') "Extrapolating wavefunction / Hamiltonian for scf reinitialization."
        call localorb_info ( info_str,use_unit,'(A)', OL_norm )

        call wf_extrapolate_eigenvectors(occ_numbers, KS_eigenvector, KS_eigenvector_complex)
      end if

      ! reinitialize lists of integration centers, Hartree potential
      ! multipole centers, etc, whose shape depends on the boundary
      ! conditions (periodic or not)
      call get_timestamps(time_bc_init, clock_time_bc_init)

      call deallocate_pbc_lists( )
      prev_n_kpoints = n_k_points
      call initialize_bc_dependent_lists( )
      if( packed_matrix_format /= PM_none ) then
        if(n_periodic .eq. 0) then
          call initialize_packed_matrix_cluster
        end if
      end if
      call get_times(time_bc_init, clock_time_bc_init, tot_time_bc_init, tot_clock_time_bc_init)

      ! Got k-points, reset symmetry arrays
      if (use_symmetry_reduced_spg)then
        call get_sym_overlap()
      endif
      ! Hartree potential real part have to be reallocated also, because the size
      ! of the atoms list can be different.
      call cleanup_hartree_potential_real

      ! If Ewald's decomposition of the Hartree potential is to be
      ! used in the non-periodic case, then perform the initialization.
      if (use_hartree_non_periodic_ewald) then
         if(output_level .eq. 'MD_light') output_priority = OL_high
         call initialize_hartree_non_periodic_ewald
         if(output_level .eq. 'MD_light') output_priority = OL_low
      end if

      if (calculate_atom_bsse) then
         !new
             call initialize_separate_core_states
             call initialize_KH_core_states
      endif

      if (use_plus_u) then
         call deallocate_plus_u
         call allocate_plus_u
      endif

      if(packed_matrix_format == PM_index .and. .not. use_local_index)then
         call reshape_matrices_bigger
      end if

      ! Pulay mixing must be reinitialized from scratch for next iteration
      call cleanup_pulay ( )
      call cleanup_broyden ( )

        ! the grid is known - preload all free-atom quantities which are
        ! needed in an scf. run


      n_full_points_old = n_full_points

      call get_timestamps(time_partition, clock_time_partition)



      if (recompute_batches_in_relaxation) then
         ! Reset the batches to repartition the entire integration grid.
         call reset_grid_partitions()
      else
         ! Just compute the new coordinates.
         call update_grid_coords()
      end if

      ! switch priority on a case by case basis.
      if(output_level .eq. 'MD_light') output_priority = OL_high
      call partition_grid()
      if(output_level .eq. 'MD_light') output_priority = OL_low


      if (use_embedding_pp) then
          call get_timestamps(time_embedding, clock_time_embedding)
          call put_localpot_on_integration_grid()
          if(use_nonlinear_core) then
             call partial_core_dens2grid_v2()
          endif
          call get_times(time_embedding, clock_time_embedding, &
          &              tot_time_embedding, tot_clock_time_embedding)
      end if


      call get_timestamps(rtime,clock_rtime)
      time_partition = rtime - time_partition
      tot_time_partition = tot_time_partition + time_partition
      call sync_timing(time_partition)
      clock_time_partition = clock_rtime - clock_time_partition
      tot_clock_time_partition = tot_clock_time_partition + clock_time_partition

      ! begin measurement of preparation of grid quantities here.
      call get_timestamps(time_free_offset, clock_time_free_offset)

      task_distribution_method = 0
      call distribute_tasks(task_distribution_method)

      if ((n_full_points_old.ne.n_full_points) .or. calculate_atom_bsse) then
         ! Size of all grid storage changed. Must deallocate and reallocate. Bleah.
         call deallocate_grid_storage( )
         call allocate_grid_storage( )
      end if

      grid_storage_initialized = .false.
      if(output_level .eq. 'MD_light') output_priority = OL_high
      call initialize_grid_storage_p1 &
           ( partition_tab, weight_tab, hartree_partition_tab, hartree_weight_tab, hartree_partition_deriv, &
             partition_deriv_delley, &
           hartree_potential, &
           free_hartree_superpos, free_rho_superpos, &
           free_rho_gradient_superpos, &
           rho, rho_gradient, &
           pot_ion_embed &
           )
      if(output_level .eq. 'MD_light') output_priority = OL_low

      if (compensate_multipole_errors.or.normalize_initial_density) then
         ! If we are interested in having an initial  density that has the exact
         ! charge as given by control.in and as integrated on the actual, 3D overlapping
         ! atom-centered integration grid, we must enforce this zero here.
         ! REMINDER: This 'if' statement is executed for *both* compensate_multipole_errors
         ! and normalize_initial_density. The functionality for 'normalize_initial_density'
         ! should be an exact subset of that for 'compensate_multipole_errors'.
         call renormalize_initial_density (rho, rho_gradient, n_spin, partition_tab, n_electrons)
      end if

      if (compensate_multipole_errors) then
         ! If we are interested in having a density that has zero charge as
         ! integrated on the actual integration grid, we must enforce this
         ! zero here.
         ! There are two arrays that we must renormalize:
         ! (1) the actual electron density, rho
         ! (2) the reference density for the Hartree potential, that is the
         !     sum of free atom densities - which is subtracted from rho before the Hartree
         !     potential is computed.
         !
         ! Reminder: If we repair the code to allow the use of a hartree_partition_table that is
         !     different from the normal partition_tab, we must renormalize here
         !     and in scf_solver
         !     using hartree_partition_tab instead AND make sure that the missing
         !     integration weights are added to both subroutines.
         ! extra call for free_rho_superpos, which has a factor 4 pi multiplied in.
         call renormalize_free_density (free_rho_superpos, free_rho_gradient_superpos, 1, &
                                   partition_tab, n_electrons)
         !
         ! What is not renormalized now is the electrostatic potential of the free atom
         ! superposition - because that is currently exactly consistent point by point with
         ! the nuclear potential.
         ! It is a good question whether this will become a problem later.
      end if

      ! update dimensions for workspaces in integrations
      ! switch output priority on a case by case basis
      if(output_level .eq. 'MD_light') output_priority = OL_high
      got_n_compute_maxes = .false.
      call get_n_compute_maxes_p1( partition_tab)
      if(output_level .eq. 'MD_light') output_priority = OL_low


      call get_timestamps(rtime, clock_rtime)
      time_free_offset =  rtime - time_free_offset
      tot_time_free_offset = tot_time_free_offset + time_free_offset
      call sync_timing(time_free_offset)
      clock_time_free_offset =  clock_rtime - clock_time_free_offset
      tot_clock_time_free_offset = tot_clock_time_free_offset + clock_time_free_offset

      call get_timestamps(time_superpos, clock_time_superpos)

      if (use_mbd_std .and. mbd_self_consistent) then
          call mbd_std_finalize()
          call mbd_std_initialize()
      end if


      call get_free_superpos_energy_p1 &
        ( partition_tab, &
          rho, rho_gradient, kinetic_density, pot_ion_embed, &
          hartree_energy_free, en_elec_free, en_xc, en_pot_xc, en_ion_ion, &
          en_ion_embed, en_density_embed, en_vdw, en_ll_vdw, en_ll_vdw_err, en_lda_c, &
          en_pbe_c, ionic_forces, ext_charge_forces &
        )


      call get_timestamps(rtime, clock_rtime)
      time_superpos = rtime-time_superpos
      tot_time_superpos = tot_time_superpos+time_superpos
      call sync_timing(time_superpos)
      clock_time_superpos = clock_rtime-clock_time_superpos
      tot_clock_time_superpos = tot_clock_time_superpos+clock_time_superpos

      ! Reset geometry dependend variables in scalapack

      if(use_scalapack) call reinitialize_scalapack()

      ! Reinitialize ELSI
      call aims_elsi_reinit_scf()

      ! When the partioning is known, we can set the index for the sparse local matrices
      ! init_local_index_communication must be called after reinitialize_scalapack()!
      if(use_local_index) then
         if(output_level .eq. 'MD_light') output_priority = OL_high
         call set_index_hamiltonian
         if(output_level .eq. 'MD_light') output_priority = OL_low
         call init_comm_sparse_local_matrix_scalapack
      endif

      ! Reinitialize storage for rho_multipole
      call cleanup_hartree_potential_storage()
      call initialize_hartree_potential_storage()

      ! Reset load balancing data completely
      if(use_load_balancing) call reset_load_balancing()

!     re-integrate the overlap matrix only, for a fixed (given) grid

      call get_timestamps(time_intg, clock_time_intg)

      if (.not. use_local_index .and. use_scalapack) then
         if (size(overlap_matrix) /= n_hamiltonian_matrix_size) then
            call aims_deallocate( overlap_matrix,                          "overlap_matrix" )
            call aims_allocate( overlap_matrix, n_hamiltonian_matrix_size, "overlap_matrix" )
         else if (.not.(allocated(overlap_matrix))) then
            call aims_allocate( overlap_matrix, n_hamiltonian_matrix_size, "overlap_matrix" )
         end if
      end if

      if(use_local_index) then

         ! Since the hamiltonian is normally not computed at all
         ! in reinitialize_scf.f90 we do the first calibration for
         ! load balancing with integrate_ovlp_matrix_p2.

         call aims_allocate( ham_ovlp_work, n_hamiltonian_matrix_size, "ham_ovlp_work" )

         if(use_load_balancing) get_batch_weights = .true.
         call integrate_ovlp_matrix_p2(partition_tab,l_shell_max,ham_ovlp_work)

         call set_sparse_local_ovlp_scalapack(ham_ovlp_work)

         if(.not. use_elsi_dm) then
            call normalize_eigenvectors_scalapack()
         endif

         get_batch_weights = .false.

         if(use_load_balancing) then
            call compute_balanced_batch_distribution(n_bp_integ)
            call init_comm_full_local_matrix(batch_perm(n_bp_integ)%n_basis_local, &
                                             batch_perm(n_bp_integ)%i_basis_local)
         endif

         call aims_deallocate( ham_ovlp_work, "ham_ovlp_work" )

      else

         call integrate_ovlp_matrix_p2(partition_tab,l_shell_max,overlap_matrix)

         if(use_lapack_fast .and. .not. use_elsi) then
            ! if this run uses the new lapack version, we must here tell it that a new overlap
            ! matrix exists. This is implicitly done by deallocating the auxiliary array that
            ! contained the Cholesky factors of the last step.
            call clear_lapack_overlap_matrix()
         endif

         if(use_scalapack) then

            if(packed_matrix_format /= PM_none) then
               call construct_overlap_scalapack(overlap_matrix)
            else
               call setup_overlap_scalapack(overlap_matrix)
            endif

            if(.not. use_elsi_dm) then
               call normalize_eigenvectors_scalapack()
            endif

         endif
      endif

      if(use_plus_u)then
         !only used for plus_u with mulliken charges
         call plus_u_init_idx
      endif

      ! OH: This is, as far as I have understood, a good place to initialize the
      ! symmetry helper. But I might very well be wrong and change my mind at any
      ! point in the near future.
      if ( use_rlsy_symmetry ) then
          call rl_check_if_reinit_needed(rlsy_h,frac_coords,lattice_vector,rlsy_need_reinit)

          if ( rlsy_need_reinit ) then
              rlsy_verbosity=1
              !skip_rlsy_symmetry=.true.
              call generate_rlsy_handle(rlsy_h,&
                  species,species_z,frac_coords,lattice_vector,&
                  basis_m,basis_l,basis_fn,basis_atom,outer_radius,&
                  r_grid_min,r_grid_inc,basis_wave_spl,basis_deriv_spl,outer_partition_radius,&
                  n_k_points_xyz,k_point_list,k_points_offset,&
                  atom_radius,free_rho_spl,free_drho_dr_spl,free_pot_es_spl,free_pot_es_at_zero,&
                  n_radial,n_angular,r_radial,r_angular,w_radial,w_angular,scale_radial,l_hartree,&
                  n_electrons,n_states,n_spin,&
                  elsi_solver,elsi_out_level,occupation_type,solver_method,&
                  elsi_elpa_n_single,occupation_width,occupation_acc,basis_threshold,&
                  rlsy_verbosity,use_unit,MPI_COMM_WORLD,skip_rlsy_symmetry)
              ! Generate the indexing so that I can access the Hamiltonian/overlap/other things
              !call rl_generate_index_map(rlsy_h,&
              !    coords_center, center_to_cell, center_to_atom, &
              !    index_hamiltonian, column_index_hamiltonian, cbasis_to_basis, &
              !    cbasis_to_center, centers_basis_integrals, frac_coords, lattice_vector,&
              !    position_in_hamiltonian)
              ! Collect the overlap matrix in case it is needed for something
              !call rl_collect_overlap(rlsy_h,overlap_matrix)
              !call rl_collect_hamiltonian(rlsy_h,hamiltonian)
              call compute_batch_rlsy_interface()
          endif
      endif

      ! If scalapack is used, the actual overlap matrix is now stored in scalapack format.
      ! Then, the original array used during the integration can be deallocated.
      if (.not. use_local_index .and. use_scalapack .and. packed_matrix_format /= PM_index ) then
         ! if packed_matrix_format index, then this deallocation only happens in scf_solver.f90 ...

         if(use_constraint .or. flag_run_mulliken .or. out_lowdin .or. out_band .or. &
              out_matrices .or. flag_KS==-2 .or. use_plus_u .or. use_out_eigenvec_ovlp .or. &
              transport_lead_calculation .or. transport_calculation .or. flag_run_tddft_real_time_propagation) then
            write(info_str,'(2X,A)') "Not deallocating overlap matrix."
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
         else
            write(info_str,'(2X,A)') "Deallocating overlap matrix."
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            call aims_deallocate( overlap_matrix,  "overlap_matrix" )
            call aims_allocate( overlap_matrix, 1, "overlap_matrix" ) ! make it defined
         endif
      endif

      ! Save the scalapack overlap in the case it is needed for post processing
      if(use_scalapack .and. (flag_run_mulliken .or. out_lowdin)) call save_overlap_scalapack

      call get_timestamps(rtime, clock_rtime)
      time_intg =  rtime - time_intg
      tot_time_intg = tot_time_intg + time_intg
      call sync_timing(time_intg)
      clock_time_intg =  clock_rtime - clock_time_intg
      tot_clock_time_intg = tot_clock_time_intg + clock_time_intg

      if (use_symmetry_reduced_spg.and.prev_n_kpoints.ne.n_k_points)then
        call reset_hamiltonian(converged)
      endif
 !--------------------------------------------------------
      if (calculate_atom_bsse) then
        write(info_str,'(2X,A)') "Recalculating the Hamiltonian for atomization BSSE correction"
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

      ! redo the Hamiltonian
      ! deallocate hamiltonian
        if (.not. use_local_index .and. use_scalapack) then
            call aims_deallocate( hamiltonian,                                 "hamiltonian" )
            call aims_allocate( hamiltonian, n_hamiltonian_matrix_size,n_spin, "hamiltonian" )
        end if

        ! Integrate the  Hamiltonian matrix
        first_integration = .true. ! new
        if(use_local_index) then

           ! TODO: Use as calibration for load balancing
           call aims_allocate( ham_ovlp_work, n_hamiltonian_matrix_size*n_spin, "ham_ovlp_work" )

           call integrate_real_hamiltonian_matrix_p2 &
                 ( hartree_potential,  &
                  rho, rho_gradient, kinetic_density, &
                  partition_tab, l_shell_max, &
                  en_xc, en_pot_xc, ham_ovlp_work, en_vdw, en_pot_vdw )

           call set_sparse_local_ham_scalapack(ham_ovlp_work)
           call aims_deallocate( ham_ovlp_work, "ham_ovlp_work" )
        else
            call integrate_real_hamiltonian_matrix_p2 &
                 ( hartree_potential,  &
                 rho, rho_gradient, kinetic_density, &
                 partition_tab, l_shell_max, &
                 en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )

        endif

        if(.not.use_local_index .and. packed_matrix_format == PM_index)then
          call remove_small_numbers_in_hamiltonian_and_ovlp(hamiltonian, overlap_matrix,info_t)
          if(info_t) call reshape_matrices
        end if

        if (.not.use_local_index .and. out_matrices) then
          if (out_matrices_format_2005) then
            call output_real_matrices &
                 ( overlap_matrix, hamiltonian )
          else
            call output_real_matrices_aij &
                 ( overlap_matrix, hamiltonian )
          end if
        end if

        if (.not.use_local_index .and. out_aitranss) then
          call output_ka_overlap(overlap_matrix)
        end if

        ! If use_scalapack is set, allocate and set the distributed overlap matrix now.
        ! Under certain circumstances we can deallocate the (non-distributed) overlap_matrix
        ! here and save some non-scalable memory
        if (.not.use_local_index .and. use_scalapack) then
          ! If scalapack is used, the hamiltonian has to be stored
          ! into the scalapack arrays

          if(n_periodic>0 .or. packed_matrix_format /= PM_none ) then
            call construct_hamiltonian_scalapack(hamiltonian)
          else
            call setup_hamiltonian_scalapack(hamiltonian)
          endif
        endif
        !     After initial integrations, solve eigenvalue problem for the first time ...
        !     ... UNLESS this is a restart run, where the initial eigenvectors are obtained from
        !     a separate file anyway



        !     we obtain the Kohn-Sham eigenvalues from a Kohn_Sham solution

        !  1st reset values : new
        KS_eigenvalue(:,:,:)=0.d0
        KS_eigenvector(:,:,:,:)=0.d0
        KS_eigenvector_complex(:,:,:,:)=0.d0
        occ_numbers(:,:,:) = 0.0d0
        chemical_potential=0.d0
        chemical_potential_spin(:)=0.d0 ! new

! use newer version of get_KS_orbitals_p0
        call advance_KS_solution(overlap_matrix, hamiltonian, n_electrons, &
             KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
             occ_numbers, chemical_potential, chemical_potential_spin)

        call evaluate_KH_core_states(   overlap_matrix, rho, &
             rho_gradient, kinetic_density, &
             partition_tab, hartree_potential ) !new

        call add_KH_states(KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex)!new

        call output_real_eigenfunctions &
             ( KS_eigenvalue, KS_eigenvector, occ_numbers)

!     determine sum of eigenvalues
        call get_ev_sum_p0 &
             ( occ_numbers, KS_eigenvalue, av_core_shift, ev_sum, &
             ev_sum_shifted )


       !continue to prepare for the actual scf cycle

        if(use_vdw_correction_hirshfeld) then
           write(info_str,'(2X,A)') &
              "Evaluating Hirshfeld Density Update"
           call localorb_info(info_str, use_unit, '(A)', OL_norm)
           call hirsh()
        endif

        converged = .false.
        hartree_delta_energy = 0.d0
        hartree_multipole_correction = 0.d0
        en_elec_delta = 0.d0

!       initialize self-consistency criteria
        deallocate (rho_change)
        allocate(rho_change(n_spin),stat=info)
        call check_allocation(info, 'rho_change                    ')

        call get_free_superpos_energy_p1 &   ! new
             ( partition_tab, &
             rho, rho_gradient, kinetic_density, pot_ion_embed, &
             hartree_energy_free, en_elec_free, en_xc, en_pot_xc, en_ion_ion, &
             en_ion_embed, en_density_embed, en_vdw, en_ll_vdw, en_ll_vdw_err, &
             en_lda_c, en_pbe_c, ionic_forces, ext_charge_forces &
             )

        ! Compute total energy of first iteration

!       if finite-width smearing was used, must get entropy correction for T=0 total energy
        call get_entropy_correction_p1(KS_eigenvalue, occ_numbers, &
             chemical_potential, chemical_potential_spin, entropy_correction)

        call get_penalty_energy(KS_eigenvector, KS_eigenvector_complex, &
        &                       occ_numbers, condition_penalty, penalty_energy)

        fock_energy=0.d0 ! new
        call get_total_energy &
        ( ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
          en_ion_embed, en_density_embed, &
          en_vdw, en_pot_vdw, en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, en_elec_free, &
          en_elec_delta, fock_energy, &
          hartree_energy_free, hartree_delta_energy, &
          hartree_multipole_correction, &
          entropy_correction, penalty_energy, total_energy &
        )

        first_integration = .false.

      endif      ! atom_bsse_ends here
 !--------------------------------------------

      if(orthonormalize_evs .and. (.not. calculate_atom_bsse)) then
         call get_timestamps(time_on_evs,clock_time_on_evs)
         if(.not. use_scalapack) then
            call orthonormalize_eigenvectors()
         elseif(.not. use_elsi_dm) then
            if(real_eigenvectors) then
               ! KS_eigenvectors is output only, if collect_eigenvectors is true
               call orthonormalize_eigenvectors_scalapack_real(KS_eigenvector)
            else
               call orthonormalize_eigenvectors_scalapack_complex(KS_eigenvector_complex)
            endif
         else
            call extrapolate_dm_scalapack()
         endif

         call get_timestamps(rtime,clock_rtime)
         time_on_evs = rtime - time_on_evs
         tot_time_on_evs = tot_time_on_evs + time_on_evs
         call sync_timing(time_on_evs)
         clock_time_on_evs = clock_rtime - clock_time_on_evs
         tot_clock_time_on_evs = tot_clock_time_on_evs + clock_time_on_evs
      endif


      if(flag_delta_rho_in_multipole .and. (.not. calculate_atom_bsse ))then
         rho_radial_multipole_old = 0.d0
         rho_multipole_old = 0.d0
         multipole_feed_back_parameter = 1.0d0
      end if

      if ((use_forces).and.(sc_iter_limit.le.1)) then
        ! too few scf iterations - must switch force computation on right away.
        write(info_str,'(2X,A)') &
        "Too few scf iterations specified."
        call localorb_info(info_str,use_unit,'(A)',OL_norm)
        write(info_str,'(2X,A)') &
        "Switching on the force computation right away."
        call localorb_info(info_str,use_unit,'(A)',OL_norm)
        call localorb_info('',use_unit)

        forces_on = .true.
        pulay_forces_on = .true.
        if (use_gga) then
          gga_forces_on = .true.
          if (use_meta_gga) then
            meta_gga_forces_on = .true.
          endif
        end if
        if (use_analytical_stress) then
          write(info_str,'(2X,A)') &
          "Switching on the stress tensor computation right away."
          call localorb_info(info_str,use_unit,'(A)',OL_norm)
          call localorb_info('',use_unit)
          AS_stress_on    = .true.
          AS_pulay_stress_on = .true.
          ! FIXME: Most probably problematic:
          ! AS_stress_on --> true for sc_iter_limit.le.2
          ! AS_pulay_stress_on --> true for sc_iter_limit.le.1
        end if

        if (use_embedding_pp.and.use_nonlinear_core) then
           nlcc_forces_on = .true.
        end if
      else
        forces_on = .false.
        gga_forces_on = .false.
        meta_gga_forces_on = .false.
        pulay_forces_on = .false.
        AS_stress_on       = .false.
        AS_pulay_stress_on = .false.
        nlcc_forces_on = .false.
      end if

!       For Hartree-Fock calculation, the initialization is done here.
!       We hereby construct the product basis, calculate its overlap integrals
!       with products of individual basis functions, and the Coulomb interaction
!       matrix for product basis functions.
!       FOR BSSE: skip this initialization and shuffle the matrix elements
      if (use_hartree_fock .and. (.not. calculate_atom_bsse )) then
         if(use_periodic_hf .and. .not. use_hf_kspace)then
            if(output_level .eq. 'MD_light') output_priority = OL_high
            call initialize_hartree_fock_p0()
            if(output_level .eq. 'MD_light') output_priority = OL_low
         else
            call initialize_hartree_fock()
         endif
      endif

      if (calculate_atom_bsse ) then
         if (use_hartree_fock ) then
            call allocate_hartree_fock_bsse()
            call bsse_nghr_mapping()
         endif
         if (use_rpa_ene ) then
            call bsse_nghr_mapping() ! this is enough for pbe+rpa calculations
         endif
      endif


      ! Finally, a setting for the Pulay mixing: Since this is
      ! a reinitialization from a given wave function, the charge density from the
      ! initialization (superposition of atomic densities) will not be included in
      ! the mixing procedure.
      if (calculate_atom_bsse) then
         first_iter_mixing= .true.
      else
         first_iter_mixing = .false.
      end if

      ! for reinitialization of convergence criteria
      ev_sum = 0.d0
      previous_forces = 0.d0
      total_energy = 0.d0

      converged = .false.

      mpb_converged = .True.
      if (solvent_method.eq.SOLVENT_MPB) then
         call deallocate_lpbsolver
         call deallocate_mpbsolver
         call allocate_lpbsolver
         if (.not.solve_lpbe_only) then
           call allocate_mpbsolver
         else
           call allocate_mpbsolver_dummy
           if (.not.allocated(ln_cosh_function)) then
             allocate(ln_cosh_function(n_full_points),stat=info2)
             call check_allocation(info2, 'ln_cosh_function                  ')
             ln_cosh_function(:) = 0d0
           end if
         end if
      !  reevaluate dielectric function and ionic exclusion function
         if (allocated(free_rho_superpos_copy)) deallocate(free_rho_superpos_copy)
         allocate(free_rho_superpos_copy(n_spin,n_full_points))
         if (allocated(free_rho_gradient_superpos_copy)) deallocate(free_rho_gradient_superpos_copy)
         allocate(free_rho_gradient_superpos_copy(3,n_spin,n_full_points))
         free_rho_gradient_superpos_copy = 0d0
         free_rho_superpos_copy = 0d0
         if (dynamic_cavity) then
           call dielecfunc_from_density(rho,dielec_func_mpb)
           call dielecfunc_gradient_from_density_gradient(dielec_func_mpb,dielec_func_gradient_mpb,&
             rho, rho_gradient) !*pi4_inv
         else
           free_rho_superpos_copy(1,:) = free_rho_superpos(:)
           do i_coords=1, 3, 1
             free_rho_gradient_superpos_copy(i_coords,1,:) = free_rho_gradient_superpos(i_coords,:)
           end do
           call dielecfunc_from_density(free_rho_superpos_copy*pi4_inv,dielec_func_mpb)
           call dielecfunc_gradient_from_density_gradient(dielec_func_mpb,dielec_func_gradient_mpb,&
             free_rho_superpos_copy*pi4_inv, free_rho_gradient_superpos_copy) !*pi4_inv
         end if
         if (dynamic_ions) then
           call alphafunc_from_density( rho, rho_gradient, dielec_func_mpb, &
           n_full_points)
         else
           free_rho_superpos_copy(1,:) = free_rho_superpos(:)
           do i_coords=1, 3, 1
             free_rho_gradient_superpos_copy(i_coords,1,:) = free_rho_gradient_superpos(i_coords,:)
           end do
           call alphafunc_from_density( free_rho_superpos_copy*pi4_inv, free_rho_gradient_superpos_copy, &
             dielec_func_mpb, n_full_points)
         end if
         call prepare_SPE_solver
         !evaluate gradient of superposition of free atom potentials
         call evaluate_v_free_gradient
         !reset convergence flags
         mpb_converged = .False.
         mpb_solver_started = .False.
         initialized_lpb_solver = .False.
         mpb_forces_on = .False.
         Gnonmf_forces_on = .False.
      else
        call allocate_lpbsolver_dummy
      end if !solvent_method.eq.SOLVENT_MPB
      ! AJL: Turn first_integration off again if we forced it on to make meta-GGAs work
      if (use_meta_gga.and.first_integration) first_integration = .false.

      call reshape_matrices_n_states(n_states_init)

      ! SCF initialization ends here
      call get_timestamps( rtime, clock_rtime )
      time_sc_loop = rtime - time_sc_loop
      clock_time_sc_loop = clock_rtime - clock_time_sc_loop
      call sync_timing(time_sc_loop)

      number_of_scf_inits = number_of_scf_inits + 1

      write(info_str,'(2X,A)') ""
      call localorb_info(info_str,use_unit,'(A)', OL_norm)
      call output_timeheader(deffmt, "End scf reinitialization - timings")
      call output_times(deffmt, "Time for scf. reinitialization", &
      &                 time_sc_loop, clock_time_sc_loop)
      call output_times(deffmt, "Boundary condition initialization", &
      &                 time_bc_init, clock_time_bc_init)
      call output_times(deffmt, "Integration", &
      &                 time_intg, clock_time_intg)
      call output_times(deffmt, "Grid partitioning", &
      &                 time_partition, clock_time_partition)
      call output_times(deffmt, "Preloading free-atom quantities on grid", &
      &                 time_free_offset, clock_time_free_offset)
      call output_times(deffmt, "Free-atom superposition energy", &
      &                 time_superpos, clock_time_superpos)
      if (orthonormalize_evs) then
        call output_times(deffmt, "K.-S. eigenvector reorthonormalization", &
        &                 time_on_evs, clock_time_on_evs)
      end if
      if (use_hartree_fock) then
         call ovlp3fn_timings(deffmt)
      endif
      write(info_str,'(A)') &
      "------------------------------------------------------------"
      call localorb_info(info_str,use_unit,'(A)', OL_norm)

      if(output_level .eq. 'MD_light' .and. MD_stepcount .gt. 1) output_priority = 1

      ! Allocatable arrays that are tracked
      if (allocated(ham_ovlp_work)) call aims_deallocate( ham_ovlp_work, "ham_ovlp_work" )

    end subroutine reinitialize_scf
!******
