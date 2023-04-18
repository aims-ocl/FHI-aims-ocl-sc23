!****s* FHI-aims/initialize_scf
!  NAME
!    initialize_scf
!  SYNOPSIS

    subroutine initialize_scf ( converged )

!  PURPOSE
!  High-level wrapper around the former "first" scf loop
!
!  Here, we "initialize the integrals" (i.e. we compute the
!  overlap matrix and set up all grids), we solve the KS eigenvalue problem
!  for the first time to produce an initial wave function, we allocate all
!  needed quantities across the grid, we store all "fixed" free-atom like
!  quantities on the grid, and we set the initial running parameters for the scf loop.
!
!  USES

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
      use hartree_non_periodic_ewald
      use hartree_fock
      use plus_u
      use scalapack_wrapper
      use separate_core_states
      use KH_core_states
      use rel_x2c_mod
      use precondition
      use SPE_solver
      use restart
      use force_occupation
      use thermodynamic_integration
      use wf_extrapolation
      use relaxation
      use load_balancing
      use hartree_potential_storage, only: initialize_hartree_potential_storage
      use geometry, only : total_initial_charge
      use pseudodata
      use vdw_correction
      use mbd_std_wrapper, only: mbd_std_initialize, mbd_self_consistent
      use hirshfeld, only: hirshfeld_initialize
      use lpb_solver_utilities, only: dielecfunc_from_density, &
          dielecfunc_gradient_from_density_gradient,alphafunc_from_density,&
          dynamic_cavity,allocate_lpbsolver_dummy,&
          allocate_lpbsolver,mpb_start_crit,start_at_scf_step,dielec_func_mpb,&
          dielec_func_gradient_mpb,dielec_func_mpb,&
          solve_lpbe_only,get_cut_f,decrement_kind,&
          correct_dielectric_decrement,initialize_eps_s_func,dynamic_ions,deps_drho,&
          d2eps_drho2,dielec_func_grad2,d2eps2_drho2,deps2_drho,dielec_func2,evaluate_v_free_gradient,&
          mpb_converged, initialized_lpb_solver,mpb_solver_started, limit_l_max_SPE
      use mpb_solver_utilities, only: allocate_mpbsolver,ln_cosh_function,reg_method, allocate_mpbsolver_dummy
      use elsi_wrapper, only: eh_scf,aims_elsi_init_scf,&
          aims_elsi_compute_mu_and_occ
      use initialize_grid_storage, only: initialize_grid_storage_p1
      use aims_memory_tracking, only: aims_allocate, aims_deallocate, &
          aims_mem_current_output
      use hartree_potential_recip, only: initialize_recip_hartree_potential
      use restart_elsi, only: elsi_restart_check
      use free_atoms, only: free_rho_spl,free_rho_spl,free_drho_dr_spl,free_pot_es_spl,free_pot_es_at_zero
      use rlsy_interface, only: rlsy_h,generate_rlsy_handle
      use generate_index_map, only: rl_generate_index_map,rl_collect_overlap,rl_collect_hamiltonian
      use batch_rlsy_interface, only: batch_irr, compute_batch_rlsy_interface

      implicit none

!  ARGUMENTS

      logical :: converged

!  INPUTS
!    o none
!  OUTPUTS
!    o converged -- if force_potential then .true. otherwise initialized to .false.
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
!  SOURCE
!

      ! imported variables


      logical :: info

      integer :: info2

      ! local variables

      integer :: task_distribution_method, max_full_points
      real*8 :: temp_time, clock_temp_time
      real*8 :: penalty_energy
      character(*), parameter :: deffmt = '2X'

      real*8, allocatable :: ham_ovlp_work(:)


      character*8  :: cdate
      character*10 :: ctime
      real*8 :: time_relax_init(4)

!      integer :: nnz
      integer :: i_state
      integer :: i_electron
      integer :: i_n, i_l, i_m, i_k, l_max

      character*100 :: info_str

      real*8 :: clock_ke_time, ke_time

      ! counters

      integer :: i_spin, allocation_info
      integer :: i_error

      real*8, allocatable, dimension(:,:) :: free_rho_superpos_copy
      real*8, allocatable, dimension(:,:,:) :: free_rho_gradient_superpos_copy
      integer :: i_coords
      complex*16 :: dummy_complex(1)

      integer :: rlsy_verbosity

      character(*), parameter :: func = 'initialize_scf'

      real*8 time0, time_work, time2, time_work2

      ! begin work

      if(output_level .ne. 'full') output_priority = OL_low  ! Even print those

!     First scf iteration is started separately - if we need to adapt all grids,
!     or starting density / potential, etc etc, then this must be done here

      write(info_str,'(A)') ''
      call localorb_info ( info_str )

      write(info_str,'(A)') "------------------------------------------------------------"
      call localorb_info( info_str )

      write(info_str,'(10X,A)') "Begin self-consistency loop: Initialization."
      call localorb_info( info_str )

      write(info_str,'(A)') ''
      call localorb_info ( info_str )

      call date_and_time(cdate, ctime)
      write(info_str,'(10X,A,A,A,A)') "Date     :  ", cdate, ", Time     :  ", ctime
      call localorb_info ( info_str )

      write(info_str,'(A)') "------------------------------------------------------------"
      call localorb_info( info_str )

      call get_timestamps( time_sc_loop, clock_time_sc_loop)

      ! minor repetition - this counter is used for relaxation or MD, to keep track of the
      ! time used for one fully converged (single-point) scf cycle. Must be initialized here
      ! as well as in reinitialize_scf()
      call get_timestamps ( time_sc_tot, clock_time_sc_tot )

      ! initialize lists of integration centers, Hartree potential
      ! multipole centers, etc, whose shape depends on the boundary
      ! conditions (periodic or not)
      time_bc_init = 0.d0
      call get_timestamps(time_bc_init, clock_time_bc_init)
      call initialize_bc_dependent_lists()
      call get_times(time_bc_init, clock_time_bc_init, tot_time_bc_init, tot_clock_time_bc_init)


      ! if we have relaxation, must allocate / initialize a few things.  Do
      ! this after initialize_bc_dependent_lists() as we need the outer radii.
      time_relax_init = 0.d0
      if (use_geo_relaxation) then
        call start_timer(time_relax_init)
        if(output_level .ne. 'full') output_priority = OL_norm
        call initialize_relaxation(coords, lattice_vector)
        if(output_level .ne. 'full') output_priority = OL_low
        call stop_timer(time_relax_init)
      end if

      ! If Ewald's decomposition of the Hartree potential is to be
      ! used in the non-periodic case, then perform the initialization.
      if (use_hartree_non_periodic_ewald) then
         call initialize_hartree_non_periodic_ewald
      end if

      if (n_periodic > 0) then
         call initialize_recip_hartree_potential
      end if

      if (use_scalapack) then
         ! If ScaLapack is really used, initialize all required arrays
         call initialize_scalapack( )
      end if

      ! Initialize load balancing data
      if(use_load_balancing) call reset_load_balancing()

      call initialize_separate_core_states

      call initialize_KH_core_states

      ! Initialize ELSI
      ! Must be after initialize_bc_dependent_lists and initialize_scalapack
      call aims_elsi_init_scf(my_scalapack_comm_work, my_blacs_ctxt, mb, mxld, &
           mxcol,my_k_point)

      ! Must be after elsi_init
      call elsi_restart_check()

      if (use_wf_extrapolation) then
         ! Can only be done after initialize_bc_dependent_lists()
         ! because it needs to know the k_point architecture.
         call initialize_wf()
      end if

      if (use_plus_u) then
         call allocate_plus_u
      endif

! DB: preparing nonlocal part of the pseudopotentials
      if (use_embedding_pp) then
          call get_timestamps(time_embedding, clock_time_embedding)
!          call baswave_pp_overlap2()
          call baswave_pp_overlap3()
!          if(use_forces) then
!             call baswave_pp_overlap()
!          endif
          call get_nonlocal_pot()
          if(use_forces) then
!             call get_pp_nonlocal_force_matrix()
          endif
          call get_times(time_embedding, clock_time_embedding, tot_time_embedding, tot_clock_time_embedding)

      end if

!     integrate T+V

!     Calculate both overlap and Hamiltonian matrix
!     * If desired, optimize all integration grids
!     * Tabulate all needed quantities on initial integration grids
!     * Compute overlap _and_ Hamiltonian matrix as test quantities for adaptive grids.


      if (ALL(specified_grid).and.ALL(angular_acc == 0.0d0)) then
         ! fixed grid: get grid storage first and matrices only afterwards

         call get_timestamps(time_partition, clock_time_partition)

         call initialize_fixed_grids()

         call partition_grid()

         call get_times(time_partition, clock_time_partition, tot_time_partition, tot_clock_time_partition)

         call get_timestamps(temp_time, clock_temp_time)

         call allocate_grid_storage()

         if (use_embedding_pp) then
             call get_timestamps(time_embedding, clock_time_embedding)
             call put_localpot_on_integration_grid()
     !        call sync_vector(whole_local_pseudpot_on_intgrid, n_full_points)
             if(use_nonlinear_core) then
                 call partial_core_dens2grid_v2()
             endif

             call get_times(time_embedding, clock_time_embedding, tot_time_embedding, tot_clock_time_embedding)
         end if


         call initialize_grid_storage_p1 &
              ( partition_tab, weight_tab, hartree_partition_tab, hartree_weight_tab, hartree_partition_deriv, &
              partition_deriv_delley, hartree_potential, free_hartree_superpos, free_rho_superpos, &
              free_rho_gradient_superpos, rho, rho_gradient, pot_ion_embed )

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
           ! If we are interested in having a scf density and Hartree potential that has the
           ! exact charge as given in control.in and as integrated on the actual integration grid,
           ! we must enforce this zero here (in addition to the preceding if, which also included
           ! normalize_initial_density).
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
           call renormalize_free_density (free_rho_superpos, free_rho_gradient_superpos, 1, partition_tab, n_electrons)
           !
           ! What is not renormalized now is the electrostatic potential of the free atom
           ! superposition - because that is currently exactly consistent point by point with
           ! the nuclear potential.
           ! It is a good question whether this will become a problem later.
        end if

        ! When using PEXSI, a copy of hartree_potential is needed to calculate
        ! the change of potential during SCF
        if(use_elsi_dm .and. elsi_solver == 3) then ! PEXSI
           if(allocated(hartree_potential_save)) then
              call aims_deallocate(hartree_potential_save,"hartree_potential_save")
           endif

           call aims_allocate(hartree_potential_save,n_full_points,"hartree_potential_save")

           hartree_potential_save = hartree_potential
        endif

        if ( use_rlsy_symmetry ) then
            rlsy_verbosity = 1
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
            call compute_batch_rlsy_interface()
        endif


         call get_n_compute_maxes_p1( partition_tab )


!        when autoselect density update method is set, the final decision about
!        packed_matrix_format is made in get_n_compute_maxes_p1
         if( packed_matrix_format /= PM_none ) then
            if(n_periodic .eq. 0) then
               call initialize_packed_matrix_cluster
            end if
         end if

         call allocate_matrices ( )

         if (force_occupation_basis) then
           call allocate_force_occ_state
         end if

         call get_times(temp_time, clock_temp_time)

         ! When the partioning is known, we can set the index for the sparse local matrices
         if(use_local_index) then
            call set_index_hamiltonian

            if (use_scalapack) then
               call init_comm_sparse_local_matrix_scalapack
            else if(.not.use_scalapack .and. .not.use_load_balancing) then
               call aims_stop("Sparse local matrices not supported when &
                              &using LAPACK eigensolver.", func)
            end if
         endif

         !TODO BL: What has to be done, do distribute the Hamiltonian right?
         if(use_local_index .and. .not.use_scalapack) then
            call set_index_hamiltonian
         endif

         call get_timestamps(time_intg, clock_time_intg)

            !OTH/CMH: Loading density
         if(restart_from_density) then
           open(file = "CMH-Density.out", unit = 27, status = 'old', form = 'unformatted', iostat = i_error)
           close(27)
           if (i_error .eq. 0) then
           open(file = "CMH-Density.out", unit = 27, status = 'old', form = 'unformatted', iostat = i_error)
           read(27) rho
           close(27)
           open(file = "CMH-Density-grad.out", unit =27, status = 'old', form = 'unformatted', iostat = i_error)
           read(27) rho_gradient
           close(27)

           !open(27, file="CMH-Density-values.out",form="formatted",status="replace")
           !write(27,*) rho
           !close(27)
           endif
         endif

         ! this has to come after grid initialization (mbdvdw allocates grid
         ! array - vdw xc potential) and before hamiltonian construction
         if (use_mbd_std) then
             call mbd_std_initialize()
             call hirshfeld_initialize()
         end if

         if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then ! allocate V,T,W,S matrices
            call rel_matrix_initialize()
         endif

         ! Integrate the initial Hamiltonian matrix

         ! Technically, this need not be integrated in the case of a restart, because we never use the
         ! results for the next electronic structure updates.
         ! There are a bunch of methods though that rely on the existence of the initial Hamiltonian matrix
         ! for setup tasks (matrix packing is one example), so it is bes to just integrate it anyway.

         if(use_local_index) then
            ! If use_load_balancing is set, this run serves also as calibrating run for the next ones
            if (use_load_balancing) get_batch_weights = .true.

            if (use_scalapack) then
               ! For ScaLAPACK, we create a temporary work matrix, integrate the real-space Hamiltonian using
               ! the *sparse* local matrix format, convert the real-space Hamiltonian into the Kohn-Sham
               ! Hamiltonian to be passed to the eigensolver, then delete the temporary work matrix.
               ! We always use the sparse local matrix format, even when our target is the full local matrix
               ! format, because the former is always available, whereas the later is being set up by this
               ! run
               call aims_allocate( ham_ovlp_work, n_hamiltonian_matrix_size*n_spin, "ham_ovlp_work" )
               call integrate_real_hamiltonian_matrix_p2 (hartree_potential, rho, rho_gradient, kinetic_density, partition_tab, l_shell_max, &
                    en_xc, en_pot_xc, ham_ovlp_work, en_vdw, en_pot_vdw)
               call set_sparse_local_ham_scalapack(ham_ovlp_work)
               call aims_deallocate( ham_ovlp_work, "ham_ovlp_work" )
            else
               ! For LAPACK, we don't have the sparse local format implemented.  Thus, we integrate the
               ! real-space Hamiltonian here solely to set up the full local matrix format.
               call integrate_real_hamiltonian_matrix_p2 (hartree_potential, rho, rho_gradient, kinetic_density, partition_tab, l_shell_max, &
                    en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )
            endif

            if(use_load_balancing) then
               ! Using the timings that we measured when integrating the real-space Hamiltonian, set up
               ! the batch distributions
               get_batch_weights = .false.
               call compute_balanced_batch_distribution(n_bp_integ)
               call init_comm_full_local_matrix(batch_perm(n_bp_integ)%n_basis_local,batch_perm(n_bp_integ)%i_basis_local)

               if (.not. use_scalapack) then
                  ! In the LAPACK case, now that we have the batch distribution set up, we integrate the
                  ! real-space Hamiltonian "for real" this time.
                  use_batch_permutation = n_bp_integ
                  call aims_deallocate(hamiltonian, "hamiltonian")
                  call aims_deallocate(overlap_matrix, "overlap_matrix")
                  call aims_allocate(hamiltonian, batch_perm(use_batch_permutation)%n_local_matrix_size, n_spin, "hamiltonian")
                  call aims_allocate(overlap_matrix, batch_perm(use_batch_permutation)%n_local_matrix_size, "overlap_matrix")
                  call integrate_real_hamiltonian_matrix_p2 (hartree_potential, rho, rho_gradient, kinetic_density, partition_tab, l_shell_max, &
                       en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )
                  use_batch_permutation = 0
               end if
            end if
         else
            call integrate_real_hamiltonian_matrix_p2 (hartree_potential, rho, rho_gradient, kinetic_density, partition_tab, l_shell_max, &
                 en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )
         endif

         ! Integrate the overlap matrix after the hamiltonian.
         ! If load balancing is in effect, we can already profit from it

         if(use_local_index) then
            if (use_load_balancing) use_batch_permutation = n_bp_integ

            if (use_scalapack) then
              ! Like before, for ScaLAPCAK we create a work array, integrate the real-space overlap matrix,
              ! create the Kohn-Sham overlap matrix, then delete the temporary array.  The difference from
              ! the last part is that we've set up the load balancing, so we can use it here.
              if(use_load_balancing) then
                call aims_allocate( ham_ovlp_work, batch_perm(n_bp_integ)%n_local_matrix_size, "ham_ovlp_work" )
              else
                call aims_allocate( ham_ovlp_work, n_hamiltonian_matrix_size, "ham_ovlp_work" )
              endif

              call integrate_ovlp_matrix_p2 ( partition_tab, l_shell_max, ham_ovlp_work )

              if (use_load_balancing) then
                call set_full_local_ovlp(ham_ovlp_work, overlap_matrix, dummy_complex)
              else
                call set_sparse_local_ovlp_scalapack(ham_ovlp_work)
              endif
              call aims_deallocate( ham_ovlp_work, "ham_ovlp_work" )
            else
              ! As before, we'll create the KS overlap matrices later, so we need only store
              ! the overlap matrix here
              if (use_load_balancing) then
                call integrate_ovlp_matrix_p2 ( partition_tab, l_shell_max, overlap_matrix )
              else
                call aims_stop("Sparse local matrices not supported when &
                               &using LAPACK eigensolver.", func)
              end if
            end if

            use_batch_permutation = 0
         else
            call integrate_ovlp_matrix_p2 ( partition_tab, l_shell_max, overlap_matrix )
         end if

         call get_timestamps(rtime, clock_rtime)
         time_intg =  rtime - time_intg
         tot_time_intg = tot_time_intg + time_intg
         call sync_timing(time_intg)

         clock_time_intg =  clock_rtime - clock_time_intg
         tot_clock_time_intg = tot_clock_time_intg + clock_time_intg

      else

!        for fixed grids the following statements have to go after get_n_compute_maxes_p1
!        here we can call them right now
         if( packed_matrix_format /= PM_none ) then
            if(n_periodic .eq. 0) then
               call initialize_packed_matrix_cluster
            end if
         end if

         call allocate_matrices ( )

         if (force_occupation_basis) then
           call allocate_force_occ_state
         end if
         ! adapted grids: get matrices first and grid storage only afterwards
         call estimate_n_compute_maxes( )

         call get_timestamps(time_intg, clock_time_intg)


         call initialize_integrals_p0 ( l_shell_max, overlap_matrix, hamiltonian )
         call get_timestamps(rtime,clock_rtime)
         time_intg =  rtime - time_intg
         tot_time_intg = tot_time_intg + time_intg
         call sync_timing(time_intg)
         clock_time_intg =  clock_rtime - clock_time_intg
         tot_clock_time_intg = tot_clock_time_intg + clock_time_intg

         call get_timestamps(time_partition, clock_time_partition)

         call partition_grid()

         call get_times(time_partition, clock_time_partition, tot_time_partition, tot_clock_time_partition)

         call get_timestamps(temp_time, clock_temp_time)
         call allocate_grid_storage()

         if (use_embedding_pp) then
             call get_timestamps(time_embedding, clock_time_embedding)
             call put_localpot_on_integration_grid()
             if(use_nonlinear_core) then
                 call partial_core_dens2grid_v2()
             endif

             call get_times(time_embedding, clock_time_embedding, tot_time_embedding, tot_clock_time_embedding)
         end if

         call initialize_grid_storage_p1 ( partition_tab, weight_tab, hartree_partition_tab, hartree_weight_tab, hartree_partition_deriv, &
              partition_deriv_delley, hartree_potential, free_hartree_superpos, free_rho_superpos, free_rho_gradient_superpos, &
              rho, rho_gradient, pot_ion_embed )

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
           call renormalize_density (rho, rho_gradient, kinetic_density, n_spin, partition_tab, n_electrons)
           ! extra call for free_rho_superpos, which has a factor 4 pi multiplied in.
           call renormalize_free_density (free_rho_superpos, free_rho_gradient_superpos, 1, partition_tab, n_electrons)
           !
           ! What is not renormalized now is the electrostatic potential of the free atom
           ! superposition - because that is currently exactly consistent point by point with
           ! the nuclear potential.
           ! It is a good question whether this will become a problem later.
         end if

         call get_n_compute_maxes_p1( partition_tab )

         call get_times(temp_time, clock_temp_time)

      end if

      if (use_embedding_pp) then
          call get_timestamps(time_embedding, clock_time_embedding)
          call add_nonlocal_pot(hamiltonian)
          call get_times(time_embedding, clock_time_embedding, tot_time_embedding, tot_clock_time_embedding)
      end if


      !   ! make a note that the matrices are integrated already at least once
      !first_integration = .false.


      if(.not.use_local_index .and. packed_matrix_format == PM_index .and. (.not. use_periodic_hf).and.(.not.use_DFPT_phonon) &
         .and. (.not. use_DFPT_dielectric).and.(.not.use_DFPT_polarizability).and.(.not.use_DFPT_polar_reduce_memory) )then
         call remove_small_numbers_in_hamiltonian_and_ovlp(hamiltonian, overlap_matrix,info)
         if(info) call reshape_matrices
      end if

      if (.not.use_local_index .and. out_matrices) then
         if (out_matrices_format_2005) then
            call output_real_matrices ( overlap_matrix, hamiltonian )
         else
            call output_real_matrices_aij ( overlap_matrix, hamiltonian )
         end if
      end if

!     >>> AB: feb 2012
      if (.not.use_local_index .and. out_aitranss) then
        call output_ka_overlap(overlap_matrix)
      end if
!     <<< done with insert: AB, feb 2012


      ! For fully-relativistic cases, The H and S' matrices that will be used for KS solvers will be generated here.
      ! With spinor integration terms (with k phase factor applied to it) generated for V, T, W, S matrices in former codes,
      ! we now assemble the X2C Hamiltonian and overlap matrices for the first step diagonalization.
      if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
      if(flag_rel.eq.REL_x2c) call x2c_1st_iteration()
      if(flag_rel.eq.REL_4c_dks) call q4c_1st_iteration()
      ! (Rundong) The present code (from the following line to 1007) has been too complex. Thus, I wrap all the relativistic
      ! parts in subroutine x2c_1st_iteration. After this, we skip to 1007 to go on. They might be merged some day.
         goto 1007
      endif

      if ( use_rlsy_symmetry ) then
          !call rl_collect_overlap(rlsy_h,overlap_matrix)
          !call rl_collect_hamiltonian(rlsy_h,hamiltonian)
      endif


      ! If use_scalapack is set, allocate and set the distributed overlap matrix now.
      ! Under certain circumstances we can deallocate the (non-distributed) overlap_matrix
      ! here and save some non-scalable memory
      if (.not.use_local_index .and. use_scalapack) then
         if(n_periodic>0 .or. packed_matrix_format /= PM_none ) then
            call construct_overlap_scalapack(overlap_matrix)
         else
            call setup_overlap_scalapack(overlap_matrix)
         endif
         if(use_constraint .or. flag_run_mulliken .or. out_lowdin .or. out_band .or. out_band_mulliken .or. out_matrices &
            .or. flag_KS==-2 .or. out_dos .or. use_plus_u .or. use_out_eigenvec_ovlp .or.transport_calculation .or.transport_lead_calculation &
            .or. flag_run_tddft_real_time_propagation .or. fo_finalise & 
            .or. use_DFPT_phonon.or.use_DFPT_phonon_reduce_memory &
            .or. use_DFPT_polarizability.or.use_DFPT_polar_reduce_memory  & 
            .or. use_DFPT_dielectric ) then
            write(info_str,'(2X,A)') "Not deallocating overlap matrix."
            call localorb_info(info_str,use_unit,'(A)')
         else
            write(info_str,'(2X,A)') "Deallocating overlap matrix."
            call localorb_info(info_str,use_unit,'(A)')
            call aims_deallocate( overlap_matrix,  "overlap_matrix")
            call aims_allocate( overlap_matrix, 1, "overlap_matrix") ! make it defined
         endif

         ! If scalapack is used, the hamiltonian has to be stored
         ! into the scalapack arrays

         if(n_periodic>0 .or. packed_matrix_format /= PM_none ) then
            call construct_hamiltonian_scalapack(hamiltonian)
         else
            call setup_hamiltonian_scalapack(hamiltonian)
         endif
      endif

      ! Save the scalapack overlap in the case it is needed for post processing
      if(use_scalapack .and. (flag_run_mulliken .or. out_lowdin .or. flag_out_ev_scalapack_hdf5 .or. out_band_mulliken)) &
           call save_overlap_scalapack

      if (use_plus_u) then
             ! only used for plus_u with mulliken charges!
             call plus_u_init_idx
      endif

!     After initial integrations, solve eigenvalue problem for the first time ...
!     ... UNLESS this is a restart run, where the initial eigenvectors are obtained from
!     a separate file anyway!

      call get_timestamps(time_diag, clock_time_diag)

      if(keep_restart_info .and. (use_scalapack.and.(.not.use_density_matrix)) .or. use_mpi)then
         call sync_logical(restart_file_exists, SYNC_OR)
         call sync_logical(restart_read, SYNC_OR)
      endif

      if ( .not. (keep_restart_info.and.restart_file_exists) ) then
        if(elsi_read_dm) then
            ! Restart with ELSI
            restart_zero_iteration = .true.
        else
            ! This is the standard case without a restart - we obtain the
            ! Kohn-Sham eigenvalues from a Kohn_Sham solution

            ! Code treating the Kohn-Sham eigenvalue problem has been moved
            ! inside advance_KS_solution.
            call advance_KS_solution(overlap_matrix, hamiltonian, n_electrons, KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
                 occ_numbers, chemical_potential, chemical_potential_spin)

            call reshape_matrices_n_states(n_states)

            restart_zero_iteration = .false.
        endif
      else if ( fo_finalise.and..not.(fo_comb_m.eq.'scalap')) then
            ! this is the case for fo-dft run
            call read_restart_info()
            restart_zero_iteration = .true.
      else if ( fo_finalise.and.(fo_comb_m.eq.'scalap') ) then
            restart_zero_iteration = .true.
            ! same as below, "solve eigenvalue problem anyway"

            ! get_KS_orbitals_p0 replaced by advance_KS_solution with exactly
            ! the same arguments. This will make no difference unless one
            ! changes the control.in file.
            call advance_KS_solution(overlap_matrix, hamiltonian, n_electrons, KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
                 occ_numbers, chemical_potential, chemical_potential_spin)

            call reshape_matrices_n_states(n_states)
      else if ( .not.(use_scalapack .and. use_density_matrix) ) then
            ! this is the normal case of a restart - just read the old eigenvectors right here.
            if (force_occupation_projector) then
              call allocate_previous_eigenvector
            end if

            call read_restart_info()

            call aims_elsi_compute_mu_and_occ(eh_scf, n_electrons, n_states, n_spin, n_k_points, k_weights, KS_eigenvalue, occ_numbers, &
                 chemical_potential)

            chemical_potential_spin(1:n_spin) = chemical_potential

            restart_zero_iteration = .true.
      else if ( (use_scalapack .and. use_density_matrix).and.(force_single_restartfile.or.restart_eigenvectors_periodic) ) then
            ! in this case, even for scalapack runs ("restart from density matrix"), the default
            ! LAPACK routine is used to avoid writing two different restart files (..and the need
            ! to rotate the density matrix representation of the KS_eigenvector.
            ! This is necessary for all calculations with use_density_matrix and
            ! rotated restart files.
            if (force_occupation_projector) then
              call allocate_previous_eigenvector
            end if

            call read_restart_info()

            call aims_elsi_compute_mu_and_occ(eh_scf, n_electrons, n_states, n_spin, n_k_points, k_weights, KS_eigenvalue, occ_numbers, &
                 chemical_potential)

            chemical_potential_spin(1:n_spin) = chemical_potential

            restart_zero_iteration = .true.

            call spread_eigenvectors_scalapack(KS_eigenvector, KS_eigenvector_complex)
            ! at this point, the correct KS_eigenvalue and KS_eigenvector are known,
            ! but for use_density_matrix the densmat needs to be re-calculated from
            ! the KS orbitals.
            ! As in the default case without rotations this happens in
            ! update_density_densmat.f90
      else if ( (use_scalapack .and. use_density_matrix).and..not.(force_single_restartfile.or.restart_eigenvectors_periodic) ) then
            ! This is the case of a restart with scalapack used and the density is obtained from the density matrix.
            ! Here, it is much more efficient to restart directly from the last density matrix,
            ! rather than reading in the Kohn-Sham orbitals from a previous run
            ! In this case, reading the initial information happens directly
            ! in update_density_densmat.f90 , where the density matrix
            ! is actually known.
            !
            ! We just set a dummy flag here.
            restart_zero_iteration = .true.

            ! ... and, we solve the eigenvalue problem anyway. Although these eigenvalues
            ! will not enter the density matrix later, this prevents unnecessary bogus
            ! output for the initialization. (Why not modify the output part??)
            ! This is also an obvious case for savings in case the restart takes too much time.

            ! get_KS_orbitals_p0 replaced by advance_KS_solution with exactly
            ! the same arguments. This will make no difference unless one
            ! changes the control.in file.
            call advance_KS_solution(overlap_matrix, hamiltonian, n_electrons, KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
                 occ_numbers, chemical_potential, chemical_potential_spin)

            call reshape_matrices_n_states(n_states)
      end if

      call evaluate_KH_core_states( overlap_matrix, rho, rho_gradient, kinetic_density, partition_tab, hartree_potential )

      call add_KH_states(KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex)

 1007 continue
                      ! write(use_unit,*)'non-rel 1st step eigenvalues:'
                      ! write(use_unit,*)KS_eigenvalue
                      ! write(use_unit,*)'occ_numbers'
                      ! do i_k=1, n_k_points
                      !   do i_n=1, n_spin
                      !     write(use_unit,"('k=',i4,3x,'spin=',i3,3x,20f10.5)")i_k,i_n,occ_numbers(:,i_n,i_k)
                      !   enddo
                      ! enddo
                      !stop

    !  if (n_atoms.eq.1 .and. n_periodic.eq.0 .and. n_electrons.le.18)
      if (use_aufbau_principle)    call aufbau_principle !exception - single atom. Set occupations according to aufbau principle

      call get_times(time_diag, clock_time_diag, tot_time_diag, tot_clock_time_diag)

      if ((.not. use_elsi_dm) .and. (.not. elsi_read_dm)) then
         call output_real_eigenfunctions ( KS_eigenvalue, KS_eigenvector, occ_numbers)
         ! if we're adjusting any s.c.f. stability related settings,
         ! this must be done right after calling output_real_eigenfunctions,
         ! in which the type of system (low gap or not) is determined
         if (adjust_scf.and.(adjust_scf_iteration.eq.0)) then
            call adjust_scf_settings
            if (.not.adjust_scf_always) then
               ! adjust only once
               adjust_scf = .false.
            end if
         end if
      end if

!     determine sum of eigenvalues
      if ((.not. use_elsi_dm) .and. (.not. elsi_read_dm)) then
         call get_ev_sum_p0 ( occ_numbers, KS_eigenvalue, av_core_shift, ev_sum, ev_sum_shifted )
      end if

!     Eigenvalue problem is solved; next, decide whether to stop with the initial eigenvalue
!     sum only, or whether to continue into the real scf cycle

      if ( use_vdw_correction_hirshfeld_sc .or. (use_mbd_std .and. mbd_self_consistent) ) then
         write(info_str,'(2X,A)') "Evaluating Hirshfeld Density Update"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm )
         call hirsh()
      endif

      if (force_potential.ne.0) then
!       non-selfconsistent loop - need only first iteration
        write(info_str,'(2X,A)') "Non-selfconsistent run - finishing after first iteration."
        call localorb_info ( info_str )

        time_free_offset = temp_time
        tot_time_free_offset = tot_time_free_offset + time_free_offset
        clock_time_free_offset = clock_temp_time
        tot_clock_time_free_offset = tot_clock_time_free_offset + clock_time_free_offset

        call output_ev_sum ( ev_sum, ev_sum_shifted )

        write(info_str,'(A)') ''
        call localorb_info ( info_str )

        converged = .true.

      else
        ! force_potential=0 - the normal case, continue to prepare for the actual scf cycle

        converged = .false.

        if ((use_forces).and.(sc_iter_limit.le.1)) then
          ! too few scf iterations - must switch force computation on right away.
          write(info_str,'(2X,A)') "Too few scf iterations specified."
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)') "Switching on the force computation right away."
          call localorb_info(info_str,use_unit,'(A)')
          call localorb_info('',use_unit)

          forces_on = .true.
          pulay_forces_on = .true.
          if (use_gga) then
            gga_forces_on = .true.
            if (use_meta_gga) then
              meta_gga_forces_on = .true.
            end if
          end if

          if (use_analytical_stress) then
            write(info_str,'(2X,A)') "Switching on the stress tensor computation right away."
            call localorb_info(info_str,use_unit,'(A)')
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
          forces_on          = .false.
          gga_forces_on      = .false.
          meta_gga_forces_on = .false.
          pulay_forces_on    = .false.
          AS_stress_on       = .false.
          AS_pulay_stress_on = .false.
          nlcc_forces_on = .false.
        end if

        hartree_delta_energy = 0.d0
        hartree_multipole_correction = 0.d0
        en_elec_delta = 0.d0

!       initialize self-consistency criteria
        if (.not.allocated(rho_change)) then
          allocate(rho_change(n_spin),stat=allocation_info)
          call check_allocation(allocation_info, 'rho_change                    ')
        end if

        call get_timestamps(time_free_offset, clock_time_free_offset)

!       distribute the MPI tasks according to the atomistic workload obtained from
!       initialize_integrals above

           task_distribution_method = 0
           call distribute_tasks(task_distribution_method)

        ! Next, distribute the spline storage
        ! Currently, this has an effect only on the Kerker multipole spline arrays for the Kerker preconditioner
        call distribute_spline_storage()

        ! if communication_type is shmem_comm, allocate shared memory segment for storing splines
        if(communication_type.eq.shmem_comm) then
           call aims_shm_init((l_pot_max+1)**2 * 2 * n_hartree_grid * n_atoms * 8, allocation_info)
           if(allocation_info /= 0) then
              call aims_stop('Preparation of Hartree potential: Error allocating shared memory.', func)
           endif
        endif

        ! other physics allocations (DO NOT CONTAIN n_full_points!!)
        call allocate_physics ( )

        call allocate_multipole_moments ( )

!       Allocate arrays for forced occupation
        if (force_occupation_projector) then
          call allocate_previous_eigenvector
        end if

        if (use_kerker_preconditioner) then
	   call prepare_preconditioner
	end if

	if (solvent_method.eq.SOLVENT_MPB) then
	  !Initialize LPB solver
	  if (limit_l_max_SPE == 0) then
	  	limit_l_max_SPE = maxval(l_hartree)
		write(info_str,'(2X,A,I4)') "Using maximum angular momentum for SPE multipole expansion of l_max = ", limit_l_max_SPE
            	call localorb_info(info_str,use_unit,'(A)')
	  end if

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

	  if (correct_dielectric_decrement.and.decrement_kind == 2) then
	    call initialize_eps_s_func()
	  end if

	  if (.not.allocated(free_rho_superpos_copy)) allocate(free_rho_superpos_copy(n_spin,n_full_points))
	  if (.not.allocated(free_rho_gradient_superpos_copy)) allocate(free_rho_gradient_superpos_copy(3,n_spin,n_full_points))

	  free_rho_gradient_superpos_copy = 0d0
	  free_rho_superpos_copy = 0d0

	  if (dynamic_cavity) then
	    call dielecfunc_from_density(rho,dielec_func_mpb)
	    call dielecfunc_gradient_from_density_gradient(dielec_func_mpb,dielec_func_gradient_mpb,rho, rho_gradient) !*pi4_inv
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
	    call alphafunc_from_density( rho, rho_gradient, dielec_func_mpb,n_full_points)
	  else
	    free_rho_superpos_copy(1,:) = free_rho_superpos(:)
	    do i_coords=1, 3, 1
	      free_rho_gradient_superpos_copy(i_coords,1,:) = free_rho_gradient_superpos(i_coords,:)
	    end do
	    call alphafunc_from_density( free_rho_superpos_copy*pi4_inv, free_rho_gradient_superpos_copy,dielec_func_mpb, n_full_points)
	  end if
! 	  if (.not. solve_lpbe_only) then
! 	    call get_cut_f(cutf_glob)
! 	  end if
          call prepare_SPE_solver
        else
          call allocate_lpbsolver_dummy
        end if!solvent_method.eq.SOLVENT_MPB

!        if (force_new_int.eq.9) then
!           call update_hartree_partition_tab_p1 &
!                ( hartree_partition_tab, hartree_partition_deriv )
!        end if

        call get_local_ylm_tab ( l_hartree )

        ! End preloading time measurement
        call get_times(time_free_offset, clock_time_free_offset, tot_time_free_offset, tot_clock_time_free_offset)


        call get_timestamps(time_superpos, clock_time_superpos)

        call get_free_superpos_energy_p1 ( partition_tab, rho, rho_gradient, kinetic_density, pot_ion_embed, &
            hartree_energy_free, en_elec_free, en_xc, en_pot_xc, en_ion_ion, &
            en_ion_embed, en_density_embed, en_vdw, en_ll_vdw, en_ll_vdw_err, &
            en_lda_c, en_pbe_c, ionic_forces, ext_charge_forces )
        if (out_nuc_nuc_repulsion) then
          write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') "| Nuc-Nuc repulsion energy      :", en_ion_ion, &
                " Ha", en_ion_ion * hartree, " eV"
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        end if

        call get_times(time_superpos, clock_time_superpos, tot_time_superpos, tot_clock_time_superpos)

        ! Compute total energy of first iteration

        call get_timestamps(time_etot, clock_time_etot)

!       if finite-width smearing was used, must get entropy correction for T=0 total energy
        if ((.not. use_elsi_dm) .and. (.not. elsi_read_dm) .and. (.not. fo_finalise)) then
           call get_entropy_correction_p1(KS_eigenvalue, occ_numbers, chemical_potential, chemical_potential_spin, entropy_correction)

           call get_penalty_energy(KS_eigenvector, KS_eigenvector_complex, occ_numbers, condition_penalty, penalty_energy)
        else
           entropy_correction = 0.0d0
           penalty_energy = 0.0d0
        end if

        if(.not. elsi_read_dm) then
           call get_total_energy ( ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, en_ion_embed, en_density_embed, &
             en_vdw, en_pot_vdw, en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, en_elec_free, en_elec_delta, fock_energy, &
             hartree_energy_free, hartree_delta_energy, hartree_multipole_correction, entropy_correction, penalty_energy, total_energy )
        endif

        ! make a note that the matrices are integrated already at least once
        first_integration = .false.


        call get_times(time_etot, clock_time_etot, tot_time_etot, tot_clock_time_etot)

	! Symmetry initialization before HF
        if ( use_symmetry_reduced_spg ) then
          call get_sym_overlap()
        endif
!       For Hartree-Fock calculation, the initialization is done here.
!       We hereby construct the product basis, calculate its overlap integrals
!       with products of individual basis functions, and the Coulomb interaction
!       matrix for product basis functions.
        if (use_hartree_fock) then
           if(use_hf_realspace)then
              call initialize_hartree_fock_p0()
           else
              call initialize_hartree_fock()
           endif
        endif

        if (out_density) then
          if (.not.use_mpi) then

                call output_density_p1( rho, partition_tab )
          else
            write(use_unit,'(1X,A)') '* Output of density is not supported in MPI runs. Aborting.'
            stop
          end if
        end if

        if (out_v_eff) then
           if (.not.use_mpi) then

              call output_potential_p1( hartree_potential )

           else
              write(use_unit,*) '* Output of potential is not supported'
              write(use_unit,*) '* in MPI runs. Aborting.'
              stop
           end if
        endif

        if (out_v_eff_new) then
!         write current effective potential on integration grid to file
!         If we are using a gradient functional, we are only writing
!         the local part of the potential, for now.

!         FIXME - THIS IS NOW ONLY THE HARTREE POTENTIAL, THE XC CONTRIBUTION IS NEVER ADDED ANY MORE.
!         MUST ADD THIS IN OUTPUT_POTENTIAL IF NECESSARY!!
           call output_potential_p1_new( hartree_potential,"initialize")
       end if

      end if

! VB: Since all matrices etc will be safely allocated here, we may here overwrite the Hamiltonian
!     if anyone required the nuclear potential matrix instead.
      if (out_nuc_matrix) then
        call integrate_v_external_matrix_p2 ( partition_tab, l_shell_max, hamiltonian )

        call output_v_external_matrix ( hamiltonian )
      end if

      if (out_t_plus_v_matrix) then
        hamiltonian = 0
        call integrate_t_plus_v_matrix_p2 ( partition_tab, l_shell_max, hamiltonian )

        call output_v_external_matrix ( hamiltonian )
      end if

      ! finally a setting for Pulay mixing:
      ! Since this is the initialization, the superposition of atomic densities
      ! charge density will be included in the mixing procedure
      first_iter_mixing = .true.

      ! This call initializes and allocates the multipole components of the difference between the
      ! sum of free atom densities and the actual density as it evolves through the s.c.f. cycle.
      ! The array is first used in update_hartree_potential and could also be allocated as late as there.
      ! Allocating this array as late as possible will keep memory free for temporary use on many platforms.
      call initialize_hartree_potential_storage()

      mpb_converged = .True.
      if (solvent_method.eq.SOLVENT_MPB) then
	!evaluate gradient of superposition of free atom potentials to regularize the full electrostatic potential
        call evaluate_v_free_gradient
        !we set mpb_converged = .false.. the SCF cycle will not converge until mpb_solver_started
	mpb_converged = .False.
	mpb_solver_started = .false.
	initialized_lpb_solver = .false.
	!see prepare_mpb_solver for how it is dealt with all the logicals
      end if

      ! SCF initialization ends here
      call get_times(time_sc_loop, clock_time_sc_loop)

      number_of_scf_inits = number_of_scf_inits + 1

      write(info_str,'(2X,A)') ""
      call localorb_info(info_str)
      call output_timeheader(deffmt, "End scf initialization - timings")
      call output_times(deffmt, "Time for scf. initialization", time_sc_loop, clock_time_sc_loop)
      call output_times(deffmt, "Boundary condition initialization", time_bc_init, clock_time_bc_init)
      if (use_geo_relaxation) then
         call output_timer("Relaxation initialization", time_relax_init(1:2), deffmt)
      end if
      call output_times(deffmt, "Integration", time_intg, clock_time_intg)
      if (.not. elsi_read_dm) then
         call output_times(deffmt, "Solution of K.-S. eqns.", time_diag, clock_time_diag)
      end if
      call output_times(deffmt, "Grid partitioning", time_partition, clock_time_partition)
      call output_times(deffmt, "Preloading free-atom quantities on grid", time_free_offset, clock_time_free_offset)
      if (force_potential.eq.0) then
         call output_times(deffmt, "Free-atom superposition energy", time_superpos, clock_time_superpos)
      end if
      if ((force_potential.eq.0) .and. (.not. elsi_read_dm)) then
         call output_times(deffmt, "Total energy evaluation", time_etot, clock_time_etot)
      end if
      if (use_hartree_fock) then
         call ovlp3fn_timings(deffmt)
      endif
      write(info_str,*)
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      call aims_mem_current_output ()
      write(info_str,'(A)') "------------------------------------------------------------"
      call localorb_info(info_str)
      if(output_level .ne. 'full') output_priority = OL_norm

    if (solvent_method.eq.SOLVENT_MPB) then
      if (allocated(free_rho_gradient_superpos_copy)) deallocate(free_rho_gradient_superpos_copy)
      if (allocated(free_rho_superpos_copy)) deallocate(free_rho_superpos_copy)
    end if

    ! Allocatable arrays that are tracked
    if (allocated(ham_ovlp_work)) call aims_deallocate( ham_ovlp_work, "ham_ovlp_work" )

    end subroutine initialize_scf
!******
