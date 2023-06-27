!****s* FHI-aims/scf_solver
!  NAME
!    scf_solver
!  SYNOPSIS

    subroutine scf_solver &
    ( converged, enough_walltime_left )

!  PURPOSE
!  High-level wrapper around a single, converged scf calculation
!  This routine can only be called after a suitable initialization routine
!  for the scf cycle. Necessary conditions for running are:
!  * all grids to be set up
!  * an existing wave function
!  * an existing overlap matrix
!  * all free-atom quantities to be tabulated correctly on the grid
!
!  and all of that with the correct array dimensions
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
      use hartree_fock
      use localized_basbas
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
      use get_dielectric_function
      use hartree_fock_p0
      use calculate_fock_matrix_p0, only: evaluate_exchange_matr_realspace_p0
      use force_occupation, only: force_occ_pr_state_periodic, force_occ_pr_state, &
                                  previous_eigenvector, previous_eigenvector_complex, &
                                  force_occ_state_periodic, get_occupation_numbers_occ_p0, &
                                  force_occ_state, force_occ_reduce_subspace, &
                                  force_occ_pr_state2
      ! start_force_occ, n_force_occ, force_occupation_projector (in dimensions.f90)
      use calculate_fock_matrix_p0
      use octree_routines
      use vdwkernel
      use load_balancing
      use hartree_potential_storage
      use numerical_stress
      use pseudodata
      use constraint
      use hartree_potential_recip
      use analytical_stress
      use energy_density
      use constants
      use external_pressure_on_system
      use mbd_rsSCS_numerical_forces
      use plot !temporary, this should go elsewhere
      use mpe_interface, only: mpe_main, &
          mpe_solved
      use fodft, only: fodft_in_potential, fodft_out_potential, fodft_select_hab
      use mbd_std_wrapper, only: mbd_std_calculate, mbd_self_consistent, &
          mbd_first_step, mbd_scf_converged, run_numerical_pulay_forces
      use hirshfeld, only: eval_hirsh_vol_pulay_deriv_dm, use_pulay, &
          numerical_hirsh_pulay_forces, evaluate_free_atom_quantities, run_hirshfeld
      use xml_write, only: xml_open, xml_close, xml_elem, tostr
      use mbd_dev_wrapper, only: mbd_dev_calculate
      use python_interface, only: run_python_hook, python_hooks
      use lpb_solver_utilities, only: prepare_mpb_solver, &
          postprocess_mpb_solver, mpb_converged, forces_mpb_force_off, &
          mpb_forces_on, dielec_func_mpb, Gnonmf_forces, Gnonmf_forces_on, &
          evaluate_sc_Gnonmf, use_Gnonmf_forces, postprocess_mpb_solver, &
          mpb_converged, forces_mpb_force_off, mpb_forces_on, dielec_func_mpb, &
          Gnonmf_forces, Gnonmf_forces_on, evaluate_sc_Gnonmf, &
          use_Gnonmf_forces, mpb_solver_started
      use elsi_wrapper, only: eh_scf, aims_elsi_occ, aims_elsi_set_output, &
          aims_elsi_init_scf, aims_elsi_finalize_scf
      use aims_memory_tracking, only: aims_mem_current_output, aims_allocate, &
          aims_deallocate
      use density_matrix_evaluation, only: output_densmat
      use batch_statistics, only:  synchronize_batch_statistics
      use json_output, only: write_scf_iteration_to_json_log, &
          write_final_output_to_json_log, write_geometry_to_json_log, &
          close_json_log_file
      use boys, only: boys_sub_flag
      use rel_x2c_mod, only: scf_iteration
      use heat_flux, only: assemble_heat_flux, print_atomic_stress
      implicit none

!  ARGUMENTS

      logical, intent(OUT) :: converged
      logical, intent(OUT) :: enough_walltime_left

!  INPUTS
!    none
!  OUTPUT
!    o  converged -- did the scf cycle converged or not.
!    o  enough_walltime_left -- was there enough wall time for calculations
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

      ! imported variables
      integer :: i

      logical :: elec_converged
      logical :: force_converged
      logical :: stress_converged
      logical :: density_is_close_to_scf
      logical :: frozen_core_converged

      ! local variables

      real*8, dimension(:,:), allocatable :: Pulay_matrix

      integer :: task_distribution_method
      character*200 :: info_str

      integer :: sc_abandon_counter = 0

      logical :: below_it_limit, out_header
      logical :: control_update_requested = .false.
      logical :: abort_requested = .false. ! Checks whether user requests abort
      logical :: abort_and_postprocess = .false. ! Check whethehr user requested
                                                 ! abort and postprocess

      logical :: require_post_prc

      logical :: precondition_check
      logical :: info
      logical :: pulay_forces_backup

      logical :: gga_forces_backup
      logical :: meta_gga_forces_backup
      logical :: nlcc_forces_backup
      logical :: Gnonmf_forces_backup
      logical :: flag_calc_nonlocal = .false.

      real*8 :: time_postprocessing_start
      real*8 :: clock_ke_time, ke_time

      character*8  :: cdate
      character*10 :: ctime
      character(*), parameter :: deffmt = '2X'

      real*8 :: alpha
      ! Fermi loop variables for Kubo-Greenwood
      integer :: fermisteps, fermiindex
      real*8 :: fermilevel

      ! END fermi loop variables

      ! VB: The following four variables are ok for now, but should
      !     move to module physics in the future, for consistency
      real*8 ::  pot_ion_embed_nuc
      real*8, dimension(3) :: embedding_force
      real*8, dimension(3, n_atoms) :: vdw_forces
      real*8, dimension(3, n_atoms) :: ll_vdw_forces
      real*8, dimension(n_spin) :: kinetic_energy
      real*8, dimension(n_spin) :: x_energy
      real*8 :: x_post_energy
      real*8 :: xc_energy
      real*8 :: c_energy, c_post_energy
      real*8 :: c_energy_lda
      real*8 :: penalty_energy
      real*8, dimension(n_spin) :: x_energy_lda

      ! matrices to be used when use_local_index is on

      real*8, allocatable :: my_hamiltonian_matrix(:,:)
      real*8, allocatable :: my_density_matrix(:)

      real*8, allocatable :: volumes_hirsh(:), volumes_hirsh_free(:)
      real*8, allocatable :: alpha_0_vv_nm(:), C6_vv_nm(:)
      real*8, allocatable :: alpha_0_vv_free(:), C6_vv_free(:)

      real*8, allocatable :: mbd_rpa_orders(:)

      ! counters

      integer :: i_spin
      integer :: i_coords
      integer :: i_atom
      integer :: i_k_point,i_force_occ,i_basis,i_k,new_k_point
      integer :: count2
      integer :: AS_counter
      !AJL, Feb2018
      integer :: i_full_points

      !for density restart
      integer :: i_error

      integer :: i_tmp

      ! SVL array for storing current occ_numbers
      real*8, dimension(:,:,:), allocatable :: occ_numbers_save
      character(*), parameter :: func = 'scf_solver'

      ! CC: Energy density: Avoid recomputation of density gradient when
      !     updating density
      logical          :: use_density_gradient_save
      logical          :: use_meta_gga_save

      real*8  en_nlcorr !SAG
      real*8, dimension(3) :: coord_current
      integer i_my_batch, i_index
      logical use_hartree_non_periodic_ewald_save
      logical first_iteration

      ! for doubly-hybrid density functional
      real*8  total_energy_old, en_xc_old

      ! Exact exchange energy and its gradient
      real*8 :: exx_ene = 0.d0
      real*8 d_exx_ene(3, n_atoms)
      logical :: EXX_forces_on = .false.

      ! Igor, variables for dft component print out
      integer     :: i_dft
      integer     :: flag_xc_old
      real*8      :: hybrid_coeff_old
      logical     :: use_meta_gga_old
      character*8 :: dft_flag
      integer :: n_states_occupied

      ! VVG
      real*8 :: ene_pairwise_ts      !<---C6/R^6 using TS alpha and C6
      real*8 :: ene_pairwise_scs     !<---C6/R^6-using SCS alpha and C6
      real*8 :: ene_mbd_scs          !<---MBD@SCS
      real*8 :: ene_mbd_ts           !<---MBD@TS
      real*8 :: ene_mbd_rsSCS        !<---MBD@rsSCS

      integer :: exc_code
      character(50) :: exc_origin
      character(150) :: exc_msg

      ! exact exchange restart
      character*40 :: exx_restart

      logical :: Gnonmf_forces_calc

      real*8 :: tmp

      logical :: use_density_matrix_save

      integer :: update_hartree_potential_p1_counter

      ! begin work

      ! SVL allocate the temporary array for occ_numbers here
      if(use_periodic_hf)&
           allocate(occ_numbers_save(n_states, n_spin, n_k_points))

      !SVL: restart
      if (keep_restart_info.and.restart_write.and.use_hf_realspace) then
         write(exx_restart,'(A,A)') &
              'exx_',trim(restart_write_file)
      endif

      number_of_loops = 0
      ! The initial loop calculates the density - using the full scf solver, but
      ! it is not an actual iteration.
      if (restart_zero_iteration) number_of_loops = -1

      below_it_limit = (number_of_loops.lt.sc_iter_limit)

      ! Don't compute forces in the stress part when starting numerical stress
      if (use_numerical_stress.and.  &
           counter_numerical_stress .eq. 0 .and. &
           .not. num_stress_finished) then
         original_force_flag             = use_forces
         original_analytical_stress_flag = use_analytical_stress
         use_forces                      = .false.
         use_analytical_stress           = .false.
      end if

      elec_converged = .false.

      if (use_forces) then
        force_converged = .false.
      else
        force_converged = .true.
      end if

      if ( use_analytical_stress ) then
        stress_converged = .false.
      else
        stress_converged = .true.
      end if
      AS_counter = 0

      ! check if there is time for at least one scf cycle
      enough_walltime_left = time_check_scf_cycle()

      out_header = .false.

      ! initialize forces to zero, just in case
      vdw_forces = 0d0
      ll_vdw_forces = 0d0

      if (use_kerker_preconditioner) then
        kerker_preconditioner_on = .true.
        ! Store actual charge mixing parameter and use special one for Kerker
        ! instead.
        stored_charge_mix_param = charge_mix_param(1)
        charge_mix_param(1) = prec_mix_param(1)
      end if

      if (force_occupation_projector.and.n_periodic>0) then
        do i_k_point = 1, n_k_points, 1
          force_occ_pr_state_periodic(1:n_force_occ,i_k_point) = &
              force_occ_pr_state(1:n_force_occ)
        end do
      end if
      if(output_level .eq. 'MD_light') output_priority = 2

     if(use_vdw_method.or.use_vdw_post.or.use_nlcorr_post)then
        allocate(root_full)
        nullify(root_full%branch)
        call SetupKernel   ! Removed at end of subroutine. SAG
     endif

     use_hartree_non_periodic_ewald_save = use_hartree_non_periodic_ewald

     first_iteration = .true.
     if (lda_update_switch) then
        use_density_matrix = .true.
     end if

     if (use_fo_potential) then
        if (allocated(local_fo_potential))  deallocate(local_fo_potential)
        allocate(local_fo_potential(n_full_points))

        local_fo_potential = 0.d0
        if (fo_potential_file) then
          ! Read in potential from file!
          call fodft_in_potential( local_fo_potential )
        end if
     else
        if (allocated(local_fo_potential))  deallocate(local_fo_potential)
        allocate(local_fo_potential(1))
     end if !(use_fo_potential)

     mbd_first_step = .true.
     mbd_scf_converged = .false.

     if (use_forces.and.solvent_method.eq.SOLVENT_MPB.and.&
         .not.forces_mpb_force_off.and.evaluate_sc_Gnonmf) then
        use_Gnonmf_forces = .true.
     else
        use_Gnonmf_forces = .false.
     end if
     Gnonmf_forces_on = .false.
     if (restart_zero_iteration) mbd_first_step = .false.

     if (use_meta_gga) then
          write(info_str,'(2X,A)') ''
          call localorb_info(info_str,use_unit,'(A)', OL_norm)
          write(info_str,'(2X,A)') &
                  "Initialising kinetic-energy density for meta-GGA"
          call localorb_info(info_str,use_unit,'(A)', OL_norm)

          call get_timestamps (ke_time, clock_ke_time)

          call calculate_kinetic_density &
            (rho, rho_gradient, kinetic_density, &
            hartree_partition_tab, partition_tab, l_shell_max, &
            KS_eigenvector, KS_eigenvector_complex, occ_numbers)

          call get_timestamps (rtime, clock_rtime)
          !call output_timeheader(deffmt, &
          !                       "Relaxation / MD: End force evaluation.", &
          !                       OL_high)
          call output_times(deffmt, "Time to calculate kinetic-energy density",&
                            rtime-ke_time, clock_rtime-clock_ke_time, OL_norm)
          !  here may be expanded to sub timings
          !write(info_str,'(2X,A)') ''
          !call localorb_info(info_str,use_unit,'(A)', OL_norm)
     endif

! AL: Necessary if we are using a non-hybrid (e.g. meta-GGA) after
!     pre-convergence. It would be tidier to have a replica infrastructure
!     for pre-convergence, whereby we could use all functionals with all
!     functionals. One for the to-do list
!
!     Pre-convergence is only available with non-hybrids, hence the
!     hybrid_coeff is hardcoded at present
     if (flag_xc_pre > 0) then
!         !flag_xc_old = flag_xc
         hybrid_coeff_old = hybrid_coeff
!         !flag_xc = flag_xc_pre
         hybrid_coeff = 0.d0
     end if

      update_hartree_potential_p1_counter = 0

!===============================================================================
!=                          BEGIN SELF-CONSISTENCY LOOP                        =
!===============================================================================

  SCF_LOOP: do while ( (.not.converged) .and.  &
  &                    below_it_limit .and. &
  &                    enough_walltime_left .and. &
  &                    (.not.restarting_scf) .and. &
  &                    (.not.abort_requested) )
        if(myid == 0) flush(use_unit)

        number_of_loops = number_of_loops + 1

        if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)scf_iteration=number_of_loops

        ! ELSI restart
        if(elsi_write_dm) then
           if(mod(number_of_loops-1,elsi_write_dm_freq) == 0) then
              elsi_write_dm_this = .true.
           else
              elsi_write_dm_this = .false.
           endif
        endif

        if (number_of_loops > 1) mbd_first_step = .false.

        if (solvent_method.eq.SOLVENT_MPB.and.forces_on) then
          ! If forces_on = true, we have to calculate the HF force corrections
          ! for the MPB solver during this SCF step
          mpb_forces_on = .true.
        else
          mpb_forces_on = .false.
        end if

        ! DFT+U ramping: The U value is increased stepwise until the
        ! defined target value is reached.
        ! This only affects the hubbard_u value, nothing more.
        if(plus_u_ramping_defined .eqv. .true.) then
            call plus_u_ramping_work(total_energy, previous_total_energy)
        endif

        ! If we want to do a Boys-localization at the beginning or throughout the SCF, do it here
        if (apply_boys_flag) then
          if (((boys_sub_flag.eq.0).and.(number_of_loops.eq.1).and.(.not.restart_zero_iteration)).or.(boys_sub_flag&
          &.eq.1).or.((boys_sub_flag .eq.0).and.(number_of_loops.eq.0).and.(restart_zero_iteration))) then
            call get_dipolematrix(KS_eigenvalue, KS_eigenvector, &
                                  KS_eigenvector_complex, occ_numbers, &
                                  chemical_potential, partition_tab, &
                                  l_shell_max)
            ! Unset flag not to trigger a second time after a restart run
            if (boys_sub_flag.ne.1) then
              boys_sub_flag = -1
            end if
          end if
        end if

        ! Before anything else is done -- output the cube files as they are
        ! right now.
        ! WHY is that here and not at the end of the loop? The advantage is that
        ! this allows to plot both the initial density.
        ! if (out_cube_nth_iteration) then
        !    if (mod(number_of_loops,cube_output_nth_iteration)==0) then
        !       call setup_cube_defaults ! needs to be called every time, to make
        !                                ! sure a new filename is generated
        !       call output_cube_files_p2
        !    endif
        ! endif

        write(info_str,'(A)') ''
        call localorb_info ( info_str,use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "------------------------------------------------------------"
        call localorb_info ( info_str,use_unit,'(A)', OL_norm  )
        if (restart_zero_iteration .and. .not. (output_level .eq. 'MD_light')) then
          write(info_str,'(10X,A,1X,I4)') "Initial SCF-cycle to compute last known "
          call localorb_info ( info_str, use_unit )
          write(info_str,'(10X,A,1X,I4)') "density from restart information."
        else
          write(info_str,'(10X,A,1X,I4)') "Begin self-consistency iteration #", number_of_loops
        end if
        call localorb_info ( info_str,use_unit,'(A)', OL_norm )
        if (output_level .eq. 'MD_light') then  ! different output for one-line-per-iteration
           if (.not. out_header) then
              if (n_spin.eq.1) then
                 call localorb_info( &
                 'Convergence:    q app. |  density  | eigen (eV) | Etot (eV) ',&
                 use_unit,'(A,$)', OL_high )
              else
                 call localorb_info( &
                 'Convergence:    q app. |  density,  spin     | eigen (eV) | Etot (eV) ', &
                 use_unit,'(A,$)', OL_high )
              end if
              if (use_forces) then
                 call localorb_info(' | forces (eV/A)', use_unit,'(A,$)', OL_high)
              else
                 write(info_str,'(A)') ' |             .'
                 call localorb_info(info_str, use_unit,'(A,$)', OL_high)
              end if
              call localorb_info(' |       CPU time |     Clock time', use_unit,'(A)', OL_high)
              out_header = .true.
           end if
           if (restart_zero_iteration) then
              write(info_str,'(2X,A)') "last SCF : "
           else
              write(info_str,'(2X,A,I4,A)') "SCF ", number_of_loops,' : '
           end if
           call localorb_info ( info_str,use_unit,'(A,$)', OL_high )
        end if
        write(info_str,'(A)') ''
        call localorb_info ( info_str,use_unit,'(A)', OL_norm )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, ", Time     :  ", ctime
        call localorb_info ( info_str,use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "------------------------------------------------------------"
        call localorb_info ( info_str,use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )

        previous_ev_sum = ev_sum
        previous_total_energy = total_energy
        previous_pot_jump = pot_jump

        if (force_occupation_projector_redsub) then
          do i_force_occ = 1, n_force_occ
            call force_occ_reduce_subspace(number_of_loops, i_force_occ, 1)
          end do
        end if

        !if (force_occupation_projector.and.n_periodic.eq.0) then
        if (force_occupation_projector.and.n_periodic.eq.0) then
          previous_eigenvector(:,:,:,:) = KS_eigenvector(:,:,:,:)
        else if (force_occupation_projector) then
          if (real_eigenvectors) then
            previous_eigenvector(:,:,:,:) = KS_eigenvector(:,:,:,:)
          else
            previous_eigenvector_complex(:,:,:,:) = KS_eigenvector_complex(:,:,:,:)
          end if
        end if

        ! If an external perturbation is requested to break the symmetry of a
        ! system, enforce and / or switch it here
        if (use_atom_ref) then

          if (number_of_loops.eq.ext_cycle) then

            write(info_str,'(2X,A)') "External electric field switched off"
            call localorb_info(info_str, use_unit)
            homogeneous_field = 0.d0

            en_ion_embed=0.d0

            do i_atom=1, n_occ_atoms

            call embedding_potential(coords(1,i_atom),pot_ion_embed_nuc,&
                                      embedding_force)
            ! needs negative sign because the embedding potential already
            ! contains the negative charge of the electron.
            en_ion_embed = en_ion_embed  &
                - species_z(species(i_atom)) * pot_ion_embed_nuc
            enddo

          endif
        endif


        if(number_of_loops > 10 )then
          multipole_feed_back_parameter = 0.d0
        end if

!-------------------- Electron Density Update and Forces --------------------->>

        ! Here is where we calculate the change in the electron density arising
        ! from the current solution to the generalized Kohn-Sham equation, as
        ! well as(most of the) contributions to the forces.  Preconditioning and
        ! the electron density mixing will occur in the next code block.
        !
        ! Density:
        ! There are two approaches to calculating the change in electron density
        ! in FHI-aims:
        ! - Orbital-based:  This is the familiar sum over probablility densities
        !   of states approach which you learned in college.  It is O(N^2) and
        !   thus not suitable for large-scale calculations.  This update method
        !   is only used for relatively small cluster calculations in FHI-aims.
        ! - Density-matrix-based: This is the method where the density matrix
        !   (essentially an outer product of KS eigenvectors weighted by the
        !   occupation number) is constructed, and the density is evaluated via
        !   a batch integration scheme similar to the Hamiltonian matrix.  For
        !   sufficiently large calculations, this is an O(N) method, and it
        !   is used for all periodic calculations as well as larger cluster
        !   calculations.  Of these two approaches, this should be considered
        !   the "default" method; as it is not too much more expensive than
        !   the orbital-based method in the worst case (small molecules) and has
        !   much more favorable scaling.
        ! Note that, in both cases, we return the *difference* of the new
        ! electron density relative to the previous electron density.  This is
        ! done as the mixers are defined in terms of differences of electron
        ! density.
        !
        ! Forces:
        ! Forces are really expensive to calculate, and so we only want to
        ! calculate them when absolutely necessary.  Whether forces should be
        ! calculated during the lifetime of the calculation is set by the
        ! use_forces flag.  To cut down on calculation time, we only start
        ! calculating forces after the electronic structure has converged,
        ! creating a new convergence criteron to keep iterating the SCF cycle
        ! until the forces have converged to some value at the current geometry.
        ! (This should not be confused with relaxation, where we modify the
        ! geometry until the SCF-converged forces for the current geometry are
        ! essentially zero.)
        !
        ! The flag which controls whether forces should be calculated as part of
        ! *this* SCF iteration is forces_on, which will be set during the
        ! convergence check code block later.  Various contributions to the
        ! force also have their own version of this flag; for example,
        ! pulay_forces_on controls whether the Pulay forces should be calculated
        ! during this SCF iteration.
        !
        ! We will not calculate the total forces until we check for SCF
        ! convergence, via the subroutine call get_total_forces().  For
        ! information about what all the contributions to the forces are, and
        ! where they are calculated, see the header for get_total_forces(),
        ! which is uncharacteristically well-documented.
        !
        ! Most force components are calculated alongside the density, hence the
        ! combination of these two procedures into one code block (and one set
        ! of timings). Exceptions abound, and new force components are
        ! continuously added as new features are added to FHI-aims, so the best
        ! approach to figure out which force components are calculated is to
        ! inspect get_total_forces() .

        call get_timestamps(time_density, clock_time_density)

        ! Allocate arrays for electron density mixing (in the next code block)
        if ((number_of_loops.eq.1).or.restart_zero_iteration) then
          if ( mixer.eq.MIX_LINEAR ) then
            call allocate_delta_rho ( )
          else if ( mixer.eq.MIX_PULAY) then
            call allocate_pulay ( )
          else if ( mixer.eq.MIX_BROYDEN) then
            call allocate_broyden ( )
          end if
        end if

        ! SVL: Save the original occupation numbers for the periodic exact
        !      exchange.
        if(use_periodic_hf) occ_numbers_save = occ_numbers

        ! VB: An ugly kludge, for now. occ_numbers are reweighted with k weights
        !     if needed.
        !     This means that occ_numbers no longer has its physical meaning.
        !     This is a temp solution until the problem is really fixed, but at
        !     least, the output is now correct.
        !
        !     The correct way would be to multiply in the k_weights for each
        !     k-point further down.
        !     THIS SHOULD BE DONE IN THE FUTURE.
        if ((.not. use_elsi_dm) .and. (.not. elsi_read_dm)) then
           call kweight_occs('scf_solver:charge', occ_numbers)
        end if

        ! Notice that rho and rho_gradient are only input here, from the
        ! previous generation.
        ! rho and rho_gradient are NOT updated by update_density!!!
        ! Instead, the routine updates the density differences delta_rho_KS and
        ! rho_diff, which are declared in the module mixing.f but they're public.
        ! rho and rho_gradient themselves will only be updated after the mixing
        ! is done.
        if (solvent_method.eq.SOLVENT_MPB) then
          !set some convergence flags for the SMPB solver
          call prepare_mpb_solver(elec_converged,number_of_loops,rho,rho_gradient)
        end if

        if(elsi_read_dm) then
           use_density_matrix_save = use_density_matrix
           use_density_matrix = .true.
        end if

        if(.not. use_density_matrix)then
          ! VB: It is important to note that this case is only ever true
          !     for non-periodic systems. Periodic systems always rely on
          !     use_density_matrix=.true. in FHI-aims.

          if(split_updates) then ! do not calculate forces
             pulay_forces_backup = pulay_forces_on
             gga_forces_backup = gga_forces_on
             meta_gga_forces_backup = meta_gga_forces_on
             nlcc_forces_backup = nlcc_forces_on
             Gnonmf_forces_backup = Gnonmf_forces_on
             gga_forces_on = .false.
             meta_gga_forces_on = .false.
             pulay_forces_on = .false.
             nlcc_forces_on = .false.
             Gnonmf_forces_on = .false.
          end if
          ! VB: Switched the input criterion forces_on to pulay_forces_on
          !     instead.
          call update_density_and_forces_orbital &
                ( KS_eigenvector, KS_eigenvector_complex, &
                KS_eigenvalue, occ_numbers, &
                partition_tab, hartree_partition_tab, rho, rho_gradient, &
                kinetic_density, hartree_potential, l_shell_max, &
                delta_rho_KS, delta_rho_gradient, delta_kinetic_density, rho_change, &
                hellman_feynman_forces, pulay_forces, &
                gga_forces, nlcc_forces, pulay_forces_on, &
                gga_forces_on, nlcc_forces_on, meta_gga_forces_on)

          ! Construct density matrix and write to file
          if(elsi_write_dm_this) then
             call output_densmat(KS_eigenvector,KS_eigenvector_complex,occ_numbers)
          endif

          if(split_updates) then
             gga_forces_on = gga_forces_backup
             meta_gga_forces_on = meta_gga_forces_backup
             pulay_forces_on = pulay_forces_backup
             nlcc_forces_on = nlcc_forces_backup
             Gnonmf_forces_on = Gnonmf_forces_backup
          end if
        end if

        if(fo_density) then
            call aims_allocate(my_density_matrix, n_hamiltonian_matrix_size, &
                               "my_density_matrix")
            call update_density_densmat &
            ( KS_eigenvector, KS_eigenvector_complex, &
            occ_numbers, &
            partition_tab,  hartree_partition_tab, rho, rho_gradient, kinetic_density, &
            l_shell_max, delta_rho_KS, delta_rho_gradient, delta_kinetic_density, &
            rho_change, my_density_matrix &
            )

        elseif(use_density_matrix .or. split_updates)then
          ! WPH: Why doesn't split_updates support load balancing?
          if(.not. split_updates) then
             if(use_local_index) then
               if(use_load_balancing) then
                 if(first_iteration) then
                   ! We use the batch permutation from the initial integration,
                   ! which should be better than the default batch permutation,
                   ! to calculate the weights for the density calculation
                   get_batch_weights = .true.
                   use_batch_permutation = n_bp_integ
                 else
                   ! Use the already-calculated batch permutation for density
                   get_batch_weights = .false.
                   use_batch_permutation = n_bp_density
                 endif
                 call aims_allocate(my_density_matrix, &
                                    batch_perm(use_batch_permutation)%n_local_matrix_size,&
                                   "my_density_matrix" )
               else
                 call aims_allocate(my_density_matrix, &
                                    n_hamiltonian_matrix_size, &
                                    "my_density_matrix")
               endif

               call update_density_densmat &
                     ( KS_eigenvector, KS_eigenvector_complex, &
                     occ_numbers, &
                     partition_tab, hartree_partition_tab, rho, rho_gradient, &
                     kinetic_density, l_shell_max, delta_rho_KS, delta_rho_gradient, &
                     delta_kinetic_density, rho_change, my_density_matrix  &
                     )

               if(use_load_balancing .and. first_iteration) then
                 ! Use the weights that we calculated to construct the batch
                 ! permutation
                 call compute_balanced_batch_distribution(n_bp_density)
               end if
               ! Return the calculation back to normal and exit.
               get_batch_weights     = .false.
               use_batch_permutation = 0
               call aims_deallocate( my_density_matrix, "my_density_matrix" )
             else
               call update_density_densmat &
                     ( KS_eigenvector, KS_eigenvector_complex, &
                     occ_numbers, &
                     partition_tab,  hartree_partition_tab, rho, rho_gradient, &
                     kinetic_density, l_shell_max, delta_rho_KS, delta_rho_gradient, &
                     delta_kinetic_density, rho_change, hamiltonian(1,1) &
                     )

             endif
          endif !.not.split_updates
          if( forces_on .and. pulay_forces_on .and. &
              ( (.not.use_analytical_stress) .or. (AS_counter.ge.1) ) ) then
            ! The output of update_forces_densmat combines the Pulay, GGA,
            ! meta-GGA, and gnonmf contributions to the forces together
            ! into one variable.
            ! Here, we use the pulay_forces variable to hold the combined sum,
            ! since it is always initialized for a force calculation
            if(use_local_index) then
              if(use_load_balancing) then
                get_batch_weights = .false.
                if (use_analytical_stress) then
                  ! When calculating the analytical stress tensor, the
                  ! AS_dde_potential variable, as calculated by
                  ! sum_up_whole_potential_p1(), is used.  This variable is
                  ! defined on the integration grid for the current MPI task,
                  ! which is a function of the batches for the current MPI task.
                  ! Thus we must use the same batch permutation for both the
                  ! forces/analytical stress tensor and the Hartree multipole
                  ! summation.
                  use_batch_permutation = n_bp_hpot
                else
                  ! We use the same batch permutation as the density, as their
                  ! workload is similar
                  use_batch_permutation = n_bp_density
                end if
                call aims_allocate(my_density_matrix, &
                                   batch_perm(use_batch_permutation)%n_local_matrix_size*n_spin,&
                                  "my_density_matrix")
              else
                call aims_allocate(my_density_matrix, &
                                   n_hamiltonian_matrix_size*n_spin,&
                                   "my_density_matrix" )
              end if
              call update_forces_densmat &
                      ( KS_eigenvector, KS_eigenvector_complex, &
                        KS_eigenvalue, occ_numbers, &
                        partition_tab, &
                        rho, rho_gradient, kinetic_density, hartree_potential, &
                        l_shell_max,  &
                        my_density_matrix, pulay_forces )

              get_batch_weights = .false.
              use_batch_permutation = 0
              call aims_deallocate( my_density_matrix, "my_density_matrix" )
            else
              call update_forces_densmat &
                      ( KS_eigenvector, KS_eigenvector_complex, &
                        KS_eigenvalue, occ_numbers, &
                        partition_tab, &
                        rho, rho_gradient, kinetic_density, hartree_potential, &
                        l_shell_max,  &
                        hamiltonian(1,1), pulay_forces )
            endif
            ! Now zero out all force contributions which are already included
            ! in pulay_forces when using the density-matrix based approach
            if (use_gga) gga_forces = 0.d0
            if (use_Gnonmf_forces) Gnonmf_forces = 0d0
          else if (use_forces) then
            ! We're not calculating force contributions in this iteration, zero
            ! out the relevant components
            pulay_forces = 0.d0
            if (use_gga) gga_forces = 0.d0
            if (use_Gnonmf_forces) Gnonmf_forces = 0d0
            if (use_embedding_pp.and.use_nonlinear_core) nlcc_forces = 0.0d0
          end if
        end if !use_density_matrix

        ! ELSI restart
        if(elsi_read_dm) then
           elsi_read_dm = .false.
           use_density_matrix = use_density_matrix_save
        end if

        if(restart_from_density) then
          open(file = "CMH-Density.out", unit = 27, status = 'old', &
               form = 'unformatted', iostat = i_error)
          close(27)
          if (i_error .eq. 0.and.(freeze_restart_density)) then
            delta_rho_KS = 0.0D0
            delta_rho_gradient = 0.0D0
            delta_kinetic_density = 0.d0
            rho_change = 0.0D0
          endif
        endif

        ! Inserted for PIMD wrapper - if code is calculating different replicas,
        ! start from superposition of free atom densities - need to do it also
        ! in linear mixing? (I think this is sufficient)
        if (need_dens_superpos) then
           delta_rho_KS = 0.0D0
           delta_rho_gradient = 0.0D0
           delta_kinetic_density = 0.d0
           rho_change = 0.0D0
           need_dens_superpos=.false. ! only do this in the first pass, then not
                                      ! anymore
        endif
        ! End insertion for PIMD wrapper

        call look_at_rho(rho, delta_rho_KS, partition_tab)

        call get_times(time_density, clock_time_density, tot_time_density, &
                       tot_clock_time_density)

!        ! Note: when using include_kinetic_density_in_mixing, currently some of the regression tests
!        ! fail with differences in energy, so the implementation isn't perfect still
!        if (use_meta_gga.and.use_density_matrix) then
!          write(info_str,'(2X,A)') ''
!          call localorb_info(info_str,use_unit,'(A)', OL_norm)
!          write(info_str,'(2X,A)') &
!               "Calculating kinetic-energy density for meta-GGA"
!          call localorb_info(info_str,use_unit,'(A)', OL_norm)
!
!          call get_timestamps (ke_time, clock_ke_time)
!
!          call calculate_kinetic_density &
!            (rho, rho_gradient, new_kinetic_density, &
!            hartree_partition_tab, partition_tab, l_shell_max, &
!            KS_eigenvector, KS_eigenvector_complex, occ_numbers)
!
!          ! AJL: I don't like having this here, but it works as proof of concept.
!          ! Updating the infrastructure to calculate delta_kinetic_density (see below)
!          ! would remove the need for this diff calculation, and in fact this entire if statement.
!          !if (include_kinetic_density_in_mixing.and.use_density_matrix) then
!            !do i_spin = 1, n_spin, 1
!              do i_full_points = 1, n_full_points, 1
!                ! FIXME: AJL, Feb2018: This isn't quite correct, yet, but Sheng Bi has done an excellent job in
!                ! helping put this together in a way that is manageable. You'll see below
!                ! that the delta_kinetic_density isn't calculated as a function of spin:
!                ! differences in both \alpha and \beta channel are with respect to just the
!                ! \alpha channel of previous iterations. This points to a flaw in the current calculation
!                ! of the kinetic_density, somehow, in the beta channel, but in a way that doesn't effect
!                ! end energies (as these all match implementations in other codes, and the gradients etc. are consistent)
!                !
!                ! It would be ideal to remove the need to calculate kinetic density difference here at all,
!                ! but instead calculate delta_kinetic_density when we calculate delta_rho_KS and delta_rho_gradient.
!                ! This must be on the to-do list, as for now at least we have a working implementation
!                if (spin_treatment.eq.0) then
!                  delta_kinetic_density(i_full_points,1) = &
!                  (new_kinetic_density(1,i_full_points) - kinetic_density(1,i_full_points))
!                else
!                  delta_kinetic_density(i_full_points,1) = &
!                  new_kinetic_density(1,i_full_points) + new_kinetic_density(2,i_full_points) &
!                  - (kinetic_density(1,i_full_points) + kinetic_density(2,i_full_points))
!                  delta_kinetic_density(i_full_points,2) = &
!                  new_kinetic_density(1,i_full_points) - new_kinetic_density(2,i_full_points) &
!                  - (kinetic_density(1,i_full_points) - kinetic_density(2,i_full_points))
!                endif
!              end do
!            !end do
!          !endif
!
!          call get_timestamps (rtime, clock_rtime)
!          !call output_timeheader(deffmt, &
!          !                       "Relaxation / MD: End force evaluation.", &
!          !                       OL_high)
!          call output_times(deffmt, "Time to calculate kinetic-energy density",&
!                            rtime-ke_time, clock_rtime-clock_ke_time, OL_norm)
!          !  here may be expanded to sub timings
!          write(info_str,'(2X,A)') ''
!          call localorb_info(info_str,use_unit,'(A)', OL_norm)
!        endif

! AJL Debugs
!        write(use_unit,*) 'IN KE = ', kinetic_density(:,1)
!        write(use_unit,*) 'IN dKE = ', delta_kinetic_density(1,:)

!---------------- Preconditioner and Electron Density Mixing ----------------->>

        ! Here is where we take the change in the electron density arising from
        ! the current solution to the generalized Kohn-Sham equation, (possibly)
        ! apply a preconditioner to it, and mix it with previous densities,
        ! generating the *actual* density which will be used for calculating
        ! subsequent quantities.

        call get_timestamps(time_mixing, clock_time_mixing)

        ! If preconditioning is suggested /before/ mixing, it should be
        ! applied to delta_rho_KS (and to delta_rho_gradient, if applicable).
        if (use_kerker_preconditioner .and. &
        &   kerker_preconditioner_on .and. &
        &   precondition_before_mixer) then
           do i_spin = 1, n_spin
              call precondition_kerker(delta_rho_KS(:, i_spin), &!delta_rho_KS
              &                        hartree_partition_tab)

              if (use_density_gradient) then
                 do i_coords = 1, 3
                    call precondition_kerker( &
                    &         delta_rho_gradient(i_coords,:, i_spin), &
                    &         hartree_partition_tab)
                 end do
                if(use_meta_gga) then
                  call precondition_kerker(delta_kinetic_density(:, i_spin), &
                  &                        hartree_partition_tab)
                endif
              end if
           end do
        end if

        ! NO mixing at all for restarting information - need pure information.
        if (restart_zero_iteration) then
          if (n_spin .eq. 1) then
            i_spin = 1
            rho(i_spin,:) = rho(i_spin,:) + delta_rho_KS(:,i_spin)
            if (use_density_gradient) then
              do i_coords = 1, 3
                rho_gradient(i_coords,i_spin,:) = &
                     rho_gradient(i_coords,i_spin,:) + &
                     delta_rho_gradient(i_coords,:,i_spin)
              end do
              if (use_meta_gga) then
                kinetic_density(i_spin,:) = kinetic_density(i_spin,:) + delta_kinetic_density(:,i_spin)
              end if
            end if
          else
            rho(1,:) = rho(1,:) + 0.5d0*(delta_rho_KS(:,1)+delta_rho_KS(:,2))
            rho(2,:) = rho(2,:) + 0.5d0*(delta_rho_KS(:,1)-delta_rho_KS(:,2))
            if (use_density_gradient) then
              do i_coords = 1, 3
                  rho_gradient(i_coords,1,:) = rho_gradient(i_coords,1,:) + &
                      0.5d0 * (delta_rho_gradient(i_coords,:,1) + delta_rho_gradient(i_coords,:,2))
                  rho_gradient(i_coords,2,:) = rho_gradient(i_coords,2,:) + &
                      0.5d0 * (delta_rho_gradient(i_coords,:,1) - delta_rho_gradient(i_coords,:,2))
              end do
              if (use_meta_gga) then
                kinetic_density(1,:) = kinetic_density(1,:) + 0.5d0*(delta_kinetic_density(:,1)+delta_kinetic_density(:,2))
                kinetic_density(2,:) = kinetic_density(2,:) + 0.5d0*(delta_kinetic_density(:,1)-delta_kinetic_density(:,2))
              end if
            end if
          end if

        else if (mixer.eq.MIX_LINEAR) then
          ! linear mixing
          call linear_mix_p1(partition_tab, hartree_partition_tab, rho, &
                             rho_gradient,kinetic_density)
        else if (mixer.eq.MIX_PULAY) then
          ! if Pulay mixing requested, take care of all storage and / or mixing

          if (first_iter_mixing.or.(number_of_loops.ge.2)) then
            if (number_of_loops.le.ini_linear_mixing) then
              ! i.e. if Pulay mixing was requested later in the run,
              ! Must still perform linear mixing at this point
              call linear_mix_p1(partition_tab, hartree_partition_tab, rho, &
                                 rho_gradient,kinetic_density)
              ! store everything that will be needed for Pulay mixing ...
              ! MUST BE CALLED AFTER LINEAR MIXING !!
              call prepare_pulay_mixing()
            else
              ! Pulay mixing of previous charge densities
              call execute_pulay_mixing(number_of_loops, partition_tab, &
                                        hartree_partition_tab, rho, &
                                        rho_gradient, kinetic_density)
            end if

          else if (number_of_loops.gt.ini_linear_mixing) then
            ! this is the first iteration in a run that is being restarted
            ! simply add the full delta_rho_KS to the present density
            ! pulay_update handles this task for us as usual, but
            ! no density differences from this iteration are stored
            rho_diff = delta_rho_KS
            if (use_density_gradient) then
              rho_gradient_diff = delta_rho_gradient
              if (use_meta_gga) then
                kinetic_density_diff = delta_kinetic_density
              end if
            end if

            call pulay_update_p1(partition_tab, rho, rho_gradient, kinetic_density)
          end if
        else if (mixer.eq.MIX_BROYDEN) then
          ! Broyden mixing of previous charge densities
          call execute_broyden_mixing(number_of_loops, partition_tab, &
                                      hartree_partition_tab, rho, rho_gradient, kinetic_density)
        end if

        if (compensate_multipole_errors) then
           ! If we are interested in having a density that has zero charge as
           ! integrated on the actual integration grid, we must enforce this
           ! zero here.
           ! The mixed density is composed of:
           ! (1) The new density from the Kohn-Sham equations, which
           !     was properly normalized on the 3D grid
           ! (2) The "sum of free atoms" density, which is NOT exactly
           !     normalized on the 3D integration grid
           ! (3) The preconditioned mixed density difference between the
           !     new Kohn-Sham density and the "sum of free atoms" density,
           !     where the preconditioner yet again breaks the normalization
           !     slightly.
           ! Hence, the most practical avenue is to renormalize here.
           !
           ! Keeping this call distinct from 'normalize_initial_density'.
           ! Since the Pulay mixer should mix exactly normalized densities
           ! already, very little renormalization should occur in the
           ! call below.
           !
           call renormalize_density(rho, rho_gradient, kinetic_density, &
                                    n_spin, partition_tab, n_electrons)
        end if

        call get_times(time_mixing, clock_time_mixing, tot_time_mixing, &
                       tot_clock_time_mixing)

        if (out_density) then
          if (.not.use_mpi) then
            call output_density_p1(rho, partition_tab)
          else
            write(info_str,'(1X,A,A)') '* Output of density is not', &
            & ' supported in MPI runs. Aborting.'
            call aims_stop(info_str, func)
          end if
        end if

! AJL Debugs
!        write(use_unit,*) 'OUT KE = ', kinetic_density(:,1)
!        write(use_unit,*) 'OUT dKE = ', delta_kinetic_density(1,:)

!---------------- Update the Electrostatic (Hartree) Potential --------------->>

        ! Update the Hartree potential using a multipole decomposition
        ! via update_hartree_potential_p1() and subsequent summation by
        ! sum_up_whole_potential_p1().  The calculation time for this section
        ! is dominated by sum_up_whole_potential_p1().

        ! just for test
        ! do i_tmp = 1, 10, 1

        call get_timestamps( time_hartree_multi, clock_time_hartree_multi)

        ! Multipole decomposition

        if (.not. use_rho_multipole_shmem) then
          call update_hartree_potential_p1 &
              ( hartree_partition_tab, free_rho_superpos, &
              rho, rho_multipole_old, &
              delta_v_hartree_part_at_zero, &
              delta_v_hartree_deriv_l0_at_zero, &
              multipole_moments, multipole_radius_sq, &
              l_hartree_max_far_distance, rho_radial_multipole_old, &
              outer_potential_radius , AS_stress_on, .true. &
              )
        else
          call update_hartree_potential_p1_without_old &
              ( hartree_partition_tab, free_rho_superpos, &
              rho, rho_multipole_old, &
              delta_v_hartree_part_at_zero, &
              delta_v_hartree_deriv_l0_at_zero, &
              multipole_moments, multipole_radius_sq, &
              l_hartree_max_far_distance, &
              outer_potential_radius , AS_stress_on, .true. &
              )
        endif


        !write(info_str,'(1X,A)') "YY multipole_moments"
        !call localorb_info(info_str,use_unit,'(A)', OL_norm)
        !do i_atom = 1, 8
        !  do i = 1, 5
        !    write(info_str,'(10E15.7)') multipole_moments((i-1)**2+1:(i)**2,i_atom)
        !    call localorb_info(info_str,use_unit,'(A)', OL_norm)
        !  end do
        !end do
        !write(info_str,*) 
        if (out_rho_multipole) then
          call output_rho_multipole
        end if

        if (out_v_hartree) then
          call output_hartree_pot ( multipole_moments )
        end if

        call get_times(time_hartree_multi, clock_time_hartree_multi, &
        &          tot_time_hartree_multi, tot_clock_time_hartree_multi)



        call output_times(deffmt, "Hartree multipole update", &
        &                 time_hartree_multi, clock_time_hartree_multi,OL_norm)

        ! just for test
        !   if(myid .eq. 0) print*, mpi_wtime()
        !   call sleep(30)
        !   if(myid .eq. 0) print*, mpi_wtime()
        !   if(myid .eq. 0) print*, "===================="
        ! end do
        ! call aims_stop("Just for test the time of update_hartree_potential_p1")

        ! update_hartree_potential_p1_counter = update_hartree_potential_p1_counter + 1
        ! if (update_hartree_potential_p1_counter > 6) call aims_stop("Just for test the time of update_hartree_potential_p1")

        ! Implicit solvation:
        ! now that the Hartree potential is updated, we can compute the
        ! reaction field, which is then added as an external potential
        ! in sum_up_whole_potential_p1
        if (solvent_method.eq.SOLVENT_MPE) then
           call mpe_main(number_of_loops)
        end if

        call get_timestamps( time_hartree_sum, clock_time_hartree_sum)

        ! On every integration grid point, tabulate the entire electrostatic
        ! potential.  Note that the exchange-correlation potential is handled
        ! separately during the integration itself.

        ! See density calculation section for comments on load-balancing
        ! calibration
        if(use_local_index .and. use_load_balancing) then
           if(first_iteration) then
              get_batch_weights = .true.
              use_batch_permutation = n_bp_integ
           else
              get_batch_weights = .false.
              use_batch_permutation = n_bp_hpot
           endif
        endif

        ! run LPBE/MPBE solver
        if (solvent_method.eq.SOLVENT_MPB.and.mpb_solver_started) then
          call run_mpb_solver(rho,rho_gradient,dielec_func_mpb)
        end if

        ! Multipole summation
        call sum_up_whole_potential_p1 &
                ( delta_v_hartree_part_at_zero, &
                delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
                partition_tab, rho, &
                free_hartree_superpos, free_rho_superpos,  &
                hartree_potential,  &
                hartree_delta_energy, en_elec_delta, &
                hartree_multipole_correction, &
                pot_ion_embed, en_density_embed, &
                multipole_forces, forces_on, multipole_radius_sq, &
                l_hartree_max_far_distance, hellman_feynman_forces, &
                energy_deriv_tress_tensor, rho_multipole_old, &
                outer_potential_radius, &
                AS_stress_on, .true., local_fo_potential )

        if(use_local_index .and. use_load_balancing .and. first_iteration) then
           call compute_balanced_batch_distribution(n_bp_hpot)
           if(use_distributed_spline_storage) then
              ! The distribution of rho_multipole over the tasks has to be
              ! reinitialized
              call reset_hartree_potential_storage &
                   ( batch_perm(n_bp_hpot)%n_my_batches, &
                     batch_perm(n_bp_hpot)%batches, &
                     batch_perm(n_bp_hpot)%partition_tab )
           endif
        endif
        get_batch_weights = .false.
        use_batch_permutation = 0

        call get_times(time_hartree_sum, clock_time_hartree_sum, &
        &              tot_time_hartree_sum, tot_clock_time_hartree_sum)

        if(use_elsi_dm .and. elsi_solver == 3) then ! PEXSI
           ! Estimate min and max delta V
           if(allocated(hartree_potential_save)) then ! SCF initialized
              hartree_potential_save = hartree_potential-hartree_potential_save

              tmp = minval(hartree_potential_save)
              call get_min_double(dv_hartree_min,tmp)

              tmp = maxval(hartree_potential_save)
              call get_max_double(dv_hartree_max,tmp)

              hartree_potential_save = hartree_potential
           else ! SCF reinitialized
              call aims_allocate(hartree_potential_save,n_full_points,&
                      "hartree_potential_save")

              hartree_potential_save = hartree_potential
           endif
        endif

        if (out_v_eff) then
           if (.not.use_mpi) then
              ! write current effective potential on integration grid to file
              call output_potential_p1 ( hartree_potential )
           else
              write(info_str,*) '* Output of potential is not supported'
              call localorb_info(info_str, use_unit)
              write(info_str,*) '* in MPI runs. Aborting.'
              call aims_stop(info_str, func)
           end if
        endif

        if (out_v_eff_new) then
           call output_potential_p1_new ( hartree_potential,"iteration")
        end if

        if(use_vdw_method.and.(.not.use_nlcorr_post) )then  !SAG
           if(associated(root_full%branch)) then
              count2 = 0
              call deallocate_tree(root_full,count2)
           endif
           call gen_master_cube()
        endif

!        ! AJL: This is saved currently just for continuity and debugging purposes
!        ! Once kinetic-mixing is accepted as standard, this loop can be removed.
!        ! Note: when using include_kinetic_density_in_mixing, currently some of the regression tests
!        ! fail with differences in energy, so the implementation isn't perfect still
!        if (use_meta_gga.and..not.include_kinetic_density_in_mixing.and.use_density_matrix) then
!          write(info_str,'(2X,A)') &
!                  "Calculating kinetic-energy density for meta-GGA"
!          call localorb_info(info_str,use_unit,'(A)', OL_norm)
!
!          call get_timestamps (ke_time, clock_ke_time)
!
!          call calculate_kinetic_density &
!          (rho, rho_gradient, kinetic_density, &
!          hartree_partition_tab, partition_tab, l_shell_max, &
!          KS_eigenvector, KS_eigenvector_complex, occ_numbers)
!
!          call get_timestamps (rtime, clock_rtime)
!          !call output_timeheader(deffmt, &
!          !                       "Relaxation / MD: End force evaluation.", &
!          !                       OL_high)
!          call output_times(deffmt, "Time to calculate kinetic-energy density",&
!                            rtime-ke_time, clock_rtime-clock_ke_time, OL_norm)
!          !  here may be expanded to sub timings
!          write(info_str,'(2X,A)') ''
!          call localorb_info(info_str,use_unit,'(A)', OL_norm)
!        endif

!------------------- Compute REAL SPACE Hamiltonian Matrix ------------------->>

        ! Construct the REAL SPACE Hamiltonian matrix via the batch integration
        ! scheme as outlined in Blum et al., Comput Phys Comm 2009 and Havu et
        ! al., J. Comput Phys 2009.
        ! An important point: this is NOT the Hamiltonian which enters into the
        ! eigensolver!  That will be constructed later, shortly before the
        ! solution of the generalized Kohn-Sham equation.
        !
        ! As of this writing (2018 January 7), there are four storage formats
        ! for the real-space Hamiltonian:
        ! - Non-packed:
        !   (use_local_index .false., packed_matrix_format PM_none)
        !   This is the naive approach where every MPI task has the full copy
        !   of the real-space Hamiltonian matrix elements <i|h|j> between all
        !   basis elements for the entire calculation.  This storage format is
        !   only used for small cluster calculations as it is MASSIVELY
        !   INEFFICIENT IN MEMORY.
        ! - Packed:
        !   (use_local_index .false., packed_matrix_format PM_index)
        !   In this storage format, every MPI task has a full copy of the
        !   real-space Hamiltonian matrix, but it is compressed via a CSR format
        !   to alleviate the memory overhead.  The header of the physics module
        !   explains this storage format far better than I (WPH) could here; the
        !   reader is encourage to look there for more information.  This is the
        !   default storage format for all periodic calculations and non-trivial
        !   cluster calculations.  While this storage format alleviates memory
        !   overhead, the non-distributed nature of this storage format will
        !   ultimately pose a memory bottleneck for sufficiently large-scale
        !   calculations.  For calculations which do not hit this memory
        !   bottleneck, this should be considered the default storage format.
        ! - Domain decomposition without load balancing:
        !   (use_local_index .true., use_load_balancing .false.)
        !   This storage format differs from the previous two formats in that
        !   the real-space Hamiltonian is distributed across MPI tasks based on
        !   the matrix elements that were computed by an MPI task during the
        !   batch integration scheme.  No MPI task has access to the full
        !   real-space Hamiltonian matrix at any point in the calculation.  This
        !   ensures that memory consumption is properly parallelized, however,
        !   there is a performance hit incurred when constructing the
        !   Hamiltonian matrix entering into the Kohn-Sham solver, as the
        !   various MPI tasks must communicate their pieces of the real-space
        !   Hamiltonian to to one another.  This storage method uses a CSR
        !   format locally and the same indexing arrays as the packed matrix
        !   format, though with differing array dimensions, so code using this
        !   storage format usually follows the same code path as the packed
        !   matrix format and should "just work" with little to no code
        !   restructuring.
        ! - Domain decomposition with load balancing:
        !   (use_local_index .true., load_balancing .true.)
        !   This storage format distributes the real-space Hamiltonian across
        !   MPI tasks using domain decomposition; however, the batches used in
        !   batch integration are also shifted around via pointers (yes,
        !   pointers in aims) to cut down on the performance cost associated
        !   with the domain decomposition.  This is the preferred approach for
        !   large-scale calculations when the packed matrix format has become
        !   too expensive in memory; however, there is some light coding that
        !   needs to be done to support this format.  First, the code must be
        !   made aware that the batches have been redistributed, and second,
        !   this matrix format is stored as a full matrix but only on the
        !   reduced set of basis elements contributing to matrix elements held
        !   by the current MPI task.

        ! Overlap matrix and grids are fixed - only recalculate the
        ! Hamiltonian matrix for given potential

        call get_timestamps(time_intg, clock_time_intg)

        if(use_local_index) then
          if(use_load_balancing) then
            ! If load balancing is in effect, there already exists a batch
            ! permutation for integrations calibrated during SCF initialization.
            ! However, this was from the overlap matrix calculation or was done
            ! for slightly different conditions, so we re-calibrate here during
            ! the first cycle.
            if (first_iteration) then
              get_batch_weights = .true.
            else
              get_batch_weights = .false.
            end if
            use_batch_permutation = n_bp_integ
          endif

          if (use_scalapack) then
            if (use_load_balancing) then
              call aims_allocate(my_hamiltonian_matrix, &
                                 batch_perm(use_batch_permutation)%n_local_matrix_size, &
                                 n_spin, "my_hamiltonian_matrix")
            else
              call aims_allocate(my_hamiltonian_matrix, &
                                 n_hamiltonian_matrix_size, n_spin, &
                                 "my_hamiltonian_matrix")
            end if
            call integrate_real_hamiltonian_matrix_p2 &
                  ( hartree_potential,  rho, rho_gradient, kinetic_density, &
                    partition_tab, l_shell_max, &
                    en_xc, en_pot_xc, my_hamiltonian_matrix, en_vdw, en_pot_vdw )
          else
            call integrate_real_hamiltonian_matrix_p2 &
                  ( hartree_potential,  rho, rho_gradient, kinetic_density, &
                    partition_tab, l_shell_max, &
                    en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )
          end if
          ! Normally this is where we'd calculate the new batch permutation, but
          ! in this case we don't do this until shortly before we solve the
          ! generalized KS equation for the first time.  I (WPH) don't know why
          ! we hold off on doing this, but there must be a good reason, possibly
          ! related to the IO subroutines?
        else
          call integrate_real_hamiltonian_matrix_p2 &
                ( hartree_potential,  rho, rho_gradient, kinetic_density, &
                  partition_tab, l_shell_max, &
                  en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )
        endif

        if(.not. use_local_index .and. packed_matrix_format == PM_index)then
          ! This here is for geometry relaxation, because then we need to
          ! reduce again. Otherwise it do not do anything.
           if(use_periodic_hf.or.use_DFPT_phonon.or.use_DFPT_dielectric.or. &
              use_DFPT_polarizability.or.use_DFPT_polar_reduce_memory)then
              info = .false.
           else
              call remove_small_numbers_in_hamiltonian_and_ovlp &
                 ( hamiltonian, overlap_matrix, info )
           endif
          if(info)then

            call reshape_matrices

            if (use_scalapack) then

              if(.not. use_elsi_dm) then
                if(n_periodic > 0 .or. packed_matrix_format /= PM_none ) then
                  call construct_overlap_scalapack(overlap_matrix)
                else
                  call setup_overlap_scalapack(overlap_matrix)
                endif
              endif

              if( use_constraint .or. flag_run_mulliken .or. &
                  out_lowdin .or. out_band .or. out_band_mulliken .or. &
              &   out_matrices .or. flag_KS==-2 .or. out_dos .or. &
                  use_wf_extrapolation .or. use_out_eigenvec_ovlp .or. &
                  transport_lead_calculation .or. transport_calculation .or. &
                  flag_run_tddft_real_time_propagation .or. fo_finalise.or. &
                  use_DFPT_dielectric) then
                write(info_str,'(2X,A)') "Not deallocating overlap matrix."
                call localorb_info(info_str,use_unit,'(A)', OL_norm)
              else
                write(info_str,'(2X,A)') "Deallocating overlap matrix."
                call localorb_info(info_str,use_unit,'(A)', OL_norm)
                call aims_deallocate( overlap_matrix,  "overlap_matrix" )
                call aims_allocate( overlap_matrix, 1, "overlap_matrix" )
              endif
            endif


          end if
        end if

        call get_times(time_intg, clock_time_intg, &
        &              tot_time_intg, tot_clock_time_intg)

!------------------- Calculate Fock Matrix (When Applicable) ----------------->>

        ! Get Hartree Fock hamiltonian.
        ! Fock term is added/updated here.

        ! Don't let the variable names in the following conditionals fool you:
        ! use_periodic_hf will also be .true. for cluster calculations using
        ! the LVL_fast resolution-of-identities method, and periodic
        ! calculations usuallly have use_hf_kspace.eq..false. by default.
        ! evaluate_exchange_matr_realspace_p0() is the subroutine that is
        ! mostly likely to be called for both periodic and non-periodic
        ! calculations.

        if (use_hartree_fock) then !.and. flag_xc_pre < 0) then
          call get_timestamps(time_fock, clock_time_fock)

          ! FK: In general, the EXX forces are switched on when the force
          ! calculation is on.
          EXX_forces_on = forces_on

          if (use_periodic_hf)then
             ! SVL: Add exact exchange energy here in case of periodic HF method
             if (number_of_loops.gt.0)then
               if (use_hf_kspace) then
                 call evaluate_exchange_matr_kspace_p0 &
                      (KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
                       occ_numbers_save)
               else
                ! FK: If we calculate the analytical stress, we do at least two
                !     iterations.  Therefore, we need to calculate the EXX
                !     forces only in the second iteration of the stress
                !     calculation.  That is the reason why we set
                !     'EXX_forces_on' to 'pulay_forces_on'.
                if (use_analytical_stress) EXX_forces_on = pulay_forces_on

                ! FK: We switch on the calculation of the EXX stress only in the
                !     second iteration of the stress calculation since we do not
                !     need its contribution earlier.  That is the reason why the
                !     argument is 'AS_pulay_stress_on' and not 'AS_stress_on'.
                if (.not.(keep_restart_info.and.restart_file_exists.and.number_of_loops.eq.1))& !SVL: restart
                  call evaluate_exchange_matr_realspace_p0 &
                       (KS_eigenvector,KS_eigenvector_complex,occ_numbers_save, &
                        exx_ene,d_exx_ene,EXX_forces_on,AS_pulay_stress_on)
               endif
               call exchange_matrix_mixing_p0(number_of_loops)

             endif

             !call get_exchange_energy_p0 &
             !     (hybrid_coeff, KS_eigenvector, KS_eigenvector_complex, &
             !      KS_eigenvalue, k_weights, occ_numbers_save, &
             !      chemical_potential,fock_energy,en_xc)

             call get_timestamps(rtime, clock_rtime)
             time_fock = rtime - time_fock
             tot_time_fock = tot_time_fock + time_fock
             call sync_timing(time_fock)
             clock_time_fock = clock_rtime - clock_time_fock
             tot_clock_time_fock = tot_clock_time_fock + clock_time_fock
          else ! .not. use_periodic_hf
            if (use_local_index) then
              ! RJ: Have still to check what exactly doesn't work here
              ! WPH (21 Feb 2018):  Apparantly, all RI methods other than LVL
              !     support non-packed matrices only, but use_local_index
              !     requires the CSR format (aka PM_index) be initialized.
              call aims_stop('use_local_index is only supported for EXX &
                             &calculations when using the LVL resolution of &
                             &identity method', func)
            endif

            en_hf=en_xc
            if (output_level == 'normal' .and. number_of_loops <= 1) &
               output_priority = 0
            if(use_screx) then
              call get_screx_hamiltonian &
              !call get_mp2_hamiltonian &
                   (number_of_loops,KS_eigenvalue,KS_eigenvector,n_electrons,&
              !      overlap_matrix,occ_numbers,hamiltonian, en_xc)
                    occ_numbers,hamiltonian, en_xc )
            else

              call get_hf_hamiltonian &
                    (number_of_loops, KS_eigenvalue, KS_eigenvector, &
                     n_electrons, occ_numbers, hamiltonian, fock_energy, en_xc)

            endif
            en_hf=en_xc-en_hf
            if (output_level == 'normal' .and. number_of_loops <= 1) &
               output_priority = 1

            call get_times(time_fock, clock_time_fock, &
                           tot_time_fock, tot_clock_time_fock)
          endif ! use_periodic_hf
        endif ! use_hartree_fock

!--------------- Solve Generalized Kohn-Sham Eigenvalue Problem -------------->>

        ! Here we construct the Hamiltonian that will enter into the
        ! Kohn-Sham solver from the real-space Hamiltonian evaluated previously,
        ! then we invoke the Kohn-Sham solver to solve the generalized
        ! Kohn-Sham eigenvalue problem.  The subroutine that performs the
        ! heavy lifting in this section is advance_KS_solution().

        if (use_embedding_pp) then
            call get_timestamps(time_embedding, clock_time_embedding)
            call add_nonlocal_pot(hamiltonian)
            call get_times(time_embedding, clock_time_embedding, &
            &              tot_time_embedding, tot_clock_time_embedding)
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

        ! in case of fo-dft no change in electron density is allowed.
        ! -> Set to converged
        if (fo_finalise) then
          elec_converged = .true.
          converged = .true.
        end if

      if (.not.fo_finalise) then
        call get_timestamps( time_diag, clock_time_diag)

        ! If ScaLAPACK is used, the Hamiltonian has to be converted from the
        ! real space format into a BLACS format.
        ! When LAPACK is used, the Hamiltonian is converted from the real space
        ! format into an upper-triangular form; this is done as part of
        ! advance_KS_solution() and is omitted here.
        if (use_scalapack) then
           if (.not.use_local_index) then
              if(n_periodic>0 .or. packed_matrix_format /= PM_none ) then
                 call construct_hamiltonian_scalapack(hamiltonian)
              else
                 call setup_hamiltonian_scalapack(hamiltonian)
              endif
           else ! use_local_index
              if(use_load_balancing) then
                 if(.not.first_iteration) then
                    ! Set the index arrays for local_index communication to the
                    ! integration batch permutation
                    call init_comm_full_local_matrix &
                         ( batch_perm(n_bp_integ)%n_basis_local, &
                         batch_perm(n_bp_integ)%i_basis_local )
                    call set_full_local_ham(my_hamiltonian_matrix)
                 else
                    ! We're already using the comms from the initial
                    ! integration batch permutation, so no need to reset them
                    call set_full_local_ham(my_hamiltonian_matrix)
                    ! Get new, hopefully better distribution
                    call reset_batch_permutation(n_bp_integ)
                    call compute_balanced_batch_distribution(n_bp_integ)
                 endif
              else
                 call set_sparse_local_ham_scalapack(my_hamiltonian_matrix)
              endif
              get_batch_weights = .false.
              use_batch_permutation = 0

              call aims_deallocate(my_hamiltonian_matrix, &
                                   "my_hamiltonian_matrix")
           end if
        end if
        ! wyj: TODO
        !print *, 'print_ham_SCF'
        !call print_ham_cpscf(ham)
        if (use_local_index .and. use_load_balancing) then
            if (n_periodic .ne. 0) call aims_stop("can't save ham when n_periodic>0", func)
            call save_ham_scalapack()
        endif

        ! Now set first_iteration to false
        first_iteration = .false.

        if(flag_rel.eq.REL_x2c)then
           call x2c_scf()
        elseif(flag_rel.eq.REL_4c_dks)then
           call q4c_scf()
        endif

        ! Solve the generalized Kohn-Sham eigenvalue problem
        call advance_KS_solution(overlap_matrix, hamiltonian, n_electrons, &
                KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
                occ_numbers, chemical_potential, chemical_potential_spin)

        call reshape_matrices_n_states(n_states)

        call evaluate_KH_core_states(overlap_matrix, rho, rho_gradient, &
                kinetic_density, partition_tab, hartree_potential)

        call add_KH_states(KS_eigenvalue, KS_eigenvector, &
                           KS_eigenvector_complex)

        if(use_periodic_hf) then !.and. flag_xc_pre < 0)then
           ! Make sure occ_numbers_save is always defined
           if(restart_zero_iteration) then
              occ_numbers_save = occ_numbers
           end if

           call get_exchange_energy_p0(hybrid_coeff, KS_eigenvector, &
                   KS_eigenvector_complex, KS_eigenvalue, k_weights, &
                   occ_numbers_save, chemical_potential, fock_energy, en_xc)
        endif

        call get_times(time_diag, clock_time_diag, &
        &              tot_time_diag, tot_clock_time_diag)

        if (output_level == 'normal' .and. number_of_loops <= 1) &
           output_priority = 0

        if (.not. use_elsi_dm) then
           call output_real_eigenfunctions &
                ( KS_eigenvalue, KS_eigenvector, occ_numbers  )
           ! if we're adjusting any s.c.f. stability related settings,
           ! this must be done right after calling output_real_eigenfunctions,
           ! in which the type of system (low gap or not) is determined
           if (adjust_scf.and.(adjust_scf_iteration.eq.number_of_loops)) then
              call adjust_scf_settings
              if (.not.adjust_scf_always) then
                 ! adjust only once
                 adjust_scf = .false.
              else
                 ! adjust for each new s.c.f. cycle.
                 ! There is a corner case, which is that for all following
                 ! s.c.f. cycles, we do not solve the eigenvalue problem during
                 ! the initialization.
                 ! This means that we can only check from the first iteration
                 ! onwards, not in the zeroth iteration:
                 if (adjust_scf_iteration.eq.0) then
                    adjust_scf_iteration = 1
                 end if
              end if
           end if
        end if

        if(n_periodic > 0  )then
           if (use_dipole_correction .or. calculate_work_function) then
              call output_energy_levels &
                   (KS_eigenvalue, occ_numbers, chemical_potential, &
                    chemical_potential_spin, converged )
           end if
        end if


        if (output_level == 'normal' .and. number_of_loops <= 1) &
           output_priority = 1

        ! determine sum of eigenvalues

        ! if (use_constraint) then
        !   call get_energy_correction &
        !        ( KS_eigenvector, occ_numbers, overlap_matrix )
        ! endif

        if (.not. use_elsi_dm) then
           call get_ev_sum_p0 &
                (occ_numbers, KS_eigenvalue, av_core_shift, ev_sum, &
                 ev_sum_shifted)
        end if

        if (use_embedding_pp) then
           call get_timestamps(time_embedding, clock_time_embedding)
           call evaluate_nonlocal_energy(KS_eigenvector, &
                                         occ_numbers)
           if(use_nonlinear_core) then
              call evaluate_nlcc_correction_energy(partition_tab, rho, &
                                                   rho_gradient)
           endif
           call get_times(time_embedding, clock_time_embedding, &
           &              tot_time_embedding, tot_clock_time_embedding)
        endif


        if ( &
            use_vdw_correction_hirshfeld_sc &
            .or. (use_mbd_std .and. mbd_self_consistent) &
        ) then
           write (info_str,'(2X,A)') "Evaluating Hirshfeld Density Update"
           call localorb_info(info_str, use_unit, '(A)', OL_norm)
           call hirsh()
        endif
      else
            write(info_str,'(2X,A)') ''
            call localorb_info(info_str,use_unit,'(A)', OL_norm )
            write(info_str,'(2X,A)') 'Final Fragment-Orbital DFT step, skipped the solving of '&
                                  // 'the eigenvalue problem!'
            call localorb_info(info_str,use_unit,'(A)', OL_norm )

            write(info_str,'(2X,A)') 'All following post-processing steps are done with the initial '&
                                  // 'density read from the fragment restart files.'
            call localorb_info(info_str,use_unit,'(A)', OL_norm )
            write(info_str,'(2X,A)') 'Please note: This can lead to unexpected behaviour in '&
                                  // 'different quantities'
            call localorb_info(info_str,use_unit,'(A)', OL_norm )
            write(info_str,'(2X,A)') 'due to the mixed occupations. Especially cube files'
            call localorb_info(info_str,use_unit,'(A)', OL_norm )
            write(info_str,'(2X,A)') 'are known to be affected.'
            call localorb_info(info_str,use_unit,'(A)', OL_norm )
      end if !fo_finalise

!---------------------------- Check SCF Convergence -------------------------->>

        ! Check convergence to determine whether we need to perform another
        ! iteration.

        ! First, check convergence of total energy, using only total energy
        ! pieces pertaining to the Kohn-Sham density of the present s.-c. loop

        call get_timestamps(time_etot, clock_time_etot)

        ! If finite-width smearing was used, must get entropy correction for T=0
        ! total energy.
        ! This routine supports k-points.

        if ((.not. use_elsi_dm) .and. (.not. fo_finalise)) then
           call get_entropy_correction_p1(KS_eigenvalue, occ_numbers, &
                chemical_potential, chemical_potential_spin, entropy_correction)

           call get_penalty_energy(KS_eigenvector, KS_eigenvector_complex, &
                occ_numbers, condition_penalty, penalty_energy)
        else
           entropy_correction = 0.0d0
           penalty_energy = 0.0d0
        end if

        call get_total_energy &
        ( ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
          en_ion_embed, en_density_embed, &
          en_vdw, en_pot_vdw, en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, &
          en_elec_free, en_elec_delta, fock_energy, &
          hartree_energy_free, hartree_delta_energy, &
          hartree_multipole_correction, &
          entropy_correction, penalty_energy, total_energy &
        )


        ! If we have reached the point in the SCF cycle where forces are being
        ! calculated, add up all contributions to the forces to get the total
        ! force, then store the forces to check for convergence
        if (forces_on) then
          call sync_forces(pulay_forces, hellman_feynman_forces, &
                           multipole_forces, gga_forces, gga_forces_on, &
                           nlcc_forces, nlcc_forces_on,Gnonmf_forces, &
                           Gnonmf_forces_on)

          call get_total_forces &
              (ionic_forces, pulay_forces, hellman_feynman_forces, &
               multipole_forces, gga_forces, gga_forces_on, &
               total_forces, vdw_forces, d_exx_ene, EXX_forces_on, &
               pseudocore_forces, nlcc_forces, nlcc_forces_on,Gnonmf_forces, &
               Gnonmf_forces_on)

          call compare_forces &
          ( pulay_forces, hellman_feynman_forces, &
            multipole_forces, previous_forces, &
            diff_forces )

          if (AS_stress_on) then
            ! We count how many calculations of the analytical stress are done.
            ! This is used to switch on the convergence check not until
            ! AS_counter>=2.
            ! See further down for more detailed explanation.
            AS_counter = AS_counter+1

            ! Syncing first
            if (use_mpi) then
              call AS_sync_analytical_stress(AS_pulay_stress_on)
            end if

            ! CC: Gets and prints the total analytical stress:
            ! FK: We need 2 iterations for the first meaningful result.
            ! FK: If vdW is used, the stress tensor w/o the vdW correction will
            !     be printed here.
            if (AS_counter .ge. 2) then
              call get_total_analytical_stress( &
                   AS_pulay_stress_on, (en_xc-en_pot_xc), fock_energy, &
                   previous_stress, diff_stress, &
                   use_vdw_correction_hirshfeld .or. use_mbd_dev .or. use_mbd_std .or. use_libmbd)
            end if
          end if
        end if

        ! For XAS calculations needed: set projection_state to actual forced
        ! state
        if (force_occupation_projector.and.start_force_occ.and.n_periodic.eq.0) then
          force_occ_pr_state2(:) = force_occ_pr_state(:)
          force_occ_pr_state(:) = force_occ_state(:)
        else if (force_occupation_projector.and.start_force_occ) then
          force_occ_pr_state_periodic(:,:) = force_occ_state_periodic(:,:)
        end if

        call get_times(time_etot, clock_time_etot, &
        &              tot_time_etot, tot_clock_time_etot)

        ev_sum_change = ev_sum - previous_ev_sum

        ! check convergence of self-consistency loop
        ! give self-consistency convergence status report

        write(info_str,'(A)') ''
        call localorb_info(info_str,use_unit,'(A)', OL_norm  )

        write(info_str,'(2X,A)') &
        "Self-consistency convergence accuracy:"
        call localorb_info(info_str,use_unit,'(A)', OL_norm  )

        if (n_spin.eq.2) then
          write(info_str,'(2X,A,1X,E11.4,1X,E11.4)') &
                "| Change of charge/spin density :", &
                (rho_change(i_spin),i_spin=1,n_spin,1)
        else
          write(info_str,'(2X,A,1X,E11.4,1X,E11.4)') &
                "| Change of charge density      :", &
                (rho_change(i_spin),i_spin=1,n_spin,1)
        end if
        call localorb_info(info_str,use_unit,'(A)', OL_norm  )

        write(info_str,'(2X,A,1X,E11.4,A)') &
             "| Change of sum of eigenvalues  :", &
             (ev_sum - previous_ev_sum)*hartree, " eV"
        call localorb_info(info_str,use_unit,'(A)', OL_norm  )

        write(info_str,'(2X,A,1X,E11.4,A)') &
             "| Change of total energy        :", &
             (total_energy - previous_total_energy) * hartree, " eV"
        call localorb_info(info_str,use_unit,'(A)', OL_norm  )

        if (forces_on) then
          write(info_str,'(2X,A,1X,E11.4,A)') &
          "| Change of forces              :", diff_forces, " eV/A"
          call localorb_info(info_str,use_unit,'(A)', OL_norm  )
        end if

        if ( AS_stress_on) then
          ! Printing the change of analytical stress makes only sense for
          ! AS_counter>=3, because we need 2 iterations for the first result.
          if (AS_counter .ge. 3) then
            write(info_str,'(2X,A,1X,E11.4,A)') &
            "| Change of analytical stress   :", diff_stress, " eV/A**3"
            call localorb_info(info_str,use_unit,'(A)', OL_norm  )
          end if
        end if
        write(info_str,'(A)') ''
        call localorb_info(info_str,use_unit,'(A)', OL_norm  )

        ! summarized convergence for really brief output if desired
        if (output_level .eq. 'MD_light') then
          if (n_spin.eq.2) then
              write(info_str,'(A,E9.2,1X,E9.2)') &
                  " | ", (rho_change(i_spin),i_spin=1,n_spin,1)
          else
              write(info_str,'(A,E9.2)') &
                  " | ", (rho_change(i_spin),i_spin=1,n_spin,1)
          end if
          call localorb_info(info_str,use_unit,'(A,$)',OL_high)
          write(info_str,'(A,E10.2)') &
                " | ", (ev_sum - previous_ev_sum)*hartree
          call localorb_info(info_str,use_unit,'(A,$)',OL_high)
          write(info_str,'(A,E9.2)') &
                " | ", (total_energy - previous_total_energy) * hartree
          call localorb_info(info_str,use_unit,'(A,$)',OL_high)
          if (forces_on) then
              write(info_str,'(A,E9.2)') &
                  " |     ", diff_forces
          else
              write(info_str,'(A)') &
                  " |             ."
          end if
          call localorb_info(info_str,use_unit,'(A,$)',OL_high)
        end if

        ! check if preconditioner is on and if it is still necessary
        if (use_kerker_preconditioner.and.kerker_preconditioner_on) then
          precondition_check = .true.
          ! check charge convergence
          do i_spin = 1, n_spin
              precondition_check = precondition_check.and.(rho_change(i_spin).lt.preconditioner_turnoff_charge)
          end do
          ! check convergence wrt energy and sum eigenvalues
          precondition_check = precondition_check &
                .and.(abs(total_energy - previous_total_energy)*hartree.lt.preconditioner_turnoff_energy) &
                .and.(abs((ev_sum - previous_ev_sum)*hartree).lt.preconditioner_turnoff_sum_ev)
          if (precondition_check) then
              kerker_preconditioner_on = .false.
              ! substitute back the non-Kerker Pulay mixing parameter
              charge_mix_param(1) = stored_charge_mix_param
              write(info_str,'(2X,A)') 'Preliminary charge convergence reached. Turning off preconditioner.'
              call localorb_info(info_str,use_unit,'(A)', OL_norm )
              call localorb_info('',use_unit,'(A)', OL_norm )
          end if
        end if


        ! Check electronic self-consistency convergence criteria
        ! Only check criteria which were set explicitly.
        ! All convergence criteria must be fulfilled.

        if (.not.elec_converged) then
          if (force_potential /= 0 .or. &
          &   (sc_accuracy_rho.gt.0.d0) .or. flag_acc_eev .or. flag_acc_etot .or. flag_acc_potjump .or. &
          &   (use_forces .and. sc_accuracy_forces > 0.d0)) then
              ! There is some convergence criterion.
              elec_converged = .true.
              ! If some electronic criterion is not met yet, elec_converged
              ! will be reset.  Otherwise (if there is none or all are met),
              ! only force_converged is still missing.

              frozen_core_converged = .true.
          end if
          ! If no convergence criterion is given (not even for the forces),
          ! iterate forever (or until number_of_loops > sc_iter).

          ! sc_accuracy_rho is always set to some default value.
          ! Only if it is deliberately set to zero or below should it be not checked.
          if (sc_accuracy_rho.gt.0.d0) then
              do i_spin = 1, n_spin, 1
                elec_converged = elec_converged .and. &
                (rho_change(i_spin) .lt. sc_accuracy_rho)

                frozen_core_converged = frozen_core_converged .and. &
                (rho_change(i_spin) .lt. frozen_core_scf_factor*sc_accuracy_rho)
              enddo
          end if

          if (flag_acc_eev) then
              elec_converged = elec_converged .and. &
              (abs(ev_sum - previous_ev_sum)*hartree .lt. sc_accuracy_eev)
          end if

          if (flag_acc_etot) then
              elec_converged = elec_converged .and. &
              (abs(total_energy - previous_total_energy) * hartree .lt. sc_accuracy_etot)
          end if

          if (flag_acc_potjump) then
             elec_converged = elec_converged .and. &
             (abs(pot_jump - previous_pot_jump) * hartree * 2.d0 .lt. sc_accuracy_potjump)
          endif
        end if

        ! If the electronic self-consistency criteria are not yet fulfilled, we
        ! have the option (configurable in control.in using keyword sc_init_iter)
        ! to stop the scf cycle anyway and restart the mixer from scratch based
        ! only on the current, not yet self-consistent orbitals. This will happen
        ! only once and is intended to look exactly like a geometry step - except
        ! that the geometry remains unchanged.
        if (.not.scf_restarted) then
           if (.not.elec_converged) then
              ! check to make sure that we only ever restart the scf
              ! cycle from scratch if we are not already practically convergedd.
                 if ( (number_of_loops.ge.sc_init_iter) .and. (number_of_loops.lt.sc_iter_limit) )then
                    if (flag_acc_rho) then
                       density_is_close_to_scf = .true.
                       do i_spin = 1, n_spin, 1
                         density_is_close_to_scf = density_is_close_to_scf .and. &
                                          (dabs(rho_change(i_spin)) .le. (sc_init_factor*sc_accuracy_rho) )
                       end do
                       restarting_scf = .not.(density_is_close_to_scf) ! true unless the density
                                                                       ! is close to scf in all
                                                                       ! spin channels
                    else
                       restarting_scf = .true.
                    end if
                 end if
           end if
        end if

        ! End of frozen core approximation
        if (.not. elec_converged .and. frozen_core_scf .and. frozen_core_converged) then
           frozen_core_scf = .false.
           restarting_scf = .true.

           if (use_scalapack) then
              call unpermute_ovlp_scalapack()
           end if

           call aims_elsi_set_output(eh_scf,0)
           call aims_elsi_finalize_scf()
           call aims_elsi_init_scf(my_scalapack_comm_work, my_blacs_ctxt, mb, &
                mxld, mxcol,my_k_point)
        end if

        ! If forces needed and electronic self-consistency cycle is converged,
        ! check force convergence subsequently.

        if (use_forces) then
          ! switch GGA forces on only in the second iteration with forces
          if (use_gga.and.forces_on) then
            gga_forces_on = .true.
            if (use_meta_gga) then
              meta_gga_forces_on = .true.
            end if
          end if
          if (use_Gnonmf_forces.and.forces_on) then
            !if we use a self-consistent Gnonmf, we get additional terms for the
            !forces, to the pulay forces and extra terms, which are calculated
            !in the next SCF step
            Gnonmf_forces_on = .true.
          end if

          !VB: The nlcc functionality is currently not supported but
          !    the check done here is misplaced. This should only
          !    happen if forces_on is true ....
          !    Please fix if this is ever picked up again (not
          !    fixing now to avoid misunderstanding on my part)
          if (use_nonlinear_core) then
             nlcc_forces_on = .true.
          end if

          ! Density matrix case: switch Pulay+GGA forces on only in the second
          ! iteration with forces
          if (forces_on) then
            pulay_forces_on = .true.
          end if
          ! ... same for stress
          if (AS_stress_on.and.AS_counter.ge.1) then
          !if (AS_stress_on) then
            AS_pulay_stress_on = .true.
          end if

          if (elec_converged.and.mpb_converged) then

            ! If forces_on was not true, this is the first converged iteration -
            ! set it to true anyway.
            if (.not.forces_on) then

              write(info_str,'(2X,A,A)') &
                  "Electronic self-consistency reached - switching on the &
                  &force computation."
              call localorb_info(info_str,use_unit,'(A)',OL_norm)
              call localorb_info('',use_unit,'(A)',OL_norm)

              forces_on = .true.

              diff_forces = 2*sc_accuracy_forces+1.d0

              ! only if we're already in the last allowed scf cycle, switch
              ! gga forces on already here!
              ! But not in the analytical stress case since we need at least
              ! 2 iterations
              if ( (number_of_loops.ge.(sc_iter_limit-1)).and.(.not.use_analytical_stress) ) then
                pulay_forces_on = .true.
                if (use_gga) then
                  gga_forces_on = .true.
                  if (use_meta_gga) then
                    meta_gga_forces_on = .true.
                  end if
                end if
                if (use_Gnonmf_forces) then
                    Gnonmf_forces_on = .true.
                end if
                if (use_nonlinear_core) then
                  nlcc_forces_on = .true.
                end if
              end if

              ! no force convergence check? Turn on GGA forces!
              ! But not in the analytical stress case since we need at least
              ! 2 iterations
              if ( ((sc_accuracy_forces.le.0d0).and.(.not.use_analytical_stress)) .or. &
                   ((sc_accuracy_forces.le.0d0).and.(AS_counter.ge.1           )) )    then
                pulay_forces_on   = .true.
                if (use_gga) then
                  gga_forces_on   = .true.
                  if (use_meta_gga) then
                    meta_gga_forces_on = .true.
                  end if
                end if
                if (use_Gnonmf_forces) then
                    Gnonmf_forces_on = .true.
                end if
                if (use_nonlinear_core) then
                  nlcc_forces_on = .true.
                end if
              end if

            else ! this is the case if forces_on
              if (sc_accuracy_forces.gt.0.d0) then
                force_converged = (diff_forces .lt. sc_accuracy_forces)
              else
                ! force convergence will not be checked with explicitly
                ! zero or negative setting of sc_accuracy_forces
                if ((.not.use_analytical_stress).or.(AS_counter.ge.2)) then
                  force_converged= .true.
                else
                  force_converged= .false.
                end if

              end if
            end if ! forces_on

            ! CC: Switch on stress:
            if ( use_analytical_stress ) then
              if ( .not. AS_stress_on ) then
                write(info_str,'(2X,A)') &
                    "Electronic self-consistency reached - switching on the &
                    &analytical stress tensor computation."
                call localorb_info(info_str,use_unit,'(A)',OL_norm)
                call localorb_info('',use_unit,'(A)',OL_norm)
                AS_stress_on = .true.
                !pulay_forces_on = .true.
                diff_stress = 2*sc_accuracy_stress+1.d0
                !AS_counter  = 1
              else
                ! We need at least two iterations for the stress.
                ! In the first interation we calculate the electrostatic
                ! contributions in sum_up_whole_potential.
                ! We then use d/de v from the first iteration in the second
                ! intertion in Pulay_force_densmat.
                ! Therefore, we check for convergence not until the second
                ! interation.
                if (AS_counter .ge. 2) then
                  if (sc_accuracy_stress .gt. 0.0d0) then
                    stress_converged = ( diff_stress .lt. sc_accuracy_stress )
                  else
                    stress_converged = .true.
                  end if
                end if
              end if
            end if
          else if ( (number_of_loops.ge.(sc_iter_limit-1)) .and. &
                    (.not.use_analytical_stress) .and. &
                    (.not.restarting_scf) ) then
            if (.not.forces_on) then
              write(info_str,'(2X,A,A)') &
              "Reaching maximum number of scf iterations ", &
              "(no convergence)."
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
              "Switching on the force computation anyway."
              call localorb_info(info_str,use_unit,'(A)')
              call localorb_info('',use_unit,'(A)',OL_norm)
            end if
            forces_on = .true.
            pulay_forces_on = .true.
            if (use_gga) then
              gga_forces_on = .true.
              if (use_meta_gga) then
                meta_gga_forces_on = .true.
              end if
            end if
            if (use_nonlinear_core) then
              nlcc_forces_on = .true.
            end if
            if (use_Gnonmf_forces) then
                Gnonmf_forces_on = .true.
            end if
          else if ( (number_of_loops.ge.(sc_iter_limit-2)) .and. &
                    (use_analytical_stress) .and. &
                    (.not.restarting_scf) ) then
            if (.not.forces_on) then
              write(info_str,'(2X,A,A)') &
              "Reaching maximum number of scf iterations ", &
              "(no convergence)."
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
              "Switching on the force computation anyway."
              call localorb_info(info_str,use_unit,'(A)')
              call localorb_info('',use_unit,'(A)',OL_norm)
            end if
            forces_on = .true.
            pulay_forces_on = .true.
            if (use_gga) then
              gga_forces_on = .true.
              if (use_meta_gga) then
                meta_gga_forces_on = .true.
              end if
            end if
            if (use_nonlinear_core) then
              nlcc_forces_on = .true.
            end if
            if (use_Gnonmf_forces) then
                Gnonmf_forces_on = .true.
            end if
            write(info_str,'(2X,A)') &
                "Switching on analytical stress tensor computation."
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            call localorb_info('',use_unit,'(A)',OL_norm)
            AS_stress_on = .true.
            !pulay_forces_on = .true.
            diff_stress = 2*sc_accuracy_stress+1.d0
            !AS_counter  = 1
          end if

          ! compute pseudopotential force contribution
          if(use_embedding_pp.and.force_converged) then
            call localorb_info('',use_unit,'(A)',OL_norm)
            write(info_str,'(2X,A)') &
                 "All-electron forces are converged. Now adding pseudopotential &
                 &contribution: "
            call localorb_info(info_str,use_unit,'(A)',OL_norm)

            call evaluate_pp_nonlocal_forces(KS_eigenvector, occ_numbers)

            if(nlcc_forces_on) then
              call evaluate_pp_nlcc_forces(rho, partition_tab, rho_gradient, kinetic_density)
            endif

            call evaluate_pp_hellman_feynman_forces_p0(rho, partition_tab, pseudocore_forces)

            call get_total_forces &
                 (ionic_forces, pulay_forces, hellman_feynman_forces, &
                 multipole_forces, gga_forces, gga_forces_on, &
                 total_forces, vdw_forces, d_exx_ene, EXX_forces_on, pseudocore_forces, &
                 nlcc_forces, nlcc_forces_on, Gnonmf_forces, Gnonmf_forces_on)
          endif
        end if ! use_forces

        converged = (elec_converged.and.force_converged.and.stress_converged.and.mpb_converged)

        !SR: probably we can remove this postprocessing, since a new flag
        !restarting_scf emerged, have to check this
        if (solvent_method.eq.SOLVENT_MPB) then
          ! Basically, only update the cavity and surface volume
          call postprocess_mpb_solver(elec_converged,converged)
        end if

        if ((number_of_loops > n_steps_xc_pre .or. converged) .and. flag_xc_pre > 0) then
            flag_xc_pre = -2
            restarting_scf = .true.
            hybrid_coeff = hybrid_coeff_old
            flag_xc = flag_old_xc
            ! AJL: Necessary as we can have situations where we don't want to
            ! use Hybrid-DFT after pre-convergence e.g. meta-GGAs
            if (hybrid_coeff .gt. 0.0) then
              use_hartree_fock = .true.
              use_periodic_hf = .true.
            end if
        end if
        if (restarting_scf) then
           ! Fake convergence. However, we must make sure below we do not
           ! write confusing output.
           converged = .true.
        end if

        ! TODO: mpe currently does not actively try to converge.
        ! Therefore, if(converged), more SCF cycles will not improve the mpe
        ! solution. Here we simply check whether or not it is solved.
        if ((solvent_method.eq.SOLVENT_MPE) .and. converged .and. (.not. mpe_solved)) then
           write(info_str, '(A)') ""
           call localorb_info(info_str)
           write(info_str, '(A)') "*******************************************************************"
           call localorb_info(info_str)
           write(info_str, '(1X,A)') '* Could not solve MPE equations in given basis.'
           call localorb_info(info_str)
           write(info_str, '(1X,A)') '* This is even after SCF convergence. Try increasing'
           call localorb_info(info_str)
           write(info_str, '(A)') ""
           call localorb_info(info_str)
           write(info_str, '(1X,A)') "          mpe_lmax_rf"
           call localorb_info(info_str)
           write(info_str, '(A)') ""
           call localorb_info(info_str)
           write(info_str, '(1X,A)') "* See user's guide for further info. Aims will abort now."
           call localorb_info(info_str)
           write(info_str, '(A)') "*******************************************************************"
           call localorb_info(info_str)
           write(info_str, '(A)') ""
           call localorb_info(info_str)
           call aims_stop_coll("Error (MPE): Adjusted R^2 below 1-tolerance. "//&
             "Could not solve MPE equations even after SCF convergence.")
        endif

!---------------------- Output of Final SCF Quantities ----------------------->>
        ! plus_u matrix control: writes the self consistent dft+u occupation matrix to a file
        ! which then can be edited manually
        if(use_plus_u.and.converged.and.plus_u_occupation_matrix_control_write) then
          call plus_u_write_occ_mat()
        endif

        if(use_plus_u.and.plus_u_matrix_release_defined) then
        ! check if occupation matrix should be calculated self-consistently again
          if(abs((total_energy - previous_total_energy) * hartree) .lt. plus_u_matrix_release) then
            ! we switch of matrix control read
            ! if the release criteria is fullfilled everything will be self-
            ! consistently from that point on
            ! However, the occupation matrix is still written out to a file
            ! (for later use)
            plus_u_occupation_matrix_control_read = .false.
            plus_u_matrix_release_defined = .false.
            write(info_str,'(2X,A,A)') &
               "Switching off DFT+U matrix control ", &
               "DFT+U occupation matrix is now calculated self-consistently"
               call localorb_info(info_str,use_unit,'(A)')
          endif
        ! the occupation matrix is now calculated self-consistently. If converged, the occupation matrix
        ! is written to a file.
        endif

        if(use_plus_u.and.converged.and.plus_u_matrix_error_defined) then
           call plus_u_matrix_error()
        endif

        if(use_plus_u.and.converged.and.plus_u_eigenvalues_defined) then
           call plus_u_eigenvalues()
        endif

!  ---------- Update electron density and perform mixing ---->>

        if (converged) then
          ! We are done - no further evaluation of density / potential needed
          ! Prepare output of final s.c.f. quantities.

          if(output_level .eq. 'MD_light') then

            ! Must output the time for current s.c.f. loop already here
            ! to avoid breaking output format
            call get_times(time_sc_loop, clock_time_sc_loop)

            write(info_str,'(A,F12.3,A,F12.3,A)') &
                   " | ", time_sc_loop, " s | ", clock_time_sc_loop, " s"
            call localorb_info(info_str,use_unit,'(A)',OL_high)


            output_priority = 1
            call get_total_energy &
            ( ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
              en_ion_embed, en_density_embed, &
              en_vdw, en_pot_vdw, en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c,&
              en_elec_free, en_elec_delta, fock_energy, &
              hartree_energy_free, hartree_delta_energy, &
              hartree_multipole_correction, &
              entropy_correction, penalty_energy, total_energy &
            )

            if(n_periodic > 0  )then
               if(use_dipole_correction .or. calculate_work_function)then
                  call output_energy_levels &
                      (KS_eigenvalue, occ_numbers, chemical_potential, &
                       chemical_potential_spin, converged)
               end if
            end if

          end if ! output level

          ! (1) output eigenvalues if not already done so
          if (output_level.eq.'normal') output_priority = 0

          if (.not. use_elsi_dm) then
             call output_real_eigenfunctions &
                  ( KS_eigenvalue, KS_eigenvector, occ_numbers)
          end if

          if(n_periodic > 0  )then
             if(use_dipole_correction .or. calculate_work_function)then
                call output_energy_levels &
                     (KS_eigenvalue, occ_numbers, chemical_potential, &
                      chemical_potential_spin, converged)
             end if
          end if


          if (output_level.eq.'normal') output_priority = 1

          write(info_str,'(A)') ''
          call localorb_info(info_str,use_unit,'(A)',OL_norm)
          if (.not.restarting_scf) then
             ! normal convergence was reached
             write(info_str,'(2X,A)') "Self-consistency cycle converged."
          else
             ! In this case, convergence was not reached, simply restarting scf.
             write(info_str,'(2X,A)') &
                  "Self-consistency cycle not yet converged - restarting mixer &
                   &to attempt better convergence."
          end if
          call localorb_info(info_str,use_unit,'(A)',OL_norm)
          write(info_str,'(A)') ''
          call localorb_info(info_str,use_unit,'(A)',OL_norm)

          !Output of eigenvalues on reciprocal k_grid
          if (out_fermisurface.and.(.not.restarting_scf)) call output_bxsf

          ! Final output of density matrix
          if(elsi_write_dm) then
             if(.not. use_elsi_dm) then
                call kweight_occs("output_densmat",occ_numbers)
             endif

             call output_densmat(KS_eigenvector,KS_eigenvector_complex,occ_numbers)

             if(.not. use_elsi_dm) then
                call de_kweight_occs("output_densmat",occ_numbers)
             endif
          endif

          if (keep_restart_info.and.restart_write.and.use_hf_realspace) then !SVL: restart
             if (use_restart_save_iterations) then
                if (.not.(mod(number_of_loops,restart_save_iterations).eq.0)) then
                   if(myid.eq.0)&
                        write(use_unit,*) 'Writing Hartree-Fock exchange &
                                          &matrix to disk for restart'
                   open(file = exx_restart, unit = 7, status = 'unknown', &
                        form = 'unformatted')
                   if(real_eigenvectors)then
                      write(7) hf_exchange_matr_real
                   else
                      write(7) hf_exchange_matr_complex
                   endif
                   close(7)
                end if
             else
                if(myid.eq.0)&
                     write(use_unit,*) 'Writing Hartree-Fock exchange matrix &
                                       &to disk for restart'
                open(file = exx_restart, unit = 7, status = 'unknown', &
                     form = 'unformatted')
                if(real_eigenvectors)then
                   write(7) hf_exchange_matr_real
                else
                   write(7) hf_exchange_matr_complex
                endif
                close(7)
             end if
          end if

          if (restarting_scf) then
            if (use_hartree_fock) then
               if(use_hf_realspace .and. exx_band_structure_version .eq. 2)then
                  call cleanup_hartree_fock_p0()
               else
                  call cleanup_hartree_fock()
                  if (use_lvl) call cleanup_localized_basbas()
               endif
            endif
          else
            if (use_hartree_fock &
                 .and. (.not.use_corr) &
                 .and. (.not.calculate_atom_bsse) &
                 .and. (.not.use_DFPT_polarizability) &
                 .and. (.not.use_DFPT_polar_reduce_memory) &
                 .and. (.not.magnetic_response)) then
               if(use_hf_realspace .and. exx_band_structure_version .eq. 2)then
                  call cleanup_hartree_fock_p0()
               else
                  call cleanup_hartree_fock()
                  if (use_lvl) call cleanup_localized_basbas()
               endif
            endif
          endif

          ! This is for the case where you have to "unshift" the eigenvalues in
          ! the end of self consistency when spin polarized calculations with
          ! the multiplicity flag are used
          if(use_hf_multiplicity.and.n_periodic.eq.0) then
            call unshift_eigenvalues(n_electrons, &
            KS_eigenvalue,KS_eigenvector, occ_numbers,  &
            chemical_potential, chemical_potential_spin)
          endif
          if(output_level .eq. 'MD_light') output_priority = 2
        else if (number_of_loops.ge.sc_iter_limit) then
          ! We are not converged, but this was the last self-consistency cycle -
          ! we do not need any more potential / density evaluations
          below_it_limit = .false.
        end if ! converged

        ! check whether or not to save restart information on this iteration,
        ! which can happen only if ...
        if (use_restart_save_iterations) then
          if ((keep_restart_info.and.restart_write).and. &
                (.not.restart_zero_iteration).and.        &
                (mod(number_of_loops,restart_save_iterations).eq.0)) then
            call write_restart_info()
            if(use_hf_realspace)then !SVL: restart
               if(myid.eq.0)&
                    write(use_unit,*) 'Writing Hartree-Fock exchange matrix to &
                                      &disk for restart'
               open(file = exx_restart, unit = 7, status = 'unknown', &
                    form = 'unformatted')
               if(real_eigenvectors)then
                  write(7) hf_exchange_matr_real
               else
                  write(7) hf_exchange_matr_complex
               endif
               close(7)
            endif
          end if
        end if

!------------------- Print Out Timings For This Iteration -------------------->>

        restart_zero_iteration = .false.

        ! Computation for the current SCF iteration ends here.

        if (pulay_forces_on .and. lda_update_switch) then
          use_density_matrix = .false.
          write(info_str,'(2X,A)') &
             "Switched to orbital based density and force update"
          call localorb_info(info_str,use_unit,'(A)')
        end if

        ! if converged and output_level MD_light, we have already printed the
        ! timing.
        if (.not. (converged.and.(output_level.eq.'MD_light')) ) then
          call get_times(time_sc_loop, clock_time_sc_loop)
        end if

        if (.not.converged) then
          ! If converged, we must already measure the times higher up to produce
          ! readable output
          if (output_level.eq.'MD_light') then
            ! finalize output line for this s.c.f. iteration with the timing
            write(info_str,'(A,F12.3,A,F12.3,A)') &
               " | ", time_sc_loop, " s | ", clock_time_sc_loop, " s"
            call localorb_info(info_str,use_unit,'(A)',OL_high)
          end if
        end if

        write(info_str,'(A,I5)') &
        & "End self-consistency iteration # ", number_of_loops
        call output_timeheader(deffmt, info_str, OL_norm)
        call output_times(deffmt, "Time for this iteration", &
        &                 time_sc_loop, clock_time_sc_loop, OL_norm)
        if (.not.forces_on) then
          call output_times(deffmt, "Charge density update", &
          &                 time_density, clock_time_density, OL_norm)
        else
          call output_times(deffmt, "Charge density & force component update",&
          &                 time_density, clock_time_density, OL_norm)
        end if
        if (use_kerker_preconditioner.and.kerker_preconditioner_on) then
          call output_times(deffmt, "Density mixing & preconditioning", &
          &                 time_mixing, clock_time_mixing, OL_norm)
        else
          call output_times(deffmt, "Density mixing", &
          &                 time_mixing, clock_time_mixing, OL_norm)
        end if
        call output_times(deffmt, "Hartree multipole update", &
        &                 time_hartree_multi, clock_time_hartree_multi,OL_norm)
        call output_times(deffmt, "Hartree multipole summation", &
        &                 time_hartree_sum, clock_time_hartree_sum, OL_norm)
        call output_times(deffmt, "Integration", &
        &                 time_intg, clock_time_intg, OL_norm)
        if(use_hartree_fock) then
          call output_times(deffmt, "Fock matrix evaluation", &
          &                 time_fock, clock_time_fock, OL_norm)
        endif
        call output_times(deffmt, "Solution of K.-S. eqns.", &
        &                 time_diag, clock_time_diag, OL_norm)
        call output_times(deffmt, "Total energy evaluation", &
        &                 time_etot, clock_time_etot, OL_norm)

        write(info_str,*)
        call localorb_info(info_str,use_unit,'(A)',OL_norm)
        call aims_mem_current_output ()
        write(info_str,'(A)') &
        "------------------------------------------------------------"
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

        ! When requested, output information about the completed SCF iteration
        ! to the JSON log file
        if (out_aims_json_log) call synchronize_batch_statistics()
        if (out_aims_json_log) call write_scf_iteration_to_json_log(AS_counter)

        ! Check the walltime limitations and keep track of maximal scf cycle
        ! time
        call update_time_max_scf_cycle(time_sc_loop)
        enough_walltime_left = time_check_scf_cycle()

        ! check if there is control.update.in
        inquire (file="control.update.in", exist=control_update_requested)

        !check if the user wants to stop the SCF
        INQUIRE(FILE="abort_scf", EXIST=abort_requested) ! abort_requested will be TRUE if the file
                                                         ! exists and FALSE otherwise
        INQUIRE(FILE="abort_and_postprocess", EXIST=abort_and_postprocess)   ! allow the user to ask for postprocessing anyway
        if (abort_and_postprocess) then
            abort_requested=.true.
            postprocess_anyway=PP_ANYWAY_EVERYTHING
        endif

        ! Check if we need to abandon the s.c.f. cycle anyway, due to excessive
        ! total energy changes !between iterations (if the mixer has mixed
        ! itself into a disaster, it makes no sense to waste X more iterations
        ! in a queue).
        if (check_sc_abandon_etot) then
           if ( ( abs(total_energy - previous_total_energy)*hartree/dble(n_atoms) ) .lt. sc_abandon_etot_threshold ) then
             sc_abandon_counter = 0
           else
             sc_abandon_counter = sc_abandon_counter + 1
           end if

           if (sc_abandon_counter.ge.sc_abandon_etot_iter) then
             if (myid.eq.0) then
               write (use_unit,'(1X,A)') "**********************************************************************"
               write (use_unit,'(1X,A)') "*** "
               write (use_unit,'(1X,A)') "*** S.C.F. convergence problem!"
               write (use_unit,'(1X,A,F20.8,A)') &
                 "*** The total energy change between s.c.f. iterations has been above ", sc_abandon_etot_threshold, " eV/atom"
               write (use_unit,'(1X,A,I8,A)') &
                 "*** for ", sc_abandon_counter, " successive s.c.f. iterations."
               write (use_unit,'(1X,A)') &
                 "*** According to our predefined criteria, the s.c.f. cycle is most likely not going to "
               write (use_unit,'(1X,A)') &
                 "*** converge any more. We will stop it to prevent an unhelpful waste of computer time."
               write (use_unit,'(1X,A)') &
                 "*** To adjust this behaviour, check out the 'sc_abandon_etot' keyword in read_control.f90 ."
               write (use_unit,'(1X,A)') "*** "
               write (use_unit,'(1X,A)') "**********************************************************************"
             end if

             ! This is one place where we must execute a hard stop right now.
             ! Outsourcing the actual stop to a more graceful place has not
             ! worked in the past - if the s.c.f. cycle did not converge in this
             ! drastic way, we simply can not go on to produce any output.  The
             ! user can still override the entire check in control.in.
             call aims_stop_coll('SCF cycle not converged.', func)
           end if
        end if

        if (abort_requested.and.myid.eq.0) then
           write(use_unit,*)
           write(use_unit,*) '*** Abort request found. ***'
           write(use_unit,*) '*** SCF cycle will be stopped. ***'
        endif
        if (control_update_requested) then
            call read_control_update ()
        end if

        ! this is the end of the self-consistency cycle.

  end do SCF_LOOP

!===============================================================================
!=                           END SELF-CONSISTENCY LOOP                         =
!===============================================================================

      ! Everything after this point is post-processing to be done before handing
      ! control back to main and starting the "real" post-processed part of the
      ! calculation

      mbd_scf_converged= .true.

      if(allocated(hartree_potential_save)) then
         call aims_deallocate(hartree_potential_save,"hartree_potential_save")
      endif

      if (python_hooks%post_scf%registered) then
          call run_python_hook(python_hooks%post_scf)
      end if

      if (out_sorted_eigenvalue_list) then
           write(info_str,*) 'Starting output of sorted eigenvalues'
           call localorb_info(info_str)
           call output_sorted_eigenvalues(KS_eigenvalue, occ_numbers, chemical_potential)
      endif

      use_hartree_non_periodic_ewald = use_hartree_non_periodic_ewald_save

      ! In case we want to do a Boys localization at the end of the SCF, we should actually do this
      ! here so that we write the transformed restart info

      if (apply_boys_flag.and.(boys_sub_flag.eq.2)) then
        call get_dipolematrix(KS_eigenvalue, KS_eigenvector, &
                              KS_eigenvector_complex, occ_numbers, &
                              chemical_potential, partition_tab, &
                              l_shell_max)
      end if
!      if(output_level .eq. 'MD_light') output_priority = 1
      ! definitely keep checkpoint information here if requested, and if not already done 10 lines above
      if (keep_restart_info.and.restart_write) then
        if (use_restart_save_iterations) then
          if (.not.(mod(number_of_loops,restart_save_iterations).eq.0)) then
            call write_restart_info()
          end if
        else
          call write_restart_info()
        end if
      end if

      if (number_of_loops .ge. 50) warn_slow_scf_convergence = .true.
      total_number_of_loops = total_number_of_loops + number_of_loops

      ! Perform any post-processing that is required after every scf cycle:
      ! * scaled ZORA, if required
      ! * output of a band structure
      ! * Direct calculation of a binding energy
      !
      !!! Please add any other postprocessing right here in the future !!!

      require_post_prc = &
          ((flag_rel == REL_zora) .and. (force_potential == 0)) &
          .or. out_band .or. out_band_mulliken .or. use_qmmm & ! add mlk
          .or. out_hirshfeld .or. out_hirshfeld_iterative &
          .or. use_vdw_correction_hirshfeld &
          .or. use_mbd_old &
          .or. use_mbd_dev &
          .or. use_libmbd &
          .or. (use_mbd_std .and. (.not. mbd_self_consistent .or. forces_on)) &
          .or. use_ll_vdwdf .or. flag_compute_kinetic .or. use_meta_gga_post &
          .or. use_nlcorr_post .or. use_vdw_post .or. flag_energy_density &
          .or. use_nlcorr_post .or. use_vdw_post.or. flag_out_dielectric &
          .or. out_mommat .or. flag_out_dclimit .or. flag_flex_fermi &
          .or. flag_out_dipmat .or. flag_out_coulmat_ovl .or. flag_out_boys &
          .or. flag_out_absorption .or. flag_out_coulmat_lvl &
          .or. flag_out_greenwood .or. flag_out_ev_scalapack_hdf5 &
          .or. flag_out_dipmat_k_k &
          .or. use_meta_gga_printout &
          .or. flag_out_locpot_atom


      ! Only require post-processing if this is not a dummy restart to
      ! just restart the scf cycle though.
      require_post_prc = require_post_prc .and. (.not.restarting_scf)

      if (require_post_prc) then

      if (flag_out_locpot_atom) then
         open(unit=88, file="On-site_ESP.dat")
         write(88,'(2X,A,X,A)') "#species", "On-site_electrostatic_potential [eV]"
         !write(88,'(A,X,A,X,A,X,A,X,A)') "#species", "free_energy_contribution", "delta_energy_contribution_real" ,"delta_energy_contribution_recip", "Total"
         do i_atom = 1, n_atoms
            elec_tot(i_atom) = elec_free_atom(i_atom) + elec_delta_atom(i_atom)-average_potential
            !write(88,'(F15.8)') average_potential
            write(88,'(2X,A,2X,F15.8)') &
            trim(species_name(species(i_atom))),  elec_tot(i_atom)*hartree
            ! write(88,'(A,X, F15.8, X, F15.8, X, F15.8)')  &
           !trim(species_name(species(i_atom))), elec_free_atom(i_atom)*hartree, elec_delta_atom(i_atom)*hartree, elec_tot(i_atom)*hartree
         enddo
      endif

      !   DIELECTRIC/Kubo-Greenwood electronic transport starts here
    if(flag_out_greenwood) then
         write(info_str,'(6X,A,1X,I4)') "Momentum Matrix/Kubo-Greenwood post processing starts"
         call localorb_info (info_str)

         if(use_batch_permutation > 0) then
            call aims_stop('*** use_batch_permutation does not work for calculation of dielectric tensor')
         endif

         if(greenwood_method .eq. 1) then
              call get_full_matrix_greenwood(KS_eigenvalue, KS_eigenvector, &
                   KS_eigenvector_complex, occ_numbers, chemical_potential,&
                   partition_tab, l_shell_max, ep1, ep2)
         else if (greenwood_method .eq. 2) then
              call get_sparse_matrix_greenwood(KS_eigenvalue, KS_eigenvector, &
                   KS_eigenvector_complex, occ_numbers, chemical_potential, &
                   partition_tab, l_shell_max, ep1, ep2)
         else if (greenwood_method .eq. 3) then
              call get_full_matrix_greenwood(KS_eigenvalue, KS_eigenvector, &
                   KS_eigenvector_complex, occ_numbers, chemical_potential,&
                   partition_tab, l_shell_max, ep1, ep2)
              call get_sparse_matrix_greenwood(KS_eigenvalue, KS_eigenvector, &
                   KS_eigenvector_complex, occ_numbers, chemical_potential, &
                   partition_tab, l_shell_max, ep1, ep2)
         endif

    endif
          !   DIELECTRIC/Kubo-Greenwood electronic transport ends here
          if(out_mommat) then
             call get_momentummatrix(KS_eigenvalue, KS_eigenvector, &
                                     KS_eigenvector_complex, occ_numbers, &
                                     chemical_potential, partition_tab,&
                                     l_shell_max)
          endif

          if(flag_out_dipmat .and. .not. calculate_perturbative_soc ) then
             call get_dipolematrix(KS_eigenvalue, KS_eigenvector, &
                                     KS_eigenvector_complex, occ_numbers, &
                                     chemical_potential, partition_tab,&
                                    l_shell_max)
          endif
          if(flag_out_boys) then
             ! since we can't have periodic systems here,
             ! this determination of Emin and Emax should be okay.
             !
             CALL determine_occupied_states(n_states_occupied)
             Emin=(KS_eigenvalue(1, 1, 1)*hartree)-0.01
             Emax=(KS_eigenvalue(n_states_occupied, 1, 1))*hartree+0.01
             write(info_str,'(1X, A, E11.4, 1X, E11.4)') "| Energies that encompass all occupied orbitals: ", Emin, Emax
             call localorb_info ( info_str )
             call get_dipolematrix(KS_eigenvalue, KS_eigenvector, &
                                     KS_eigenvector_complex, occ_numbers, &
                                     chemical_potential, partition_tab,&
                                    l_shell_max)
          endif
          if(flag_out_dipmat_k_k) then
             call get_dipolemat_k_k(KS_eigenvalue, KS_eigenvector, &
                                     KS_eigenvector_complex, occ_numbers, &
                                     chemical_potential, partition_tab,&
                                    l_shell_max, k_k_method)
          endif

          if(flag_out_coulmat_ovl) then
             call get_coulelement_ovl(partition_tab, l_shell_max,&
                          KS_eigenvalue, KS_eigenvector,KS_eigenvector_complex&
                          )
          endif

          if(flag_out_coulmat_lvl) then
            call get_coulelement_lvl_v0(KS_eigenvalue,KS_eigenvector,&
                                        KS_eigenvector_complex,.False.)
          endif
          if(flag_out_ev_scalapack_hdf5) then
            call get_ev_scalapack_hdf5 (KS_eigenvalue, eigenvec_complex, &
                                   ovlp_complex_stored, l_row, l_col, mxld, &
                                   mxcol,my_k_point, chemical_potential)
          endif
        if(output_level .eq. 'normal') output_priority = 0

        enough_walltime_left =  time_check_postprocessing().and.enough_walltime_left

        if (enough_walltime_left) then

          call cpu_time(time_postprocessing_start)

          if ( (flag_rel.eq.REL_zora) .and. (force_potential.eq.0) .and. (.not. use_geo_relaxation) ) then
            ! require full scaled ZORA corrections to eigenvalues

            write(info_str,'(A)') ''
            call localorb_info ( info_str, use_unit )

            write(info_str,'(A)') "------------------------------------------------------------"
            call localorb_info ( info_str, use_unit )

            write(info_str,'(10X,A)') &
              "Post-processing: scaled ZORA corrections to eigenvalues and total energy."
            call localorb_info ( info_str, use_unit )
            write(info_str,'(A)') ''
            call localorb_info ( info_str, use_unit )

            call date_and_time(cdate, ctime)
            write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, ", Time     :  ", ctime
            call localorb_info ( info_str, use_unit )

            write(info_str,'(A)') "------------------------------------------------------------"
            call localorb_info ( info_str, use_unit )

            call get_timestamps(time_scaled_zora, clock_time_scaled_zora)

            if(flag_KH_post_correction)then


              call relativistic_correction_term_p1 &
                    ( KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers, partition_tab,  &
                    rho, rho_gradient, kinetic_density, hartree_potential, l_shell_max )

            else

              call allocate_scaled_zora_transform

              call integrate_scaled_zora_transf_p2( &
                  rho, rho_gradient, kinetic_density, hartree_potential,    &
                  partition_tab, l_shell_max )


              call evaluate_scaled_zora_transf_p1( &
                  KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue)

              call deallocate_scaled_zora_transform

            end if

            if(force_occupation_projector .and. start_force_occ) then
               if(mu_method == 0) then ! zeroin
                  call get_occupation_numbers_occ_p0(KS_eigenvalue,n_electrons,&
                          .true.,occ_numbers,chemical_potential)
               else ! bisection
                  call aims_elsi_occ(n_electrons,n_states,n_spin,n_k_points,&
                          k_weights,KS_eigenvalue,occ_numbers,&
                          chemical_potential)
               endif
            else if (fixed_spin_moment .or. use_hf_multiplicity) then
               ! occupation numbers must be updated per spin channel, separate
               ! chemical potentials for each channel will be returned.
               if(mu_method == 0) then ! zeroin
                  call get_occupation_numbers_fsm(KS_eigenvalue,occ_numbers,&
                          chemical_potential_spin)
               else ! bisection
                  do i_spin = 1,n_spin
                     call aims_elsi_occ(fixed_spin_moment_electrons(i_spin),&
                             n_states,1,n_k_points,k_weights,&
                             KS_eigenvalue(1:n_states,i_spin,1:n_k_points),&
                             occ_numbers(1:n_states,i_spin,1:n_k_points),&
                             chemical_potential_spin(i_spin))
                  enddo
               endif
            else
               if(mu_method == 0) then ! zeroin
                  call get_occupation_numbers_p0(KS_eigenvalue,n_electrons,&
                          .true.,occ_numbers,chemical_potential)
               else ! bisection
                  call aims_elsi_occ(n_electrons,n_states,n_spin,n_k_points,&
                          k_weights,KS_eigenvalue,occ_numbers,&
                          chemical_potential)
               endif

               ! Set chemical_potential_spin here (same for both channels)
               chemical_potential_spin(:) = chemical_potential
            endif


            call output_real_eigenfunctions &
              ( KS_eigenvalue, KS_eigenvector, occ_numbers)


            if(n_periodic > 0  )then
               if(use_dipole_correction .or. calculate_work_function)then
                  call output_energy_levels &
                  ( KS_eigenvalue, occ_numbers, chemical_potential, chemical_potential_spin, converged )
               end if
            end if


            call get_ev_sum_p0 &
            ( occ_numbers, KS_eigenvalue, av_core_shift, ev_sum, &
              ev_sum_shifted )

            call get_entropy_correction_p1(KS_eigenvalue, occ_numbers, &
               chemical_potential, chemical_potential_spin, entropy_correction)

            if (output_level.eq.'MD_light') output_priority = 1
            call get_total_energy &
              ( ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
              en_ion_embed, en_density_embed, &
              en_vdw, en_pot_vdw, en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, en_elec_free, &
              en_elec_delta, fock_energy, &
              hartree_energy_free, hartree_delta_energy, &
              hartree_multipole_correction, &
              entropy_correction, penalty_energy, total_energy &
              )
            if (output_level.eq.'MD_light') output_priority = 2

            call get_times(time_scaled_zora, clock_time_scaled_zora, &
            &              tot_time_scaled_zora, tot_clock_time_scaled_zora)

            write(info_str,'(2X,A)') &
              "End evaluation of scaled ZORA corrections."
            call localorb_info(info_str, use_unit)
            call output_times(deffmt, "Time for scaled ZORA corrections", &
            &                 time_scaled_zora, clock_time_scaled_zora,OL_norm)
            write(info_str,'(A)') &
              "------------------------------------------------------------"
            call localorb_info(info_str, use_unit)

          end if ! scaled ZORA corrections


          call cpu_time(rtime)
          call update_time_max_postprocessing(rtime - time_postprocessing_start)

!    QM/MM forces on multipoles
          if (use_qmmm) then
            write(info_str,'(6X,A,1X,I4)') "Calculating Hellman-Feynman forces on external charges ..."
            call localorb_info ( info_str, use_unit )
            call evaluate_ext_charge_hellman_feynman_forces_p0(rho, partition_tab, ext_charge_forces)
            call sync_ext_charge_forces(ext_charge_forces)
          endif

        if ( &
            use_vdw_correction_hirshfeld .or. use_mbd_old &
            .or. use_mbd_dev .or. use_mbd_std .or. use_libmbd &
        ) then
           call get_timestamps ( time_vdw, clock_time_vdw )
        end if

!    If Hirshfeld analysis requested, fill in charge density on grid points that are not covered by "normal" partition_tab
          if ( &
              out_hirshfeld .or. out_hirshfeld_iterative .or. use_vdw_correction_hirshfeld &
              .or. use_mbd_old .or. use_mbd_dev .or. use_mbd_std .or. use_libmbd &
          ) then
            ! Prior to Hirshfeld analysis, must make sure that the density was calculated on
            ! all grid points - previously, the density was only updated on the grid points
            ! with a non-zero partition table.
            ! These subroutines must be kept exacly in sync with the usual density update routines.
            ! try to get new maxes....
            call get_n_max_for_densities( partition_tab, hartree_partition_tab)

            call kweight_occs('scf_solver:missing', occ_numbers)

            if (use_density_matrix) then
              if (use_local_index .or. transport_calculation .or. transport_lead_calculation) then
                call aims_allocate( my_density_matrix, n_hamiltonian_matrix_size, "my_density_matrix" )
                call update_missing_density_densmat &
                       ( KS_eigenvector, KS_eigenvector_complex, &
                       occ_numbers, partition_tab, hartree_partition_tab, rho, &
                       l_shell_max, my_density_matrix )
                call aims_deallocate( my_density_matrix, "my_density_matrix" )
              else
                call update_missing_density_densmat &
                       ( KS_eigenvector, KS_eigenvector_complex, &
                        occ_numbers, partition_tab,  hartree_partition_tab, rho, &
                        l_shell_max, hamiltonian(1,1) )
              end if
            else
              call update_missing_density_orbital &
                    ( KS_eigenvector, KS_eigenvector_complex, &
                    KS_eigenvalue, occ_numbers, &
                    partition_tab, hartree_partition_tab, rho, &
                    l_shell_max )
            end if

            call de_kweight_occs('scf_solver:missing', occ_numbers)

          end if

!         Van der Waals correction based on Hirshfeld volume partitioning
          if (use_vdw_correction_hirshfeld) then
!            if(output_level .eq. 'MD_light') output_priority = 2
            write(info_str,'(2X,A)') &
              "Evaluating non-empirical van der Waals correction (Tkatchenko/Scheffler 2009)."
            call localorb_info ( info_str,use_unit,'(A)', OL_norm )
            call hirshfeld_analysis( )
            vdw_forces = 0.0
            call calc_vdw_hirshfeld(en_vdw,vdw_forces)

            if (output_level.eq.'MD_light') output_priority = 1
            call get_total_energy &
              ( ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
              en_ion_embed, en_density_embed, &
              en_vdw, en_pot_vdw, en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, en_elec_free, &
              en_elec_delta, fock_energy, &
              hartree_energy_free, hartree_delta_energy, &
              hartree_multipole_correction, &
              entropy_correction, penalty_energy, total_energy &
              )
            if (output_level.eq.'MD_light') output_priority = 2

            if (forces_on) then
              call get_total_forces &
                  (ionic_forces, pulay_forces, hellman_feynman_forces, &
                  multipole_forces, gga_forces, gga_forces_on, &
                  total_forces, vdw_forces, d_exx_ene, EXX_forces_on, pseudocore_forces, &
                  nlcc_forces, nlcc_forces_on, Gnonmf_forces, Gnonmf_forces_on)
              ! FK: If vdW is used, the full stress tensor with vdW correction will be printed here.
              if (use_analytical_stress) then
                call get_total_analytical_stress( AS_pulay_stress_on, (en_xc - en_pot_xc), fock_energy, &
                     previous_stress, diff_stress, .false. )
              end if
            end if

          else if(use_mbd_old) then
            ! VB: Verified in read_control.f90 that vdw_correction_hirshfeld
            !     and use_mbd_old are mutually exlusive

            if(.NOT.forces_on) then
              call hirshfeld_analysis ( )

              en_vdw=0.d0
              vdw_forces =0.d0
              call mbd_rsSCS_energy_and_forces (en_vdw,vdw_forces)

              call get_total_energy &
                ( ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
                en_ion_embed, en_density_embed, &
                en_vdw, en_pot_vdw, en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, en_elec_free, &
                en_elec_delta, fock_energy, &
                hartree_energy_free, hartree_delta_energy, &
                hartree_multipole_correction, &
                entropy_correction, penalty_energy, total_energy &
                )

            else
              call hirshfeld_analysis ( )
              en_vdw=0.d0
              vdw_forces =0.d0
              call mbd_rsSCS_energy_and_forces (en_vdw,vdw_forces)

              if (output_level.eq.'MD_light') output_priority = 1

              call get_total_energy &
                ( ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
                en_ion_embed, en_density_embed, &
                en_vdw, en_pot_vdw, en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, en_elec_free, &
                en_elec_delta, fock_energy, &
                hartree_energy_free, hartree_delta_energy, &
                hartree_multipole_correction, &
                entropy_correction, penalty_energy, total_energy &
                )
              if (output_level.eq.'MD_light') output_priority = 2


              call get_total_forces &
                  (ionic_forces, pulay_forces, hellman_feynman_forces, &
                  multipole_forces, gga_forces, gga_forces_on, &
                  total_forces, vdw_forces, d_exx_ene, EXX_forces_on, pseudocore_forces, &
                  nlcc_forces, nlcc_forces_on, Gnonmf_forces, Gnonmf_forces_on)

            endif ! check force_on

          else if (use_mbd_std) then
              call hirshfeld_analysis()
              if (forces_on .and. use_pulay) then
                  if (use_local_index) then
                      call aims_allocate( my_density_matrix, n_hamiltonian_matrix_size*n_spin, "my_density_matrix" )
                      call eval_hirsh_vol_pulay_deriv_dm( &
                          KS_eigenvector, KS_eigenvector_complex, occ_numbers, &
                          my_density_matrix &
                      )
                      call aims_deallocate( my_density_matrix, "my_density_matrix" )
                  else
                      if (run_numerical_pulay_forces) then
                          call numerical_hirsh_pulay_forces( &
                              KS_eigenvector, KS_eigenvector_complex, &
                              occ_numbers, hamiltonian &
                          )
                      end if
                      call eval_hirsh_vol_pulay_deriv_dm( &
                          KS_eigenvector, KS_eigenvector_complex, occ_numbers, &
                          hamiltonian &
                      )
                  endif
              end if
              call mbd_std_calculate(rho, en_vdw, vdw_forces, as_vdw_stress)
              if (.not. mbd_self_consistent) then
                  call get_total_energy( &
                      ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
                      en_ion_embed, en_density_embed, en_vdw, en_pot_vdw, &
                      en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, en_elec_free, &
                      en_elec_delta, fock_energy, hartree_energy_free, &
                      hartree_delta_energy, hartree_multipole_correction, &
                      entropy_correction, penalty_energy, total_energy &
                  )
              end if
              if (forces_on) then
                  call get_total_forces( &
                      ionic_forces, pulay_forces, hellman_feynman_forces, &
                      multipole_forces, gga_forces, gga_forces_on, &
                      total_forces, vdw_forces, d_exx_ene, EXX_forces_on, &
                      pseudocore_forces, nlcc_forces, nlcc_forces_on, &
                      Gnonmf_forces, Gnonmf_forces_on &
                  )
                  if (use_analytical_stress) then
                      call get_total_analytical_stress( &
                          as_pulay_stress_on, en_xc-en_pot_xc, fock_energy, &
                          previous_stress, diff_stress, .false. &
                      )
                  end if
              end if

          else if (use_mbd_dev) then
              call hirshfeld_analysis()
              vdw_forces(:, :) = 0.d0
              call mbd_dev_calculate( &
                  hirshfeld_volume(:n_occ_atoms), &
                  en_vdw, &
                  vdw_forces(:, :n_occ_atoms), as_vdw_stress)
              if (output_level == 'MD_light') output_priority = 1
              call get_total_energy( &
                  ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
                  en_ion_embed, en_density_embed, en_vdw, en_pot_vdw, &
                  en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, en_elec_free, &
                  en_elec_delta, fock_energy, hartree_energy_free, &
                  hartree_delta_energy, hartree_multipole_correction, &
                  entropy_correction, penalty_energy, total_energy)
              if (output_level == 'MD_light') output_priority = 2
              if (forces_on) then
                  call get_total_forces( &
                      ionic_forces, pulay_forces, hellman_feynman_forces, &
                      multipole_forces, gga_forces, gga_forces_on, &
                      total_forces, vdw_forces, d_exx_ene, EXX_forces_on, &
                      pseudocore_forces, nlcc_forces, nlcc_forces_on, &
                      Gnonmf_forces, Gnonmf_forces_on &
                  )
                  if (use_analytical_stress) then
                      call get_total_analytical_stress( &
                          AS_pulay_stress_on, (en_xc-en_pot_xc), &
                          fock_energy, &
                          previous_stress, diff_stress, .false.)
                  end if
              end if

          else if (use_libmbd) then
              call localorb_info('')
              call localorb_info('  Entering Libmbd.')
              call xml_open('libmbd')
              call localorb_info('  Libmbd: VdW method: ' // mbd_input%method)
              call xml_elem('method', mbd_input%method)
              select case (mbd_input%method)
              case ('mbd-rsscs', 'ts')
                  call localorb_info('  Libmbd: Performing Hirshfeld analysis...')
                  call localorb_info('')
                  call hirshfeld_analysis()
                  call mbd_calc%update_vdw_params_from_ratios(hirshfeld_volume)
              case ('mbd-nl')
                  call localorb_info('  Libmbd: Evaluating kinetic energy density...')
                  call localorb_info('')
                  call kweight_occs('scf_solver:mbd-nl', occ_numbers)
                  call calculate_kinetic_density( &
                      rho, rho_gradient, kinetic_density, &
                      hartree_partition_tab, partition_tab, l_shell_max, &
                      KS_eigenvector, KS_eigenvector_complex, occ_numbers &
                  )
                  call de_kweight_occs('scf_solver:mbd-nl', occ_numbers)
                  call localorb_info('')
                  allocate (volumes_hirsh_free(n_species))
                  allocate (alpha_0_vv_free(n_species), C6_vv_free(n_species))
                  call evaluate_free_atom_quantities( &
                      volumes_hirsh_free, alpha_0_vv_free, C6_vv_free &
                  )
                  allocate (volumes_hirsh(n_atoms))
                  allocate (alpha_0_vv_nm(n_atoms), C6_vv_nm(n_atoms))
                  call localorb_info('  Libmbd: Performing Hirshfeld analysis...')
                  call run_hirshfeld( &
                      rho, volumes_hirsh, rho_grad_grid=rho_gradient, &
                      kinetic_density_grid=kinetic_density, &
                      alpha_0_vv_nm=alpha_0_vv_nm, C6_vv_nm=C6_vv_nm &
                  )
                  call mbd_calc%update_vdw_params_nl( &
                      alpha_0_vv_nm/alpha_0_vv_free(species), &
                      C6_vv_nm/C6_vv_free(species) &
                  )
              end select
              call localorb_info('  Libmbd: Calculating MBD energy...')
              call mbd_calc%update_coords(coords)
              if (n_periodic > 0) call mbd_calc%update_lattice_vectors(lattice_vector)
              call mbd_calc%evaluate_vdw_method(en_vdw)
              call mbd_calc%get_exception(exc_code, exc_origin, exc_msg)
              if (exc_code > 0) then
                  call localorb_info('  Error in MBD routine')
                  call aims_stop(exc_msg, exc_origin)
              end if
              call localorb_info('  Libmbd: Evaluated energy: ' // tostr(en_vdw))
              call xml_elem('energy', en_vdw)
              if (output_level == 'MD_light') output_priority = 1
              call get_total_energy( &
                  ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
                  en_ion_embed, en_density_embed, en_vdw, en_pot_vdw, &
                  en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, en_elec_free, &
                  en_elec_delta, fock_energy, hartree_energy_free, &
                  hartree_delta_energy, hartree_multipole_correction, &
                  entropy_correction, penalty_energy, total_energy &
              )
              if (output_level == 'MD_light') output_priority = 2
              if (forces_on) then
                  call mbd_calc%get_gradients(vdw_forces)
                  vdw_forces = -vdw_forces
                  call get_total_forces( &
                      ionic_forces, pulay_forces, hellman_feynman_forces, &
                      multipole_forces, gga_forces, gga_forces_on, &
                      total_forces, vdw_forces, d_exx_ene, EXX_forces_on, &
                      pseudocore_forces, nlcc_forces, nlcc_forces_on, &
                      Gnonmf_forces, Gnonmf_forces_on &
                  )
                  if (use_analytical_stress) then
                      call mbd_calc%get_lattice_derivs(as_vdw_stress)
                      as_vdw_stress = ( &
                          matmul(as_vdw_stress, transpose(lattice_vector)) &
                          - matmul(vdw_forces, transpose(coords)) &
                      )
                      call get_total_analytical_stress( &
                          AS_pulay_stress_on, en_xc-en_pot_xc, &
                          fock_energy, previous_stress, diff_stress, .false. &
                      )
                  end if
              end if
              if (mbd_input%rpa_orders) then
                  call mbd_calc%get_rpa_orders(mbd_rpa_orders)
                  call xml_elem('rpa_orders', mbd_rpa_orders)
                  deallocate (mbd_rpa_orders)
              end if
              call xml_close()
              call localorb_info('')
              call localorb_info('  Leaving Libmbd.')
              call localorb_info('')

          else if (out_hirshfeld_always) then

            if (out_hirshfeld) then
               call hirshfeld_analysis ( )
            end if
            if (out_hirshfeld_iterative) then
               call hirshfeld_analysis_iterative ( )
            end if

          endif ! check vdw_correction_hirshfeld or mbd flag

          if ( &
              use_vdw_correction_hirshfeld .or. use_mbd_old &
              .or. use_mbd_dev .or. use_mbd_std .or. use_libmbd &
          ) then
             ! Timings for vdw_correction_hirshfeld.
             ! For certain solids, the time can matter, but was not assessed.
             ! The timings are now written out. At a later stage, the
             ! output priority could be switched to 'OL_low'.
             call get_times(time_vdw, clock_time_vdw, &
             &              tot_time_vdw, tot_clock_time_vdw)
             write(info_str,'(2X,A)') &
                   "Timings for VdW correction (Tkatchenko-Scheffler 2009 or later):"
             call localorb_info(info_str,use_unit,'(A)',OL_norm)
             write(info_str,'(2X,A,F12.3,A,F12.3,A)') &
                   "| CPU time: ", time_vdw, " s, Wall clock time: ", clock_time_vdw, " s."
             call localorb_info(info_str,use_unit,'(A)',OL_norm)
          end if

          !   LL VDW (2010). Nonlocal correlation post-processed correction. SAG
          if(use_nlcorr_post)then

             write(use_unit,*)" VDW post-processing starts"

             if(associated(root_full%branch)) then
                count2 = 0
                call deallocate_tree(root_full,count2)
             endif

             call gen_master_cube()

             call integrate_nlcorr(en_nlcorr)

             vdw_total_energy = total_energy + en_nlcorr

             if(myid.eq.0) then
                write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
                     "Total energy with nonlocal correlation:", vdw_total_energy, &
                     "Ha,", vdw_total_energy*hartree, " eV"
             endif

          endif



          !LL VDW (2010). pbe_vdw or revpbe_vdw post-processed correction. SAG
          if(use_vdw_post)then
             write(info_str,'(2X,A)') "VDW post-processing starts"
             call localorb_info(info_str,use_unit,'(A)', OL_norm)

             if(.not.associated(root_full%branch)) then
                call gen_master_cube()
             endif

             write(info_str,'(2X,A)') "calculating post-processed xc energy"
             call localorb_info(info_str,use_unit,'(A)', OL_norm)
             call integrate_post_xc_energy &
                  ( partition_tab, &
                  rho,rho_gradient,  &
                  kinetic_density,&
                  en_post_xc, &
                  x_post_energy, &
                  c_post_energy )

             vdw_total_energy = total_energy - en_xc + en_post_xc

            if(myid.eq.0) then
               if(flag_post_xc.eq.3)then
                  write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
                       "PBE_VDW corrected total energy", vdw_total_energy, &
                       " Ha,", vdw_total_energy*hartree, " eV"
                  write(use_unit,*)" "
               endif
               if(flag_post_xc.eq.4)then
                  write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
                       "REVPBE_VDW corrected total energy", vdw_total_energy, &
                       " Ha,", vdw_total_energy*hartree, " eV"
                  write(use_unit,*)" "
               endif
            endif

         endif




!    LL van der Waals functional
          if (use_ll_vdwdf) then
            write(info_str,'(6X,A,1X,I4)') "Calculating charge density projected on even grids ..."
            call localorb_info ( info_str,use_unit,'(A)', OL_norm  )
            call kweight_occs('scf_solver:ll_vdwdf', occ_numbers)
            call obtain_charge_density( )
            call de_kweight_occs('scf_solver:ll_vdwdf', occ_numbers)
            ll_vdw_forces = 0.0
            write(info_str,'(6X,A,1X,I4)') "Calculating LL van der Waals energy term ..."
            call localorb_info ( info_str,use_unit,'(A)', OL_norm  )
            call calc_ll_vdw(en_ll_vdw,en_ll_vdw_err,en_lda_c,en_pbe_c)
            write(info_str,'(6X,A,1X,I4)') "Synchronizing LL van der Waals energy term and errors ..."
            call localorb_info ( info_str,use_unit,'(A)', OL_norm  )
            call sync_ll_vdwdf(en_ll_vdw,en_ll_vdw_err,ll_vdw_forces,en_lda_c,en_pbe_c)
!            do i_atom=1, n_atoms
!              total_forces(:,i_atom) = total_forces(:,i_atom) + ll_vdw_forces(:,i_atom)
!            end do
            if (output_level.eq.'MD_light') output_priority = 1
            call get_total_energy &
              ( ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
              en_ion_embed, en_density_embed, &
              en_vdw, en_pot_vdw, en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, en_elec_free, &
              en_elec_delta, fock_energy, &
              hartree_energy_free, hartree_delta_energy, &
              hartree_multipole_correction, &
              entropy_correction, penalty_energy, total_energy &
              )
            if (output_level.eq.'MD_light') output_priority = 2
         endif

! Computation of meta_gga
          if(use_meta_gga_post) then
            write(info_str,'(2X,A)') "Meta-GGA post processing starts"
            call localorb_info ( info_str, use_unit )

            call kweight_occs('scf_solver:meta', occ_numbers)

!      First calculate the kinetic density
            call localorb_info( &
            "Calculating meta-GGA energy as post processing, using DM", use_unit,'(2X,A)', OL_norm )

            call calculate_kinetic_density &
              (rho, rho_gradient, kinetic_density, &
              hartree_partition_tab, &
              partition_tab, l_shell_max, &
              KS_eigenvector, KS_eigenvector_complex, occ_numbers)

!      Then integrate the post proc. xc energy
            call integrate_post_xc_energy &
               ( partition_tab, &
               rho,rho_gradient,  &
               kinetic_density, &
               en_post_xc, &
               x_post_energy, &
               c_post_energy)

            if (use_hybrid_meta) then
               call prepare_corr_energy_calc()
               call qpe_calculation()
            endif

            call de_kweight_occs('scf_solver:meta', occ_numbers)

!       Finally subtract old xc energy and add new
            meta_gga_total_energy = total_energy - en_xc + en_post_xc
            if (use_hybrid_meta) then
               meta_gga_total_energy=meta_gga_total_energy+en_hf
            endif

            if(myid.eq.0) then
              write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
              "| Meta-gga total energy", meta_gga_total_energy, &
              " Ha,", meta_gga_total_energy*hartree, " eV"
!              write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
!              "  DFT/HF total energy",  (total_energy+en_hf), &
!              " Ha,", (total_energy+en_hf)*hartree, " eV"
!
              write(use_unit,*)"------------------------------------------", &
              "----------------------------------"
              write(use_unit,*)

            endif

         endif
         !if(use_meta_gga_printout) then
         !  call kweight_occs('scf_solver:meta', occ_numbers)
         !  call calculate_kinetic_density &
         !    (rho, rho_gradient, kinetic_density, &
         !    hartree_partition_tab, &
         !    partition_tab, l_shell_max, &
         !    KS_eigenvector, KS_eigenvector_complex, occ_numbers)
         !  call de_kweight_occs('scf_solver:meta', occ_numbers)
         !  use_meta_gga_old = use_meta_gga
         !  use_meta_gga = .true.
         !end if
! Computation of the kinetic energy matrix
          if ( (flag_compute_kinetic).or.( (flag_chetty_martin_energy_density) ) ) then
            write(info_str,'(6X,A)') "Calculating kinetic energy matrix..."
            call localorb_info ( info_str, use_unit )
            if ( n_periodic .gt. 0 )  then
              write(info_str,'(6X,A)') " *** Warning: Kinetic energy computation seems not to be working for PBC"
              call localorb_info ( info_str, use_unit )
            end if
            call kweight_occs('scf_solver:kinetic', occ_numbers)
            call get_kinetic_energy(rho, partition_tab, l_shell_max, &
                                    KS_eigenvector, KS_eigenvector_complex, &
                                    occ_numbers, kinetic_energy &
                                  )
            call de_kweight_occs('scf_solver:kinetic', occ_numbers)
          endif
       end if ! enough wall time
        if(output_level .eq. 'normal') output_priority = 1
     end if ! post-processing tasks

      if (out_dipole_moment) then
        call output_dipole_moment ( )
      end if

      if (out_quadrupole_moment) then
        call output_quadrupole_moment ( )
      end if

!     If band structure output was requested based on the standard k-point
!     grid used in the s.c.f. cycle, perform this output here - after all modifications
!     to the eigenvalues are done.
      if (out_band_during_scf) then
         call output_bands_during_scf
      end if


!     Write final self-consistency data, if required:

      if (force_potential.eq.0) then
        if (out_v_hartree) then
          call output_hartree_pot (multipole_moments)
        end if

        if (out_rho_multipole) then
          call output_rho_multipole
        end if

        if (out_density) then

          call output_density_p1( rho, partition_tab )

        end if

        if (out_v_eff) then
           if (.not.use_mpi) then
              !             write current effective potential on integration grid to file

              call output_potential_p1 ( hartree_potential )


           else
              write(use_unit,*) '* Output of potential is not supported'
              write(use_unit,*) '* in MPI runs. Aborting.'
              stop
           end if
        end if
        if (out_v_eff_new) then
           call output_potential_p1_new ( hartree_potential, "consistent")
        end if

     end if

     if (use_fo_potential.and.(.not.fo_file_exists)) then
        call fodft_out_potential( hartree_potential, 0 )
     elseif (use_fo_potential.and.(fo_file_exists)) then
        call fodft_out_potential( hartree_potential, 1, local_fo_potential)
     end if

    ! Assembly of atomic stress
    ! CC: Evaluate stress per atom when AIMS is used as a library for MD
    !     For ``real'' GK MD calculations, J is computed in
    !     output_molecular_dynamics using actual velocities
    if (compute_heat_flux .and. .not. use_molecular_dynamics) then
      call assemble_heat_flux()
    end if

      !     write energy and forces in a compact form in standard-output file
      if (force_potential.eq.0) then
        call output_energy_and_forces()
        ! FlK: if atomic stress was computed, print to standard output file
        if (compute_heat_flux .and. .not. use_molecular_dynamics) then
          call print_atomic_stress()
        end if
      end if

      ! FK: Output summary for external pressure
      if (flag_external_pressure) then
        if (use_numerical_stress) then
          if (num_stress_finished) then
            call external_pressure_summary(total_energy, external_pressure, cell_volume, numerical_stress_tensor(1:3,1:3))
          end if
        else if (use_analytical_stress) then
          call external_pressure_summary(total_energy, external_pressure, cell_volume, analytical_stress_tensor(1:3,1:3))
        else
          call external_pressure_summary(total_energy, external_pressure, cell_volume)
        end if
      end if

      ! write density of states if requested
      if (out_dos) then
        if (.not.pert_dos_on) then
          ! pert_dos_on would mean that a post-processing calculation
          ! of the density of states with a denser k-point grid was requested. If that is the case,
          ! the present information on eigenvalues and eigenvectors would be overwritten.
          ! So, we here only plot density of states information (per scf cycle) if it can be
          ! derived directly from the current eigenvectors, without modifying them.
          ! Look for the "post-processed" density of states with a denser k-grid at the end
          ! of main.f90 instead.
          if (calculate_perturbative_soc) then
            ! The SOC-perturbed version will be calculated later, here calculate
            ! the non-SOC-perturbed version
            call output_KS_dos( KS_eigenvalue, occ_numbers, chemical_potential, n_electrons, &
                 "KS_DOS_total.dat.no_soc", "KS_DOS_total_raw.dat.no_soc" )
          else
            call output_KS_dos( KS_eigenvalue, occ_numbers, chemical_potential, n_electrons, &
                 "KS_DOS_total.dat", "KS_DOS_total_raw.dat" )
          end if
        end if
      end if
      ! write density of states if requested
      if (out_dos_tetrahedron) then
        if (.not.pert_dos_on) then
          ! pert_dos_on would mean that a post-processing calculation
          ! of the density of states with a denser k-point grid was requested. If that is the case,
          ! the present information on eigenvalues and eigenvectors would be overwritten.
          ! So, we here only plot density of states information (per scf cycle) if it can be
          ! derived directly from the current eigenvectors, without modifying them.
          ! Look for the "post-processed" density of states with a denser k-grid at the end
          ! of main.f90 instead.
          if (calculate_perturbative_soc) then
            ! The SOC-perturbed version will be calculated later, here calculate
            ! the non-SOC-perturbed version
            call output_KS_dos_tetrahedron( KS_eigenvalue, occ_numbers, chemical_potential, n_electrons, &
                 "KS_DOS_total_tetrahedron.dat.no_soc", "KS_DOS_total_raw_tetrahedron.dat.no_soc" )
          else
            call output_KS_dos_tetrahedron( KS_eigenvalue, occ_numbers, chemical_potential, n_electrons, &
                 "KS_DOS_total_tetrahedron.dat", "KS_DOS_total_raw_tetrahedron.dat" )
          end if
        end if
      end if


      !    write the atom-projected KS density of states
      if(out_species_dos) then
        call get_species_projected_dos('KS',1, n_states,KS_eigenvalue)
      endif

      ! do double-hybrid density functional calculations, IGOR
      if ((use_dftpt2 .or. use_dftpt2_qpe) .and. .not. use_meta_gga_post) then
          total_energy_old = total_energy
          en_xc_old        = en_xc
          if(use_periodic_hf) en_hf = exx_ene * hybrid_coeff
          dftpt2_total_energy_dft  = total_energy - en_xc + 2.0d0 * en_hf
          call dftpt2_dft_part_energy ( partition_tab, &
                                        rho, rho_gradient, &
                                        kinetic_density, &
                                        xc_energy &
                                        )
          !en_hf  = en_hf/hybrid_coeff*dftpt2_Ex_hf
          en_xc  = xc_energy
          dftpt2_total_energy_dft = dftpt2_total_energy_dft + en_xc
          if (myid .eq. 0) then
              write(use_unit,'(2X,A)') "---------------------------------------- "
              write(use_unit,'(2X,A)') "Double-Hybrid Density Functional is used "
              write(use_unit,'(2X,A)') "---------------------------------------- "
              if (flag_dftpt2_dft_part .eq. 30) then
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_total(B3LYP)        :", total_energy_old, " Ha     ", (total_energy_old)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_xc(B3LYP)           :", en_xc_old, " Ha     ", (en_xc_old)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_x_hf(B3LYP)         :", en_hf, " Ha     ", (en_hf)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_xc(XYG3/DFT)        :", en_xc, " Ha     ", (en_xc)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_total(XYG3/DFT)     :", dftpt2_total_energy_dft, " Ha     ", (dftpt2_total_energy_dft)*hartree, " eV"
                  write(use_unit,'(2X,A)') "------------------------------------ "
              elseif (flag_dftpt2_dft_part .eq. 31) then
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_total(PBE0)         :", total_energy_old, " Ha     ", (total_energy_old)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_xc(PBE0)            :", en_xc_old, " Ha     ", (en_xc_old)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_x_hf(PBE0)          :", en_hf, " Ha     ", (en_hf)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_xc(xDH-PBE0/DFT)    :", en_xc, " Ha     ", (en_xc)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_total(xDH-PBE0/DFT) :", dftpt2_total_energy_dft, " Ha     ", (dftpt2_total_energy_dft)*hartree, " eV"
                  write(use_unit,'(2X,A)') "------------------------------------ "
              elseif (flag_dftpt2_dft_part .eq. 32) then
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_total(B3LYP)        :", total_energy_old, " Ha     ", (total_energy_old)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_xc(B3LYP)           :", en_xc_old, " Ha     ", (en_xc_old)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_x_hf(B3LYP)         :", en_hf, " Ha     ", (en_hf)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_xc(XYGJOS/DFT)      :", en_xc, " Ha     ", (en_xc)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_total(XYGJOS/DFT)   :", dftpt2_total_energy_dft, " Ha     ", (dftpt2_total_energy_dft)*hartree, " eV"
                  write(use_unit,'(2X,A)') "------------------------------------ "
              elseif (flag_dftpt2_dft_part .eq. 33) then
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_total(PBE0)         :", total_energy_old, " Ha     ", (total_energy_old)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_xc(PBE0)            :", en_xc_old, " Ha     ", (en_xc_old)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_x_hf(PBE0)          :", en_hf, " Ha     ", (en_hf)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_xc(ZRPS/DFT)        :", en_xc, " Ha     ", (en_xc)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_total(ZRPS/DFT)     :", dftpt2_total_energy_dft, " Ha     ", (dftpt2_total_energy_dft)*hartree, " eV"
                  write(use_unit,'(2X,A)') "------------------------------------ "
              elseif (flag_dftpt2_dft_part .eq. 34) then
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_total(B3LYP)        :", total_energy_old, " Ha     ", (total_energy_old)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_xc(B3LYP)           :", en_xc_old, " Ha     ", (en_xc_old)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_x_hf(B3LYP)         :", en_hf, " Ha     ", (en_hf)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_xc(XYG5/DFT)        :", en_xc, " Ha     ", (en_xc)*hartree, " eV"
                  write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
                   "E_total(XYG5/DFT)     :", dftpt2_total_energy_dft, " Ha     ", (dftpt2_total_energy_dft)*hartree, " eV"
                  write(use_unit,'(2X,A)') "------------------------------------ "
              endif
          endif
          !en_hf  = en_hf*hybrid_coeff/dftpt2_Ex_hf
          en_xc  = en_xc_old
      else
          dftpt2_total_energy_dft = total_energy
      endif
      ! ============================================================
      ! Igor, to print out exchange-correlation components for DFT
      ! methods specified in control.f90
      ! ============================================================
      ! for meta-GGAs, preparing kinetic_density
      if(use_meta_gga_printout) then
        call kweight_occs('scf_solver:meta', occ_numbers)
        call calculate_kinetic_density &
          (rho, rho_gradient, kinetic_density, &
          hartree_partition_tab, &
          partition_tab, l_shell_max, &
          KS_eigenvector, KS_eigenvector_complex, occ_numbers)
        call de_kweight_occs('scf_solver:meta', occ_numbers)
      end if
      if (n_dft_methods .gt. 0) then
          write(info_str,'(2X,A)') &
              "----------------------------------------------------------- "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)') &
              "Start evaluation of XC contributions for various DFT methods"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)') &
              "----------------------------------------------------------- "
          call localorb_info(info_str,use_unit,'(A)')
          do i_dft = 1, n_dft_methods
              if (dft_print_out(i_dft) .eq. 6) then
                  dft_flag = 'PBE'
              else if (dft_print_out(i_dft) .eq. 9) then
                  dft_flag = 'BLYP'
              else if (dft_print_out(i_dft) .eq. 16) then
                  dft_flag = 'SVWN'
              else if (dft_print_out(i_dft) .eq. 22) then
                  dft_flag = 'xPBE'
              else if (dft_print_out(i_dft) .eq. 51) then
                  dft_flag = 'TPSS'
              else if (dft_print_out(i_dft) .eq. 59) then
                  dft_flag = 'SCAN'
              else
                  dft_flag = 'UNKNOWN'
              endif
              write(info_str,'(2X,A,A8)') "XC contributuion for ", dft_flag
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A,A16,A16)') &
                 "|              ", 'current', 'original'
              call localorb_info(info_str,use_unit,'(A)')
              flag_xc_old         = flag_xc
              flag_xc             = dft_print_out(i_dft)
              write(info_str,'(2X,A,I16,I16)') &
                 "| flag_xc      ", flag_xc, flag_xc_old
              call localorb_info(info_str,use_unit,'(A)')
              hybrid_coeff_old    = hybrid_coeff
              hybrid_coeff        = 0.0d-1
              write(info_str,'(2X,A,F16.8,F16.8)') &
                 "| hybrid_coeff ", hybrid_coeff, hybrid_coeff_old
              call localorb_info(info_str,use_unit,'(A)')
              ! to evaluate meta-gga components
              if(dft_flag=='TPSS'.or.dft_flag=='SCAN') then
                use_meta_gga_old = use_meta_gga
                use_meta_gga = .true.
              end if
              call integrate_xc_energy ( &
                                         partition_tab,&
                                         rho,rho_gradient,&
                                         kinetic_density,&
                                         xc_energy,&
                                         x_energy,&
                                         c_energy,&
                                         x_energy_lda,&
                                         c_energy_lda&
                                         )
              if (n_spin.eq.2) then
                write(info_str,'(2X,A,A8,A,F20.9,A,F20.9,A)') &
                   "| X ",dft_flag,"               : ", x_energy(1)+x_energy(2), " Ha     ", &
                                         (x_energy(1)+x_energy(2))*hartree, " eV"
              else
                write(info_str,'(2X,A,A8,A,F20.9,A,F20.9,A)') &
                   "| X ",dft_flag,"               : ", x_energy(1), &
                   " Ha     ", (x_energy(1))*hartree, " eV"
              endif
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A,A8,A,F20.9,A,F20.9,A)') &
                   "| C ",dft_flag,"               : ", c_energy, " Ha     ", c_energy*hartree, " eV"
              call localorb_info(info_str,use_unit,'(A)')
              ! now recover the original DFT setting
              flag_xc = flag_xc_old
              hybrid_coeff = hybrid_coeff_old
              if(dft_flag=='TPSS'.or.dft_flag=='SCAN') then
                  use_meta_gga = use_meta_gga_old
              end if
          enddo
          write(info_str,'(2X,A)') &
              "----------------------------------------------------------- "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)') &
              "End evaluation of XC contributions for various DFT methods"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)') &
              "----------------------------------------------------------- "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)') " "
          call localorb_info(info_str,use_unit,'(A)')
      endif

      if(.not.use_embedding_pp) then
      if ((flag_xc.ne.3) .and. (flag_xc.ne.20)) then
          if(myid.eq.0 .and. OL_norm >= output_priority) then
              write(use_unit,'(2X,A)') "------------------------------------ "
              write(use_unit,'(2X,A)') "Start decomposition of the XC Energy "
              write(use_unit,'(2X,A)') "------------------------------------ "
          endif

          call integrate_xc_energy ( &
                                     partition_tab,&
                                     rho,rho_gradient,&
                                     kinetic_density,&
                                     xc_energy,&
                                     x_energy,&
                                     c_energy,&
                                     x_energy_lda,&
                                     c_energy_lda&
                                     )

        if(myid.eq.0 .and. OL_norm >= output_priority) then
          ! The HF contributions report 0 Ha without this update:
          if(use_periodic_hf) en_hf = -exx_ene * hybrid_coeff
          write(use_unit,'(2X,A)') "X and C from original XC functional choice "
          write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
           "Hartree-Fock Energy : ", en_hf, " Ha     ", (en_hf)*hartree, " eV"
          if (n_spin.eq.2) then
            write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
               "X Energy            : ", x_energy(1)+x_energy(2)+en_hf, " Ha     ", &
                                     (x_energy(1)+x_energy(2)+en_hf)*hartree, " eV"
          else
            write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
               "X Energy            : ", x_energy(1)+en_hf, &
               " Ha     ", (x_energy(1)+en_hf)*hartree, " eV"
          endif ! n_spin
          write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
           "C Energy            : ", c_energy, " Ha     ", c_energy*hartree, " eV"
          write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
           "Total XC Energy     : ", xc_energy+en_hf, " Ha     ", (xc_energy+en_hf)*hartree, " eV"

          ! AJL, Apr 2017. Tried to simplify this but still too many if statements!
          ! No output if restarting_scf (sc_init_iter reached)
          if (use_meta_gga_post .and. .not. restarting_scf) then
            write(use_unit,'(2X,A)') "------------------------------------ "
            write(use_unit,'(2X,A)') "Post-processing meta-GGA X and C from self-consistent density "
            write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
             "Hartree-Fock Energy : ", en_hf, " Ha     ", (en_hf)*hartree, " eV"
            write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
             "X Energy meta-GGA   : ", x_post_energy+en_hf, &
                 " Ha     ", (x_post_energy+en_hf)*hartree, " eV"
            write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
             "C Energy meta-GGA   : ", c_post_energy, " Ha     ", c_post_energy*hartree, " eV"
            write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
             "Total XC Energy     : ", en_post_xc+en_hf, " Ha     ", (en_post_xc+en_hf)*hartree, " eV"
          endif ! use_meta-gga

          write(use_unit,'(2X,A)') "------------------------------------ "
          write(use_unit,'(2X,A)') "LDA X and C from self-consistent density "
          if (n_spin.eq.2) then
            write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
               "X Energy LDA        : ", x_energy_lda(1)+x_energy_lda(2), " Ha     ", &
                                     (x_energy_lda(1)+x_energy_lda(2))*hartree, " eV"
          else
            write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
              "X Energy LDA        : ", x_energy_lda(1), " Ha     ", x_energy_lda(1)*hartree, " eV"
          endif
          write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
           "C Energy LDA        : ", c_energy_lda, " Ha     ", c_energy_lda*hartree, " eV"
          write(use_unit,'(A)') "  ------------------------------------ "
          write(use_unit,'(A)') "  End decomposition of the XC Energy "
          write(use_unit,'(A)') "  ------------------------------------ "
          write(use_unit,'(A)') ""
          write(use_unit,'(A)') "------------------------------------------------------------"
        endif
     endif
     endif
      ! ------------------------------------------------------------------------------------------------------------------------------------
      ! CC: Energy density calculations (incl. post-processing) and output
      if (flag_energy_density) then

        call localorb_info ( '', use_unit )
        write(info_str,'(2X,A)') "Computation and output of the energy density:"
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,2A)') "----------------------------------------", &
           "-----------------------------------------------------------------"
        call localorb_info ( info_str, use_unit )

        ! CC: Calculate eigenvalue related energy density as required
        !     in the Harris-Foulkes-like energy density
        if ( flag_harris_foulkes_energy_density ) then

          write(info_str,'(2X,2A)') "= Calculating the KS-state specific ", &
             "energy density as required for the Harris-Foulkes energy density:"
          call localorb_info ( info_str, use_unit )
          write(info_str,'(2X,3A)') "--------------------------------------", &
             "-------------------------------------------------------------", &
             "------"
          call localorb_info ( info_str, use_unit )

          ! CC: Allocate the necessary arrays
          if (.not.allocated(ed_eev_times_kweights)) allocate(ed_eev_times_kweights(n_states,n_spin,n_k_points))
          if (.not.allocated(ed_eev_energy_density)) allocate(ed_eev_energy_density(n_full_points,n_spin))


          ! CC: UGLY! We reweight the occupation number of each orbital with its respective
          !           eigenvalue and then construct the associated density.
          !           N.B. To ensure the respective occ_numbers to be positive (as required in the
          !           subroutines used to construct the density), we weigh with -1 * KS_eigenvalue
          ! CC: Update: We now shift the complete spectrum by ed_E_shift so that the lowest eigenvalue
          !             becomes positive - ed_construct_pseudo_occupation returns ed_eev_shift. After calculating
          !             the eev ed, we the shift that pointwise by ed_eev_shift*rho
          call kweight_occs(   'scf_solver:eev_density', occ_numbers)
          call ed_construct_pseudo_occupation( occ_numbers, KS_eigenvalue, ed_eev_times_kweights, ed_eev_shift )
          call de_kweight_occs('scf_solver:eev_density', occ_numbers)

          if(.not. use_density_matrix)then

            ! do not calculate forces
            pulay_forces_backup = pulay_forces_on
            gga_forces_backup = gga_forces_on
            meta_gga_forces_backup = meta_gga_forces_on
            nlcc_forces_backup = nlcc_forces_on
            Gnonmf_forces_backup = Gnonmf_forces_on
            gga_forces_on = .false.
            meta_gga_forces_on = .false.
            pulay_forces_on = .false.
            nlcc_forces_on = .false.
            Gnonmf_forces_on = .false.

            ! Avoid calculation of density gradients
            use_density_gradient_save = use_density_gradient
            use_density_gradient      = .false.
            use_meta_gga_save         = use_meta_gga
            use_meta_gga              = .false.

            ! VERY UGLY: Dummy 0 array
            if(.not.allocated(ed_null1))      allocate(ed_null1(n_spin,n_full_points))
            if(.not.allocated(ed_null2))      allocate(ed_null2(3,n_spin,n_full_points))
            if(.not.allocated(ed_null3))      allocate(ed_null3(3,n_spin,n_full_points))
            if(.not.allocated(ed_null4))      allocate(ed_null4(n_spin,n_full_points))
            if(.not.allocated(ed_null5))      allocate(ed_null5(n_spin,n_full_points))
            ed_null1(:,:)         = 0.0d0
            ed_null2(:,:,:)       = 0.0d0
            ed_null3(:,:,:)       = 0.0d0
            ed_null4(:,:)         = 0.0d0
            ed_null5(:,:)         = 0.0d0

            call update_density_and_forces_orbital &
                  ( KS_eigenvector, KS_eigenvector_complex, &
                  KS_eigenvalue, ed_eev_times_kweights, &
                  partition_tab, hartree_partition_tab, ed_null1(:,:), ed_null2(:,:,:) , &
                  ed_null4(:,:), hartree_potential, l_shell_max, &
                  ed_eev_energy_density, ed_null3(:,:,:), ed_null5(:,:), rho_change, &
                  hellman_feynman_forces, pulay_forces, &
                  gga_forces, nlcc_forces, pulay_forces_on, gga_forces_on, &
                  nlcc_forces_on, meta_gga_forces_on)

            if(allocated(ed_null1))    deallocate(ed_null1)
            if(allocated(ed_null2))    deallocate(ed_null2)
            if(allocated(ed_null3))    deallocate(ed_null3)
            if(allocated(ed_null4))    deallocate(ed_null4)
            if(allocated(ed_null5))    deallocate(ed_null5)

            ! Restore use_density_gradient
            gga_forces_on = gga_forces_backup
            meta_gga_forces_on = meta_gga_forces_backup
            pulay_forces_on = pulay_forces_backup
!            use_density_gradient_save = use_density_gradient !? This doesn't seem right - should this be the other way round as we  are restoring the saved flag?
!            use_density_gradient      = .false. ! As above, this seems the wrong way round?
            use_density_gradient = use_density_gradient_save
            ! AJL, March2018
            use_meta_gga = use_meta_gga_save
            nlcc_forces_on = nlcc_forces_backup
            Gnonmf_forces_on = Gnonmf_forces_backup

          else  ! use_density_matrix = .true.

               !if(squeeze_memory) then

               !  write(info_str,'(6X,A)') " *** squeeze_memory ==  .true. not yet tested in computation of ed_eev_energy_density ..."
               !  call localorb_info ( info_str,use_unit,'(A)' )
               !  !!! call aims_stop()

               !  !!! !! CODE SHOULD WORK - but not yet tested
               !  !!! !! CC: So far, squeeze_memory must not be set for the computation of ed_eev_energy_density
               !  !!! !!     This has no formal reason and might get changed in future (FIXME)
               !  !!! !!     For future reference, I hence keep the respective code here.
               !  !!!
               !  !!! !! CC: Begin code that has to be adapted
               !  !!! ! call update_density_densmat &
               !  !!! !       ( KS_eigenvector, KS_eigenvector_complex, &
               !  !!! !       occ_numbers, &
               !  !!! !       partition_tab, hartree_partition_tab, rho, rho_gradient, &
               !  !!! !       l_shell_max, delta_rho_KS, delta_rho_gradient, rho_change,&
               !  !!! !       ham_ovlp_work  &
               !  !!! !       )
               !  !!! !! CC: End code that has to be adapted
               !
               !  ! Avoid calculation of density gradients
               !  use_density_gradient_save = use_density_gradient
               !  use_density_gradient      = .false.

               !  ! VERY UGLY: Dummy 0 array
               !  if(.not.allocated(ed_null1))      allocate(ed_null1(n_spin,n_full_points))
               !  if(.not.allocated(ed_null2))      allocate(ed_null2(3,n_spin,n_full_points))
               !  if(.not.allocated(ed_null3))      allocate(ed_null3(3,n_spin,n_full_points))
               !  ed_null1(:,:)         = 0.0d0
               !  ed_null2(:,:,:)       = 0.0d0
               !  ed_null3(:,:,:)       = 0.0d0

               !  call update_density_densmat &
               !        ( KS_eigenvector, KS_eigenvector_complex, &
               !        ed_eev_times_kweights, &
               !        partition_tab,  hartree_partition_tab, ed_null1(:,:), ed_null2(:,:,:) , &
               !        l_shell_max, ed_eev_energy_density, ed_null3(:,:,:), rho_change,&
               !        ham_ovlp_work )

               !  if(allocated(ed_null1))    deallocate(ed_null1)
               !  if(allocated(ed_null2))    deallocate(ed_null2)
               !  if(allocated(ed_null3))    deallocate(ed_null3)

               !  ! Restore use_density_gradient
               !  use_density_gradient_save = use_density_gradient
               !  use_density_gradient      = .false.

               !else

                 ! Avoid calculation of density gradients
                 use_density_gradient_save = use_density_gradient
                 use_density_gradient      = .false.
                 ! AJL March 2018
                 use_meta_gga_save = use_meta_gga
                 use_meta_gga = .false.

                 ! VERY UGLY: Dummy 0 array
                 if(.not.allocated(ed_null1))      allocate(ed_null1(n_spin,n_full_points))
                 if(.not.allocated(ed_null2))      allocate(ed_null2(3,n_spin,n_full_points))
                 if(.not.allocated(ed_null3))      allocate(ed_null3(3,n_spin,n_full_points))
                 if(.not.allocated(ed_null4))      allocate(ed_null4(n_spin,n_full_points))
                 if(.not.allocated(ed_null5))      allocate(ed_null5(n_spin,n_full_points))

                 ed_null1(:,:)         = 0.0d0
                 ed_null2(:,:,:)       = 0.0d0
                 ed_null3(:,:,:)       = 0.0d0
                 ed_null4(:,:)         = 0.0d0
                 ed_null5(:,:)         = 0.0d0

                 call update_density_densmat &
                       ( KS_eigenvector, KS_eigenvector_complex, &
                       ed_eev_times_kweights, &
                       partition_tab,  hartree_partition_tab, ed_null1(:,:), ed_null2(:,:,:) , &
                       ed_null4(:,:), l_shell_max, ed_eev_energy_density, ed_null3(:,:,:), &
                       ed_null5(:,:), rho_change, hamiltonian(1,1) )

                 if(allocated(ed_null1))    deallocate(ed_null1)
                 if(allocated(ed_null2))    deallocate(ed_null2)
                 if(allocated(ed_null3))    deallocate(ed_null3)
                 if(allocated(ed_null4))    deallocate(ed_null4)
                 if(allocated(ed_null5))    deallocate(ed_null5)

                 ! Restore use_density_gradient
!                 use_density_gradient_save = use_density_gradient ! As above, this doesn't restore? Is this a bug?
!                 use_density_gradient      = .false.
                 use_density_gradient = use_density_gradient_save
                 ! AJL March 2018
                 use_meta_gga = use_meta_gga_save

               !endif

          end if

          ! Re-shift eev energy density
          if (ed_eev_shift .gt. 0.0d0) then
            write(info_str,'(2X,A,F20.8,A)') &
               "= Shifting eev-energy density by ", &
               -1.0d0 * ed_eev_shift * hartree, &
               " eV to correct for the previous shifting of the eigenvalues:"
            call localorb_info ( info_str, use_unit )
            call ed_eev_shift_back(ed_eev_energy_density,ed_eev_shift,rho)
          end if

          ! Deallocate
          if (allocated(ed_eev_times_kweights)) deallocate(ed_eev_times_kweights)
          write(info_str,'(2X,3A)') "--------------------------------------", &
             "-------------------------------------------------------------", &
             "------"
          call localorb_info ( info_str, use_unit )

        end if

        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A)') "= Postprocessing energy density:"
        call localorb_info ( info_str, use_unit )
        call ed_write_chetty_martin_values(partition_tab)
        call localorb_info ( '', use_unit )
      end if
      ! End Energy Density related procdures------------------------------------------------------------------------------------------------

      ! Write eigenvectors for the purpose of a wave function extrapolation
      ! over several MD steps.
      if (use_molecular_dynamics .and. use_wf_extrapolation) then
        call wf_save_eigenvectors(occ_numbers, KS_eigenvector, KS_eigenvector_complex)
      end if

      if (restart_from_density) then
         open(unit = 27, file="CMH-Density.out",form="unformatted",status="replace")
         write(27) rho
         close(27)

        ! !open(unit = 28, file="CMH-Density-values.out",form="formatted",status="replace")
         !write(28,*) rho
         !close(28)
         open(unit = 29, file="CMH-Density-grad.out",form="unformatted",status="replace")
         write(29) rho_gradient
         close(29)
         !open(unit = 29, file="CMH-Hartree_Potential.out",form="formatted",status="new")
         !write(29,*) hartree_potential
         !close(29)
         !open(unit = 30, file="CMH-partition_tab.out", form="formatted",status="new")
         !write(30,*) partition_tab
         !close(30)
       endif


      if(flag_output_effective_potential)then

         call  output_effective_potential( hartree_potential, rho, rho_gradient, kinetic_density, &
                                           free_hartree_superpos, &
                                           free_rho_superpos, free_rho_gradient_superpos, partition_tab   )

      end if

      if(flag_neutral_excitation) then
         call optical_response()
      end if

!         if (force_occupation_basis.or.force_occupation_projector) then
!           do i_basis = 1, n_basis, 1
!             i_fn = basis_fn(i_basis)
!             write(info_str,*) i_basis,basis_atom(i_basis),basisfn_type(i_fn),basisfn_n(i_fn),basis_l(i_basis),&
!               basis_m(i_basis)
!             call localorb_info(info_str, use_unit)
!           end do
!           do i_force_occ = 1, n_force_occ, 1
!             do i_basis = 1, n_basis, 1
!               write(info_str,*) "composition",i_force_occ,force_occ_state_periodic(i_force_occ,1),i_basis,KS_eigenvector_complex(i_basis,force_occ_state_periodic(i_force_occ,1),&
!                             force_occ_spin(i_force_occ),1)
!               call localorb_info(info_str, use_unit)
!             end do
!           end do
!         end if

      if (use_molecular_dynamics .or. use_geo_relaxation) then
        call get_timestamps (rtime, clock_rtime)
        if(converged) then
          call output_timeheader(deffmt, "Relaxation / MD: End force evaluation.", OL_high)
          call output_times(deffmt, "Time for this force evaluation", rtime-time_sc_tot, clock_rtime-clock_time_sc_tot, OL_high)
          !  here may be expanded to sub timings
        end if
      end if

      if(output_level .eq. 'MD_light') output_priority = 1

      if(transport_lead_calculation)then

         call transport_find_local_min_occupated_state(occ_numbers, KS_eigenvalue,KS_eigenvector, &
              KS_eigenvector_complex, overlap_matrix )

!         call transport_average_potential( hartree_potential, rho, rho_gradient, partition_tab )
!         call transport_evaluate_reference_energy_level (delta_v_hartree_part_at_zero)



         if(myid==0) call transport_write_boundary_conditions( overlap_matrix, hamiltonian,chemical_potential )

      elseif(transport_calculation)then

         call transport_find_local_min_occupated_state(occ_numbers, KS_eigenvalue,KS_eigenvector, &
              KS_eigenvector_complex, overlap_matrix )



 !        call transport_average_potential( hartree_potential, rho, rho_gradient, partition_tab )
 !        call transport_evaluate_reference_energy_level (delta_v_hartree_part_at_zero)

         call cleanup_multipole_moments
         call deallocate_grid_storage
         call cleanup_pulay
         call cleanup_broyden

         if(use_scalapack .or. use_elpa)then

            call construct_ham_and_ovl_transport_scalapack( hamiltonian, overlap_matrix )
!            call deallocate_physics()
            call read_lead_information_scalapack( chemical_potential )
         else

            call read_lead_information( overlap_matrix, hamiltonian,chemical_potential )
            ! if(myid==0) call read_lead_information( overlap_matrix, hamiltonian,chemical_potential )
         end if
      end if

      if (fo_finalise) then
         call fodft_select_hab(hamiltonian, overlap_matrix)
      end if

      ! Final warning - if we were not converged, this must be brought to the user's attention
     ! if (.not.below_it_limit) then
      if (.not.converged) then
        write(info_str,'(A)') ''
        call localorb_info(info_str, use_unit)
        write(info_str,'(1X,A)') &
        "***********************************************************"
        call localorb_info(info_str, use_unit)
        write(info_str,'(1X,A)') &
        "* "
        call localorb_info(info_str, use_unit)
        write(info_str,'(1X,A)') &
        "* WARNING! SELF-CONSISTENCY CYCLE DID NOT CONVERGE"
        call localorb_info(info_str, use_unit)
        write(info_str,'(1X,A)') &
        "* USING YOUR PRESELECTED ACCURACY CRITERIA."
        call localorb_info(info_str, use_unit)
        write(info_str,'(1X,A)') &
        "* DO NOT RELY ON ANY FINAL DATA WITHOUT FURTHER CHECKS."
        call localorb_info(info_str, use_unit)
        write(info_str,'(1X,A)') &
        "* "
        call localorb_info(info_str, use_unit)
        write(info_str,'(1X,A)') &
        "***********************************************************"
        call localorb_info(info_str, use_unit)

        if (postprocess_anyway == PP_ANYWAY_NOTHING) then
           write(info_str,'(1X,A)') &
           "* NO POSTPROCESSING WILL BE DONE"
           call localorb_info(info_str, use_unit)
           write(info_str,'(1X,A)') &
          "* IF YOU WANT POSTPROCESSING ANYWAY, SET postprocess_anyway TO .true."
           call localorb_info(info_str, use_unit)
           write(info_str,'(1X,A)') &
           "***********************************************************"
           call localorb_info(info_str, use_unit)

           if (out_aims_json_log) then
             call write_geometry_to_json_log("final")
             call write_final_output_to_json_log(.false.)
             call close_json_log_file()
           end if
           call final_deallocations ( )
           call final_timings_and_memory(.false.)
           call aims_stop_coll('SCF cycle not converged.', func)
        end if

        write(info_str,'(A)') ''
        call localorb_info(info_str, use_unit)
      end if

      if(allocated(occ_numbers_save)) deallocate(occ_numbers_save)
!       if (allocated(scf_gradient_basis_wave)) deallocate(scf_gradient_basis_wave)

      if(use_vdw_method.or.use_vdw_post.or.use_nlcorr_post)then
         call RemoveKernelFromMemory  !SAG.
      endif

      ! Reset the index arrays for local_index communication to the integration
      ! batch permutation, since some post-processed functionality assumes that
      ! these arrays will be allocated thusly
      if (use_local_index .and. use_load_balancing) then
        call init_comm_full_local_matrix &
             ( batch_perm(n_bp_integ)%n_basis_local, &
               batch_perm(n_bp_integ)%i_basis_local )
      end if

      ! Allocatable arrays that are tracked
      if (allocated(my_hamiltonian_matrix)) call aims_deallocate( my_hamiltonian_matrix, "my_hamiltonian_matrix" )
      if (allocated(my_density_matrix))     call aims_deallocate( my_density_matrix,         "my_density_matrix" )

    end subroutine scf_solver
!******
