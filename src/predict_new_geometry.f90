! VB 2006
!
! Subroutine predict_new_geometry is a high-level wrapper around
! (possibly) different relaxation algorithms that move the geometry
! of the system by one step.
!
! Historically, the "geometry" module is used by the DFT solver to house
! only the quantities needed for a single-shot run. A separate module
! "relaxation" will keep track of the relaxation history of the system,
! including all parameters relevant to the various relaxation algorithms.
!
! FIXME: (1) Numerical stress tensor:
!            For every relaxation step that uses the stress tensor one scf cycle could
!            be saved, because the zero-point energy is calculated twice (9% speedup)
!            This also applies for calculations without relaxation. Here the first
!            distorted geometry for the numerical stress could be calculated already
!            in initialize_scf.

  subroutine predict_new_geometry ( converged_scf, final_geo_found, valid_geo, cont_bsse )

    use localorb_io
    use dimensions
    use constants
    use runtime_choices
    use timing
    use physics
    use geometry
    use species_data
    use relaxation
    use molecular_dynamics
    use pi_molecular_dynamics       ! XZL: added for PIMD
    use thermodynamic_integration
    use numerical_stress
    use analytical_stress
    use external_pressure_on_system
    use mpi_utilities
    use potconstrain
    use bravais , only: get_cell_volume
    use symmetry_constrained_relaxation
    ! TARP: Done to fix collisions, or if the volume decreases significantly
    use trust_radius_method
    use geometry
    implicit none

    ! imported variables

    logical :: converged_scf
    logical :: final_geo_found
    logical :: valid_geo
    logical :: cont_bsse  ! for atom_bsse

    ! Local variables

    logical :: converged_scf_bead
    logical :: walltime_left_bead

    real*8 :: free_energy, pot_energy

    real*8, dimension(3,n_atoms) :: converted_forces
    real*8, dimension(3,n_atoms) :: converted_coords

    real*8, dimension(3,3)       :: lattice_vector_orig
    real*8, dimension(3,3)       :: lattice_vector_tmp
    real*8, dimension(3,3)       :: lv_tmp2
    real*8                       :: lattice_sum

    real*8, dimension(3,n_atoms) :: periodic_unit_cell_translations_tmp

    real*8, dimension(3,3)       :: stress_tensor_for_relaxation

    character*80                 :: info_str


    !for atom_bsse: might need cleaning up later
    integer  current_bsse_atom
    !   bsse time keeping
    real*8 :: time_bsse_per_atom
    real*8  temp_time_bsse
    real*8  temp_clock_time_bsse

    logical :: abort_requested = .false.

    ! counters

    integer :: i_atom, i_bead, j_bead   ! XZL: i_bead added for PIMD
    integer :: i_coord, i_periodic
    integer :: i,j,g,k,l

    ! check variable

    integer :: MD_stop = 1

    !! CC & MOL: symmetry-constrained relaxation
    integer :: n_mobile_lv_components_tmp, n_mobile_atom_coords_tmp, n_mobile_coords_tmp
    real*8  :: volume, SCR_scaling
    ! temporary arrays to store transformed lattice/coords/forces
    ! real*8, dimension(3,n_periodic) :: SCR_lv, SCR_lv_tmp
    real*8, dimension(3,n_periodic) :: SCR_lv, SCR_lv_tmp
    real*8, dimension(3,n_atoms)    :: SCR_fcoords, SCR_coords_tmp, SCR_fcoords_tmp
    real*8, dimension(3,n_atoms)    :: SCR_coords_sr, SCR_forces_sr
    real*8, dimension(3,n_periodic) :: SCR_lv_sr, SCR_forces_lv_sr
    real*8, dimension(3,n_atoms)    :: SCR_forces_tmp
    real*8, dimension(3,n_periodic)    :: SCR_forces_lv_tmp
    ! 1d arrays
    real*8, dimension(3*n_atoms)           :: SCR_fcoords_1d, SCR_forces_1d
    real*8, dimension(3*n_periodic)        :: SCR_lv_1d, SCR_forces_lv_1d

    ! TARP: Collision check for symmetry-constrained relaxation
    real*8, dimension(SCR_n_params) :: X, F
    real*8, dimension(SCR_n_params) :: DX_step
    ! real*8, dimension(3) :: SCR_scaling
    real*8 :: this_dist, new_volume, volume_change, i_trust_rad
    ! allocatable arrays like transformation matrices and symm-reduced arrays
    ! are allocated and created in the relaxation module and used here

   ! Begin work

    call localorb_info('',use_unit)
    call localorb_info("------------------------------------------------------------",use_unit,'(A)')
    if (use_geo_relaxation .and. (num_stress_finished .or. (.not. use_numerical_stress))) then
      call localorb_info("Geometry optimization: Attempting to predict improved coordinates.",use_unit,'(2X,A)' )
    else if (use_molecular_dynamics .and. (num_stress_finished .or. (.not. use_numerical_stress))) then
      call localorb_info("Molecular dynamics: Attempting to update all nuclear coordinates.",use_unit,'(2X,A)' )
    endif
    call localorb_info('',use_unit)

    ! Initialize check for validity of geometry
    valid_geo = .false.

    ! initialize continuation flag for loop over various atoms
    cont_bsse = .false.

    ! Memorize lattice vectors for one definitive check of whether
    ! or not unit cell shape changed, at end of subroutine.
    ! lattice_vector_orig should not be altered at any point
    ! in this subroutine.
    if (n_periodic > 0) then
       lattice_vector_orig = lattice_vector
    endif

    if (.not.converged_scf) then
      ! no good forces/energies available from last iteration - aborting relaxation

      valid_geo = .false.

      call localorb_info("******************************************************************",use_unit,'(1X,A)' )
      call localorb_info("* Attention: preceding scf. cycle did not converge.",use_unit,'(1X,A)' )
      call localorb_info("* Cannot provide updated geometry/numerical stress tensor ",use_unit,'(1X,A)' )
      call localorb_info("* based on non-converged forces/energies - ",use_unit,'(1X,A)' )
      call localorb_info("* aborting relaxation.",use_unit,'(1X,A)' )
      call localorb_info("******************************************************************",use_unit,'(1X,A)' )
      call localorb_info('',use_unit)

    elseif (use_numerical_stress .and. (.not. num_stress_finished)) then
       valid_geo = .true.

       if (counter_numerical_stress == 0) then

         ! We freshly enter numerical stress
         call localorb_info('',use_unit)
         call localorb_info("------------------------------------------------------------",use_unit,'(A)')
         call localorb_info("Calculation of numerical stress starts. ",use_unit,'(2X,A)' )
         call localorb_info("------------------------------------------------------------",use_unit,'(A)')
         call localorb_info('',use_unit)

         ! check whether we are in 3 dimensions
         if (n_periodic .ne. 3 ) then
            call localorb_info("* Attention: Calculation of numerical stress can be only done in 3D",use_unit,'(2X,A)' )
            call localorb_info("* Please adjust your input. FHI-aims will stop here. ",use_unit,'(2X,A)' )
            call localorb_info('',use_unit)
            stop
         endif

         ! reduce output level and remember skip_scf flag for time saving options later
         original_output_level  = output_level
!         output_level           = 'MD_light'
         original_skip_scf_flag = skip_scf

         ! Nothing happend so far to the geometry. Store it!
         call allocate_numerical_stress ()
         call remember_original_geometry ()

       endif

       if (counter_numerical_stress <  12) then

         ! If the preceeding scf calculation contributed to the stress,
         ! then save the free energy
         if (counter_numerical_stress > 0 .and. .not. skip_scf) then

            ! Temporary change counter in order to have the right index picture
            ! when storing the free energies
            free_energy              = total_energy + 2.d0*entropy_correction
            counter_numerical_stress = counter_numerical_stress-1
            call energy_contribution_to_numerical_stress (free_energy, counter_numerical_stress)
            ! Restore the current counter
            counter_numerical_stress = counter_numerical_stress+1

         endif
         skip_scf = original_skip_scf_flag

         ! Get distorted structure
         call apply_strain_transformation (counter_numerical_stress)

         call get_indices_numerical_stress (counter_numerical_stress, i ,j, g)

         call localorb_info&
         ("Components of the stress tensor (for mathematical background see comments in numerical_stress.f90).",use_unit,'(2X,A)' )
         call localorb_info&
         ("For a distortion R -> (1+e)R with 3x3 strain-matrix e compute derivative of free energy by  ",use_unit,'(2X,A)' )
         write(info_str,'(2X,A,I2.1,A,I2.1,A,F11.8)') &
                         "strain component: ", i, ",", j, "; displacement: ", (-1)**g*delta_numerical_stress
         call localorb_info(info_str,use_unit,'(A)')

         ! set skip_scf flag using the strain transformations
         if (numerical_stress_save_scf .and. use_relaxation_constraints) then
            lattice_vector_tmp(:,:) = 1d0
            do k = 1, 3
               do l = 1, 3
                  if (constrain_lv_components(k,l)) lattice_vector_tmp(k,l) = 0d0
               end do
            end do

            lattice_sum = 0d0
            do  i_periodic = 1, 3
               call apply_strain_transformation_on_vector(i, j, g, lattice_vector_tmp(:,i_periodic),lv_tmp2(:,i_periodic))
               lv_tmp2(:,i_periodic) = (lv_tmp2(:,i_periodic)-lattice_vector_tmp(:,i_periodic))*lattice_vector_tmp(:,i_periodic)
               lattice_sum = lattice_sum + sum(lv_tmp2(:,i_periodic))
            end do

            if (abs(lattice_sum) < 1d-8) then
               call localorb_info&
               ("According to the relaxation constraints, this displacement does not need to be calculated.",use_unit,'(2X,A)')
               call localorb_info("Skipping the next scf cycle.",use_unit,'(2X,A)')
               skip_scf = .true.
            endif
         endif

         counter_numerical_stress =  counter_numerical_stress+1

      else if (counter_numerical_stress == 12) then

         skip_scf = original_skip_scf_flag

         ! There is one component from the previous geometry
         ! missing to the numerical stress (counter_numerical_stress 11)

         free_energy              = total_energy + 2.d0*entropy_correction
         ! Temporary change counter in order to have the right index picture
         ! when storing the free energies
         counter_numerical_stress = counter_numerical_stress-1
         call energy_contribution_to_numerical_stress (free_energy, counter_numerical_stress)
         ! Restore the current counter
         counter_numerical_stress = counter_numerical_stress+1

         ! Finished calculation of all components of numerical stress - put it together
         call update_numerical_stress_tensor ()

         ! RESTORE VARIABLES before leaving the numerical stress

         ! Reset counter for the calculation of the stress tensor
         ! for the next geometry relaxation step
         counter_numerical_stress = 0
         num_stress_finished = .true.
         output_level =  original_output_level
         ! Set back quantities which changed during stress-calculation.
         ! The rest should be recalculated in intialize_bc_dependent_lists
         call restore_zero_point_geometry ()
         use_forces = original_force_flag
         use_analytical_stress = original_analytical_stress_flag

         ! DERIVED QUANTIES:
         ! Calculate Pressure: strictly speaking this only makes sense in the hydrostatic case
         ! where the stress tensor is a multiple of unit.
         ! Here we calculate the normed trace of the stress tensor, which corresponds
         ! to the arithmetic mean of the diagonal elements (times minus 1).
         numerical_pressure = (  numerical_stress_tensor(1,1) &
                               + numerical_stress_tensor(2,2) &
                               + numerical_stress_tensor(3,3) ) /3d0 *(-1)

         ! OUTPUT of numerical stress and related quantities
         call localorb_info('',use_unit)
         call localorb_info("Calculation of numerical stress completed",use_unit,'(2X,A)' )
         call localorb_info('',use_unit)
         call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
         call localorb_info("|                       Numerical stress tensor                     |",use_unit,'(2X,A)')
         call localorb_info("|                    Cartesian components [eV/A**3]                 |",use_unit,'(2X,A)')
         call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
         call localorb_info("|                x                y                z                |",use_unit,'(2X,A)')
         call localorb_info("|                                                                   |",use_unit,'(2X,A)')

         write(info_str,'(2X,A1,A,3(2X,F15.8),11X,A1)') "|", "  x  ", &
                           (numerical_stress_tensor(1,j)*hartree / (bohr**3),j=1,3,1), "|"
         call localorb_info(info_str,use_unit,'(A)')
         write(info_str,'(2X,A1,A,3(2X,F15.8),11X,A1)') "|", "  y  ", &
                           (numerical_stress_tensor(2,j)*hartree / (bohr**3),j=1,3,1), "|"
         call localorb_info(info_str,use_unit,'(A)')
         write(info_str,'(2X,A1,A,3(2X,F15.8),11X,A1)') "|", "  z  ", &
                           (numerical_stress_tensor(3,j)*hartree / (bohr**3),j=1,3,1), "|"
         call localorb_info(info_str,use_unit,'(A)')

         call localorb_info("|                                                                   |",use_unit,'(2X,A)')
         write(info_str,'(2X,A1,A,2X,F15.8,2X,A,26X,A1)') "|", "  Pressure:", &
                           numerical_pressure *hartree / (bohr**3), " [eV/A**3] " ,"|"
         call localorb_info(info_str,use_unit,'(A)')
         call localorb_info("|                                                                   |",use_unit,'(2X,A)')
         call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')

         call check_pressure(numerical_stress_tensor)

         ! ZERO POINT in the next SCF
         call localorb_info('',use_unit)
         call localorb_info("Calculation of original geometry",use_unit,'(2X,A)' )
         call localorb_info("------------------------------------------------------------",use_unit,'(A)')
         call localorb_info('',use_unit)

         ! Deallocate the arrays which stored the original geometric quantities during stress calculation
         call deallocate_numerical_stress ()

      endif ! counter_numerical_stress

    ! CC & MOL: Standard relaxation algorithm (not symmetry constrained)
    else if ((use_geo_relaxation) .and. (.not.use_symm_const_geo_relaxation)) then

      ! Calculate and output the forces on the lattice - cleaned from the atomic forces
      if (relax_unit_cell .ne. 0 ) then

         ! Set stress tensor which is used for relaxation.
         if (stress_for_relaxation.eq.RELAX_NUMERICAL) then
            stress_tensor_for_relaxation(1:3,1:3) = numerical_stress_tensor(1:3,1:3)
         else if (stress_for_relaxation.eq.RELAX_ANALYTICAL) then
            stress_tensor_for_relaxation(1:3,1:3) = analytical_stress_tensor(1:3,1:3)
            if (AS_components.eq.9) then
               call AS_symmetrize_stress(stress_tensor_for_relaxation)
            endif
         endif

         if (flag_cpu_consistency) &
            call check_cpu_consistency_matrix (stress_tensor_for_relaxation, &
                                     "stress_tensor_initial", 3, 3)

         ! Add external pressure to stress
         if (flag_external_pressure) then
            call add_p_to_stress(stress_tensor_for_relaxation(1:3,1:3), &
                                 external_pressure)
            if (flag_cpu_consistency) &
              call check_cpu_consistency_matrix (stress_tensor_for_relaxation, &
                                     "stress_tensor_pressure", 3, 3)
         endif

         ! the energy gradient does not go into the optimizer. For the time being we
         ! still calculate and output it to see how big the atomic contributions are
         ! -> could/should be removed
         call  get_derivatives_on_lattice (stress_tensor_for_relaxation, energy_gradient_on_lattice)

         if (flag_cpu_consistency) &
            call check_cpu_consistency_matrix(total_forces,"forces",3,n_atoms)

         call  extract_lv_forces &
               (forces_lv, stress_tensor_for_relaxation, total_forces)

         ! Output the energy gradient with respect to the lattice vectors
         call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
         call localorb_info("| Generalized derivatives on lattice vectors [eV/A]                 |",use_unit,'(2X,A)' )
         call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
         do i_periodic = 1, n_periodic, 1
            write(info_str,'(2X,A1,A,3(2X,F15.8),2X,A1)') "|", "lattice_vector", &
                 (energy_gradient_on_lattice(i_coord,i_periodic)*hartree/(bohr) ,i_coord=1,3,1), "|"
            call localorb_info(info_str,use_unit,'(A)')
         enddo
         call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
         call localorb_info("| Forces on lattice vectors cleaned from atomic contributions [eV/A]|",use_unit,'(2X,A)' )
         call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
         do i_periodic = 1, n_periodic, 1
            write(info_str,'(2X,A1,A,3(2X,F15.8),2X,A1)') "|", "lattice_vector", &
                 (forces_lv(i_coord,i_periodic)*hartree/(bohr) ,i_coord=1,3,1), "|"
            call localorb_info(info_str,use_unit,'(A)')
         enddo

         ! If 'shape' unit cell relaxation is required, calculated the force-projections in the direction
         ! of the lattice vectors and output them
         if (relax_unit_cell == 2 ) then

            call  project_lv_forces_4_shape_relaxation (forces_lv)

            call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
            call localorb_info("| Projected forces for shape relaxation  [eV/A]                     |",use_unit,'(2X,A)' )
            call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
            do i_periodic = 1, n_periodic, 1
               write(info_str,'(2X,A1,A,3(2X,F15.8),2X,A1)') "|", "lattice_vector", &
                    (forces_lv(i_coord,i_periodic)*hartree/(bohr) ,i_coord=1,3,1), "|"
               call localorb_info(info_str,use_unit,'(A)')
            enddo

         endif
         call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
      endif!relax_unit_cell .ne. 0


      ! MR: trying to add a constraining potential
      if (use_potential_constrain) then
         call apply_constrain(pot_energy)
         if (myid.eq.0) write(use_unit,*) '| Energy from potential constrain:', pot_energy
      endif

      ! Explicitly clean unitary translations and rotations,
      ! and remove forces on constrained atoms.
      ! This might change the forces.

      call clean_force_components(total_forces, forces_lv)
      call check_geometry_convergence(total_forces, forces_lv, final_geo_found)

      INQUIRE(FILE="abort_opt", EXIST=abort_requested)   ! abort_requested will be TRUE if the file

      if ( (.not.final_geo_found).and.(relaxation_step_number .lt. max_relaxation_steps).and. &
         &  .not.abort_requested ) then

         valid_geo = .true.
         relaxation_step_number = relaxation_step_number + 1

         write(info_str,'(2X,A,I6,A)') "Relaxation step number ", relaxation_step_number, &
              ": Predicting new coordinates."
         call localorb_info(info_str)
         call localorb_info('')

         ! (re-)compute free energy of the present structure
         ! Note that it is the free energy which is consistent with
         ! the forces, not the total energy

         free_energy = total_energy + 2.d0*entropy_correction

         ! add energy from constrain potential
         if (use_potential_constrain) then
            free_energy=free_energy+pot_energy
         endif
         ! Add external pressure * cell_volume to energy
         if (flag_external_pressure) then
            call add_pV_to_energy(free_energy, external_pressure, cell_volume)
         endif

         if (external_forces_on) then
            ! External work from last step
            free_energy = free_energy + external_work
         endif

         ! Pass dummy variable lattice_vector,  forces_lv  also for the cluster case
         ! since want to use the relaxation routine in the periodic/cluster case on equal footing
         if (n_periodic > 0) then

            if (flag_cpu_consistency) &
               call check_cpu_consistency_matrix(lattice_vector, &
                                       "lattice_vector_init", 3, n_periodic)
            lattice_vector_tmp = lattice_vector

            ! FK: Pass periodic_unit_cell_translations to relaxation routine to keep track of it
            !     for the individual geometries. Otherwise, rejecting a relaxation step might lead to
            !     inconsistencies between geometry and periodic shift (atoms jump out of the unit cell).
            !     In the cluster case, we just pass the dummy variable.
            periodic_unit_cell_translations_tmp = periodic_unit_cell_translations
         endif

         if (flag_cpu_consistency) &
            call check_cpu_consistency_matrix (coords,"coords",3,n_atoms)

         call get_next_relaxation_step &
                (free_energy, coords, total_forces,&
                 lattice_vector_tmp, forces_lv, periodic_unit_cell_translations_tmp, &
                 KS_eigenvector,  KS_eigenvector_complex, occ_numbers, valid_geo)

         if (external_forces_on) then

            ! W = -F * delta s
            do i_atom = 1, n_atoms, 1
               do i_coord = 1, 3, 1
                  external_work = external_work &
                  - external_forces(i_coord,i_atom) &
                     * (coords(i_coord,i_atom) - stored_coords(i_coord,i_atom))
               end do
            end do
            write(info_str,'(2X,A,F13.6,A)') "External work evaluated: ", external_work, " H"

         endif


         if (n_periodic > 0) then
            if (any(lattice_vector_tmp /= lattice_vector)) then
               lattice_vector = lattice_vector_tmp
               if (flag_cpu_consistency) &
                  call check_cpu_consistency_matrix(lattice_vector, &
                                          "lattice_vector_out", 3, n_periodic)
               call initialize_bravais_quantities()
            endif

            if (any(periodic_unit_cell_translations_tmp /= periodic_unit_cell_translations)) then
               periodic_unit_cell_translations = periodic_unit_cell_translations_tmp
            endif
         endif

         ! Coordinated have been updated - accounting:
         if (flag_external_pressure) &
            call external_pressure_notice(external_pressure)
         write(info_str,'(2X,A)') "Updated atomic structure: "
         call localorb_info( info_str )

         call output_structure( )
         write(info_str,'(A)') &
              "------------------------------------------------------------"
         call localorb_info(info_str)

         ! Overwrites the lower triangle of the Hessian with values
         ! from the upper triangle. This was originally needed for the
         ! BFGS->TRM transition (see 432ca0ef). Currently the
         ! necessity of this subroutine is unclear.
         call symmetrize_hessian()

         call write_new_geometry_file(coords, lattice_vector)

      else

         if (abort_requested) then
           write(info_str,'(1X,A)') &
                 "* Abort request for optimization found *"
            call localorb_info(info_str)
            write(info_str,'(1X,A)') &
                 "* Aborting optimization."
            call localorb_info(info_str)
            write(info_str,'(A)') ''
            valid_geo = .false.

         else if (.not.final_geo_found) then

            write(info_str,'(1X,A)') &
                 "* Maximum number of relaxation steps exceeded."
            call localorb_info(info_str)
            write(info_str,'(1X,A)') &
                 "* Aborting optimization."
            call localorb_info(info_str)
            write(info_str,'(A)') ''
            valid_geo = .false.

         endif

         write(info_str,'(A)') &
              "------------------------------------------------------------"
         call localorb_info( info_str )
         if (flag_external_pressure) &
            call external_pressure_notice(external_pressure)
         write(info_str,'(2X,A)') "Final atomic structure: "
         call localorb_info( info_str )

         call output_structure( )

         write(info_str,'(A)') &
              "------------------------------------------------------------"
         call localorb_info( info_str )

      endif

      ! VA: Geometry has changed, calculate numerical stress again
      num_stress_finished = .false.

    ! CC & MOL: EXPERIMENTAL: Symmetry constrained relaxation algorithm
    else if ((use_geo_relaxation) .and. (use_symm_const_geo_relaxation)) then
      ! CC & MOL:
      ! This implements a symmetry-constrained relaxation where the coordinates
      ! (both of lattice and atoms) as well as the respective forces on them
      ! are projected (we call it transformed) to a reduced parameter space defined
      ! in the geometry.in by keywords like `symmetry_n_params` etc.
      ! We use the existing arrays of dimension (3,n_atoms) and (3,n_periodic),
      ! fill them with the symmetry-reduced parameters and fill them up with 0.
      ! The relaxation uses the same index vectors (e.g. mobile_atom_coords etc.)
      ! but reduced in the same way to contain only as many entries as the number
      ! of free parameters. (Also the Hessian has this reduced size.)
      !
      ! Not compatible or not yet supported: shape relaxation, external
      ! pressure/forces/potential, potential constrain --> Feel free to implement
      ! it. Compare with the `else if` clause before - the one for the standard
      ! relaxation algorithm - for details on which parts were deleted here.

      ! Calculate and output the forces on the lattice - cleaned from the atomic forces
      if (relax_unit_cell .ne. 0 ) then
        ! Set stress tensor which is used for relaxation.
        if (stress_for_relaxation.eq.RELAX_NUMERICAL) then
           stress_tensor_for_relaxation(1:3,1:3) = numerical_stress_tensor(1:3,1:3)
        else if (stress_for_relaxation.eq.RELAX_ANALYTICAL) then
           stress_tensor_for_relaxation(1:3,1:3) = analytical_stress_tensor(1:3,1:3)
           if (AS_components.eq.9) then
              call AS_symmetrize_stress(stress_tensor_for_relaxation)
           endif
        endif

        if (flag_cpu_consistency) &
           call check_cpu_consistency_matrix (stress_tensor_for_relaxation, &
                                    "stress_tensor_initial", 3, 3)

        ! --------- calculate non-symmetrized forces on lattice and output -------------------
        ! the energy gradient does not go into the optimizer. For the time being we
        ! still calculate and output it to see how big the atomic contributions are
        ! -> could/should be removed
        call  get_derivatives_on_lattice (stress_tensor_for_relaxation, energy_gradient_on_lattice)

        if (flag_cpu_consistency) &
           call check_cpu_consistency_matrix(total_forces,"forces",3,n_atoms)

        call  extract_lv_forces &
              (forces_lv, stress_tensor_for_relaxation, total_forces)
        ! Output the energy gradient with respect to the lattice vectors
        call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
        call localorb_info("| Generalized derivatives on lattice vectors [eV/A]                 |",use_unit,'(2X,A)' )
        call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
        do i_periodic = 1, n_periodic, 1
           write(info_str,'(2X,A1,A,3(2X,F15.8),2X,A1)') "|", "lattice_vector", &
                (energy_gradient_on_lattice(i_coord,i_periodic)*hartree/(bohr) ,i_coord=1,3,1), "|"
           call localorb_info(info_str,use_unit,'(A)')
        enddo
        call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
        call localorb_info("| Forces on lattice vectors cleaned from atomic contributions [eV/A]|",use_unit,'(2X,A)' )
        call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
        do i_periodic = 1, n_periodic, 1
           write(info_str,'(2X,A1,A,3(2X,F15.8),2X,A1)') "|", "lattice_vector", &
                (forces_lv(i_coord,i_periodic)*hartree/(bohr) ,i_coord=1,3,1), "|"
           call localorb_info(info_str,use_unit,'(A)')
        enddo
        call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')

        ! --------- calculate symmetrized forces on lattice and output -------------------
        ! MOL: Symmetrize coords, lattice_vector and forces on atoms
        !      for extraction of forces on lattice vectors
        !      If there are no parameters for the atomic coordinates: all forces
        !      symmetrized are zero --> still do.
      end if
      ! save original in tmp variables
      SCR_forces_lv_tmp = forces_lv
      SCR_forces_tmp = total_forces
      lattice_vector_tmp = lattice_vector
      SCR_coords_tmp = coords
      call cart2frac(lattice_vector, coords, SCR_fcoords)
      SCR_fcoords_tmp = SCR_fcoords
      ! Symmetrize lattice vectors
      if (SCR_n_params_lv > 0) call SCR_symmetrize_lv(lattice_vector)
      SCR_lv = lattice_vector
      if (flag_cpu_consistency) &
         call check_cpu_consistency_matrix (lattice_vector,"lattice_vector_init",3,n_periodic)

      ! Symmetrize atom coordinates
      if (SCR_n_params_coords > 0) then
        call cart2frac(lattice_vector, coords, SCR_fcoords)
        call SCR_symmetrize(SCR_fcoords)
        call frac2cart(lattice_vector, SCR_fcoords, coords)
      endif
      if (flag_cpu_consistency) &
         call check_cpu_consistency_matrix (coords,"coords",3,n_atoms)

      ! Symmetrize forces on atoms
      ! TARP: inverse the transformations since position is contravariant and forces are covariant
      if (SCR_n_params_coords > 0) then
        call SCR_generate_force_Jacobian(lattice_vector)
        call SCR_symmetrize_forces(total_forces)
      else
        total_forces = 0.d0
      end if

      if (flag_cpu_consistency) &
         call check_cpu_consistency_matrix (total_forces,"forces",3,n_atoms)
      call localorb_info("  Symmetrized total atomic cartesian forces [eV/Ang]:", use_unit)
      do i_atom = 1, n_atoms, 1
        write (info_str,'(2X,A,I4,1X,3(E20.10,1X))') &
             "|",i_atom, total_forces(:,i_atom) * hartree / bohr
        call localorb_info(info_str, use_unit)
      enddo
      if (relax_unit_cell .ne. 0) then
        !MOL
        if (flag_cpu_consistency) &
           call check_cpu_consistency_matrix (total_forces,"forces",3,n_atoms)

        ! Extract forces_lv symmetrized
        call  extract_lv_forces &
              (forces_lv, stress_tensor_for_relaxation, total_forces)

        call SCR_symmetrize_forces_lv(forces_lv)

        call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
        call localorb_info("| Symmetrized forces on lattice vectors cleaned from symmetrized     ",use_unit,'(2X,A)' )
        call localorb_info("| atomic contributions [eV/A]|",use_unit,'(2X,A)' )
        call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
        do i_periodic = 1, n_periodic, 1
           write(info_str,'(2X,A1,A,3(2X,F15.8),2X,A1)') "|", "lattice_vector", &
                (forces_lv(i_coord,i_periodic)*hartree/(bohr) ,i_coord=1,3,1), "|"
           call localorb_info(info_str,use_unit,'(A)')
        enddo
        call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
      else ! no unit cell relaxation
        forces_lv = 0.d0
      endif !relax_unit_cell .ne. 0

      if (n_periodic > 0) then
         if (flag_cpu_consistency) &
            call check_cpu_consistency_matrix(lattice_vector, "lattice_vector_init", 3, n_periodic)
         ! FK: Pass periodic_unit_cell_translations to relaxation routine to keep track of it
         !     for the individual geometries. Otherwise, rejecting a relaxation step might lead to
         !     inconsistencies between geometry and periodic shift (atoms jump out of the unit cell).
         !     In the cluster case, we just pass the dummy variable.
         periodic_unit_cell_translations_tmp = periodic_unit_cell_translations
      endif

      if (flag_cpu_consistency) &
          call check_cpu_consistency_matrix (coords,"coords",3,n_atoms)

      ! Explicitly clean unitary translations and rotations, and remove forces on constrained atoms.
      ! This might change the forces.
      call clean_force_components(total_forces, forces_lv)  ! clean cartesian

      ! ------- CHECK CONVERGENCE USING SYMMETRIZED CARTESIAN FORCES -------------
      call check_geometry_convergence(total_forces, forces_lv, final_geo_found)

      ! Transform to parameter space (actually they're symmetrized already...)
      call SCR_transform(SCR_fcoords, SCR_coords_sr)
      call SCR_transform_lv(lattice_vector, SCR_lv_sr)

      call SCR_transform_forces(total_forces, SCR_forces_sr)
      call SCR_transform_forces_lv(forces_lv, SCR_forces_lv_sr)

      if (SCR_n_params_coords > 0) then
      call localorb_info("  Symmetry-reduced parameters:", use_unit)
        do i=1,SCR_n_params_coords
          write (info_str, '(2X,A,I3,1X,E30.15,1X)') '| atomic ', i, &
                SCR_coords_sr( modulo((i-1),3)+1, (i-1)/3 +1 )
          call localorb_info(info_str, use_unit)
        enddo
      endif
      if (SCR_n_params_lv > 0) then
        do i=1,SCR_n_params_lv
          write (info_str, '(2X,A,I3,1X,E30.15,1X)') '| lattice', i, &
                SCR_lv_sr( modulo((i-1),3)+1, (i-1)/3+1 ) * bohr
          call localorb_info(info_str, use_unit)
        enddo
      endif
      call localorb_info("  Symmetry-reduced forces:", use_unit)
      if (SCR_n_params_coords > 0) then
        do i=1,SCR_n_params_coords
          write (info_str, '(2X,A,I3,1X,E30.15,1X)') '| atomic ', i, &
                SCR_forces_sr( modulo((i-1),3)+1, (i-1)/3 +1 ) * hartree
          call localorb_info(info_str, use_unit)
        enddo
      endif
      if (SCR_n_params_lv > 0) then
        do i=1,SCR_n_params_lv
          write (info_str, '(2X,A,I3,1X,E30.15,1X)') '| lattice', i, &
                SCR_forces_lv_sr( modulo((i-1),3)+1, (i-1)/3+1 ) * hartree / bohr
          call localorb_info(info_str, use_unit)
        enddo
      endif

      INQUIRE(FILE="abort_opt", EXIST=abort_requested)   ! abort_requested will be TRUE if the file

      if ( (.not.final_geo_found).and.(relaxation_step_number .lt. max_relaxation_steps).and. &
        &  .not.abort_requested ) then

        valid_geo = .true.
        relaxation_step_number = relaxation_step_number + 1

        write(info_str,'(2X,A,I6,A)') "Relaxation step number ", relaxation_step_number, &
             ": Predicting new coordinates."
        call localorb_info(info_str)
        call localorb_info('')

        ! (re-)compute free energy of the present structure
        ! Note that it is the free energy which is consistent with
        ! the forces, not the total energy

        free_energy = total_energy + 2.d0*entropy_correction

        ! Scale atom parameters with average cell length to make it same order of magnitude
        ! Otherwise TRM will make steps too large because it works for lengths in Cartesian space
        call get_cell_volume(lattice_vector, volume)
        SCR_scaling = volume**(1./3)
        SCR_coords_sr = SCR_coords_sr * SCR_scaling
        SCR_forces_sr = SCR_forces_sr / SCR_scaling
        call coords2X(SCR_coords_sr, SCR_lv_sr, X)
        call coords2X(SCR_forces_sr, SCR_forces_lv_sr, F)

        call get_next_relaxation_step &
               (free_energy, SCR_coords_sr, SCR_forces_sr,&
                SCR_lv_sr, SCR_forces_lv_sr, periodic_unit_cell_translations_tmp, &
                KS_eigenvector,  KS_eigenvector_complex, occ_numbers, valid_geo)

        SCR_coords_sr = SCR_coords_sr / SCR_scaling
        SCR_forces_sr = SCR_forces_sr * SCR_scaling

        ! Back transform to full geometry (no backtrafo of forces needed)
        if (SCR_n_params_coords > 0) then
          call SCR_backtransform(SCR_coords_sr, SCR_fcoords, .false.)
        else
          SCR_fcoords = SCR_fcoords_tmp
        endif
        if (SCR_n_params_lv > 0) then
          call SCR_backtransform_lv(SCR_lv_sr, SCR_lv)
        else
          SCR_lv = lattice_vector_tmp
        endif
        call frac2cart(SCR_lv, SCR_fcoords, coords)

        ! TARP: Check for collisions and Volume SCR_scaling
        call min_atomic_dist(SCR_lv, coords, this_dist)
        i_trust_rad = 1.d0
        do while ( this_dist < 0.69d0/bohr .and. ( i_trust_rad < 20 ) )
            i_trust_rad = i_trust_rad + 1.d0
            write(info_str, "(2X,'| Use smaller step to avoid collisions.')")
            call localorb_info(info_str)
            call trm_trusted_step(F, hessian, trust_radius / i_trust_rad, DX_step)
            call X2coords(SCR_coords_sr, SCR_lv_sr,  X)
            call X2new_coords(SCR_coords_sr, SCR_lv_sr, X, DX_step)
            SCR_coords_sr = SCR_coords_sr / SCR_scaling
            if (SCR_n_params_coords > 0) then
              call SCR_backtransform(SCR_coords_sr, SCR_fcoords, .false.)
            end if
            if (SCR_n_params_lv > 0) then
              call SCR_backtransform_lv(SCR_lv_sr, SCR_lv)
            end if
            call frac2cart(SCR_lv, SCR_fcoords, coords)
            call min_atomic_dist(SCR_lv, coords, this_dist)
        end do

        call get_cell_volume(SCR_lv, new_volume)
        volume_change = ( volume-new_volume ) / volume
        do while ( ( volume_change > 0.3d0 ) .and. ( i_trust_rad < 20 ) )
            i_trust_rad = i_trust_rad + 1.d0
            write(info_str, "(2X,A,F4.1,A)") '** The unit cell volume has decreased by ', 100*volume_change, ' percent.'
            call localorb_info(info_str)
            write(info_str, "(2X,A,F4.1,A)") '** This is more than the threshold of ', 30.0, ' percent.'
            call localorb_info(info_str)
            write(info_str, "(2X,'** Use smaller step to avoid collapse of unit cell.')")
            call localorb_info(info_str)

            call trm_trusted_step(F, hessian, trust_radius / i_trust_rad, DX_step)
            call X2coords(SCR_coords_sr, SCR_lv_sr,  X)
            call X2new_coords(SCR_coords_sr, SCR_lv_sr, X, DX_step)
            SCR_coords_sr = SCR_coords_sr / SCR_scaling
            if (SCR_n_params_coords > 0) then
              call SCR_backtransform(SCR_coords_sr, SCR_fcoords, .false.)
            end if
            if (SCR_n_params_lv > 0) then
              call SCR_backtransform_lv(SCR_lv_sr, SCR_lv)
            end if
            call frac2cart(SCR_lv, SCR_fcoords, coords)
            call min_atomic_dist(SCR_lv, coords, this_dist)
        end do

        lattice_vector_tmp = SCR_lv

        if (n_periodic > 0) then
           if (any(lattice_vector_tmp /= lattice_vector)) then
              lattice_vector = lattice_vector_tmp
              if (flag_cpu_consistency) &
                 call check_cpu_consistency_matrix(lattice_vector, &
                                         "lattice_vector_out", 3, n_periodic)
              call initialize_bravais_quantities()
           endif

           if (any(periodic_unit_cell_translations_tmp /= periodic_unit_cell_translations)) then
              periodic_unit_cell_translations = periodic_unit_cell_translations_tmp
           endif
        endif

        ! MOL: for symm-const relaxation do not shift
        periodic_unit_cell_translations = 0.

        write(info_str,'(2X,A)') "Updated atomic structure: "
        call localorb_info( info_str )

        call output_structure( )
        write(info_str,'(A)') &
             "------------------------------------------------------------"
        call localorb_info(info_str)

        ! Overwrites the lower triangle of the Hessian with values
        ! from the upper triangle. This was originally needed for the
        ! BFGS->TRM transition (see 432ca0ef). Currently the
        ! necessity of this subroutine is unclear.
        call symmetrize_hessian()

        ! TARP: Preserve the initial symmetry in the geometry.in.next_step_file for restarts
        if (SCR_n_params_coords > 0) then
          call SCR_backtransform(SCR_coords_sr, SCR_fcoords, .true.)
        endif
        call frac2cart(SCR_lv, SCR_fcoords, SCR_coords_tmp)
        call write_new_geometry_file(SCR_coords_tmp, lattice_vector)
      else

         if (abort_requested) then
           write(info_str,'(1X,A)') &
                 "* Abort request for optimization found *"
            call localorb_info(info_str)
            write(info_str,'(1X,A)') &
                 "* Aborting optimization."
            call localorb_info(info_str)
            write(info_str,'(A)') ''
            valid_geo = .false.

         else if (.not.final_geo_found) then

            write(info_str,'(1X,A)') &
                 "* Maximum number of relaxation steps exceeded."
            call localorb_info(info_str)
            write(info_str,'(1X,A)') &
                 "* Aborting optimization."
            call localorb_info(info_str)
            write(info_str,'(A)') ''
            valid_geo = .false.

         endif

         write(info_str,'(A)') &
              "------------------------------------------------------------"
         call localorb_info( info_str )
         if (flag_external_pressure) &
            call external_pressure_notice(external_pressure)
         write(info_str,'(2X,A)') "Final atomic structure: "
         call localorb_info( info_str )

         call output_structure( )

         write(info_str,'(A)') &
              "------------------------------------------------------------"
         call localorb_info( info_str )

      endif

      ! VA: Geometry has changed, calculate numerical stress again
      num_stress_finished = .false.



    else if (use_molecular_dynamics) then

      ! explicitly clean unitary translations and rotations in subroutine, this might change the forces slightly
      if (skip_SCF) call initialize_cleaning_forces()

      ! MR: trying to add a constraining potential
      if (use_potential_constrain) then
         call apply_constrain(pot_energy)
         if (myid.eq.0) write(use_unit,*) '| Energy from potential constrain:', pot_energy
      endif

      call clean_force_components(total_forces, forces_lv)
      ! CC: * updating the schedule before the time step instead of afterwards
      !       is inconsistent and can lead to problems when classical fields that
      !       that depend on the schedule are used
      !
      !if (MD_use_schedule) then
      !  if (tsystem.ge.MD_time) call change_MD_schedule_step
      !endif
      !
      !     * Also, it makes much more sense to check if the trajectory
      !       should end before yet another SCF is run.
      !
      !! check if file MD_stop exists
      !if (check_MD_stop) then
      !open(300,file = "MD_stop", status = 'old', iostat = MD_stop)
      !endif
      !! final_geo_found determines the point when the MD simulation is over
      !final_geo_found = (tsystem.ge.MD_time) .or. (use_MD_max_steps.and.(MD_stepcount.ge.MD_maxsteps)) .or. (MD_stop.eq.0)

      ! apparently the SCF cycle has converged if the code got this far
      valid_geo = .true.

      ! (re-)compute free energy of the present structure
      ! Note that it is the free energy which is consistent with
      ! the forces, not the total energy
      ! Check classical_field.f90 as well
      ! entropy already included in TDI total E
      if (use_thermodynamic_integration .or. use_reffree_AS) then
        free_energy = total_energy
      else
        free_energy = total_energy + 2.d0*entropy_correction
      endif

      ! one more force has been computed
      MD_force_evaluations = MD_force_evaluations + 1

      if ( (MD_stepcount.eq.0) .and. (.not.MD_successful_restart_read) ) then
         ! initialize only if it was not already done so by an MD restart file
         ! initial MD call, with all required initializations and the first MD step
         call initialize_MD(free_energy)

      else
         ! next MD step
         call MD_step(free_energy)
         if (.not. use_harmonic_pot_only) then
           if (MD_restart_binary) then
             call write_MD_restart_binary
           else
             call write_MD_restart
           endif
         endif
      endif

      call output_molecular_dynamics(free_energy)

      !CC: Update lambda for the next step
      !    and write restart info for the NEW lambda value
      if (use_thermodynamic_integration .or. use_reffree_AS) then
        call TDI_update_lambda()
        if (.not. use_harmonic_pot_only) call write_TDI_restart
      endif

      ! CC: see comments above
      if (MD_use_schedule) then
        if (tsystem.ge.MD_time) call change_MD_schedule_step()
      endif

      ! check if file MD_stop exists
      if (check_MD_stop) then
      open(300,file = "MD_stop", status = 'old', iostat = MD_stop)
      endif
      ! final_geo_found determines the point when the MD simulation is over
      final_geo_found = (tsystem.ge.MD_time) .or. (use_MD_max_steps.and.(MD_stepcount.ge.MD_maxsteps)) .or. (MD_stop.eq.0)

      ! Geometry has changed, calculate numerical stress again
      num_stress_finished = .false.

    else if(use_pathint) then        ! XZL: added for PIMD

     PIMD_beadcount = 1

     do i_bead = 1, n_beads

       if(PIMD_beadcount.eq.1) then

         converged_scf_bead = converged_scf

         call clean_force_components(total_forces, forces_lv)

         free_energy = total_energy + 2d0 * entropy_correction

         total_forces_beads(:,PIMD_beadcount,:) = total_forces(:,:)

         PIMD_beadcount = PIMD_beadcount + 1

       else

         coords(:,:) = coords_beads(:,i_bead,:)

         if (skip_SCF) then

           call reinitialize_noscf ( converged_scf_bead )

         else

           call reinitialize_scf ( converged_scf_bead )

         endif

         ! run scf cycle for updated geometry

         if (skip_SCF) then

           call noscf_solver ( converged_scf_bead , walltime_left_bead)

         else

           call scf_solver ( converged_scf_bead , walltime_left_bead)

         endif

         ! Calculate additional classical force field contributions

         call classical_field ( )

         ! Use history of earlier geometry steps to predict the
         ! next geometry

         call clean_force_components(total_forces,forces_lv)

         total_forces_beads(:,PIMD_beadcount,:) = total_forces(:,:)

         forces_lv_beads (:,PIMD_beadcount,:) = forces_lv(:,:)

         free_energy_beads(PIMD_beadcount) = total_energy + 2d0 * entropy_correction

         PIMD_beadcount = PIMD_beadcount + 1

       endif

     end do

     valid_geo = .true.

     if(PIMD_stepcount.eq.0) then

       call initialize_PIMD(free_energy)

     else

       call PIMD_step(free_energy)

       call write_PIMD_restart

     endif

     do i_atom = 1, n_atoms
       coords(:,i_atom) = coords_beads(:,1,i_atom)
     end do

     call output_pi_molecular_dynamics(free_energy)

     ! final_geo_found determines the point when the MD simulation is over
     final_geo_found = (tsystem_beads.ge.PIMD_time) .or.    &
     &                 (use_PIMD_max_steps.and.(PIMD_stepcount.ge.PIMD_maxsteps))

     ! Geometry has changed, calculate numerical stress again
     num_stress_finished = .false.

     PIMD_beadcount = 1

    else if (calculate_atom_bsse) then

        if (current_atom_num_bsse < n_atoms) then
          call get_timestamps(temp_time_bsse, temp_clock_time_bsse )
          time_bsse = temp_time_bsse
          clock_time_bsse = temp_clock_time_bsse
          current_atom_num_bsse=current_atom_num_bsse+1
          current_bsse_atom=current_atom_num_bsse
          !if (current_atom_num_bsse==0) then
          !   call cpu_time(rtime)
          !   time_bsse_per_atom= rtime
          !else
          !   call cpu_time(rtime)
          !   time_bsse_per_atom= rtime-time_bsse_per_atom
          !   write(use_unit,'(2X,A,F12.3,A)') &
          !  "Total time taken for the last BSSE step              : ", &
          !   time_bsse_per_atom, " s"
          !endif
          write(info_str,'(2X,A)') "Geometry for the next BSSE step "
          call localorb_info( info_str )
          call localorb_info("-----------------------------------------------------------", use_unit,'(2X,A)')
          n_occ_atoms=1! originally n_occ_atoms,n_empty_atoms was found while parsing geometry.in
          n_empty_atoms=n_atoms-1
          ! set charges to zero
          charge=0.d0
          ! set hf_version==HF_EIGEN
            hf_version = HF_EIGEN

          call output_bsse_structure(current_bsse_atom, n_electrons)
          call localorb_info("-----------------------------------------------------------", use_unit,'(2X,A)')
          !  to test part of the loop  :
          !if (current_atom_num_bsse< 2) then
               cont_bsse= .true.
          !endif


        endif! for looping over current atoms

    endif

   ! If unit cell changed in any way, write basic unit cell properties
   ! to make sure computational quantities such as volume are correctly
   ! udated.
   ! Change unit cell related parameters here.
   if (n_periodic > 0) then
      if (any(lattice_vector_orig /= lattice_vector)) then
         ! unit cell geometry changed.
         ! Notice that the variable unit_cell_shape_changed will be set to .false.
         ! after the first call to evaluate_hartree_recip_coef in hartree_potential_recip.f90 .
         unit_cell_shape_changed = .true.
         call output_bravais_quantities()
         if (Ewald_radius_automatic) then
            call estimate_Ewald_radius(n_atoms, cell_volume, Ewald_radius)
         endif
      else
         unit_cell_shape_changed = .false.
      endif
   endif

   ! Last operation of this subroutine:
   ! check cpu consistency of the updated geometry.
   ! if we perform a supercell calculation:
   if (n_periodic > 0 .and. flag_cpu_consistency) then
      ! check for consistent lattice vectors
      call check_cpu_consistency_matrix(lattice_vector, &
         "lattice for partition grid", 3, n_periodic)

      ! check for consistent map_to_center_cell_matrix
      call check_cpu_consistency_matrix(map_to_center_cell_matrix, &
         "map_to_center_cell_matrix for partition grid", 3, n_periodic)

      ! check for consistent coords
      call check_cpu_consistency_matrix (coords,"coords",3,n_atoms)

   endif


 end subroutine predict_new_geometry
