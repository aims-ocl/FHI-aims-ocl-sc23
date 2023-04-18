
!****s* FHI-aims/aims
!  NAME
!    aims
!  SYNOPSIS

    subroutine aims(mpi_comm_input,in_unit,mpi_switch)

      !  PURPOSE
      !    This subroutine is the main (high-level) file for the FHI-aims executable.
      !    All tasks branch out from this subroutine:
      !    * Initialization
      !    * Reading of all input data
      !    * Preparation, initialization and execution of the first full scf cycle
      !    * Loop(s) over further scf cycles for updated geometries
      !     (structure optimization or molecular dynamics)
      !    * Post-processing tasks (Charge analyses, perturbative corrections, ...)
      !
      !    There should be as few "use" statements as possible here. Any meaningful work
      !    should be done in high-level subroutines that are well separated from one
      !    another.
      !
      !  USES
      
            use dimensions
            use localorb_io
            use xml_write, only: xml_close_file
            use applicable_citations
            use mpi_utilities
            use runtime_choices
            use aims_memory_tracking, only: aims_mem_initialize
            use timing
            use tddft_real_time_propagation
            use ipc, only: ipc_end_session
            use ipc_hamiltonian_and_ovl_transfer, only: initiate_hamiltonian_and_ovl_ipc_transfer
            use generate_aims_uuid, only: write_aims_uuid
            use aims_gpu, only: initialize_gpu, finalize_gpu
            use json_output, only: open_json_log_file, close_json_log_file, &
                write_geometry_to_json_log, write_final_output_to_json_log
            use WannierCenters, only: WCC_calc
            use opencl_util, only: read_opencl_settings, mpi_platform_relative_id, mpi_per_node
      
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
      
      
      !  Declaration of variables
      
            implicit none
      
            ! Local variables
      
            ! "converged_scf" indicates whether or not the scf cycle for a single
            ! geometry was converged. It is now passed around directly.
            logical :: converged_scf, converged_cpscf
      
            ! "final_geo_found" indicates whether a possible geometry relaxation
            ! is converged (i.e. whether the magnitude of the forces is below
            ! the prescribed threshold), or whether an MD simulation is over.
            logical :: final_geo_found = .false.
      
            ! "valid_geo" indicates whether a new valid geometry could be found
            ! during relaxation
            logical :: valid_geo
      
            ! "walltime_left" indicates that the program does not have to quit
            ! due to walltime restrictions
            logical :: walltime_left
      
            ! indicates whether a bsse run should be continued
            logical :: cont_bsse
      
            character*150 info_str
            character(LEN=80) ::  uuid_str
      
            ! FHI-aims library version input variables
            integer,intent(in) :: mpi_comm_input
            integer,intent(in) :: in_unit
            logical,intent(in) :: mpi_switch
            integer :: Z2_i_plane
            ! CC: Interface to make out_plot_band callable with optional 
            !     arguments despite being a subroutine and no module
            interface
              subroutine out_plot_band (Z2_i_plane )
                    integer,intent(in),optional :: Z2_i_plane
              end subroutine out_plot_band 
            end interface
              
      
            logical :: scf_reinit
            integer :: n_loops
            integer :: info !Gets the return from allocate statements
            ! FHI-aims library version input variables
            !set the global communicator
            use_mpi = mpi_switch
            mpi_comm_global = mpi_comm_input
            call get_my_task()
            !set the default output unit
            default_unit = in_unit
            use_unit     = in_unit
            !Begin work
      
            ! call crts_init()
      
            ! Write Preamble at the beginning of the the output file
            ! E.g. greet the user, uuid identifier, initial_mpi_report, etc.
            call write_preamble()
      
            ! Initialize list of applicable references to be written out at end of code.
            call init_citations
      
            ! Initialize known debug modules
            call init_debug
      
            ! register FHI-aims CPC citation for final citations
            call cite_reference("FHI-aims-CPC")
      
            ! Set default values for variables in various modules
            call set_aims_defaults ( )
      
            ! Start the memory tracking infrastructure
            call aims_mem_initialize ( )
      
            ! if(myid==0) write(use_unit,*) 'setting stacksize to maximum now...'
            !call unlimit_stack()
      
            ! do tests which can be useful to prevent bad performance and/or
            ! erroneous results
            call prepare_system
      
            ! Read all input files:
            ! * Set up initial dimensions
            ! * allocate initial data structures
            ! * read all input information from files control.in and geometry.in
      
            call read_input_data ( )
      
            ! call hip_and_hipblas_init()
            ! Initialize GPU, if the user requested it
            if (use_gpu) call initialize_gpu ()

            ! read opencl settings
            mpi_platform_relative_id = mod(myid, mpi_per_node)
            ! call read_opencl_settings()
      
            ! Set up all fixed computational quantities:
            ! * Basic integration grid specifications for all species
            ! * Free-atom DFT solution for each species on a logarithmic grid
            ! * Radial basis functions for each species, tabulated on log. grids
            ! * Coefficients that prepare the computation of Ylm functions and
            !   their derivatives from cartesian products x^l*y^m*z^n
            call prepare_scf ( )
      
            ! If this is a pure parser executable OR if the 'dry_run' flag was set, stop now.
            ! In that case, the sole purpose of the run was to check the consistency of
            ! the input files 'control.in' and 'geometry.in'.
            call stop_if_parser()
      
      !-----Start the SCF cycles----------------------------------------------------------------
      
            ! Initial iteration of the scf loop, the one that generates
            ! the grids and the ovp. matrix (initialize_integrals), the first
            ! wave function, all free-atom derived quantities, etc.
            if (skip_SCF .or. use_pimd_wrapper) then
              call initialize_noscf ( converged_scf )
            else
              call initialize_scf ( converged_scf )
            end if
      
            if (use_rrs_pbc) then
              call parse_rrs_pbc()
            end if
      
            ! Open JSON log file
            if (out_aims_json_log) then
              call open_json_log_file()
              call write_geometry_to_json_log("initial")
            end if
      
            walltime_left   = .true. ! Set .true. for the first initialization to enter loop
            n_loops = 0
            scf_reinit = .false.
      
            SCF_INITS_LOOP: do while (((scf_reinit.or.n_loops==0).and. walltime_left))
            ! n_loops == 0: corresponds to first initialization
            ! n_loops  > 0: SCF is reinitialized (E.g.: new geometry, bad SCF convergence, ...)
      
              if (n_loops > 0) then
                ! reinitializes all geometry-dependent quantities in scf, but
                ! leaves the basic integration grids and the wave function intact
                if (skip_SCF) then
                  call reinitialize_noscf ( converged_scf )
                else
                  call reinitialize_scf ( converged_scf )
                end if
              endif
      
              if (skip_SCF) then
                call noscf_solver ( converged_scf , walltime_left)
      
              else if (use_pimd_wrapper .and. n_loops == 0) then
                ! MR: if doing the pimd wrapper, we want to skip the normal scf.
                call run_driver(converged_scf, walltime_left)
      
              else if (skip_scf_for_postscf .and. n_loops == 0) then
                ! IYZ: skip scf for postscf calculations
                ! SK: The if line above looks like a dev hack ... is it still needed?
                walltime_left = .true.
                converged_scf = .true.
      
              else
                ! Start the SCF cycle
                call scf_solver ( converged_scf , walltime_left)
              end if
      
              ! RRS-PBC scheme, igor
              if (use_rrs_pbc .and. n_loops == 0) call run_rrs_pbc()
      
              if(calculate_atom_bsse) then
                 if (use_rpa_ene) then
                   if (n_loops == 0) call prepare_corr_energy_calc ()
                   call qpe_calculation ( )
                 elseif(use_mp2 .or. use_dftpt2 .or. use_os_mp2) then
                   if (n_loops == 0) call prepare_corr_energy_calc ()
                   call mp2_calculation ( )
                 endif
              endif
      
              ! Do the DFPT-related calculations if requested
              ! E.g. vibrations, phonons, polarizability, dielectric, friction
              if ((.not. use_pimd_wrapper .and. n_loops == 0) .or. n_loops > 0) then
                call DFPT_wrapper(converged_cpscf,converged_scf)
              endif
      
              ! Calculate additional classical force field contributions
              call classical_field ( )
      
              ! Decide here whether to reinitialize the SCF
              if (restarting_scf) then
                scf_reinit = .true.
                call localorb_info( &
                "No geometry change - only restart of scf mixer after initial iterations.",use_unit,'(2X,A)' )
              else if (use_geo_relaxation.or.use_molecular_dynamics.or.use_pathint.or.&
                       use_numerical_stress.or.calculate_atom_bsse) then
                call predict_new_geometry ( converged_scf, final_geo_found, valid_geo, cont_bsse )
                scf_reinit = (valid_geo.and.(.not.final_geo_found)).or.cont_bsse
              else
                scf_reinit = .false.
              endif
      
              n_loops = n_loops + 1
            end do SCF_INITS_LOOP
      
      !-----post-production begins here--------------------------------------------------------
      
            if ( use_out_eigenvec_ovlp ) then
              ! output eigenvector and overlap matrix (for periodic system)
              ! into unformatted files
              call out_eigenvec_ovlp_p1 ()
            end if
      
            ! TDDFT real-time propagation
            if ( flag_run_tddft_real_time_propagation ) then
              !call reinitialize_scf ( converged_scf )
              call tddft_real_time_propagation_init()
              call tddft_real_time_propagation_run()
              call tddft_real_time_propagation_end()
            end if
      
            ! check whether anything remains to be done!!!
            if (.not.walltime_left.and.(.not.final_geo_found)) then
              call report_out_of_walltime_error()
              call final_timings_and_memory(.false.)
              if (out_aims_json_log) then
                call write_geometry_to_json_log("final")
                call write_final_output_to_json_log(.false.)
                call close_json_log_file()
              end if
              call deallocate_mpi()
              call finalize_mpi()
              stop
            end if
      
            if ( flag_run_mulliken .and. .not. calculate_perturbative_soc) then
              ! write Mulliken charge analysis to standard output,
              ! and full list of Mulliken partial charges to separate file
              call mulliken_analysis ( converged_scf )
            end if
      
            if ( out_lowdin .and. flag_run_lowdin ) then
              call lowdin_analysis ( converged_scf )
            end if
      
            if (out_hirshfeld) then
              ! Perform Hirshfeld analysis of partial charges
              call hirshfeld_analysis ( )
            end if
      
            if (out_hirshfeld_iterative) then
              ! Perform Hirshfeld-I analysis of partial charges
              call hirshfeld_analysis_iterative ( )
            end if
      
            if (out_esp.or.out_esp_full) then
              ! Fitting of ESP charges to Long Range Potential
              call get_esp_charges ( )
            end if
      
            if (use_hf_post) then
               call initialize_hartree_fock()
               call hf_postproc()
            end if
            if (n_periodic==0 .and. use_corr .and. (.not.calculate_atom_bsse)) then
              ! if (use_dftpt2 .and. use_hse) use_dftpt2_and_hse=.true.
              call prepare_corr_energy_calc ()
      !        if (use_qpe .and. (.not. neutral_excitation_bse)) then
              if (use_qpe) then
                 call qpe_calculation ( )
              endif
              if (use_os_mp2_qpe .or. use_os_rpa_qpe .or. use_sic_rpa_qpe) then
                 call qpe_osPT2 ( )
              endif
              if(use_scgw.or.use_scgw0)then
                 call self_consistent_gw ( )
              endif
              if(use_mp2 .or. use_dftpt2 .or. use_os_mp2) then
                 call mp2_calculation ( )
              endif
              if(use_ci) then
                 call ci_calculation ( )
              endif
              if (use_cc) then
                 Call CC_main()
              endif
              if(use_bare_ci) then
                 call evaluate_bare_ci ( )
              endif
              if (use_dftpt2_and_lrc) then
                 call prepare_lrc_corr_energy_calc ()
                 call mp2_calculation ( )
              endif
           end if
      
           if(neutral_excitation_bse) then
              call bse()
           endif
      
           if (output_ks_coulomb_integral) then
              call output_coulomb_integrals_mo()
           endif
      
           if(use_dmft_gw) then
            call self_consistent_DMFT ( )
      
           elseif(use_dmft_pbe0) then
            call self_consistent_DMFT_PBE0 ( )
           endif
      
            if (calculate_perturbative_soc) then
              call calculate_second_variational_soc ()
            end if
      
            if ( flag_run_mulliken .and. calculate_perturbative_soc) then
              ! write Mulliken charge analysis to standard output,
              ! and full list of Mulliken partial charges to separate file
              call mulliken_analysis ( converged_scf )
            end if
      
            if(flag_out_dielectric .or. flag_out_absorption) then
              call dielectric_analysis ( converged_scf )
            endif
      
           if(n_periodic > 0 .and. (use_rpa_ene .or. use_periodic_gw)) then
            call post_scf_correlation_treatment_p0()
           else if(n_periodic > 0 .and. (use_mp2 .or. use_dftpt2) .and. use_mp2_blacs) then
            call post_cpt2_blacs()
           else if(n_periodic > 0 .and. (use_mp2 .or. use_dftpt2 .or. use_os_mp2)) then
            call post_cpt2_lapack()
           endif
      
            if (magnetic_response) call MR_main()
      
            if (out_cube) then
              call cube_output ( )
            end if
            call mpi_barrier(mpi_comm_global,info)
      
            ! CC: output routine to generate files for DGrid
            !  interface developed by: Alim Ormeci - MPI CPfS, Dresden (August 2018)
            if ( out_dgrid ) then
               call output_dgrid()
            end if
      
            ! At this point, all (periodic) calculations that rely on an existing
            ! list of k-points are completed. The remaining routines might reset
            ! the k-point list for their own purposes.
      
            ! plot density of states if a denser k-point grid was requested for the plot
            ! (post-processing only)
            if (pert_dos_on) then
              ! pert_dos_on can only be true if out_dos is true or
              ! if postscf_eigenvalues is true. This is ensured in read_control.f90.
              ! Please check read_control.f90 to understand and possibly modify these restrictions.
              call output_KS_dos_pert()
            end if
      
            call mpi_barrier(mpi_comm_global,info)
      
            call initiate_hamiltonian_and_ovl_ipc_transfer()
      
            ! plot band structure
            ! plot band structure
            if (out_band_mulliken) then
              if (use_scalapack) then
                if (calculate_perturbative_soc) then
                  write(info_str, *) "Output mulliken decomposition of band k points for SOC in ScaLAPACK"
                  call localorb_info(info_str)
                  call out_plot_band_soc_mulliken_scalapack()
                else
                  write(info_str, *) "Output mulliken decomposition of band k points in ScaLAPACK"
                  call localorb_info(info_str)
                  call out_plot_band_mulliken_scalapack ( )
                end if
              else
                if (calculate_perturbative_soc) then
                  write(info_str, *) "Output mulliken decomposition of band k points for SOC in LAPACK"
                  call localorb_info(info_str)
                  call out_plot_band_soc_mulliken()
                else
                  write(info_str, *) "Output mulliken decomposition of band k points in LAPACK"
                  call localorb_info(info_str)
                  call out_plot_band_mulliken()
                end if
              end if
            end if
      
            if (out_band .and. (.not. use_gw)) then
              if (use_scalapack) then
                call out_plot_band_scalapack ( )
              else
                if(use_hf_kspace .or. use_hf_realspace) then
                   if (exx_band_structure_version.eq.2) then
                      call out_plot_band_hf_k_space ( )
                   else
                      call out_plot_band ( )
                   endif
                else
                  call out_plot_band ( )
                endif
              end if
            end if
      
            ! CC: Wannier Centers / adapted from Z2 implementation of Cmera 
            if (.not.use_local_index  .and. out_z2_invariant .and. (.not. use_gw)) then
      
              call localorb_info('')
              write(info_str,'(2X,A)') "Starting with the calculation of the Wannier Centers:"
              call localorb_info(info_str)
      
              if (use_scalapack) then
                write(info_str, *) "**** SCALAPACK MISSING SO FAR ************** "
                call localorb_info(info_str)
                write(info_str, *) "** Cannot perform the computation of Wannier centers. Sorry."
                call localorb_info(info_str)
              else
                WCC_calc = .true.
                do Z2_i_plane = 1, Z2_n_planes
                  write(info_str,'(2X,A,I2,A,I2,A)') " |-> Using LAPACK Code for plane" , &
                    & Z2_i_plane, " of ",Z2_n_planes," planes"
                  call localorb_info(info_str)
                  call out_plot_band ( Z2_i_plane )
                end do
              end if
            end if
      
            if (calculate_atom_bsse) then
               call atom_bsse_results
            endif
      
            if (out_aims_json_log) then
              call write_geometry_to_json_log("final")
              call write_final_output_to_json_log(converged_scf)
              call close_json_log_file()
            end if
      
            call mpi_barrier(mpi_comm_global,info)
      
            call final_deallocations ( )
      
            call mpi_barrier(mpi_comm_global,info)
      
            call final_energy_output ( )
      
            if (out_xml) then
                call xml_close_file ()
            end if
      
            call print_references ( )
            call cleanup_citationlist ( )
      
            !Make sure 'Have a nice day' is only told if the scf is converged
            call final_timings_and_memory(converged_scf)
      
            if (use_gpu) call finalize_gpu()
      
            grid_storage_initialized = .false.
            converged_scf            = .false.
            first_integration        = .true.
            got_n_compute_maxes      = .false.
      ! end work - that's it
      
            CALL ipc_end_session()
      
      ! AJL/Jan13 Can this be moved somewhere tidier?
      ! Needed in case number of atoms changes in QM/MM
            call deallocate_mpi ( )
      
            ! call athread_halt() 
      
          end subroutine aims
      !******
      