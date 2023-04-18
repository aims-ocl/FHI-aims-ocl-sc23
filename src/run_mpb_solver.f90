!****s* FHI-aims/mpb_solver_utilities
!  NAME
!   mpb_solver_utilities
!  AUTHOR
!   Stefan Ringe
!  REFERENCE
!   This routine is part of the MPBE solver for the modeling of-
!   electrolytes with FHI-aims. Large parts of this code are still-
!   experimental, that is why we highly recommend you to contact
!   the main authors of the corresponding paper (below) before starting
!   to calculate anything. We are highly motivated to help and-
!   cooperate with interested FHI-aims users.
!  SEE ALSO
!    Ringe, S., Oberhofer, H., Hille, C., Matera, S., Reuter, K., "Function-Space oriented Solution
!    Scheme for the-
!    Size-Modified Poisson-Boltzmann Equation in Full-Potential DFT", JCTC 12, 4052-4066 (2016).
!    DOI: 10.1021/acs.jctc.6b00435
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2016).
!  SYNOPSIS

subroutine run_mpb_solver(rho, rho_gradient, dielec_func)

    !purpose: organizes all flags of the MPBE/LPBE solver, for communication with SCF cycle

    use runtime_choices
    use SPE_solver
    use species_data, only:l_shell_max
    use physics, only:hartree_partition_tab, partition_tab
    use lpb_solver_utilities
    use mpi_utilities
    use physics, only:free_rho_superpos
    use constants, only:pi4_inv
    implicit none

    real*8 :: rho(n_spin, n_full_points)
    real*8 :: rho_gradient(3, n_spin, n_full_points)
    real*8 :: dielec_func(n_full_points)
    integer :: iter_steps
    !	SR: DEBUG FLAGS
    ! 	if (myid==0) then
    ! 	write(use_unit,*) 'flagsss',not_converge_rho_mpb,lpb_solver_going_to_start,mpb_solver_started,&
        ! 	  mpb_converged,converged_scf_mpb,converged_scf_el_mpb
    ! 	end if


    ! 	 if (cube_dielec_func_outputted) then
    ! 	    !output dielectric function in when conv crit is hit
    ! 	    if ((mpb_converge_crit.eq.'step'.and.start_at_scf_step-1.eq.number_of_scf_loops)) then
    ! 	      !only one step left, but we do not need the LPB hartree saving stuff in sum_up_whole_potential
    ! 	      lpb_solver_going_to_start = .false.
    ! 	    else if ((mpb_converge_crit.eq.'converge'.and.converged_scf_mpb).or.&
        ! 	      (mpb_converge_crit.eq.'step'.and.start_at_scf_step.eq.number_of_scf_loops)) then
    ! 	      !if conv crit is reached, output dielec_func and stop
    ! 	      mpb_converged = .True.
    ! 	      mpb_solver_started = .False.
    ! 	      lpb_solver_going_to_start = .False.
    ! 	      converged_scf_mpb = .True.
    ! 	      converged_scf_el_mpb = .True.
    ! 	      !sum_up_whole_potential will now skip all lpb part, but let the outer SCF
    ! 	      !stop in the next step
    ! 	      return
    ! 	    end if

    ! 		if (solve_lpb_with_constant_dielec.and.mpb_solver_started) then
    ! 		      call solve_SPE(.false.,sum(rho,DIM=1)- pi4_inv * free_rho_superpos, delta_v_MPB(1), delta_v_MPB_gradient(1,1), &
        ! 			.False.,0.0d0,0.0d0,.True.,.True.,.False.,.false.,&
        ! 			.false.,iter_steps,.false.,delta_v_MPB_laplacian,0,.false.)
    ! 		      call output_delta_v_MPB_step(delta_v_MPB,1666)
    ! 		      stop
    ! 		end if
    ! !
    ! ! 	 else if (not_converge_rho_mpb) then
    ! 	 if (not_converge_rho_mpb) then
    !
    ! 	    if (mpb_solver_started .and..not.restart_read) then
    ! 	      !as soon as lpb_solver is started, stop after one run (no response of rho to solvent)
    ! 	      mpb_converged = .True.
    ! 	      converged_scf_mpb = .True.
    ! 	      converged_scf_el_mpb = .True.
    ! 	      if (.not.solve_lpbe_only) then
    ! 		call run_newton_mpb
    ! 	      else
    ! 		call lpb_solver(.false.,0d0,0d0,0d0,0d0)
    ! 	      end if
    ! 	      return
    ! 	    else if (restart_read) then
    ! 	      !as soon as lpb_solver is started, stop after one run (no response of rho to solvent)
    ! 	      mpb_converged = .True.
    ! 	      converged_scf_mpb = .True.
    ! 	      converged_scf_el_mpb = .True.
    ! 	      if (.not.solve_lpbe_only) then
    ! 		call run_newton_mpb
    ! 	      else
    ! 		call lpb_solver(.false.,0d0,0d0,0d0,0d0)
    ! 	      end if
    ! 	      return
    ! 	    end if

    if (mpb_solver_started) then

        !		SR: DEBUG: plot electronic hartree potential
        ! 		i_full_points = 0
        ! 		do i_my_batch = 1, n_my_batches, 1
        !
        ! 		      do i_index = 1, batches(i_my_batch)%size, 1
        ! 			i_full_points = i_full_points + 1
        ! 			call  tab_single_atom_centered_coords_p0( 1, batches(i_my_batch) % points(i_index) % coords(:), &
            ! 			    dist_tab_sq, dir_tab )
        ! 			hartree_potential_nuc(i_full_points) = 1.0/sqrt(dist_tab_sq)
        ! 		      end do
        ! 		end do
        ! 		call output_delta_v_MPB(hartree_potential+hartree_potential_nuc)
        ! 		stop
        ! 		END DEBUG

        ! 		if (.not.mpb_converged) then
        ! 		  mpb_converged = .True.
        ! 		end if

        ! go on and run lpb_solver, because not converged, yet...
        if (solve_lpb_with_constant_dielec) then
            call solve_SPE(.false., rho(1, 1), delta_v_MPB(1), delta_v_MPB_gradient(1, 1), &
                           .False., 0.0d0, 0.0d0, .True., .True., .false., &
                           .false., iter_steps, .false., delta_v_MPB_laplacian, 0, .false.)
        else
            if (.not. solve_lpbe_only) then
                call run_newton_mpb
            else
                call lpb_solver(.false., 0d0, 0d0)
            end if
        end if

    end if

    return

end subroutine run_mpb_solver
