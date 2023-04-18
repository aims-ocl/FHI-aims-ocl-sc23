!****s* FHI-aims/run_newton_mpb
!  NAME
!   run_newton_mpb
!  AUTHOR
!   Stefan Ringe
!  REFERENCE
!   This routine is part of the MPBE solver for the modeling of
!   electrolytes with FHI-aims. Large parts of this code are still
!   experimental, that is why we highly recommend you to contact
!   the main authors of the corresponding paper (below) before starting
!   to calculate anything. We are highly motivated to help and
!   cooperate with interested FHI-aims users.
!  SEE ALSO
!    Ringe, S., Oberhofer, H., Hille, C., Matera, S., Reuter, K., "Function-Space oriented Solution Scheme for the
!    Size-Modified Poisson-Boltzmann Equation in Full-Potential DFT", JCTC, 2016 12, 4052-4066 (2016).
!    DOI: 10.1021/acs.jctc.6b00435
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2016).
!  SYNOPSIS
subroutine lpb_solver(use_mpbe_params,h_function_calc,q_times_epsinf_mpb_calc)
!  PURPOSE
!    Solves a kind of linearized Poisson-Boltzmann equation (LPBE) via the Multipole Expansion Relaxation Method (MERM) using
!    a self-consistency cycle and multipole expansions in each step. The multipole expansions are performed via the SPE_solver.
!  USES

    use physics, only: partition_tab,free_rho_superpos,&
        hartree_potential, hartree_partition_tab,rho,rho_gradient,free_hartree_superpos,&
        free_rho_gradient_superpos
    use lpb_solver_utilities, only: &
        dielecfunc_from_density,alphafunc_from_density,dielecfunc_gradient_from_density_gradient,&
        evaluate_q_times_epsinf,evaluate_SPE_source,evaluate_rmsd_merm,lpb_converged,&
        initialized_lpb_solver,delta_rho_in_merm,&
        dynamic_cavity,epsinf_mpb,tau_merm,eta_merm, delta_v_MPB_gradient, delta_v_MPB,&
        qSPE_old , qSPE_plus_delta_rho,kappa_debye,&
        dielec_func_mpb, delta_v_MPB_laplacian,&
        v_free_gradient_public,dielec_func_gradient_mpb,z_mpb,kBT_mpb, &
        dielec_func2,dielec_func_grad2,&
        kappainf_mpb,correct_dielectric_decrement,decrement_slope,start_dielec_decr_corr_acc,&
        solve_lpbe_only,MERM_in_SPE_solver,qSPE,use_mpbe_free_energy,&
        dynamic_ions,phi_zero_mpb,alpha_func_mpb, rhomin_mpb, rhomax_mpb,rhomin_kappa,rhomax_kappa,&
        d_alpha_ion,xi_alpha_ion,forces_mpb_force_off,mpb_forces_on
    use mpb_solver_utilities, only: ln_cosh_function,evaluate_ln_cosh_function
    use SPE_solver, only: solve_SPE, prepare_SPE_solver
    use constants, only : pi4_inv
    use localorb_io
    use mpi_tasks, only: aims_stop, myid, mpi_wtime, mpi_comm_global
    use synchronize_mpi_basic, only: sync_real_number
    use grids
    use dimensions
    use species_data
    use geometry, only: species

    implicit none

    !MERM_in_SPE_solver:	the sclpbe cycle is performed withing the multipole representations
    !			without full update of potential on whole aims grid which gives additional speed
    !			especially for c^s=0
    !

    !LOGICAL
    character*200 :: info_str

    integer :: i_part_step, i_iterate_lpb
    integer :: i_spin, i_full_points, i_coords
    integer :: total_iter_steps
    integer :: info

    real*8 :: rmsd_merm !root mean square change in potential
    real*8 :: delta_rho_lpb(n_full_points) !delta density rho - rho_free
    real*8, dimension(n_full_points) :: delta_v_MPB_old
    real*8, dimension(n_full_points) :: q_times_epsinf
    real*8, dimension(n_full_points) :: qSPE_delta
    real*8, dimension(n_full_points) :: rho_spinless
    real*8 :: grid_coord(3)

    logical :: epsinfinv_factor = .False. !if a 1/eps_inf factor should be included in the calculation
    logical :: use_partition = .False. !enable this, if the MPE solver should be called with specific eps and kappa
    !of the nuclei potential, to use a specific eps and q set: use_partition = .true. and:
    !	  call reinitialize_green_mpb(precondition_q)
    ! 	  call solve_SPE(rho, potential, potential_gradient, &
        ! 	    hartree_partition_tab,  l_shell_max,&
        ! 	    use_partition, precondition_q, epsinf_mpb,.False.,.False.,&
        ! 	    calculate_electronic_potential,partition_tab,.true.)

    ! Timing
    real*8 time0, time_work, time_all

    !MPBE stuff
    logical :: use_mpbe_params
    real*8, dimension(n_full_points) :: q_times_epsinf_mpb_calc !if use_mpbe_params=true, this is the part of the right hand side of the LPBE, which is updated in each Newton step
    real*8, dimension(n_full_points) :: h_function_calc !hfunction of Newton method for MPBE


    if (solve_lpbe_only) then
        write(info_str,'(2X,A)') 'Starting LPBE Solver'
        call localorb_info(info_str,use_unit,'(A)',OL_norm)
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !TIMING
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!
    time0 = mpi_wtime()
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !EVALUATE DENSITY AND DIELECTRIC FUNCTION
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We do not need a spin-independent density, so sum over all spin channels
    rho_spinless = 0.0d0
    do i_spin = 1, n_spin, 1
        rho_spinless(:) = rho_spinless(:) + rho(i_spin,:)
    end do
    delta_rho_lpb = rho_spinless - pi4_inv * free_rho_superpos
    !if we want to update the eps in each step from the current electron density, update it here
    if (dynamic_cavity.and.solve_lpbe_only) then
        call dielecfunc_from_density(rho,dielec_func_mpb)
        call dielecfunc_gradient_from_density_gradient(dielec_func_mpb,dielec_func_gradient_mpb,&
            rho, rho_gradient)
    end if
    if (dynamic_ions.and.solve_lpbe_only) then
        call alphafunc_from_density( rho, rho_gradient, dielec_func_mpb, n_full_points)
    end if
    !call output_delta_v_step(sum(rho,DIM=1),10)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !INITIALIZATION
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    total_iter_steps = 0
    !
    !use already initialized greens function (in initialize_scf)
    use_partition = .False.

    if (.not.initialized_lpb_solver) then

        lpb_converged = .False.
        qSPE = 0.0d0
        qSPE_old = 0.0d0
        q_times_epsinf = 0.0d0
        delta_v_MPB_old = 0.0d0
        delta_v_MPB = 0.0d0
        delta_v_MPB_gradient = 0.0d0

        !if .not.MERM_in_SPE_solver:
        !	initialize solver with q_iter of LPBE and delta_v_0 = 0
        !else
        !	initialize solver with q_iter of MPBE and delta_v_0 = delta_v_lpbe
        if (use_mpbe_params.and..not.MERM_in_SPE_solver) then
            call evaluate_q_times_epsinf(rho_spinless,  &
                v_free_gradient_public,&
                dielec_func_mpb,dielec_func_gradient_mpb,  q_times_epsinf,&
                delta_rho_lpb,q_times_epsinf_mpb_calc,use_mpbe_params)
            call evaluate_SPE_source(q_times_epsinf, dielec_func_mpb,&
                dielec_func_gradient_mpb, &
                delta_v_MPB, delta_v_MPB_gradient, qSPE,  &
                .False.,0.0d0,use_mpbe_params,h_function_calc)
        else
            call evaluate_q_times_epsinf(rho_spinless,  &
                v_free_gradient_public,&
                dielec_func_mpb,dielec_func_gradient_mpb,  q_times_epsinf,&
                delta_rho_lpb,q_times_epsinf_mpb_calc,use_mpbe_params)
        end if

        if (.not.MERM_in_SPE_solver) then
            !we solve here SPE the first time, to get an initial value
            if (.not.use_mpbe_params) qSPE = q_times_epsinf
            call solve_SPE(.false., qSPE, delta_v_MPB(1), delta_v_MPB_gradient(1,1), &
                use_partition,0.0d0,0.0d0, .False.,epsinfinv_factor,&
                .true.,.false.,total_iter_steps,.false.,&
                delta_v_MPB_laplacian,0,use_mpbe_params)
        end if

    else !if already initialized
        call evaluate_q_times_epsinf(rho_spinless,  &
            v_free_gradient_public,&
            dielec_func_mpb,dielec_func_gradient_mpb, q_times_epsinf,&
            delta_rho_lpb,q_times_epsinf_mpb_calc,use_mpbe_params)
        lpb_converged = .False.
    end if !initialized_lpb_solver

    rmsd_merm = 0d0
    i_iterate_lpb = 0


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !RUN SOLVER
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (MERM_in_SPE_solver) then

        call solve_SPE(.false.,q_times_epsinf, &
            delta_v_MPB, delta_v_MPB_gradient, &
            use_partition, 0d0, 0d0,.False.,epsinfinv_factor,&
            .true.,.true.,total_iter_steps,.false.,&
            delta_v_MPB_laplacian,0,use_mpbe_params)

        initialized_lpb_solver = .True.

        !we take as input the effective density calculated with the potential of the previous step
        call evaluate_SPE_source(q_times_epsinf, dielec_func_mpb,&
            dielec_func_gradient_mpb, &
            delta_v_MPB, delta_v_MPB_gradient, qSPE,  &
            .False.,0.0d0,use_mpbe_params,h_function_calc)
        ! 	end if
        qSPE_plus_delta_rho(:) = (rho_spinless(:) - pi4_inv * free_rho_superpos(:)) -&
            qSPE(:)/epsinf_mpb!*1d0/dielec_func_mpb(:)


    else !ifnot MERM_in_SPE_solver

        if (solve_lpbe_only) then
            write(info_str,*) ''
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            write(info_str,'(2X,A)') ' |SCLPBE Initialization finished.'
            call localorb_info(info_str,use_unit,'(A)',OL_norm)


            write(info_str,*) ''
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            write(info_str,'(2X,A)') ' |Solving LPBE...'
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            write(info_str,*) '--------------------------'
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            write(info_str,*) ''
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
        end if

        do while (.not.(lpb_converged)) !.or.i_iterate_lpb.lt.2)
            ! 	  call system_clock(counti,count_rate)
            i_iterate_lpb = i_iterate_lpb + 1



            if (solve_lpbe_only) then
                write(info_str,'(2X,A,I4)') ' |Solving LPB equation with step: ', i_iterate_lpb
                call localorb_info(info_str,use_unit,'(A)',OL_norm)
            end if

            delta_v_MPB_old = delta_v_MPB

            if (correct_dielectric_decrement.and.solve_lpbe_only.and.&
                    (rmsd_merm.lt.start_dielec_decr_corr_acc).and.rmsd_merm.ne.0d0) then
                !we have to update the dielectric function in each step dependent on the current ion distribution
                if (dynamic_cavity) then
                    call dielecfunc_from_density(rho,dielec_func_mpb)
                    call dielecfunc_gradient_from_density_gradient(dielec_func_mpb,dielec_func_gradient_mpb,&
                        rho, rho_gradient)
                else
                    call dielecfunc_from_density(free_rho_superpos*pi4_inv,dielec_func_mpb)
                    call dielecfunc_gradient_from_density_gradient(dielec_func_mpb,dielec_func_gradient_mpb,&
                        free_rho_superpos*pi4_inv, free_rho_gradient_superpos)
                end if
            end if


            qSPE_old = qSPE

            !call output_delta_v_MPB_step(q_times_epsinf,2000+i_iterate_lpb)

            !   ! 	  !obtain rho_iter_n(u_{n-1},grad u_{n-1})
            call evaluate_SPE_source(q_times_epsinf, dielec_func_mpb,&
                dielec_func_gradient_mpb, &
                delta_v_MPB, delta_v_MPB_gradient, qSPE,  &
                .False.,0.0d0,use_mpbe_params,h_function_calc)


            !mix iterative density with previous step

            qSPE_delta = qSPE - qSPE_old
            qSPE = qSPE_old + eta_merm * qSPE_delta

            if (delta_rho_in_merm) then
                !we want to integrate the change in rho_iter
                qSPE = qSPE - qSPE_old
            end if


            !  	   !solve SPE with iterative density qSPE. output: solution of SPE delta_v_MPB and gradient delta_v_MPB_gradient
            call solve_SPE(.false.,qSPE(1), delta_v_MPB(1), delta_v_MPB_gradient(1,1), &
                use_partition, 0d0, 0d0,.False.,epsinfinv_factor,&
                .true.,.false.,total_iter_steps,&
                .false.,delta_v_MPB_laplacian,0,use_mpbe_params)



            if (delta_rho_in_merm) then
                delta_v_MPB = delta_v_MPB_old + delta_v_MPB
            end if

            call evaluate_rmsd_merm(delta_v_MPB(1), delta_v_MPB_old(1), rmsd_merm, tau_merm, lpb_converged)

            if (solve_lpbe_only) then
                write(info_str,'(2X,A,ES14.3)') '|Current RMSD: ', rmsd_merm
                call localorb_info(info_str,use_unit,'(A)',OL_norm)
            end if

            if (rmsd_merm.gt.1d6) then
                if (myid==0) write(use_unit,'(6(A,F10.6,$))') 'The following parameters were applied, nmin =',rhomin_mpb, ', nmax =',rhomax_mpb, ' ,nmin^alpha = ',&
                    rhomin_kappa, ' ,nmax^alpha = ',rhomax_kappa, ' ,d_alpha = ',d_alpha_ion,' ,xi_alpha = ',xi_alpha_ion
                call aims_stop('RMSD was above 1d6, SCLPB cycle may not converge. Reconsider your grid and angular momentum settings.')
            end if


            total_iter_steps = total_iter_steps + 1


            !call output_delta_v_step(delta_v_MPB,100+total_iter_steps)

        end do !sclpb loop

        initialized_lpb_solver = .true.

        qSPE_plus_delta_rho(:) = (rho_spinless(:) - pi4_inv * free_rho_superpos(:)) -&
            qSPE(:)/epsinf_mpb!*1.d0/dielec_func_mpb(:)

    end if !solve in mpe


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !FINAL OUTPUTS AND EVALUATIONS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !       call output_delta_v_MPB_step(-pi4_inv*alpha_func_mpb*kappainf_mpb**2*z_mpb/kBT_mpb*&
        !  		1d0/(1d0+phi_zero_mpb*(alpha_func_mpb-1d0))*(delta_v_MPB+free_hartree_superpos),201)
    ! 	call output_delta_v_MPB_step(delta_v_MPB,301)
    ! !call output_delta_v_step(dielec_func_mpb,9501)
    !       call output_delta_v_MPB_step(sqrt(sum(delta_v_MPB_gradient**2,DIM=1)),801)
    !       call output_delta_v_MPB_step(delta_rho_lpb,133)
    !       call output_delta_v_MPB_step(sum(rho,DIM=1),1777)
    !       call output_delta_v_MPB_step(free_rho_superpos*pi4_inv,1778)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (solve_lpbe_only) then
        !	if we do not use the Newton solver in run_newton.f90, where ln_cosh_function is calculated, do it here:
        call evaluate_ln_cosh_function(free_hartree_superpos,delta_v_MPB,ln_cosh_function)
    end if

    !      if (use_forces.and..not.forces_mpb_force_off.and.solve_lpbe_only) then
    !	call evaluate_dnion_dalpha(free_hartree_superpos+delta_v_MPB)
    !      end if

    if (solve_lpbe_only) then
        write(info_str,*) '--------------------------'
        call localorb_info(info_str,use_unit,'(A)',OL_norm)
        write(info_str,*) '|SCLPB STATS: '
        call localorb_info(info_str,use_unit,'(A)',OL_norm)
    end if
    write(info_str,'(2X,A,I4)') '|LPBE Total steps: ', total_iter_steps
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

    time_work = mpi_wtime()-time0
    call mpi_barrier(mpi_comm_global,info)
    time_all = mpi_wtime()-time0
    call sync_real_number(time_work)
    call sync_real_number(time_all)
    if (solve_lpbe_only) then
        write(info_str,*) ' |Total time:  real work ', &
            time_work,' s, elapsed ',time_all,' s'
        call localorb_info(info_str, use_unit, "(A)", OL_norm)
        write(info_str,*) '--------------------------'
        call localorb_info(info_str,use_unit,'(A)',OL_norm)
    end if


end subroutine lpb_solver
