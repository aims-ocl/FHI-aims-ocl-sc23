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
!    Size-Modified Poisson-Boltzmann Equation in Full-Potential DFT", JCTC 12, 4052-4066 (2016).
!    DOI: 10.1021/acs.jctc.6b00435
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2016).
!  SYNOPSIS
subroutine run_newton_mpb
!  PURPOSE
!    Newton iterations to solve the MPBE. In each Newton step, the lpb_solver is called to solve an LPB like equation
!    using the MERM method.
!  SOURCE

    use lpb_solver_utilities
    use mpb_solver_utilities
    use localorb_io
    use physics, only: rho, rho_gradient, free_rho_superpos, free_rho_gradient_superpos, free_hartree_superpos !SR:laplacian of nfree,free_rho_laplace_superpos
    use synchronize_mpi_basic, only: sync_real_number
    use mpi_tasks, only: aims_stop
    implicit none

    real*8, dimension(n_full_points) :: delta_v_MPB_lpbe_temp,delta_v_MPB_mpbe_temp
    real*8, dimension(3,n_full_points) :: delta_v_MPB_gradient_lpbe_temp,delta_v_MPB_gradient_mpbe_temp
    real*8, dimension(n_full_points) :: qSPE_plus_delta_rho_lpbe_temp,qSPE_plus_delta_rho_mpbe_temp
    real*8, dimension(n_full_points) :: rho_spinless
    real*8, dimension(n_full_points) :: delta_rho_mpb
    integer :: i_iterate_mpb,total_newton_steps,i_spin,i_full_points,iter_steps
    character*400 :: info_str

    real*8, dimension(n_full_points) :: ion_dens

    logical :: dielectric_decrement_converged
    real*8 :: rmsd_dd
    integer :: i_iterate_dd, i_coords
    real*8, dimension(n_full_points) :: delta_v_MPB_old
    !with current LPBE solution, if false initialize with MPBE solution of step before
    real*8 :: charge ! this is the integrated ionic charge

    write(info_str,'(2X,A)') 'Starting Newton Solver'
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!INITIALIZATION OF VARIABLES. !!!!!!!!!!!!!!!!
    !!!!UPDATE ALL VARIABLES THAT DEPEND ON THE UPDATED ELECTRON DENSITY
    !!!!THIS IS THE DIELECTRIC FUNCTION, ION EXCLUSION FUNCTION AND PARTITION FUNCTION (IF PARTITIONED POTENTIAL IS USED FOR REGULARIZATION)
    rho_spinless = 0.0d0
    do i_spin = 1, n_spin, 1
        rho_spinless(:) = rho_spinless(:) + rho(i_spin,:)
    end do
    if (reg_method .ne. 'lpbe') then
        delta_rho_mpb = rho_spinless - pi4_inv * free_rho_superpos
    end if
    !  delta_rho_mpb = 0d0

    if (dynamic_cavity) then
        call dielecfunc_from_density(rho,dielec_func_mpb)
        call dielecfunc_gradient_from_density_gradient(dielec_func_mpb,dielec_func_gradient_mpb,&
            rho, rho_gradient)
    end if


    if (dynamic_ions) then
        call alphafunc_from_density( rho, rho_gradient, dielec_func_mpb, n_full_points)
    end if



    !   if (reg_method.eq.'vacandconst') then
    !   ! we need the solution of the LPBE with constant coefficients (eps = epsinf_mpb, kappa = kappainf_mpb)
    !     if (dynamic_cavity) then
    !        !call evaluate_rho_laplacian(rho_laplacian)
    ! !        call evaluate_partition_func_mpb(dielec_func_mpb,dielec_func_gradient_mpb,&
        ! !  	  free_rho_gradient_superpos,free_rho_laplace_superpos,&
        ! !  	  deps_drho,d2eps_drho2,&
        ! !  	  partition_func_mpb,&
        ! !  	  partition_func_mpb_gradient,partition_func_mpb_laplacian)
    !        call evaluate_partition_func_mpb(dielec_func2,dielec_func_grad2,&
        !  	  free_rho_gradient_superpos,free_rho_laplace_superpos,&
        !  	  deps2_drho,d2eps2_drho2,&
        !  	  partition_func_mpb,&
        !  	  partition_func_mpb_gradient,partition_func_mpb_laplacian)
    !     end if
    !     call solve_SPE(.false.,-delta_rho_mpb, &
        ! 	 delta_v_MPB_const(1), delta_v_MPB_const_gradient(1,1), &
        !  	.False.,0.0d0,0.0d0,.True.,.True.,.False.,.false.,&
        !  	.false.,iter_steps,.true.,delta_v_MPB_laplacian,0,.false.)
    !     !calculate the regularizion potential and gradient
    !     delta_v_MPB_vacandconst = partition_func_mpb*delta_v_MPB_const + (1d0-partition_func_mpb)*delta_v_MPB_vac
    !     do i_coords = 1, 3, 1
    !       delta_v_MPB_vacandconst_gradient(i_coords,:) = partition_func_mpb_gradient(i_coords,:)*delta_v_MPB_const +&
        ! 	partition_func_mpb*delta_v_MPB_const_gradient(i_coords,:) + (1d0-partition_func_mpb)*delta_v_MPB_vac_gradient(i_coords,:) -&
        ! 	partition_func_mpb_gradient(i_coords,:)*delta_v_MPB_vac
    !     end do
    !   end if

    !   call output_delta_v_MPB_step(free_rho_laplace_superpos,2)
    !   call output_delta_v_MPB_step(delta_v_MPB_vac,3)
    !   call output_delta_v_MPB_step(delta_v_MPB_vacandconst,4)
    !   call output_delta_v_MPB_step(partition_func_mpb,5)
    !   call output_delta_v_MPB_step(delta_v_MPB_const,6)
    !   call output_delta_v_MPB_step(dielec_func_mpb,8)
    !   call output_delta_v_MPB_step(dielec_func2,9)
    !   call output_delta_v_MPB_step(-pi4*delta_rho_mpb/epsinf_mpb + kappa_debye**2*delta_v_MPB_const,10)
    !   call output_delta_v_MPB_step(sqrt(sum(partition_func_mpb_gradient**2,DIM=1)),11)
    !

    !   stop
    !!!!!!!!!END INITIALIZATION. !!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!INITIALIZATION OF FIRST NEWTON WITH LPBE. !!!!!!!!!!!!!!!!

    if (.not.initialized_mpb_solver.or.reg_method.eq.'lpbe') then

        !if the MPBE solver is called the first time, we need to calculate the LPBE solution
        !for initialization. if we regularize with the current LPBE solution, we have to
        !recalculate it every time the Newton solver is called
        !in case of reg_method.eq.'lpbe' we thus need to save the MPBE potential
        delta_v_MPB_mpbe_temp = delta_v_MPB
        delta_v_MPB_gradient_mpbe_temp = delta_v_MPB_gradient
        qSPE_plus_delta_rho_mpbe_temp = qSPE_plus_delta_rho
        !we need the LPBE solution for initialization/regularization
        write(info_str,'(2X,A)') '|Solving LPBE for Initialization/Regularization.'
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

        call lpb_solver(.false.,0d0,0d0)
        !we have just calculated LPBE solution
        !save all functions we need to add them again in end of routine
        delta_v_MPB_lpbe_temp = delta_v_MPB
        delta_v_MPB_gradient_lpbe_temp = delta_v_MPB_gradient
        qSPE_plus_delta_rho_lpbe_temp = qSPE_plus_delta_rho
        !if we regularize with free pot, we want delta_v_MPB to be the LPBE solution for
        !initialization, if we regularize with the LPBE potential, we recover the MPBE solution of the previous step
        !       if (reg_method.eq.'lpbe') then
        ! 	delta_v_MPB = delta_v_MPB_mpbe_temp
        ! 	delta_v_MPB_gradient = delta_v_MPB_gradient_mpbe_temp
        ! 	rho_mpb = rho_mpb_mpbe_temp
        !       end if
        !       if (reg_method.eq.'vacandconst') then
        ! 	!we have to calculate the actual delta v = v - vreg, with vreg = w*vconst + (1-w)*vvac
        ! 	delta_v_MPB_lpbe_temp = delta_v_MPB_lpbe_temp - delta_v_MPB_vacandconst
        ! 	delta_v_MPB_gradient_lpbe_temp = delta_v_MPB_gradient_lpbe_temp - delta_v_MPB_vacandconst_gradient
        !       end if
        write(info_str,'(2X,A)') '|Newton Solver Initialization Finished '
        call localorb_info(info_str,use_unit,'(A)',OL_norm)
    end if
    if (reg_method.eq.'lpbe') then
        !update total hartree potential
        pot_lpbe = free_hartree_superpos+delta_v_MPB_lpbe_temp
        pot_lpbe_gradient = v_free_gradient_public+delta_v_MPB_gradient_lpbe_temp
    end if

    !!!!!!!!!END INITIALIZATION. !!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !   call output_delta_v_MPB_step(delta_v_MPB_lpbe_temp,17)
    !   call output_delta_v_MPB_step(sqrt(sum(delta_v_MPB_gradient_lpbe_temp**2,DIM=1)),18)
    !!!!!!!!!NEWTON LOOP. !!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    dielectric_decrement_converged = .false.
    i_iterate_dd = 0

    do while (.not.(dielectric_decrement_converged))
        !if we want to correct for the dielectric decrement we run another scf loop around the Newton method (eps_s then depends on ionic concentration)
        !if not just to the Newton


        if (correct_dielectric_decrement .and. dielectric_decrement_method == 'SCF') then !.and.&
            delta_v_MPB_old = delta_v_MPB

            i_iterate_dd = i_iterate_dd + 1

            write(info_str,'(2X,A,I4)') '|SC-Dielectric-Decrement, step:',i_iterate_dd
            call localorb_info(info_str,use_unit,'(A)',OL_norm)

            !       (rmsd_newton.lt.start_dielec_decr_corr_acc).and.rmsd_newton.ne.0d0) then
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
            !       if (initialized_mpb_solver) then
            ! 	  call output_delta_v_MPB_step(dielec_func_mpb,i_iterate_dd)
            ! 	  call output_delta_v_MPB_step(sqrt(sum(dielec_func_gradient_mpb**2,DIM=1)),i_iterate_dd+200)
            !       endif
        else
            dielectric_decrement_converged = .true.
        end if

        !now heres the actual Newton loop
        newton_converged = .false.
        i_iterate_mpb = 0
        rmsd_newton = 0d0

        !Newton solver. In each iteration step we solve a LPB kind of equation. The solving is performed by the MERM as carried out by lpb_solver.f90.
        do while (.not.(newton_converged))

            i_iterate_mpb = i_iterate_mpb + 1

            !at first step we will initialize with delta_mpbe = 0
            !all following steps we will initialize with the delta_mpbe from step before
            !in the case of MERM_in_SPE_solver, we want to take all the splines and parameters from the LPBE solution, so we
            !have already initialized everything
            if (i_iterate_mpb==1.and.reg_method.eq.'lpbe') initialized_lpb_solver = .False.

            write(info_str,'(2X,A,I4)') '|Newton Solver of MPBE, step:',i_iterate_mpb
            call localorb_info(info_str,use_unit,'(A)',OL_norm)

            if (i_iterate_mpb==1) then
                !initializations of the Newton solver
                if (.not.initialized_mpb_solver) then
                    !we run the Newton solver the first time
                    !we want to initialize with the LPBE solution calculated above
                    if (reg_method.eq.'lpbe') then
                        !regularization with LPBE, residue is set to zero
                        delta_mpbe = 0d0 !delta_v_MPB_lpbe_temp
                        delta_mpbe_gradient = 0d0 !delta_v_MPB_gradient_lpbe_temp
                        ! 	  else if (reg_method.eq.'vacandconst') then
                        ! 	    !regularization with free atom potential, residue is the lpbe residue
                        ! 	    delta_mpbe = delta_v_MPB_lpbe_temp - delta_v_MPB_vacandconst
                        ! 	    delta_mpbe_gradient = delta_v_MPB_gradient_lpbe_temp
                    else
                        delta_mpbe = delta_v_MPB_lpbe_temp
                        delta_mpbe_gradient = delta_v_MPB_gradient_lpbe_temp
                    end if
                else
                    ! 	  if (reg_method.eq.'vacandconst') then
                    ! 	    !delta_v_MPB = delta v = v - vfree
                    ! 	    !delta_mpbe = v - vreg, where vref = w*vconst + (1-w)*vvac
                    ! 	    continue
                    ! 	  else
                    !we set the potential to the MPBE potential from previous SCF step
                    delta_mpbe = delta_v_MPB
                    delta_mpbe_gradient = delta_v_MPB_gradient
                    ! 	  end if
                end if
            end if


            delta_mpbe_old = delta_mpbe


            if (reg_method.eq.'lpbe') then
                !update h_function
                !in this case, we need different h and f functions
                call evaluate_h_and_f_function(pot_lpbe,pot_lpbe_gradient,delta_mpbe,&
                    delta_mpbe_gradient, h_function,&
                    f_function,ion_conc,ion_conc_deriv)
                !evaluate right hand side q_times_epsinf_mpb
                call evaluate_iterative_mpbe_function(dielec_func_mpb,&
                    dielec_func_gradient_mpb, pot_lpbe,delta_mpbe,delta_mpbe_gradient,h_function,f_function,&
                    q_times_epsinf_mpb,v_free_gradient_public,free_rho_superpos,delta_rho_mpb)
                !       else if (reg_method.eq.'vacandconst') then
                ! 	!update h_function
                ! 	if (correct_dielectric_decrement) then
                ! 	  call evaluate_h_and_f_function(free_hartree_superpos,delta_mpbe + delta_v_MPB_vacandconst,h_function,&
                    ! 	    f_function,ion_conc,ion_conc_deriv,.true.)
                ! 	else
                ! 	  call evaluate_h_and_f_function(free_hartree_superpos,delta_mpbe + delta_v_MPB_vacandconst,h_function,&
                    ! 	    f_function,ion_conc,ion_conc_deriv,.false.)
                ! 	end if
                ! 	call evaluate_iterative_mpbe_function(dielec_func_mpb,&
                    ! 	  dielec_func_gradient_mpb, free_hartree_superpos,delta_mpbe,delta_mpbe_gradient,h_function,f_function,&
                    ! 	  q_times_epsinf_mpb,v_free_gradient_public,free_rho_superpos,delta_rho_mpb)
            else !(reg_method.eq.'lpbe')
                !update h_function
                call evaluate_h_and_f_function(free_hartree_superpos, v_free_gradient_public, delta_mpbe,&
                    delta_mpbe_gradient,h_function,&
                    f_function,ion_conc,ion_conc_deriv)
                !evaluate right hand side q_times_epsinf_mpb
                call evaluate_iterative_mpbe_function(dielec_func_mpb,&
                    dielec_func_gradient_mpb, free_hartree_superpos,delta_mpbe,delta_mpbe_gradient,h_function,f_function,&
                    q_times_epsinf_mpb,v_free_gradient_public,free_rho_superpos,delta_rho_mpb)
            end if
            !    call output_delta_v_MPB_step(f_function,555)
            !    call output_delta_v_MPB_step(h_function,666)

            !copy working arrays in newton methods to working arrays in lpbe solver (in general this should have not effect, just for clarity)
            delta_v_MPB = delta_mpbe
            delta_v_MPB_gradient = delta_mpbe_gradient

            !       call output_delta_v_MPB_step(q_times_epsinf_mpb,i_iterate_mpb+550)

            !solve LPBE with Newton parameters
            call lpb_solver(.true.,h_function,q_times_epsinf_mpb)

            !copy the calculated arrays to the name used in this loop
            delta_mpbe = delta_v_MPB
            delta_mpbe_gradient = delta_v_MPB_gradient
            qSPE_plus_delta_rho_mpbe = qSPE_plus_delta_rho

            call evaluate_rmsd_newton(delta_mpbe(1), delta_mpbe_old(1), rmsd_newton, tau_newton, newton_converged)

            !IMPORTANT OUTPUTS
            !ionic charge density
            !     call output_delta_v_MPB_step(-pi4_inv*f_function,8000)
            !     call output_delta_v_MPB_step(dielec_func_mpb,9000)
            !     call output_delta_v_MPB_step(delta_v_MPB+free_hartree_superpos,9200)


            write(info_str,'(2X,A,ES14.3)') '|Current RMSD Newton: ', rmsd_newton
            call localorb_info(info_str,use_unit,'(A)',OL_norm)

        end do !end newton



        total_newton_steps = i_iterate_mpb

        write(info_str,'(2X,A,I4)') '|Newton converged in Total steps: ', total_newton_steps
        call localorb_info(info_str,use_unit,'(A)',OL_norm)



        if (correct_dielectric_decrement .and. dielectric_decrement_method == 'SCF') then
            call evaluate_rmsd_newton(delta_v_MPB,delta_v_MPB_old,rmsd_dd,tau_newton,dielectric_decrement_converged)
            write(info_str,'(2X,A,ES14.3)') '|Current RMSD SC-Dielectric-Decrement: ', rmsd_dd
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
        end if
    end do !end dielectric decrement

    !!!!!!!!!END NEWTON LOOP. !!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    call output_delta_v_step(dielec_func_mpb,100)
    !call output_delta_v_step(sqrt(sum(dielec_func_gradient_mpb**2,DIM=1)),200)


    !copy arrays to eventual location
    !delta_v_MPB = v - v_free = (pot_lpbe - v_free) + delta_mpbe
    if (reg_method.eq.'lpbe') then
        !in this case we have to add the LPBE solution, which we saved before
        delta_v_MPB = delta_v_MPB_lpbe_temp + delta_mpbe
        delta_v_MPB_gradient = delta_v_MPB_gradient_lpbe_temp + delta_mpbe_gradient
        qSPE_plus_delta_rho = qSPE_plus_delta_rho_lpbe_temp + qSPE_plus_delta_rho_mpbe
        call aims_stop('|Regularization with LPBE is not working yet for energies')
        !else: delta_v_MPB and delta_v_MPB_gradient are already the solution and gradient of the MPBE
        !   else if (reg_method.eq.'vacandconst') then
        !     !add other delta potentials
        !     !to obtain our delta v = v - vfree which we usually use in aims
        !     delta_v_MPB = delta_v_MPB + delta_v_MPB_vacandconst
        !     do i_coords = 1, 3, 1
        !       delta_v_MPB_gradient(i_coords,:) = delta_v_MPB_gradient(i_coords,:) +  delta_v_MPB_vacandconst_gradient(i_coords,:)
        !     end do
    end if



    !evaluate ln_cosh function for free energy
    call evaluate_ln_cosh_function(free_hartree_superpos,delta_v_MPB,ln_cosh_function)
    call evaluate_cosh_cosh_function(free_hartree_superpos,delta_v_MPB,cosh_cosh_function)

    !call evaluate_rdf(4.d-1,200,4.d-1,sum(ion_conc,dim=1),"rdf_ion_conc.dat")
    !  call output_delta_v_step(pi4_inv/2d0*kappainf_mpb**2*&
    !	(alpha_func_mpb*exp(-z_mpb*(free_hartree_superpos+delta_v_MPB)/kBT_mpb)/(1d0-phi_zero_mpb+phi_zero_mpb*alpha_func_mpb*cosh(z_mpb*(free_hartree_superpos+delta_v_MPB)/kBT_mpb))+&
        !	alpha_func_anion_mpb*exp(z_mpb*(free_hartree_superpos+delta_v_MPB)/kBT_mpb)/(1d0-phi_zero_mpb+phi_zero_mpb*alpha_func_mpb*cosh(z_mpb*(free_hartree_superpos+delta_v_MPB)/kBT_mpb))),1111)
    !  call output_delta_v_step(c_mpb *&
        !alpha_func_mpb*exp(-z_mpb*(free_hartree_superpos+delta_v_MPB)/kBT_mpb)/(1d0-phi_zero_mpb+phi_zero_mpb*alpha_func_mpb*cosh(z_mpb*(free_hartree_superpos+delta_v_MPB)/kBT_mpb)),2222)
    !  call output_delta_v_step(c_mpb*&
        !	alpha_func_anion_mpb*exp(z_mpb*(free_hartree_superpos+delta_v_MPB)/kBT_mpb)/(1d0-phi_zero_mpb+phi_zero_mpb*alpha_func_mpb*cosh(z_mpb*(free_hartree_superpos+delta_v_MPB)/kBT_mpb)),3333)
    !    if (use_forces.and..not.forces_mpb_force_off) then
    !    call evaluate_dnion_dalpha(free_hartree_superpos+delta_v_MPB)
    !  end if

    !FEM OUTPUTS
    ! call output_delta_v_MPB_step(delta_mpbe,900)
    !   call output_delta_v_MPB_step(-pi4_inv*f_function,200)
    ! call output_delta_v_step(alpha_func_mpb,100)
    ! if (use_separate_alpha) then
    !   call output_delta_v_MPB_step(alpha_func_anion_mpb,101)
    ! end if
    ! call output_delta_v_step(dielec_func_mpb,200)
    ! call output_delta_v_MPB_step(sum(rho,DIM=1),300)
    ! call output_delta_v_MPB_step(delta_v_MPB,300)
    ! call output_delta_v_MPB_step(free_hartree_superpos,400)
    ! call output_delta_v_MPB_step(sqrt(sum(delta_v_MPB_gradient**2,DIM=1)),800)
    ! stop
    !  write(use_unit,*) 'outputting'
    !  call output_delta_v_step(f_function,1)
    !  call output_delta_v_step(dielec_func_mpb,2)
    initialized_mpb_solver = .True.

    charge = 0d0
    do i_full_points = 1, n_full_points, 1
        charge = charge +&
        (ion_conc(1,i_full_points)-ion_conc(2,i_full_points))*partition_tab(i_full_points)
    end do
    call sync_real_number(charge)
    write(info_str,'(X)')
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A,F10.8)') '|Integrated ionic charge = ',-charge
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

end subroutine run_newton_mpb
