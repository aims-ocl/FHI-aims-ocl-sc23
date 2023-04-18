!****s* FHI-aims/mpb_solver_utilities
!  NAME
!   mpb_solver_utilities
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

module mpb_solver_utilities

!  PURPOSE
!
!  Module mpb_solver_utilities contains variables & routines specific to the size-
!  modification of the PBE and the nonlinear character, as it is also the Newton method
!  infrastructure and all free energy formulas
!
!  VARIABLES
!
!  //GENERAL CONVERGENCE//
!  newton_converged : true if Newton method is converged
!  initialized_mpb_solver : false if the Newton method has not been called during the SCF, true if already called at least once
!  rmsd_newton : root mean square deviation of previous delta v with current one
!  tau_newton : accuracy criterium of delta v that has to be fulfilled
!
!  //SYSTEM SPECIFICATIONS//
!  freeen_MPBE_nonelstat_solv_energy : nonelectrostatic free energy part of the total energy
!
!  //SPECIFIC CALCULATIONS//
!  KS_mpb_corr1 : KS Hamiltonian correction due to the electron density dependence of the dielectric function eps
!  KS_mpb_corr2 : KS Hamiltonian correction due to the electron density dependence of the exclusion function alpha
!
!  //TECHNICAL VARIABLES//
!  reg_method : lpbe or vacandconst  or none. "none" is the only proved working choice
!	and means that v^reg = vfree. vacandconst: regularize with v^reg = w*vvac+(1-w)*v^SPE, where
!	vvac is the vacuum Hartree potential and v^SPE the far field solution. By that one can perfectly localize
!	the delta v potential in the transition region but we observed worse convergence with lmax. to try comment in all relevant
!  	parts where "vacandconst" appears in mpb_solver_utilities, run_newton_mpb and sum_up_whole_potential and
!	initialize_scf, lpb_solver_utilities. 'lpbe' is not working, yet, but can be in principle realily implemented. it
! 	is however not so neat since we regularize with a potential that already has a multipole error
!  ffunction_warning_nan, ffunction_warning_nan_done : for very large electrostatic potential we get numerical problems
!       in the evaluation of the f and h functions, so we use some limit approximations, this warning reminds the user of that
!
!  //VARIABLES ON FHI-AIMS GRID//
!  pot_lpbe : LPBE electrostatic potential for reg_method = 'lpbe'
!  pot_lpbe_gradient : corresponding gradient
!  delta_mpbe, delta_mpbe_gradient : delta v potential. local variable which in the end is the same as the global delta_v_MPB
!  delta_mpbe_old : delta_mpbe from previous iteration
!  h_function, f_function : h and f functions, where f = -pi4 * nion^MPB and h_function = h^2 from paper
!  q_times_epsinf_mpb : q * epsinf_mpb, where q is the q from the Newton method
!  qSPE_plus_delta_rho_mpbe : copy of qSPE_plus_delta_rho withing Newton method (only for reg_method = 'lpbe')
!  ln_cosh_function : eq. (18) in the paper, "ln(...)" part without prefactor
!  cosh_cosh_function : eq. (17) in the paper, fraction that contains two cosh function without prefactor
!  cut_ln_cosh : value of the argument of the ln(a-cosh(x)) argument x from which on function is always  x+ln(0.5*phi/(1-phi))
!  cut_sinh_cosh : value at which argument of sinh(x)/(1-phi+phi*cosh(x)) arg x from which on which function is always. 1/(phi*alpha)
!  cut_cosh_cosh : value of the argument of the cosh(x)/(1-phi+phi*cosh(x) arg x from which on function is always 1/(phi*alpha)
!  USES

    use dimensions
    use constants, only: pi4, pi4_inv
    use lpb_solver_utilities, only: z_mpb, kBT_mpb,phi_zero_mpb,epsinf_mpb,kappainf_mpb,&
        kappa_debye,c_mpb, alpha_func_mpb, solve_lpbe_only, solve_debye_only,&
        v_free_gradient_public, use_separate_alpha,&
        alpha_func_anion_mpb,a_mpb, solve_pbe_only, dalpha_drho, deps_dcion,&
        correct_dielectric_decrement, dielectric_decrement_method,dcion_dalpha,&
        dcion_dgradv, atomic_MERM
    !comment in for reg_method=='vacandconst':
    !,pot_helmh_vac_gradient, pot_helmh_const_gradient, pot_helmh_vac, pot_helmh_const
    use localorb_io
    use physics, only: partition_tab



    implicit none

    logical :: newton_converged
    character*15 :: reg_method
    logical :: initialized_mpb_solver = .False.
    real*8 :: ffunction_warning_nan = 0d0
    logical :: ffunction_warning_nan_done = .False.
    ! 	  logical :: ffunction_warning_nan_special
    logical :: KS_mpb_corr1
    logical :: KS_mpb_corr2

    real*8 :: rmsd_newton
    real*8 :: tau_newton = 1d-10

    real*8 :: freeen_MPBE_nonelstat_solv_energy !nonelstat. solvation energy according to MPBE free energy expression

    real*8, dimension(:), allocatable :: pot_lpbe
    real*8, dimension(:,:), allocatable :: pot_lpbe_gradient
    real*8, dimension(:), allocatable :: delta_mpbe
    real*8, dimension(:,:), allocatable :: delta_mpbe_gradient
    real*8, dimension(:), allocatable :: q_times_epsinf_mpb
    real*8, dimension(:), allocatable :: delta_mpbe_old
    real*8, dimension(:), allocatable :: h_function
    real*8, dimension(:), allocatable :: f_function
    real*8, dimension(:), allocatable :: qSPE_plus_delta_rho_mpbe
    real*8, dimension(:), allocatable :: ln_cosh_function
    real*8, dimension(:), allocatable :: cosh_cosh_function
    real*8, dimension(:), allocatable :: rho_laplacian

    !reg_method=='vacandconst'
    ! 	  real*8, dimension(:), allocatable :: pot_helmh_vacandconst
    ! 	  real*8, dimension(:,:), allocatable :: pot_helmh_vacandconst_gradient
    ! 	  real*8, dimension(:), allocatable :: partition_func_mpb
    ! 	  real*8, dimension(:,:), allocatable :: partition_func_mpb_gradient
    ! 	  real*8, dimension(:), allocatable :: partition_func_mpb_laplacian

    real*8 :: cut_ln_cosh = 100
    real*8 :: cut_sinh_cosh = 100
    real*8 :: cut_cosh_cosh = 100


contains


    !****s* FHI-aims/allocate_mpbsolver
    !  NAME
    !  allocate_mpbsolver
    !  SYNOPSIS
    subroutine allocate_mpbsolver()
    !  PURPOSE
    !	allocate all variables of mpb_solver_utilities

        use runtime_choices
        use mpi_tasks
        implicit none

        integer:: info

        !reg_method=='vacandconst'
        ! 	  if (.not.allocated(pot_helmh_vacandconst_gradient)) then
        ! 	    allocate(pot_helmh_vacandconst_gradient(3,n_full_points),stat=info)
        ! 	    call check_allocation(info, 'pot_helmh_vacandconst_gradient                  ')
        ! 	    pot_helmh_vacandconst_gradient(:,:) = 0d0
        ! 	  end if
        !
        ! 	  if (.not.allocated(pot_helmh_vacandconst)) then
        ! 	    allocate(pot_helmh_vacandconst(n_full_points),stat=info)
        ! 	    call check_allocation(info, 'pot_helmh_vacandconst                  ')
        ! 	    pot_helmh_vacandconst(:) = 0d0
        ! 	  end if


        ! 	  if (.not.allocated(partition_func_mpb_gradient)) then
        ! 	    allocate(partition_func_mpb_gradient(3,n_full_points),stat=info)
        ! 	    call check_allocation(info, 'partition_func_mpb_gradient                  ')
        ! 	    partition_func_mpb_gradient(:,:) = 0d0
        ! 	  end if
        !
        ! 	  if (.not.allocated(partition_func_mpb_laplacian)) then
        ! 	    allocate(partition_func_mpb_laplacian(n_full_points),stat=info)
        ! 	    call check_allocation(info, 'partition_func_mpb_laplacian                  ')
        ! 	    partition_func_mpb_laplacian(:) = 0d0
        ! 	  end if
        !
        ! 	  if (.not.allocated(partition_func_mpb)) then
        ! 	    allocate(partition_func_mpb(n_full_points),stat=info)
        ! 	    call check_allocation(info, 'partition_func_mpb                  ')
        ! 	    partition_func_mpb(:) = 0d0
        ! 	  end if

        if (.not.allocated(qSPE_plus_delta_rho_mpbe)) then
            allocate(qSPE_plus_delta_rho_mpbe(n_full_points),stat=info)
            call check_allocation(info, 'qSPE_plus_delta_rho_mpbe                  ')
            qSPE_plus_delta_rho_mpbe(:) = 0d0
        end if

        if (.not.allocated(delta_mpbe)) then
            allocate(delta_mpbe(n_full_points),stat=info)
            call check_allocation(info, 'delta_mpbe                  ')
            delta_mpbe(:) = 0d0
        end if

        if (.not.allocated(delta_mpbe_gradient)) then
            allocate(delta_mpbe_gradient(3,n_full_points),stat=info)
            call check_allocation(info, 'delta_mpbe_gradient                  ')
            delta_mpbe_gradient(:,:) = 0d0
        end if

        if (.not.allocated(pot_lpbe)) then
            allocate(pot_lpbe(n_full_points),stat=info)
            call check_allocation(info, 'pot_lpbe                  ')
            pot_lpbe(:) = 0d0
        end if

        if (.not.allocated(pot_lpbe_gradient)) then
            allocate(pot_lpbe_gradient(3,n_full_points),stat=info)
            call check_allocation(info, 'pot_lpbe_gradient                  ')
            pot_lpbe_gradient(:,:) = 0d0
        end if

        if (.not.allocated(q_times_epsinf_mpb)) then
            allocate(q_times_epsinf_mpb(n_full_points),stat=info)
            call check_allocation(info, 'q_times_epsinf_mpb                  ')
            q_times_epsinf_mpb(:) = 0d0
        end if

        if (.not.allocated(delta_mpbe_old)) then
            allocate(delta_mpbe_old(n_full_points),stat=info)
            call check_allocation(info, 'delta_mpbe_old                  ')
            delta_mpbe_old(:) = 0d0
        end if

        if (.not.allocated(h_function)) then
            allocate(h_function(n_full_points),stat=info)
            call check_allocation(info, 'h_function                  ')
            h_function(:) = 0d0
        end if

        if (.not.allocated(f_function)) then
            allocate(f_function(n_full_points),stat=info)
            call check_allocation(info, 'f_function                  ')
            f_function(:) = 0d0
        end if

        if (.not.allocated(ln_cosh_function)) then
            allocate(ln_cosh_function(n_full_points),stat=info)
            call check_allocation(info, 'ln_cosh_function                  ')
            ln_cosh_function(:) = 0d0
        end if

        if (.not.allocated(cosh_cosh_function)) then
            allocate(cosh_cosh_function(n_full_points),stat=info)
            call check_allocation(info, 'cosh_cosh_function                  ')
            cosh_cosh_function(:) = 0d0
        end if

    end subroutine allocate_mpbsolver

    subroutine allocate_mpbsolver_dummy()
    !  PURPOSE
    !	allocate all variables of mpb_solver_utilities with smallest array sizes,
    !       since they are sometimes passed even in the case we do not need the MPB solver
    !       we need to still allocate them

        use runtime_choices
        use mpi_tasks
        implicit none

        integer:: info

        !reg_method=='vacandconst'
        ! 	  if (.not.allocated(pot_helmh_vacandconst_gradient)) then
        ! 	    allocate(pot_helmh_vacandconst_gradient(3,n_full_points),stat=info)
        ! 	    call check_allocation(info, 'pot_helmh_vacandconst_gradient                  ')
        ! 	    pot_helmh_vacandconst_gradient(:,:) = 0d0
        ! 	  end if
        !
        ! 	  if (.not.allocated(pot_helmh_vacandconst)) then
        ! 	    allocate(pot_helmh_vacandconst(n_full_points),stat=info)
        ! 	    call check_allocation(info, 'pot_helmh_vacandconst                  ')
        ! 	    pot_helmh_vacandconst(:) = 0d0
        ! 	  end if


        ! 	  if (.not.allocated(partition_func_mpb_gradient)) then
        ! 	    allocate(partition_func_mpb_gradient(3,n_full_points),stat=info)
        ! 	    call check_allocation(info, 'partition_func_mpb_gradient                  ')
        ! 	    partition_func_mpb_gradient(:,:) = 0d0
        ! 	  end if
        !
        ! 	  if (.not.allocated(partition_func_mpb_laplacian)) then
        ! 	    allocate(partition_func_mpb_laplacian(n_full_points),stat=info)
        ! 	    call check_allocation(info, 'partition_func_mpb_laplacian                  ')
        ! 	    partition_func_mpb_laplacian(:) = 0d0
        ! 	  end if
        !
        ! 	  if (.not.allocated(partition_func_mpb)) then
        ! 	    allocate(partition_func_mpb(n_full_points),stat=info)
        ! 	    call check_allocation(info, 'partition_func_mpb                  ')
        ! 	    partition_func_mpb(:) = 0d0
        ! 	  end if

        if (.not.allocated(qSPE_plus_delta_rho_mpbe)) then
            allocate(qSPE_plus_delta_rho_mpbe(1),stat=info)
            call check_allocation(info, 'qSPE_plus_delta_rho_mpbe                  ')
            qSPE_plus_delta_rho_mpbe(:) = 0d0
        end if

        if (.not.allocated(delta_mpbe)) then
            allocate(delta_mpbe(1),stat=info)
            call check_allocation(info, 'delta_mpbe                  ')
            delta_mpbe(:) = 0d0
        end if

        if (.not.allocated(delta_mpbe_gradient)) then
            allocate(delta_mpbe_gradient(1,1),stat=info)
            call check_allocation(info, 'delta_mpbe_gradient                  ')
            delta_mpbe_gradient(:,:) = 0d0
        end if

        if (.not.allocated(pot_lpbe)) then
            allocate(pot_lpbe(1),stat=info)
            call check_allocation(info, 'pot_lpbe                  ')
            pot_lpbe(:) = 0d0
        end if

        if (.not.allocated(pot_lpbe_gradient)) then
            allocate(pot_lpbe_gradient(3,n_full_points),stat=info)
            call check_allocation(info, 'pot_lpbe_gradient                  ')
            pot_lpbe_gradient(:,:) = 0d0
        end if

        if (.not.allocated(q_times_epsinf_mpb)) then
            allocate(q_times_epsinf_mpb(1),stat=info)
            call check_allocation(info, 'q_times_epsinf_mpb                  ')
            q_times_epsinf_mpb(:) = 0d0
        end if

        if (.not.allocated(delta_mpbe_old)) then
            allocate(delta_mpbe_old(1),stat=info)
            call check_allocation(info, 'delta_mpbe_old                  ')
            delta_mpbe_old(:) = 0d0
        end if

        if (.not.allocated(h_function)) then
            allocate(h_function(1),stat=info)
            call check_allocation(info, 'h_function                  ')
            h_function(:) = 0d0
        end if

        if (.not.allocated(f_function)) then
            allocate(f_function(1),stat=info)
            call check_allocation(info, 'f_function                  ')
            f_function(:) = 0d0
        end if

        if (.not.allocated(cosh_cosh_function)) then
            allocate(cosh_cosh_function(1),stat=info)
            call check_allocation(info, 'cosh_cosh_function                  ')
            cosh_cosh_function(:) = 0d0
        end if

    end subroutine allocate_mpbsolver_dummy

    !****s* FHI-aims/deallocate_mpbsolver
    !  NAME
    !  deallocate_mpbsolver
    !  SYNOPSIS
    subroutine deallocate_mpbsolver()
    !  PURPOSE
    !	deallocate all variables of mpb_solver_utilities

        use runtime_choices

        implicit none

        integer:: info

        ! 	  if (allocated(partition_func_mpb_gradient)) then
        ! 	    deallocate(partition_func_mpb_gradient)
        ! 	  end if
        !
        ! 	  if (allocated(partition_func_mpb_laplacian)) then
        ! 	    deallocate(partition_func_mpb_laplacian)
        ! 	  end if
        !
        ! 	  if (allocated(partition_func_mpb)) then
        ! 	    deallocate(partition_func_mpb)
        ! 	  end if


        if (allocated(qSPE_plus_delta_rho_mpbe)) then
            deallocate(qSPE_plus_delta_rho_mpbe)
        end if

        if (allocated(delta_mpbe)) then
            deallocate(delta_mpbe)
        end if

        if (allocated(delta_mpbe_gradient)) then
            deallocate(delta_mpbe_gradient)
        end if

        if (allocated(delta_mpbe_old)) then
            deallocate(delta_mpbe_old)
        end if

        if (allocated(pot_lpbe)) then
            deallocate(pot_lpbe)
        end if

        if (allocated(pot_lpbe_gradient)) then
            deallocate(pot_lpbe_gradient)
        end if

        if (allocated(q_times_epsinf_mpb)) then
            deallocate(q_times_epsinf_mpb)
        end if


        if (allocated(h_function)) then
            deallocate(h_function)
        end if

        if (allocated(f_function)) then
            deallocate(f_function)
        end if


        if (allocated(ln_cosh_function)) then
            deallocate(ln_cosh_function)
        end if

        if (allocated(cosh_cosh_function)) then
            deallocate(cosh_cosh_function)
        end if

    end subroutine deallocate_mpbsolver

    !****s* FHI-aims/evaluate_cosh_cosh_function
    !  NAME
    !  evaluate_cosh_cosh_function
    !  SYNOPSIS
    subroutine evaluate_cosh_cosh_function(pot_regul,delta_mpbe_calc,cosh_cosh_function_calc)
    !  PURPOSE
    ! evaluate cosh(beta z v) / (1-phi0+phi0*alpha[n_el]*cosh(beta z v)) + appropriate limits at high v
        use mpi_utilities

        integer :: i_full_points
        real*8, dimension(n_full_points),intent(in) :: pot_regul
        real*8, dimension(n_full_points),intent(in) :: delta_mpbe_calc
        real*8, dimension(n_full_points), intent(out) :: cosh_cosh_function_calc
        real*8 :: current_potential

        character*600 :: info_str
        real*8 :: cosh_func
        real*8 :: abs_dalpha_drho

        !the function is only evaluated at points, where |dalpha_drho| .gt.0.0!!!! because
        !these are the points, where we actually need it


        !SR: we can use a limiting expression here also. I however commented out, because usually
        ! we do not need it. if someone needs it please check this carefully

        cosh_cosh_function_calc = 0d0

        do i_full_points = 1, n_full_points, 1
            abs_dalpha_drho = abs(dalpha_drho(i_full_points))
            if (abs_dalpha_drho.gt.0d0.and.partition_tab(i_full_points).gt.0d0) then
                ! 	      if (alpha_func_mpb(i_full_points).lt.1d-10) then
                ! 		write(info_str,'(2X,2A)') '|WARNING: Calculating the MPBE free energy cosh-cosh-function as limiting expression,',&
                    ! 		'although the limit may not be reached for this particular point, since the alpha-function is very small.'
                ! 		call localorb_info(info_str,use_unit,'(A)',OL_norm)
                ! 	      end if
                current_potential = pot_regul(i_full_points) + delta_mpbe_calc(i_full_points)
                ! 	      if (abs(z_mpb*current_potential/kBT_mpb).gt.cut_cosh_cosh.and..not.solve_pbe_only) then
                ! 		if (current_potential.lt.0d0) then
                ! 		  cosh_cosh_function_calc(i_full_points) =&
                    ! 		    (1d0/(phi_zero_mpb*alpha_func_mpb(i_full_points)))
                ! 		else if (current_potential.gt.0d0) then
                ! 		  cosh_cosh_function_calc(i_full_points) =&
                    ! 		    (1d0/(phi_zero_mpb*alpha_func_mpb(i_full_points)))
                ! 		end if
                ! 	      else
                cosh_func = cosh(z_mpb*current_potential/kBT_mpb)
                if (.not.solve_pbe_only) then
                    cosh_cosh_function_calc(i_full_points) = &
                        cosh_func/(1-phi_zero_mpb+phi_zero_mpb*alpha_func_mpb(i_full_points)*&
                        cosh_func)
                else
                    cosh_cosh_function_calc(i_full_points) = cosh_func
                end if
                ! 	      end if
            end if
        end do


    end subroutine evaluate_cosh_cosh_function

    !****s* FHI-aims/evaluate_ln_cosh_function
    !  NAME
    !  evaluate_ln_cosh_function
    !  SYNOPSIS
    subroutine evaluate_ln_cosh_function(pot_regul,delta_mpbe_calc,ln_cosh_function_calc)
    !  PURPOSE
    ! evaluate ln(1+phi0/(1-phi0)*alpha[n_el]*cosh(beta z v)) + appropriate limits at high v

        use mpi_utilities

        integer :: i_full_points
        real*8, dimension(n_full_points),intent(in) :: pot_regul
        real*8, dimension(n_full_points),intent(in) :: delta_mpbe_calc
        real*8, dimension(n_full_points), intent(out) :: ln_cosh_function_calc
        real*8 :: current_potential
        real*8 :: part_cat, part_an,exp_cat,exp_an,Afunc

        character*600 :: info_str

        !the function is only evaluated at points, where alpha_func_mpb.gt.0.0!!!! because
        !these are the points, where we actually need it

        ln_cosh_function_calc = 0d0

        do i_full_points = 1, n_full_points, 1
            if (alpha_func_mpb(i_full_points).gt.0d0.and.partition_tab(i_full_points).gt.0d0.and..not.use_separate_alpha) then
                if (alpha_func_mpb(i_full_points).lt.1d-12) then
                    write(info_str,'(2X,2A)') '|WARNING: Calculating the MPBE free energy ln-cosh-function as limiting expression,',&
                        'although the limit may not be reached for this particular point, since the alpha-function is very small.'
                    call localorb_info(info_str,use_unit,'(A)',OL_norm)
                end if
                current_potential = pot_regul(i_full_points) + delta_mpbe_calc(i_full_points)
                if (abs(z_mpb*current_potential/kBT_mpb).gt.cut_ln_cosh.and..not.solve_pbe_only) then
                    if (current_potential.lt.0d0) then
                        ln_cosh_function_calc(i_full_points) =&
                            (-z_mpb*current_potential/kBT_mpb+log(0.5d0*alpha_func_mpb(i_full_points)*phi_zero_mpb/(1d0-phi_zero_mpb)))
                    else if (current_potential.gt.0d0) then
                        ln_cosh_function_calc(i_full_points) =&
                            (z_mpb*current_potential/kBT_mpb+log(0.5d0*alpha_func_mpb(i_full_points)*phi_zero_mpb/(1d0-phi_zero_mpb)))
                    end if
                else
                    if (.not.solve_lpbe_only) then
                        if (.not.solve_pbe_only) then
                            ln_cosh_function_calc(i_full_points) = &
                                log(1d0+phi_zero_mpb/(1-phi_zero_mpb)*alpha_func_mpb(i_full_points)*cosh(z_mpb*current_potential/kBT_mpb))
                        else !we will use 1/a^3*ln_cosh_function which is for a->0:
                            ln_cosh_function_calc(i_full_points) = &
                                2d0*c_mpb*alpha_func_mpb(i_full_points)*cosh(z_mpb*current_potential/kBT_mpb)
                        end if
                    else
                        if (.not.solve_debye_only) then
                            ln_cosh_function_calc(i_full_points) = &
                                log(1d0+phi_zero_mpb/(1-phi_zero_mpb)*alpha_func_mpb(i_full_points))
                        else !we will use 1/a^3*ln_cosh_function which is for a->0:
                            ln_cosh_function_calc(i_full_points) = &
                                2d0*c_mpb*alpha_func_mpb(i_full_points)
                        end if
                    end if
                end if
            end if !
            if (use_separate_alpha) then
                !if we have both a cationic and anionic alpha function
                if (partition_tab(i_full_points).gt.0d0) then
                    current_potential = pot_regul(i_full_points) + delta_mpbe_calc(i_full_points)
                    exp_cat = 0d0
                    exp_an = 0d0
                    if (alpha_func_anion_mpb(i_full_points).gt.0d0) then
                        exp_an = exp(z_mpb*current_potential/kBT_mpb)
                    end if
                    if (alpha_func_mpb(i_full_points).gt.0d0) then
                        exp_cat = exp(-z_mpb*current_potential/kBT_mpb)
                    end if
                    if (solve_pbe_only) then
                        ln_cosh_function_calc(i_full_points) = &
                            c_mpb*(alpha_func_mpb(i_full_points)*exp_cat+alpha_func_anion_mpb(i_full_points)*exp_an)
                    else
                        Afunc = 1d0
                        ln_cosh_function_calc(i_full_points) = 0d0
                        if (alpha_func_mpb(i_full_points).gt.0d0) then
                            part_cat = 1d0-phi_zero_mpb+phi_zero_mpb*alpha_func_mpb(i_full_points)*cosh(z_mpb*current_potential/kBT_mpb)
                            Afunc = Afunc -phi_zero_mpb/2d0*alpha_func_mpb(i_full_points)*exp_cat/part_cat
                            ln_cosh_function_calc(i_full_points) = ln_cosh_function_calc(i_full_points)+&
                                c_mpb*a_mpb**3*alpha_func_mpb(i_full_points)*exp_cat/part_cat*log(part_cat/(1d0-phi_zero_mpb))
                        end if
                        if (alpha_func_anion_mpb(i_full_points).gt.0d0) then
                            part_an = 1d0-phi_zero_mpb+phi_zero_mpb*alpha_func_anion_mpb(i_full_points)*cosh(z_mpb*current_potential/kBT_mpb)
                            Afunc = Afunc -phi_zero_mpb/2d0*alpha_func_anion_mpb(i_full_points)*exp_an/part_an
                            ln_cosh_function_calc(i_full_points) = ln_cosh_function_calc(i_full_points)+&
                            c_mpb*a_mpb**3*alpha_func_anion_mpb(i_full_points)*exp_an/part_an*log(part_an/(1d0-phi_zero_mpb))
                        end if
                        ln_cosh_function_calc(i_full_points) =  ln_cosh_function_calc(i_full_points)&
                            -Afunc*log(Afunc)
                    end if !pbe_only
                end if !not zero anyways
            end if !use_separate_alpha
        end do


    end subroutine evaluate_ln_cosh_function

    !****s* FHI-aims/evaluate_h_and_f_function
    !  NAME
    !  evaluate_h_and_f_function
    !  SYNOPSIS
    subroutine evaluate_h_and_f_function(pot_regul,pot_regul_gradient,&
            delta_mpbe_calc,delta_mpbe_gradient_calc,&
            h_function_calc,f_function_calc,&
            ion_conc_calc,ion_conc_deriv_calc)
    !  PURPOSE
    ! evaluate h_function = h^2 from eq. (34), needed for Newton method
    ! evaluate f_function = -pi4 * nion^mpb, with nion^mpb defined by eq. (6)
    ! and respective limits at high v

        ! reference for dielectric decrement: Nakayama and Andelmann J. Chem. Phys, 044706 (2015) -- Nakayama2015

        !evaluates the "kappa" function of the LPBE solved in each Newton step, called here h_function

        ! 	  use lpb_solver_utilities, only: cutf_glob
        use mpi_utilities
        use synchronize_mpi_basic, only: sync_real_number

        implicit none

        integer :: i_full_points
        real*8, dimension(n_full_points),intent(in) :: pot_regul
        real*8, dimension(3,n_full_points), intent(in) :: pot_regul_gradient
        real*8, dimension(n_full_points),intent(in) :: delta_mpbe_calc
        real*8, dimension(3,n_full_points),intent(in) :: delta_mpbe_gradient_calc
        real*8, dimension(n_full_points),intent(out) :: h_function_calc
        real*8, dimension(n_full_points),intent(out) :: f_function_calc
        real*8, dimension(2,n_full_points),intent(out) :: ion_conc_calc
        real*8, dimension(2,n_full_points), intent(out) :: ion_conc_deriv_calc
        ! eq(6) in Nakayama2015:
        real*8 :: rho_cat, rho_an, drho_cat_dv, drho_an_dv
        real*8, dimension(3) :: drho_cat_dgradv, drho_an_dgradv

        real*8 :: denom, current_potential,denom_anion,exp_cat,exp_an,sinh_term, cosh_term
        real*8, dimension(3) :: current_potential_gradient
        real*8 :: current_potential_gradient_squared
        real*8 :: deriv_func

        real*8, dimension(n_full_points) :: canion, ccation

        character*600 :: info_str

        !  	 if (dielectric_decrement_newton) then
        ! 	    write(use_unit,*) 'correct values here and everything exchange kappa_mpb with alpha_mpb'
        ! 	    stop
        ! 	  end if

        h_function_calc = 0d0
        f_function_calc = 0d0
        canion = 0d0
        ccation = 0d0
        ion_conc_calc=0d0
        ion_conc_deriv_calc=0d0


        do i_full_points = 1, n_full_points, 1
            if (.not.use_separate_alpha) then
                if (alpha_func_mpb(i_full_points).gt.0d0.and.(partition_tab(i_full_points).gt.0d0 .or. atomic_MERM)) then
                    !only for all points outside of cavity where we actually have ions
                    !if large potentials insert here, we may get NaN here
                    !if that happens we can use cutf_glob for f_function_calc:
                    ! 	    if (abs(current_potential).lt.cutf_glob) then
                    ! 	      denom = 1d0-phi_zero_mpb+phi_zero_mpb*cosh(z_mpb*current_potential/kBT_mpb)
                    ! 	      f_function_calc(i_full_points) = kappa_func(i_full_points)**2*&
                        ! 		sinh(z_mpb*current_potential/kBT_mpb)/&
                        ! 		denom
                    ! 	    else if (current_potential.lt.cutf_glob) then
                    ! 	      f_function_calc(i_full_points) = -kappa_func(i_full_points)**2*1/phi_zero_mpb
                    ! 	    else if (current_potential.gt.cutf_glob) then
                    ! 	      f_function_calc(i_full_points) = kappa_func(i_full_points)**2*1/phi_zero_mpb
                    ! 	    end if
                    current_potential = pot_regul(i_full_points) + delta_mpbe_calc(i_full_points)
                    current_potential_gradient = pot_regul_gradient(:,i_full_points) + delta_mpbe_gradient_calc(:,i_full_points)
                    if (abs(z_mpb*current_potential/kBT_mpb).gt.cut_sinh_cosh .and..not. correct_dielectric_decrement) then
                        !argument of sinh and cosh functions is VERY large leading to enourmos numbers (e.g. sinh(arg)>1e38).
                        !we are then however in regions where both functions have already reached a constant number
                        !this cutoff of 90 is always save! the worst that can happen is phi_zero_mpb = 1d-30. under this phi_zero_mpb, the
                        !functions will reach the limit 1/phi (f_function) or 0 (h_function) at around 70.
                        !SR: there is a problem with the new corrected expressions. the new expressions converge differenty fast to the limiting
                        !    value dependent on the value of the alpha function. the smaller values of alpha appear, the slower it converges
                        if (current_potential.lt.0d0) then
                            if (ffunction_warning_nan==0d0) then
                                ffunction_warning_nan = ffunction_warning_nan + 1d0
                            end if
                            ion_conc_calc(1,i_full_points) = c_mpb * 2/phi_zero_mpb
                            ion_conc_calc(2,i_full_points) = 0d0
                            h_function_calc(i_full_points) = 0d0
                            f_function_calc(i_full_points) = -kappainf_mpb**2/phi_zero_mpb
                        else if (current_potential.gt.0d0) then
                            if (ffunction_warning_nan==0d0) then
                                ffunction_warning_nan = ffunction_warning_nan + 1d0
                            end if
                            ion_conc_calc(2,i_full_points) = c_mpb * 2/phi_zero_mpb
                            ion_conc_calc(1,i_full_points) = 0d0
                            h_function_calc(i_full_points) = 0d0
                            f_function_calc(i_full_points) = kappainf_mpb**2/phi_zero_mpb
                        end if
                    else
                        !normal case, argument is in region where we can safely evaluate the hyperbolic function
                        !for a= 0, we do not need the limit considerations

                        if (correct_dielectric_decrement .and. dielectric_decrement_method == 'Newton' ) then
                            current_potential_gradient_squared = sum(current_potential_gradient**2,DIM=1)
                            exp_cat = exp(-z_mpb*current_potential/kBT_mpb + pi4_inv/(2d0*kBT_mpb)*&
                                deps_dcion(1,i_full_points)*current_potential_gradient_squared)
                            exp_an = exp(z_mpb*current_potential/kBT_mpb + 1d0/(2.*kBT_mpb)*&
                                deps_dcion(2,i_full_points)*current_potential_gradient_squared)

                        else
                            exp_cat = exp(-z_mpb*current_potential/kBT_mpb)
                            exp_an = exp(z_mpb*current_potential/kBT_mpb)
                            cosh_term = cosh(z_mpb*current_potential/kBT_mpb)
                            sinh_term = sinh(z_mpb*current_potential/kBT_mpb)
                        end if

                        if (.not.solve_pbe_only) then

                            if (correct_dielectric_decrement .and. dielectric_decrement_method == 'Newton' ) then
                                !we devide deps_dcion by 2, since we want the change per species not per salt molecule
                                rho_cat = c_mpb * alpha_func_mpb(i_full_points) * exp_cat
                                rho_an = c_mpb * alpha_func_mpb(i_full_points) * exp_an

                                denom = 1d0-phi_zero_mpb+a_mpb**3*(rho_cat+rho_an)

                                ion_conc_calc(1,i_full_points) =  rho_cat/denom
                                ion_conc_calc(2,i_full_points) =  rho_an/denom
                                !
                                !calculate dc_ion/dv
                                !SR: Attention here! We can only analytically calculate this derivative if eps(c_ion) is a linear function
                                ! of c_ion! For decrement_kind = 2 or 3, the concentrations depend on themselves!
                                drho_cat_dv = -c_mpb* alpha_func_mpb(i_full_points)*z_mpb/kBT_mpb*exp_cat
                                drho_an_dv = c_mpb* alpha_func_mpb(i_full_points)*z_mpb/kBT_mpb*exp_an
                                !dc_ion/dv ("d" means partial derivative)
                                ion_conc_deriv_calc(1,i_full_points) = (drho_cat_dv*denom-rho_cat*a_mpb**3*(drho_cat_dv+drho_an_dv))/denom**2
                                ion_conc_deriv_calc(2,i_full_points) = (drho_an_dv*denom-rho_an*a_mpb**3*(drho_cat_dv+drho_an_dv))/denom**2
                                !in order to get the full functional derivative delta c_ion / delta v, we also need the derivative
                                !with respect to the gradient of v
                                drho_cat_dgradv = rho_cat * pi4_inv/kBT_mpb*deps_dcion(1,i_full_points)*current_potential_gradient
                                drho_an_dgradv = rho_an * pi4_inv/kBT_mpb*deps_dcion(2,i_full_points)*current_potential_gradient
                                !dcion_dgradv = dcion/drhocat * drhocat/dgradv
                                dcion_dgradv(1,:,i_full_points) = &
                                    (drho_cat_dgradv*denom-rho_cat*a_mpb**3*(drho_cat_dgradv+drho_an_dgradv))/denom**2
                                dcion_dgradv(2,:,i_full_points) = &
                                    (drho_an_dgradv*denom-rho_an*a_mpb**3*(drho_cat_dgradv+drho_an_dgradv))/denom**2
                                !another derivative, we will need dc_ion/dalpha
                                dcion_dalpha(1,i_full_points) = c_mpb/denom * (exp_cat &
                                    - rho_cat/denom * a_mpb**3 * (exp_cat+exp_an))
                                dcion_dalpha(2,i_full_points) = c_mpb/denom * (exp_an &
                                    - rho_an/denom * a_mpb**3 * (exp_cat+exp_an))
                                !
                                !CORRECT THIS!!!:
                                h_function_calc(i_full_points) = 0d0
                                f_function_calc(i_full_points) = -kappainf_mpb**2/2d0*&
                                    (ion_conc_calc(1,i_full_points)+ion_conc_calc(2,i_full_points))
                            else

                                denom = 1d0-phi_zero_mpb+phi_zero_mpb*alpha_func_mpb(i_full_points)*cosh_term

                                ion_conc_calc(1,i_full_points) =  c_mpb * alpha_func_mpb(i_full_points)* exp_cat/denom
                                ion_conc_calc(2,i_full_points) =  c_mpb * alpha_func_mpb(i_full_points)* exp_an/denom
                                if (correct_dielectric_decrement .and. dielectric_decrement_method == 'SCF' ) then
                                    !   !dc_ion/dv
                                    ion_conc_deriv_calc(1,i_full_points) = c_mpb*&
                                        z_mpb/kBT_mpb*(-exp_cat*(1d0-phi_zero_mpb)-phi_zero_mpb)/denom**2
                                    ion_conc_deriv_calc(2,i_full_points) = c_mpb*&
                                        z_mpb/kBT_mpb*(exp_an*(1d0-phi_zero_mpb)+phi_zero_mpb)/denom**2
                                end if
                                h_function_calc(i_full_points) = alpha_func_mpb(i_full_points)*kappainf_mpb**2*z_mpb/kBT_mpb*&
                                    (phi_zero_mpb*alpha_func_mpb(i_full_points)-(phi_zero_mpb-1d0)*cosh_term)/denom**2
                                f_function_calc(i_full_points) = alpha_func_mpb(i_full_points)*kappainf_mpb**2*sinh_term/denom
                            end if

                            !
                            ! end if

                        else !use expressions for a=0
                            if (correct_dielectric_decrement .and. dielectric_decrement_method == 'Newton' ) then
                                rho_cat = c_mpb * alpha_func_mpb(i_full_points) * exp_cat
                                rho_an = c_mpb * alpha_func_mpb(i_full_points) * exp_an
                                ion_conc_calc(1,i_full_points) =  rho_cat
                                ion_conc_calc(2,i_full_points) =  rho_an
                                !
                                !calculate dc_ion/dv
                                !SR: Attention here! We can only analytically calculate this derivative if eps(c_ion) is a linear function
                                ! of c_ion! For decrement_kind = 2 or 3, the concentrations depend on themselves!
                                drho_cat_dv = -c_mpb* alpha_func_mpb(i_full_points)*z_mpb/kBT_mpb*exp_cat
                                drho_an_dv = c_mpb* alpha_func_mpb(i_full_points)*z_mpb/kBT_mpb*exp_an
                                !dc_ion/dv
                                ion_conc_deriv_calc(1,i_full_points) = drho_cat_dv
                                ion_conc_deriv_calc(2,i_full_points) = drho_an_dv
                                !CORRECT THIS!!!:
                                h_function_calc(i_full_points) = 0d0
                                f_function_calc(i_full_points) = -kappainf_mpb**2/2d0*(ion_conc_calc(1,i_full_points)+ion_conc_calc(2,i_full_points))
                            else
                                h_function_calc(i_full_points) = alpha_func_mpb(i_full_points)*kappainf_mpb**2*z_mpb/kBT_mpb*cosh_term
                                f_function_calc(i_full_points) = alpha_func_mpb(i_full_points)*kappainf_mpb**2*sinh_term
                                ion_conc_calc(1,i_full_points) =  c_mpb * alpha_func_mpb(i_full_points)*exp_cat
                                ion_conc_calc(2,i_full_points) =  c_mpb * alpha_func_mpb(i_full_points)*exp_an
                            end if
                        end if !solve_pbe_only
                    end if !use limit expression
                end if !if alpha>0 and partition_tab>0
            else !if use_separate_alpha
                if (.not.(alpha_func_mpb(i_full_points).eq.0d0.and.alpha_func_anion_mpb(i_full_points).eq.0d0).and.&
                        (partition_tab(i_full_points).gt.0d0 .or. atomic_MERM)) then
                    current_potential = pot_regul(i_full_points) + delta_mpbe_calc(i_full_points)
                    exp_cat = exp(-z_mpb*current_potential/kBT_mpb)
                    exp_an = exp(z_mpb*current_potential/kBT_mpb)
                    cosh_term = cosh(z_mpb*current_potential/kBT_mpb)
                    sinh_term = sinh(z_mpb*current_potential/kBT_mpb)
                    if (.not.solve_pbe_only) then
                        denom = 1d0-phi_zero_mpb+phi_zero_mpb*&
                            alpha_func_mpb(i_full_points)*cosh_term
                        denom_anion = 1d0-phi_zero_mpb+phi_zero_mpb*&
                            alpha_func_anion_mpb(i_full_points)*cosh_term
                        ! 		      write(use_unit,*) 'vals',current_potential,exp_cat,sinh_term
                        h_function_calc(i_full_points) = -kappainf_mpb**2/2d0*z_mpb/kBT_mpb*&
                            (exp_cat*alpha_func_mpb(i_full_points)*&
                            (-denom-phi_zero_mpb*alpha_func_mpb(i_full_points)*sinh_term)/denom**2+&
                            exp_an*alpha_func_anion_mpb(i_full_points)*&
                            (-denom_anion+phi_zero_mpb*alpha_func_anion_mpb(i_full_points)*sinh_term)/denom_anion**2)
                        f_function_calc(i_full_points) = -kappainf_mpb**2/2d0*&
                            (alpha_func_mpb(i_full_points)*exp_cat/denom -&
                            alpha_func_anion_mpb(i_full_points)*exp_an/denom_anion)
                        canion(i_full_points) = c_mpb*&
                            alpha_func_anion_mpb(i_full_points)*exp_an/denom_anion
                        ccation(i_full_points) = c_mpb*&
                            alpha_func_mpb(i_full_points)*exp_cat/denom
                        ion_conc_calc(1,i_full_points) =  c_mpb * alpha_func_mpb(i_full_points)* &
                            exp_cat/denom
                        ion_conc_calc(2,i_full_points) =  c_mpb * alpha_func_anion_mpb(i_full_points)*&
                            exp_an/denom_anion
                    else
                        h_function_calc(i_full_points) = kappainf_mpb**2/2d0*z_mpb/kBT_mpb*&
                            (alpha_func_mpb(i_full_points)*exp_cat+&
                            alpha_func_anion_mpb(i_full_points)*exp_an)
                        f_function_calc(i_full_points) = -kappainf_mpb**2/2d0*&
                            (alpha_func_mpb(i_full_points)*exp_cat-&
                            alpha_func_anion_mpb(i_full_points)*exp_an)
                        canion(i_full_points) = c_mpb*&
                            alpha_func_anion_mpb(i_full_points)*exp_an
                        ccation(i_full_points) = c_mpb*&
                            alpha_func_mpb(i_full_points)*exp_cat
                        ion_conc_calc(1,i_full_points) =  c_mpb * alpha_func_mpb(i_full_points)* &
                            exp_cat
                        ion_conc_calc(2,i_full_points) =  c_mpb * alpha_func_anion_mpb(i_full_points)*&
                            exp_an
                    end if !solve_pbe_only
                end if!relevant point
            end if !use_separate_alpha
        end do
        call sync_real_number(ffunction_warning_nan)

        if (.not.ffunction_warning_nan_done.and.(ffunction_warning_nan>0d0)) then
            write(info_str,'(2X,3A)') '|WARNING: Calculating the f- and h-functions in the MPBE solver with ',&
                'a very high potential, most probably due to your dielectric function settings. This can lead to NaN in the functions, therefore using a ',&
                'limit approximation for some points.'
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            ffunction_warning_nan_done = .True.
        end if

    end subroutine evaluate_h_and_f_function

    !   subroutine evaluate_partition_func_mpb(dielec_func_calc,dielec_func_grad_calc,&
        !     rho_gradient_calc,rho_laplacian_calc,&
        !     dielec_func_deriv_calc,dielec_func_2nd_deriv_calc,&
        !     partition_func_mpb_calc,&
        !     partition_func_mpb_gradient_calc,partition_func_mpb_laplacian_calc)
    !
    ! 	  implicit none
    !
    ! 	  !IN
    ! 	  real*8, dimension(n_full_points) :: dielec_func_calc
    ! 	  real*8, dimension(3,n_full_points) :: dielec_func_grad_calc
    ! 	  real*8, dimension(3,n_full_points) :: rho_gradient_calc
    ! 	  real*8, dimension(n_full_points) :: rho_laplacian_calc
    ! 	  real*8, dimension(n_full_points) :: dielec_func_deriv_calc
    ! 	  real*8, dimension(n_full_points) :: dielec_func_2nd_deriv_calc
    ! 	  !OUT
    ! 	  real*8, dimension(n_full_points) :: partition_func_mpb_calc
    ! 	  real*8, dimension(3,n_full_points) :: partition_func_mpb_gradient_calc
    ! 	  real*8, dimension(n_full_points) :: partition_func_mpb_laplacian_calc
    !
    ! 	  integer :: i_full_points,i_coords
    !
    !   ! 	if ((d-1d0).gt.1e-20) then
    ! 	  do i_full_points = 1, n_full_points, 1
    ! 	    if (partition_tab(i_full_points).gt.0d0) then
    ! 	      partition_func_mpb_calc(i_full_points) = 1d0/(epsinf_mpb-1d0)*&
        ! 		(dielec_func_calc(i_full_points)-1d0)
    ! 	      do i_coords = 1, 3, 1
    ! 		partition_func_mpb_gradient_calc(i_coords,i_full_points) = 1d0/(epsinf_mpb-1d0)*&
        ! 		  dielec_func_grad_calc(i_coords,i_full_points)
    ! 	      end do
    ! 	      partition_func_mpb_laplacian_calc(i_full_points) = 1d0/(epsinf_mpb-1d0) *&
        ! ! 		sum(rho_gradient_calc(:,i_full_points)**2,DIM=1)*dielec_func_2nd_deriv_calc(i_full_points)+&
        !   	      (sum(rho_gradient_calc(:,i_full_points)**2,DIM=1)*dielec_func_2nd_deriv_calc(i_full_points)+&
        !   	      rho_laplacian_calc(i_full_points) * dielec_func_deriv_calc(i_full_points))
    ! 	    end if
    ! 	  end do
    !   ! 	else
    !   ! 	  partition_tab_d = 0d0
    !   ! 	end if
    !
    ! ! 	  stop
    !
    ! 	  return
    !   end subroutine evaluate_partition_func_mpb

    !****s* FHI-aims/evaluate_iterative_mpbe_function
    !  NAME
    !  evaluate_iterative_mpbe_function
    !  SYNOPSIS
    subroutine evaluate_iterative_mpbe_function(dielec_func_calc,dielec_func_grad_calc, &
            pot_regul,delta_mpbe_calc, delta_mpbe_grad_calc, h_function_calc,f_function_calc,&
            q_times_epsinf_mpb_calc,v_free_gradient_calc,free_rho_calc,delta_rho_calc)
    !  PURPOSE
    ! evaluates q_times_epsinf_mpb_calc which is q * epsinf, with q defined by eq. (35)

        !evaluates the right hand site of the LPBE solved in each Newton step

        implicit none

    !INPUTS
        integer :: i_full_points, i_coords
        real*8, dimension(n_full_points), intent(in) :: dielec_func_calc
        real*8, dimension(3,n_full_points), intent(in) :: dielec_func_grad_calc
        real*8, dimension(n_full_points),intent(in) :: pot_regul
        real*8, dimension(n_full_points),intent(in) :: delta_mpbe_calc
        real*8, dimension(3,n_full_points), intent(in) :: delta_mpbe_grad_calc
        real*8, dimension(n_full_points),intent(in) :: h_function_calc
        real*8, dimension(n_full_points),intent(in) :: f_function_calc
        real*8, dimension(3,n_full_points),intent(in) :: v_free_gradient_calc
        real*8, dimension(n_full_points),intent(in) :: free_rho_calc
        real*8, dimension(n_full_points),intent(in) :: delta_rho_calc

    !OUTPUTS
        real*8, dimension(n_full_points),intent(out) :: q_times_epsinf_mpb_calc


        q_times_epsinf_mpb_calc = 0d0
        do i_full_points = 1, n_full_points, 1
            if (partition_tab(i_full_points).gt.0d0 .or. atomic_MERM) then
                !right hand side of (nabla (eps nabla ) -h ) v_mpbe = -4pi q_times_epsinf_mpb
                !  	    q_times_epsinf_mpb = pi4_inv*epsinf_mpb/dielec_func_calc(i_full_points)*&
                    !  	      (-kappa_func_calc(i_full_points)**2*z_mpb/kBT_mpb*pot_regul +&
                    !  	      f_function_calc(i_full_points)-h_function_calc(i_full_points)*delta_mpbe_calc(i_full_points))

                !  	    !now right hand side of: (-Delta+kappa**2)v_mpbe = 4pi rho_iter_mpbe_calc
                !  	    do i_coords = 1, 3, 1
                !  	      rho_iter_mpbe_calc(i_full_points) = rho_iter_mpbe_calc(i_full_points) + q_times_epsinf_mpb +&
                !  		epsinf_mpb*pi4_inv*(1d0/dielec_func_calc(i_full_points) *&
                    !  		dielec_func_grad_calc(i_coords,i_full_points)*delta_mpbe_grad_calc(i_coords,i_full_points)+ &
                    !  		(kappa_debye**2-h_function_calc(i_full_points)/dielec_func_calc(i_full_points))*delta_mpbe_calc(i_full_points))
                !  	    end do

                !right hand side of (nabla (eps nabla ) -h ) v_mpbe = -4pi q_times_epsinf_mpb
                if (reg_method == 'lpbe') then
                    q_times_epsinf_mpb_calc(i_full_points) = -pi4_inv*epsinf_mpb/dielec_func_calc(i_full_points)*&
                        (-alpha_func_mpb(i_full_points)*kappainf_mpb**2*z_mpb/kBT_mpb*pot_regul(i_full_points) +&
                        f_function_calc(i_full_points)-h_function_calc(i_full_points)*delta_mpbe_calc(i_full_points))
                    !  	      else if (reg_method == 'vacandconst') then
                    !  	    !in the comments I write the components of eps_s*q:
                    !  		q_times_epsinf_mpb_calc(i_full_points) = -pi4_inv*epsinf_mpb/dielec_func_calc(i_full_points)*&
                        !  	    !eps_s*nfree[1/eps-1]
                    !  		  (-free_rho_calc(i_full_points)*(1d0-dielec_func_calc(i_full_points))-&
                        !  		  pi4*delta_rho_calc(i_full_points)*(1d0-dielec_func_calc(i_full_points)+&
                        !  		  dielec_func_calc(i_full_points)*partition_func_mpb(i_full_points)*&
                        !  		  (1d0-1d0/epsinf_mpb))+&
                        !  	    !eps_s/eps*n_ion
                    !  		  f_function_calc(i_full_points)&
                        !  	    !eps_s/eps*(-1/4pi*h^2*delta v)
                    !  		  -h_function_calc(i_full_points)*&
                        !  		  delta_mpbe_calc(i_full_points)&
                        !  	    !eps_s/eps*(-1/4pi*eps(vconst-vvac)*Delta w)
                    !  		  -partition_func_mpb_laplacian(i_full_points)*dielec_func_calc(i_full_points)*&
                        !  		  (pot_helmh_const(i_full_points)-pot_helmh_vac(i_full_points))&
                        !  	    !eps_s/eps*(-1/4pi*eps*w*kappa^2*deltavconst)
                    !  		  -dielec_func_calc(i_full_points)*partition_func_mpb(i_full_points)*&
                        !  		  kappa_debye**2*pot_helmh_const(i_full_points))
                    !  		do i_coords = 1, 3, 1
                    !  		  q_times_epsinf_mpb_calc(i_full_points) = q_times_epsinf_mpb_calc(i_full_points) +&
                        !  	    !eps_s/eps*(1/4pi*grad eps*grad vfree)
                    !  		    pi4_inv*epsinf_mpb/dielec_func_calc(i_full_points)*&
                        !  		    dielec_func_grad_calc(i_coords,i_full_points)*&
                        !  		    v_free_gradient_calc(i_coords,i_full_points) +&
                        !  	    !eps_s/eps*(1/4pi*grad eps*grad vvac)
                    !  		    pi4_inv*epsinf_mpb/dielec_func_calc(i_full_points)*&
                        !  		    dielec_func_grad_calc(i_coords,i_full_points)*&
                        !  		    pot_helmh_vac_gradient(i_coords,i_full_points)+&
                        !  	    !w eps_s/eps*(1/4pi*grad eps*(grad vconst-grad vvac))
                    ! 		    partition_func_mpb(i_full_points)*&
                        ! 		    pi4_inv*epsinf_mpb/dielec_func_calc(i_full_points)*&
                        ! 		    dielec_func_grad_calc(i_coords,i_full_points)*&
                        ! 		    (pot_helmh_const_gradient(i_coords,i_full_points) - pot_helmh_vac_gradient(i_coords,i_full_points))+&
                        ! 	    !eps_s/eps * (1/4pi *eps  grad w * (grad vconst - grad vvac)
                    ! 		    pi4_inv*epsinf_mpb/dielec_func_calc(i_full_points)*&
                        ! 		    2d0*dielec_func_calc(i_full_points)*&
                        ! 		    partition_func_mpb_gradient(i_coords,i_full_points)*&
                        ! 		    (pot_helmh_const_gradient(i_coords,i_full_points) - pot_helmh_vac_gradient(i_coords,i_full_points))+&
                    ! 	    !eps_s/eps*(1/4pi*grad eps*grad w *(vconst - vvac)
                    ! 		    pi4_inv*epsinf_mpb/dielec_func_calc(i_full_points)*&
                        ! 		    dielec_func_grad_calc(i_coords,i_full_points)*&
                        ! 		    partition_func_mpb_gradient(i_coords,i_full_points)*&
                        ! 		    (pot_helmh_const(i_full_points) - pot_helmh_vac(i_full_points))
                    ! 		end do
                else
                    !in the comments I write the components of eps_s*q:
                    q_times_epsinf_mpb_calc(i_full_points) = -pi4_inv*epsinf_mpb/dielec_func_calc(i_full_points)*&
                        !eps_s/eps*nfree*(1-eps) + eps_s/eps*delta n = eps_s*[nfree/eps - nfree + deltan/eps] = eps_s*[n/eps-nfree]
                        (-free_rho_calc(i_full_points)*(1d0-dielec_func_calc(i_full_points))-&
                        pi4*delta_rho_calc(i_full_points)+&
                        !eps_s/eps*n_ion
                        f_function_calc(i_full_points)&
                        !eps_s/eps*(-1/4pi*h^2*delta v)
                        -h_function_calc(i_full_points)*&
                        delta_mpbe_calc(i_full_points))
                    do i_coords = 1, 3, 1
                        q_times_epsinf_mpb_calc(i_full_points) = q_times_epsinf_mpb_calc(i_full_points) +&
                            !eps_s/eps*(1/4pi*grad eps*grad vfree)
                            pi4_inv*epsinf_mpb/dielec_func_calc(i_full_points)*&
                            dielec_func_grad_calc(i_coords,i_full_points)*&
                            v_free_gradient_calc(i_coords,i_full_points)
                    end do
                end if
            end if !partition_tab
        end do



        ! 	    rho_iter_mpb = -epsinf_mpb/dielec_func_mpb * (&
            ! 	      free_rho*(1d0-dielec_func_mpb)+&
            ! 	      delta_rho_lpb-pi4_inv*kappa_func_mpb**2*z_mpb/kBT_mpb*&
            ! 		  free_hartree_superpos) + &
            ! 	      epsinf_mpb * pi4_inv * (kappa_debye**2-  &
            ! 	      (kappa_func_mpb / sqrt(dielec_func_mpb) * sqrt(z_mpb/kBT_mpb))**2) *&
            ! 	      pot_helmh
    end subroutine evaluate_iterative_mpbe_function

    !****s* FHI-aims/evaluate_rmsd_newton
    !  NAME
    !  evaluate_rmsd_newton
    !  SYNOPSIS
    subroutine evaluate_rmsd_newton(pot_new_calc, pot_old_calc, rmsd_calc,  tau_calc, converged_calc)
    !  PURPOSE
    ! evaluate root mean square change of potential during Newton iteration
        use mpi_utilities
        use synchronize_mpi

        implicit none

        integer :: i_full_points,i_myid,info
        real*8 :: rmsd_calc!, norm
        real*8 :: rmsd_mpi(n_tasks)
        real*8 :: norm_mpi(n_tasks)

        real*8, dimension(n_full_points), intent(in) :: pot_new_calc
        real*8, dimension(n_full_points), intent(in) :: pot_old_calc

        real*8 :: tau_calc

        logical, intent(out) :: converged_calc



        rmsd_calc = 0.0d0
        ! 	norm = 0.0d0

        i_full_points = 0

        do i_myid = 0, n_tasks-1, 1
            if (myid == i_myid) then
                rmsd_calc = 0
                do i_full_points = 1, n_full_points, 1
                    rmsd_calc = rmsd_calc + (pot_new_calc(i_full_points)- pot_old_calc(i_full_points))**2 &
                        * partition_tab(i_full_points)
                    ! 			norm = norm + partition_tab(i_full_points)
                end do
            end if

            call mpi_barrier(mpi_comm_global,info)
        end do

        call sync_real_number(rmsd_calc) !sums all rmsd_calc on different processors and saves them to rmsd
        ! 	call sync_real_number(norm)

        ! 	rmsd = sqrt(rmsd/norm) !/n_full_points
        rmsd_calc = sqrt(rmsd_calc) !/n_full_points

        if (rmsd_calc.lt.tau_calc) then
            converged_calc = .True.
        else
            converged_calc = .False.
        end if


        return

    end subroutine evaluate_rmsd_newton

    !subroutine get_gradient(   )
    !
    ! This subroutine calculates the gradient of the \delta v-potential w.r.t. the nuclear coordinate R_at.
    ! The nuclear coordinates refer of course to non-atom-centered coordinates.
    !
    ! USES

    ! INPUT
    ! o the maximum multipole order -- l_max_hartree
    ! o the current normalized direction in atom-centered spherical coordinates -- curr_ndir
    ! o the spherical harmonics -- ylm_tab
    ! o the gradient of the spherical harmonics -- grad_ylm
    ! o the multipole moment \delta v_{at,lm} -- multmom
    ! o the derivative of the multipole moment w.r.t. r_at -- deriv_multmom
    !
    !    integer*8
    !    real*8 curr_ndir(3)
    !    real*8 ylm_tab( n_atoms )
    !    real*8 grad_ylm( 3, n_atoms )
    !    real*8 multmom( 3, 2*(l_max_hartree + 1)**2 )
    !    real*8 deriv_multmom( 2*(l_max_hartree +1)**2 )

    ! OUTPUT
    ! o the gradient of \delta v

    !    real*8 grad_delta_v( 3, n_atoms )

    !   call tab_atom_centered_coordinates( coord_current,dist_tab,i_r,dir_tab )
    !  curr_ndir(:) = dir_tab( :, i_atom )

    ! call tab_single_wave_ylm_p0( current_center,&
        !           trigonom_tab,basis_l_max,l_ylm_max,&
        !            ylm_tab )
    !
    !do

    ! end do

    !end subroutine get_gradient

    !!subroutine get_gradient(  )
    !
    ! real*8 :: coord_current(3)
    ! real*8 :: dist_tab( n_atoms )
    ! real*8 :: i_r( n_atoms )
    ! real*8 :: dir_tab( 3,n_atoms )
    ! real*8 :: curr_ndir(3)
    !    real*8 ::



end module mpb_solver_utilities
