!****s* FHI-aims/SPE_solver
!  NAME
!   SPE_solver
!  AUTHOR
!   Stefan Ringe & FHI-aims team
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

module SPE_solver
!  PURPOSE
!    Solves Screened Poisson Equation (SPE) of the form:
!	(laplacian - kappa_SPE**2) * psi = -4 pi R_in
!	R_in is given by input, psi by output
!    The equation is solved via Green's function integration using a multipole expansion of the Green's function
!    (see paper above)
!    In the case of the Newton-MERM scheme to solve the MPBE, we solve the SPE given by eq. (36) with
!    R_in = q - pi4_inv * L_1 * delta v_{n+1}, and psi = delta v_{n+1} and kappa_SPE = kappa_debye
!
!    There are in principle 2 methods to use this module in combination with the MERM scheme
!    1. simply to integrate the SPE via multipole expansion (MERM_in_SPE_solver = .false.)
!    2. by performing the iterations of the MERM scheme directly in this routine (MERM_in_SPE_solver = .true.). This has the additional
!       advantage of higher speed (especially for c=0) by avoiding the full update of delta v_{n+1} on the full grid in each iterative
!       step
!  REMARKS
!    The routines have been copied from precondition_kerker and modified to suit the MPBE solver.
!  SOURCE

    !   logical :: kerker_preconditioner_on
    !   logical :: preconditioner_first_warnings = .true.
    !   integer, private :: n_spline_grid_dim, prec_l_dim
    !
    !   real*8,  dimension(:,:,:)   , allocatable, private :: green_I, green_K  ! radial Green functions of modified Helmholtz Eqn, fct of r, l, and species
    !   real*8,  dimension(:,:,:)   , allocatable, private :: R_multipole       ! angular and radially resolved multipole components
    !   real*8,  dimension(:,:,:,:) , allocatable, private :: R_multipole_spl   ! the spline functions of above multipole components - FIXME: use same mem as hartree pot!!
    !   real*8,  dimension(:,:)     , allocatable, private :: R_prec_multipole  ! the multipole version of the output residual
    !   real*8,  dimension(:,:)     , allocatable, private :: aux_R_spline            ! temporarily contains the splines for one angular momentum shell on one atom: rad grid
    !   real*8,  dimension(:,:)     , allocatable, private :: aux_R_prec_spline       ! same for logarithmic grid
    !   real*8,  dimension(:,:)     , allocatable, private :: angular_integral_log    ! temporarily contains the values of R_spl(alpha,l,m,r)
    !   real*8,  dimension(:)       , allocatable, private :: integral_zero_r         ! running variables for Green function integration
    !   real*8,  dimension(:)       , allocatable, private :: integral_r_infty
    !   real*8,  dimension(:,:,:)   , allocatable, private :: current_R_multipole_spl ! to store the multipole splines of one given atom for calculation of final answer
    !   real*8,  dimension(:)       , allocatable, private :: max_lm_spl_rad_sq       ! stores the maximum relevant radius of the logarithmically splined R_prec(l,m)
    !   real*8,  dimension(:)       , allocatable, private :: aux_R_prec_result       ! stores the spline result for R_prec in the final calculation
    !   real*8,  dimension(:,:)     , allocatable, private :: inner_radial_spl        ! stores the splines inside the two innermost radial points in real space
    !   integer, dimension(:,:)     , allocatable, private :: index_lm                ! (m,l)->1D translation array
    !   integer, dimension(:)       , allocatable, private :: prec_max_l_species
    !   integer, dimension(:)       , allocatable, private :: prec_l_dim_species

    use lpb_solver_utilities, only: lpb_converged,&
        initialized_lpb_solver,debug_lpb_solver,&
        solve_lpb_with_constant_dielec,kappainf_mpb,&
        kBT_mpb,z_mpb,&
        l_dim_SPE_public, l_dim_species_SPE_public, &
        kappa_debye,c_mpb,&
        evaluate_g_gradient,dielec_func_mpb,&
        evaluate_rmsd_merm,tau_merm,eta_merm,&
        dielec_func_gradient_mpb, limit_l_max_SPE_ff,&
        solve_lpbe_only, alpha_func_mpb,phi_zero_mpb,&
        mpb_forces_on, delta_v_MPB_gradient_atom,&
        forces_mpb_force_off, epsinf_mpb_inv, epsinf_mpb,&
        limit_l_max_SPE, kappa_SPE, multipole_radius_free_SPE,&
        atomic_MERM,sqrt_z_mpb_div_kBT_mpb
    !forces_hf_mpb

    use constants, only: pi4_inv

    !use physics, only: rho

    !variables for MPB solvent calculations, depending on choice of precondition_l_max_mpb
    logical, private :: SPE_first_warnings = .true.
    integer :: n_spline_grid_dim
    integer, private :: limit_l_dim_SPE
    integer, private :: limit_l_dim_SPE_ff
    integer :: n_max_spline_mpb = 4
    
    logical, dimension(:), allocatable :: atom_converged

    integer, dimension(:,:), allocatable :: index_atoms
    integer, dimension(:), allocatable ::index_atoms_size
    real*8,  dimension(:,:)     , allocatable, private :: aux_R_spline_mpb            ! temporarily contains the splines for one angular momentum shell on one atom: rad grid
    real*8,  dimension(:,:)     , allocatable, private :: aux_R_prec_spline_mpb       ! same for logarithmic grid
    real*8,  dimension(:)       , allocatable, private :: dylm_dtheta_tab
    real*8,  dimension(:)       , allocatable, private :: dylm_dtheta_tab_ff
    real*8,  dimension(:)       , allocatable, private :: scaled_dylm_dphi_tab
    real*8,  dimension(:)       , allocatable, private :: scaled_dylm_dphi_tab_ff
    real*8,  dimension(:,:)     , allocatable, private :: R_prec_multipole_mpb  ! the multipole version of the output residual
    real*8,  dimension(:,:,:)   , allocatable, private :: R_multipole_mpb       ! angular and radially resolved multipole components
    real*8,  dimension(:,:,:)   , allocatable, private :: R_multipole_mpb_old       ! angular and radially resodelta_v_MPB_at_nuclved multipole components
    real*8,  dimension(:,:,:)   , allocatable, private :: R_in_multipole       ! angular and radially resolved multipole components
    real*8,  dimension(:,:,:,:) , allocatable :: R_multipole_spl_mpb   ! the spline functions of above multipole components - FIXME: use same mem as hartree pot!!
    real*8,  dimension(:,:)     , allocatable, private :: angular_integral_log_mpb    ! temporarily contains the values of R_spl(alpha,l,m,r)
    real*8,  dimension(:,:)     , allocatable, private :: inner_parabolic_radial_splines    ! temporarily contains the values of R_spl(alpha,l,m,r)
    real*8,  dimension(:,:,:)   , allocatable :: current_R_multipole_spl_mpb ! to store the multipole splines of one given atom for calculation of final answer
    real*8,  dimension(:,:,:)   , allocatable, private :: green_I_mpb, green_K_mpb  ! radial Green functions of modified Helmholtz Eqn, fct of r, l, and species
    integer, dimension(:)       , allocatable :: limit_l_max_SPE_species
    integer, dimension(:)       , allocatable, private :: limit_l_max_SPE_species_ff
    integer, dimension(:)       , allocatable, private :: limit_l_dim_SPE_species_ff
    real*8,  dimension(:)       , allocatable, private :: integral_zero_r_mpb         ! running variables for Green function integration
    real*8,  dimension(:)       , allocatable, private :: integral_r_infty_mpb
    real*8,  dimension(:)       , allocatable :: max_lm_spl_rad_sq_mpb       ! stores the maximum relevant radius of the logarithmically splined R_prec(l,m)
    real*8,  dimension(:)       , allocatable :: aux_R_prec_result_mpb       ! stores the spline result for R_prec in the final calculation
    real*8,  dimension(:)       , allocatable, private :: aux_R_prec_result_mpb_ff
    real*8,  dimension(:)       , allocatable, private :: aux_R_prec_result_at_zero       ! stores the spline result for R_prec in the final calculation
    real*8,  dimension(:)       , allocatable, private :: aux_R_prec_result_at_zero_ff       ! stores the spline result for R_prec in the final calculation
    real*8,  dimension(:)       , allocatable, private :: aux_R_prec_deriv_result       ! stores the spline result for the radial derivative of R_prec in the final calculation (MPB solvent calculation)
    real*8,  dimension(:)       , allocatable, private :: aux_R_prec_deriv_result_ff       ! stores the spline result for the radial derivative of R_prec in the final calculation (MPB solvent calculation)
    real*8,  dimension(:)       , allocatable, private :: aux_R_prec_deriv_result_at_zero
    real*8,  dimension(:)       , allocatable, private :: aux_R_prec_deriv_result_at_zero_ff

    !   real*8, dimension(:,:), allocatable, private :: deriv_2nd
    !   real*8, dimension(:,:), allocatable, private :: deriv_1st
    !   real*8, dimension(:,:), allocatable, private :: multipole_moments
    real*8,  dimension(:)       , allocatable, private :: aux_R_prec_2nd_deriv_result       ! stores the spline result for the radial derivative of R_prec in the final calculation (MPB solvent calculation)
    real*8,  dimension(:)       , allocatable, private :: aux_R_prec_2nd_deriv_result_ff       ! stores the spline result for the radial derivative of R_prec in the final calculation (MPB solvent calculation)
    real*8,  dimension(:)       , allocatable, private :: aux_R_prec_gradient_result       ! stores the spline result for the gradient of R_prec*Ylm in the final calculation (MPB solvent calculation)
    real*8,  dimension(:,:)       , allocatable, private :: R_prec_gradient_at_zero       ! stores the spline result for the gradient of R_prec*Ylm in the final calculation (MPB solvent calculation)
    real*8,  dimension(:,:)     , allocatable, private :: inner_radial_spl_mpb        ! stores the splines inside the two innermost radial points in real space
    integer, dimension(:)       , allocatable :: limit_l_dim_SPE_species
    integer, dimension(:)       , allocatable :: limit_l_dim_SPE_species_current
    integer, dimension(:,:)     , allocatable :: index_lm_mpb                ! (m,l)->1D translation array
    real*8,  dimension(:)       , allocatable, private :: R_multipole_mpb_at_zero !monopole potential at r = 0
    real*8,  dimension(:,:)     , allocatable, private :: R_multipole_mpb_deriv_at_zero
    real*8, dimension(:) , allocatable, private :: R_prec_at_zero
    real*8, dimension(:,:,:), allocatable, private :: rho_multipole_mpb
    real*8, dimension(:,:,:), allocatable, private :: rho_multipole_mpb_in
!     real*8, dimension(:,:,:), allocatable, private :: rho_multipole_mpb_old
    real*8, dimension(:), allocatable, private :: R_prec_old
    real*8, dimension(:,:), allocatable, private :: R_prec_old_atoms
    real*8, dimension(:,:), allocatable, private :: R_prec_gradient_old
    real*8, dimension(:,:,:), allocatable, private :: R_prec_gradient_old_atoms
    
    real*8, dimension(:,:), allocatable, private :: multipole_moment_ff_glob


    real*8, dimension(:,:,:),   allocatable :: delta_v_MPB_multipole_component_at_zero
    real*8, dimension(:,:,:),   allocatable :: delta_v_MPB_multipole_deriv_at_zero
    real*8, dimension(:,:,:),   allocatable :: delta_v_MPB_gradient_at_zero


    real*8 :: aux_R_prec_laplacian_result


    integer :: sum_saved_zero = 0

contains

    !******

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !copied and changed Modules for solving the modified Helmholtz equation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !------------------------------------------------------------------------------
    !****s* SPE_solver/prepare_SPE_solver
    !  NAME
    !    prepare_SPE_solver
    !  SYNOPSIS
    subroutine prepare_SPE_solver
    !  PURPOSE
    !    This takes care of initialize all global variables for the SPE solver.
    !    Prepares the SPE solver to solve the modified helmholtz equation / SPE  for MPB
    !  USES
        use constants, only: pi
        use runtime_choices
        use dimensions
        use grids
        use localorb_io
        use species_data
        use mpi_tasks, only: aims_stop
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  SOURCE
        implicit none
        integer :: i_species, i_grid, i_m, i_l, i_index
        character*200 :: info_str
        real*8 :: kappa_times_rgrid
        integer :: i_full_points,  i_my_batch, current_atom


        !   if (avoid_grad_eps) then
        !     kappa_SPE = kappa_debye/SQRT(epsinf_mpb)*SQRT(eps_ref)
        !   else
        kappa_SPE = kappa_debye
        !   end if

        multipole_radius_free_SPE = multipole_radius_free !+ 30d0


        if (debug_lpb_solver) then
            write(info_str,*) 'precond', kappa_SPE
            call localorb_info(info_str)
        end if
        !   kappa_SPE = 1e-10

        ! ---------- Allocation and initialization ----------------------------------
        call deallocate_SPE_solver    ! Just in case the spline distribution & other things have failed ...

        ! go through the various species to check on the size of the angular momentum expansion
        if (.not.allocated(limit_l_max_SPE_species))            allocate(limit_l_max_SPE_species(n_species))
        if (.not.allocated(limit_l_dim_SPE_species))            allocate(limit_l_dim_SPE_species(n_species))
        if (.not.allocated(limit_l_dim_SPE_species_current))    allocate(limit_l_dim_SPE_species_current(n_species))
        if (.not.allocated(limit_l_max_SPE_species_ff))         allocate(limit_l_max_SPE_species_ff(n_species))
        if (.not.allocated(limit_l_dim_SPE_species_ff))         allocate(limit_l_dim_SPE_species_ff(n_species))


        if (n_centers_hartree_potential.ne.n_atoms) then
            call aims_stop('* WARNING!!!! n_centers_hartree_potential .ne. n_atoms...MPB integrator not working for periodic boundary conditions yet ')
        end if
        limit_l_max_SPE_species(:) = limit_l_max_SPE    ! initialize species-dependent l-values

        write(info_str,'(2X,A)') '| Actual lmax for multipole expansion to solve the SPE needed for the MPBE solver'
        call localorb_info(info_str)

        i_l = 0
        do i_species = 1, n_species
            if (l_hartree(i_species).lt.limit_l_max_SPE) then
                write(info_str,'(2X,2A)') '* WARNING!!!! l_hartree < limit_l_max_SPE for species ', species_name(i_species)
                call localorb_info(info_str)
                !         write(info_str,'(2X,A)')  '* The preconditioner would crash under this condition.'
                !         call localorb_info(info_str)
                !         write(info_str,'(2X,A)')  '* Setting limit_l_max_SPE = l_hartree for this species to prevent difficulties.'
                !         call localorb_info(info_str)
                limit_l_max_SPE_species(i_species) = l_hartree(i_species)
            end if
            i_l = max(i_l,limit_l_max_SPE_species(i_species))
            write(info_str,'(3X,3A,I4)') '| species ',species_name(i_species),': ', limit_l_max_SPE_species(i_species)
            call localorb_info(info_str)
        end do

        if (i_l.lt.limit_l_max_SPE) then
            write(info_str,'(2X,A)') '* The above warnings cause AIMS to reduce limit_l_max_SPE to an overall'
            call localorb_info(info_str)
            write(info_str,'(2X,A,I4,A)') '* value of ',i_l,' for overall consisitency.'
            call localorb_info(info_str)
            limit_l_max_SPE = i_l
        end if

        limit_l_dim_SPE_species(:) = (limit_l_max_SPE_species(:)+1)**2

        limit_l_dim_SPE_species_current = limit_l_dim_SPE_species
        l_dim_species_SPE_public = limit_l_dim_SPE_species
        ! R_multipole is only used in this routine, but it does contain a runtime choice of the max number of angular momentum shells
        ! same goes for R_multipole_spl
        limit_l_dim_SPE = (limit_l_max_SPE+1)**2                        ! total number of (l,m) pairs - this limit is easier to handle

        l_dim_SPE_public = limit_l_dim_SPE

        if (debug_lpb_solver) then
            write(use_unit,*) 'limit_l_dim_SPE, limit_l_max_SPE',limit_l_dim_SPE, limit_l_max_SPE
            write(use_unit,*) 'limit_l_dim_SPE_species',limit_l_dim_SPE_species
        end if

        if (.not.allocated(aux_R_spline_mpb))       allocate(aux_R_spline_mpb(n_max_spline_mpb, n_max_radial+1))
        if (.not.allocated(atom_converged)) allocate(atom_converged(n_spline_atoms))
        if (.not.allocated(index_atoms)) allocate(index_atoms(n_spline_atoms,n_full_points))
        if (.not.allocated(index_atoms_size)) allocate(index_atoms_size(n_spline_atoms))
        if (.not.allocated(aux_R_prec_spline_mpb))  allocate(aux_R_prec_spline_mpb(n_max_spline_mpb, n_max_grid))
        if (.not.allocated(R_multipole_mpb))        allocate(R_multipole_mpb(n_max_radial+1, limit_l_dim_SPE, n_atoms))
        if (.not.allocated(R_multipole_mpb_old))        allocate(R_multipole_mpb_old(n_max_radial+1, limit_l_dim_SPE, n_atoms))
        if (.not.allocated(R_in_multipole))        allocate(R_in_multipole(n_max_radial+1, limit_l_dim_SPE, n_atoms))
        if (.not.allocated(rho_multipole_mpb))        allocate(rho_multipole_mpb(n_max_radial+2, limit_l_dim_SPE, n_atoms))
        if (.not.allocated(rho_multipole_mpb_in))        allocate(rho_multipole_mpb_in(n_max_radial+2, limit_l_dim_SPE, n_atoms))
!         if (.not.allocated(rho_multipole_mpb_old))        allocate(rho_multipole_mpb_old(n_max_radial+2, limit_l_dim_SPE, n_atoms))

        ! this array is used twice in a row with different meanings each time, but one meaning after the other
        ! hence, over-allocation is not exactly necessary
        n_spline_grid_dim = max(n_max_grid,n_max_radial+1)

        if (.not.allocated(R_multipole_spl_mpb)) allocate(R_multipole_spl_mpb(limit_l_dim_SPE, n_max_spline_mpb, n_spline_grid_dim, n_spline_atoms))
        if (.not.allocated(R_multipole_mpb_at_zero)) allocate(R_multipole_mpb_at_zero(n_spline_atoms))
        if (.not.allocated(green_I_mpb))                 allocate(green_I_mpb(n_max_grid, limit_l_max_SPE+1, n_species))
        if (.not.allocated(green_K_mpb))                 allocate(green_K_mpb(n_max_grid, limit_l_max_SPE+1, n_species))
        if (.not.allocated(angular_integral_log_mpb))    allocate(angular_integral_log_mpb(limit_l_dim_SPE, n_max_grid ))
        if (.not.allocated(inner_parabolic_radial_splines))    allocate(inner_parabolic_radial_splines(limit_l_dim_SPE, n_max_radial ))
        if (.not.allocated(integral_zero_r_mpb))         allocate(integral_zero_r_mpb(limit_l_dim_SPE))
        if (.not.allocated(integral_r_infty_mpb))        allocate(integral_r_infty_mpb(limit_l_dim_SPE))
        if (.not.allocated(R_prec_multipole_mpb))        allocate(R_prec_multipole_mpb(n_max_grid, limit_l_dim_SPE))
        if (.not.allocated(current_R_multipole_spl_mpb)) allocate(current_R_multipole_spl_mpb(limit_l_dim_SPE, n_max_spline_mpb, n_spline_grid_dim))

        if (.not.allocated(max_lm_spl_rad_sq_mpb))       allocate(max_lm_spl_rad_sq_mpb(limit_l_dim_SPE))
        if (.not.allocated(aux_R_prec_result_mpb))       allocate(aux_R_prec_result_mpb(limit_l_dim_SPE))
        if (.not.allocated(aux_R_prec_result_at_zero))       allocate(aux_R_prec_result_at_zero(limit_l_dim_SPE))
        if (.not.allocated(aux_R_prec_deriv_result))       allocate(aux_R_prec_deriv_result(limit_l_dim_SPE))
        !   if (.not.allocated(deriv_2nd))       allocate(deriv_2nd(limit_l_dim_SPE,n_full_points))
        !   if (.not.allocated(deriv_1st))       allocate(deriv_1st(limit_l_dim_SPE,n_full_points))
        !   if (.not.allocated(multipole_moments))       allocate(multipole_moments(limit_l_dim_SPE,n_full_points))
        if (.not.allocated(aux_R_prec_2nd_deriv_result))       allocate(aux_R_prec_2nd_deriv_result(limit_l_dim_SPE))

        if (.not.allocated(aux_R_prec_gradient_result))       allocate(aux_R_prec_gradient_result(3))
        if (.not.allocated(index_lm_mpb))                allocate(index_lm_mpb(-limit_l_max_SPE:limit_l_max_SPE,0:limit_l_max_SPE))
        if (.not.allocated(inner_radial_spl_mpb))        allocate(inner_radial_spl_mpb(3,limit_l_dim_SPE))
        if (.not.allocated(dylm_dtheta_tab))       allocate(dylm_dtheta_tab(limit_l_dim_SPE))
        if (.not.allocated(scaled_dylm_dphi_tab))      allocate(scaled_dylm_dphi_tab(limit_l_dim_SPE))

        if (.not.allocated(R_prec_old)) allocate(R_prec_old(n_full_points))
        if (.not.allocated(R_prec_gradient_old)) allocate(R_prec_gradient_old(3,n_full_points))
        if (.not.allocated(R_prec_old_atoms)) allocate(R_prec_old_atoms(n_spline_atoms,n_full_points))
        if (.not.allocated(R_prec_gradient_old_atoms)) allocate(R_prec_gradient_old_atoms(3,n_spline_atoms,n_full_points))
        
        if (.not.allocated(R_prec_at_zero)) allocate(R_prec_at_zero(n_spline_atoms))


        !the far field lmax must always be lower than the normal lmax
        if (limit_l_max_SPE_ff.gt.limit_l_max_SPE) then
            write(info_str,'(3X,A)') "WARNING in SPE_solver: Far field lmax > near field lmax. This does not work in the&
                &current implementation and also makes no sense therefore settings lmax_ff = lmax"
            call localorb_info(info_str)
        end if
        limit_l_max_SPE_ff = min(limit_l_max_SPE,limit_l_max_SPE_ff)
        limit_l_dim_SPE_ff = (limit_l_max_SPE_ff+1)**2

        write(info_str,'(2X,A)') '| Actual lmax used in the far field'
        call localorb_info(info_str)
        do i_species = 1, n_species, 1
            limit_l_max_SPE_species_ff(i_species) = min(limit_l_max_SPE_ff,&
                limit_l_max_SPE_species(i_species))
            limit_l_dim_SPE_species_ff(i_species) = (limit_l_max_SPE_species_ff(i_species)+1)**2
            write(info_str,'(3X,3A,I4)') '| species ',species_name(i_species),': ', limit_l_max_SPE_species_ff(i_species)
        end do

        if (.not.allocated(multipole_moment_ff_glob)) allocate(multipole_moment_ff_glob(n_spline_atoms,limit_l_dim_SPE_ff))
        if (.not.allocated(aux_R_prec_result_mpb_ff)) allocate(aux_R_prec_result_mpb_ff(limit_l_dim_SPE_ff))
        if (.not.allocated(aux_R_prec_2nd_deriv_result_ff))       allocate(aux_R_prec_2nd_deriv_result_ff(limit_l_dim_SPE_ff))
        if (.not.allocated(scaled_dylm_dphi_tab_ff))      allocate(scaled_dylm_dphi_tab_ff(limit_l_dim_SPE_ff))
        if (.not.allocated(dylm_dtheta_tab_ff))       allocate(dylm_dtheta_tab_ff(limit_l_dim_SPE_ff))
        if (.not.allocated(aux_R_prec_deriv_result_ff))       allocate(aux_R_prec_deriv_result_ff(limit_l_dim_SPE_ff))
        if (.not.allocated(aux_R_prec_result_at_zero_ff))       allocate(aux_R_prec_result_at_zero_ff(limit_l_dim_SPE_ff))

        !FORCES
        if (use_forces.and..not.forces_mpb_force_off) then
            if (.not.allocated(R_multipole_mpb_deriv_at_zero)) allocate(R_multipole_mpb_deriv_at_zero(3,n_spline_atoms))
            if (.not.allocated(R_prec_gradient_at_zero))       allocate(R_prec_gradient_at_zero(3,n_spline_atoms))
            if (.not.allocated(aux_R_prec_deriv_result_at_zero))       allocate(aux_R_prec_deriv_result_at_zero(limit_l_dim_SPE))
            if (.not.allocated(aux_R_prec_deriv_result_at_zero_ff))       allocate(aux_R_prec_deriv_result_at_zero_ff(limit_l_dim_SPE_ff))

            if (.not.allocated(delta_v_MPB_multipole_component_at_zero)) then
                allocate (delta_v_MPB_multipole_component_at_zero(l_dim_SPE_public,n_spline_atoms,n_spline_atoms))
                delta_v_MPB_multipole_component_at_zero(:,:,:) = 0d0
            end if
            if (.not.allocated(delta_v_MPB_multipole_deriv_at_zero)) then
                allocate (delta_v_MPB_multipole_deriv_at_zero(l_dim_SPE_public,n_spline_atoms,n_spline_atoms))
                delta_v_MPB_multipole_deriv_at_zero(:,:,:) = 0d0
            end if
            if (.not.allocated(delta_v_MPB_gradient_at_zero)) then
                allocate (delta_v_MPB_gradient_at_zero(3,n_spline_atoms,n_spline_atoms))
                delta_v_MPB_gradient_at_zero(:,:,:) = 0d0
            end if
        end if

        ! map each (m,l) to an index and keep in storage for later use when splining
        i_index = 0
        index_lm_mpb = 0
        do i_l = 0, limit_l_max_SPE, 1
            do i_m = -i_l, i_l
                i_index = i_index + 1
                index_lm_mpb(i_m, i_l) = i_index
            enddo
        enddo

        ! initialize the radial Green functions for this particular value of q0; might possibly change throughout calculation ??????
        do i_species = 1, n_species
            do i_grid = 1, n_grid(i_species)

                ! calculate all l-values at once - for the radial argument (q_0 r)

                ! FIXME: investigate whether or not is is possible to switch arguments 1 and 2 in green_I and green_K for speed reasons?
                kappa_times_rgrid = kappa_SPE*r_grid(i_grid,i_species)
                call bessel_I_mpb(kappa_times_rgrid,limit_l_max_SPE, green_I_mpb(i_grid, :, i_species))
                call bessel_K_mpb(kappa_times_rgrid,limit_l_max_SPE, green_K_mpb(i_grid, :, i_species))
                ! now calculate modified spherical bessel function of 1st and 3rd (=modified spherical hankel function) kind by
                !i_n = sqrt(pi/(2precond r)*I_n+1/2 and k_n = sqrt(pi/(2precond r)*K_n+1/2
                !the monopole components are e.g. given by:
                !	k_0(precond r) = pi/(2precond r)*exp(-precond r)
                !	i_0(precond r) = sinh(precond r) / (precond r)
                green_I_mpb(i_grid,:,i_species) = green_I_mpb(i_grid,:,i_species)*sqrt(pi/(2.*kappa_times_rgrid))
                green_K_mpb(i_grid,:,i_species) = green_K_mpb(i_grid,:,i_species)*sqrt(pi/(2.*kappa_times_rgrid))
                !         green_I_mpb(i_grid,:,i_species) = green_I_mpb(i_grid,:,i_species)/sqrt(kappa_SPE*r_grid(i_grid,i_species))
                !         green_K_mpb(i_grid,:,i_species) = green_K_mpb(i_grid,:,i_species)/sqrt(kappa_SPE*r_grid(i_grid,i_species))
                ! 	if (any(green_I_mpb(i_grid,:,i_species).gt.1e-4)) then
                ! 	  write(use_unit,*) i_species,i_grid,'bessel1', green_I_mpb(i_grid,:,i_species)
                ! 	  write(use_unit,*) i_species,i_grid,'besselK', green_K_mpb(i_grid,:,i_species)
                ! 	end if
            end do         ! radial loop, initialization of bessel green function
        end do            ! species loop, initialization of bessel fct Green function

    
        !last but not least, devide the grid into atoms
        !SR: there should be an Aims array doing this?
        if (atomic_MERM) then
            i_full_points = 0
            index_atoms_size = 0
            index_atoms = 0
            do i_my_batch = 1, n_my_batches, 1                          ! loop through all total batches
                do i_index = 1, batches(i_my_batch)%size, 1          ! loop through all points in batch
                    i_full_points = i_full_points + 1
                    current_atom    = batches(i_my_batch) % points(i_index) % index_atom  
                    index_atoms_size(current_atom) = index_atoms_size(current_atom) + 1
                    index_atoms(current_atom,index_atoms_size(current_atom)) = i_full_points
                end do
            end do
        end if
    end subroutine prepare_SPE_solver
    !******

    !------------------------------------------------------------------------------
    !****s* SPE_solver/deallocate_SPE_solver
    !  NAME
    !    deallocate_SPE_solver
    !  SYNOPSIS
    subroutine deallocate_SPE_solver
    !  PURPOSE
    !    Clean up all of the variables allocated above
    !  USES
        use runtime_choices
        use dimensions
    !  SOURCE
        implicit none

        if (allocated(multipole_moment_ff_glob)) deallocate(multipole_moment_ff_glob)
        if (allocated(limit_l_max_SPE_species_ff)) deallocate(limit_l_max_SPE_species_ff)
        if (allocated(limit_l_dim_SPE_species_ff)) deallocate(limit_l_dim_SPE_species_ff)
        if (allocated(inner_radial_spl_mpb))        deallocate(inner_radial_spl_mpb)
        if (allocated(limit_l_max_SPE_species))      deallocate(limit_l_max_SPE_species)
        if (allocated(limit_l_dim_SPE_species))      deallocate(limit_l_dim_SPE_species)
        if (allocated(limit_l_dim_SPE_species_current))      deallocate(limit_l_dim_SPE_species_current)
        if (allocated(aux_R_prec_gradient_result))       deallocate(aux_R_prec_gradient_result)
        if (allocated(dylm_dtheta_tab))       deallocate(dylm_dtheta_tab)
        if (allocated(scaled_dylm_dphi_tab))       deallocate(scaled_dylm_dphi_tab)
        if (allocated(R_prec_multipole_mpb))       deallocate(R_prec_multipole_mpb)
        if (allocated(aux_R_prec_result_mpb))       deallocate(aux_R_prec_result_mpb)
        if (allocated(aux_R_prec_result_mpb_ff))       deallocate(aux_R_prec_result_mpb_ff)
        if (allocated(aux_R_prec_deriv_result)) deallocate(aux_R_prec_deriv_result)
        if (allocated(aux_R_prec_deriv_result_ff)) deallocate(aux_R_prec_deriv_result_ff)
        if (allocated(aux_R_prec_2nd_deriv_result)) deallocate(aux_R_prec_2nd_deriv_result)
        if (allocated(aux_R_prec_2nd_deriv_result_ff)) deallocate(aux_R_prec_2nd_deriv_result_ff)
        if (allocated(aux_R_prec_result_at_zero))       deallocate(aux_R_prec_result_at_zero)
        if (allocated(aux_R_prec_result_at_zero_ff))       deallocate(aux_R_prec_result_at_zero_ff)

        if (allocated(R_multipole_mpb))       deallocate(R_multipole_mpb)
        if (allocated(R_multipole_mpb_old))       deallocate(R_multipole_mpb_old)
        if (allocated(R_in_multipole))       deallocate(R_in_multipole)
        if (allocated(rho_multipole_mpb))       deallocate(rho_multipole_mpb)
        if (allocated(rho_multipole_mpb_in))       deallocate(rho_multipole_mpb_in)
!         if (allocated(rho_multipole_mpb_old))       deallocate(rho_multipole_mpb_old)

        if (allocated(R_multipole_spl_mpb))       deallocate(R_multipole_spl_mpb)
        if (allocated(angular_integral_log_mpb))       deallocate(angular_integral_log_mpb)
        if (allocated(inner_parabolic_radial_splines)) deallocate(inner_parabolic_radial_splines)
        if (allocated(current_R_multipole_spl_mpb))       deallocate(current_R_multipole_spl_mpb)
        if (allocated(green_I_mpb))       deallocate(green_I_mpb)
        if (allocated(green_K_mpb))       deallocate(green_K_mpb)
        if (allocated(integral_zero_r_mpb))       deallocate(integral_zero_r_mpb)
        if (allocated(integral_r_infty_mpb))       deallocate(integral_r_infty_mpb)
        if (allocated(max_lm_spl_rad_sq_mpb))       deallocate(max_lm_spl_rad_sq_mpb)
        if (allocated(inner_radial_spl_mpb))       deallocate(inner_radial_spl_mpb)
        if (allocated(aux_R_spline_mpb))            deallocate(aux_R_spline_mpb)
        if (allocated(atom_converged)) deallocate(atom_converged)
        if (allocated(index_atoms)) deallocate(index_atoms)
        if (allocated(index_atoms_size)) deallocate(index_atoms_size)
        if (allocated(aux_R_prec_spline_mpb))       deallocate(aux_R_prec_spline_mpb)
        if (allocated(index_lm_mpb))                deallocate(index_lm_mpb)
        if (allocated(R_multipole_mpb_at_zero)) deallocate(R_multipole_mpb_at_zero)
        if (allocated(R_prec_at_zero)) deallocate(R_prec_at_zero)
        if (allocated(R_prec_old)) deallocate(R_prec_old)
        if (allocated(R_prec_gradient_old)) deallocate(R_prec_gradient_old)
        if (allocated(R_prec_old_atoms)) deallocate(R_prec_old_atoms)
        if (allocated(R_prec_gradient_old_atoms)) deallocate(R_prec_gradient_old_atoms)
        
        !forces
        if (use_forces.and.forces_mpb_force_off) then
            if (allocated(R_prec_gradient_at_zero)) deallocate(R_prec_gradient_at_zero)
            if (allocated(R_multipole_mpb_deriv_at_zero)) deallocate(R_multipole_mpb_deriv_at_zero)
            if (allocated(aux_R_prec_deriv_result_at_zero)) deallocate(aux_R_prec_deriv_result_at_zero)
            if (allocated(aux_R_prec_deriv_result_at_zero_ff)) deallocate(aux_R_prec_deriv_result_at_zero_ff)
            if (allocated(delta_v_MPB_multipole_component_at_zero)) then
                deallocate(delta_v_MPB_multipole_component_at_zero)
            end if

            if (allocated(delta_v_MPB_multipole_deriv_at_zero)) then
                deallocate(delta_v_MPB_multipole_deriv_at_zero)
            end if
            if (allocated(delta_v_MPB_gradient_at_zero)) then
                deallocate(delta_v_MPB_gradient_at_zero)
            end if
        end if


        !******
    end subroutine deallocate_SPE_solver

    !------------------------------------------------------------------------------
    !****s* SPE_solver/reinitialize_green_mpb
    !  NAME
    !    reinitialize_green_mpb
    !  SYNOPSIS
    subroutine reinitialize_green_mpb(kappa_SPE_step)
    !  PURPOSE
    !    Reinitializes Bessel and Hankel functions for specific value of kappa_SPE_step, so that if different kappa's should be used
    !    in a single FHI-aims calculation not the whole SPE_solver has to be reinitialized
        use constants, only: pi
        use runtime_choices
        use dimensions
        use grids
        use localorb_io
        use species_data

        implicit none

        real*8 :: kappa_SPE_step
        integer :: i_species, i_grid
        real*8 :: kappa_times_rgrid

        green_I_mpb = 0.0d0
        green_K_mpb = 0.0d0

        ! initialize the radial Green functions for this particular value of q0; might possibly change throughout calculation ??????
        do i_species = 1, n_species
            do i_grid = 1, n_grid(i_species)

                ! calculate all l-values at once - for the radial argument (q_0 r)
                kappa_times_rgrid = kappa_SPE_step*r_grid(i_grid,i_species)
                call bessel_I_mpb(kappa_times_rgrid,limit_l_max_SPE, green_I_mpb(i_grid, :, i_species))
                call bessel_K_mpb(kappa_times_rgrid,limit_l_max_SPE, green_K_mpb(i_grid, :, i_species))
                ! now calculate modified spherical bessel function of 1st and 3rd kind by i_n = sqrt(pi/(2precond r)*I_n+1/2 and k_n = sqrt(pi/(2precond r)*K_n+1/2
                green_I_mpb(i_grid,:,i_species) = green_I_mpb(i_grid,:,i_species)*sqrt(pi/(2.*kappa_times_rgrid))
                green_K_mpb(i_grid,:,i_species) = green_K_mpb(i_grid,:,i_species)*sqrt(pi/(2.*kappa_times_rgrid))
            end do         ! radial loop, initialization of bessel green function
        end do            ! species loop, initialization of bessel fct Green function

    end subroutine reinitialize_green_mpb

    !------------------------------------------------------------------------------
    !****s* SPE_solver/partition_SPE
    !  NAME
    !    partition_SPE
    !  SYNOPSIS
    subroutine partition_SPE(partition_steps, kappa_SPE_step,dielec_func_current)
    !  PURPOSE
    ! partitions the dielectric function/ kappa space into partition_steps pieces, at which the Helmholtz equation is solved

        implicit none

        integer :: partition_steps
        real*8 :: dielec_interval, kappa_interval
        real*8 :: kappa_func_current
        real*8, dimension(partition_steps) :: kappa_SPE_step
        real*8, dimension(partition_steps) :: dielec_func_current
        integer :: i_part_step

        kappa_SPE_step = 0.0d0
        dielec_interval = (epsinf_mpb-1)/(DBLE(partition_steps)-1)
        kappa_interval = (kappainf_mpb)/(DBLE(partition_steps)-1)

        do i_part_step = 1, partition_steps, 1
            dielec_func_current(i_part_step) = 1 + dielec_interval * (i_part_step-1)
            kappa_func_current = kappa_interval * (i_part_step-1)
            kappa_SPE_step(i_part_step) = kappa_func_current*&
                sqrt(z_mpb/(kBT_mpb*dielec_func_current(i_part_step)))
        end do

        kappa_SPE_step(1) = 1d-10

    end subroutine partition_SPE

    !------------------------------------------------------------------------------
    !****s* SPE_solver/evaluate_density_mpe
    !  NAME
    !    evaluate_density_mpe
    !  SYNOPSIS
    subroutine evaluate_density_mpe(R_in_copy,dielec_func_current,set_rho_mp,&
            use_partition,full_iteration,take_part,R_out, R_out_gradient)
    !  PURPOSE
    ! evaluate the mp moments of the density R_multipole_mpb and saves it also to rho_multipole_mpb
    ! this is the integral over the angular grid points of each shell
    ! solves angular integral of eq. (42) & (43)
    !  USES
        use lpb_solver_utilities
        use mpb_solver_utilities, only: h_function
        use dimensions
        use grids
        use physics, only: hartree_partition_tab, hartree_weight_tab
        use synchronize_mpi
        use geometry, only: species

        implicit none

        logical :: full_iteration
        integer :: take_part
        real*8, intent(in) :: dielec_func_current
        real*8, dimension(n_full_points), intent(in) :: R_in_copy
        real*8, dimension(limit_l_dim_SPE) :: ylm_tab
        real*8, dimension(n_full_points), intent(in) :: R_out
        real*8, dimension(3,n_full_points), intent(in) :: R_out_gradient
        logical, intent(in) :: use_partition, set_rho_mp
        real*8 :: temp_rho, kappa_current
        integer :: i_full_points, i_my_batch, i_index, i_coords
        integer :: current_atom, current_radial, current_angular
        real*8 :: hartree_partition_tab_choice

        integer :: i_atom
        
        real*8 :: current_rout
        real*8, dimension(3) :: current_rout_grad

        i_full_points = 0
        if (.not.full_iteration) then
            rho_multipole_mpb = 0d0
            R_multipole_mpb = 0d0
        end if

!         if (atomic_MERM) then
!             n_batches_iter = n_grid_batches
!             !each CPU has to run over the whole grid
!         else
!             n_batches_iter = n_my_batches
!         end if
            
        do i_my_batch = 1, n_my_batches, 1                          ! loop through all total batches

!             if (atomic_MERM) then
!                 n_index_iter = batch_sizes(i_batch)
!                 !each CPU has to run over the whole grid
!             else 
!                 n_index_iter = batches(i_my_batch)%size
!             end if
        
            do i_index = 1, batches(i_my_batch)%size, 1          ! loop through all points in batch
            
                ! get the technical details of current grid point
                current_atom    = batches(i_my_batch) % points(i_index) % index_atom
                !       if (current_atom==2) then !only c atom center
                current_radial  = batches(i_my_batch) % points(i_index) % index_radial
                current_angular = batches(i_my_batch) % points(i_index) % index_angular
                i_full_points   = i_full_points + 1
                
                if (atomic_MERM .and. atom_converged(current_atom)) then
                    cycle
                end if

                if (atomic_MERM) then
                    hartree_partition_tab_choice = hartree_weight_tab(i_full_points) !=w_angular(current_angular, current_radial, species(current_atom)) * pi4
                else
                    hartree_partition_tab_choice = hartree_partition_tab(i_full_points)
                end if


                ! treat grid point in question to obtain multipole expansion
                if (hartree_partition_tab(i_full_points).gt.0.d0 .or. atomic_MERM) then ! - well, only if there is a non-zero summation weight of course

                    kappa_current = sqrt(alpha_func_mpb(i_full_points)/(1d0+phi_zero_mpb*(alpha_func_mpb(i_full_points)-1d0)))*kappainf_mpb / sqrt(dielec_func_mpb(i_full_points)) *sqrt_z_mpb_div_kBT_mpb

                    if (partition_tab_mpb(0,i_full_points,take_part,0,.false.,kappa_current)==1) then
                        current_rout = R_out(i_full_points)
                        current_rout_grad = R_out_gradient(:,i_full_points)
                        ! look up the sperical harmonics for this point ....
                        ylm_tab (1:limit_l_dim_SPE_species(species(current_atom))) = &
                            local_ylm_tab(1:limit_l_dim_SPE_species(species(current_atom)),current_angular, &      ! local_ylm_tab contains information for all lm values serially
                            lebedev_grid_index(current_radial,species(current_atom)))

                        if (take_part==0) then !normal case, just multipole expand input function
                            temp_rho = 0.0d0
                            if (use_partition) then
                                temp_rho = temp_rho + R_in_copy(i_full_points) / dielec_func_current!already scaled by position dependent eps
                            else
                                temp_rho = temp_rho + R_in_copy(i_full_points) * epsinf_mpb_inv !have to scale, negativ because charge density is required here
                            end if

                            ! ... to calculate the multipole expansion contribution by implied looping over all (l,m)
                            R_multipole_mpb(current_radial, 1:limit_l_dim_SPE_species(species(current_atom)), current_atom) = &
                                R_multipole_mpb(current_radial, 1:limit_l_dim_SPE_species(species(current_atom)), current_atom) + &          
                                ylm_tab(1:limit_l_dim_SPE_species(species(current_atom)))*&
                                hartree_partition_tab(i_full_points)*&
                                temp_rho

                        else if (take_part==1) then !LPBE case, multipole expand iterative density
                            R_multipole_mpb(current_radial, 1:limit_l_dim_SPE_species(species(current_atom)), current_atom) = &
                                R_multipole_mpb(current_radial, 1:limit_l_dim_SPE_species(species(current_atom)), current_atom) + &          
                                ylm_tab(1:limit_l_dim_SPE_species(species(current_atom)))*current_rout*&
                                pi4_inv*(kappa_debye**2-kappa_current**2)*&
                                hartree_partition_tab_choice
                            do i_coords = 1, 3, 1
                                R_multipole_mpb(current_radial, 1:limit_l_dim_SPE_species(species(current_atom)), current_atom) = &
                                    R_multipole_mpb(current_radial, 1:limit_l_dim_SPE_species(species(current_atom)), current_atom) + &          
                                    ylm_tab(1:limit_l_dim_SPE_species(species(current_atom)))* current_rout_grad(i_coords)*&
                                    pi4_inv*dielec_func_gradient_mpb(i_coords,i_full_points)/dielec_func_mpb(i_full_points)*&
                                    hartree_partition_tab_choice
                            end do
                        else if (take_part==2) then !MPBE case, multipole expand iterative density
                            R_multipole_mpb(current_radial, 1:limit_l_dim_SPE_species(species(current_atom)), current_atom) = &
                                R_multipole_mpb(current_radial, 1:limit_l_dim_SPE_species(species(current_atom)), current_atom) + &          
                                ylm_tab(1:limit_l_dim_SPE_species(species(current_atom)))*current_rout*&
                                pi4_inv*(kappa_debye**2-h_function(i_full_points)/dielec_func_mpb(i_full_points))*&
                                hartree_partition_tab_choice
                            do i_coords = 1, 3, 1
                                R_multipole_mpb(current_radial, 1:limit_l_dim_SPE_species(species(current_atom)), current_atom) = &
                                    R_multipole_mpb(current_radial, 1:limit_l_dim_SPE_species(species(current_atom)), current_atom) + &          
                                    ylm_tab(1:limit_l_dim_SPE_species(species(current_atom)))* current_rout_grad(i_coords)*&
                                    pi4_inv*dielec_func_gradient_mpb(i_coords,i_full_points)/dielec_func_mpb(i_full_points)*&
                                    hartree_partition_tab_choice
                            end do
                        end if !take part
                        !		!here we add the missing multipole moments of the delta density from the Hartree potential evaluation
                        !		do i_lm = 1, limit_l_dim_SPE_species(species(current_atom))
                        !		  !there can be problems here. lmax of hartree and kerker must be the same and also
                        !		  !hartree must be evaluated for EVERY atom!
                        !		  if (use_partition) then
                        !		    R_multipole_mpb(current_radial, i_lm, current_atom) = R_multipole_mpb(current_radial, i_lm, current_atom) &
                            !		      - 1.0d0/dielec_func_mpb(i_full_points)*&
                            !			rho_multipole_save_for_lpb(i_lm,current_radial, current_atom )
                        !		  else
                        !		    R_multipole_mpb(current_radial, i_lm, current_atom) = R_multipole_mpb(current_radial, i_lm, current_atom) &
                            !		      - 1.0d0/dielec_func_mpb(i_full_points)*&
                            !			rho_multipole_save_for_lpb(i_lm,current_radial, current_atom )
                        !		  end if
                        !		end do
                    end if !if partition_tab_mpb ==1
                end if       ! (hartree_partition_tab(i_full_points).gt.0.d0) ?
            end do          ! index = 1, batches(i_mybatch)%size
        end do                ! i_batch = 1, n_grid_batches, 1

        ! synchronize multipole expansion across all threads	
        call sync_kerker_multipole(R_multipole_mpb, n_max_radial+1, limit_l_dim_SPE)
        if (set_rho_mp) then
        !we save the angular integrals in order to later reconstruct the multipole error on the source term q + L_1 * delta_v
        !which we need to evaluate the multipole correction on the multipole error
             if (atomic_MERM) then
                !in this case we only save the atoms that were still calculate. we neglect the already converged atoms
                !since those are already stored correctly in rho_multipole_mpb
                do i_atom = 1, n_spline_atoms, 1
                    if (.not. atom_converged(i_atom)) then
                        rho_multipole_mpb(2:,:,i_atom) = R_multipole_mpb(:,:,i_atom)
                    end if
                end do
             else
                 rho_multipole_mpb(2:,:,:) = R_multipole_mpb
             end if
        end if
        
        ! and save the angular integrals, in case we simply evaluated the radial integral of "q", since we need to add this
        ! on top of the angular integral of L_1 * delta_v in later iterations
       
        if (full_iteration.and.take_part==0) then
            R_in_multipole = R_multipole_mpb
        end if
        if (set_rho_mp) then
           ! call sync_kerker_multipole(rho_multipole_mpb, n_max_radial+2, limit_l_dim_SPE)
            if (full_iteration.and.take_part==0) then
                rho_multipole_mpb_in = rho_multipole_mpb
            end if
        end if

        !     do i_lm = 1, limit_l_dim_SPE_species(species(1)), 1
        !       write(loopstr,'(I3.3)') i_lm
        !       filewrite = loopstr//"_mp_pot.dat"
        !       open(100,file=filewrite)
        ! !       write(100,*) multipole_moment_ff(1,i_lm)
        !       do i_radial = 1, n_radial(species(1)), 1
        ! 	write(100,*) r_radial(i_radial,species(1)),R_multipole_mpb(i_radial,i_lm,1)
        !       end do
        ! !       close(100)
        ! ! !       call output_delta_v_MPB_step(multipole_output(:,i_lm,i_atom),i_lm+100*(i_atom-1))
        !     end do
        !     stop


        !output monopole moments
        !    if (myid==0) open(111,file='mp_comp.dat')
        !    do i_atom = 1, n_atoms, 1
        !     do i_radial = 1, n_radial(species(i_atom)), 1
        !       if (myid==0) write(111,'(I,2(E,1X))') i_atom, r_radial(i_radial,species(i_atom)), R_multipole_mpb(i_radial,1,i_atom)
        !     end do
        !    end do
        !
        !    if (myid==0) open(111,file='mp_comp.dat')
        !    i_atom = 2
        !     do i_radial = 1, n_radial(species(i_atom)), 1
        !       if (myid==0) write(111,'(6(E,1X))') r_radial(i_radial,species(i_atom)), &
            ! 	(R_multipole_mpb(i_radial,i_lm,i_atom),i_lm=1,5,1)
        !     end do
        !
        !   stop

    end subroutine evaluate_density_mpe

    !------------------------------------------------------------------------------
    !****s* SPE_solver/evaluate_radial_integral
    !  NAME
    !    evaluate_radial_integral
    !  SYNOPSIS
    subroutine evaluate_radial_integral(set_rho_mp, use_partition, use_analytic_far_field,&
            kappa_SPE_step, multipole_moment_ff,&
            full_iteration, evaluate_only_mpe)
    !  PURPOSE
    ! Performs the radial integration of eq. (41)

        use dimensions, only: n_spline_atoms, n_atoms, n_max_radial, &
            use_distributed_spline_storage
        use runtime_choices, only: extra_adding_to_hartree_potential_distance
        use grids
        use spline, only: cubic_spline, spline_vector
        use synchronize_mpi
        use geometry, only: species
        use mpi_utilities, only: task_list, spline_atom_storage

        implicit none

        logical, intent(in) :: set_rho_mp, use_partition, use_analytic_far_field,&
            full_iteration, evaluate_only_mpe
        real*8, intent(in) :: kappa_SPE_step
        real :: drho
        integer :: i_atom, n_grid_limit, i_l, i_m, i_radial, n_inner_grid,&
            i_spline, i_grid
        real*8 :: delta, delta_2, delta_3, r1, r2, f1, f2, denominator,&
            radius_sq, radius, alpha, i_r_outer
        integer :: current_spl_atom
        real*8, dimension(n_spline_atoms,limit_l_dim_SPE_ff) :: multipole_moment_ff


        ! spline the above multipole expansion for use with the logarithmic grid
        !       code is the same as in update_hartree_potential_p1; except that two array arguments in the splines are switched as suggested in that routine
        !
        ! this part of the code runs over all atoms - we are no longer on the original grid because there now is a multipole
        ! expansion that has a basis which is independent of the grid points.
        ! This also means that each thread must know whether or not the atom under consideration is in its task list ...

        ! The same loop also calculates the atomic multipole expansion of the preconditioned residual akin to routine integrate_hartree_log_grid_p1

        ! loop over all atoms
        real*8 :: R_multipole_spl_mpb_temp(limit_l_dim_SPE, n_max_spline_mpb, n_spline_grid_dim, n_spline_atoms)
        real*8 :: R_multipole_mpb_at_zero_temp(n_spline_atoms)
        real*8, dimension(n_spline_atoms,limit_l_dim_SPE_ff) :: multipole_moment_ff_glob_temp, multipole_moment_ff_temp
        real*8 :: R_multipole_mpb_deriv_at_zero_temp(3,n_spline_atoms)

        if (atomic_MERM) then
            !need to save the data from previous iterations (already converged atoms), here
            R_multipole_spl_mpb_temp = R_multipole_spl_mpb
            multipole_moment_ff_glob_temp = multipole_moment_ff_glob
            multipole_moment_ff_temp = multipole_moment_ff
            R_multipole_mpb_at_zero_temp = R_multipole_mpb_at_zero
            R_multipole_mpb_deriv_at_zero_temp = R_multipole_mpb_deriv_at_zero
        end if

        if (full_iteration) then
            R_multipole_mpb_at_zero(:) = 0.d0
            R_multipole_spl_mpb = 0d0
            multipole_moment_ff_glob = 0d0
            multipole_moment_ff = 0d0
        else
            R_multipole_mpb_at_zero(:) = 0.d0
        end if

        if (mpb_forces_on) then
            R_multipole_mpb_deriv_at_zero(:,:) = 0.d0
        end if


        do i_atom = 1, n_atoms   
        
            if (atomic_MERM .and. atom_converged(i_atom)) then
                cycle
            end if
        
            ! set end of multipole expansion to zero explicitly, - to ensure that there is no diverging charge density
            ! outside the known grid quantities - this should be done on all atoms; this is done again on the finally splined multipole density later
            R_multipole_mpb(n_radial(species(i_atom))+1,:,i_atom) = 0.d0

            if (set_rho_mp) then
                rho_multipole_mpb(n_radial(species(i_atom))+2,:,i_atom) = 0.d0

                ! At zero (only for rho_multipole_mpb)
                !this must be done to obtain a smooth spline interpolation of the density at the zero
                !in Kerker preconditioner this is later done for the integration weight, but not for the splines
                !so we save the corrected spline interpolation as rho_multipole_mpb
                drho = (rho_multipole_mpb(2,1,i_atom) - rho_multipole_mpb(3,1,i_atom)) &
                    /    (r_radial(1,species(i_atom)) - r_radial(2,species(i_atom)))

                rho_multipole_mpb(1,1,i_atom) = rho_multipole_mpb(2,1,i_atom) &
                    - drho * r_radial(1,species(i_atom))

                rho_multipole_mpb(1,2:limit_l_dim_SPE_species(species(i_atom)),i_atom) = 0.d0
            end if

            ! atomic task list distribution - why work with atomic tasks
            if (myid.eq.task_list(i_atom)) then

                ! spline_atom_storage(n_atoms) is defined on each thread and contains the order of atoms which are treated by this
                ! particular thread - in case that distributed spline_storage is required, the (atomic) index for a spline
                ! in the spline storage array is given by spline_atom_storage(i_atom)
                current_spl_atom = spline_atom_storage(i_atom)

                ! do the splining itself: separately, for each atom, and each angular momentum shell.
                ! the treatment of all radial shells is done in the splining routine

                do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                    do i_m = - i_l, i_l, 1

                        call cubic_spline &
                            (R_multipole_mpb(:,index_lm_mpb(i_m, i_l),i_atom), &
                            n_radial(species(i_atom))+1, &
                            aux_R_spline_mpb )

                        ! copy the array of spline coefficients to its eventual location ...
                        do i_radial = 1, n_radial(species(i_atom)), 1
                            do i_spline = 1, n_max_spline_mpb, 1
                                if (evaluate_only_mpe) then
                                    R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), i_spline, i_radial, current_spl_atom) &
                                        = aux_R_spline_mpb(i_spline,i_radial)
                                else
                                    R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), i_spline, i_radial+1, current_spl_atom) &
                                        = aux_R_spline_mpb(i_spline,i_radial)
                                end if
                            enddo
                        enddo
                    end do
                end do              ! end splining loop


                ! doctor the atomic splines such that they don't yield infinite residuals in the far field - just in case ...
                ! find outermost radial grid point that is possibly non-zero
                i_radial = n_radial(species(i_atom))

                if (evaluate_only_mpe) then
                    do while ( ( r_radial(i_radial,species(i_atom)) .ge. &
                            multipole_radius_free_SPE(species(i_atom)) ) &
                            .and.(i_radial.gt.1) )
                        R_multipole_spl_mpb(:,:,i_radial-1,current_spl_atom) = 0.d0    !SR: changes this to i_radial+1 set whole spline to zero for non-valid points
                        i_radial = i_radial - 1                                        ! next innermost radial point
                    enddo
                else
                    do while ( ( r_radial(i_radial,species(i_atom)) .ge. &
                            multipole_radius_free_SPE(species(i_atom)) ) &
                            .and.(i_radial.gt.1) )
                        R_multipole_spl_mpb(:,:,i_radial,current_spl_atom) = 0.d0    !SR: changes this to i_radial+1 set whole spline to zero for non-valid points
                        i_radial = i_radial - 1                                      ! next innermost radial point
                    enddo
                end if

                ! Outermost atom radius in units of the radial integration grid  - this is for the free atom ...
                ! FH: I did not go through the following spline-doctoring procedure in detail, but copied it from update_hartree_potential_p1 -
                !     trusting that it was tested there in detail.
                i_r_outer = invert_radial_grid &
                    ( multipole_radius_free_SPE(species(i_atom)), &
                    n_radial(species(i_atom)), &
                    scale_radial(species(i_atom)) ) + 1

                ! difference from outermost density point to outermost free atom point
                delta = i_r_outer - i_radial
                delta_2 = delta*delta
                delta_3 = delta_2*delta

                ! doctor the spline coefficients at the outermost finite value
                ! i_radial to go smoothly to zero at multipole_radius_free
                !SR: changed i_radial to i_radial+1
                if (evaluate_only_mpe) then
                    do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                        do i_m = - i_l, i_l, 1
                            R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), 3, i_radial-1, current_spl_atom ) = &
                                - 3.d0 / delta_2 * &
                                R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), 1, i_radial-1, current_spl_atom ) &
                                - 2.d0 / delta * &
                                R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), 2, i_radial-1, current_spl_atom )
                            R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), 4, i_radial-1, current_spl_atom ) = &
                                2.d0 / delta_3 * &
                                R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), 1, i_radial-1, current_spl_atom ) &
                                + 1.d0 / delta_2 * &
                                R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), 2, i_radial-1, current_spl_atom )
                        enddo
                    enddo
                else
                    do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                        do i_m = - i_l, i_l, 1
                            R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), 3, i_radial, current_spl_atom ) = &
                                - 3.d0 / delta_2 * &
                                R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), 1, i_radial, current_spl_atom ) &
                                - 2.d0 / delta * &
                                R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), 2, i_radial, current_spl_atom )
                            R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), 4, i_radial, current_spl_atom ) = &
                                2.d0 / delta_3 * &
                                R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), 1, i_radial, current_spl_atom ) &
                                + 1.d0 / delta_2 * &
                                R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), 2, i_radial, current_spl_atom )
                        enddo
                    enddo
                end if

                if (.not.evaluate_only_mpe) then


                    ! figure out parabolic splines for the inside of r_radial(2,species(i_atom))
                    ! These are REAL-SPACE splines that are INDEPENDENT of any mapping of an index/radius onto a
                    ! nonlinear grid and are there to take out a divergence of the radial grid mapping near the
                    ! origin.

                    ! all the splines are written as parabolas to make evaluation a lot simpler
                    inner_radial_spl_mpb = 0d0
                    ! calculate the inner spline for l = 0, m = 0: linear, going through f(r1) and f(r2)
                    r1 = r_radial(1,species(i_atom))
                    r2 = r_radial(2,species(i_atom))
                    f1 = R_multipole_mpb(1,index_lm_mpb(0,0),i_atom)
                    f2 = R_multipole_mpb(2,index_lm_mpb(0,0),i_atom)
                    inner_radial_spl_mpb(2,index_lm_mpb(0,0)) = (f2-f1)/(r2-r1)     ! linear term
                    inner_radial_spl_mpb(1,index_lm_mpb(0,0)) = f1-r1*inner_radial_spl_mpb(2, index_lm_mpb(0,0)) ! constant term

                    ! calculate the inner spline for l > 0, m - parabolic, going through f(r1), f(r2), and (0,0)
                    denominator = (r1*r1*r2-r2*r2*r1)
                    do i_l = 1, limit_l_max_SPE_species(species(i_atom)), 1
                        do i_m = -i_l, i_l
                            f1 = R_multipole_mpb(1, index_lm_mpb(i_m, i_l), i_atom)
                            f2 = R_multipole_mpb(2, index_lm_mpb(i_m, i_l), i_atom)
                            ! constant term = 0
                            inner_radial_spl_mpb(2, index_lm_mpb(i_m, i_l)) = (r1*r1*f2-r2*r2*f1)/denominator  ! linear term
                            inner_radial_spl_mpb(3, index_lm_mpb(i_m, i_l)) = (r2*f1-f2*r1)/denominator        ! parabolic term
                        end do
                    end do

                    ! integrate the Green function integral (12b) in Delley, 1990; with the modified Green functions I_{l+1/2}(q_0 r) and K_{l+1/2}(q_0 r)
                    ! for each atom; for now on the logarithmic grid and using the splines that were just computed.
                    ! what follows is also done in the routine integrate_hartree_log_grid

                    ! FIXME: why do we even have to store all the splines? In case that this is a memory bottle-neck, there might be a large scale
                    !        workaround by working on each (l,m) value separately, store the temporary splines and then continue to evaluate the
                    !        completed R_prec(alpha,l,m,r) ... which would then have to be synchronized again, of course

                    ! FIXME: how are the speed, memory requirements, and accuracy of the overall result affected by the calculation of the
                    !        splines vs the use of the regular radial grids for the integrations? The latter should be much faster but should also
                    !        yield a lot less accurate results. Does that matter for a preconditioning correction to the density????


                    angular_integral_log_mpb = 0d0
                    alpha = log(r_grid_inc(species(i_atom)))   ! prefactor for the calculation of the grid radius & the integration weight

                    ! calculate the indexes for integration limits
                    ! n_inner_grid = the point where the second radial shell starts, i.e. where we go to the regular splines.
                    n_inner_grid = invert_log_grid(r_radial(2,species(i_atom)), &
                        r_grid_min(species(i_atom)),r_grid_inc(species(i_atom)))

                    ! n_grid_limit = the point after which there is no density anyway, i.e. the free atom multipole radius
                    n_grid_limit = invert_log_grid(multipole_radius_free_SPE(species(i_atom))+extra_adding_to_hartree_potential_distance, &
                        r_grid_min(species(i_atom)),r_grid_inc(species(i_atom)))

                    ! first spline evaluation loop:
                    ! everything inside the innermost two radial shells.
                    do i_grid = 1, n_inner_grid

                        ! this is the actual radius known to the parabolic splines !
                        radius = r_grid(i_grid,species(i_atom))
                        radius_sq = radius*radius

                        ! evaluate real-space parabolic spline for all (l,m)
                        angular_integral_log_mpb(:, i_grid) = &
                            inner_radial_spl_mpb(1,:)+inner_radial_spl_mpb(2,:)*radius+inner_radial_spl_mpb(3,:)*radius_sq

                        ! transform log grid onto uniform index space for integration
                        angular_integral_log_mpb(:,i_grid) = alpha*(radius**3)*angular_integral_log_mpb(:,i_grid)
                    end do




                    ! second spline evaluation loop on logarithmic grid:
                    ! everything for which we know to have well-defined splines.
                    do i_grid = n_inner_grid, n_grid_limit + 1


                        ! calculate index of grid point on radial grid, where the splines reside
                        radius = invert_radial_grid( r_grid(i_grid,species(i_atom)), n_radial(species(i_atom)), &
                            scale_radial(species(i_atom)))+1d0

                        ! evaluate all (l,m) splines for this particular grid point    - arguments in the declaration of spline_vector ...
                        call spline_vector( radius,                         & !  real*8  :: r_output
                            R_multipole_spl_mpb(:,:,:,current_spl_atom),    & !  real*8  :: spl_param(n_l_dim,4,n_grid_dim) - for this particular atom only
                            n_max_radial+1,                                 & !  integer :: n_grid_dim
                            limit_l_dim_SPE,                                & !  integer :: n_l_dim
                            n_radial(species(i_atom))+1,                    & !  integer :: n_points
                            limit_l_dim_SPE_species(species(i_atom)),       & !  integer :: n_vector -- the number of points actually evaluated
                            angular_integral_log_mpb(:,i_grid))               !  real*8  :: out_result(n_vector)

                        ! transform log grid onto uniform index space for integration
                        angular_integral_log_mpb(:,i_grid) = alpha*(r_grid(i_grid,species(i_atom))**3)*angular_integral_log_mpb(:,i_grid)
                    end do   ! end evaluation of splines across the whole grid, for atom current_spl_atom

                    !SR: correct the angular_integral_log_mpb. i observed this factor to be necessary to get the right results...
                    if (use_partition) then
                        angular_integral_log_mpb = angular_integral_log_mpb * 8*kappa_SPE_step
                    else
                        angular_integral_log_mpb = angular_integral_log_mpb * 8*kappa_SPE
                    end if

                    ! we now know, for i_atom, the values of the splined multipole residual R_in on the logarithmic grid
                    !                                            - let's calculate R_prec_multipole_mpb on the same logarithmic grid.
                    ! This is done using the Adams-Moulton linear multistep integrator, using up to 4 terms.
                    !    see http://en.wikipedia.org/wiki/Linear_multistep_method or Abramowitz/Stegun p. 896 for details.

                    ! as in integrate_hartree_log_grid, do so in two parts
                    ! part (1) contains the integral from 0 up to the grid radius r, with I_{l+1/2} in the integrand and 4 pi K_{l+1/2} as the prefactor
                    ! part (2) is the other way around, i.e. from r to "infinity"    with K_{l+1/2} in the integrand and 4 pi I_{l+1/2} as the prefactor

                    ! first integral: from zero to radius r, for all possible r

                    ! first three terms: do one-by-one with increasing order of integration method.


                    ! 	integrate_full_grid = .false. !tested for small molecules, no effect in solvation energy
                    ! 	!it means for an h-atom, e.g., to integrate until 13.33 instead of 99.0, but at 13.33 the electron
                    ! 	!density is already zero anyway, so it does not matter
                    !
                    ! 	if (integrate_full_grid) then
                    !
                    ! 	  integral_zero_r_mpb  = 0d0
                    ! 	  do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                    ! 	    do i_m = - i_l, i_l, 1
                    !
                    ! 	      ! TERM 1 : Integral_1 = h*f_1; but h = 1
                    ! 	      integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) =  &
                    ! 		    angular_integral_log_mpb(index_lm_mpb(i_m, i_l),1) * green_I_mpb(1,i_l+1,species(i_atom))
                    ! 	      R_prec_multipole_mpb(1,index_lm_mpb(i_m,i_l)) =  &
                        ! 		    integral_zero_r_mpb(index_lm_mpb(i_m, i_l) ) * green_K_mpb(1,i_l+1,species(i_atom))
                    !
                    ! 	      ! TERM 2 : Integral_2 = Integral_1 + h(f_2+f_1)/2
                    ! 	      integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) =  integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) &
                    ! 		  + ( 1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),2) * green_I_mpb(2,i_l+1,species(i_atom)) &
                        ! 		    + 1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),1) * green_I_mpb(1,i_l+1,species(i_atom)))/2d0
                    ! 	      R_prec_multipole_mpb(2,index_lm_mpb(i_m,i_l)) =  &
                        ! 		    integral_zero_r_mpb(index_lm_mpb(i_m, i_l) ) * green_K_mpb(2,i_l+1,species(i_atom))
                    !
                    ! 	      ! TERM 3 : Integral_3 = Integral_2 + h(5f_3 + 8f_2 - f_1)/12
                    ! 	      integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) =  integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) &
                    ! 		  + ( 5d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),3) * green_I_mpb(3,i_l+1,species(i_atom)) &
                        ! 		    + 8d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),2) * green_I_mpb(2,i_l+1,species(i_atom)) &
                        ! 		    - 1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),1) * green_I_mpb(1,i_l+1,species(i_atom)))/12d0
                    ! 	      R_prec_multipole_mpb(3,index_lm_mpb(i_m,i_l)) =  &
                        ! 		    integral_zero_r_mpb(index_lm_mpb(i_m, i_l) ) * green_K_mpb(3,i_l+1,species(i_atom))
                    ! 	    end do
                    ! 	  end do
                    !
                    ! 	  ! TERM i_grid > 4 : Integral_i = Integral_(i-1) + h[9 f_i + 19 f_(i-1) - 5 f_(i-2) + f_(i-3)]/24
                    ! 	  do i_grid = 4, n_grid(species(i_atom)), 1
                    !
                    ! ! 	    if ( r_grid(i_grid,species(i_atom)) .lt. &
                        ! ! 		  (multipole_radius_free(species(i_atom))+extra_adding_to_hartree_potential_distance)) then
                    ! 		!       radial integration weight on the logarithmic grid alpha*r times usual radial
                    ! 		!       integration weight from integral r^2 dr
                    !
                    ! 		! loop over all angular momentum components
                    ! 		do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                    ! 		  do i_m = - i_l, i_l, 1
                    !
                    ! 		      integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) =   integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) &
                        ! 			+ (  9d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid  )*green_I_mpb(i_grid  ,i_l+1,species(i_atom)) &
                        ! 			  + 19d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid-1)*green_I_mpb(i_grid-1,i_l+1,species(i_atom)) &
                        ! 			  -  5d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid-2)*green_I_mpb(i_grid-2,i_l+1,species(i_atom)) &
                        ! 			  +  1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid-3)*green_I_mpb(i_grid-3,i_l+1,species(i_atom)))/24d0
                    !
                    ! 		      R_prec_multipole_mpb(i_grid,index_lm_mpb(i_m,i_l)) =  &                        ! spice with outside constant Green function
                    ! 			  integral_zero_r_mpb(index_lm_mpb(i_m, i_l) ) * green_K_mpb(i_grid,i_l+1,species(i_atom))
                    !
                    ! 		  end do
                    ! 		end do
                    ! ! 	    else
                    ! ! 		do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                    ! ! 		  do i_m = - i_l, i_l, 1
                    ! ! 		      R_prec_multipole_mpb(i_grid,index_lm_mpb(i_m,i_l)) =  &
                    ! ! 			  integral_zero_r_mpb(index_lm_mpb(i_m, i_l) ) * green_K_mpb(i_grid,i_l+1,species(i_atom))
                    ! ! 		  end do
                    ! ! 		end do
                    ! ! 	    end if
                    ! 	  end do        ! end radial loop for calculation of first integral
                    !
                    ! 	  ! second integral: from radius r to infinity - well, do it backwards and start by figuring out where the
                    ! 	  !    integrand starts being non-zero, i.e. where we have to start worrying about it.
                    ! 	  n_grid_limit = n_grid(species(i_atom))
                    ! 	  do while ( r_grid(n_grid_limit,species(i_atom)) .gt. &
                        ! 		    (multipole_radius_free_SPE(species(i_atom))+extra_adding_to_hartree_potential_distance))
                    ! 	    n_grid_limit = n_grid_limit - 1
                    ! 	  end do
                    !
                    ! 	  ! start integrating from the outside in, again using the Adams-Moulton linear multistep integrator.
                    ! 	  ! the first three terms warrant special treatment, similar to the above integral.
                    ! 	  integral_r_infty_mpb = 0d0
                    ! 	  do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                    ! 	    do i_m = - i_l, i_l, 1
                    !
                    ! 		! TERM 1 : Integral_N = h*f_N; but h = 1
                    ! 		integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) = &
                    ! 		    angular_integral_log_mpb(index_lm_mpb(i_m, i_l),n_grid(species(i_atom))) * green_K_mpb(n_grid(species(i_atom)),i_l+1,species(i_atom))
                    ! 		R_prec_multipole_mpb(n_grid(species(i_atom)),index_lm_mpb(i_m,i_l)) =   &
                    ! 		    R_prec_multipole_mpb(n_grid(species(i_atom)),index_lm_mpb(i_m,i_l)) &
                    ! 		    + integral_r_infty_mpb(index_lm_mpb(i_m, i_l) ) * green_I_mpb(n_grid(species(i_atom)),i_l+1,species(i_atom))
                    !
                    ! 		! TERM 2 : Integral_(N-1) = Integral_N + h(f_(N-1)+f_N)/2
                    ! 		integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) = integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) + &
                    ! 		    ( 1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),n_grid(species(i_atom))-1)*green_K_mpb(n_grid(species(i_atom))-1,i_l+1,species(i_atom)) &
                        ! 		    + 1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),n_grid(species(i_atom))  )*green_K_mpb(n_grid(species(i_atom))  ,i_l+1,species(i_atom)))/2d0
                    ! 		R_prec_multipole_mpb(n_grid(species(i_atom))-1,index_lm_mpb(i_m,i_l)) =   &
                    ! 		    R_prec_multipole_mpb(n_grid(species(i_atom))-1,index_lm_mpb(i_m,i_l)) &
                    ! 		    + integral_r_infty_mpb(index_lm_mpb(i_m, i_l) ) * green_I_mpb(n_grid(species(i_atom))-1,i_l+1,species(i_atom))
                    !
                    ! 		! TERM 3 : Integral_(N-2) = Integral_(N-1) + h(5f_(N-2) + 8f_(N-1) - f_N)/12
                    ! 		integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) = integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) + &
                    ! 		    ( 5d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),n_grid(species(i_atom))-2)*green_K_mpb(n_grid(species(i_atom))-2,i_l+1,species(i_atom)) &
                        ! 		    + 8d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),n_grid(species(i_atom))-1)*green_K_mpb(n_grid(species(i_atom))-1,i_l+1,species(i_atom)) &
                        ! 		    - 1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),n_grid(species(i_atom))  )*green_K_mpb(n_grid(species(i_atom))  ,i_l+1,species(i_atom)))/12d0
                    ! 		R_prec_multipole_mpb(n_grid(species(i_atom))-2,index_lm_mpb(i_m,i_l)) =   &
                    ! 		    R_prec_multipole_mpb(n_grid(species(i_atom))-2,index_lm_mpb(i_m,i_l)) &
                    ! 		    + integral_r_infty_mpb(index_lm_mpb(i_m, i_l) ) * green_I_mpb(n_grid(species(i_atom))-2,i_l+1,species(i_atom))
                    ! 	    end do
                    ! 	  end do
                    !
                    ! 	  ! all remaining terms
                    ! 	  ! Integral_i = Integral_(i+1) + h[9 f_i + 19 f_(i+1) - 5 f_(i+2) + f_(i+3)]/24
                    ! 	  do i_grid = n_grid(species(i_atom))-3, 1, -1
                    ! 	    do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                    ! 		do i_m = - i_l, i_l, 1
                    ! 		integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) = integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) + &
                    ! 		    (  9d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid  )*green_K_mpb(i_grid  ,i_l+1,species(i_atom)) &
                        ! 		    + 19d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid+1)*green_K_mpb(i_grid+1,i_l+1,species(i_atom)) &
                        ! 		    -  5d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid+2)*green_K_mpb(i_grid+2,i_l+1,species(i_atom)) &
                        ! 		    +  1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid+3)*green_K_mpb(i_grid+3,i_l+1,species(i_atom)))/24d0
                    ! 		  R_prec_multipole_mpb(i_grid,index_lm_mpb(i_m,i_l)) =   &
                    ! 			R_prec_multipole_mpb(i_grid,index_lm_mpb(i_m,i_l)) &
                    ! 			+ integral_r_infty_mpb(index_lm_mpb(i_m, i_l) ) * green_I_mpb(i_grid,i_l+1,species(i_atom))
                    ! 		end do
                    ! 	    end do
                    ! 	  end do ! end calculation of second integral
                    !
                    ! 	else !integrate just to multipole_radius+extra_adding_to_hartree_potential_distance

                    integral_zero_r_mpb  = 0d0
                    do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                        do i_m = - i_l, i_l, 1

                            ! TERM 1 : Integral_1 = h*f_1; but h = 1
                            integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) =  &
                                angular_integral_log_mpb(index_lm_mpb(i_m, i_l),1) * green_I_mpb(1,i_l+1,species(i_atom))
                            R_prec_multipole_mpb(1,index_lm_mpb(i_m,i_l)) =  &
                                integral_zero_r_mpb(index_lm_mpb(i_m, i_l) ) * green_K_mpb(1,i_l+1,species(i_atom))

                            ! TERM 2 : Integral_2 = Integral_1 + h(f_2+f_1)/2
                            integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) =  integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) &
                                + ( 1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),2) * green_I_mpb(2,i_l+1,species(i_atom)) &
                                + 1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),1) * green_I_mpb(1,i_l+1,species(i_atom)))/2d0
                            R_prec_multipole_mpb(2,index_lm_mpb(i_m,i_l)) =  &
                                integral_zero_r_mpb(index_lm_mpb(i_m, i_l) ) * green_K_mpb(2,i_l+1,species(i_atom))

                            ! TERM 3 : Integral_3 = Integral_2 + h(5f_3 + 8f_2 - f_1)/12
                            integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) =  integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) &
                                + ( 5d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),3) * green_I_mpb(3,i_l+1,species(i_atom)) &
                                + 8d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),2) * green_I_mpb(2,i_l+1,species(i_atom)) &
                                - 1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),1) * green_I_mpb(1,i_l+1,species(i_atom)))/12d0
                            R_prec_multipole_mpb(3,index_lm_mpb(i_m,i_l)) =  &
                                integral_zero_r_mpb(index_lm_mpb(i_m, i_l) ) * green_K_mpb(3,i_l+1,species(i_atom))
                        end do
                    end do


                    ! TERM i_grid > 4 : Integral_i = Integral_(i-1) + h[9 f_i + 19 f_(i-1) - 5 f_(i-2) + f_(i-3)]/24
                    do i_grid = 4, n_grid(species(i_atom)), 1

                        if ( r_grid(i_grid,species(i_atom)) .lt. &
                                (multipole_radius_free_SPE(species(i_atom))+extra_adding_to_hartree_potential_distance)) then
                            !       radial integration weight on the logarithmic grid alpha*r times usual radial
                            !       integration weight from integral r^2 dr

                            ! loop over all angular momentum components
                            do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                                do i_m = - i_l, i_l, 1

                                    integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) =   integral_zero_r_mpb(index_lm_mpb(i_m, i_l)) &
                                        + (  9d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid  )*green_I_mpb(i_grid  ,i_l+1,species(i_atom)) &
                                        + 19d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid-1)*green_I_mpb(i_grid-1,i_l+1,species(i_atom)) &
                                        -  5d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid-2)*green_I_mpb(i_grid-2,i_l+1,species(i_atom)) &
                                        +  1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid-3)*green_I_mpb(i_grid-3,i_l+1,species(i_atom)))/24d0

                                    R_prec_multipole_mpb(i_grid,index_lm_mpb(i_m,i_l)) =  &                        ! spice with outside constant Green function
                                        integral_zero_r_mpb(index_lm_mpb(i_m, i_l) ) * green_K_mpb(i_grid,i_l+1,species(i_atom))

                                end do
                            end do
                        else
                            !integral_zero_r_mpb will not change anymore, because multipole moments are zero for larger radii (no more contributions to integral)
                            do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                                do i_m = - i_l, i_l, 1
                                    R_prec_multipole_mpb(i_grid,index_lm_mpb(i_m,i_l)) =  &
                                        integral_zero_r_mpb(index_lm_mpb(i_m, i_l) ) * green_K_mpb(i_grid,i_l+1,species(i_atom))
                                end do
                            end do
                            if (use_analytic_far_field) then
                                do i_l = 0, limit_l_max_SPE_species_ff(species(i_atom)), 1
                                    do i_m = - i_l, i_l, 1
                                        multipole_moment_ff(current_spl_atom,index_lm_mpb(i_m,i_l)) =&
                                            integral_zero_r_mpb(index_lm_mpb(i_m, i_l))
                                    end do
                                end do
                            end if

                        end if
                    end do        ! end radial loop for calculation of first integral


                    ! second integral: from radius r to infinity - well, do it backwards and start by figuring out where the
                    !    integrand starts being non-zero, i.e. where we have to start worrying about it.
                    n_grid_limit = n_grid(species(i_atom))

                    do while ( r_grid(n_grid_limit,species(i_atom)) .gt. &
                            (multipole_radius_free_SPE(species(i_atom))+extra_adding_to_hartree_potential_distance))
                        n_grid_limit = n_grid_limit - 1
                    end do

                    ! start integrating from the outside in, again using the Adams-Moulton linear multistep integrator.
                    ! the first three terms warrant special treatment, similar to the above integral.
                    integral_r_infty_mpb = 0d0
                    do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                        do i_m = - i_l, i_l, 1

                            ! TERM 1 : Integral_N = h*f_N; but h = 1
                            integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) = &
                                angular_integral_log_mpb(index_lm_mpb(i_m, i_l),n_grid_limit) * green_K_mpb(n_grid_limit,i_l+1,species(i_atom))
                            R_prec_multipole_mpb(n_grid_limit,index_lm_mpb(i_m,i_l)) =   &
                                R_prec_multipole_mpb(n_grid_limit,index_lm_mpb(i_m,i_l)) &
                                + integral_r_infty_mpb(index_lm_mpb(i_m, i_l) ) * green_I_mpb(n_grid_limit,i_l+1,species(i_atom))

                            ! TERM 2 : Integral_(N-1) = Integral_N + h(f_(N-1)+f_N)/2
                            integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) = integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) + &
                                ( 1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),n_grid_limit-1)*green_K_mpb(n_grid_limit-1,i_l+1,species(i_atom)) &
                                + 1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),n_grid_limit  )*green_K_mpb(n_grid_limit  ,i_l+1,species(i_atom)))/2d0
                            R_prec_multipole_mpb(n_grid_limit-1,index_lm_mpb(i_m,i_l)) =   &
                                R_prec_multipole_mpb(n_grid_limit-1,index_lm_mpb(i_m,i_l)) &
                                + integral_r_infty_mpb(index_lm_mpb(i_m, i_l) ) * green_I_mpb(n_grid_limit-1,i_l+1,species(i_atom))

                            ! TERM 3 : Integral_(N-2) = Integral_(N-1) + h(5f_(N-2) + 8f_(N-1) - f_N)/12
                            integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) = integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) + &
                                ( 5d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),n_grid_limit-2)*green_K_mpb(n_grid_limit-2,i_l+1,species(i_atom)) &
                                + 8d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),n_grid_limit-1)*green_K_mpb(n_grid_limit-1,i_l+1,species(i_atom)) &
                                - 1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),n_grid_limit  )*green_K_mpb(n_grid_limit  ,i_l+1,species(i_atom)))/12d0
                            R_prec_multipole_mpb(n_grid_limit-2,index_lm_mpb(i_m,i_l)) =   &
                                R_prec_multipole_mpb(n_grid_limit-2,index_lm_mpb(i_m,i_l)) &
                                + integral_r_infty_mpb(index_lm_mpb(i_m, i_l) ) * green_I_mpb(n_grid_limit-2,i_l+1,species(i_atom))
                        end do
                    end do

                    ! all remaining terms
                    ! Integral_i = Integral_(i+1) + h[9 f_i + 19 f_(i+1) - 5 f_(i+2) + f_(i+3)]/24
                    do i_grid = n_grid_limit-3, 1, -1
                        do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                            do i_m = - i_l, i_l, 1
                                integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) = integral_r_infty_mpb(index_lm_mpb(i_m, i_l)) + &
                                    (  9d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid  )*green_K_mpb(i_grid  ,i_l+1,species(i_atom)) &
                                    + 19d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid+1)*green_K_mpb(i_grid+1,i_l+1,species(i_atom)) &
                                    -  5d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid+2)*green_K_mpb(i_grid+2,i_l+1,species(i_atom)) &
                                    +  1d0*angular_integral_log_mpb(index_lm_mpb(i_m, i_l),i_grid+3)*green_K_mpb(i_grid+3,i_l+1,species(i_atom)))/24d0
                                R_prec_multipole_mpb(i_grid,index_lm_mpb(i_m,i_l)) =   &
                                    R_prec_multipole_mpb(i_grid,index_lm_mpb(i_m,i_l)) &
                                    + integral_r_infty_mpb(index_lm_mpb(i_m, i_l) ) * green_I_mpb(i_grid,i_l+1,species(i_atom))
                            end do
                        end do
                    end do ! end calculation of second integral
                    ! 	end if !if full integration or just integration to multipole_radius



                    !calculate monopole potential at r = 0

                    n_grid_limit = n_grid(species(i_atom))
                    do while ( r_grid(n_grid_limit,species(i_atom)) .gt. &
                            (multipole_radius_free_SPE(species(i_atom))+extra_adding_to_hartree_potential_distance))
                        n_grid_limit = n_grid_limit - 1
                    end do



                    !evaluate multipole moments and derivatives at radius = zero
                    !in case of the moments, only the monopoles, for the derivatives only the dipoles give non-zero
                    !contributions
                    i_l = 1
                    do  i_grid = 1,  n_grid_limit,1 !n_grid_limit, 1 !n_grid(species(i_atom)),1

                        R_multipole_mpb_at_zero(current_spl_atom) = &
                            R_multipole_mpb_at_zero(current_spl_atom) &
                            +  angular_integral_log_mpb(index_lm_mpb(0, 0),i_grid ) * green_K_mpb(i_grid ,1,species(i_atom))

                        if (mpb_forces_on) then
                            do i_m = -i_l, i_l, 1
                                R_multipole_mpb_deriv_at_zero(i_m+2,current_spl_atom) = &
                                    R_multipole_mpb_deriv_at_zero(i_m+2,current_spl_atom) +&
                                    angular_integral_log_mpb(index_lm_mpb(i_m,i_l),i_grid) *&
                                    green_K_mpb(i_grid,i_l+1,species(i_atom)) !SR replaced green_K_mpb(i_grid,3,species(i_atom))
                                !		!multiply the dipole moments with kappa/3. The 1/3 stems from the derivative
                                !		!of the Bessel I-function at zero, the kappa from the inner derivative
                                !		  if (use_partition) then
                                !		    R_multipole_mpb_deriv_at_zero(i_m+2,current_spl_atom) = &
                                    !		      dble(1d0/3d0)*kappa_SPE_step*&
                                    !		      R_multipole_mpb_deriv_at_zero(i_m+2,current_spl_atom)
                                !		  else
                                !		    R_multipole_mpb_deriv_at_zero(i_m+2,current_spl_atom) = &
                                    !		      dble(1d0/3d0)*kappa_SPE*&
                                    !		      R_multipole_mpb_deriv_at_zero(i_m+2,current_spl_atom)
                                !		  end if
                                !			green_K_mpb(i_grid,2,species(i_atom)) !SR replaced green_K_mpb(i_grid,3,species(i_atom))

                            end do
                        end if
                    end do
                    !multiply the dipole moments with kappa/3. The 1/3 stems from the derivative
                    !of the Bessel I-function at zero, the kappa from the inner derivative
                    if (mpb_forces_on) then
                        do i_m = -i_l, i_l, 1
                            if (use_partition) then
                                R_multipole_mpb_deriv_at_zero(i_m+2,current_spl_atom) = &
                                    dble(1d0/3d0)*kappa_SPE_step*&
                                    R_multipole_mpb_deriv_at_zero(i_m+2,current_spl_atom)
                            else
                                R_multipole_mpb_deriv_at_zero(i_m+2,current_spl_atom) = &
                                    dble(1d0/3d0)*kappa_SPE*&
                                    R_multipole_mpb_deriv_at_zero(i_m+2,current_spl_atom)
                            end if
                        end do
                    end if

                    !	forces
                    ! 	  write(use_unit,*) 'R_multipole_mpb_deriv_at_zero'
                    ! 	    do i_grid=1,3,1
                    ! 	    write(use_unit,*) 'R_m', current_spl_atom, R_multipole_mpb_deriv_at_zero(current_spl_atom,i_grid)
                    ! 	    end do


                    !		forces
                    ! 	  write(use_unit,*) '1st part, R_m1', R_multipole_mpb_deriv_at_zero(:,1)
                    !           write(use_unit,*) '1st part, R_m2', R_multipole_mpb_deriv_at_zero(:,2)
                    ! 	  write(use_unit,*) '1st part, R_m3', R_multipole_mpb_deriv_at_zero(:,3)
                    !
                    ! 	  write(use_unit,*) 'kappa_SPE',kappa_SPE


                    ! 	green_I_mpb(1 ,1,species(i_atom)) is exactly one


                    ! set end of multipole expansion to zero explicitly, - to ensure that there is no diverging charge density
                    ! outside the known grid quantities - this should be done on all atoms; this is done again on the finally splined multipole density later
                    !         R_prec_multipole_mpb(n_grid(species(i_atom))+1,:) = 0.d0 !BUG solved, changed n_radial to n_grid

                    ! 	do i_grid=1, n_grid(species(i_atom)),1
                    ! 	  write(use_unit,*) 'integralout',i_grid,R_prec_multipole_mpb(i_grid,1)
                    ! 	end do
                    ! 	stop
                    ! spline_atom_storage(n_atoms) is defined on each thread and contains the order of atoms which are treated by this
                    ! particular thread - in case that distributed spline_storage is required, the (atomic) index for a spline
                    ! in the spline storage array is given by spline_atom_storage(i_atom)

                    ! do the splining itself: separately, for each atom, and each angular momentum shell.
                    ! the treatment of all radial shells is done in the splining routine


                    do i_l = 0, limit_l_max_SPE_species(species(i_atom)), 1
                        do i_m = - i_l, i_l, 1
                            !DEBUG: plot multipole moments
                            ! 	      write(use_unit,*) 'multipolemoment ', i_atom, i_l, i_m, sum(abs(R_prec_multipole_mpb(:,index_lm_mpb(i_m, i_l))))

                            aux_R_prec_spline_mpb = 0d0
                            call cubic_spline &
                                (R_prec_multipole_mpb(:,index_lm_mpb(i_m, i_l)), &
                                n_grid(species(i_atom)), &
                                aux_R_prec_spline_mpb )


                            ! copy the array of spline coefficients to its eventual location ...
                            do i_radial = 1, n_grid(species(i_atom)), 1
                                do i_spline = 1, n_max_spline_mpb, 1
                                    R_multipole_spl_mpb( index_lm_mpb(i_m, i_l), i_spline, i_radial, current_spl_atom) &
                                        = aux_R_prec_spline_mpb(i_spline,i_radial)
                                end do
                            end do
                        end do
                    end do               ! end splinig loop

                end if !if evaluate_only_mpe
            end if                 ! atomic task list distribution in spline loop
            !      write(use_unit,*) 'nloggrid', n_grid(species(i_atom))
            !      do i_grid = 1, n_grid(species(i_atom)), 1
            !       write(use_unit,*) 'Rprec',i_grid,r_grid(i_grid,species(i_atom)),R_prec_multipole_mpb(i_grid,1)
            !      end do
        end do                    ! spline and integration loop over all atoms

        !   if (.not.evaluate_only_mpe) then
        call sync_vector(R_multipole_mpb_at_zero,n_spline_atoms)
        if (use_analytic_far_field) then
            call sync_matrix(multipole_moment_ff,n_spline_atoms,limit_l_dim_SPE_ff)
        end if
        if (full_iteration) then
!             if (atomic_MERM) then
!                 multipole_moment_ff_glob(i_atom_MERM,:) = multipole_moment_ff(i_atom_MERM,:)
!             else
                multipole_moment_ff_glob = multipole_moment_ff
!             end if
        end if

        if (mpb_forces_on) then
            !  	call sync_vector(delta_v_MPB_multipole_component_at_zero,l_dim_SPE_public*n_atoms**2)
            !  	call sync_vector(delta_v_MPB_multipole_deriv_at_zero,l_dim_SPE_public*n_atoms**2)
            !  	call sync_vector(R_prec_gradient_at_zero,3*n_spline_atoms)
            call sync_matrix(R_multipole_mpb_deriv_at_zero,3,n_spline_atoms)
            ! !  	if (myid==1) write(use_unit,*) 'rprecatzero', sqrt(sum(R_prec_gradient_at_zero**2,DIM=1))
            !      end if
        end if



        ! ! 	SR DEBUG: output the splines...
        !     do i_lm = 1, limit_l_dim_SPE_species(species(1)), 1
        !       write(loopstr,'(I3.3)') i_lm
        !       filewrite = loopstr//"_mp_pot.dat"
        !       open(100,file=filewrite)
        ! !       write(100,*) multipole_moment_ff(1,i_lm)
        !       do i_radial = 1, n_radial(species(1)), 1
        ! 	write(100,*) r_radial(i_radial,species(1)),R_prec_multipole_mpb(i_radial,i_lm)
        !       end do
        ! !       close(100)
        ! ! !       call output_delta_v_MPB_step(multipole_output(:,i_lm,i_atom),i_lm+100*(i_atom-1))
        !     end do
        !     stop

        !grid was atomvisely distributed on cpus (atom1 = cpu1, atom2 = cpu2,...), now synchronize the vector

        ! if every cpu keeps its own version of the splines, distribute them to
        ! all the CPU's before doing much more ...
        if (.not.use_distributed_spline_storage) &
            call sync_kerker_multipole_splines(R_multipole_spl_mpb,limit_l_dim_SPE,n_max_spline_mpb,&
            n_spline_grid_dim,n_spline_atoms)

        if (atomic_MERM) then
            !add to the arrays the already calculated values from the converged atoms
            do i_atom = 1, n_spline_atoms, 1
                if (atom_converged(i_atom)) then
                    R_multipole_spl_mpb(:,:,:,i_atom) =&
                            R_multipole_spl_mpb_temp(:,:,:,i_atom)
                    multipole_moment_ff_glob(i_atom,:) = multipole_moment_ff_glob_temp(i_atom,:)
                    multipole_moment_ff(i_atom,:) = multipole_moment_ff_temp(i_atom,:)
                    R_multipole_mpb_at_zero(i_atom) = R_multipole_mpb_at_zero_temp(i_atom)
                    R_multipole_mpb_deriv_at_zero(:,i_atom) =&
                        R_multipole_mpb_deriv_at_zero_temp(:,i_atom)
                end if
            end do
        end if

    end subroutine evaluate_radial_integral

    !------------------------------------------------------------------------------
    !****f* SPE_solver/partition_tab_mpb
    !  NAME
    !    partition_tab_mpb
    !  SYNOPSIS
    integer function partition_tab_mpb(routine,i_full_points,take_part,full_iteration,only_outer,kappa_current)
    ! PURPOSE
    !if we solve the MPBE/LPBE iteratively, we do not have to update the potential on the whole grid_point
    !but only on the grid points, where the part of the input density which depends on
    !the potential (L_1 * delta_v) is actually non-zero. this routine
    !gives 1, if current point is in calculation domain and zero if not...
    !routine: = 0 if called in evaluate_density_mpe, = 1 if called in sum_up_mpe_potential
        use mpb_solver_utilities, only: h_function
        use synchronize_mpi
        implicit none

        integer, intent(in) :: routine, i_full_points
        logical, intent(in) :: only_outer
        integer, intent(in) :: take_part,full_iteration
        real*8, intent(in) :: kappa_current
        real*8 :: eps_cut = 1d-30
        real*8 :: kappa_cut = 1d-30
        real*8 :: ion_cut = 1d-8
        real*8 :: c_cut = 1d-20
        partition_tab_mpb = 0

        if (routine==0) then

            if (take_part==0.or.take_part==1) then
                !select the points for which the part of rho_iter that depends on the potential is zero anyway
                !(independent from the value of the potential)
                if (take_part==0.or.(take_part==1.and.(&
                        !1st case: c_mpb = 0. only points in dielectric transition region are non-zero
                        (c_mpb.le.c_cut.and.(dielec_func_mpb(i_full_points)-1d0).gt.eps_cut.and.&
                        (epsinf_mpb-dielec_func_mpb(i_full_points)).gt.eps_cut).or.&
                        !2nd case: c_mpb > 0. points that lie outside of the dielectric AND kappa cavity are zero
                        !	and can be neglected
                        (c_mpb.gt.c_cut.and.&
                        (sum(dielec_func_gradient_mpb(:,i_full_points)**2).gt.eps_cut.or.&
                        abs(kappa_current**2-kappa_debye**2).gt.kappa_cut))))) then
                    partition_tab_mpb = 1
                end if
            else if (take_part==2) then
                !this is for the MPBE solver:
                if (take_part==2.and.(&
                        (c_mpb.le.c_cut.and.(dielec_func_mpb(i_full_points)-1d0).gt.eps_cut.and.&
                        (epsinf_mpb-dielec_func_mpb(i_full_points)).gt.eps_cut).or.&
                        !if c>0, we can throw away points for which grad eps * grad v =0 AND hfunction-kappa =0
                        (c_mpb.gt.c_cut.and.(&
                        ((dielec_func_mpb(i_full_points)-1d0).gt.eps_cut.and.&
                        (epsinf_mpb-dielec_func_mpb(i_full_points)).gt.eps_cut).or.&
                        (abs(kappa_debye**2-h_function(i_full_points)/dielec_func_mpb(i_full_points)).gt.eps_cut)&
                        )))) then
                    partition_tab_mpb = 1
                end if
            end if
            !   !select the points for which the part of rho_iter that depends on the potential is zero anyway
            !   !(independent from the value of the potential)
            !     if (take_part==0.or.(take_part==1.and.(&
                !       !1st case: c_mpb = 0. only points in dielectric transition region are non-zero
            !       (c_mpb.le.c_cut.and.(dielec_func_mpb(i_full_points)-1d0).gt.eps_cut.and.&
                !       (epsinf_mpb-dielec_func_mpb(i_full_points)).gt.eps_cut).or.&
                !       !2nd case: c_mpb > 0. points that lie outside of the dielectric AND kappa cavity are zero
            !       !	and can be neglected
            !       (c_mpb.gt.c_cut.and.&
                !       (sum(dielec_func_gradient_mpb(:,i_full_points)**2).gt.eps_cut.or.&
                !       abs(kappa_current**2-kappa_debye**2).gt.kappa_cut)))).or.&
                !       !this is for the MPBE solver:
            !       (take_part==2.and.(&
                !       (c_mpb.le.c_cut.and.(dielec_func_mpb(i_full_points)-1d0).gt.eps_cut.and.&
                !       (epsinf_mpb-dielec_func_mpb(i_full_points)).gt.eps_cut).or.&
                !       !if c>0, we can throw away points for which grad eps * grad v =0 AND hfunction-kappa =0
            !       (c_mpb.gt.c_cut.and.(&
                !       ((dielec_func_mpb(i_full_points)-1d0).gt.eps_cut.and.&
                !       (epsinf_mpb-dielec_func_mpb(i_full_points)).gt.eps_cut).or.&
                !       (abs(kappa_debye**2-h_function(i_full_points)/dielec_func_mpb(i_full_points)).gt.eps_cut)&
                !       ))))) then
            !
            !       partition_tab_mpb = 1
            !     end if
        else !routine != 0
            !   if (myid==0)  write(use_unit,*) 'in ptab', epsinf_mpb,dielec_func_mpb(i_full_points),kappainf_mpb,kappa_func_mpb(i_full_points),&
                !       abs(kappa_current**2-kappa_debye**2),&
                !
            !       (sum(dielec_func_gradient_mpb(:,i_full_points)**2).gt.eps_cut.or.&
            !       abs(kappa_current**2-kappa_debye**2).gt.kappa_cut),&
                !
            !       (sum(dielec_func_gradient_mpb(:,i_full_points)**2).le.eps_cut.and.&
                !       abs(kappa_current**2-kappa_debye**2).le.kappa_cut),&
                !
            !       (sum(dielec_func_gradient_mpb(:,i_full_points)**2).le.eps_cut.and.&
                !       abs(kappa_current**2-kappa_debye**2).gt.kappa_cut.and.epsinf_mpb-dielec_func_mpb(i_full_points).le.eps_cut)

            !only_outer: recalculate the potential on the grid points which we have thrown away so far
            !non only_outer: calculate the potential on the points where the part of q_iter that depends on the potential
            !is zero anyway
            if (full_iteration==0.or.full_iteration==1) then
                if (full_iteration==0.or.&
                        !LPBE case
                        !  |calculate all non-zero ("inside")  points (during iterative cycle)
                        (full_iteration==1.and..not.only_outer.and.&
                        ((c_mpb.le.c_cut.and.&
                        (dielec_func_mpb(i_full_points)-1d0).gt.eps_cut.and.&
                        (epsinf_mpb-dielec_func_mpb(i_full_points)).gt.eps_cut).or.&
                        (c_mpb.gt.c_cut.and.&
                        (sum(dielec_func_gradient_mpb(:,i_full_points)**2).gt.eps_cut.or.&
                        abs(kappa_current**2-kappa_debye**2).gt.kappa_cut)))).or.&
                        !  |only calculate "outer" points (that have not been calculated yet)
                        (full_iteration==1.and.only_outer.and.&
                        ((c_mpb.le.c_cut.and.((dielec_func_mpb(i_full_points)-1d0).le.eps_cut.or.&
                        (epsinf_mpb-dielec_func_mpb(i_full_points)).le.eps_cut)).or.&
                        (c_mpb.gt.c_cut.and.&
                        (sum(dielec_func_gradient_mpb(:,i_full_points)**2).le.eps_cut.and.&
                        abs(kappa_current**2-kappa_debye**2).le.kappa_cut))))) then
                    partition_tab_mpb = 1
                end if
            else if (full_iteration==2) then
                !MPBE case
                !  |calculate all non-zero ("inside") points (during iterative cycle)
                if ((full_iteration==2.and..not.only_outer.and.&
                        ((c_mpb.le.c_cut.and.&
                        (dielec_func_mpb(i_full_points)-1d0).gt.eps_cut.and.&
                        (epsinf_mpb-dielec_func_mpb(i_full_points)).gt.eps_cut).or.&
                        (c_mpb.gt.c_cut.and.(((dielec_func_mpb(i_full_points)-1d0).gt.eps_cut.and.&
                        (epsinf_mpb-dielec_func_mpb(i_full_points)).gt.eps_cut).or.&
                        (abs(kappa_debye**2-h_function(i_full_points)/dielec_func_mpb(i_full_points)).gt.ion_cut))))).or.&
                        !  |only calculate "outer" points (that have not been calculated yet)
                        (full_iteration==2.and.only_outer.and.&
                        ((c_mpb.le.c_cut.and.((dielec_func_mpb(i_full_points)-1d0).le.eps_cut.or.&
                        (epsinf_mpb-dielec_func_mpb(i_full_points)).le.eps_cut)).or.&
                        (c_mpb.gt.c_cut.and.(((dielec_func_mpb(i_full_points)-1d0).le.eps_cut.or.&
                        (epsinf_mpb-dielec_func_mpb(i_full_points)).le.eps_cut).and.&
                        abs(kappa_debye**2-h_function(i_full_points)/dielec_func_mpb(i_full_points)).le.ion_cut))))&
                        ) then
                    partition_tab_mpb = 1
                end if
            end if
            !     if (full_iteration==0.or.&
            !       !LPBE case
            !       !  |calculate all non-zero ("inside")  points (during iterative cycle)
            !       (full_iteration==1.and..not.only_outer.and.&
            !       ((c_mpb.le.c_cut.and.&
            !       (dielec_func_mpb(i_full_points)-1d0).gt.eps_cut.and.&
            !       (epsinf_mpb-dielec_func_mpb(i_full_points)).gt.eps_cut).or.&
            !       (c_mpb.gt.c_cut.and.&
            !       (sum(dielec_func_gradient_mpb(:,i_full_points)**2).gt.eps_cut.or.&
            !       abs(kappa_current**2-kappa_debye**2).gt.kappa_cut)))).or.&
            !       !  |only calculate "outer" points (that have not been calculated yet)
            !       (full_iteration==1.and.only_outer.and.&
            !       ((c_mpb.le.c_cut.and.((dielec_func_mpb(i_full_points)-1d0).le.eps_cut.or.&
            !       (epsinf_mpb-dielec_func_mpb(i_full_points)).le.eps_cut)).or.&
            !       (c_mpb.gt.c_cut.and.&
            !       (sum(dielec_func_gradient_mpb(:,i_full_points)**2).le.eps_cut.and.&
            !       abs(kappa_current**2-kappa_debye**2).le.kappa_cut)))).or.&
            !       !MPBE case
            !       !  |calculate all non-zero ("inside") points (during iterative cycle)
            !       (full_iteration==2.and..not.only_outer.and.&
            !       ((c_mpb.le.c_cut.and.&
            !       (dielec_func_mpb(i_full_points)-1d0).gt.eps_cut.and.&
            !       (epsinf_mpb-dielec_func_mpb(i_full_points)).gt.eps_cut).or.&
            !       (c_mpb.gt.c_cut.and.(((dielec_func_mpb(i_full_points)-1d0).gt.eps_cut.and.&
            !       (epsinf_mpb-dielec_func_mpb(i_full_points)).gt.eps_cut).or.&
            !       (abs(kappa_debye**2-h_function(i_full_points)/dielec_func_mpb(i_full_points)).gt.ion_cut))))).or.&
            !       !  |only calculate "outer" points (that have not been calculated yet)
            !       (full_iteration==2.and.only_outer.and.&
            !       ((c_mpb.le.c_cut.and.((dielec_func_mpb(i_full_points)-1d0).le.eps_cut.or.&
            !       (epsinf_mpb-dielec_func_mpb(i_full_points)).le.eps_cut)).or.&
            !       (c_mpb.gt.c_cut.and.(((dielec_func_mpb(i_full_points)-1d0).le.eps_cut.or.&
            !       (epsinf_mpb-dielec_func_mpb(i_full_points)).le.eps_cut).and.&
            !       abs(kappa_debye**2-h_function(i_full_points)/dielec_func_mpb(i_full_points)).le.ion_cut))))&
            !       ) then
            ! 	partition_tab_mpb = 1
            !       end if
        end if !routine == 0

        return
    end function

    !------------------------------------------------------------------------------
    !****s* SPE_solver/sum_up_mpe_potential
    !  NAME
    !    sum_up_mpe_potential
    !  SYNOPSIS
    subroutine sum_up_mpe_potential(epsinfinv_factor,use_partition, kerker_gradient, &
            use_analytic_far_field, set_mult_rad_to_addtomrad, &
            addtomprad,multipole_moment_ff,kappa_SPE_step,&
            multipole_thresh, R_prec, R_prec_gradient, full_iteration,only_outer,evaluate_laplacian,&
            R_prec_laplacian,evaluate_only_mpe,lmax_cut,mpb_forces_on)
    !  PURPOSE
    ! do the actual multipole summation of eq. (40)
        use constants, only: pi
        use runtime_choices, only: extra_adding_to_hartree_potential_distance
        use grids
        use dimensions
        use synchronize_mpi
        use pbc_lists
        use physics, only: partition_tab
        use lpb_solver_utilities, only: changed_SPE_ff_cutoffs
        use spline, only: spline_vector, spline_deriv_vector, &
            spline_2nd_deriv_vector
        use localorb_io, only: localorb_info, use_unit
        use geometry, only: species

        implicit none

        !INPUTS
        logical, intent(in) :: epsinfinv_factor, use_partition, kerker_gradient, &
            use_analytic_far_field, set_mult_rad_to_addtomrad, &
            evaluate_laplacian, mpb_forces_on
        logical, intent(in) :: evaluate_only_mpe !< if true evaluate only mp expansion of input density not of integral
        logical :: only_outer
        integer :: full_iteration
        real*8, dimension(n_species), intent(in) :: addtomprad
        ! radial Green functions of modified Helmholtz Eqn, fct of r, l, and species
        real*8, intent(in) :: multipole_thresh, kappa_SPE_step

        !OUTPUTS
        real*8, dimension(n_full_points), intent(inout) :: R_prec
        
        real*8, dimension(3,n_full_points), intent(inout) :: R_prec_gradient
        real*8, dimension(n_full_points), intent(out) :: R_prec_laplacian

        integer :: i_center, i_l, i_m, current_spl_atom, current_center, i_full_points,&
            i_my_batch, i_index, i_coords, i_grid, prev_atom, i_spline, i_atom,&
            current_spl_atom_calc
        character*1000 :: info_str
        real*8, dimension(3) :: coord_current, dir_tab, dir_tab_in
        real*8 :: radius_rad, log_weight, radial_weight, kappa_times_rgrid, delta, max_spl_rad_sq,&
            dist_tab_sq,  trigonom_tab(4), ddot, radius
        real*8, dimension(limit_l_dim_SPE) :: ylm_tab

        !far field
        real*8, dimension(limit_l_dim_SPE_ff) :: ylm_tab_ff
        real*8,  dimension(limit_l_max_SPE_ff+1) :: green_K_current
        real*8, dimension(limit_l_max_SPE_ff+1) :: green_K_deriv
        real*8, dimension(n_spline_atoms,limit_l_dim_SPE_ff) :: multipole_moment_ff

        !SR: if the multipoles should be cut off already in the near field, we can use this
        !cutoff radius beyond which only multipoles up to lmax_cut will be considered in the summation
        integer :: lmax_cut, current_atom,i,j
        real*8 :: radius_cut

        real*8 :: kappa_current, current_rprec
        
        real*8 :: sqrt_pi_kappa_times_rgrid_radius_sq, sqrt_pi_kappa_times_rgrid
        real*8 :: sqrt_dist_tab_sq  !< compute the sqrt of dist_tab_sq only once and use it everywhere else
        
        real*8, dimension(n_full_points) :: R_temp
        real*8, dimension(3,n_full_points) :: R_temp_gradient
        
        ! ========== BEGIN INITIALIZATIONS ==========
        if (debug_lpb_solver) then
            if (myid==0) write(use_unit,*) 'limit_l_max_SPE_ff', limit_l_max_SPE_ff
        end if

        ! ---------- begin force initializations ----------
        if (mpb_forces_on) then
            if (atomic_MERM) then
                do i_atom  = 1, n_spline_atoms, 1
                    if (.not. atom_converged(i_atom)) then
                        delta_v_MPB_multipole_component_at_zero(:,i_atom,:) = 0d0
                        delta_v_MPB_multipole_deriv_at_zero(:,i_atom,:) = 0d0
                    end if
                end do
            else
                delta_v_MPB_multipole_component_at_zero = 0d0
                delta_v_MPB_multipole_deriv_at_zero = 0d0
            end if
        end if

        if (use_forces.and..not.forces_mpb_force_off.and..not.only_outer) then
            !SR: until now we are calculating delta_v_MPB_gradient_atom in every SCF step. we however need this in principle
            ! only in the step before the force calculations, for the pulay force calculations which are
            !before the call of sum_up_whole_potential we need these derivatives
            if (atomic_MERM) then
                do i_atom = 1, n_spline_atoms, 1
                    if (.not. atom_converged(i_atom)) then
                        delta_v_MPB_gradient_atom(i_atom,:,:) = 0d0
                    end if
                end do
            else
                delta_v_MPB_gradient_atom = 0d0
            end if
        end if
        ! ---------- end force initializations ----------

        radius_cut = 5d0

        ! FINALLY: sum up the solved and splined preconditioned residual from all different atoms
        !   over all grid points.
        !  need to pay attention to the memory location of the various splines....
        ! intialize variables: R_prec has not been used so far
        R_prec_at_zero = 0d0
        current_spl_atom = 0
        current_spl_atom_calc = 0

        R_prec_laplacian = 0d0

        if (full_iteration==0) then
            green_K_current = 0d0
            green_K_deriv = 0d0
            aux_R_prec_result_mpb_ff = 0d0
        else
            if (atomic_MERM) then
                do i_atom = 1, n_spline_atoms, 1
                    if (.not. atom_converged(i_atom)) then
                        multipole_moment_ff(i_atom,:) = multipole_moment_ff_glob(i_atom,:)
                    end if
                end do
            else
                multipole_moment_ff = multipole_moment_ff_glob
            end if
        end if

        if (atomic_MERM) then
            if (.not.only_outer) then ! .or. (only_outer .and. atomic_MERM)) then
                do i_atom=1, n_spline_atoms, 1
                    if (.not. atom_converged(i_atom)) then
                        do j = 1, index_atoms_size(i_atom), 1 
                            i = index_atoms(i_atom,j)
                            R_prec(i) = 0d0
                            R_prec_gradient(:,i) = 0d0
                        end do
                    end if
                end do
            else
                !save the data we have already calculated
                R_temp = R_prec
                R_temp_gradient = R_prec_gradient
                !initialize the resulting arrays
                R_prec = 0d0
                R_prec_gradient = 0d0
            end if
        else if (.not. only_outer) then
            R_prec = 0d0
            R_prec_gradient = 0d0
        end if

        if (mpb_forces_on) then
            R_prec_gradient_at_zero = 0d0
        end if

        ! ========== END OF INITIALIZATIONS ==========
        
        !short explanation on the loops here for atomic_MERM = .True.:
        !if only_outer = .False. 
        !    - we are in the middle of the MERM iterations
        !    - we evaluate the potential on all atomic grids of the unit cell atoms
        !    - each grid point gets only contributions from the atom center
        !if only_outer = .True.
        !    - this is the last MERM step
        !    - we need to update the whole potential, meaning that we need to evaluate all contributions from all centers on the FHI-aims grid points
        !    - we already calculated the on-side potential before, so we don't need to do this again, therefore only evaluate contributions where current_center.ne.current_spl_atom

        do i_center = 1, n_centers_hartree_potential, 1 

            prev_atom        = current_spl_atom_calc
            current_center   = centers_hartree_potential(i_center)
            current_spl_atom = center_to_atom(current_center)       ! this is the actual atom number, .ne. current_center for PBC's
            ! obtain the atomic spline for current center from where it is actually stored ...
            if (.not.only_outer .and. atomic_MERM .and. atom_converged(current_spl_atom)) then
                cycle
            end if

            !in the atomic_MERM case we want to some over images only in the last iteration step
            !this basically makes the loop over centers a loop over the unit cell atoms only!!
            if (atomic_MERM.and..not.only_outer.and.current_center.ne.current_spl_atom) then
               cycle
            end if

            !in this atomic_MERM case, current_spl_atom_calc is the spline atom which is really
            !calculated, so that prev_atom references to the correct atom
            current_spl_atom_calc = current_spl_atom
            
            if (current_spl_atom.ne.prev_atom) then
                ! distribute the information about that particular atom to all the threads
                ! simultaneously, each one is going to use it in order to compute its part
                ! of the new residual ...

                current_R_multipole_spl_mpb = 0d0

                call get_R_prec_spline_mpi(current_spl_atom,limit_l_dim_SPE,current_R_multipole_spl_mpb,R_multipole_spl_mpb, n_spline_grid_dim)

                ! ========== CUTOFF PARAMS ==========
                !
                ! Determine the cutoff distance for the far field case. If this was not set in the input file as
                ! SPE_cut_and_lmax_ff (which then leads to changed_SPE_ff_cutoffs being .true.), determine a suitable
                ! value from current_R_multipole_spl_mpb by starting far away from the atom and then moving closer to
                ! the nucleus while checking whether we encounter a value greater than the multipole_thresh.
                !
                ! This distance should be a function of l and m, therefore the loops below. The largest distance that
                ! is encountered in this procedure is then stored in max_spl_rad_sq for further use.
                !
                ! max_lm_spl_rad_sq_mpb = max radius for i_center, where multipole moment has significant value
                ! max_spl_rad_sq        = max of all max_lm_spl_rad_sq_mpb
                ! multipole_radius_free = innermost point where the free-atom charge density is still zero
                ! extra_adding_to_hartree_potential_distance
                if (.not. changed_SPE_ff_cutoffs) then
                    max_spl_rad_sq = 0d0

                    do i_l = 0, limit_l_max_SPE_species(species(current_spl_atom))
                        do i_m = -i_l, i_l

                            ! start with the furthest point on each grid and move inwards until you find non-zero splines.
                            if (evaluate_only_mpe) then
                                i_grid = n_radial(species(current_spl_atom)) + 1
                            else
                                i_grid = n_grid(species(current_spl_atom)) + 1
                            end if

                            ! check through splines ... see if they are equal to zero, decrease i_grid if they are
                            ! the check is representative when we square the multipole coefficients.
                            ! FIXME: one might also consider a threshold for delta below which 'there is no far field'
                            delta = 0d0
                            do while ((delta .le. multipole_thresh) .and. (i_grid .gt. 1))   ! FIXME: this should of course be a threshold 1E-12
                                i_grid = i_grid - 1
                                do i_spline = 1, n_max_spline_mpb
                                    ! FIXME: this would be good enough if we only used the i_spline = 1 (or 0?) value !!!
                                    delta = delta + current_R_multipole_spl_mpb(index_lm_mpb(i_m, i_l), i_spline, i_grid)**2
                                end do
                            end do

                            ! set maximal radius for the atomic preconditioned residual multipole decomposition
                            if (evaluate_only_mpe) then
                                max_lm_spl_rad_sq_mpb(index_lm_mpb(i_m, i_l)) = r_radial(i_grid, species(current_spl_atom))**2
                            else
                                max_lm_spl_rad_sq_mpb(index_lm_mpb(i_m, i_l)) = r_grid(i_grid, species(current_spl_atom))**2
                            end if

                            if (max_lm_spl_rad_sq_mpb(index_lm_mpb(i_m, i_l)) .gt. max_spl_rad_sq) then
                                max_spl_rad_sq = max_lm_spl_rad_sq_mpb(index_lm_mpb(i_m, i_l))
                            end if
                        end do
                    end do
                end if

                ! ========== END CUTOFF PARAMS ==========

                if (debug_lpb_solver) then
                    write(use_unit,*) 'max_spl_rad_sq',max_spl_rad_sq, (multipole_radius_free_SPE(species(current_spl_atom)) &
                        +extra_adding_to_hartree_potential_distance)**2
                end if

                ! output warning, atom number etc if the maximal spline radius is exactly equal to the outermost
                ! grid radius
                if ((max_spl_rad_sq.gt.(multipole_radius_free_SPE(species(current_spl_atom)) &
                        +extra_adding_to_hartree_potential_distance)**2).and. &
                        SPE_first_warnings) then

                    ! FIXME: VB suggests that you should stop right here, but I don't want to do this until the routine is known to be working.

                    write(info_str,'(2X,A)'   ) '* WARNING: Multipole expansion of preconditioned residual '
                    call localorb_info(info_str,use_unit,'(A)')
                    write(info_str,'(2X,A,I4)') '* extends too far at atom ', current_spl_atom
                    call localorb_info(info_str,use_unit,'(A)')
                end if

                if (set_mult_rad_to_addtomrad) then
                    max_spl_rad_sq = 0.0d0
                end if
            end if !current_spl_atom.ne.prev_atom

            ! FIXME: This is done without the long-range contribution in periodic cases, does that matter?
            ! what are the physical implementations of the 'Far-field' ? Does the filtering of low-frequencies
            ! enhance the regions which a given atom influences? - see comments in header of this file
            ! loop over all BATCHES

            i_full_points = 0

            do i_my_batch = 1, n_my_batches, 1

                do i_index = 1, batches(i_my_batch)%size, 1
                    i_full_points   = i_full_points + 1

                    !calculate atom to which this grid point belongs to
                    current_atom    = batches(i_my_batch) % points(i_index) % index_atom
                    !in case we perform atom-wise iterations, we need to update the potential only on the grid of the
                    !respective atom, since the following multipole expansion evaluate_density_mpe will also perform
                    !atom-centered angular integrations

!                     if (atomic_MERM .and. (current_atom .ne. i_atom_MERM) .and. .not.&
!                         only_outer) then
                    if (atomic_MERM .and. (current_atom .ne. current_spl_atom) .and. .not. only_outer) then
                        !for atomic_MERM = .True.:
                        !in case of any but the last iteration step, we need to simply evaluate the potential on a single 
                        !   atomic grid (NOT partitioned!)
                        !in case of the last iteration step (only_outer = .True.), we need to evaluate the atomic potential on the full
                        !   FHI-aims integration grid of all atoms (partitioned!)
                        cycle
                    end if

                    !SR: actually in case of the atomic_MERM method, we would need to update the
                    !potential on the full grid not only the partitioned one. The reason for that is
                    !that partition_tab_mpb can slightly change during the SCF cycle, since n_el
                    !changes and therefore we need the correct R_prec_atom_wise to initialize the
                    !MERM of the next SCF cycle. We leave this out here anyways, because we believe
                    !the benefit of updating the potential on the full grid is not worth the small
                    !benefit in better initialization
                    if (partition_tab(i_full_points).gt.0.d0 .or. (atomic_MERM .and. .not.only_outer)) then
                    
                        kappa_current = sqrt(alpha_func_mpb(i_full_points)/(1d0+phi_zero_mpb*(alpha_func_mpb(i_full_points)-1d0)))*&
                        kappainf_mpb/ sqrt(dielec_func_mpb(i_full_points)) * sqrt_z_mpb_div_kBT_mpb

                        !if we currently consider an atom of the unit cell and its contribution on the atomic grid belonging to this atom
                        !we can just take the potential of the last step, since we already evaluated this 
                        if (atomic_MERM .and. (current_center.eq.current_spl_atom).and.(current_atom .eq. current_spl_atom) .and.&
                            only_outer .and..not.partition_tab_mpb(1,i_full_points,0,full_iteration,only_outer,kappa_current)==1) then
                            !we have already calculated the contribution from current_center to R_prec(i_full_points) in the last iteration of the respective atom
                            R_prec(i_full_points) = R_prec(i_full_points) + R_temp(i_full_points)
                            R_prec_gradient(:,i_full_points) = R_prec_gradient(:,i_full_points) + R_temp_gradient(:,i_full_points)
                            cycle
                        end if

                        !the reverse to the case before in which we actually have to evaluate the potential
                        if (partition_tab_mpb(1,i_full_points,0,full_iteration,only_outer,kappa_current)==1 .or.&
                            (atomic_MERM .and. only_outer .and. ((current_center.ne.current_spl_atom).or.(current_atom .ne. current_spl_atom)))) then ! .and.&
                            !(.not.(atomic_MERM.and.only_outer))).or.&
                            !(atomic_MERM .and. only_outer .and.&
                            !!if we are using the atom-wise MERM, we need to calculate the potential
                            !!1.) at all points belonging to other atoms than i_atom_MERM
                            !!2.) at all points belonging to i_atom_MERM but the ones already
                            !!calculated
                            !(current_atom .ne. current_spl_atom .or.&
                            !(current_atom .eq. current_spl_atom .and. partition_tab_mpb(1,i_full_points,0,full_iteration,only_outer,kappa_current)==1)))) then
                            ! 		    write(use_unit,*) 'flags',full_iteration,only_outer
                            
                            coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
                            

                            call  tab_single_atom_centered_coords_p0( current_center, coord_current, &
                                dist_tab_sq, dir_tab )

                            ! ========== SUMMATION EQUATION ==========
                            ! the sum to be done here is the same as Delley (1990) Eqn (12c), give or take integration factors:
                            ! R_prec(r) = sum_{alpha,l,m} Y_lm (r-R_alpha) R_prec_{alpha,l,m}(|r-R_alpha|)

                            ! from the current_center and coord_current, compute the distance squared, dist_tab_sq,
                            ! and the vector from the center to the grid point, dir_tab
                            call tab_single_atom_centered_coords_p0( current_center, coord_current, &
                                dist_tab_sq, dir_tab )

                            radius = SQRT(dist_tab_sq)  ! this is the actual distance from point to atom: spline this ...
                            sqrt_dist_tab_sq = radius   ! store this for later use: avoid use of SQRT

                            ! trigonom_tab contains trig functions of polar angles in dir_tab:
                            call tab_single_trigonom_p0( dir_tab, trigonom_tab )

                            ! these can be used in the ylm tab calculation routine to get all necessary spherical harmonics:
                            ! INVESTIGATE/FIXME: this might be faster if we use the proper (and already known!) (l,m) dependence
                            call tab_single_wave_ylm_p0(current_center,trigonom_tab, &
                                limit_l_max_SPE_species,limit_l_max_SPE, ylm_tab )

                            if (radius.gt.radius_cut.and.lmax_cut.gt.0) then
                                limit_l_dim_SPE_species_current = (lmax_cut + 1)**2
                            else
                                limit_l_dim_SPE_species_current = limit_l_dim_SPE_species
                            end if

                            ! in the following if statement, obtain spline-interpolated values of the multipole
                            ! components of the preconditioned residual, splined on the logarithmic integration grid

                            ! is current grid point within the realm of current_spl_atom?
                            if (dist_tab_sq .le. max_spl_rad_sq+addtomprad(species(current_spl_atom))) then

                                ! DEBUG: keep the actual radius for a moment, need it for output in a few lines
                                !                 delta = radius

                                ! need to invert on the logarithmic grid to find the index of the point in question and
                                ! where we should be looking between the two grid points (for all l,m values at once)

                                ! INVESTIGATE/FIXME: this might be faster if we use the proper (and already known!)
                                ! (l,m) dependence !

                                if (evaluate_only_mpe) then
                                    radius_rad = invert_radial_grid(radius,         &
                                        n_radial(species_center(current_center)),   &
                                        scale_radial(species_center(current_center)))
                                    call spline_vector                                          &
                                        ( radius_rad,                                           &
                                        current_R_multipole_spl_mpb,                            &
                                        n_spline_grid_dim,                                      &
                                        limit_l_dim_SPE,                                        &   ! FIXME: this is only necessary for some maximum l, should be faster
                                        n_radial(species_center(current_center)),               &
                                        limit_l_dim_SPE_species(species_center(current_center)),&
                                        aux_R_prec_result_mpb)
                                    radius = invert_log_grid(radius,            &
                                        r_grid_min(species(current_spl_atom)),  &
                                        r_grid_inc(species(current_spl_atom)))
                                else
                                    radius = invert_log_grid(radius,            &
                                        r_grid_min(species(current_spl_atom)),  &
                                        r_grid_inc(species(current_spl_atom)))
                                    call spline_vector                                          &
                                        ( radius,                                               &
                                        current_R_multipole_spl_mpb,                            &
                                        n_spline_grid_dim,                                      &
                                        limit_l_dim_SPE,                                        &   ! FIXME: this is only necessary for some maximum l, should be faster
                                        n_grid(species_center(current_center)),                 &
                                        limit_l_dim_SPE_species(species_center(current_center)),&
                                        aux_R_prec_result_mpb)
                                end if

                                ! multipole_moments(1:limit_l_dim_SPE_species(species_center(current_center)),i_full_points) = &
                                !     aux_R_prec_result_mpb(1:limit_l_dim_SPE_species(species_center(current_center)))

                                ! add electronic contribution from this particular source to R_prec
                                if (full_iteration==1.or.full_iteration==2) then
                                    !current multipole moment
                                    current_rprec = ddot ( limit_l_dim_SPE_species_current(species_center(current_center)), &
                                        aux_R_prec_result_mpb(:), 1, ylm_tab(:), 1 )
                                    R_prec(i_full_points) = R_prec(i_full_points) + current_rprec
                                        
                                else
                                    current_rprec = ddot ( limit_l_dim_SPE_species_current(species_center(current_center)), &
                                        aux_R_prec_result_mpb(:), 1, ylm_tab(:), 1 )
                                    R_prec(i_full_points) = R_prec(i_full_points)+ current_rprec

                                end if !full_iteration


                                if (kerker_gradient) then
                                    !explicitly calculate the gradient

                                    radius_rad = invert_radial_grid(sqrt_dist_tab_sq,   &
                                        n_radial(species_center(current_center)),       &
                                        scale_radial(species_center(current_center)))

                                    call tab_single_radial_weights_v2                       &
                                        ( current_spl_atom, sqrt_dist_tab_sq, radius_rad,   &
                                        log_weight, radial_weight )

                                    aux_R_prec_deriv_result = 0.0d0

                                    if (evaluate_only_mpe) then
                                        call spline_deriv_vector                                    &
                                            ( radius_rad,                                           &
                                            current_R_multipole_spl_mpb,                            &
                                            n_spline_grid_dim,                                      &
                                            limit_l_dim_SPE,                                        &
                                            n_radial(species_center(current_center)),               &
                                            limit_l_dim_SPE_species(species_center(current_center)),&
                                            aux_R_prec_deriv_result)
                                        aux_R_prec_deriv_result(1:limit_l_dim_SPE_species(species_center(current_center))) =&
                                            aux_R_prec_deriv_result(1:limit_l_dim_SPE_species(species_center(current_center))) &
                                            * radial_weight
                                    else
                                        call spline_deriv_vector                                    &
                                            ( radius,                                               &
                                            current_R_multipole_spl_mpb,                            &
                                            n_spline_grid_dim,                                      &
                                            limit_l_dim_SPE,                                        &
                                            n_grid(species_center(current_center)),                 &
                                            limit_l_dim_SPE_species(species_center(current_center)),&
                                            aux_R_prec_deriv_result)
                                        aux_R_prec_deriv_result(1:limit_l_dim_SPE_species(species_center(current_center))) =&
                                            aux_R_prec_deriv_result(1:limit_l_dim_SPE_species(species_center(current_center))) &
                                            * log_weight
                                    end if

                                    if (evaluate_laplacian) then
                                        !evaluate second derivative for laplacian
                                        aux_R_prec_2nd_deriv_result = 0d0

                                        if (evaluate_only_mpe) then
                                            call spline_2nd_deriv_vector                                &
                                                ( radius_rad,                                           &
                                                current_R_multipole_spl_mpb,                            &
                                                n_spline_grid_dim,                                      &
                                                limit_l_dim_SPE,                                        &
                                                n_radial(species_center(current_center)),               &
                                                limit_l_dim_SPE_species(species_center(current_center)),&
                                                aux_R_prec_2nd_deriv_result)
                                            aux_R_prec_2nd_deriv_result(1:limit_l_dim_SPE_species(species_center(current_center))) =&
                                                aux_R_prec_2nd_deriv_result(1:limit_l_dim_SPE_species(species_center(current_center))) &
                                                * radial_weight**2&
                                                -aux_R_prec_deriv_result(1:limit_l_dim_SPE_species(species_center(current_center)))*&
                                                1d0/sqrt_dist_tab_sq
                                        else
                                            call spline_2nd_deriv_vector                                &
                                                ( radius,                                               &
                                                current_R_multipole_spl_mpb,                            &
                                                n_spline_grid_dim,                                      &
                                                limit_l_dim_SPE,                                        &
                                                n_grid(species_center(current_center)),                 &
                                                limit_l_dim_SPE_species(species_center(current_center)),&
                                                aux_R_prec_2nd_deriv_result)
                                            aux_R_prec_2nd_deriv_result(1:limit_l_dim_SPE_species(species_center(current_center))) =&
                                                aux_R_prec_2nd_deriv_result(1:limit_l_dim_SPE_species(species_center(current_center))) &
                                                * log_weight**2&
                                                -aux_R_prec_deriv_result(1:limit_l_dim_SPE_species(species_center(current_center)))*&
                                                1d0/sqrt_dist_tab_sq
                                        end if
                                    end if

                                    dir_tab_in = dir_tab / sqrt_dist_tab_sq

                                    call tab_single_gradient_ylm_p0( trigonom_tab,  limit_l_max_SPE_species, &
                                        limit_l_max_SPE, current_center, ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab)
                                    aux_R_prec_gradient_result = 0.0d0
                                    call evaluate_g_gradient(limit_l_dim_SPE_species(species_center(current_center)), &
                                        limit_l_dim_SPE_species_current(species_center(current_center)),&
                                        1.0d0/sqrt_dist_tab_sq, dir_tab_in, trigonom_tab, limit_l_max_SPE, ylm_tab, dylm_dtheta_tab, &
                                        scaled_dylm_dphi_tab, aux_R_prec_result_mpb, aux_R_prec_deriv_result, aux_R_prec_2nd_deriv_result, &
                                        aux_R_prec_gradient_result,limit_l_max_SPE_species(species(current_spl_atom)),&
                                        aux_R_prec_laplacian_result,index_lm_mpb,evaluate_laplacian)

                                    do i_coords = 1, 3, 1
                                        R_prec_gradient(i_coords,i_full_points) = R_prec_gradient(i_coords,i_full_points) + &
                                            aux_R_prec_gradient_result(i_coords)
                                        !SR: at this point we would need to add a contribution in
                                        !case of (atomic_MERM .and. (.not. only_outer .or. (only_outer .and. current_spl_atom .eq. current_atom))), 
                                        !in order to get the correct R_prec_old in the next SCF step on the full atomic grid (not only on the one 
                                        !where partition_tab_mpb .ne. 0), since the region where partition_tab_mpb = 0 can change among SCF cycles. 
                                        !Seems to work however also without this modification
                                    end do

                                    if (use_forces.and..not.forces_mpb_force_off) then
                                        delta_v_MPB_gradient_atom(current_center,1:3,i_full_points) =&
                                            delta_v_MPB_gradient_atom(current_center,1:3,i_full_points) + &
                                            aux_R_prec_gradient_result(:)
                                    end if

                                    if (full_iteration==1.or.full_iteration==2) then

                                        if (evaluate_laplacian) then
                                            R_prec_laplacian(i_full_points) = R_prec_laplacian(i_full_points) + &
                                                aux_R_prec_laplacian_result
                                        end if

                                    else !if .not.full_iteration

                                        if (evaluate_laplacian) then
                                            R_prec_laplacian(i_full_points) = R_prec_laplacian(i_full_points) + &
                                                aux_R_prec_laplacian_result
                                        end if
                                    end if !full_iteration

                                end if !kerker_gradienit

                            else if (use_analytic_far_field) then ! far field calculation
                                !calculate Hankel function at current radial shell
                                ! calculate all l-values at once - for the radial argument (q_0 r)
                                call tab_single_wave_ylm_p0(current_center,trigonom_tab, &
                                    limit_l_max_SPE_species_ff,limit_l_max_SPE_ff, ylm_tab_ff )

                                if (use_partition) then
                                    kappa_times_rgrid = kappa_SPE_step*radius
                                    if (kerker_gradient) then
                                        call bessel_K_deriv(kappa_times_rgrid,kappa_SPE_step,&
                                            limit_l_max_SPE_ff, green_K_deriv)
                                    end if
                                else
                                    kappa_times_rgrid = kappa_SPE*radius
                                    if (kerker_gradient) then
                                        call bessel_K_deriv(kappa_times_rgrid,kappa_SPE,&
                                            limit_l_max_SPE_ff, green_K_deriv)
                                    end if
                                end if
                                call bessel_K_mpb(kappa_times_rgrid,limit_l_max_SPE_ff, green_K_current)
                                sqrt_pi_kappa_times_rgrid = sqrt(pi/(2.*kappa_times_rgrid))
                                sqrt_pi_kappa_times_rgrid_radius_sq = sqrt(pi/(8.*kappa_times_rgrid*radius**2))

                                if (kerker_gradient) then
                                    green_K_deriv = green_K_deriv*sqrt_pi_kappa_times_rgrid-&
                                        green_K_current*sqrt_pi_kappa_times_rgrid_radius_sq
                                end if

                                green_K_current = green_K_current*sqrt_pi_kappa_times_rgrid

                                do i_l = 0,limit_l_max_SPE_species_ff(species(current_spl_atom)),1
                                    do i_m = -i_l, i_l, 1
                                        aux_R_prec_result_mpb_ff(index_lm_mpb(i_m,i_l)) = green_K_current(i_l+1) * &
                                            multipole_moment_ff(current_spl_atom,index_lm_mpb(i_m,i_l))
                                    end do
                                end do
                                current_rprec = ddot ( limit_l_dim_SPE_species_ff(species_center(current_center)), &
                                    aux_R_prec_result_mpb_ff(:), 1, ylm_tab_ff(1:limit_l_dim_SPE_ff), 1 )
                                R_prec(i_full_points) = R_prec(i_full_points) + current_rprec
                                    

                                if (kerker_gradient) then

                                    radius_rad = invert_radial_grid(radius,                                &
                                        n_radial(species_center(current_center)), &
                                        scale_radial(species_center(current_center)))

                                    do i_l = 0, limit_l_max_SPE_species_ff(species(current_spl_atom)), 1
                                        do i_m = -i_l, i_l, 1
                                            aux_R_prec_deriv_result_ff(index_lm_mpb(i_m,i_l)) = green_K_deriv(i_l+1) * &
                                                multipole_moment_ff(current_spl_atom,index_lm_mpb(i_m,i_l))
                                        end do
                                    end do

                                    dir_tab_in = dir_tab / radius

                                    call tab_single_gradient_ylm_p0( trigonom_tab,  limit_l_max_SPE_species_ff, &
                                        limit_l_max_SPE_ff, current_center, &
                                        ylm_tab_ff, dylm_dtheta_tab_ff, scaled_dylm_dphi_tab_ff)

                                    call evaluate_g_gradient(limit_l_dim_SPE_species_ff(species_center(current_center)), &
                                        limit_l_dim_SPE_species_ff(species_center(current_center)), &
                                        1.0d0/radius, dir_tab_in, trigonom_tab, &
                                        limit_l_max_SPE_ff, ylm_tab_ff(1:limit_l_dim_SPE_ff),&
                                        dylm_dtheta_tab_ff(1:limit_l_dim_SPE_ff), &
                                        scaled_dylm_dphi_tab_ff(1:limit_l_dim_SPE_ff),&
                                        aux_R_prec_result_mpb_ff(1:limit_l_dim_SPE_ff),&
                                        aux_R_prec_deriv_result_ff(1:limit_l_dim_SPE_ff), &
                                        aux_R_prec_2nd_deriv_result_ff(1:limit_l_dim_SPE_ff), &
                                        aux_R_prec_gradient_result,&
                                        limit_l_max_SPE_species_ff(species(current_spl_atom)),&
                                        aux_R_prec_laplacian_result,index_lm_mpb,evaluate_laplacian)

                                    do i_coords = 1, 3, 1
                                        R_prec_gradient(i_coords,i_full_points) = R_prec_gradient(i_coords,i_full_points) + &
                                        aux_R_prec_gradient_result(i_coords)
!                                         if (atomic_MERM .and. (.not. only_outer .or. (only_outer .and. current_spl_atom .eq. current_atom))) then
                                    end do

                                    if (use_forces.and..not.forces_mpb_force_off) then
                                        delta_v_MPB_gradient_atom(current_center,1:3,i_full_points) =&
                                            delta_v_MPB_gradient_atom(current_center,1:3,i_full_points) +&
                                            aux_R_prec_gradient_result(:)
                                    end if

                                end if  ! kerker gradient

                            end if  ! use_ANALYTIC_farfield, in mp radius

                        end if  ! if not full iteration or non-zero integrand
                    end if  ! partition_tab.ne.0

                end do  ! loop over points in a given batch

            end do  ! loop over batches

            !Now calculate potential at center of atoms for nuclei nuclei interaction correction

            ! if (epsinfinv_factor) then
            !     write(use_unit,*) 'check kerker_mpb. there should not be a epsinfinv_factor while calculating corrections'
            !     stop
            ! end if

            ! check whether this is the last iteration step if MERM is performed inside solve_SPE
            if (full_iteration==0.or.((full_iteration==1.or.full_iteration==2).and.only_outer)) then

                ! evaluate electrostatic potential and gradient at position of nuclei
                call evaluate_potential_at_zero(use_partition, kappa_SPE_step,max_spl_rad_sq,&
                    current_center, addtomprad,i_center, current_spl_atom,&
                    multipole_moment_ff, use_analytic_far_field,mpb_forces_on,radius_cut,lmax_cut)

            end if

        end do          ! loop over relevant Hartree centers

    end subroutine sum_up_mpe_potential

    !------------------------------------------------------------------------------
    !****s* SPE_solver/evaluate_potential_at_zero
    !  NAME
    !    evaluate_potential_at_zero
    !  SYNOPSIS
    subroutine evaluate_potential_at_zero(use_partition, kappa_SPE_step,max_spl_rad_sq,&
            current_center, addtomprad,i_center, current_spl_atom, &
            multipole_moment_ff,use_analytic_far_field,mpb_forces_on,radius_cut,lmax_cut)
    !  PURPOSE
    ! evaluate delta v and gradient of delta v at the center of the atom "current_center"
        use constants, only: pi
        use grids
        use dimensions
        use synchronize_mpi
        use pbc_lists
        use lpb_solver_utilities, only: scf_mul_vec_1
        use spline, only: spline_vector, spline_deriv_vector
        use geometry, only: coords, species

        implicit none

        !IN
        integer :: i_coords
        logical, intent(in) :: use_partition, use_analytic_far_field, mpb_forces_on
        real*8 :: max_spl_rad_sq, kappa_SPE_step
        integer, intent(in) :: current_center, i_center, current_spl_atom
        real*8, dimension(n_species), intent(in) :: addtomprad

        integer :: i_center_2, i_lm, current_spl_atom_2, current_center_2, i_l, i_m
        real*8 :: radius, kappa_times_rgrid, dist_tab_sq, log_weight, radial_weight, radius_rad
        real*8 :: coord_current(3), dir_tab(3), dir_tab_in(3)
        real*8 :: ddot, trigonom_tab(4)
        real*8, dimension(limit_l_dim_SPE) :: ylm_tab

        !far field
        real*8, dimension(limit_l_dim_SPE_ff) :: ylm_tab_ff
        real*8,  dimension(limit_l_max_SPE_ff+1) :: green_K_current
        real*8, dimension(limit_l_max_SPE_ff+1) :: green_K_deriv
        real*8, dimension(n_spline_atoms,limit_l_dim_SPE_ff),intent(in) :: multipole_moment_ff

        real*8 :: radius_cut
        integer :: lmax_cut
        real*8 :: sqrt_dist_tab_sq

        do i_center_2 = 1, n_centers_hartree_potential, 1

            !TODO: check, is the position here right? changed position from lower line to here

            current_center_2   = centers_hartree_potential(i_center_2)
            current_spl_atom_2 = center_to_atom(current_center_2)                     ! this is the actual atom number, .ne. current_center for PBC's

            coord_current(:) = coords(:,current_spl_atom_2)

            ! 	write(use_unit,*) 'check', current_spl_atom, coord_current

            call  tab_single_atom_centered_coords_p0( current_center, coord_current, &
                dist_tab_sq, dir_tab )



            ! trigonom_tab contains trig functions of polar angles in dir_tab:
            call tab_single_trigonom_p0( dir_tab, trigonom_tab )

            ! these can be used in the ylm tab calculation routine to get all necessary spherical harmonics:
            ! INVESTIGATE/FIXME: this might be faster if we use the proper (and already known!) (l,m) dependence
            call tab_single_wave_ylm_p0(current_center,trigonom_tab, &
                limit_l_max_SPE_species,limit_l_max_SPE, ylm_tab )
            ! obtain spline-interpolated values of the multipole components
            ! of the preconditioned residual, splined on the logarithmic
            ! integration grid
            radius = SQRT(dist_tab_sq)    ! this is the actual distance from point to atom: spline this ...


            if (dist_tab_sq.le.max_spl_rad_sq+addtomprad(species(current_spl_atom))) then            ! is current grid point within the realm of current_spl_atom?

                if (radius.gt.radius_cut.and.lmax_cut.gt.0) then
                    limit_l_dim_SPE_species_current = (lmax_cut + 1)**2
                else
                    limit_l_dim_SPE_species_current = limit_l_dim_SPE_species
                end if

                if (i_center_2.eq.i_center) then

                    !dist_tab_sq = corr_cutoff**2!
                    dist_tab_sq = r_grid_min(species(current_spl_atom))**2 + 1e-15
                    radius = SQRT(dist_tab_sq) !like for the hartree potential in sum_up_whole_potential

                end if

                dir_tab_in = dir_tab / radius
                sqrt_dist_tab_sq = radius
                radius = invert_log_grid(radius,                                &
                    r_grid_min(species(current_spl_atom)), &
                    r_grid_inc(species(current_spl_atom)))


                !spline value of multipole moment of atom1 at atom2
                call spline_vector                             &
                    ( radius,                                 &
                    current_R_multipole_spl_mpb,                &
                    n_spline_grid_dim,                      &
                    limit_l_dim_SPE,                             &   ! FIXME: this is only necessary for some maximum l, should be faster
                    n_grid(species_center(current_center)), &
                    limit_l_dim_SPE_species(species_center(current_center)), &
                    aux_R_prec_result_at_zero)



                if (i_center_2.eq.i_center) then !if nuclei interaction with its own potential
                    !use extrapolation for all moments from dipole on and use explicitly integrated value R_multipole_mpb_at_zero for
                    !the monopole moment
                    aux_R_prec_result_at_zero(1) = R_multipole_mpb_at_zero(current_spl_atom)

                    do i_lm = 2, limit_l_dim_SPE_species(species_center(current_center)), 1

                        aux_R_prec_result_at_zero(i_lm) = 2* current_R_multipole_spl_mpb(i_lm,1,1) &
                            - current_R_multipole_spl_mpb(i_lm,1,2)
                    end do

                end if

                if (mpb_forces_on) then
                    delta_v_MPB_multipole_component_at_zero(:,current_spl_atom,current_spl_atom_2) = aux_R_prec_result_at_zero(:)
                endif

                R_prec_at_zero(current_spl_atom_2) = R_prec_at_zero(current_spl_atom_2) + &
                    ddot ( limit_l_dim_SPE_species(species_center(current_center)), &
                    aux_R_prec_result_at_zero(:), 1, ylm_tab(:), 1 )
                !
                ! 		  end if !(use_partition)
                if ( mpb_forces_on ) then
                    radius_rad = invert_radial_grid(sqrt_dist_tab_sq,                                &
                        n_radial(species_center(current_center)), &
                        scale_radial(species_center(current_center)))

                    call tab_single_radial_weights_v2 &
                        ( current_spl_atom, sqrt_dist_tab_sq, radius_rad, &
                        log_weight, radial_weight )

                    aux_R_prec_deriv_result_at_zero = 0d0

                    !                    call tab_single_radial_weights_v2( current_spl_atom,dist_tab_sq,radius,log_weight,radial_weight )
                    call spline_deriv_vector(radius,current_R_multipole_spl_mpb,n_spline_grid_dim,limit_l_dim_SPE,&
                        n_grid(species_center(current_center)),&
                        limit_l_dim_SPE_species(species_center(current_center)),&
                        aux_R_prec_deriv_result_at_zero)

                    aux_R_prec_deriv_result_at_zero(1:limit_l_dim_SPE_species(species_center(current_center))) =&
                        aux_R_prec_deriv_result_at_zero(1:limit_l_dim_SPE_species(species_center(current_center)))*log_weight

                    !SR: do we have to multiply with log_weight after setting the components?
                    if ( i_center_2 .eq. i_center ) then
                        aux_R_prec_deriv_result_at_zero(2) = R_multipole_mpb_deriv_at_zero(1,current_spl_atom)
                        aux_R_prec_deriv_result_at_zero(3) = R_multipole_mpb_deriv_at_zero(2,current_spl_atom)
                        aux_R_prec_deriv_result_at_zero(4) = R_multipole_mpb_deriv_at_zero(3,current_spl_atom)
                    end if

                    delta_v_MPB_multipole_deriv_at_zero(:,current_spl_atom,current_spl_atom_2) = aux_R_prec_deriv_result_at_zero(:)

                    call tab_single_gradient_ylm_p0( trigonom_tab,limit_l_max_SPE_species, &
                        limit_l_max_SPE,current_center,ylm_tab,dylm_dtheta_tab,scaled_dylm_dphi_tab)

                    aux_R_prec_gradient_result = 0d0

                    call evaluate_g_gradient( limit_l_dim_SPE_species(species_center(current_center)), &
                        limit_l_dim_SPE_species_current(species_center(current_center)),&
                        1.0d0/sqrt_dist_tab_sq, dir_tab_in, trigonom_tab, limit_l_max_SPE, &
                        ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab, aux_R_prec_result_at_zero, &
                        aux_R_prec_deriv_result_at_zero, aux_R_prec_2nd_deriv_result, &
                        aux_R_prec_gradient_result,limit_l_max_SPE_species(species(current_spl_atom)),&
                        aux_R_prec_laplacian_result,index_lm_mpb,.FALSE.)

                    do i_coords = 1, 3, 1
                        R_prec_gradient_at_zero(i_coords,current_spl_atom_2) = &
                            R_prec_gradient_at_zero(i_coords,current_spl_atom_2) + aux_R_prec_gradient_result(i_coords)
                    end do

                end if !if MPB forces


                if (mpb_forces_on) then
                    !this is the far field gradient at the position of the nuclei
                    !indeed we only need the far field part of this, as calculated below,
                    !but just to be save that not the far field radii cutoffs
                    !of SPE_solver and aims interfere, we also add the correct vale here
                    delta_v_MPB_gradient_at_zero(:,current_spl_atom,current_spl_atom_2) = aux_R_prec_gradient_result
                end if


            else if (use_analytic_far_field) then ! far field calculation
                call tab_single_wave_ylm_p0(current_center,trigonom_tab, &
                    limit_l_max_SPE_species_ff,limit_l_max_SPE_ff, ylm_tab_ff )

                if (use_partition) then
                    kappa_times_rgrid = kappa_SPE_step*radius
                    if (mpb_forces_on) then
                        call bessel_K_deriv(kappa_times_rgrid,kappa_SPE_step,&
                            limit_l_max_SPE_ff, green_K_deriv)
                    end if
                else
                    kappa_times_rgrid = kappa_SPE*radius
                    if (mpb_forces_on) then
                        call bessel_K_deriv(kappa_times_rgrid,kappa_SPE,&
                            limit_l_max_SPE_ff, green_K_deriv)
                    end if
                end if
                call bessel_K_mpb(kappa_times_rgrid,limit_l_max_SPE_ff, green_K_current)
                if (mpb_forces_on) then
                    green_K_deriv = green_K_deriv*sqrt(pi/(2.*kappa_times_rgrid))-&
                        green_K_current*sqrt(pi/(8.*kappa_times_rgrid*radius**2))
                end if
                green_K_current = green_K_current*sqrt(pi/(2.*kappa_times_rgrid))


                do i_l = 0,limit_l_max_SPE_species_ff(species(current_spl_atom)),1
                    do i_m = -i_l, i_l, 1
                        aux_R_prec_result_at_zero_ff(index_lm_mpb(i_m,i_l)) = green_K_current(i_l+1)* &
                            multipole_moment_ff(current_spl_atom,index_lm_mpb(i_m,i_l))
                    end do
                end do
                R_prec_at_zero(current_spl_atom_2) = R_prec_at_zero(current_spl_atom_2) +&
                    ddot ( limit_l_dim_SPE_species_ff(species_center(current_center)), &
                    aux_R_prec_result_at_zero_ff(1:limit_l_dim_SPE_ff), 1, ylm_tab_ff(1:limit_l_dim_SPE_ff), 1 )


                if (mpb_forces_on) then

                    radius_rad = invert_radial_grid(radius,                                &
                        n_radial(species_center(current_center)), &
                        scale_radial(species_center(current_center)))

                    do i_l = 0, limit_l_max_SPE_species_ff(species(current_spl_atom)), 1
                        do i_m = -i_l, i_l, 1
                            aux_R_prec_deriv_result_at_zero_ff(index_lm_mpb(i_m,i_l)) = green_K_deriv(i_l+1)* &
                                multipole_moment_ff(current_spl_atom,index_lm_mpb(i_m,i_l))
                        end do
                    end do

                    dir_tab_in = dir_tab / radius
                    call tab_single_gradient_ylm_p0( trigonom_tab,  limit_l_max_SPE_species_ff, &
                        limit_l_max_SPE_ff, current_center, &
                        ylm_tab_ff, dylm_dtheta_tab_ff, scaled_dylm_dphi_tab_ff)

                    call evaluate_g_gradient(limit_l_dim_SPE_species_ff(species_center(current_center)), &
                        limit_l_dim_SPE_species_ff(species_center(current_center)), &
                        1.0d0/radius, dir_tab_in, trigonom_tab, &
                        limit_l_max_SPE_ff, ylm_tab_ff(1:limit_l_dim_SPE_ff),&
                        dylm_dtheta_tab_ff(1:limit_l_dim_SPE_ff), &
                        scaled_dylm_dphi_tab_ff(1:limit_l_dim_SPE_ff),&
                        aux_R_prec_result_at_zero_ff(1:limit_l_dim_SPE_ff),&
                        aux_R_prec_deriv_result_at_zero_ff(1:limit_l_dim_SPE_ff), &
                        aux_R_prec_2nd_deriv_result_ff(1:limit_l_dim_SPE_ff), &
                        aux_R_prec_gradient_result,&
                        limit_l_max_SPE_species_ff(species(current_spl_atom)),&
                        aux_R_prec_laplacian_result,index_lm_mpb,.false.)
                    do i_coords = 1, 3, 1
                        R_prec_gradient_at_zero(i_coords,current_spl_atom_2) = &
                            R_prec_gradient_at_zero(i_coords,current_spl_atom_2) + &
                            aux_R_prec_gradient_result(i_coords)
                    end do
                    delta_v_MPB_gradient_at_zero(:,current_spl_atom,current_spl_atom_2) = aux_R_prec_gradient_result
                end if !mpb_forces_on

            end if !far field/near field

        end do

    end subroutine evaluate_potential_at_zero

    !------------------------------------------------------------------------------
    !****s* SPE_solver/solve_SPE
    !  NAME
    !    solve_SPE
    !  SYNOPSIS
    subroutine solve_SPE(evaluate_only_mpe, R_in, R_out, R_out_gradient, &
            use_partition,kappa_SPE_step,dielec_func_current, change_rho_sign, &
            epsinfinv_factor,set_rho_mp,full_iteration,&
            total_steps,evaluate_laplacian, R_out_laplacian,lmax_cut,use_mpbe_params)
    !  INPUTS
    !   o R_in (input)   - on input: q - pi4_inv * L1 delta v_{n+1}
    !  OUTPUT
    !   o R_out (output)  - on output: delta v_{n+1} = delta_v_MPB
    !   o R_out_gradient (output)  - on output: delta_v_MPB_gradient
    !  SOURCE
    !  PURPOSE
    !    utilize precondition_kerker to solve the modified Helmholtz equation or &
        !   screened Poisson equation (SPE) and obtain
    !   the solution + the gradient
    !
    !  USES
        use localorb_io
        use dimensions, only: n_species, n_atoms, n_full_points, n_spline_atoms
        use grids, only: n_radial,r_radial
        use geometry, only: species
        use runtime_choices, only: extra_adding_to_hartree_potential_distance
        use synchronize_mpi_basic, only : sync_logical_vector
        use lpb_solver_utilities, only: multipole_radius_SPE_sq, changed_SPE_ff_cutoffs
        use mpi_tasks, only: aims_stop, myid
        implicit none

    !  ARGUMENTS

        !INPUTS:
        logical :: kerker_gradient, change_rho_sign, set_rho_mp, epsinfinv_factor, full_iteration
        logical :: evaluate_only_mpe !if true evaluate only mp expansion of input density not of integral
        !use_partition ! if true, divide eps/kappa space into pieces and sum up solutions later
        !epsinfinv_factor ! if true: normal calculation, if false: for SCLPB method
        logical :: use_mpbe_params !shall we solve MPBE or LPBE?
        real*8 :: kappa_SPE_step
        
        !variables concerning partitioning of kappa/eps space into equal pieces. this can be used
        !in case convergence with ions is far from optimal, then we can start from the ion-free case
        !and gradually increase the ionic strength
        logical :: use_partition
        !if use_partition, this is the dielectric permittivity to use:
        real*8 :: dielec_func_current
        
        ! input source term from previous Newton/MERM cycle
        real*8, dimension(n_full_points), intent(in)      :: R_in                

        !OUTPUTS
        ! output potential solving the SPE or LPB-type of PDE
        real*8, dimension(n_full_points), intent(out)      :: R_out                
        real*8, dimension(3,n_full_points), intent(out)      :: R_out_gradient                ! gradient of output residual from previous 
        
        !if needed, the laplacian of the potential:
        logical :: evaluate_laplacian
        real*8, dimension(n_full_points) :: R_out_laplacian

        integer, intent(out) :: total_steps

        !INTERNAL VARS
        logical :: set_mult_rad_to_addtomrad = .false. ! if true addtomprad is the actual radius,
        logical :: use_analytic_far_field = .true.

        integer :: i_species, i_atom
        real*8 ::  multipole_thresh
        real*8, dimension(n_full_points)      :: R_in_copy                ! input residual from previous scf cycle
        real*8, dimension(n_species) :: addtomprad!100.0
        real*8,  dimension(n_spline_atoms,limit_l_dim_SPE_ff) :: multipole_moment_ff
        character*1000 :: info_str


        real*8, dimension(:), allocatable :: R_out_laplacian_tot
        
        !if full_iteration
        integer :: i_iterate_lpb,i,j
        real*8 :: rmsd

        integer :: lmax_cut
        
        real*8, dimension(n_full_points) :: R_temp 
        real*8, dimension(n_full_points) :: R_temp_old 
        real*8, dimension(n_spline_atoms) :: rmsd_atom
        
        R_temp= 0d0
        R_temp_old = 0d0

        !these are the summed up atomic components for atomic_MERM = True
        atom_converged = .False.
        
        if (evaluate_laplacian .and. atomic_MERM) then
            if (.not. allocated(R_out_laplacian_tot)) allocate(R_out_laplacian_tot(n_full_points))
            R_out_laplacian_tot = 0d0
        end if
        
        ! FIXME: make explicit parameter if necessary - threshold to determine extend of multipole expansion
        multipole_thresh = 1d-15 !12

        if (debug_lpb_solver) then
            if (myid==0) write(use_unit,*) 'multipole_radius_free_SPE', multipole_radius_free_SPE
        end if
        if (use_analytic_far_field) then
            set_mult_rad_to_addtomrad = .true.
            do i_species = 1, n_species, 1
                if (changed_SPE_ff_cutoffs) then
                    addtomprad(i_species) = multipole_radius_SPE_sq
                else
                    addtomprad(i_species) =  (multipole_radius_free_SPE(i_species) &
                        +extra_adding_to_hartree_potential_distance)**2
                end if
            end do
        end if


        ! ! 	BEGIN DEBUG: Greens functions output
        !   do i_species = 1, n_species
        !      do i_grid = 1, n_grid(i_species)
        !
        !         ! calculate all l-values at once - for the radial argument (q_0 r)
        !
        !         ! FIXME: investigate whether or not is is possible to switch arguments 1 and 2 in green_I and green_K for speed reasons?
        ! 	kappa_times_rgrid = kappa_SPE*r_grid(i_grid,i_species)
        !         call bessel_K_deriv(kappa_times_rgrid,kappa_SPE,&
            ! 	  limit_l_max_SPE, green_K_deriv(i_grid,:,i_species))
        ! 	! now calculate modified spherical bessel function of 1st and 3rd kind by i_n = sqrt(pi/(2precond r)*I_n+1/2 and k_n = sqrt(pi/(2precond r)*K_n+1/2
        !         green_K_deriv(i_grid,:,i_species) = green_K_deriv(i_grid,:,i_species)*&
            ! 	  sqrt(pi/(2.*kappa_times_rgrid))-&
            ! 	  green_K_mpb(i_grid,:,i_species)*1./sqrt(pi/(2.*kappa_times_rgrid))*&
            ! 	  sqrt(pi/(8.*kappa_times_rgrid*r_grid(i_grid,i_species)**2))
        !      end do         ! radial loop, initialization of bessel green function
        !       if (i_species == 1) then
        ! 	do i_l = 0, limit_l_max_SPE_species(i_species), 1
        ! 	  do i_grid = 1, n_grid(i_species), 1
        ! 	    if (myid==0) then
        ! 	      write(use_unit,*) i_l, i_grid, r_grid(i_grid,i_species),&
            ! 		green_K_mpb(i_grid,i_l+1,i_species)
        ! 	    end if
        ! 	  end do
        ! 	end do
        !       end if
        ! !       if (myid ==0) then
        ! ! 	write(use_unit,*) 'change'
        ! !       end if
        ! !       if (i_species==1) then
        ! ! 	do i_l = 0, limit_l_max_SPE_species(i_species), 1
        ! ! 	  do i_grid = 1, n_grid(i_species), 1
        ! ! 	    if (myid==0) then
        ! ! 	      write(use_unit,*) i_l, i_grid, r_grid(i_grid,i_species),&
            ! ! 		green_K_deriv(i_grid,i_l+1,i_species)
        ! ! 	    end if
        ! ! 	  end do
        ! ! 	end do
        ! !       end if
        !   end do            ! species loop, initialization of bessel fct Green function
        !
        !    call mpi_barrier(mpi_comm_global,info)
        !    stop
        ! 	END DEBUG

        if (debug_lpb_solver.and.myid==0) then
            do i_atom = 1, n_atoms, 1
                write(use_unit,*) 'multipole_radius_free_SPE',multipole_radius_free_SPE(species(i_atom))
                write(use_unit,*) 'r_radial',r_radial(n_radial(species(i_atom)),species(i_atom))
            end do
        end if

        !   do i_atom = 1, n_atoms, 1
        !     multipole_radius_free_SPE(species(i_atom)) = r_radial(n_radial(species(i_atom))-1,species(i_atom))
        !     write(use_unit,*) 'iatom',i_atom,r_radial(n_radial(species(i_atom))-1,species(i_atom))
        !   end do

        !DEBUG: mp radius, end of radial grid
        if (debug_lpb_solver.and.myid==0) then
            write(use_unit,*) 'multipole radius LPB', multipole_radius_free_SPE
            write(use_unit,*) 'rradial', shape(r_radial)
        end if

        kerker_gradient = .True.
        !   only_electronic = .false. !neglect the nuclei potential completely

        if (full_iteration) then
            if (.not.initialized_lpb_solver) then
                !for the case this is called the first time, we do not have any splines
                !in memory, yet, thus initialize with zero
                R_multipole_mpb_old = 0d0
                R_multipole_spl_mpb  = 0d0
                multipole_moment_ff_glob = 0d0
                multipole_moment_ff = 0d0
                R_prec_old = 0d0
                R_prec_gradient_old = 0d0
                R_prec_old_atoms = 0d0
                R_prec_gradient_old_atoms = 0d0
!                 if (set_rho_mp) then
!                     rho_multipole_mpb_old = 0d0
!                 end if
            end if
            R_in_multipole = 0d0
            if (set_rho_mp) then
                rho_multipole_mpb_in = 0d0
                rho_multipole_mpb = 0d0
            end if
            R_out = 0d0
            R_out_gradient = 0d0
        else
            multipole_moment_ff = 0d0
            R_multipole_mpb      = 0d0                                         ! better safe than sorry ...
            R_multipole_spl_mpb  = 0d0
            R_prec_multipole_mpb = 0d0
            rho_multipole_mpb = 0d0
        end if


        ! multipole expansion is done according to Delley (1990), modifying Eqn (11) according to this particular purpose:
        ! R_multipole(s) = Int d\Omega ylm(r-r_atom) hartree_partition_tab R_in(r-r_atom)
        ! i.e.  R_in = partition_function * R_total already contains the partitioned residual for each atom,
        !       this routine will only have to play with the atomic contributions and add up all
        !       points in a given atomic integration shell later.

        !write(use_unit,*) '	|Solve LPB, partition:', use_partition
        if (debug_lpb_solver) then
            if (use_partition.and.myid==0.and.solve_lpbe_only) then
                write(info_str,'(2X,A,ES14.3,2X,A,F6.3)') '|MPE with lambda = ', kappa_SPE_step, 'dielec_func_current = ', dielec_func_current
                call localorb_info(info_str)
            else if (myid==0.and.solve_lpbe_only) then
                write(info_str,'(2X,A,ES14.3,2X,A,F6.3)') '|MPE with lambda = ', kappa_SPE, 'dielec_func_current = ', epsinf_mpb
                call localorb_info(info_str)
            end if
        end if
        if (change_rho_sign) then !if not already changed before calling routine, change to electronic charge density
            R_in_copy = -R_in
        else
            R_in_copy = R_in
        end if

        if (evaluate_only_mpe) then
            !1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!MPE OF INPUT DENSITY, ANGULAR INTEGRATION!!!!!!!!!!!!!!!!!!!!
            call evaluate_density_mpe(R_in_copy,dielec_func_current,&
                set_rho_mp,use_partition,.false.,0, R_out, R_out_gradient)
            !set on output: R_multipole_mpb and rho_multipole_mpb (results of angular integrals)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!WE USE THIS ROUTINE JUST TO SPLINE THE MP MOMENTS AND MANIPULATE THE SPLINES, NO RADIAL INTEGRATION!!!!!!!!!!!!!!!!!!!
            call evaluate_radial_integral(set_rho_mp, use_partition, use_analytic_far_field,&
                kappa_SPE_step, multipole_moment_ff,full_iteration, evaluate_only_mpe)
            !set on output: R_multipole_spl_mpb (splined mp moments), R_multipole_mpb_at_zero,
            !		multipole_moment_ff (radial multipole functions)
            !modified: rho_multipole_mpb
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!PERFORM THE ACTUAL MP SUMMATION TO GET POTENTIAL AND GRADIENT ON AIMS GRID
            call sum_up_mpe_potential(epsinfinv_factor,use_partition, kerker_gradient, &
                use_analytic_far_field, set_mult_rad_to_addtomrad, &
                addtomprad,multipole_moment_ff,kappa_SPE_step,&
                multipole_thresh, R_out, R_out_gradient,0,.false.,evaluate_laplacian,&
                R_out_laplacian,evaluate_only_mpe,lmax_cut,mpb_forces_on)
            !set on output: R_out (solution to SPE), R_out_gradient (gradient of R_out)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else !if .not.evaluate_only_mpe
            if (full_iteration) then !we solve the LPBE
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !	the following scheme is applied:
                !		in the beginning of each lpb loop, we need to calculate the constant part of rho_iter, which is independent of u but on rho: R_in_multipole
                !>>> 1. if lpb is called first time: (initialized_lpb_solver = .false.)
                !------------------------------------------------------------------------
                !ooo i_iterate_lpb = 1 !!!!!!!!!!!!REMOVE UNNESSESARY INITIALIZATION RUN
                !- evaluate_density_mpe (with rho_scaled) -> R_in_multipole = R_multipole_mpb
                !- R_multipole_mpb_old = R_multipole_mpb, R_prec_old = R_out
                !- evaluate_radial_integral -> R_multipole_spl_mpb, R_multipole_at_zero, multipole_moment_ff_glob
                !ooo i_iterate_lpb > 1
                !- sum_up_mpe_potential -> R_out, R_out_gradient
                !- evaluate_density_mpe (with rho_iter) -> R_multipole_mpb
                !- R_multipole_mpb += R_in_multipole, linear mixing with R_multipole_mpb_old, evaluate rmsd
                !(- if converged: sum_up_mpe_potential for other points -> R_out, R_out_gradient)
                !- R_multipole_mpb_old = R_multipole_mpb, R_prec_old = R_out
                !- evaluate_radial_integral -> R_multipole_spl_mpb, R_multipole_at_zero, multipole_moment_ff_glob
                !>>> 2. if lpb is called later times: (initialized_lpb_solver = .true.)
                !------------------------------------------------------------
                !o i_iterate_lpb = 1
                !- evaluate_density_mpe (with rho_scaled) -> R_in_multipole = R_multipole_mpb
                !- R_out = R_prec_old (initialize with last scf step)
                !- evaluate_density_mpe (with rho_iter) -> R_multipole_mpb
                !- R_multipole_mpb += R_in_multipole, linear mixing with R_multipole_mpb_old, evaluate rmsd
                !(- if converged: sum_up_mpe_potential for other points -> R_out, R_out_gradient)
                !- R_multipole_mpb_old = R_multipole_mpb, R_prec_old = R_out
                !- evaluate_radial_integral -> R_multipole_spl_mpb, R_multipole_at_zero, multipole_moment_ff_glob
                !i_iterate_lpb > 1
                !	=> same as above
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                total_steps = 0
                    
                !perform iterations for each atom separately

                !some initializations
                lpb_converged = .False.
                if (.not.initialized_lpb_solver) then
                    !for the case this is called the first time, we do not have any splines
                    !in memory, yet, thus initialize with zero
!                             if (atomic_MERM) then
!                                 !do this only for the current atom, all previous atoms have already been calculated and stored here
!                                 R_multipole_mpb_old(:,:,i_atom) = 0d0
!                                 multipole_moment_ff_glob(i_atom,:) = 0d0
!                             else
                    R_multipole_mpb_old = 0d0
                    multipole_moment_ff_glob = 0d0
!                             end if
                    R_multipole_spl_mpb  = 0d0
                    multipole_moment_ff = 0d0
                    R_prec_old = 0d0
                    R_prec_old_atoms = 0d0
                    R_prec_gradient_old = 0d0
                    R_prec_gradient_old_atoms = 0d0
!                         if (set_rho_mp) then
!                             rho_multipole_mpb_old = 0d0
!                         end if
                end if
                R_in_multipole = 0d0
                if (set_rho_mp) then
                    rho_multipole_mpb_in = 0d0
                    rho_multipole_mpb = 0d0
                end if
                R_out = 0d0
                R_out_gradient = 0d0

                i_iterate_lpb = 1
                rmsd = 1d3
                rmsd_atom = rmsd
                if (myid==0.and.solve_lpbe_only) then
!                     if (i_atom==1) then
                    write(use_unit,'(2X,A)') "|Starting iterative solver using spatially partial update of potential in each step"
!                     end if
!                     if (atomic_MERM) then
!                         write(use_unit,'(2X,A,I6)') "| Atom: ", i_atom
!                     end if
                end if
                do while (.not.(lpb_converged))

                    R_prec_multipole_mpb = 0d0
                    R_multipole_mpb      = 0d0

                    !for the first iteration step, perform multipole expansion of rho_in (input density), since rho has changed
                    if (i_iterate_lpb==1) then
                        call evaluate_density_mpe(R_in_copy,dielec_func_current,set_rho_mp,&
                            use_partition,full_iteration,0, R_out, R_out_gradient)

                        !set on output: R_multipole_mpb and rho_multipole_mpb and saves it as rho_multipole_mpb_in and R_in_multipole
                    else !if i_iterate_lpb.gt.1
                        if (use_mpbe_params) then
                            call sum_up_mpe_potential(epsinfinv_factor,use_partition, kerker_gradient, &
                                use_analytic_far_field, set_mult_rad_to_addtomrad, &
                                addtomprad,multipole_moment_ff,kappa_SPE_step,&
                                multipole_thresh, R_out, R_out_gradient,2,.false.,&
                                evaluate_laplacian, R_out_laplacian,evaluate_only_mpe,lmax_cut,mpb_forces_on)
                        else
                            call sum_up_mpe_potential(epsinfinv_factor,use_partition, kerker_gradient, &
                                use_analytic_far_field, set_mult_rad_to_addtomrad, &
                                addtomprad,multipole_moment_ff,kappa_SPE_step,&
                                multipole_thresh, R_out, R_out_gradient,1,.false.,&
                                evaluate_laplacian, R_out_laplacian,evaluate_only_mpe,lmax_cut,mpb_forces_on)
                        end if
                        !set on output: R_out (solution to Helmholtz equation), R_out_gradient (gradient of
                        !solution to Helmholtz equation
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    end if !i_iterate_lpb.eq.1

                    if (i_iterate_lpb.gt.1.or.initialized_lpb_solver) then
                        if (initialized_lpb_solver.and.i_iterate_lpb.eq.1) then
                            ! in this case we use the potential of the end of the last MERM cycle to initialize
                            R_out = R_prec_old
                            R_out_gradient = R_prec_gradient_old
                            R_multipole_mpb = 0d0
                            if (set_rho_mp) then
                                rho_multipole_mpb = 0d0
                            end if
                        end if
                        if (use_mpbe_params) then
                            call evaluate_density_mpe(R_in_copy,dielec_func_current,set_rho_mp,&
                                use_partition,full_iteration,2, R_out, R_out_gradient)
                        else
                            call evaluate_density_mpe(R_in_copy,dielec_func_current,set_rho_mp,&
                                use_partition,full_iteration,1, R_out, R_out_gradient)
                        end if

                    end if

                    if (i_iterate_lpb.gt.1.or.(initialized_lpb_solver.and.i_iterate_lpb.eq.1)) then
                        !calculate full R_multipole_mpb by adding R_in_multipole
                        R_multipole_mpb = R_multipole_mpb + R_in_multipole

                        !linear mixing
                        R_multipole_mpb = R_multipole_mpb_old + eta_merm*(R_multipole_mpb-R_multipole_mpb_old)

                        if (set_rho_mp) then
                            !in case of evaluating the MP expanded source term, we need to add
                            !rho_in. But only if we did not use an old rho_multipole_mpb of an
                            !already converged atom. In that case, we already have the correct
                            !potential
                            if (atomic_MERM) then 
                                !only add rho_multipole_mpb_in for those atoms for which rho_multipole_mpb was
                                !actually calculated in this iteration step
                                !for the already converged atoms, rho_multipole_mpb should not be touched!
                                do i_atom = 1, n_spline_atoms, 1
                                    if (.not. atom_converged(i_atom)) then
                                        rho_multipole_mpb(:,:,i_atom) = rho_multipole_mpb(:,:,i_atom) +&
                                            rho_multipole_mpb_in(:,:,i_atom)
                                    end if
                                end do
                            else
                                rho_multipole_mpb = rho_multipole_mpb + rho_multipole_mpb_in
                            end if
                        end if
                    end if

                    if (i_iterate_lpb.gt.1) then
                        !we always have to calculate at least one time the potential from the current density to check the response
                        !of the solvent -> i_iterate_lpb must be larger than 1
                        ! 	if (myid==0) then
                        ! 	  write(use_unit,*) 'evaluating rmsd', i_iterate_lpb
                        ! 	end if
                        !SR: check this convergence criterion here for atomic_MERM
                        ! we should switch to an atom-wise testing of rmsd. only continue evaluating the atoms which 
                        ! have not converged, yet.
                        if (.not. atomic_MERM) then
                            call evaluate_rmsd_merm(R_out(1), R_prec_old(1), rmsd, tau_merm, lpb_converged)
                        else 
                            !we check for each atom to be converged. afterwards continue iterating only the 
                            !not yet converged atomic potentials
                            do i_atom=1, n_spline_atoms, 1
                                !if (myid.eq.task_list(i_atom)) then
                                    if (.not. atom_converged(i_atom)) then
                                        !get 2 arrays which contain the data of the respective atom only and zero's otherwise 
                                        !in order to easily apply the routine evaluate_rmsd_merm
                                        R_temp = 0d0
                                        R_temp_old = 0d0
                                        do j = 1, index_atoms_size(i_atom), 1 
                                            i = index_atoms(i_atom,j)
                                            R_temp(i) = R_out(i)
                                            R_temp_old(i) = R_prec_old(i)
                                        end do
                                        call evaluate_rmsd_merm(R_temp(1), R_temp_old(1), rmsd_atom(i_atom), tau_merm, atom_converged(i_atom))
                                        if (atom_converged(i_atom)) then
                                            if (myid==0) then
                                                write(use_unit,'(2X,A,I4,A,I4,A)') "|Atom ", i_atom, " converged in ", i_iterate_lpb, " steps."
                                            end if
                                        end if                                    
                                    end if
                                !end if
                            end do
                            !call sync_logical_vector(atom_converged,n_spline_atoms,SYNC_OR)
                        end if
                        if (atomic_MERM) then
                            rmsd = maxval(rmsd_atom)
                        end if
                                    
                        if ((rmsd.lt.tau_merm .and. .not. atomic_MERM) .or.&
                             (all(atom_converged) .and. atomic_MERM)) then
                            lpb_converged = .True.
                            if (solve_lpbe_only) then
                                write(info_str,'(2X,A)') '|Self-consistency cycle converged '
                                call localorb_info(info_str,use_unit,'(A)',OL_norm)
                                write(info_str,'(2X,A,ES14.3)') '|Final RMSD: ', rmsd
                                call localorb_info(info_str,use_unit,'(A)',OL_norm)
                            end if
                            !in this case
                            !we will evaluate the splines only in the outer part for one last time
                            !if not wished, we will integrate the cLPBE one more time and in the end evaluate the splines on the full grid
                            !we have to sum up the multipoles to get the new potential and gradient
                            !we have to sum up the multipoles to get the new potential and gradient
                            !3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            !!!PERFORM THE ACTUAL MP SUMMATION TO GET POTENTIAL AND GRADIENT ON AIMS GRID
                            if (atomic_MERM) then
                                R_prec_old = R_out
                                R_prec_gradient_old = R_out_gradient
                            end if
                            
                            if (use_mpbe_params) then
                                call sum_up_mpe_potential(epsinfinv_factor,use_partition, kerker_gradient, &
                                    use_analytic_far_field, set_mult_rad_to_addtomrad, &
                                    addtomprad,multipole_moment_ff,kappa_SPE_step,&
                                    multipole_thresh, R_out, R_out_gradient,2,.true.,evaluate_laplacian, &
                                    R_out_laplacian,evaluate_only_mpe,lmax_cut,mpb_forces_on)
                            else
                                call sum_up_mpe_potential(epsinfinv_factor,use_partition, kerker_gradient, &
                                    use_analytic_far_field, set_mult_rad_to_addtomrad, &
                                    addtomprad,multipole_moment_ff,kappa_SPE_step,&
                                    multipole_thresh, R_out, R_out_gradient,1,.true.,evaluate_laplacian, &
                                    R_out_laplacian,evaluate_only_mpe,lmax_cut,mpb_forces_on)
                            end if
                            
                            if (.not. atomic_MERM) then
                                R_prec_old = R_out
                                R_prec_gradient_old = R_out_gradient
                            end if                                


                            !set on output: R_out (solution to Helmholtz equation), R_out_gradient (gradient of
                            !solution to Helmholtz equation
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            total_steps = total_steps + i_iterate_lpb-1
                        else
                            lpb_converged = .false.
                        end if
                        if (.not. lpb_converged) then
                            if (rmsd.gt.1d6) then
                                call aims_stop('RMSD was over 1d6, SCLPBE cycle may not converge. Reconsider your grid and angular momentum settings.')
                            end if
                            if (solve_lpbe_only.and..not.atomic_MERM) then
                                write(info_str,'(2X,A,ES14.3)') '|Current RMSD: ', rmsd
                                call localorb_info(info_str,use_unit,'(A)',OL_norm)
                            end if
                        end if
                    end if !i_iterate_lpb.gt.1
                    if (.not.lpb_converged) then
                        if (atomic_MERM) then
                            do i_atom = 1, n_atoms, 1
                                if (.not. atom_converged(i_atom)) then
                                    R_multipole_mpb_old(:,:,i_atom) = R_multipole_mpb(:,:,i_atom)
                                end if
                            end do
                        else
                            R_multipole_mpb_old = R_multipole_mpb
                        end if
                        R_prec_old = R_out
                        R_prec_gradient_old = R_out_gradient
                        if (i_iterate_lpb.eq.1.and.solve_lpbe_only) then
                            write(info_str,'(2X,A)') '|Multipole density initialization finished.'
                            call localorb_info(info_str,use_unit,'(A)',OL_norm)
                        end if
                        if (solve_lpbe_only .and. .not.atomic_MERM) then
                            write(info_str,'(2X,A,I4)') '|Solving LPB equation with step: ', i_iterate_lpb
                            call localorb_info(info_str,use_unit,'(A)',OL_norm)
                        end if

                        !!!SPLINE INTERPOLATION ON FINE LOG GRID AND RADIAL INTEGRATION!!!!!!!!!!!!!!!!!!!!
                        call evaluate_radial_integral(set_rho_mp, use_partition, use_analytic_far_field,&
                            kappa_SPE_step,multipole_moment_ff,full_iteration, evaluate_only_mpe)
                        !set on output: R_multipole_spl_mpb (splined mp moments), R_multipole_mpb_at_zero,
                        !		multipole_moment_ff
                        i_iterate_lpb = i_iterate_lpb + 1
                        SPE_first_warnings = .false.
                    end if
                end do
                total_steps = i_iterate_lpb
!                         if (atomic_MERM) then
!                             !add up the atomic components
!                             R_out_tot = R_out_tot + R_out
!                             R_out_gradient_tot = R_out_gradient_tot + R_out_gradient
!                             R_prec_at_zero_tot = R_prec_at_zero_tot + R_prec_at_zero
!                             R_prec_gradient_at_zero_tot = R_prec_gradient_at_zero_tot + R_prec_gradient_at_zero
!                             if (evaluate_laplacian) then
!                                 R_out_laplacian_tot = R_out_laplacian_tot + R_out_laplacian
!                             end if
!                         end if
!                          if (.not. atomic_MERM .or. (i_atom .eq. n_atoms)) then
!                     if (atomic_MERM) then
!                         if (evaluate_laplacian) then
!                             R_out_laplacian = R_out_laplacian_tot
!                         end if
!                     end if
                return
!                         end if
!                     end if !if current atom is calculated by this CPU
!                 end do !loop over atoms
                
!                 if (atomic_MERM) then
!                     ! we need to sync a few vectors/matrices
!                     call sync_vector(R_out,n_full_points)
!                     call sync_matrix(R_out_gradient,3,n_full_points)
!                     call sync_vector(R_prec_at_zero,n_spline_atoms)
!                     call sync_matrix(R_prec_gradient_at_zero,3,n_spline_atoms)
!                     if (evaluate_laplacian) then
!                         call sync_vector(R_out_laplacian, n_full_points)
!                     end if
!                     if (use_analytic_far_field) then
!                         call sync_matrix(multipole_moment_ff,n_spline_atoms,limit_l_dim_SPE_ff)
!                     end if           
!                 end if

            else !we just solve the LPBE with constant coefficients

                !1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!MPE OF INPUT DENSITY, ANGULAR INTEGRATION!!!!!!!!!!!!!!!!!!!!
                !eq. (42) get multipole moments of source term
                call evaluate_density_mpe(R_in_copy,dielec_func_current,set_rho_mp,&
                    use_partition,full_iteration,0, R_out, R_out_gradient)
                !set on output: R_multipole_mpb and rho_multipole_mpb
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!SPLINE INTERPOLATION ON FINE LOG GRID AND RADIAL INTEGRATION!!!!!!!!!!!!!!!!!!!!
                !eq. (41)
                call evaluate_radial_integral(set_rho_mp, use_partition,use_analytic_far_field,&
                    kappa_SPE_step,&
                    multipole_moment_ff,full_iteration, evaluate_only_mpe)
                !set on output: R_multipole_spl_mpb (splined mp moments), R_multipole_mpb_at_zero,
                !		multipole_moment_ff
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!PERFORM THE ACTUAL MP SUMMATION TO GET POTENTIAL AND GRADIENT ON AIMS GRID
                !eq. (40)
                call sum_up_mpe_potential(epsinfinv_factor,use_partition, kerker_gradient, &
                    use_analytic_far_field, set_mult_rad_to_addtomrad, &
                    addtomprad,multipole_moment_ff,kappa_SPE_step,&
                    multipole_thresh, R_out, R_out_gradient,0,.false.,evaluate_laplacian,&
                    R_out_laplacian,evaluate_only_mpe,lmax_cut,mpb_forces_on)
                !set on output: R_out (solution to SPE), R_out_gradient (gradient of
                !solution to SPE
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            end if !iterative or just one step
        end if !if only mpe of input density or of integral (i.e. solve Helmholtz equation)
        ! warnings should only be output once
        SPE_first_warnings = .false.
        

    end subroutine solve_SPE


    !------------------------------------------------------------------------------
    !****s* SPE_solver/evaluate_pot_lpb_at_zero
    !  NAME
    !    evaluate_pot_lpb_at_zero
    !  SYNOPSIS
    subroutine evaluate_pot_lpb_at_zero( lpb_zero_energy)
    !  PURPOSE
    ! calculates interaction energy between nuclear charge at electrostatic potential

        use grids           ! information on integration grids and such
        use dimensions      ! here each thread learns about its local memory structure
        use constants
        use spline
        use species_data
        use runtime_choices
        use pbc_lists
        use geometry, only: coords, species

        implicit none

        real*8 :: lpb_zero_energy
        integer :: i_center, current_center, current_spl_atom
        real*8 :: coord_current(3)
        real*8 :: delta_mpb_at_zero !mpb potential - hartree potential at position of nuclei

        lpb_zero_energy = 0.0d0

        !free_pot_es_at_zero depends on species(i_atom)


        do i_center = 1, n_centers_hartree_potential, 1
            current_center   = centers_hartree_potential(i_center)
            current_spl_atom = center_to_atom(current_center)
            coord_current(:) = coords(:,current_spl_atom)

            !this is the total lpb potential at the specified cutoff minus z/rcut (self interaction nuclei)

            delta_mpb_at_zero = species_z(species(current_center)) * R_prec_at_zero(current_spl_atom)

            lpb_zero_energy = lpb_zero_energy + delta_mpb_at_zero

            !      write(use_unit,*) 'myid, iatom: ', myid, current_spl_atom
            !      write(use_unit,*) 'center		', i_center, '	elec at free pot at zero:', - free_pot_es_at_zero(species(current_spl_atom))
            !      write(use_unit,*) 'center		', i_center, '	elec solvent tot at zero: ', R_prec_at_zero(current_spl_atom)
            !      write(use_unit,*) 'center		', i_center, '	elec delta solvent at zero: ', delta_mpb_at_zero
            !      write(use_unit,*) 'center		', i_center, '	nuclei solvent interaction: ', nuclei_solvent_int
        end do
        !   write(use_unit,*) 'R_prec_at_zero: ', myid,current_spl_atom,R_prec_at_zero
        !   write(use_unit,*) 'free_pot_es_at_zero',myid,free_pot_es_at_zero

        !i checked this energy for 1cpu compared to n cpu calculations for an n2 molecule
        !the difference is 1e-12 Ha, is this a bug or just uncertainty?
        !the error in the solvation energy MPI-serial is 10e-6 Hartree, bug?

        return

    end subroutine

    !******

    !------------------------------------------------------------------------------
    !****s* preconditioner/bessel_i_mpb
    !  NAME
    !    bessel_i_mpb
    !  SYNOPSIS

    subroutine bessel_I_mpb(r,Lmax,y)

    !  PURPOSE
    ! compute the modified spherical bessel functions I_(1/2)(r), I_(3/2)(r), ... , I_(Lmax+1/2)(r) and
    ! store in the array y(0:Lmax) in that order. Routine *should* be stable for the purposes of AIMS.
    ! method: two different cases, see below
    !  USES
        use constants
        implicit none

    !  ARGUMENTS

        integer Lmax
        double precision r

    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  INPUTS
    !   o r    - argument of the functions
    !   o Lmax - maximal angular momentum shell for which to compute
    !  OUTPUT
    !   o y(0:Lmax) - values of the Bessel functions
    !  SOURCE


        integer Ndiff, i, j
        double precision y(Lmax+1),r_thresh,y_next,y_this,y_prev,buf
        double precision term_thresh, thresh_used
        parameter (Ndiff = 100)
        parameter (r_thresh = 5d-1)
        parameter (term_thresh = 1d-100)
        ! two cases:

        if (r.ge.r_thresh) then
            ! large r >= r_thresh
            ! start Ndiff functions above Lmax and iterate downwards with the recursion formula
            !     I(n-1) = I(n+1) + 2 n I(n)/r (Abramowitz, Stegun, Eqn 10.2.18)
            ! this gets you within a constant factor of the actual result - the latter is obtained
            ! by comparing with I(0) = sqrt(2 pi/r) sinh(r) which evaluates quite nicely
            !
            ! in my empirical investigation, this fixed-point type evaluation converged for r >~ 0.1;
            ! r_thresh is set well above that value
            y_next = 1d-10
            y_this = 0d0
            do i = Lmax+Ndiff, Lmax, -1
                y_prev = y_next + 2*(dble(i+1)+5d-1)*y_this/r
                y_next = y_this
                y_this = y_prev
            end do
            y(Lmax+1) = y_next + 2*(dble(Lmax+1)+5d-1)*y_this/r
            if (Lmax.gt.0) then
                y(Lmax  ) = y_this + 2*(dble(Lmax)+5d-1)*y(Lmax+1)/r
                if (Lmax.gt.1) then
                    do i = Lmax-2, 0, -1
                        y(i+1) = y(i+3) + 2*(dble(i+1)+5d-1)*y(i+2)/r
                    end do
                end if
            end if
            buf = sinh(r)*sqrt(2d0/(pi*r))/y(1)
            !This is I_i+1/2, i.e. the modified Bessel function
            y(:) = buf*y(:)
        else
            ! small r <= r_thresh
            ! direct implementation of Abramowitz, Stegun; Eqn (10.2.5).
            ! the summation is cut off when the terms in the series are smaller than
            ! either r^6 or term_thresh (defined above), whichever is smaller.
            thresh_used = r**6
            if (thresh_used.gt.term_thresh) thresh_used = term_thresh
            do i = 0, Lmax
                buf = 1d0
                y(i+1) = buf
                j    = 1
                do while (buf.gt.term_thresh)
                    buf = buf*r*r/(2d0*dble(j)*dble(2*i+2*j+1))
                    y(i+1) = y(i+1) + buf
                    j = j + 1
                end do
                buf = 1d0
                do j = 0, i
                    buf = buf*(2d0*dble(j)+1d0)
                end do
                !This is I_i+1/2, i.e. the modified Bessel function
                y(i+1) = r**(dble(i)+5d-1)*y(i+1)*sqrt(2/pi)/buf
            end do
        end if
    end subroutine bessel_I_mpb
    !******

    !****s* preconditioner/bessel_K
    !  NAME
    !    bessel_K
    !  SYNOPSIS
    subroutine bessel_K_mpb(r,Lmax,y)
    !  PURPOSE
    ! compute the modified bessel functions I_(1/2)(r), I_(3/2)(r), ... , I_(Lmax+1/2)(r) and
    ! store in the array y(0:Lmax) in that order. Routine *should* be stable for the purposes of AIMS.
    ! calculation is based on I_(1/2)(r) and I_(3/2)(r) which is then paired with the usual upwards iteration
    !  USES
        use constants
        implicit none
    ! ARGUMENTS


        integer Lmax
        double precision r,y(Lmax+1)

    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  INPUTS
    !    o r    - argument of the functions
    !    o Lmax - maximal angular momentum shell for which to compute
    !  OUTPUT
    !    y(0:Lmax) - values of the Bessel functions
    !  SOURCE

        integer i

        y(1) = sqrt(pi/(2d0*r))*exp(-r)
        if (Lmax.ge.1) then
            y(2) = (1d0+1d0/r)*y(1)
            if (Lmax.ge.2) then
                do i = 1, Lmax-1
                    y(i+2) = y(i)+2d0*(dble(i)+5d-1)*y(i+1)/r
                end do
            end if
        end if

    end subroutine bessel_K_mpb

    !****s* preconditioner/bessel_K
    !  NAME
    !    bessel_K
    !  SYNOPSIS
    subroutine bessel_K_deriv(r,lambda,Lmax,dy)
    !  PURPOSE
    ! compute the modified bessel functions I_(1/2)(r), I_(3/2)(r), ... , I_(Lmax+1/2)(r) and
    ! store in the array y(0:Lmax) in that order. Routine *should* be stable for the purposes of AIMS.
    ! calculation is based on I_(1/2)(r) and I_(3/2)(r) which is then paired with the usual upwards iteration
    !  USES
        use constants
        implicit none
    ! ARGUMENTS


        integer Lmax
        double precision r,y(Lmax+1),dy(Lmax+1)
        double precision lambda

    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  INPUTS
    !    o r    - argument of the functions
    !    o Lmax - maximal angular momentum shell for which to compute
    !  OUTPUT
    !    y(0:Lmax) - values of the Bessel functions
    !  SOURCE

        integer i

        y(1) = sqrt(pi/(2d0*r))*exp(-r)
        dy(1) = -y(1)*(lambda+1d0/(2d0*r/lambda))
        if (Lmax.ge.1) then
            y(2) = (1d0+1d0/r)*y(1)
            dy(2) = (1d0+1d0/r)*dy(1)-lambda/r**2*y(1)
            !-y(2)*lambda*(3d0/(2d0*r)+r/(1d0+r))
            if (Lmax.ge.2) then
                do i = 1, Lmax-1
                    y(i+2) = y(i)+2d0*(dble(i)+5d-1)*y(i+1)/r
                    dy(i+2) = dy(i)+lambda*(dble(i)+5d-1)*&
                        (2d0*r/lambda*dy(i+1)-2d0*y(i+1))/r**2
                end do
            end if
        end if

    end subroutine bessel_K_deriv

    subroutine get_rho_multipole_spl_mpb(rho_multipole_spl, spl_atom)

    !  PURPOSE
    !    Delivers the spline coefficients for rho_multipole
    !  USES

        use species_data
        use dimensions
        use grids
        use pbc_lists       ! for atom-atom distance tab center_to_atom - which should work for PBC as well ...
        use spline
        use geometry, only: species

        implicit none

    !  ARGUMENTS

        real*8, intent(out) :: rho_multipole_spl(limit_l_dim_SPE, n_max_spline_mpb, n_max_radial+2)
        integer, intent(in) :: spl_atom

    !  INPUTS
    !    o spl_atom -- the atom for which the spline coefficients are needed
    !  OUTPUTS
    !    o rho_multipole_spl -- the spline coefficients
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE


        integer i_radial, l_h_dim, n_rad, j, i_lm

        real*8 :: i_r_outer
        real*8 :: delta, delta_2, delta_3


        l_h_dim = limit_l_dim_SPE_species(species(spl_atom))
        n_rad   = n_radial(species(spl_atom))

        ! First, simply produce spline coefficients from the charge density as
        ! given on the radial integration grid.

        do i_lm = 1, l_h_dim, 1
            do i_radial = 1, n_rad+2, 1
                rho_multipole_spl(i_lm,1,i_radial) = rho_multipole_mpb(i_radial,i_lm,spl_atom)
            end do
        end do

        ! Spline interpolation
        call cubic_spline_v2(rho_multipole_spl, limit_l_dim_SPE, n_max_spline_mpb, n_max_radial+2, &
            n_rad+2, l_h_dim)

        ! Splines are now tabulated up to r = infinity in principle.
        ! NOW, "doctor" all splines in the far field:
        ! We know that no density must occur outside the radius of the free atom,
        ! because the partition table is zero there.
        ! Therefore, extrapolate from last radial shell inside free atom radius to
        ! become zero at multipole_radius_free

        ! find outermost radial grid point that is possibly non-zero
        i_radial = n_rad
        do while ( ( r_radial(i_radial,species(spl_atom)) .ge. &
                multipole_radius_free_SPE(species(spl_atom)) ) &
                .and.(i_radial.gt.1) )

            rho_multipole_spl(1:l_h_dim,:,i_radial+1) = 0.d0

            i_radial = i_radial - 1
        enddo

        ! Outermost atom radius in units of the radial integration grid
        i_r_outer = invert_radial_grid &
            ( multipole_radius_free_SPE(species(spl_atom)), &
            n_rad, &
            scale_radial(species(spl_atom)) )

        delta = dble(i_r_outer - i_radial)
        delta_2 = delta*delta
        delta_3 = delta_2*delta

        ! This is an ugly hack because the element i_radial+1 in rho_multipole_spl
        ! now corresponds to radial shell r_radial(i_radial). What a mess.
        i_radial = i_radial + 1

        ! doctor the spline coefficients at the outermost finite value
        ! i_radial to go smoothly to zero at multipole_radius_free

        do j = 1, l_h_dim

            rho_multipole_spl( j, 3, i_radial) = &
                - 3.d0 / delta_2 * rho_multipole_spl( j, 1, i_radial) &
                - 2.d0 / delta   * rho_multipole_spl( j, 2, i_radial)

            rho_multipole_spl( j, 4, i_radial) = &
                2.d0 / delta_3 * rho_multipole_spl( j, 1, i_radial) &
                + 1.d0 / delta_2 * rho_multipole_spl( j, 2, i_radial)

        enddo

    end subroutine get_rho_multipole_spl_mpb

end module SPE_solver


