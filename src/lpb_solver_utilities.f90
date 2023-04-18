!****s* FHI-aims/lpb_solver_utilities
!  NAME
!   lpb_solver_utilities
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

module lpb_solver_utilities

!  PURPOSE
!
!
!  Module lpb_solver_utilities contains all variables & routines of the
!  MPBE solver which are not explicitely connected to the size-
!  modification and nonlinearity extensions to the LPBE.
!
!  VARIABLES
!
!  //GENERAL CONVERGENCE//
!	mpb_solver_started : determines if the MPB solver is executed already during each SCF step
!	lpb_converged : determines if LPB solver has been converged
!	mpb_converged : initially = false. prevents SCF cycle to stop until we actually run the MPB solver and
!		mpb_solver_started = true. after which time mpb_converged = true
!	mpb_start_crit : convergence criterium. if clc = 'converge', the SCF is first converged in vacuum before the MPB solver is started in each
!		SCF step. if clc = 'step', the MPB solver is started at the SCF step to be specified with
!	start_at_scf_step : # of SCF steps at which MPB solver is started
!
!  //SYSTEM SPECIFICATIONS//
!	solve_pbe_only : if true the PBE without size-modification is solved
!	solve_lpbe_only : solve LPBE as given in the main paper above
!	solve_debye_only : solve lpbe with a = 0
!	correct_dielectric_decrement : this is to correct for the dependency of  the dielectric function on the spatially dependent ion concentration.
!		so we make eps_s spatially dependent: eps_s(kappa_func**2*z/kT*phi) very experimental, if anything wants to be done with this please contact authors
! 	use_mpbe_free_energy : the LPBE potential is often a good approximation to the true MPBE potential. if this flag is true, use therefore the
!		more accurate MPBE free energy expression containing important entropic corrections but the LPBE potential instead of the MPBE
! 	use_separate_alpha  : use 2 alpha functions, one for anions and one for cations, experimental, please refer to authors
!	solve_lpb_with_constant_dielec : instead of solving the MPBE, solve the LPBE with spatially constant dielectric function and ionic strength experimental, only for tests
!	freeen_LPBE_nonelstat_solv_energy : nonelectrostatic LPBE solvation energy stored globally here
!	dieleckind : kind of dielectric function. kind = 0: dielectric function of  Gygi & Fattebert (not recommended since gradient is not zero outside
!		of transition region, kind = 1: Marzari/Andreussi function see paper above
!	alphakind : alphakind = 0 takes sharp step function for alpha function, alphakind=1  parameterizes in terms of the Marzari/Andreussi dielectric function (see paper)
!	/PARAMETERS DEFINING THE IONS/
!	offset_kappa_mpb : shift of alpha function (alphakind=0)
!	d_alpha_ion : shift alpha function in log density units (alphakind=1) shift of the kappa function. if positive in the direction of the solvent.
! 		1d0 means rho_max(eps) = rho_min(kappa). 0.5 means ln(rho_min(kappa)) =  (ln(rho_min(eps))+ln(rho_max(eps)))/2
!	xi_alpha_ion : scale alpha function in log density units (alphakind=1) symmetrically rescales the ln(rho_min(kappa))-ln(rho_max(kappa)) distance.
!	d_alpha_ion_anion, xi_alpha_ion_anion : if use_separate_alpha: define the anionic alpha function. normal alpha function is then automatically the cationic alpha function
!	dalpha_anion_drho : derivative of dalpha_anion with respect to nel
!	rhomin_kappa, rhomax_kappa,log_rhomin_kappa,log_rhomax_kappa : density cutoffs for alpha function (alphakind=1)
!	rhomin_kappa_anion, rhomax_kappa_anion,log_rhomin_kappa_anion,log_rhomax_kappa_anion : density cutoffs for anionic alpha function
!	z_mpb : ionic charge in MPB calculation
!	a_mpb : ion size
!	c_mpb : ionic bulk concentration
!	phi_zero_mpb : volume fraction occupied by ions
!	T_mpb, kBT_mpb : temperature and kB*T of ions
!	kappainf_mpb : sqrt(8 pi c z)
!	kappa_div_sqrt_eps : sqrt(8 pi c z / epsinf)
!	kappa_debye : sqrt(8 pi c z^2 / epsinf / kT)
!!	/PARAMETERS DEFINING THE SOLVENT/
!	alphaplusgamma_mpb, beta_mpb : parameter defining Gnonmf
!	epsinf_mpb,log_epsinf_mpb : dielectric function value in the bulk &	log
!	epsinf_mpb_inv : inverse of epsinf_mpb
!	rho0_eps, beta_eps : parameter defining Gygi & Fattebert dielectric function
!	rhomin_mpb,rhomax_mpb,log_rhomin_mpb,log_rhomax_mpb : parameter defining Andreussi dielectric function
!	rhomin_mpb2,rhomax_mpb2,log_rhomin_mpb2,log_rhomax_mpb2 : define 2nd dielectric function if needed for any purpose
!	/DIELECTRIC DECREMENT/
!   dielectric_decrement_method : "Newton": integrate dd iterations into Newton solver
!                                 "SCF"   : perform an additional self-consistent iteration above the Newton solver to converge changed dielectric function
!	-> experimental, NOT yet to be used. if interested contact authors
!
!  //SPECIFIC CALCULATIONS//
!	surface_and_volume_calc : if true evaluate the quantum surface and volume of the solvation cavity needed to evaluate Gnonmf
!	dynamic_cavity : if false use free electron density to evaluate dielectric function epsilon. if true, epsilon is updated in each SCF step by the current electron density
!	dynamic_ions : same as dynamic_cavity for the Exclusion function alpha
!	not_converge_rho_mpb : if true converge SCF in vacuum and run MPBE only once to evaluate non-self-consistent energies and potential
!	mpbe_no_ion : true if c_mpb = 0 (1d-30), so the ionic strength is zero
!	evaluate_sc_Gnonmf : self-consistent evaluation of G^non-mf, means we add the corresponding functional dervative of the Gnonmf energy term
!		with respect to the electron density to the KS Hamiltonian. there may be convergence issues related to using this flag (see above paper)
!	surface_mpb, volume_mpb : quantum surface and volume of solvation cavity
!	delta_Gnonmf : finite difference parameter for calculation of Gnonmf (see paper)
!	en_pot_Gnonmf : correction term in total energy due to addition of term to KS Hamiltonian
!	en_Gnonmf : Gnonmf free energy
!	en_eps_rho, en_alpha_rho : correction terms in total energy due to addition of operators to KS Hamiltonian due to eps[n_el] and alpha[n_el]
!
!  //TECHNICAL VARIABLES//
!	precondition_with_ions : enables the user to define a different kappa_debye screening parameter for the SPE as a preconditioner. this flag is
!		highly experimental though and changing this should not yield better results
!	debug_lpb_solver : writes out several important variables of the SPE_solver which can be useful to determine numerical problems
!	delta_rho_in_merm : propagate the difference SPE source term instead of the source term itself during the iterations. NOT WORKING, YET!
!	zero_cutoff_mpb : number of radial shells within which dielectric function is set to 1 explicitely. NOT recommended to use dielectric function which
!		gradient is not exactly zero outside the transition region since this gives worse convergence
!	MERM_in_SPE_solver : perform MERM inside the SPE_solver infrastructure without updating the electrostatic potential
!		in each iterative step on the full grid. gives much faster method for pure solvent calculations, a bit faster with ions.
!	eta_merm : linear mixing parameter of SPE source functions
!	tau_merm : accuracy of delta-potential during MERM
!	l_dim_SPE_public : order of multipole expansion (lmax+1)**2, with lmax maximum of all species
!	l_dim_species_SPE_public : order of multipole expansion (lmax+1)**2, tabulated for each species
!	/FAR FIELD/
!	limit_l_max_SPE_ff : maximum angular momentum for multipole expansion in the far field
!	multipole_radius_SPE_sq : radius**2 at which far field is switched on
!	changed_SPE_ff_cutoffs : true if the far field cutoff radius or lmax has been changed from default
!	limit_l_max_SPE : this is the maximum angular momentum for the multipole expansion to solve the SPE. if any of the species
!		lmax'es l_hartree(i_species) is lower than this, these will be taken for the expansion, if they are higher limit_l_max_SPE will be taken.
!		if limit_l_max_SPE = 0, automatically the l_hartree(i_species) will be taken
!       multipole_radius_free_SPE : this is the radius where splines are doctored in SPE_solver and also this controls where the far field starts, if not
!               explicitely defined by multipole_radius_SPE_sq. until now this is always equal to multipole_radius_free, but leaving this
!		variable for later flexibility
!
!  //VARIABLES ON THE FHI-AIMS GRID//
!	delta_v_MPB : delta v potential. vfree + delta_v_MPB = v (full electrostatic potential)
!	delta_v_MPB_gradient : gradient of delta_v_MPB
!	delta_v_MPB_laplacian : laplacian of delta_v_MPB
!	dielec_func_mpb, dielec_func_gradient : dielectric function & gradient
!	dielec_func2, dielec_func_grad2 : 2nd dielectric function for any purpose
!	deps_drho,d2eps_drho2 : 1st & 2nd derivative of eps with respect to nel
!	deps2_drho,d2eps2_drho2 : same for dielec_func2
!	alpha_func_mpb, alpha_func_gradient_mpb : alpha function & gradient
!	alpha_func_anion_mpb, alpha_func_anion_gradient_mpb : alpha function & gradient for anions
!	dnion_dalpha : derivative of nion with respect to alpha
!	dalpha_drho : derivative of alpha with respect to nel
!	v_free_gradient_public : gradient of free atom potential vfree
!	v_free_gradient_public_atom : v_free_gradient_public saved atom-wise
!	qSPE : source term of SPE = (q + pi4_inv * L_1 * delta v)*epsinf_mpb
!	qSPE_old : qSPE from previous iteration
!	qSPE_plus_delta_rho: qSPE + delta n_el, as needed for multipole correction
!
!  //FORCES//
!	mpb_forces_on : if MPB forces are to be calculated in current SCF step
!	forces_mpb_force_off : if MPB forces are explicitely turned off from the user (just for testing)
!	delta_v_MPB_gradient_atom : atom wise tabulated delta_v_MPB_gradient
!	Gnonmf_forces : forces coming from Gnonmf term
!	use_Gnonmf_forces : Gnonmf_forces are calculated at some point
!	Gnonmf_forces_on : Gnonmf_forces are calculated in this SCF step
!
!  USES
    use dimensions
    use constants, only: pi4, pi4_inv
    use localorb_io, only: use_unit
    use physics, only : partition_tab, free_hartree_superpos
    use runtime_choices, only: restart_read

    implicit none


! 	    //GENERAL CONVERGENCE//
    logical :: mpb_solver_started = .False.
    logical :: lpb_converged = .false.
    logical :: initialized_lpb_solver = .False.
    logical :: mpb_converged
    character*10 :: mpb_start_crit = 'converge'
    integer :: start_at_scf_step = 0

! 	    //SYSTEM SPECIFICATIONS//
    logical :: solve_pbe_only
    logical :: correct_dielectric_decrement
    logical :: solve_lpbe_only
    logical :: solve_debye_only
    logical :: use_mpbe_free_energy
    logical :: mpbe_no_ion
    logical :: use_separate_alpha
    logical :: solve_lpb_with_constant_dielec = .false.
    integer :: dieleckind = 1
    integer :: alphakind = 1
    integer :: z_mpb
    real*8 :: sqrt_z_mpb_div_kBT_mpb
    character*10 :: dielectric_decrement_method = 'Newton'
    real*8 :: freeen_LPBE_nonelstat_solv_energy
    real*8 :: d_alpha_ion, xi_alpha_ion, offset_kappa_mpb = 0d0
    real*8 :: a_mpb, c_mpb,beta_eps
    real*8 :: phi_zero_mpb, kBT_mpb, T_mpb
    real*8 :: alphaplusgamma_mpb, beta_mpb
    real*8 :: epsinf_mpb, epsinf_mpb_inv, log_epsinf_mpb
    real*8 :: rho0_eps,eta_eps,rhomin_mpb,rhomax_mpb,log_rhomin_mpb,log_rhomax_mpb
    real*8 :: kappa_div_sqrt_eps, kappainf_mpb, kappa_debye
    real*8 :: rhomin_mpb2 = 1d-4 !0.0001 !1d-5 !1d-6
    real*8 :: rhomax_mpb2 = 1d-1 !0.005 !7d-3 !1d-2
    real*8 :: log_rhomin_mpb2
    real*8 :: log_rhomax_mpb2
    real*8 :: d_alpha_ion_anion !if use_separate_alpha = True, we use two alpha functions, one
    !for the cations and one for the anion, these are the parameters for the anionic alpha function
    real*8 :: xi_alpha_ion_anion
    real*8 :: rhomin_kappa_anion
    real*8 :: rhomax_kappa_anion
    real*8 :: log_rhomin_kappa_anion,log_rhomax_kappa_anion
    real*8 :: rhomin_kappa,rhomax_kappa,log_rhomin_kappa,log_rhomax_kappa

    !DIELECTRIC DECREMENT
    integer :: decrement_kind = 3
    real*8 :: rhofrac = 0.1
    !if decrement_kind = 0: take a linear curve obtained from Barthels experimental data: eps = 78.36 +decrement_slope * c_ion
    !(J. Barthel, R. Buchner, and M. Muensterer,DECHEMA Chemistry Data Series, Vol. 12, 1995)
    !if decrement_kind = 1: take the curve from Buchner 1999: eps = 78.45 - 15.45 * c + 3.76 * c^3/2 (eps = 78.45 + decr_a * c + decr_b * c^(3/2))
    !(R. buchner, G.T. Hefter and P.M.May, J.Phys.Chem. A 103, 1-9 (1999))
    !if decrement_kind = 2: take the curve from Buchner up to c = 4.5M (highest verified experimental value). after that we suggest a monotonic exponential
    !decrease of the dielectric permittivity of the form exp_a*exp(-exp_b*c) + exp_c. exp_c we recommend to set to the value of the Buchner curve at 6.1M (saturated NaCl) solution
    !so for infinite concentrations the curve converges to eps_buchner(6.1M). the curve is smoothly attached to the Buchner curve (same derivative and function value at c=4.5M)
    !if decrement_kind = 3: take linear curve but with different values for cations and anions as suggested by Nakayama and Andelman, JCP 142, 044706 (2015).
    !eps = 78.3 - gamma_AD_plus * c_+ - gamma_AD_min * c_-
    real*8 :: decr_a = -15.45
    real*8 :: decr_b = 3.76
    real*8 :: gamma_AD_plus = 8 !for Na+
    real*8 :: gamma_AD_min = 3 !for Cl-
    real*8 :: decrement_slope = -11.5061 !how does eps change with concentration of ions, linear behaviour assumed, this is slope...
    real*8 :: start_dielec_decr_corr_acc = 1d-5 !choose accuracy in potential that must be reached before starting add dielectric decrement effects
    real*8 :: epsmin = 2d0 !minimum epsilon for very high ionic concentrations
    real*8 :: eps_cut_decr = 4.5
    real*8 :: exp_b, exp_a
    real*8 :: exp_c = 20 !40.
    real*8, dimension(:,:), allocatable :: ion_conc !first index is type of ions, 1=positive ions, 2=negative ions
    real*8, dimension(:,:), allocatable :: ion_conc_deriv
    real*8, dimension(:,:), allocatable :: deps_dcion !partial derivative of dielectric function with respect to ion concentration
    real*8, dimension(:,:), allocatable :: dcion_dalpha !partial derivative of ionic concentrations with respect to ion exclusion function
    real*8, dimension(:,:,:), allocatable :: dcion_dgradv !partial derivative of ionic concentrations with respect to gradient of v

! 	    //SPECIFIC CALCULATIONS//
    logical :: surface_and_volume_calc = .True.
    logical :: dynamic_cavity = .false.
    logical :: dynamic_ions = .false.
    logical :: not_converge_rho_mpb = .false.
    logical :: evaluate_sc_Gnonmf = .true.
    real*8 :: surface_mpb = 0d0
    real*8 :: volume_mpb = 0d0
    real*8 :: delta_Gnonmf = 1d-8 !FD parameter for Gnonmf calculation
    real*8 :: en_pot_Gnonmf = 0d0 !Gnon-mf energy term to correct for KS corrected Hamiltonian
    real*8 :: en_Gnonmf = 0d0 !Gnon-mf free energy
    real*8 :: en_eps_rho = 0d0, en_alpha_rho = 0d0


! 	    //TECHNICAL VARIABLES//
    logical :: atomic_MERM = .false. !if MERM iterations should be performed for each atom separately
    logical :: precondition_with_ions = .true.
    logical :: debug_lpb_solver = .false.
    logical :: delta_rho_in_merm = .false.
    logical :: MERM_in_SPE_solver = .true.
    integer :: zero_cutoff_mpb = 3
    integer :: limit_l_max_SPE_ff
    real*8 :: eta_merm, tau_merm
    integer :: l_dim_SPE_public
    integer, dimension(:), allocatable :: l_dim_species_SPE_public
    real*8 :: multipole_radius_SPE_sq
    logical :: changed_SPE_ff_cutoffs = .False.
    integer :: limit_l_max_SPE
    real*8  :: kappa_SPE = 0.0 ! is defined in kerker_mpb
    real*8, dimension(:), allocatable :: multipole_radius_free_SPE



! 	    //VARIABLES ON THE FHI-AIMS GRID//
    real*8, dimension(:), allocatable :: delta_v_MPB ! elstat lpn potential from solving the Helmholtz equation
    real*8, dimension(:,:), allocatable :: delta_v_MPB_gradient ! lpb gradient of elstat potential from solving the Helmholtz equation
    real*8, dimension(:), allocatable :: delta_v_MPB_laplacian

    real*8, dimension(:), allocatable :: dielec_func_mpb
    real*8, dimension(:), allocatable :: alpha_func_mpb
    real*8, dimension(:), allocatable :: alpha_func_anion_mpb
    real*8, dimension(:,:), allocatable :: dielec_func_gradient_mpb
    real*8, dimension(:,:), allocatable :: alpha_func_gradient_mpb
    real*8, dimension(:,:), allocatable :: alpha_func_anion_gradient_mpb
    real*8, dimension(:), allocatable :: dielec_func2
    real*8, dimension(:,:), allocatable :: dielec_func_grad2
    real*8, dimension(:), allocatable :: deps2_drho,d2eps2_drho2
    real*8, dimension(:), allocatable :: deps_drho, d2eps_drho2

    real*8, dimension(:,:), allocatable :: v_free_gradient_public
    real*8, dimension(:,:,:), allocatable :: v_free_gradient_public_atom
    real*8, dimension(:), allocatable :: qSPE_plus_delta_rho
    real*8, dimension(:), allocatable :: qSPE_old
    real*8, dimension(:), allocatable :: qSPE

    real*8, dimension(:), allocatable :: dnion_dalpha,dalpha_drho, dalpha_anion_drho

! 	  //FORCES//
    logical :: mpb_forces_on = .false., forces_mpb_force_off = .false. !mpb forces are calculated
    logical :: use_Gnonmf_forces = .false., Gnonmf_forces_on = .false.
    ! 	  real*8, dimension(:,:), allocatable :: forces_hf_mpb
    real*8, dimension(:,:,:),   allocatable :: delta_v_MPB_gradient_atom
    real*8, dimension(:,:), allocatable :: Gnonmf_forces


    ! 	  real*8, dimension(:), allocatable :: delta_v_MPB_const
    ! 	  real*8, dimension(:,:), allocatable :: delta_v_MPB_const_gradient
    !           real*8, dimension(:), allocatable :: delta_v_MPB_vac
    ! 	  real*8, dimension(:,:), allocatable :: delta_v_MPB_vac_gradient

contains

    !****s* FHI-aims/allocate_lpbsolver
    !  NAME
    !  allocate_lpbsolver
    !  SYNOPSIS
    subroutine allocate_lpbsolver()
    !  PURPOSE
    !	allocate all variables of lpb_solver_utilities

        use runtime_choices
        use mpi_tasks
        implicit none

        integer:: info

        if (.not.allocated(multipole_radius_free_SPE)) then
            allocate(multipole_radius_free_SPE(n_species),stat=info)
            call check_allocation(info, 'multipole_radius_free_SPE                  ')
            multipole_radius_free_SPE(:) = 0d0
        end if

        if (.not.allocated(Gnonmf_forces)) then
            allocate(Gnonmf_forces(3,n_atoms),stat=info)
            call check_allocation(info, 'Gnonmf_forces                  ')
            Gnonmf_forces(:,:) = 0d0
        end if

        if (.not.allocated(v_free_gradient_public_atom)) then
            allocate(v_free_gradient_public_atom(n_atoms,3,n_full_points),stat=info)
            call check_allocation(info, 'v_free_gradient_public_atom                  ')
            v_free_gradient_public_atom(:,:,:) = 0d0
        end if

        if (.not.allocated(d2eps2_drho2)) then
            allocate(d2eps2_drho2(n_full_points),stat=info)
            call check_allocation(info, 'd2eps2_drho2                  ')
            d2eps2_drho2(:) = 0d0
        end if

        if (.not.allocated(deps2_drho)) then
            allocate(deps2_drho(n_full_points),stat=info)
            call check_allocation(info, 'deps2_drho                  ')
            deps2_drho(:) = 0d0
        end if

        !reg_method=='vacandconst'
        ! 	  if (.not.allocated(delta_v_MPB_vac)) then
        ! 	    allocate(delta_v_MPB_vac(n_full_points),stat=info)
        ! 	    call check_allocation(info, 'delta_v_MPB_vac                  ')
        ! 	    delta_v_MPB_vac(:) = 0d0
        ! 	  end if
        !
        ! 	  if (.not.allocated(delta_v_MPB_vac_gradient)) then
        ! 	    allocate(delta_v_MPB_vac_gradient(3,n_full_points),stat=info)
        ! 	    call check_allocation(info, 'delta_v_MPB_vac_gradient                  ')
        ! 	    delta_v_MPB_vac_gradient(:,:) = 0d0
        ! 	  end if

        ! 	      if (.not.allocated(delta_v_MPB_const_gradient)) then
        ! 		allocate(delta_v_MPB_const_gradient(3,n_full_points),stat=info)
        ! 		call check_allocation(info, 'delta_v_MPB_const_gradient                  ')
        ! 		delta_v_MPB_const_gradient(:,:) = 0d0
        ! 	      end if


        ! 	      if (.not.allocated(delta_v_MPB_const)) then
        ! 		allocate(delta_v_MPB_const(n_full_points),stat=info)
        ! 		call check_allocation(info, 'delta_v_MPB_const                  ')
        ! 		delta_v_MPB_const(:) = 0d0
        ! 	      end if

        if (.not.allocated(ion_conc)) then
            allocate(ion_conc(2,n_full_points),stat=info)
            call check_allocation(info, 'ion_conc                  ')
            ion_conc(:,:) = 0d0
        end if

        if (.not.allocated(deps_dcion)) then
            allocate(deps_dcion(2,n_full_points),stat=info)
            call check_allocation(info, 'deps_dcion                  ')
            deps_dcion(:,:) = 0d0
        end if

        if (.not.allocated(dcion_dalpha)) then
            allocate(dcion_dalpha(2,n_full_points),stat=info)
            call check_allocation(info, 'dcion_dalpha                  ')
            dcion_dalpha(:,:) = 0d0
        end if

        if (.not.allocated(dcion_dgradv)) then
            allocate(dcion_dgradv(2,3,n_full_points),stat=info)
            call check_allocation(info, 'dcion_dgradv                  ')
            dcion_dgradv(:,:,:) = 0d0
        end if

        if (.not.allocated(ion_conc_deriv)) then
            allocate(ion_conc_deriv(2,n_full_points),stat=info)
            call check_allocation(info, 'ion_conc_deriv                  ')
            ion_conc_deriv(:,:) = 0d0
        end if

        if (.not.allocated(qSPE)) then
            allocate(qSPE(n_full_points),stat=info)
            call check_allocation(info, 'qSPE                  ')
            qSPE(:) = 0d0
        end if

        if (.not.allocated(dielec_func2)) then
            allocate(dielec_func2(n_full_points),stat=info)
            call check_allocation(info, 'dielec_func2                  ')
            dielec_func2(:) = 0d0
        end if

        if (.not.allocated(d2eps_drho2)) then
            allocate(d2eps_drho2(n_full_points),stat=info)
            call check_allocation(info, 'd2eps_drho2                  ')
            d2eps_drho2(:) = 0d0
        end if

        if (.not.allocated(deps_drho)) then
            allocate(deps_drho(n_full_points),stat=info)
            call check_allocation(info, 'deps_drho                  ')
            deps_drho(:) = 0d0
        end if

        if (.not.allocated(dnion_dalpha)) then
            allocate(dnion_dalpha(n_full_points),stat=info)
            call check_allocation(info, 'dnion_dalpha                  ')
            dnion_dalpha(:) = 0d0
        end if


        if (.not.allocated(dalpha_drho)) then
            allocate(dalpha_drho(n_full_points),stat=info)
            call check_allocation(info, 'dalpha_drho                  ')
            dalpha_drho(:) = 0d0
        end if

        if (.not.allocated(dalpha_anion_drho)) then
            allocate(dalpha_anion_drho(n_full_points),stat=info)
            call check_allocation(info, 'dalpha_anion_drho                  ')
            dalpha_anion_drho(:) = 0d0
        end if


        if (.not.allocated(dielec_func_grad2)) then
            allocate(dielec_func_grad2(3,n_full_points),stat=info)
            call check_allocation(info, 'dielec_func_grad2                  ')
            dielec_func_grad2(:,:) = 0d0
        end if

        if (.not.allocated(delta_v_MPB)) then
            allocate ( delta_v_MPB(n_full_points),stat=info )
            call check_allocation(info, 'delta_v_MPB                  ')
            delta_v_MPB(:) = 0.0d0
        end if
        if (.not.allocated(qSPE_plus_delta_rho)) then
            allocate ( qSPE_plus_delta_rho(n_full_points),stat=info )
            call check_allocation(info, 'qSPE_plus_delta_rho                  ')
            qSPE_plus_delta_rho(:) = 0.0d0
        end if
        if (.not.allocated(l_dim_species_SPE_public)) then
            allocate ( l_dim_species_SPE_public(n_species),stat=info )
            call check_allocation(info, 'l_dim_species_SPE_public                  ')
            l_dim_species_SPE_public(:) = 0.0d0
        end if
        if (.not.allocated(dielec_func_mpb)) then
            allocate (dielec_func_mpb(n_full_points),stat=info)
            call check_allocation(info, 'dielec_func                           ')
            dielec_func_mpb(:) = 0d0
        end if
        if (.not.allocated(alpha_func_mpb)) then
            allocate (alpha_func_mpb(n_full_points),stat=info)
            alpha_func_mpb(:) = 0d0
            call check_allocation(info, 'alpha_func_mpb                           ')
        end if
        if (.not.allocated(alpha_func_anion_mpb).and.use_separate_alpha) then
            allocate (alpha_func_anion_mpb(n_full_points),stat=info)
            alpha_func_anion_mpb(:) = 0d0
            call check_allocation(info, 'alpha_func_anion_mpb                           ')
        end if
        if (.not.allocated(qSPE_old)) then
            allocate ( qSPE_old(n_full_points),stat=info )
            call check_allocation(info, 'qSPE_old                  ')
            qSPE_old(:) = 0d0
        end if
        if (.not.allocated(delta_v_MPB_gradient)) then
            allocate ( delta_v_MPB_gradient(3,n_full_points),stat=info )
            call check_allocation(info, 'delta_v_MPB_gradient                  ')
            delta_v_MPB_gradient(:,:) = 0d0
        end if
        if (.not.allocated(delta_v_MPB_laplacian)) then
            allocate ( delta_v_MPB_laplacian(n_full_points),stat=info )
            call check_allocation(info, 'delta_v_MPB_laplacian                  ')
            delta_v_MPB_laplacian(:) = 0d0
        end if
        if (.not.allocated(v_free_gradient_public)) then
            allocate ( v_free_gradient_public(3,n_full_points),stat=info )
            call check_allocation(info, 'v_free_gradient_public                  ')
            v_free_gradient_public(:,:) = 0d0
        end if
        if (.not.allocated(dielec_func_gradient_mpb)) then
            allocate ( dielec_func_gradient_mpb(3,n_full_points),stat=info )
            call check_allocation(info, 'dielec_func_gradient_mpb                  ')
            dielec_func_gradient_mpb(:,:) = 0d0
        end if
        if (.not.allocated(alpha_func_gradient_mpb)) then
            allocate ( alpha_func_gradient_mpb(3,n_full_points),stat=info )
            call check_allocation(info, 'alpha_func_gradient_mpb                  ')
            alpha_func_gradient_mpb(:,:) = 0d0
        end if
        if (.not.allocated(alpha_func_anion_gradient_mpb)) then
            allocate ( alpha_func_anion_gradient_mpb(3,n_full_points),stat=info )
            call check_allocation(info, 'alpha_func_anion_gradient_mpb                  ')
            alpha_func_anion_gradient_mpb(:,:) = 0d0
        end if
        if (use_forces.and..not.forces_mpb_force_off) then
            ! 	      if (.not.allocated(forces_hf_mpb)) then
            ! 		allocate (forces_hf_mpb(n_atoms,3),stat=info)
            ! 		call check_allocation(info, 'forces_hf_mpb                           ')
            ! 		forces_hf_mpb(:,:) = 0d0
            ! 	      end if
            if (.not.allocated(delta_v_MPB_gradient_atom)) then
                allocate (delta_v_MPB_gradient_atom(n_atoms,3,n_full_points),stat=info)
                call check_allocation(info, 'delta_v_MPB_gradient_atom                           ')
                delta_v_MPB_gradient_atom(:,:,:) = 0d0
            end if
        end if

        ! ============================
        ! Initializations of variables
        ! ============================
        log_rhomin_mpb2 = log(rhomin_mpb2)
        log_rhomax_mpb2 = log(rhomax_mpb2)

    end subroutine allocate_lpbsolver

    !****s* FHI-aims/allocate_lpbsolver_dummy
    !  NAME
    !  allocate_lpbsolver
    !  SYNOPSIS
    subroutine allocate_lpbsolver_dummy()
    !  PURPOSE
    !	allocate all variables of lpb_solver_utilities, if they are not
    !	actually used, but we still need to do that because some of them
    !       are passed in routines

        use runtime_choices
        use mpi_tasks
        implicit none

        integer:: info

        !    	  if (.not.allocated(multipole_radius_free_SPE)) then
        !	    allocate(multipole_radius_free_SPE(n_species),stat=info)
        !	    call check_allocation(info, 'multipole_radius_free_SPE                  ')
        !	    multipole_radius_free_SPE(:) = 0d0
        !	  end if

        if (.not.allocated(Gnonmf_forces)) then
            allocate(Gnonmf_forces(3,n_atoms),stat=info)
            call check_allocation(info, 'Gnonmf_forces                  ')
            Gnonmf_forces(:,:) = 0d0
        end if

        !	  if (.not.allocated(v_free_gradient_public_atom)) then
        !	    allocate(v_free_gradient_public_atom(n_atoms,3,n_full_points),stat=info)
        !	    call check_allocation(info, 'v_free_gradient_public_atom                  ')
        !	    v_free_gradient_public_atom(:,:,:) = 0d0
        !	  end if
        !
        if (.not.allocated(d2eps2_drho2)) then
            allocate(d2eps2_drho2(1),stat=info)
            call check_allocation(info, 'd2eps2_drho2                  ')
            d2eps2_drho2(1) = 0d0
        end if
        !
        if (.not.allocated(deps2_drho)) then
            allocate(deps2_drho(1),stat=info)
            call check_allocation(info, 'deps2_drho                  ')
            deps2_drho(1) = 0d0
        end if
        !
        !	  !reg_method=='vacandconst'
        !! 	  if (.not.allocated(delta_v_MPB_vac)) then
        !! 	    allocate(delta_v_MPB_vac(n_full_points),stat=info)
        !! 	    call check_allocation(info, 'delta_v_MPB_vac                  ')
        !! 	    delta_v_MPB_vac(:) = 0d0
        !! 	  end if
        !!
        !! 	  if (.not.allocated(delta_v_MPB_vac_gradient)) then
        !! 	    allocate(delta_v_MPB_vac_gradient(3,n_full_points),stat=info)
        !! 	    call check_allocation(info, 'delta_v_MPB_vac_gradient                  ')
        !! 	    delta_v_MPB_vac_gradient(:,:) = 0d0
        !! 	  end if
        !
        !! 	      if (.not.allocated(delta_v_MPB_const_gradient)) then
        !! 		allocate(delta_v_MPB_const_gradient(3,n_full_points),stat=info)
        !! 		call check_allocation(info, 'delta_v_MPB_const_gradient                  ')
        !! 		delta_v_MPB_const_gradient(:,:) = 0d0
        !! 	      end if
        !
        !
        !! 	      if (.not.allocated(delta_v_MPB_const)) then
        !! 		allocate(delta_v_MPB_const(n_full_points),stat=info)
        !! 		call check_allocation(info, 'delta_v_MPB_const                  ')
        !! 		delta_v_MPB_const(:) = 0d0
        !! 	      end if
        !
        !	      if (.not.allocated(ion_conc)) then
        !		allocate(ion_conc(2,n_full_points),stat=info)
        !		call check_allocation(info, 'ion_conc                  ')
        !		ion_conc(:,:) = 0d0
        !	      end if
        !
        !	      if (.not.allocated(ion_conc_deriv)) then
        !		allocate(ion_conc_deriv(2,n_full_points),stat=info)
        !		call check_allocation(info, 'ion_conc_deriv                  ')
        !		ion_conc_deriv(:,:) = 0d0
        !	      end if
        !
        !	      if (.not.allocated(qSPE)) then
        !		allocate(qSPE(n_full_points),stat=info)
        !		call check_allocation(info, 'qSPE                  ')
        !		qSPE(:) = 0d0
        !	      end if
        !
        if (.not.allocated(dielec_func2)) then
            allocate(dielec_func2(1),stat=info)
            call check_allocation(info, 'dielec_func2                  ')
            dielec_func2(1) = 0d0
        end if
        !
        if (.not.allocated(d2eps_drho2)) then
            allocate(d2eps_drho2(1),stat=info)
            call check_allocation(info, 'd2eps_drho2                  ')
            d2eps_drho2(1) = 0d0
        end if
        !
        !	      if (.not.allocated(deps_drho)) then
        !		allocate(deps_drho(n_full_points),stat=info)
        !		call check_allocation(info, 'deps_drho                  ')
        !		deps_drho(:) = 0d0
        !	      end if
        !
        !	      if (.not.allocated(dnion_dalpha)) then
        !		allocate(dnion_dalpha(n_full_points),stat=info)
        !		call check_allocation(info, 'dnion_dalpha                  ')
        !		dnion_dalpha(:) = 0d0
        !	      end if
        !
        !
        !	      if (.not.allocated(dalpha_drho)) then
        !		allocate(dalpha_drho(n_full_points),stat=info)
        !		call check_allocation(info, 'dalpha_drho                  ')
        !		dalpha_drho(:) = 0d0
        !	      end if
        !
        if (.not.allocated(dielec_func_grad2)) then
            allocate(dielec_func_grad2(1,1),stat=info)
            call check_allocation(info, 'dielec_func_grad2                  ')
            dielec_func_grad2(1,1) = 0d0
        end if
        !
        !	      if (.not.allocated(delta_v_MPB)) then
        !		allocate ( delta_v_MPB(n_full_points),stat=info )
        !		call check_allocation(info, 'delta_v_MPB                  ')
        !		delta_v_MPB(:) = 0.0d0
        !	      end if
        !	      if (.not.allocated(qSPE_plus_delta_rho)) then
        !		allocate ( qSPE_plus_delta_rho(n_full_points),stat=info )
        !		call check_allocation(info, 'qSPE_plus_delta_rho                  ')
        !		qSPE_plus_delta_rho(:) = 0.0d0
        !	      end if
        !	      if (.not.allocated(l_dim_species_SPE_public)) then
        !		allocate ( l_dim_species_SPE_public(n_species),stat=info )
        !		call check_allocation(info, 'l_dim_species_SPE_public                  ')
        !		l_dim_species_SPE_public(:) = 0.0d0
        !	      end if
        if (.not.allocated(dielec_func_mpb)) then
            allocate (dielec_func_mpb(1),stat=info)
            call check_allocation(info, 'dielec_func                           ')
            dielec_func_mpb(1) = 0d0
        end if
        !	      if (.not.allocated(alpha_func_mpb)) then
        !	      allocate (alpha_func_mpb(n_full_points),stat=info)
        !	      alpha_func_mpb(:) = 0d0
        !	      call check_allocation(info, 'alpha_func_mpb                           ')
        !	    end if
        !	      if (.not.allocated(alpha_func_anion_mpb).and.use_separate_alpha) then
        !		allocate (alpha_func_anion_mpb(n_full_points),stat=info)
        !		alpha_func_anion_mpb(:) = 0d0
        !		call check_allocation(info, 'alpha_func_anion_mpb                           ')
        !	      end if
        !	      if (.not.allocated(qSPE_old)) then
        !		allocate ( qSPE_old(n_full_points),stat=info )
        !		call check_allocation(info, 'qSPE_old                  ')
        !		qSPE_old(:) = 0d0
        !	      end if
        !	      if (.not.allocated(delta_v_MPB_gradient)) then
        !		allocate ( delta_v_MPB_gradient(3,n_full_points),stat=info )
        !		call check_allocation(info, 'delta_v_MPB_gradient                  ')
        !		delta_v_MPB_gradient(:,:) = 0d0
        !	      end if
        !	      if (.not.allocated(delta_v_MPB_laplacian)) then
        !		allocate ( delta_v_MPB_laplacian(n_full_points),stat=info )
        !		call check_allocation(info, 'delta_v_MPB_laplacian                  ')
        !		delta_v_MPB_laplacian(:) = 0d0
        !	      end if
        if (.not.allocated(v_free_gradient_public)) then
            allocate ( v_free_gradient_public(1,1),stat=info )
            call check_allocation(info, 'v_free_gradient_public                  ')
            v_free_gradient_public(1,1) = 0d0
        end if
        if (.not.allocated(dielec_func_gradient_mpb)) then
            allocate ( dielec_func_gradient_mpb(1,1),stat=info )
            call check_allocation(info, 'dielec_func_gradient_mpb                  ')
            dielec_func_gradient_mpb(1,1) = 0d0
        end if
        !	      if (.not.allocated(alpha_func_gradient_mpb)) then
        !		allocate ( alpha_func_gradient_mpb(3,n_full_points),stat=info )
        !		call check_allocation(info, 'alpha_func_gradient_mpb                  ')
        !		alpha_func_gradient_mpb(:,:) = 0d0
        !	      end if
        !	      if (.not.allocated(alpha_func_anion_gradient_mpb)) then
        !		allocate ( alpha_func_anion_gradient_mpb(3,n_full_points),stat=info )
        !		call check_allocation(info, 'alpha_func_anion_gradient_mpb                  ')
        !		alpha_func_anion_gradient_mpb(:,:) = 0d0
        !	      end if
        !	    if (use_forces.and..not.forces_mpb_force_off) then
        !! 	      if (.not.allocated(forces_hf_mpb)) then
        !! 		allocate (forces_hf_mpb(n_atoms,3),stat=info)
        !! 		call check_allocation(info, 'forces_hf_mpb                           ')
        !! 		forces_hf_mpb(:,:) = 0d0
        !! 	      end if
        !	      if (.not.allocated(delta_v_MPB_gradient_atom)) then
        !		allocate (delta_v_MPB_gradient_atom(n_atoms,3,n_full_points),stat=info)
        !		call check_allocation(info, 'delta_v_MPB_gradient_atom                           ')
        !		delta_v_MPB_gradient_atom(:,:,:) = 0d0
        !	      end if
        !	    end if

    end subroutine allocate_lpbsolver_dummy

    !****s* FHI-aims/deallocate_lpbsolver
    !  NAME
    !  deallocate_lpbsolver
    !  SYNOPSIS
    subroutine deallocate_lpbsolver()
    !  PURPOSE
    !	deallocate all variables of lpb_solver_utilities

        use runtime_choices

        implicit none

        integer:: info

        if (allocated(multipole_radius_free_SPE)) then
            deallocate(multipole_radius_free_SPE)
        end if

        if (allocated(Gnonmf_forces)) then
            deallocate(Gnonmf_forces)
        end if

        if (allocated(deps2_drho)) then
            deallocate(deps2_drho)
        end if

        if (allocated(d2eps2_drho2)) then
            deallocate(d2eps2_drho2)
        end if

        if (allocated(qSPE)) then
            deallocate(qSPE)
        end if

        if (allocated(dielec_func2)) then
            deallocate(dielec_func2)
        end if

        if (allocated(deps_drho)) then
            deallocate(deps_drho)
        end if

        if (allocated(d2eps_drho2)) then
            deallocate(d2eps_drho2)
        end if

        if (allocated(dielec_func_grad2)) then
            deallocate(dielec_func_grad2)
        end if

        if (allocated(v_free_gradient_public_atom)) then
            deallocate(v_free_gradient_public_atom)
        end if

        if (allocated(ion_conc)) then
            deallocate(ion_conc)
        end if

        if (allocated(deps_dcion)) then
            deallocate(deps_dcion)
        end if

        if (allocated(dcion_dalpha)) then
            deallocate(dcion_dalpha)
        end if

        if (allocated(dcion_dgradv)) then
            deallocate(dcion_dgradv)
        end if

        if (allocated(ion_conc_deriv)) then
            deallocate(ion_conc_deriv)
        end if

        if (allocated(delta_v_MPB)) then
            deallocate(delta_v_MPB)
        end if

        if (allocated(qSPE_plus_delta_rho)) then
            deallocate(qSPE_plus_delta_rho)
        end if

        if (allocated(dielec_func_mpb)) then
            deallocate(dielec_func_mpb)
        end if

        if (allocated(dielec_func_gradient_mpb)) then
            deallocate(dielec_func_gradient_mpb)
        end if

        if (allocated(alpha_func_mpb)) then
            deallocate(alpha_func_mpb)
        end if

        if (allocated(alpha_func_anion_mpb)) then
            deallocate(alpha_func_anion_mpb)
        end if

        if (allocated(alpha_func_gradient_mpb)) then
            deallocate(alpha_func_gradient_mpb)
        end if
        if (allocated(alpha_func_anion_gradient_mpb)) then
            deallocate(alpha_func_anion_gradient_mpb)
        end if

        if (allocated(l_dim_species_SPE_public)) then
            deallocate(l_dim_species_SPE_public)
        end if

        if (allocated(delta_v_MPB_gradient)) then
            deallocate(delta_v_MPB_gradient)
        end if

        if (allocated(delta_v_MPB_laplacian)) then
            deallocate(delta_v_MPB_laplacian)
        end if


        if (allocated(qSPE_old)) then
            deallocate(qSPE_old)
        end if

        if (allocated(v_free_gradient_public)) then
            deallocate(v_free_gradient_public)
        end if

        ! 	  if (allocated(forces_hf_mpb)) then
        ! 	      deallocate(forces_hf_mpb)
        ! 	  end if
        if (allocated(delta_v_MPB_gradient_atom)) then
            deallocate(delta_v_MPB_gradient_atom)
        end if

        if (allocated(dnion_dalpha)) then
            deallocate(dnion_dalpha)
        end if

        if (allocated(dalpha_drho)) then
            deallocate(dalpha_drho)
        end if

        if (allocated(dalpha_anion_drho)) then
            deallocate(dalpha_anion_drho)
        end if

    end subroutine

    !****s* FHI-aims/prepare_mpb_solver
    !  NAME
    !  prepare_mpb_solver
    !  SYNOPSIS
    subroutine prepare_mpb_solver &
            (elec_converged,number_of_loops,rho,rho_gradient)
    !  PURPOSE
    ! This routine reads in the current convergence data from the scf_solver and controls the start
    ! of the MPBE solver.

        implicit none

        logical, intent(in) :: elec_converged
        integer, intent(in) :: number_of_loops
        real*8, dimension(n_spin,n_full_points), intent(in) :: rho
        real*8, dimension(3,n_spin,n_full_points), intent(in) :: rho_gradient

        !the procedure is as follows:
        !1st case: mpb_start_crit.eq.'step'
        !at x'th SCF step:
        !- mpb_solver_started = True
        !- run MPBE solver
        !- check convergence (elec_converged, converged). hopefully not reached, yet
        !- mpb_converged = True
        !at x+n'th SCF step: (mpb_converged = mpb_solver_started = True)
        !- check if elec_converged, if yes evaluate forces
        !2nd case: mpb_start_crit.eq.'converge'
        !at the SCF step where in the end convergence is reached (elec_converged = True)
        !- in the end of the step we have elec_converged = True,
        !  BUT converged = False, since still mpb_converged = False
        !at the following step
        !- mpb_solver_started = True
        !- run MPBE solver
        !- set elec_converged = False, since we will restart the SCF convergence
        !- set mpb_converged = True
        !following steps: like x+n'th SCF step in 1st case

        if (.not.mpb_solver_started) then
            !if the MPB solver has not been started, yet, check if it should be started now
            if ((mpb_start_crit.eq.'step'.and.start_at_scf_step.eq.number_of_loops)&
                    .or.(mpb_start_crit.eq.'converge'.and.elec_converged)) then !.and..not.&
                    !cube_dielec_func_outputted) then
                !start crit reached, in this scf step, save all necessary arrays (like hartree
                !potential and hartree potential gradient) and in the one after that start solver
                mpb_solver_started = .true.
            end if
            if (restart_read) then
                !start MPBE solver directly already in initialization
                mpb_solver_started = .true.
            end if
        end if

    end subroutine prepare_mpb_solver

    !****s* FHI-aims/postprocess_mpb_solver
    !  NAME
    !  postprocess_mpb_solver
    !  SYNOPSIS
    subroutine postprocess_mpb_solver &
            (elec_converged,converged)
    !  PURPOSE
    !  This routine resets the "converged" boolean if the SCF solver should be continued and if
    !  the MPB solver has not been started, yet.
        use mpi_utilities


        implicit none

        logical, intent(inout) :: elec_converged,converged

        if (elec_converged.and.mpb_start_crit.eq.'converge'.and.&
                .not.mpb_converged.and.mpb_solver_started.and..not.not_converge_rho_mpb) then
            !if our SCF converged in vacuum and we want to start the MPB solver in the next
            !step reset the elec_converged flag
            elec_converged = .false.
        end if


        if (mpb_solver_started) then
            !MPB solver started in the last SCF stop.
            !we can therefore let the SCF cycle stop as soon as the normal
            !convergence criteria are reached
            !if elec_converged, now also forces will be calculated
            mpb_converged = .true.
        end if

        if (elec_converged.and.not_converge_rho_mpb.and.mpb_solver_started) then
            !we have already calculated the potential, we do not need one more step
            converged = .true.
            elec_converged = .true.
        end if

        !if restart_read_only file selected
        !call run_mpb_solver(rho,rho_gradient,dielec_func_mpb)

    end subroutine postprocess_mpb_solver

    !****s* lpb_solver_utilities/evaluate_Gnonmf
    !  NAME
    !  evaluate_Gnonmf
    !  SYNOPSIS
    subroutine evaluate_Gnonmf(rho, rho_gradient, &
            local_Gnonmf_derivs, Gnonmf_gradient_deriv, &
            i_full_points)
    !  PURPOSE
    ! Evaluates the derivatives needed to calculate the functional derivative of
    ! Gnonmf with respect to the electron density

        use runtime_choices
        use mpi_utilities
        real*8 rho(n_spin)
        real*8 rho_gradient(3,n_spin)
        real*8 local_Gnonmf_derivs !derivatives of Gnonmf x [phi_i phi_j]
        real*8 Gnonmf_gradient_deriv(3) !derivatives of Gnonmf x nabla [phi_i phi_j]
        integer, intent(in) :: i_full_points

        real*8 :: dtheta_deps, deps_drho_min, deps_drho_max, log_rho, part1a, part1, total_rho,&
            deps_drho_cur
        integer :: i_coords
        real*8 :: rho_spinless
        real*8 :: rho_gradient_spinless(3)
        real*8 :: rho_gradient_abs
        integer :: info
        rho_spinless = sum(rho)
        rho_gradient_spinless = sum(rho_gradient,DIM=2)
        rho_gradient_abs = sqrt(sum(rho_gradient_spinless**2))

        dtheta_deps = - 1d0/(epsinf_mpb-1d0)

        part1a = log_rhomax_mpb - log_rhomin_mpb

        !evaluate deps_drho at rho - delta
        total_rho = rho_spinless-delta_Gnonmf/2d0
        if (total_rho.gt.rhomin_mpb.and.total_rho.lt.rhomax_mpb) then
            log_rho = log(total_rho)
            part1 = -log_epsinf_mpb/part1a*(1-cos(pi4*0.5d0*(log_rhomax_mpb-log_rho)/part1a))
            deps_drho_min = exp(t_dielec(log_rho,rhomin_mpb,rhomax_mpb,epsinf_mpb))*&
                part1*1.0d0/total_rho
        else
            deps_drho_min = 0d0
        end if
        !evaluate deps_drho at rho + delta
        total_rho = rho_spinless+delta_Gnonmf/2d0
        if (total_rho.gt.rhomin_mpb.and.total_rho.lt.rhomax_mpb) then
            log_rho = log(total_rho)
            part1 = -log_epsinf_mpb/part1a*(1-cos(pi4*0.5d0*(log_rhomax_mpb-log_rho)/part1a))
            deps_drho_max = exp(t_dielec(log_rho,rhomin_mpb,rhomax_mpb,epsinf_mpb))*&
                part1*1.0d0/total_rho
        else
            deps_drho_max = 0d0
        end if
        !evaluate deps_drho at rho
        total_rho = rho_spinless
        if (total_rho.gt.rhomin_mpb.and.total_rho.lt.rhomax_mpb) then
            log_rho = log(total_rho)
            part1 = -log_epsinf_mpb/part1a*(1-cos(pi4*0.5d0*(log_rhomax_mpb-log_rho)/part1a))
            deps_drho_cur = exp(t_dielec(log_rho,rhomin_mpb,rhomax_mpb,epsinf_mpb))*&
                part1*1.0d0/total_rho
        else
            deps_drho_cur = 0d0
        end if

        local_Gnonmf_derivs = &
            0.0000339882737736419 * beta_mpb * dtheta_deps * deps_drho_cur !volume term
        local_Gnonmf_derivs = & !surface term
        local_Gnonmf_derivs + &
            0.0000006423049116300181 * alphaplusgamma_mpb/delta_Gnonmf * &
            rho_gradient_abs*&
            (deps_drho_max - deps_drho_min) * dtheta_deps
        Gnonmf_gradient_deriv = 0d0
        do i_coords = 1, 3, 1
            if (abs(rho_gradient_spinless(i_coords)).gt.0d0) then
                Gnonmf_gradient_deriv(i_coords) = &
                    0.0000006423049116300181*alphaplusgamma_mpb/delta_Gnonmf * &
                    rho_gradient_spinless(i_coords) / rho_gradient_abs * &
                    (volume_switch(rho_spinless+delta_Gnonmf/2d0) - &
                    volume_switch(rho_spinless-delta_Gnonmf/2d0))
            end if
        end do

        !	    Gnonmf_gradient_deriv = Gnonmf_gradient_deriv *50.
        !	    local_Gnonmf_derivs = local_Gnonmf_derivs *50.

    end subroutine evaluate_Gnonmf

    !****s* FHI-aims/evaluate_Gnonmf_energy_shell
    !  NAME
    !  evaluate_Gnonmf_energy_shell
    !  SYNOPSIS
    subroutine evaluate_Gnonmf_energy_shell &
            ( n_points, partition, local_Gnonmf_derivs, &
            Gnonmf_gradient_deriv, local_rho, local_rho_gradient, &
            en_pot_Gnonmf &
            )
    !  PURPOSE
    ! Evaluate the energy resulting from the Kohn-Sham Hamiltonian correction due to electron density
    ! dependence of Gnonmf

        use dimensions
        use energy_density
        use runtime_choices

        implicit none

        !  imported variables

        ! input
        integer :: n_points
        real*8 partition(n_points)
        real*8, dimension(n_points) :: local_Gnonmf_derivs
        real*8, dimension(3,n_points) :: &
            Gnonmf_gradient_deriv
        real*8, dimension(n_spin,n_points) :: local_rho
        real*8, dimension(3,n_spin,n_points) :: local_rho_gradient
        ! output
        real*8 :: en_pot_Gnonmf
        !  local variables
        real*8, dimension(n_points) :: energy_term
        ! counters
        integer :: i_point, i_spin, i_coord
        ! functions
        real*8, external :: ddot

        !  begin work

        energy_term = 0.d0

        do i_point = 1, n_points, 1

            do i_spin = 1, n_spin, 1
                energy_term(i_point) = energy_term(i_point) + &
                    local_Gnonmf_derivs(i_point) * &
                    local_rho(i_spin,i_point)
            enddo

        enddo

        do i_point = 1, n_points, 1

            do i_spin = 1, n_spin, 1
                do i_coord = 1, 3, 1
                    energy_term(i_point) = energy_term(i_point) + &
                        Gnonmf_gradient_deriv(i_coord,i_point) * &
                        local_rho_gradient(i_coord,i_spin,i_point)
                enddo
            enddo

        enddo

        en_pot_Gnonmf = en_pot_Gnonmf +&
            ddot(n_points,energy_term,1,partition,1)

    end subroutine evaluate_Gnonmf_energy_shell

    !****s* lpb_solver_utilities/evaluate_Gnonmf_forces
    !  NAME
    !    evaluate_Gnonmf_forces
    !  SYNOPSIS

    subroutine evaluate_Gnonmf_forces &
            (i_atom, KS_orbital, nuclear_gradient, orb_grad_dot_rho_grad_Gnonmf, &
            orb_hess_dot_rho_grad_Gnonmf, &
            partition_tab, max_occ_number, n_points, &
            global_atom, Gnonmf_forces )

    !  PURPOSE
    !  Subroutine evaluates the part of the forces arising from the dependence of the Gnonmf
    !  term in the free energy on nabla n_el. the terms arising from the dependence on n_el
    !  are included in the Pulay forces.
    !
    !  USES

        use dimensions
        use basis
        use species_data
        use geometry
        use runtime_choices
        use mpi_utilities
        use synchronize_mpi
        implicit none

        !  ARGUMENTS

        integer :: max_occ_number
        integer, intent(in) :: n_points
        real*8, dimension(max_occ_number, n_points),               intent(in)    :: KS_orbital
        real*8, dimension(n_points),                               intent(in)    :: partition_tab
        real*8, dimension(n_states, n_max_batch_size, 3),          intent(in)    :: nuclear_gradient
        real*8, dimension(n_states ,n_points),                     intent(in)    :: orb_grad_dot_rho_grad_Gnonmf
        real*8, dimension(n_states, n_max_batch_size, 3),          intent(in)    :: orb_hess_dot_rho_grad_Gnonmf
        integer, dimension(n_atoms),                               intent(in)    :: global_atom
        real*8, dimension(3, n_atoms),                             intent(inout) :: Gnonmf_forces
        ! Contributions for meta-GGAs. AJL

    !  INPUTS
    !   o max_occ_number -- maximum number of states with non-zore occupation
    !   o n_points -- number of grid points in the batch
    !   o KS_orbital -- Kohn Sman orbitals
    !   o partition_tab -- values of partition functions
    !   o nuclear_gradient -- gradients respects atom nucleus position
    !   o orb_grad_dot_rho_grad_Gnonmf -- orbital gradient times density gradient
    !   o orb_hess_dot_rho_grad_Gnonmf -- orbital hessian matrix times density gradient
    !   o global_atom -- atom which forces is currently calculated
    !
    !  OUTPUT
    !   o Gnonmf_forces -- Gnonmf force component is added here.
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




        integer :: compute_force_atom

        real*8, dimension(n_points) :: Gnonmf_k_term

        ! functions
        real*8, external :: ddot

        ! counters

        integer :: i_coord
        integer :: i_atom
        integer :: i_point

        integer :: info

        ! begin work
        compute_force_atom = global_atom(i_atom)
        do i_coord = 1, 3, 1
            do i_point = 1, n_points, 1
                Gnonmf_k_term(i_point) = &
                    ddot(max_occ_number, nuclear_gradient(1,i_point,i_coord),1,&
                    orb_grad_dot_rho_grad_Gnonmf(1,i_point),1) + &
                    ddot(max_occ_number, KS_orbital(1,i_point),1,&
                    orb_hess_dot_rho_grad_Gnonmf(1,i_point,i_coord),1)
            end do

            ! NOTE that the gga forces still lack a factor of 4, which is added at the very end of update_density
            Gnonmf_forces(i_coord, compute_force_atom) = &
                Gnonmf_forces(i_coord, compute_force_atom) + &
                ddot(n_points, partition_tab,1,Gnonmf_k_term,1)
        end do

        ! end work

    end subroutine evaluate_Gnonmf_forces
    !---------------------------------------------------------------------
    !******

    !****s* FHI-aims/evaluate_cavity_surface_and_volume
    !  NAME
    !  evaluate_cavity_surface_and_volume
    !  SYNOPSIS
    subroutine evaluate_cavity_surface_and_volume(rho,rho_gradient)
    !  PURPOSE
    ! Evaluates the quantum surface and volume of the solvation cavity

        use grids
        use mpi_utilities
        use synchronize_mpi_basic, only: sync_real_number

        implicit none

        real*8, dimension(n_spin,n_full_points) :: rho
        real*8, dimension(3,n_spin,n_full_points) :: rho_gradient
        integer :: i_spin,i_full_points, i_coords,i_my_batch,i_index
        real*8 :: rho_gradient_spinless(3),rho_spinless, rho_gradient_abs
        !   real*8, dimension(3) :: grid_coord2

        surface_mpb = 0d0
        volume_mpb = 0d0

        !   open(100,file='surface.asc')

        i_full_points = 0

        do i_my_batch = 1, n_my_batches, 1
            do i_index = 1, batches(i_my_batch)%size, 1
                i_full_points = i_full_points + 1
                if (partition_tab(i_full_points).gt.0d0) then

                    rho_spinless = sum(rho(:,i_full_points))
                    rho_gradient_spinless = sum(rho_gradient(:,:,i_full_points),DIM=2)
                    rho_gradient_abs = sqrt(sum(rho_gradient_spinless**2))
                    !       grid_coord2(:) = batches(i_my_batch) % points(i_index) % coords(:)
                    ! 		rho_current = 0d0
                    ! 		rho_gradient_current = 0d0
                    ! 		do i_spin = 1, n_spin, 1
                    ! 		  rho_gradient_current_spin = 0d0
                    ! 		  rho_current = rho_current + rho(i_spin,i_full_points)
                    ! 		  do i_coords = 1, 3, 1
                    ! 		    rho_gradient_current_spin = &
                        ! 		      rho_gradient_current_spin +&
                        ! 		      rho_gradient(i_coords,i_spin,i_full_points)**2
                    ! 		  end do
                    ! 		  rho_gradient_current_spin = sqrt(rho_gradient_current_spin)
                    ! 		  rho_gradient_current = rho_gradient_current + rho_gradient_current_spin
                    ! 		end do

                    !       if (abs(rho_current-0.0001).lt.0.2d-4) then
                    ! 	write(100,'(3(F20.10,1X))') &
                        ! 	    (grid_coord2(i_coords), i_coords=1,3,1)
                    !       end if

                    surface_mpb = surface_mpb + &
                        (volume_switch(rho_spinless+delta_Gnonmf/2d0)-volume_switch(rho_spinless-delta_Gnonmf/2d0))*&
                        rho_gradient_abs/delta_Gnonmf * partition_tab(i_full_points)

                    volume_mpb = volume_mpb + &
                        volume_switch(rho_spinless) * partition_tab(i_full_points)
                end if
            end do
        end do

        call sync_real_number(surface_mpb)
        call sync_real_number(volume_mpb)

        !   close(100)

    end subroutine evaluate_cavity_surface_and_volume

    !****s* FHI-aims/volume_switch
    !  NAME
    !  volume_switch
    !  SYNOPSIS
    real*8 function volume_switch(rho)
    !  PURPOSE
    ! theta-switch function needed in the evaluation of quantum surface and volume
        implicit none
        real*8 :: rho
        real*8 :: rho_in(1)
        real*8 :: dielec(1)
        rho_in(1) = rho
        !calculate dielectric function for current rho
        call dielecfunc_from_density_points( rho_in, 1, dielec )
        !calculate switch function
        volume_switch = (epsinf_mpb-dielec(1))/(epsinf_mpb-1d0)
        return
    end function volume_switch

    !real*8 function t_dielec(x,rhomin_t,rhomax_t,inf_t)
    !
    !	implicit none
    !
    !	real*8 :: x, innerlogterm, rhomin_t, rhomax_t,inf_t
    !
    !	innerlogterm =  0.5d0 * pi4 *(log(rhomax_t)-x)/(log(rhomax_t)-log(rhomin_t))
    !	t_dielec = 2.0d0 * pi4_inv *log(inf_t)*(innerlogterm-sin(innerlogterm)) !
    !
    !	return
    !end function t_dielec

    real*8 function volume_sphere_switch(coord,center,radius)
    !  PURPOSE
    ! switch function needed in the evaluation of ion density included in volume interval
        implicit none
        real*8 :: coord(3), center(3)
        real*8 :: radius
        !          if (radius>radius+delta/2.) then
        !              volume_sphere_switch = 0.0
        !    elif (radius<radius-delta/2.) then
        !
        !           exp(t_dielec(log(radius),radius-delta/2.,radius+delta/2.,2.))
        !          dielecfuncfrom_density_points(radius,1,dielec)
        if (sqrt(sum((coord-center)**2))> radius) then
            volume_sphere_switch = 0.0
        else
            volume_sphere_switch = 1.0
        end if
        return
    end function volume_sphere_switch

    subroutine evaluate_rdf(delta_rdf,rdf_nradii,rdf_deltaradii,function_for_rdf,filename)

        use dimensions
        use physics, only: partition_tab,rho
        use grids
        use relaxation
        use constants, only: pi
        use synchronize_mpi_basic, only: sync_real_number
        use mpi_utilities, only: myid
        integer :: i_my_batch,i_index,i_full_points, i_rdf
        real*8 :: center_of_mass(3), relative_coords(3, n_atoms)
        integer :: rdf_nradii
        real*8, dimension(:), allocatable :: rdf_radii
        real*8 :: coords_point(3)
        real*8 :: delta_rdf
        real*8 :: rdf_current, rdf_radius
        real*8, dimension(:), allocatable :: rdf
        real*8 :: volume_shell
        real*8 :: total_number
        real*8, dimension(n_full_points) :: function_for_rdf
        integer :: npoints
        real*8 :: delta_points
        real*8 :: rdf_deltaradii
        character(LEN=*) :: filename
        !rdf file


        !first find COM of molecule
        call get_relative_coords(coords, center_of_mass, relative_coords)

        allocate(rdf_radii(rdf_nradii))
        allocate(rdf(rdf_nradii))
        rdf=0d0
        do i_rdf=1, rdf_nradii, 1
            rdf_radius=i_rdf*rdf_deltaradii
            rdf_radii(i_rdf)=rdf_radius !save radii if needed for later use
            rdf_current=0d0
            volume_shell = 0d0
            ! 	!create points of icosahedron
            ! ! 	call create_points_on_sphere(rdf_radius,center_of_mass,rdf_points, distance,npoints)
            ! 	if not allocated(rho_interp) allocate(rho_interp(npoints))
            ! 	!interpolate electron density at the rdf_points wanted
            ! 	call interpolate_density(rdf_points,rho_interp)
            ! 	!interpolate ionic charge density at the rdf_points
            !
            ! 	!integrate over all just created points
            ! 	call integrate_over_sphere(rdf_points,int_dielec,int_ion_dens)
            ! 	if allocated(rdf_points) deallocate(rdf_points)
            ! 	if allocated(rho_interp) deallocate(rho_interp)
            total_number=0d0
            i_full_points = 0
            do i_my_batch = 1, n_my_batches, 1
                do i_index = 1, batches(i_my_batch)%size, 1
                    i_full_points = i_full_points + 1
                    if (partition_tab(i_full_points).gt.0d0) then
                        coords_point = batches(i_my_batch) % points(i_index) % coords(:)
                        rdf_current = rdf_current + &
                            !volume_sphere_switch gives 1 if current point is inside sphere of radius
                            !rdf_radius+-delta_rdf
                            (volume_sphere_switch(coords_point,center_of_mass,rdf_radius+delta_rdf/2d0)-volume_sphere_switch(coords_point,center_of_mass,rdf_radius-delta_rdf/2d0))*&
                            function_for_rdf(i_full_points)*&
                            partition_tab(i_full_points)
                        !evaluate the volume of the shell
                        volume_shell = volume_shell + &
                            (volume_sphere_switch(coords_point,center_of_mass,rdf_radius+delta_rdf/2d0)-volume_sphere_switch(coords_point,center_of_mass,rdf_radius-delta_rdf/2d0))*&
                            partition_tab(i_full_points)
                        total_number = total_number + &
                            function_for_rdf(i_full_points)*&
                            partition_tab(i_full_points)
                    end if
                end do
            end do
            call sync_real_number(rdf_current)
            call sync_real_number(volume_shell)
            call sync_real_number(total_number)
            !volume_shell=((rdf_radius+delta_rdf)**3-(rdf_radius-delta_rdf/2.)**3)*4./3.*pi
            rdf(i_rdf)=rdf_current/volume_shell
        end do

        !call sync_vector(rdf,rdf_nradii)
        write(use_unit,*) 'totalpoints',rdf_nradii
        write(use_unit,*) 'allradii', rdf_radii
        if (myid==0) then
            open(100,file=filename)
            do i_rdf=1, rdf_nradii, 1
                if (myid==0)  write(100,'(3(F20.10,1X))') rdf_radii(i_rdf), rdf(i_rdf),total_number
                !        if (myid==0)  write(100,*) '%',center_of_mass(1),center_of_mass(2),center_of_mass(3)
            end do
            close(100)
        end if

        !call output_delta_v_step(sum(rho,DIM=1),1)
        !stop

    end subroutine evaluate_rdf
    !****s* FHI-aims/evaluate_g_gradient
    !  NAME
    !  evaluate_g_gradient
    !  SYNOPSIS
    subroutine evaluate_g_gradient &
            ( prec_l_dim,prec_l_dim_current, &
            one_over_dist_tab, dir_tab, trigonom_tab,  l_ylm_max,&
            ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab,  &
            radial_wave, radial_wave_deriv, radial_wave_2nd_deriv,&
            gradient,prec_l,laplacian, index_lm,evaluate_laplacian)

    !  PURPOSE
    ! Evaluates the gradient of delta_v_MPB, given as input by the radial multipole moments
    ! radial_wave and its derivative radial_wave_deriv.

        !the laplacian seems to not work at c>=0.1M? have to check that again...

        use constants, only: pi
        use dimensions
        use basis
        use geometry
        use grids
        use spline
        use pbc_lists
        implicit none

    !  ARGUMENTS

        logical :: evaluate_laplacian
        integer :: prec_l_dim, prec_l_dim_current, i_lm, l_ylm_max, prec_l, i_l, i_m
        integer :: index_lm(-l_ylm_max:l_ylm_max,0:l_ylm_max)
        real*8 one_over_dist_tab
        real*8 dir_tab(3)
        real*8 trigonom_tab(4)
        real*8 ylm_tab ((l_ylm_max+1)**2)
        real*8 dylm_dtheta_tab ((l_ylm_max+1)**2)
        real*8 scaled_dylm_dphi_tab ((l_ylm_max+1)**2)

        real*8 radial_wave(prec_l_dim)
        real*8 radial_div_dist(prec_l_dim)
        real*8 radial_wave_deriv(prec_l_dim)
        real*8 radial_wave_2nd_deriv(prec_l_dim)

        real*8 gradient(3)

        !  local variables

        real*8 dpsi_dr(prec_l_dim)
        real*8 metric_dpsi_dtheta(prec_l_dim)
        real*8 metric_dpsi_dphi(prec_l_dim)

        real*8 e_theta(3)

        real*8 laplacian, laplace_mp

    !  INPUTS
    !  o one_over_dist_tab -- 1/(distance to atoms)
    !  o dir_tab -- normalized direction to relevant atoms
    !  o trigonom_tab -- values of trigonometric functions
    !  o radial_wave -- radial basis functions
    !  o radial_wave_deriv -- derivative of radial basis functions
    !  o prec_l_dim -- number of basis functions

    !  OUTPUT
    !   o gradient -- gradient stores the 3d cartesian gradient of each basis function at the
    !                 present integration point
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

        !     begin work

        !     tabulate unit vector e_theta in local spherical coordinate system
        !     ( e_r is already given through dir_tab,
        !       e_phi is simply [-sin(phi),cos(phi),0] in trigonom_tab -
        !       no need to compute e_r, e_phi again.)

        !     cos(phi) * cos(theta)

        if (evaluate_laplacian) then
            laplacian = 0d0
        end if

        e_theta(1) =  &
            trigonom_tab(4) *  &
            trigonom_tab(2)

        !     sin(phi) * cos(theta)
        e_theta(2) =  &
            trigonom_tab(3) *  &
            trigonom_tab(2)

        !       - sin(theta)
        e_theta(3) =  &
            - trigonom_tab(1)

        radial_div_dist = radial_wave * one_over_dist_tab


        call scf_mul_vec_1 ( &
            dpsi_dr, prec_l_dim, &
            ylm_tab, &
            radial_wave_deriv &
            )

        call scf_mul_vec_1 ( &
            metric_dpsi_dtheta, prec_l_dim, &
            dylm_dtheta_tab, &
            radial_div_dist &
            )
        !scaled_dylm_dphi_tab = 1/sin(theta)*dy_(lm)/d(phi)
        call scf_mul_vec_1 ( &
            metric_dpsi_dphi, prec_l_dim, &
            scaled_dylm_dphi_tab, &
            radial_div_dist &
            )

        gradient(1:3) = 0.0d0


        if(sum(abs(dir_tab))< 1d-10)then
            !dm_lm/dr * Y_lm
            !m_l = 1
            gradient(1) = &
                -radial_wave_deriv(4)/2*sqrt(3.0/pi) !* 4*pi4_inv
            !m_l = -1
            gradient(2) = &
                radial_wave_deriv(2)/2*sqrt(3.0/pi) !* 4*pi4_inv
            !m_l = 0
            gradient(3) = &
                radial_wave_deriv(3)/2*sqrt(3.0/pi) !* 4*pi4_inv
            !
            ! 		do i_lm = 1, prec_l_dim_current, 1
            !
            ! 		    gradient(1) = gradient(1) +&
                ! 			e_theta(1) *  &
                ! 			metric_dpsi_dtheta(i_lm) - &
                ! 			trigonom_tab(3) *  &
                ! 			metric_dpsi_dphi(i_lm)
            !
            ! 		    gradient(2) = gradient(2) +&
                ! 			e_theta(2) *  &
                ! 			metric_dpsi_dtheta(i_lm) + &
                ! 			trigonom_tab(4) *  &
                ! 			metric_dpsi_dphi(i_lm)
            !
            ! 		    gradient(3) = gradient(3) +&
                ! 			e_theta(3) *  &
                ! 			metric_dpsi_dtheta(i_lm)
            !
            ! 		end do
            !
        else

            do i_lm = 1, prec_l_dim_current, 1

                gradient(1) = gradient(1) +&
                    dir_tab(1)*dpsi_dr(i_lm) + &
                    e_theta(1) *  &
                    metric_dpsi_dtheta(i_lm) - &
                    trigonom_tab(3) *  &
                    metric_dpsi_dphi(i_lm)

                gradient(2) = gradient(2) +&
                    dir_tab(2)*dpsi_dr(i_lm) + &
                    e_theta(2) *  &
                    metric_dpsi_dtheta(i_lm) + &
                    trigonom_tab(4) *  &
                    metric_dpsi_dphi(i_lm)

                gradient(3) = gradient(3) +&
                    dir_tab(3)*dpsi_dr(i_lm) + &
                    e_theta(3) *  &
                    metric_dpsi_dtheta(i_lm)

            end do
        end if

        if (evaluate_laplacian) then
            do i_l = 0, prec_l
                do i_m = -i_l, i_l
                    !trigonom_tab(1,2,3,4): sin(theta),cos(theta),sin(phi),cos(phi)

                    laplace_mp = 2d0*one_over_dist_tab*radial_wave_deriv(index_lm(i_m,i_l))
                    laplace_mp = laplace_mp + radial_wave_2nd_deriv(index_lm(i_m,i_l))

                    laplacian = laplacian +&
                        ylm_tab(index_lm(i_m,i_l))*(&
                        laplace_mp-one_over_dist_tab**2 *i_l*(i_l+1)*radial_wave(index_lm(i_m,i_l)))

                end do
            end do
        end if

    end subroutine evaluate_g_gradient

    !****s* FHI-aims/get_cut_f
    !  NAME
    !  get_cut_f
    !  SYNOPSIS
    subroutine get_cut_f(cutf)
    !  PURPOSE
    ! Evaluate radius at which approximations to the functions sinh/(a+bcosh) can be used

        !determine distance where sinh/(a+bcosh) - 1/b < 1e-15
        implicit none

        real*8 :: i_rho, func, cutf, func_start
        integer :: steps

        func = 0.0d0
        steps = 0

        if ((phi_zero_mpb.lt.1.0/3.0).or.(phi_zero_mpb.gt.1.0d0)) then
            !3 extrema of derivative, taking the inflection point at:
            func_start = ACOSH((phi_zero_mpb**2+2*phi_zero_mpb-1)/((phi_zero_mpb-1)*&
                phi_zero_mpb))/(z_mpb/(kBT_mpb))
            !as the starting point for the limes evaluation, go right from there until func = 1/phi_zero_mpb
            i_rho = func_start
            do while (abs(func - 1.0d0/phi_zero_mpb).gt.1e-30.or.steps.lt.1)
                steps = steps + 1
                i_rho = i_rho + 0.1 * func_start
                func = (SINH(z_mpb*i_rho/(kBT_mpb))/(1-phi_zero_mpb+phi_zero_mpb*&
                    COSH(z_mpb*i_rho/(kBT_mpb))))
            end do
        else
            !only minimum at i_rho = 0, estimate point on interval from 0 to 1/phi_zero_mpb
            func_start = 3.0d0*10**(-2.0d0*LOG10(z_mpb/(kBT_mpb))) * (z_mpb/(kBT_mpb))
            i_rho = func_start
            do while (abs(func - 1.0d0/phi_zero_mpb).gt.1e-30.or.steps.lt.1)
                steps = steps + 1
                i_rho = i_rho + 0.1 * func_start
                func = (SINH(z_mpb*i_rho/(kBT_mpb))/(1-phi_zero_mpb+phi_zero_mpb*&
                    COSH(z_mpb*i_rho/(kBT_mpb))))
            end do
        end if

        cutf = i_rho

    end subroutine get_cut_f

    !****s* FHI-aims/dielecfunc_from_density
    !  NAME
    !  dielecfunc_from_density
    !  SYNOPSIS
    subroutine dielecfunc_from_density( dens_array, dielec_func)
    !calculates a smooth dielectrc function at a number of grid points
    !  PURPOSE
    ! Calculates the dielectric function from the density nel on the FHI-aims grid

        use constants, only: bohr
        use grids

        implicit none

        integer :: i_atom, i_my_batch, i_index
        integer :: current_radial, i_full_points, i_spin

        real*8, dimension(n_spin,n_full_points), intent(in) :: dens_array
        real*8, dimension(n_full_points), intent(out) :: dielec_func

        real*8 :: total_rho
        real*8 :: total_rho_abs
        real*8 :: epsinf_mpb_current

        real*8 :: c_current, c_current_molarity
        real*8 :: test_func(n_full_points)

        real*8 :: pot_current

        real*8 :: au_to_molarity

        au_to_molarity = 1d0/(6.02214129D23 * 1.0D-27 * bohr**3)

        if (correct_dielectric_decrement) then
            deps_dcion = 0d0
        end if

        i_full_points = 0
        do i_my_batch = 1, n_my_batches, 1

            do i_index = 1, batches(i_my_batch)%size, 1


                ! grid point index in the final preconditioned residual
                i_full_points   = i_full_points + 1
                ! obtain distance^2 between integration point and current_center as well as a direction vector
                if (partition_tab(i_full_points).gt.0d0 .or. atomic_MERM) then
                    current_radial  = batches(i_my_batch) % points(i_index) % index_radial

                    total_rho = 0.0d0
                    do i_spin = 1, n_spin, 1
                        total_rho = total_rho + dens_array(i_spin,i_full_points)
                    end do
                    total_rho_abs = abs(total_rho)
                    ! 	      if (current_radial.lt.zero_cutoff_mpb) then
                    !set dielectric function to 1 explicitly near nuclei
                    ! 		dielec_func(i_full_points) = 1.0d0
                    ! 	      else
                    if (dieleckind.eq.0) then
                        ! use dielectric function from Fatterbert and Gygi, J. Comput. Chem. 2001

                        dielec_func(i_full_points) = 1 + (epsinf_mpb-1)/2 * (1 + (1-(total_rho_abs/rho0_eps)&
                            **(2*beta_eps))/(1+(total_rho_abs/rho0_eps)**(2*beta_eps)))

                    else if (dieleckind.eq.1) then

                        !use dielectric function from Andreussi, Dabo, Marzari, J.Chem.Phys 2012, Eq. (41),(42)
                        if (correct_dielectric_decrement) then
                            !include dielectric decrement effect
                            if (solve_lpbe_only) then
                                pot_current=(free_hartree_superpos(i_full_points)-delta_v_MPB(i_full_points))
                                c_current = kappainf_mpb**2 * alpha_func_mpb(i_full_points)*1d0/kBT_mpb*pot_current
                                c_current_molarity = c_current * au_to_molarity
                            else
                                !we devide the ionic concentrations by 2 because the experimental curve is plotted against the concentration of NaCl molecules!
                                c_current = (ion_conc(1,i_full_points) + ion_conc(2,i_full_points))/2d0 !*kappa_func_mpb(i_full_points)/kappainf_mpb
                                c_current_molarity = c_current * au_to_molarity
                            end if
                            if (decrement_kind == 0) then
                                epsinf_mpb_current = max(epsinf_mpb + decrement_slope * c_current_molarity,epsmin)
                                deps_dcion(:,i_full_points) = decrement_slope * au_to_molarity/2d0
                            else if (decrement_kind == 1.or.(decrement_kind == 2 .and.c_current_molarity.le.4.5)) then
                                epsinf_mpb_current = max(epsinf_mpb + decr_a * c_current_molarity + decr_b*c_current_molarity**(3d0/2d0),epsmin)
                                deps_dcion(:,i_full_points) = (decr_a * au_to_molarity + 3d0/2d0*decr_b*c_current**(1d0/2d0)*au_to_molarity**(3d0/2d0))/2d0
                            else if (decrement_kind == 2.and.c_current_molarity.gt.4.5) then
                                epsinf_mpb_current = exp_a*exp(-exp_b*c_current_molarity)+exp_c
                                deps_dcion(:,i_full_points) = -exp_a*exp_b*au_to_molarity*exp(-exp_b*c_current_molarity)/2d0
                            else if (decrement_kind == 3) then
                                epsinf_mpb_current = max(epsinf_mpb - gamma_AD_plus * ion_conc(1,i_full_points)*au_to_molarity - gamma_AD_min * ion_conc(2,i_full_points)*au_to_molarity,epsmin)
                                deps_dcion(1,i_full_points) = - gamma_AD_plus * au_to_molarity
                                deps_dcion(2,i_full_points) = - gamma_AD_min * au_to_molarity
                            end if
                            if (c_current_molarity.gt.6d0) then
                                write(use_unit,*) 'conc>6!', c_current_molarity
                                ! 				  stop
                            end if
                        else
                            epsinf_mpb_current = epsinf_mpb
                        end if

                        if (total_rho_abs.lt.rhomin_mpb) then
                            dielec_func(i_full_points) = epsinf_mpb_current
                        else if (total_rho_abs.gt.rhomax_mpb) then
                            dielec_func(i_full_points) = 1.0d0
                        else
                            dielec_func(i_full_points) = exp(t_dielec(log(total_rho_abs),rhomin_mpb,rhomax_mpb,epsinf_mpb_current))
                        end if

                        ! 		    if (abs(total_rho).lt.rhomin_mpb*rhofrac) then  0.007
                        !
                        ! 		      dielec_func2(i_full_points) = epsinf_mpb
                        ! 		    else if (abs(total_rho).gt.rhomin_mpb) then
                        ! 		      dielec_func2(i_full_points) = 1.0d0
                        ! 		    else
                        ! 		      dielec_func2(i_full_points) = exp(t_dielec(log(abs(total_rho)),rhomin_mpb*rhofrac,rhomin_mpb,epsinf_mpb))
                        ! 		    end if

                        if (total_rho_abs.lt.rhomin_mpb2) then
                            dielec_func2(i_full_points) = epsinf_mpb
                        else if (total_rho_abs.gt.rhomax_mpb2) then
                            dielec_func2(i_full_points) = 1.0d0
                        else
                            dielec_func2(i_full_points) = exp(t_dielec(log(total_rho_abs),rhomin_mpb2,rhomax_mpb2,epsinf_mpb))
                        end if

                    end if
                end if !partition_tab

            end do
        end do

        return

    end subroutine dielecfunc_from_density

    !****s* FHI-aims/t_dielec
    !  NAME
    !  t_dielec
    !  SYNOPSIS
    real*8 function t_dielec(x,rhomin_t,rhomax_t,inf_t)

        implicit none

        real*8 :: x, innerlogterm, rhomin_t, rhomax_t,inf_t

        innerlogterm =  0.5d0 * pi4 *(log(rhomax_t)-x)/(log(rhomax_t)-log(rhomin_t))
        t_dielec = 2.0d0 * pi4_inv *log(inf_t)*(innerlogterm-sin(innerlogterm)) !

        return
    end function t_dielec

    !****s* FHI-aims/dielecfunc_from_density_points
    !  NAME
    !  dielecfunc_from_density_points
    !  SYNOPSIS
    subroutine dielecfunc_from_density_points( dens_array, n_dens_points, dielec_func )
    !  PURPOSE
    ! Calculates the dielectric function from the electron density nel on any arbitrary grid

        !calculates a smooth dielectrc function at a number of grid points
        !if denspoints = n_full_points, it is assumed that we are on the aims integration grid and we can use abs_rho_spinless, which is defined on the aims grid!

        implicit none

        integer :: n_dens_points, i_dens_points

        real*8, dimension(n_dens_points), intent(in) :: dens_array
        real*8, dimension(n_dens_points), intent(out) :: dielec_func
        real*8 :: total_rho, total_rho_abs

        do i_dens_points = 1, n_dens_points, 1
            total_rho = dens_array(i_dens_points)
            if (n_dens_points.eq.n_full_points) then
                total_rho_abs = abs(total_rho) !abs_rho_spinless(i_dens_points)
            else
                total_rho_abs = abs(total_rho)
            end if
            if (dieleckind.eq.0) then
                ! use dielectric function from Fatterbert and Gygi, J. Comput. Chem. 2001

                dielec_func(i_dens_points) = 1 + (epsinf_mpb-1)/2 * (1 + (1-(abs(total_rho)/rho0_eps)&
                    **(2*beta_eps))/(1+(abs(total_rho)/rho0_eps)**(2*beta_eps)))

            else if (dieleckind.eq.1) then
                !use dielectric function from Andreussi, Dabo, Marzari, J.Chem.Phys 2012, Eq. (41),(42)

                if (total_rho_abs.lt.rhomin_mpb) then
                    dielec_func(i_dens_points) = epsinf_mpb
                else if (total_rho_abs.gt.rhomax_mpb) then
                    dielec_func(i_dens_points) = 1.0d0
                else
                    dielec_func(i_dens_points) = &
                        exp(t_dielec(log(abs(total_rho)),rhomin_mpb,rhomax_mpb,epsinf_mpb))
                end if


            end if
        end do


        return

    end subroutine dielecfunc_from_density_points

    !****s* FHI-aims/initialize_eps_s_func
    !  NAME
    !  initialize_eps_s_func
    !  SYNOPSIS
    subroutine initialize_eps_s_func()
    !  PURPOSE
    ! Needed for dielectric decrement. Initializes the eps->c^{s,bulk} curves taken from experiment
    !initializes the eps_s concentration dependency function parameter for decrement_kind==2
        real*8 :: buchner, buchnerderiv

        buchner = epsinf_mpb + decr_a * eps_cut_decr + decr_b*eps_cut_decr**(3d0/2d0)
        buchnerderiv = decr_a + decr_b * (3d0/2d0) * eps_cut_decr**(1d0/2d0)
        exp_b = -buchnerderiv/(buchner - exp_c)
        exp_a = -buchnerderiv/exp_b*exp(exp_b*eps_cut_decr)

    end subroutine initialize_eps_s_func

    !****s* FHI-aims/dielecfunc_gradient_from_density_gradient
    !  NAME
    !  dielecfunc_gradient_from_density_gradient
    !  SYNOPSIS
    subroutine dielecfunc_gradient_from_density_gradient(dielec_func, dielec_func_grad, rho,rho_grad)
    !  PURPOSE
    ! Evaluates the gradient of the dielectric function and all derivatives with respect to electron density
    ! on the FHI-aims grid
        use constants, only: bohr
        use grids
        use mpi_utilities

        implicit none

        integer ::  current_radial, i_index, i_my_batch
        real*8 :: numer_grad, denom_grad, innerterm, x, part1, part2,part1a
        real*8, dimension(n_spin,n_full_points), intent(in) :: rho
        real*8, dimension(n_full_points), intent(in) :: dielec_func
        real*8, dimension(3, n_spin, n_full_points), intent(in) :: rho_grad
        real*8, dimension(3, n_full_points), intent(out) :: dielec_func_grad
        integer :: i_coords, i_full_points
        real*8 :: t_dielec_deriv
        real*8, dimension(n_full_points) :: gradient_abs

        real*8 :: total_rho, total_rho_abs
        real*8, dimension(3) :: total_rho_grad
        integer :: i_spin
        real*8 :: deps_deps0

        real*8 :: c_current, c_current_molarity, epsinf_mpb_current,log_epsinf_mpb_current

        real*8 :: deps0_dv, dc_dalpha, dccat_dalpha, dcan_dalpha, dc_dv, dccat_dv,dcan_dv
        real*8 :: au_to_molarity, deps0_dalpha

        real*8, dimension(3) :: deps0_dgradv

        au_to_molarity = 1d0/(6.02214129D23 * 1.0D-27 * bohr**3)

        i_full_points = 0

        do i_my_batch = 1, n_my_batches, 1

            do i_index = 1, batches(i_my_batch)%size, 1


                i_full_points = i_full_points + 1
                current_radial  = batches(i_my_batch) % points(i_index) % index_radial

                if (partition_tab(i_full_points).gt.0d0 .or. atomic_MERM) then


                    total_rho = 0.0d0
                    total_rho_grad = 0.0d0
                    do i_spin = 1, n_spin, 1
                        total_rho = total_rho + rho(i_spin,i_full_points)
                        total_rho_grad = total_rho_grad + rho_grad(:,i_spin,i_full_points)
                    end do

                    total_rho_abs = abs(total_rho) !abs_rho_spinless(i_full_points)

                    ! 		if (current_radial.lt.zero_cutoff_mpb) then
                    ! 		  !set dielectric function to 1 explicitly near nuclei
                    ! 		  dielec_func_grad(:,i_full_points) = 0.0d0
                    ! 		  dielec_func_grad2(:,i_full_points) = 0d0
                    ! 		else
                    if (dieleckind.eq.0) then

                        part1 = 2.0d0*beta_eps*(total_rho_abs/rho0_eps)**(2*beta_eps-1.0d0)

                        part2 = (total_rho_abs/rho0_eps)**(2*beta_eps)

                        dielec_func_grad(:,i_full_points) = (epsinf_mpb-1.0d0)/2.0d0*(-(part1*(1-part2))/&
                            (rho0_eps*(1.0+part2)**2)-(part1)/(rho0_eps*(1+part2)))*total_rho_grad

                        ! 	    if (gradient_abs(i_full_points).ne.gradient_abs(i_full_points)) then
                        ! 	      write(use_unit,*) 'NaN found in dielecfunc_gradient_from_density_gradient', part1, part2
                        ! 	      stop
                        ! 	    end if
                    else if (dieleckind.eq.1) then

                        !calculate here the spatially dependent eps_s
                        if (correct_dielectric_decrement) then
                            if (solve_lpbe_only) then
                                !SR: should there not be a plus before the last term??
                                c_current = kappainf_mpb**2 * alpha_func_mpb(i_full_points)/kBT_mpb*(free_hartree_superpos(i_full_points)-delta_v_MPB(i_full_points))
                            else
                                c_current = (ion_conc(1,i_full_points) + ion_conc(2,i_full_points))/2d0 !*kappa_func_mpb(i_full_points)/kappainf_mpb
                            end if
                            c_current_molarity = c_current*au_to_molarity

                            if (decrement_kind == 0) then
                                epsinf_mpb_current = max(epsinf_mpb + decrement_slope * c_current_molarity,epsmin)
                            else if (decrement_kind == 1.or. (decrement_kind == 2.and.c_current_molarity.le.eps_cut_decr)) then
                                epsinf_mpb_current = max(epsinf_mpb + decr_a * c_current_molarity + decr_b*c_current_molarity**(3d0/2d0),epsmin)
                            else if (decrement_kind == 2.and.c_current_molarity.gt.eps_cut_decr) then
                                epsinf_mpb_current = exp_a*exp(-exp_b*c_current_molarity)+exp_c
                            else if (decrement_kind == 3) then
                                epsinf_mpb_current = max(epsinf_mpb - gamma_AD_plus * ion_conc(1,i_full_points)*au_to_molarity - gamma_AD_min * ion_conc(2,i_full_points)*au_to_molarity,epsmin)
                            end if
                        else
                            epsinf_mpb_current = epsinf_mpb
                        end if
                        log_epsinf_mpb_current = log(epsinf_mpb_current)
                        !tabulate deps0_dv
                        if (correct_dielectric_decrement.and.epsinf_mpb_current.gt.epsmin) then
                            !grad eps in the area where deps/dnel = 0
                            !grad eps = deps/deps0 gradeps0 = deps/deps0(=1) deps0/dv gradv (deps/dn = 0)
                            if (solve_lpbe_only) then
                                if (decrement_kind == 0) then
                                    deps0_dv = decrement_slope*au_to_molarity*&
                                        kappainf_mpb**2 * alpha_func_mpb(i_full_points)/kBT_mpb
                                else if (decrement_kind == 1.or. (decrement_kind == 2.and.c_current_molarity.le.eps_cut_decr)) then
                                    deps0_dv = (decr_a*au_to_molarity*&
                                        kappainf_mpb**2 * alpha_func_mpb(i_full_points)/kBT_mpb+&
                                        decr_b*(3d0/2d0)*c_current_molarity**(1d0/2d0))
                                else if (decrement_kind == 2.and.c_current_molarity.gt.eps_cut_decr) then
                                    deps0_dv = -exp_a*exp_b*au_to_molarity*&
                                        kappainf_mpb**2 * alpha_func_mpb(i_full_points)/kBT_mpb*exp(-exp_b*c_current_molarity)
                                end if
                            else
                                !we here use the functional derivative delta c_ion / delta v = dc_ion/dv - nabla dcion/d(nabla v)
                                !where delta means functional and "d" means partial derivatives
                                dccat_dv = au_to_molarity*ion_conc_deriv(1,i_full_points)
                                dcan_dv = au_to_molarity*ion_conc_deriv(2,i_full_points)
                                dc_dv = (dccat_dv + dcan_dv)/2d0
                                if (decrement_kind == 0) then
                                    !first part:
                                    deps0_dv = &
                                        !deps0/dc
                                        decrement_slope*&
                                        !dc/dv
                                        dc_dv
                                    deps0_dgradv = &
                                        decrement_slope*&
                                        (dcion_dgradv(1,:,i_full_points)+dcion_dgradv(2,:,i_full_points))/2d0*au_to_molarity
                                else if (decrement_kind == 1 .or. (decrement_kind == 2.and.c_current_molarity.le.eps_cut_decr)) then
                                    deps0_dv = &
                                        !deps/dc
                                        (decr_a+decr_b*(3d0/2d0)*c_current_molarity**(1d0/2d0))*&
                                        !dc/dv
                                        dc_dv
                                    deps0_dgradv = &
                                        (decr_a+decr_b*(3d0/2d0)*c_current_molarity**(1d0/2d0))*&
                                        (dcion_dgradv(1,:,i_full_points)+dcion_dgradv(2,:,i_full_points))/2d0*au_to_molarity
                                else if (decrement_kind == 2 .and.c_current_molarity.gt.eps_cut_decr) then
                                    deps0_dv = &
                                        !deps/dc
                                        (-exp_a*exp_b*exp(-exp_b*c_current_molarity))*&
                                        !dc/dv
                                        dc_dv
                                    deps0_dgradv = &
                                        (-exp_a*exp_b*exp(-exp_b*c_current_molarity))*&
                                        (dcion_dgradv(1,:,i_full_points)+dcion_dgradv(2,:,i_full_points))/2d0*au_to_molarity
                                else if (decrement_kind == 3) then
                                    deps0_dv = &
                                        !deps0/dccat*dccat/dv + deps0/dcan*dcan/dv
                                        -gamma_AD_plus*dccat_dv-gamma_AD_min*dcan_dv
                                    !deps0/dccat*dccat/dgradv + deps0/dcan*dcan/dgradv
                                    deps0_dgradv = &
                                        -gamma_AD_plus*dcion_dgradv(1,:,i_full_points)*au_to_molarity-&
                                        gamma_AD_min*dcion_dgradv(2,:,i_full_points)*au_to_molarity
                                end if
                            end if
                        else !(correct_dielectric_decrement.and.epsinf_mpb_current.gt.epsmin)
                            deps0_dv = 0d0
                        end if

                        if (total_rho_abs.lt.rhomin_mpb) then
                            !only eps0 depends on position (deps/dnel=0)...
                            !deps0_dv*grad v
                            dielec_func_grad(:,i_full_points) = deps0_dv* &
                                (v_free_gradient_public(:,i_full_points)+delta_v_MPB_gradient(:,i_full_points))

                            !                 !we need to also add deps0/d(grad v) * grad(grad v) = deps0/d(grad v) Delta v
                            !                 !Delta v = Delta delta_v - 4pi n_free
                            !                 !Delta delta_v = -4 pi * qSPE  + kappa^2 delta_v
                            !                 !add this as deps/
                            !                 dielec_func_grad(:,i_full_points) = dielec_func_grad(:,i_full_points) + &
                                !                     deps0_dgradv*(&
                                !                     !delta_v_MPB_laplacian(i_full_points)&
                                !                     -pi4*qSPE(i_full_points)/epsinf_mpb + kappa_debye**2 * delta_v_MPB(i_full_points)&
                                !                     -free_rho_superpos(i_full_points))

                            deps_drho(i_full_points) = 0d0
                            d2eps_drho2(i_full_points) = 0d0

                        else if (total_rho_abs.gt.rhomax_mpb) then
                            !inside the solvation cavity the gradient is always zero
                            dielec_func_grad(:,i_full_points) = 0.0d0
                            deps_drho(i_full_points) = 0d0
                            d2eps_drho2(i_full_points) = 0d0

                        else
                            x = log(total_rho_abs)

                            part1a = log_rhomax_mpb - log_rhomin_mpb
                            part1 = -log_epsinf_mpb_current/part1a*(1-cos(pi4*0.5d0*(log_rhomax_mpb-x)/part1a))
                            part2 = log_epsinf_mpb_current/part1a*sin(pi4*0.5d0*(log_rhomax_mpb-x)/part1a)*pi4*0.5d0/part1a/total_rho

                            !grad eps = deps/dn grad n
                            dielec_func_grad(:,i_full_points) = exp(t_dielec(x,rhomin_mpb,rhomax_mpb,epsinf_mpb_current))*&
                                part1*1.0d0/total_rho*total_rho_grad

                            deps_drho(i_full_points) = exp(t_dielec(x,rhomin_mpb,rhomax_mpb,epsinf_mpb_current))*&
                                part1*1.0d0/total_rho

                            d2eps_drho2(i_full_points) = deps_drho(i_full_points)/total_rho*part1 - &
                                exp(t_dielec(x,rhomin_mpb,rhomax_mpb,epsinf_mpb_current))/total_rho**2*part1+&
                                exp(t_dielec(x,rhomin_mpb,rhomax_mpb,epsinf_mpb_current))/total_rho*part2


                            if (correct_dielectric_decrement.and.epsinf_mpb_current.gt.epsmin) then
                                !we write: grad eps = deps/deps0 gradeps0 (to be done) + deps/drho graddrho (this we have already)
                                !deps/deps0  = t'*eps0^(t'-1); with t' = t/ln(eps0)
                                deps_deps0 = t_dielec(x,rhomin_mpb,rhomax_mpb,epsinf_mpb_current)/log_epsinf_mpb_current*&
                                    epsinf_mpb_current**(t_dielec(x,rhomin_mpb,rhomax_mpb,epsinf_mpb_current)/log_epsinf_mpb_current-1d0)
                                !grad eps += deps/deps0 * (deps0/dphi gradphi...)
                                !only if we are in the linear region of eps vs. conc
                                dielec_func_grad(:,i_full_points) = dielec_func_grad(:,i_full_points)+&
                                    deps_deps0 * deps0_dv*&
                                    (v_free_gradient_public(:,i_full_points)+delta_v_MPB_gradient(:,i_full_points))

                                if (alphakind==1) then
                                    ! 			      part1 = log(rhomax_mpb) - log(rhomin_mpb)
                                    ! 			      part1 = -log(epsinf_mpb)/part1*(1-cos(pi4*0.5d0*(log(rhomax_mpb)-x)/part1))
                                    !grad eps += deps/deps0 (...deps0/dalpha gradalpha); gradalpha = dalpha/drho gradrho
                                    !only if we are in the linear region of eps vs. conc
                                    if (solve_lpbe_only) then
                                        dc_dalpha = 2d0*kappainf_mpb*alpha_func_mpb(i_full_points)/kBT_mpb*&
                                            (free_hartree_superpos(i_full_points)+delta_v_MPB(i_full_points))*au_to_molarity
                                    else
                                        dc_dalpha = au_to_molarity*sum(dcion_dalpha(:,i_full_points))/2.
                                    end if

                                    !deps0_alpha = deps0_dc*dc_dalpha
                                    if (decrement_kind == 0) then
                                        deps0_dalpha = decrement_slope*dc_dalpha
                                    else if (decrement_kind == 1.or. (decrement_kind == 2.and.c_current_molarity.le.eps_cut_decr)) then
                                        deps0_dalpha = (decr_a+decr_b*(3d0/2d0)*c_current_molarity**(1d0/2d0))*dc_dalpha
                                    else if (decrement_kind == 2.and.c_current_molarity.gt.eps_cut_decr) then
                                        deps0_dalpha = -exp_a*exp_b*exp(-exp_b*c_current_molarity)*dc_dalpha
                                    else if (decrement_kind == 3) then
                                        deps0_dalpha = -gamma_AD_plus*dcion_dalpha(1,i_full_points)*au_to_molarity&
                                            -gamma_AD_min*dcion_dalpha(2,i_full_points)*au_to_molarity
                                    end if
                                    dielec_func_grad(:,i_full_points) = dielec_func_grad(:,i_full_points)+&
                                        deps_deps0* deps0_dalpha*alpha_func_gradient_mpb(:,i_full_points)

                                    ! 			    testfunc(i_full_points) = exp(t_dielec(x,rhomin_mpb,rhomax_mpb,epsinf_mpb))*&
                                        ! 			      part1*1.0d0/total_rho*(kappainf_mpb/(epsinf_mpb-1.0d0))
                                end if

                            end if !(correct_dielectric_decrement.and.epsinf_mpb_current.gt.epsmin)
                        end if

                        ! 		    if (c_current / (6.02214129D23 * 1.0D-27 * bohr**3).gt.4d0) then
                        ! 		      write(use_unit,*) 'Stop here since dielectric decrease may not be linear anymore'
                        ! 		      stop
                        ! 		    end if
                        ! 		    if (abs(total_rho).lt.rhomin_mpb*rhofrac.or.abs(total_rho).gt.rhomin_mpb) then
                        ! 		      dielec_func2_gradient(:,i_full_points) = 0.0d0
                        ! 		    else
                        ! 		      x = log(abs(total_rho))
                        !
                        ! 	! 	      part1 = log(rhomax_mpb) - log(rhomin_mpb)
                        ! 	! 	      dielec_func_grad(:,i_full_points) = &
                            ! 	! 		- 2.0d0*log(epsinf_mpb)*sin(pi*(log(rhomax_mpb)-x)/part1)**2*1.0d0/(abs(rho(i_full_points))*part1)*&
                            ! 	! 		rho_grad(:,i_full_points)
                        !
                        ! 		      part1 = log_rhomin_mpb - log(rhomin_mpb*rhofrac)
                        !
                        ! 		      part1 = -log_epsinf_mpb/part1*(1-cos(pi4*0.5d0*(log_rhomin_mpb-x)/part1))
                        !
                        ! 		      dielec_func2_gradient(:,i_full_points) = exp(t_dielec(x,rhomin_mpb*rhofrac,rhomin_mpb,epsinf_mpb))*&
                            ! 			part1*1.0d0/total_rho*total_rho_grad
                        !
                        ! 		    end if
                        gradient_abs(i_full_points) = sqrt(sum(dielec_func_grad(:,i_full_points)**2))

                        if (total_rho_abs.lt.rhomin_mpb2.or.total_rho_abs.gt.rhomax_mpb2) then
                            dielec_func_grad2(:,i_full_points) = 0.0d0
                            deps2_drho(i_full_points) = 0d0
                            d2eps2_drho2(i_full_points) = 0d0
                        else
                            x = log(total_rho_abs)

                            part1a = log_rhomax_mpb2 - log_rhomin_mpb2
                            part1 = -log_epsinf_mpb_current/part1a*(1-cos(pi4*0.5d0*(log_rhomax_mpb2-x)/part1a))
                            part2 = log_epsinf_mpb_current/part1a*sin(pi4*0.5d0*(log_rhomax_mpb2-x)/part1a)*pi4*0.5d0/part1a/total_rho

                            dielec_func_grad2(:,i_full_points) = exp(t_dielec(x,rhomin_mpb2,rhomax_mpb2,epsinf_mpb))*&
                                part1*1.0d0/total_rho*total_rho_grad

                            deps2_drho(i_full_points) = exp(t_dielec(x,rhomin_mpb2,rhomax_mpb2,epsinf_mpb_current))*&
                                part1*1.0d0/total_rho

                            d2eps2_drho2(i_full_points) = deps2_drho(i_full_points)/total_rho*part1 - &
                                exp(t_dielec(x,rhomin_mpb2,rhomax_mpb2,epsinf_mpb_current))/total_rho**2*part1+&
                                exp(t_dielec(x,rhomin_mpb2,rhomax_mpb2,epsinf_mpb_current))/total_rho*part2

                        end if

                    end if !(dieleckind.eq.0

                    ! 		end if !current_radial.lt.zero_cutoff_mpb
                end if !partition_tab
            end do
        end do


        return

    end subroutine dielecfunc_gradient_from_density_gradient

    !****s* FHI-aims/dielecfunc_gradient_from_density_gradient_points
    !  NAME
    !  dielecfunc_gradient_from_density_gradient_points
    !  SYNOPSIS
    subroutine dielecfunc_gradient_from_density_gradient_points(dielec_func, dielec_func_grad, rho,rho_grad, denspoints)
    !  PURPOSE
    !  same as dielecfunc_gradient_from_density_gradient but on any arbitrary grid

        !if denspoints = n_full_points, it is assumed that we are on the aims integration grid and we can use abs_rho_spinless, which is defined on the aims grid!

        implicit none

        integer :: denspoints
        real*8 :: numer_grad, denom_grad, innerterm, x, part1, part2
        real*8, dimension(denspoints), intent(in) :: rho, dielec_func
        real*8, dimension(3, denspoints), intent(in) :: rho_grad
        real*8, dimension(3, denspoints), intent(out) :: dielec_func_grad
        integer :: i_coords, i_full_points
        real*8 :: t_dielec_deriv
        real*8, dimension(denspoints) :: gradient_abs
        real*8 :: total_rho_abs

        if (dieleckind.eq.0) then
            do i_full_points = 1, denspoints, 1
                total_rho_abs = abs(rho(i_full_points))

                part1 = 2.0d0*beta_eps*(total_rho_abs/rho0_eps)**(2*beta_eps-1.0d0)

                part2 = (total_rho_abs/rho0_eps)**(2*beta_eps)

                dielec_func_grad(:,i_full_points) = (epsinf_mpb-1.0d0)/2.0d0*(-(part1*(1-part2))/&
                    (rho0_eps*(1.0+part2)**2)-(part1)/(rho0_eps*(1+part2)))*rho_grad(:,i_full_points)

                ! 	    if (gradient_abs(i_full_points).ne.gradient_abs(i_full_points)) then
                ! 	      write(use_unit,*) 'NaN found in dielecfunc_gradient_from_density_gradient', part1, part2
                ! 	      stop
                ! 	    end if
            end do
        else if (dieleckind.eq.1) then

            !calculate here directly grad{log(dielec_func)}
            do i_full_points = 1, denspoints, 1
                total_rho_abs = abs(rho(i_full_points))
                if (total_rho_abs.lt.rhomin_mpb.or.total_rho_abs.gt.rhomax_mpb) then
                    dielec_func_grad(:,i_full_points) = 0.0d0
                else
                    x = log(total_rho_abs)

                    ! 	      part1 = log(rhomax_mpb) - log(rhomin_mpb)
                    ! 	      dielec_func_grad(:,i_full_points) = &
                        ! 		- 2.0d0*log(epsinf_mpb)*sin(pi*(log(rhomax_mpb)-x)/part1)**2*1.0d0/(abs(rho(i_full_points))*part1)*&
                        ! 		rho_grad(:,i_full_points)

                    part1 = log_rhomax_mpb - log_rhomin_mpb

                    part1 = -log_epsinf_mpb/part1*(1-cos(pi4*0.5d0*(log_rhomax_mpb-x)/part1))

                    dielec_func_grad(:,i_full_points) = exp(t_dielec(x,rhomin_mpb,rhomax_mpb,epsinf_mpb))*&
                        part1*1.0d0/rho(i_full_points)*rho_grad(:,i_full_points)

                end if
                gradient_abs(i_full_points) = sqrt(sum(dielec_func_grad(:,i_full_points)**2))
            end do

        end if

        return

    end subroutine dielecfunc_gradient_from_density_gradient_points

    !****s* FHI-aims/alphafunc_from_density
    !  NAME
    !  alphafunc_from_density
    !  SYNOPSIS
    subroutine alphafunc_from_density( dens_array, dens_array_gradient, dielec_func_array, denspoints)
    !  PURPOSE
    ! Evaluates the exclusion function alpha, the function kappa and all respective derivatives with respect
    ! to electron density and the gradients.

        implicit none

        integer :: denspoints
        integer ::  i,i_spin

        real*8, dimension(n_spin,denspoints), intent(in) :: dens_array
        real*8, dimension(3,n_spin,denspoints), intent(in) :: dens_array_gradient
        real*8, dimension(denspoints), intent(in) :: dielec_func_array
        real*8, dimension(denspoints) :: kappa_func
        real*8, dimension(3,denspoints) :: kappa_func_gradient
        real*8 :: dielec_func_in

        real*8, dimension(denspoints) :: kappa_func_anion
        real*8, dimension(3,denspoints) :: kappa_func_anion_gradient

        real*8 :: rho0_kappa, total_rho
        real*8 :: total_rho_grad(3)

        real*8 :: part1

        real*8 :: shift_dist !distance that corresponds to the rescaling of the logarithmic distance between rhomin and rhomax

        real*8 :: dkappa_drho

        real*8 :: total_rho_abs, log_total_rho_abs

        if (alphakind.eq.0) then !set kappa as a sharp step function

            do i = 1, denspoints, 1
                dielec_func_in = dielec_func_array(i)
                if (dielec_func_in.le.(epsinf_mpb-1.0d0)*offset_kappa_mpb+1.0d0) then
                    kappa_func(i) = 0.0d0
                else
                    kappa_func(i) = kappainf_mpb
                end if
            end do

            kappa_func_gradient = 0d0

        else if (alphakind.eq.1) then !set kappa as a smooth step function, similar to dielecfunc_from_density

            if (dieleckind.eq.0) then

                rho0_kappa = abs((1-offset_kappa_mpb)/offset_kappa_mpb)**(1.0/(2.0*beta_eps)) * rho0_eps

                ! 	    rho0_kappa = 3**(1/(2*beta_eps))*rho0_eps*(sqrt(epsinf_mpb)/sqrt((epsinf_mpb-4)))**(1/beta_eps)
                ! this rho0_kappa leads to the smooth transition at the points where the dielectric function has just 1/4
                ! of its bulk value from kappainf_mpb (bulk) to zero inside cavity


                do i = 1, denspoints, 1
                    total_rho = 0.0d0
                    do i_spin = 1, n_spin, 1
                        total_rho = total_rho + dens_array(i_spin,i)
                    end do
                    kappa_func(i) = kappainf_mpb/2 * (1 + (1-(abs(total_rho)/rho0_kappa)&
                        **(2*beta_eps))/(1+(abs(total_rho)/rho0_kappa)**(2*beta_eps)))
                end do

            else if (dieleckind.eq.1) then

                !kappa function has same behaviour as the dielectric function kind = 1. as long as max and min rho0_eps
                !are identical

                do i = 1, denspoints, 1
                    total_rho = 0.0d0
                    total_rho_grad = 0d0
                    do i_spin = 1, n_spin, 1
                        total_rho = total_rho + dens_array(i_spin,i)
                        total_rho_grad(:) = total_rho_grad(:) + dens_array_gradient(:,i_spin,i)
                    end do
                    total_rho_abs = abs(total_rho)
                    log_total_rho_abs = log(total_rho_abs)
                    if (total_rho_abs.lt.rhomin_kappa) then
                        kappa_func(i) = kappainf_mpb
                        dalpha_drho(i) = 0d0
                        kappa_func_gradient(:,i) = 0d0
                    else if (total_rho_abs.gt.rhomax_kappa) then
                        kappa_func(i) = 0.0d0
                        dalpha_drho(i) = 0d0
                        kappa_func_gradient(:,i) = 0d0
                    else
                        kappa_func(i) = (exp(t_dielec(log_total_rho_abs,rhomin_kappa,rhomax_kappa,epsinf_mpb))-1.0d0)&
                            *(kappainf_mpb/(epsinf_mpb-1.0d0))
                        part1 = log_rhomax_kappa - log_rhomin_kappa
                        part1 = -log_epsinf_mpb/part1*(1-cos(pi4*0.5d0*(log_rhomax_kappa-log_total_rho_abs)/part1))

                        dkappa_drho = exp(t_dielec(log_total_rho_abs,rhomin_kappa,rhomax_kappa,epsinf_mpb))*&
                            part1*1.0d0/total_rho*(kappainf_mpb/(epsinf_mpb-1.0d0))

                        kappa_func_gradient(:,i) = dkappa_drho*total_rho_grad(:)

                        !evaluate d(kappa_func**2*z_mpb/kBT_mpb)/drho
                        dalpha_drho(i) = 1d0/kappainf_mpb*dkappa_drho

                    end if
                    if (use_separate_alpha) then
                        if (total_rho_abs.lt.rhomin_kappa_anion) then
                            kappa_func_anion(i) = kappainf_mpb
                            kappa_func_anion_gradient(:,i) = 0d0
                            dalpha_anion_drho(i) = 0d0
                        else if (total_rho_abs.gt.rhomax_kappa_anion) then
                            kappa_func_anion(i) = 0.0d0
                            kappa_func_anion_gradient(:,i) = 0d0
                            dalpha_anion_drho(i) = 0d0
                        else
                            kappa_func_anion(i) = (exp(t_dielec(log_total_rho_abs,rhomin_kappa_anion,rhomax_kappa_anion,epsinf_mpb))-1.0d0)&
                                *(kappainf_mpb/(epsinf_mpb-1.0d0))
                            part1 = log_rhomax_kappa_anion - log_rhomin_kappa_anion
                            part1 = -log_epsinf_mpb/part1*(1-cos(pi4*0.5d0*(log_rhomax_kappa_anion-log_total_rho_abs)/part1))

                            dkappa_drho = exp(t_dielec(log_total_rho_abs,rhomin_kappa_anion,rhomax_kappa_anion,epsinf_mpb))*&
                                part1*1.0d0/total_rho*(kappainf_mpb/(epsinf_mpb-1.0d0))

                            kappa_func_anion_gradient(:,i) = dkappa_drho*total_rho_grad(:)
                            dalpha_anion_drho(i) = 1d0/kappainf_mpb*dkappa_drho
                        end if
                    end if
                end do

            end if

        else

            write(use_unit,*) 'chose other alphakind'
            stop

        end if

        alpha_func_mpb = kappa_func/kappainf_mpb
        alpha_func_gradient_mpb = kappa_func_gradient/kappainf_mpb

        if (use_separate_alpha) then
            alpha_func_anion_mpb = kappa_func_anion/kappainf_mpb
            alpha_func_anion_gradient_mpb = kappa_func_anion_gradient/kappainf_mpb
        end if

        return

    end subroutine alphafunc_from_density

    !****s* FHI-aims/kappafunc_from_density_points
    !  NAME
    !  kappafunc_from_density_points
    !  SYNOPSIS
    subroutine kappafunc_from_density_points( dens_array, dens_array_gradient, dielec_func_array, denspoints, kappa_func,&
            kappa_func_gradient,cutmin,cutmax)
    !  PURPOSE
    ! same as kappafunc_from_density but on arbitrary grid
        implicit none

        integer :: denspoints
        integer ::  i

        real*8, dimension(denspoints), intent(in) :: dens_array
        real*8, dimension(3,denspoints), intent(in) :: dens_array_gradient
        real*8, dimension(denspoints), intent(in) :: dielec_func_array
        real*8, dimension(denspoints), intent(out) :: kappa_func
        real*8, dimension(3,denspoints), intent(out) :: kappa_func_gradient
        real*8 :: dielec_func_in, cutmin,cutmax

        real*8 :: rho0_kappa, total_rho
        real*8 :: total_rho_grad(3)

        real*8 :: x, part1

        real*8 :: shift_dist !distance that corresponds to the rescaling of the logarithmic distance between rhomin and rhomax

        real*8 :: dkappa_drho

        real*8 :: total_rho_abs

        if (alphakind.eq.0) then !set kappa as a sharp step function

            do i = 1, denspoints, 1
                dielec_func_in = dielec_func_array(i)
                if (dielec_func_in.le.(epsinf_mpb-1.0d0)*offset_kappa_mpb+1.0d0) then
                    kappa_func(i) = 0.0d0
                else
                    kappa_func(i) = kappainf_mpb
                end if
            end do

            kappa_func_gradient = 0d0

        else if (alphakind.eq.1) then !set kappa as a smooth step function, similar to dielecfunc_from_density

            if (dieleckind.eq.0) then

                rho0_kappa = abs((1-offset_kappa_mpb)/offset_kappa_mpb)**(1.0/(2.0*beta_eps)) * rho0_eps

                ! 	    rho0_kappa = 3**(1/(2*beta_eps))*rho0_eps*(sqrt(epsinf_mpb)/sqrt((epsinf_mpb-4)))**(1/beta_eps)
                ! this rho0_kappa leads to the smooth transition at the points where the dielectric function has just 1/4
                ! of its bulk value from kappainf_mpb (bulk) to zero inside cavity


                do i = 1, denspoints, 1
                    total_rho = dens_array(i)
                    kappa_func(i) = kappainf_mpb/2 * (1 + (1-(abs(total_rho)/rho0_kappa)&
                        **(2*beta_eps))/(1+(abs(total_rho)/rho0_kappa)**(2*beta_eps)))
                end do

            else if (dieleckind.eq.1) then

                !kappa function has same behaviour as the dielectric function kind = 1. as long as max and min rho0_eps
                !are identical

                do i = 1, denspoints, 1
                    total_rho = dens_array(i)
                    total_rho_grad = dens_array_gradient(:,i)
                    if (denspoints.eq.n_full_points) then
                        total_rho_abs = abs(total_rho) !abs_rho_spinless(i)
                    else
                        total_rho_abs = abs(total_rho)
                    end if
                    if (total_rho_abs.lt.cutmin) then
                        kappa_func(i) = kappainf_mpb
                        kappa_func_gradient(:,i) = 0d0
                    else if (total_rho_abs.gt.cutmax) then
                        kappa_func(i) = 0.0d0
                        kappa_func_gradient(:,i) = 0d0
                    else
                        kappa_func(i) = (exp(t_dielec(log(abs(total_rho)),cutmin,cutmax,epsinf_mpb))-1.0d0)&
                            *(kappainf_mpb/(epsinf_mpb-1.0d0))
                        x = log(total_rho)
                        part1 = log(cutmax) - log(cutmin)
                        part1 = -log_epsinf_mpb/part1*(1-cos(pi4*0.5d0*(log(cutmax)-x)/part1))

                        dkappa_drho = exp(t_dielec(x,cutmin,cutmax,epsinf_mpb))*&
                            part1*1.0d0/total_rho*(kappainf_mpb/(epsinf_mpb-1.0d0))

                        kappa_func_gradient(:,i) = dkappa_drho*total_rho_grad(:)

                    end if
                end do

            end if

        else

            write(use_unit,*) 'chose other alphakind'
            stop

        end if


        return

    end subroutine kappafunc_from_density_points

    !****s* FHI-aims/evaluate_q_times_epsinf
    !  NAME
    !  evaluate_q_times_epsinf
    !  SYNOPSIS
    subroutine evaluate_q_times_epsinf(rho,  &
            v_free_gradient_public, &
            dielec_func, dielec_func_gradient_mpb, q_times_epsinf,&
            delta_rho_lpb,q_times_epsinf_mpb_calc,use_mpbe_params)
    !  PURPOSE
    ! Evaluates q * epsinf_mpb, the part of the source term of the SPE that remains fixed during
    ! the MERM
        use physics, only: free_rho_superpos

        implicit none

        logical :: use_mpbe_params !do we solve MPBE or LPBE?
        real*8, dimension(n_full_points), intent(out) :: q_times_epsinf
        real*8, dimension(n_full_points), intent(in) :: rho
        real*8, dimension(n_full_points), intent(in) :: dielec_func
        real*8, dimension(3,n_full_points), intent(in) :: dielec_func_gradient_mpb
        real*8, dimension(3,n_full_points), intent(in) :: v_free_gradient_public
        real*8, dimension(n_full_points), intent(in) :: q_times_epsinf_mpb_calc
        real*8 :: grad_dot_product
        real*8 :: delta_rho_lpb(n_full_points)
        real*8 :: part_gradient_comp

        integer :: i_full_points, i_coords


        q_times_epsinf = 0d0
        if (use_mpbe_params) then
            !we already evaluated eps_s*q in run_newton
            q_times_epsinf = q_times_epsinf_mpb_calc
        else
            !evaluate eps_s*q^LPB from the paper
            do i_full_points = 1, n_full_points, 1
                if (partition_tab(i_full_points).gt.0d0) then
                    grad_dot_product = 0.0d0
                    do i_coords = 1, 3, 1
                        grad_dot_product =  grad_dot_product +&
                            dielec_func_gradient_mpb(i_coords,i_full_points)*&
                            v_free_gradient_public(i_coords,i_full_points)
                    end do
                    q_times_epsinf(i_full_points) = epsinf_mpb/dielec_func(i_full_points) * (&
                        (1.0d0-dielec_func(i_full_points))* pi4_inv * free_rho_superpos(i_full_points) +&
                        delta_rho_lpb(i_full_points) + &
                        pi4_inv*(grad_dot_product-&
                        alpha_func_mpb(i_full_points)*1d0/(1d0+phi_zero_mpb*(alpha_func_mpb(i_full_points)-1d0))*&
                        kappainf_mpb**2*z_mpb/kBT_mpb*&
                        free_hartree_superpos(i_full_points)))
                end if
            end do
        end if !use_mpbe_params
        return

    end subroutine evaluate_q_times_epsinf

    !****s* FHI-aims/evaluate_SPE_source
    !  NAME
    !  evaluate_SPE_source
    !  SYNOPSIS
    subroutine evaluate_SPE_source(q_times_epsinf, dielec_func, dielec_func_grad, &
            delta_v_MPB, delta_v_MPB_gradient, rho_iter, &
            use_mpb_ions,cutf,use_mpbe_params,h_function_calc)
    !  PURPOSE
    ! Evaluates qSPE = q_times_epsinf - pi4_inv * (L_1 * delta_v_MPB) * epsinf_mpb
    !purpose: evaluate L_1*delta v and add it to q to get total source term of SPE
    !output: rho_iter = eps_s *    (q - 1/4pi*L_1 * delta v)
    !	- in kerker_mpb this is again devided by eps_s to get the source term "source"
    !	- this was multiplied by -1/4pi to get the source term of the integral >>int exp(-kappa r')/(r-r')* "source"<< which is then multipol expanded
    !		, with "source" = q - 1/4pi*L_1*delta v

        use constants, only : pi4_inv

        implicit none


        real*8 :: rho_iter_inner, kappa_current, rho_rescale
        integer :: i_full_points,  i_coords
        real*8 :: cutf

        !INPUTS
        real*8, dimension(n_full_points), intent(in) :: dielec_func
        real*8, dimension(n_full_points), intent(in) :: delta_v_MPB
        real*8, dimension(3,n_full_points), intent(in) :: delta_v_MPB_gradient
        real*8, dimension(3,n_full_points), intent(in) :: dielec_func_grad
        real*8, dimension(n_full_points), intent(in) :: q_times_epsinf
        logical, intent(in) :: use_mpb_ions, use_mpbe_params
        real*8, dimension(n_full_points), intent(in) :: h_function_calc !part that is set in each Newton step

        real*8, dimension(n_full_points) :: rho_iter_inner_ar
        real*8, dimension(n_full_points) :: rho_rescale_ar

        !
        !OUTPUTS
        real*8, dimension( n_full_points), intent(out) :: rho_iter

        rho_iter = 0d0


        do i_full_points = 1, n_full_points, 1

            if (partition_tab(i_full_points).gt.0d0) then

                if (use_mpbe_params) then
                    !MPB case: evaluate eps_s*(q - 1/4pi*L_1 * delta v)
                    !1st part of eps_s*1/4pi*L_1*delta v : eps_s*1/4pi*(kappa^2-h^2/eps) delta v
                    rho_iter(i_full_points) = q_times_epsinf(i_full_points) + epsinf_mpb*pi4_inv*( &
                        (kappa_debye**2-h_function_calc(i_full_points)/dielec_func(i_full_points))*delta_v_MPB(i_full_points))
                    do i_coords = 1, 3, 1
                        !2nd part of eps_s*1/4pi*L_1*delta v : eps_s*1/4pi*(nabla ln(eps)) (nabla delta v)
                        rho_iter(i_full_points) = rho_iter(i_full_points) + &
                        epsinf_mpb*pi4_inv*(1d0/dielec_func(i_full_points) *&
                            dielec_func_grad(i_coords,i_full_points)*delta_v_MPB_gradient(i_coords,i_full_points))
                    end do
                else
                    !LPB case: evaluate eps_s*(q^LPB -1/4pi L_1^LPB*delta v)
                    rho_iter_inner = 0.0d0

                    if (dieleckind.eq.0) then
                        ! 	      if (sum(dielec_func_grad(:, i_full_points)**2).gt.1.0d-30) then !only if dielec_func_grad is not zero anyways
                        do i_coords = 1, 3, 1
                            rho_iter_inner = rho_iter_inner +  pi4_inv *&
                                epsinf_mpb/dielec_func(i_full_points) * dielec_func_grad(i_coords, i_full_points)*&
                                delta_v_MPB_gradient(i_coords,i_full_points)
                        end do
                        ! 	      end if
                    else if (dieleckind.eq.1) then
                        ! 	      if (sum(dielec_func_grad(:, i_full_points)**2).gt.1.0d-30) then
                        do i_coords = 1, 3, 1
                            rho_iter_inner = rho_iter_inner +  pi4_inv *&
                                epsinf_mpb/dielec_func(i_full_points) * dielec_func_grad(i_coords, i_full_points)*&
                                ! 		      delta_v_MPB(i_full_points)
                                delta_v_MPB_gradient(i_coords,i_full_points)
                        end do
                        ! 	      end if
                    end if

                    rho_rescale = 0.0d0

                    kappa_current = kappainf_mpb / sqrt(dielec_func(i_full_points)) * sqrt(z_mpb/kBT_mpb)
                    rho_rescale = epsinf_mpb * pi4_inv * (kappa_debye**2- alpha_func_mpb(i_full_points)*&
                        1d0/(1d0+phi_zero_mpb*(alpha_func_mpb(i_full_points)-1d0))*kappa_current**2) * delta_v_MPB(i_full_points)

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    rho_iter(i_full_points) = q_times_epsinf(i_full_points) + & !normal scaled density
                        rho_iter_inner + & !terms for nonlinear dielectric/kappa (1.0d0-exp(-i_iterate_lpb_outer/2.0d0)) *
                        rho_rescale
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    rho_iter_inner_ar(i_full_points) = rho_iter_inner
                    rho_rescale_ar(i_full_points) = rho_rescale
                end if !use_mpbe_params
            end if !partition_tab
        end do


    end subroutine evaluate_SPE_source

    !****s* FHI-aims/integrate_effective_density
    !  NAME
    !  integrate_effective_density
    !  SYNOPSIS
    subroutine integrate_effective_density(rho_solve, integral)
    !  PURPOSE
    ! Explicit integration of qSPE to get the potential. Only for testing purposes. This is superexpensive
    ! and makes not at all sense to use it. If one however runs into problems with the SPE_solver, one coul
    ! try to use this routine to find errors.

        !does the same integration that the kerker_mpb does. but does not use multipole expansion,
        !but explicitly evaluates the two non-local integral

        !this is VERY inefficient and just for testing purposes

        ! 	use dimensions
        ! 	use species_data
        use grids

        integer :: i_my_batch2, i_index2, i_my_batch, i_index, i_full_points, i_full_points2,&
            i_coords, i_spin

        real*8, dimension(n_spin, n_full_points) :: rho_solve
        real*8, dimension(n_full_points), intent(out) :: integral

        real*8 :: greens_func

        real*8 :: lambda_mpb,dist, dist_inv

        real*8, dimension(3) :: coord_current
        real*8, dimension(3) :: coord_current2

        integral = 0.0d0

        i_full_points = 0

        lambda_mpb = kappa_div_sqrt_eps*sqrt(z_mpb/kBT_mpb)

        write(use_unit,*) 'Starting with explicit integration of modified Helmholtz equation'

        do i_my_batch = 1, n_my_batches, 1
            do i_index = 1, batches(i_my_batch)%size, 1
                i_full_points   = i_full_points + 1

                coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)

                if (mod(i_full_points,100)==0) then
                    write(use_unit,*) 'Finished ', float(i_full_points)/float(n_full_points)*100., ' %'
                end if

                i_full_points2 = 0

                do i_my_batch2 = 1, n_my_batches, 1
                    do i_index2 = 1, batches(i_my_batch2)%size, 1


                        i_full_points2   = i_full_points2 + 1

                        if ((partition_tab(i_full_points2).gt.0.0d0).and.&
                                .not.(i_my_batch == i_my_batch2 .and. i_index == i_index2)) then
                            !
                            coord_current2(:) = batches(i_my_batch2) % points(i_index2) % coords(:)

                            dist = 0.0d0
                            do i_coords = 1, 3, 1
                                dist = dist + (coord_current(i_coords)-coord_current2(i_coords))**2
                            end do
                            dist = sqrt(dist)
                            dist_inv = 1.0d0/dist

                            greens_func = 1.0d0/epsinf_mpb * exp(-lambda_mpb * dist) * dist_inv

                            do i_spin = 1, n_spin, 1
                                integral(i_full_points) = integral(i_full_points) - &
                                    rho_solve(i_spin,i_full_points2) * greens_func * &
                                    partition_tab(i_full_points2)
                            end do
                        end if

                    end do
                end do
            end do
        end do


    end subroutine integrate_effective_density

    !****s* FHI-aims/mix_delta_v_MPB
    !  NAME
    !  mix_delta_v_MPB
    !  SYNOPSIS
    subroutine mix_delta_v_merm(delta_v_MPB, delta_v_MPB_old, density_mixing_lpb)
    !  PURPOSE
    ! Linear mixing of delta_v_MPB of previous step and current
        implicit none

        integer :: i_full_points

        real*8, dimension(n_full_points), intent(inout) :: delta_v_MPB
        real*8, dimension(n_full_points), intent(in) :: delta_v_MPB_old
        real*8 :: density_mixing_lpb

        ! 	density_mixing_lpb = 0.3d0 !defines percentage of new potential that should be used

        do i_full_points = 1, n_full_points
            delta_v_MPB(i_full_points) = density_mixing_lpb * delta_v_MPB(i_full_points) +&
                (1.0d0 - density_mixing_lpb) * delta_v_MPB_old(i_full_points)
        end do

        return

    end subroutine mix_delta_v_merm

    !****s* FHI-aims/evaluate_rmsd_merm
    !  NAME
    !  evaluate_rmsd_merm
    !  SYNOPSIS
    subroutine evaluate_rmsd_merm(delta_v_MPB, delta_v_MPB_old, rmsd,  tau_merm, lpb_converged)
    !  PURPOSE
    ! Evaluate root mean square change of delta_v_MPB during iteration

        use mpi_utilities
        use synchronize_mpi

        implicit none

        integer :: i_full_points,i_myid,info
        real*8 :: rmsd!, norm
        real*8 :: rmsd_mpi(n_tasks)
        real*8 :: norm_mpi(n_tasks)

        real*8, dimension(n_full_points), intent(in) :: delta_v_MPB
        real*8, dimension(n_full_points), intent(in) :: delta_v_MPB_old

        real*8 :: tau_merm

        logical, intent(out) :: lpb_converged

        rmsd = 0.0d0
        ! 	norm = 0.0d0

        i_full_points = 0

        do i_myid = 0, n_tasks-1, 1
            if (myid == i_myid) then
                rmsd = 0
                do i_full_points = 1, n_full_points, 1
                    rmsd = rmsd + (delta_v_MPB(i_full_points)- delta_v_MPB_old(i_full_points))**2 * partition_tab(i_full_points)
                    ! 			norm = norm + partition_tab(i_full_points)
                end do
            end if

            call mpi_barrier(mpi_comm_global,info)
        end do

        call sync_real_number(rmsd) !sums all rmsd on different processors and saves them to rmsd
        ! 	call sync_real_number(norm)

        ! 	rmsd = sqrt(rmsd/norm) !/n_full_points
        rmsd = sqrt(rmsd) !/n_full_points


        if (rmsd.lt.tau_merm) then
            lpb_converged = .True.
        else
            lpb_converged = .False.
        end if

        if (rmsd.ne.rmsd) then
            if (myid==0) write(use_unit,*) 'RMSD: ', rmsd
            if (myid==0) write(use_unit,*)  'The following parameters were applied, nmin =',rhomin_mpb, ', nmax =',rhomax_mpb, ' ,nmin^alpha = ',&
                rhomin_kappa, ' ,nmax^alpha = ',rhomax_kappa, ' ,d_alpha = ',d_alpha_ion,' ,xi_alpha = ',xi_alpha_ion
            call aims_stop('NaN in LPBE solver, RMSD calculation. Reconsider your dielectric functions settings.' )
        end if

        return

    end subroutine evaluate_rmsd_merm

    ! !****s* FHI-aims/evaluate_dnion_dalpha
    ! !  NAME
    ! !  evaluate_dnion_dalpha
    ! !  SYNOPSIS
    ! subroutine evaluate_dnion_dalpha(new_potential)
    ! !  PURPOSE
    ! ! Evaluate derivative of ion density nion with respect to alpha
    ! !
    ! !calculate dnion/dalpha only for alpha>0.0, since we only need this
    ! !	function multiplied by dalpha/drho which is exactly zero
    ! !	for dalpha/drho = 0
    !   integer :: i_full_points
    !   real*8, dimension(n_full_points) :: new_potential
    !   real*8 :: current_pot, current_alpha
    !
    !   dnion_dalpha = 0d0
    !
    !   do i_full_points = 1, n_full_points, 1
    !     if ((partition_tab(i_full_points).gt.0d0).and.&
        !       (abs(dalpha_drho(i_full_points)).gt.0d0)) then
    !       current_pot = new_potential(i_full_points)
    !       if (solve_lpbe_only) then
    ! 	if (solve_debye_only) then
    ! 	  dnion_dalpha(i_full_points) = &
        ! 	    -pi4_inv*kappa_debye**2 * epsinf_mpb*&
        ! 	    current_pot
    ! 	else
    ! 	  current_alpha = alpha_func_mpb(i_full_points)
    ! 	  dnion_dalpha(i_full_points) = &
        ! 	    -pi4_inv*kappa_debye**2 * epsinf_mpb*&
        ! 	    (1d0+phi_zero_mpb)/&
        ! 	    (1d0+phi_zero_mpb*(current_alpha-1d0))**2*&
        ! 	    current_pot
    ! 	end if
    !       else
    ! 	if (solve_pbe_only) then
    ! 	  dnion_dalpha(i_full_points) =&
        ! 	    -2d0*z_mpb*c_mpb*sinh(z_mpb/kBT_mpb*current_pot)
    ! 	else
    ! 	  current_alpha = alpha_func_mpb(i_full_points)
    ! 	  dnion_dalpha(i_full_points) =&
        ! 	    -2d0*z_mpb*c_mpb*(1d0-phi_zero_mpb)*&
        ! 	    sinh(z_mpb/kBT_mpb*current_pot)/&
        ! 	    (1d0-phi_zero_mpb+phi_zero_mpb*current_alpha*&
        ! 	    cosh(z_mpb/kBT_mpb*current_pot))**2
    ! 	end if
    !       end if
    !     end if
    !   end do
    !   return
    ! end subroutine evaluate_dnion_dalpha

    ! !****s* FHI-aims/evaluate_pulay_weight_mpb
    ! !  NAME
    ! !  evaluate_pulay_weight_mpb
    ! !  SYNOPSIS
    ! subroutine evaluate_pulay_weight_mpb(i_point_to_fullpoint,n_points,pulay_weight_mpb)
    ! !  PURPOSE
    ! ! Old routine
    ! !
    !   !evaluate weight for additional pulay force caused by ions. only calculated for partition_tab>0 and&
        !   !alpha_func > 0 (and not mpbe_no_ion)
    !   !n_points are all point of current batch, for which partition_tab .gt. 0d0
    !   integer :: n_points,i_point
    !   integer, dimension(n_points), intent(in):: i_point_to_fullpoint
    !   real*8, dimension(n_points), intent(out) :: pulay_weight_mpb
    !
    !   pulay_weight_mpb = 0d0
    !
    !   do i_point = 1, n_points, 1
    !       if (abs(dalpha_drho(i_point_to_fullpoint(i_point))).gt.0d0) then
    ! 	pulay_weight_mpb(i_point) = &
        ! 	  dnion_dalpha(i_point_to_fullpoint(i_point)) * &
        ! 	  dalpha_drho(i_point_to_fullpoint(i_point))
    !       end if
    !   end do
    !
    !   return
    !
    ! end subroutine evaluate_pulay_weight_mpb

    !****s* FHI-aims/scf_mul_vec_1
    !  NAME
    !  scf_mul_vec_1
    !  SYNOPSIS
    subroutine scf_mul_vec_1 ( wave, n_mul, ylm, vec_2 )

    !  PURPOSE
    !  write an explicitly vectorizable multiplication to avoid an index mess
    !
    !  USES
        implicit none
    !  ARGUMENTS


        integer :: n_mul, i_mul
        real*8 :: wave(n_mul)
        real*8 :: ylm(n_mul)
        real*8 :: vec_2(n_mul)

    !  INPUTS
    !    o n_mul -- dimentsion of the vectors
    !    o ylm -- vector which is multiplied
    !    o factor -- factor which is multiplied
    !  OUTPUT
    !    o wave -- results of the multiplication
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

        do i_mul = 1, n_mul, 1
            wave(i_mul) = ylm(i_mul) * vec_2(i_mul)
        end do

    end subroutine scf_mul_vec_1

    !****s* FHI-aims/evaluate_v_free_gradient
    !  NAME
    !  evaluate_v_free_gradient
    !  SYNOPSIS
    subroutine evaluate_v_free_gradient
    !  PURPOSE
    ! Evaluate the gradient of the superposition of free atom potentials

        use dimensions
        use pbc_lists
        use runtime_choices
        use species_data
        use physics, only: multipole_radius_sq !???? why do we need this, this should be the end of the free atom grid...
        use geometry
        use grids
        use spline
        use free_atoms
        use load_balancing
        use mpi_tasks, only: check_allocation

        integer :: i_full_points,i_index,i_center,current_center,&
            current_spl_atom,i_batch,info,n_my_batches_work, i_atom
        real*8, dimension(3) :: dir_tab,dir_tab_in,dir_tab_out,coord_current
        real*8 :: dist_tab_in,dist_tab_out,dist_tab_sq,i_r,i_r_log,radial_weight,&
            log_weight,d_v_hartree_free_d_r
        real*8, dimension(:),allocatable :: adap_outer_radius_sq
        type (batch_of_points), pointer :: batches_work(:)
        integer n_bp
        real*8, dimension(n_atoms) :: multipole_radius_sq_local !this is the free atom multipole radius, instead for the case that
        !we do periodic calculation, for this case this has to be corrected, see integrate_hartree_log_grid

        do i_atom = 1, n_atoms, 1
            multipole_radius_sq_local(i_atom) = multipole_radius_free(species(i_atom))**2
        end do

        n_bp = use_batch_permutation
        if(use_batch_permutation > 0) then
            n_my_batches_work = batch_perm(n_bp)%n_my_batches
            batches_work => batch_perm(n_bp)%batches
        else
            n_my_batches_work = n_my_batches
            batches_work => batches
        endif

        allocate(adap_outer_radius_sq(n_atoms),stat=info)
        call check_allocation(info, 'adap_outer_radius_sq          ')

        adap_outer_radius_sq = 1e8

        v_free_gradient_public = 0d0
        v_free_gradient_public_atom = 0d0
        !determine superposition of free atom potentials and gradient and atomic components and save into variables
        do i_center = 1, n_centers_hartree_potential, 1
            current_center   = centers_hartree_potential(i_center)
            current_spl_atom = center_to_atom(current_center)
            i_full_points = 0

            do i_batch = 1, n_my_batches_work
                ! loop over one batch
                do i_index = 1, batches_work(i_batch)%size, 1
                    ! i_full_points is the index that indicates where we are in the entire grid (for external quanities like rho, potential, ...)
                    i_full_points = i_full_points + 1
                    ! Only execute if partition_tab is .gt. zero, else
                    ! we can run into 1/0 !!!
                    if (partition_tab(i_full_points).gt.0.d0) then
                        if (.not.empty(current_spl_atom)) then
                            coord_current(:) = batches_work(i_batch) % points(i_index) % coords(:)

                            call tab_single_atom_centered_coords_p0 &
                                ( current_center, &
                                coord_current,  &
                                dist_tab_sq,  &
                                dir_tab )
                            if (dist_tab_sq.lt.multipole_radius_sq_local(current_spl_atom)) then
                                call tab_single_atom_centered_coords_radial_log_p0 &
                                    ( current_center, dist_tab_sq, dir_tab,  &
                                    dist_tab_in, i_r, i_r_log, dir_tab_in )
                                call tab_single_radial_weights_v2 &
                                    ( current_spl_atom, dist_tab_in, i_r, &
                                    log_weight, radial_weight )
                                if(n_periodic == 0 .and.(.not. force_new_functional))then
                                    d_v_hartree_free_d_r =  &
                                        val_spline_deriv( i_r_log, &
                                        free_pot_es_spl(1,1,species(current_spl_atom)),  &
                                        n_grid(species(current_spl_atom))) * log_weight &
                                        - species_z(species(current_spl_atom)) / (dist_tab_sq)
                                else
                                    d_v_hartree_free_d_r = &
                                        val_spline_deriv( i_r_log, &
                                        free_pot_es_spl(1,1,species(current_spl_atom)),  &
                                        n_grid(species(current_spl_atom))) * log_weight
                                end if
                                dir_tab_out = dir_tab_in
                            else if (dist_tab_sq .lt. adap_outer_radius_sq(current_spl_atom)) then
                                dist_tab_out = sqrt(dist_tab_sq)
                                dir_tab_out = dir_tab / dist_tab_out
                                if(.not. force_new_functional)then
                                    d_v_hartree_free_d_r =  &
                                        - species_z(species(current_spl_atom))/dist_tab_sq
                                else
                                    d_v_hartree_free_d_r =  0.d0
                                end if
                            end if !far or short distance
                            v_free_gradient_public(:,i_full_points) =  v_free_gradient_public(:,i_full_points) +&
                                dir_tab_out(:) *d_v_hartree_free_d_r !* dist_tab_out
                            if (use_forces.and..not.forces_mpb_force_off) then
                                v_free_gradient_public_atom(current_center,:,i_full_points) =&
                                    v_free_gradient_public_atom(current_center,:,i_full_points) +&
                                    dir_tab_out(:) *d_v_hartree_free_d_r !* dist_tab_out
                            end if
                        end if   !empty(current_spl_atom)
                    end if !partition tab
                end do !index%batch
            end do !batches
        end do !center

    end subroutine evaluate_v_free_gradient

    !****s* lpb_solver_utilities/evaluate_mpb_energies
    !  NAME
    !  evaluate_mpb_energies
    !  SYNOPSIS
    subroutine evaluate_mpb_energies(hartree_delta_energy, en_elec_delta,&
            hartree_multipole_correction,hartree_multipole_error,&
            pot_current,total_rho, rho_multipole_current,&
            free_rho_superpos_current,v_hartree_free_current,&
            partition_tab_current,&
            !changed at output: energies:
            freeen_MPBE_nonelstat_solv_energy,&
            freeen_LPBE_nonelstat_solv_energy,&
            freen_solv_plus_ion,&
            freeen_sol_plus_solv_plus_ion,&
            freen_est,&
            f_function_current,&
            h_function_current,&
            ln_cosh_function_current,&
            cosh_cosh_function_current,&
            KS_mpb_corr1,KS_mpb_corr2,i_full_points)
    !  PURPOSE
    ! Evaluate different energy contributions in the MPB solver, as also the
    ! Kohn-Sham Hamiltonian corrections.
    !
        use constants, only: pi
        use mpi_utilities
    !  IN
        real*8 :: pot_current, total_rho, pot_old, total_rho_mpb,free_rho_superpos_current,&
            v_hartree_free_current,rho_multipole_current,partition_tab_current
        real*8 :: ion_freeen_1, ion_freeen_2
        real*8 :: freeen_MPBE_elstat, freeen_LPBE_elstat
        real*8 :: f_function_current,h_function_current
        real*8 :: ln_cosh_function_current, cosh_cosh_function_current
        integer :: i_full_points
        logical :: KS_mpb_corr1,KS_mpb_corr2
    !  INOUT
        real*8 :: hartree_delta_energy, en_elec_delta, hartree_multipole_correction,&
            hartree_multipole_error
        real*8 :: freeen_MPBE_nonelstat_solv_energy
        real*8 :: freeen_LPBE_nonelstat_solv_energy
        real*8 :: freen_solv_plus_ion
        real*8 :: freeen_sol_plus_solv_plus_ion
        real*8 :: freen_est
        real*8 :: Afunc, Nfunc, Nfunc_anion,dAfunc_drho, dNfunc_drho,&
            dNfunc_anion_drho, dAfunc_times_ln_Afunc_drho,&
            dNfunc_inv_times_ln_Nfunc,dNfunc_anion_inv_times_ln_Nfunc_anion,&
            exp_cat, exp_an, dAfunc_drho_term1,dAfunc_drho_term2
        real*8 :: abs_dalpha_drho_current, abs_dalpha_anion_drho_current
        real*8 :: cosh_func
        real*8 :: KS_corr_eps, KS_corr_alpha
        real*8 :: prefunc1,prefunc2, func1, func2

        !     //initialization

        abs_dalpha_drho_current = abs(dalpha_drho(i_full_points))
        abs_dalpha_anion_drho_current = abs(dalpha_anion_drho(i_full_points))
        KS_corr_alpha = 0d0
        KS_corr_eps = 0d0
        freeen_MPBE_elstat = 0d0
        freeen_LPBE_elstat = 0d0
        ion_freeen_2 = 0d0

        pot_old = pot_current

        !    ////////////////////////////////////////////////////////////
        !      ///////////////1. KS HAMILTONIAN CORRECTION
        !    ////////////////////////////////////////////////////////////

        if (abs(deps_drho(i_full_points)).gt.0d0.and.KS_mpb_corr1) then
            !we correct KS Hamiltonian...
            !now come the corrections to the Kohn-Sham operator, since epsilon and kappa functions depend on
            !the electron density
            !1st term: density dependence of eps
            KS_corr_eps =  - 1d0/(8d0*pi) * deps_drho(i_full_points) * &
                sum((v_free_gradient_public(:,i_full_points)+delta_v_MPB_gradient(:,i_full_points))**2)
        end if
        if (.not.use_separate_alpha) then
            if (abs_dalpha_drho_current.gt.0d0.and.KS_mpb_corr2.and..not.mpbe_no_ion) then
                !2nd term: density dependence of kappa
                if (.not.(solve_lpbe_only.and..not.use_mpbe_free_energy)) then
                    ! 	    	     !MPBE free energy
                    KS_corr_alpha = - kBT_mpb*2d0*c_mpb * dalpha_drho(i_full_points) * &
                        cosh_cosh_function_current
                else if (solve_lpbe_only.and..not.use_mpbe_free_energy) then
                    !d/dn(alpha c beta z^2 v^2) = dalpha/dn c beta z^2 v^2
                    !LPBE free energy
                    if (solve_debye_only) then
                        KS_corr_alpha = -&
                            kBT_mpb*2d0*c_mpb/(1d0+phi_zero_mpb/(1d0-phi_zero_mpb)*alpha_func_mpb(i_full_points))*&
                            1d0/(1d0-phi_zero_mpb)*dalpha_drho(i_full_points)-&
                            dalpha_drho(i_full_points)*c_mpb*z_mpb**2/kBT_mpb&
                            /(1d0+phi_zero_mpb*(alpha_func_mpb(i_full_points)-1d0))*&
                            pot_old**2 + &
                            alpha_func_mpb(i_full_points)*c_mpb*z_mpb**2/kBT_mpb*phi_zero_mpb*dalpha_drho(i_full_points)&
                            /(1d0+phi_zero_mpb*(alpha_func_mpb(i_full_points)-1d0))**2*&
                            pot_old**2
                    end if !solve_debye_only
                end if !solve_lpbe_only
            end if     !dalpha_drho>0
        else !use_separate_alpha case:
            !we need a bunch of terms for this...
            if (KS_mpb_corr2.and..not.mpbe_no_ion) then
                if (.not.solve_pbe_only) then
                    Afunc = 1d0
                    !1st we need the exponential functions only in special cases
                    exp_cat = 0d0
                    exp_an = 0d0
                    if (abs_dalpha_drho_current.gt.0d0.or.alpha_func_mpb(i_full_points).gt.0d0) then
                        exp_cat = exp(-z_mpb*pot_old/kBT_mpb)
                    end if
                    if (abs_dalpha_anion_drho_current.gt.0d0.or.alpha_func_anion_mpb(i_full_points).gt.0d0) then
                        exp_an = exp(z_mpb*pot_old/kBT_mpb)
                    end if

                    if (alpha_func_mpb(i_full_points).gt.0d0) then
                        cosh_func = cosh(z_mpb*pot_old/kBT_mpb)
                        Nfunc = 1d0-phi_zero_mpb+phi_zero_mpb*alpha_func_mpb(i_full_points)*cosh_func
                        Afunc = Afunc -phi_zero_mpb/2d0*alpha_func_mpb(i_full_points) * exp_cat/Nfunc
                    else
                        Nfunc = 1d0-phi_zero_mpb
                    end if
                    if (alpha_func_anion_mpb(i_full_points).gt.0d0) then
                        cosh_func = cosh(z_mpb*pot_old/kBT_mpb)
                        Nfunc_anion = 1d0-phi_zero_mpb+phi_zero_mpb*alpha_func_anion_mpb(i_full_points)*cosh_func
                        Afunc = Afunc -phi_zero_mpb/2d0*alpha_func_anion_mpb(i_full_points) * exp_an/Nfunc_anion
                    else
                        Nfunc_anion = 1d0-phi_zero_mpb
                    end if
                    if (abs_dalpha_drho_current.gt.0d0) then
                        cosh_func = cosh(z_mpb*pot_old/kBT_mpb)
                        dNfunc_drho = dalpha_drho(i_full_points) * phi_zero_mpb * cosh_func
                        dNfunc_inv_times_ln_Nfunc = dNfunc_drho * (1d0-log(Nfunc/(1-phi_zero_mpb)))/Nfunc**2
                    else
                        dNfunc_drho = 0d0
                        dNfunc_inv_times_ln_Nfunc = 0d0
                    end if
                    if (abs_dalpha_anion_drho_current.gt.0d0) then
                        cosh_func = cosh(z_mpb*pot_old/kBT_mpb)
                        dNfunc_anion_drho = dalpha_anion_drho(i_full_points) * phi_zero_mpb * cosh_func
                        dNfunc_anion_inv_times_ln_Nfunc_anion = dNfunc_anion_drho * (1d0-log(Nfunc_anion/(1-phi_zero_mpb)))/Nfunc_anion**2
                    else
                        dNfunc_anion_drho = 0d0
                        dNfunc_anion_inv_times_ln_Nfunc_anion = 0d0
                    end if
                    !this term should be now save to evaluate numerically for any grid point
                    dAfunc_drho_term1 = -phi_zero_mpb/2d0*(dalpha_drho(i_full_points)*Nfunc-alpha_func_mpb(i_full_points)*dNfunc_drho)/Nfunc**2
                    dAfunc_drho_term2 = -phi_zero_mpb/2d0*(dalpha_anion_drho(i_full_points)*Nfunc_anion-alpha_func_anion_mpb(i_full_points)*dNfunc_anion_drho)/Nfunc_anion**2

                    dAfunc_drho = 0d0
                    !saveguard, the exponential could explode inside the ionic cavity where the 2 terms are zero anyway
                    if (abs(dAfunc_drho_term1).gt.0d0) then
                        !better save then sorry, reevaluate the exponential, but we actually do not need it I think
                        dAfunc_drho = dAfunc_drho + exp_cat * dAfunc_drho_term1
                    end if
                    if (abs(dAfunc_drho_term2).gt.0d0) then
                        dAfunc_drho = dAfunc_drho + exp_an * dAfunc_drho_term2
                    end if

                    dAfunc_times_ln_Afunc_drho = (log(Afunc)+1d0)*dAfunc_drho
                    !now we have all the ingredients and we can calculate the actual corrections

                    KS_corr_alpha = kBT_mpb/a_mpb**3 * dAfunc_times_ln_Afunc_drho
                    if (abs_dalpha_drho_current.gt.0d0) then
                        KS_corr_alpha = KS_corr_alpha -kBT_mpb * c_mpb * exp_cat/Nfunc * dalpha_drho(i_full_points) * log(Nfunc/(1d0-phi_zero_mpb))
                    end if
                    if (abs_dalpha_anion_drho_current.gt.0d0) then
                        KS_corr_alpha = KS_corr_alpha -kBT_mpb * c_mpb * exp_an/Nfunc_anion * dalpha_anion_drho(i_full_points) * log(Nfunc_anion/(1d0-phi_zero_mpb))
                    end if
                    if (alpha_func_mpb(i_full_points).gt.0d0) then
                        KS_corr_alpha = KS_corr_alpha -kBT_mpb * c_mpb * exp_cat * alpha_func_mpb(i_full_points)*dNfunc_inv_times_ln_Nfunc
                    end if
                    if (alpha_func_anion_mpb(i_full_points).gt.0d0) then
                        KS_corr_alpha = KS_corr_alpha -kBT_mpb * c_mpb * exp_an * alpha_func_anion_mpb(i_full_points)*dNfunc_anion_inv_times_ln_Nfunc_anion
                    end if

                else !solve_pbe_only case
                    !this is much simpler...
                    KS_corr_alpha = 0d0
                    if (abs_dalpha_drho_current.gt.0d0) then
                        exp_cat = exp(-z_mpb*pot_old/kBT_mpb)
                        KS_corr_alpha = KS_corr_alpha -kBT_mpb * c_mpb * dalpha_drho(i_full_points) * exp_cat
                    end if
                    if (abs_dalpha_anion_drho_current.gt.0d0) then
                        exp_an = exp(z_mpb*pot_old/kBT_mpb)
                        KS_corr_alpha = KS_corr_alpha -kBT_mpb * c_mpb * dalpha_anion_drho(i_full_points) * exp_an
                    end if
                end if
                !!finished all corrections
            end if !KS_mpb_corr2.and..not.mpbe_no_ion
        end if !use_separate_alpha

        !add corrections to KS Hamiltonian
        pot_current = pot_current + KS_corr_eps
        pot_current = pot_current + KS_corr_alpha

        !     //ENERGY CORRECTIONS DUE TO KS CORRECTION

        if (KS_mpb_corr1) then
            !correct out the energies that correspond to the addition of
            !operators containing derivatives deps/drho and dkappa/drho to the hamiltonian.
            !they were added because of the
            !electron density dependence of the dielectric function and the kappa function
            en_eps_rho = en_eps_rho + KS_corr_eps*total_rho * partition_tab_current
        end if

        if (.not.mpbe_no_ion) then !only if we want to do DFT with ions
            ! 		 !electrostatic energy corrections, ionic interactions with Hartree potential
            if (KS_mpb_corr2) then
                !this will be substracted from the sum of eigenvalues in get_total_energy (belongs to kinetic energy!)
                en_alpha_rho = en_alpha_rho + KS_corr_alpha*total_rho*partition_tab_current
            end if
        end if


        !    ////////////////////////////////////////////////////////////
        !      ///////////////2. ENERGY CORRECTIONS
        !    ////////////////////////////////////////////////////////////

        !   //STANDARD AIMS TERMS INCLUDING MODIFIED MULTIPOLE CORRECTION

        hartree_delta_energy = hartree_delta_energy + &
            ( (pot_old) * &
            ( rho_multipole_current+ qSPE_plus_delta_rho(i_full_points) ) &
            - pi4_inv * free_rho_superpos_current * &
            v_hartree_free_current ) * partition_tab_current

        en_elec_delta = en_elec_delta + &
            ( (pot_old) * &
            ( rho_multipole_current + qSPE_plus_delta_rho(i_full_points) ) &
            - pi4_inv * free_rho_superpos_current * &
            v_hartree_free_current ) * partition_tab_current

        hartree_multipole_correction = &
            hartree_multipole_correction + &
            ( (pot_old) * &
            ( total_rho - rho_multipole_current - qSPE_plus_delta_rho(i_full_points) ) ) &
            *partition_tab_current


        !   //ADDITIONAL ENERGY CONTRIBUTIONS DUE TO IONS

        if (.not.mpbe_no_ion) then !only if we want to do DFT with ions
            !       //ELECTROSTATIC ENERGIES DUE TO IONS ARE ADDED TO HARTREE_DELTA_ENERGY

            !we add -1/8pi * int f(v) * v = 1/2*int nion*v to the total energy (electrostatic repulsion ion-v)
            if (alpha_func_mpb(i_full_points).gt.0d0) then
                if (solve_lpbe_only.and..not.use_mpbe_free_energy) then
                    !use LPBE electrostatics. here the ion-v interaction is exactly cancelled  by the nonelstatics!
                    freeen_LPBE_elstat = &
                        -pi4_inv/2d0*alpha_func_mpb(i_full_points)/(1d0-phi_zero_mpb*(alpha_func_mpb(i_full_points)-1d0))*&
                        kappainf_mpb**2*z_mpb/kBT_mpb*pot_old**2*partition_tab_current
                    hartree_delta_energy = hartree_delta_energy-2d0*freeen_LPBE_elstat
                else
                    !use MPBE electrostatics
                    freeen_MPBE_elstat = -pi4_inv/2d0*f_function_current*pot_old*&
                        partition_tab_current
                    hartree_delta_energy = hartree_delta_energy-2d0*freeen_MPBE_elstat
                end if
            end if

            !       //NONELECTROSTATIC (ENTROPIC) CORRECTIONS DUE TO IONS

            if (.not.(solve_lpbe_only.and..not.use_mpbe_free_energy)) then
                !nonelstat free energies according to MPBE free energy expression
                !free energy of ion and solvent WITHOUT solute (i.e. for v,n= 0):
                if (.not.solve_pbe_only) then
                    ion_freeen_1 = - kBT_mpb/a_mpb**3*&
                        log(1d0/(1d0-phi_zero_mpb))*partition_tab_current
                else !a=0 case
                    ion_freeen_1 = - 2d0*c_mpb*kBT_mpb*partition_tab_current
                end if
                if (.not.solve_pbe_only) then
                    !SR: we observed that in some cases the ion effect can be estimated by
                    !the following formula (neutral molecules) -- only experimental
                    !the ion effect can be estimated by int (1-alpha(r))*ion_freeen_1:
                    !old version: only valid if alpha is step function:
                    !freen_est = freen_est + (1d0-alpha_func_mpb(i_full_points))*&
                        !    kBT_mpb/a_mpb**3*&
                        !    log(1d0/(1d0-phi_zero_mpb))*partition_tab_current
                    if (use_separate_alpha) then
                        prefunc1 = (1d0-phi_zero_mpb+phi_zero_mpb*alpha_func_mpb(i_full_points))
                        prefunc2 = (1d0-phi_zero_mpb+phi_zero_mpb*alpha_func_anion_mpb(i_full_points))
                        func1 = alpha_func_mpb(i_full_points)/prefunc1
                        func2 = alpha_func_anion_mpb(i_full_points)/prefunc2
                        Afunc = 1d0-c_mpb*a_mpb**3*(func1+func2)
                        freen_est = freen_est + (-kBT_mpb*c_mpb*(&
                            func1*log(prefunc1/(1d0-phi_zero_mpb))+&
                            func2*log(prefunc2/(1d0-phi_zero_mpb)))+&
                            kBT_mpb/a_mpb**3*Afunc*log(Afunc)+kBT_mpb/a_mpb**3*log(1d0/(1d0-phi_zero_mpb)))*partition_tab_current
                    else
                        freen_est = freen_est - kBT_mpb/a_mpb**3 *&
                            log(1d0+phi_zero_mpb*(alpha_func_mpb(i_full_points)-1d0))*&
                            partition_tab_current
                    end if
                else !a=0 case
                    !old version, same as new in this case
                    !freen_est = freen_est + (1d0-alpha_func_mpb(i_full_points))*&
                        !    2d0*c_mpb*kBT_mpb*partition_tab_current
                    if (use_separate_alpha) then
                        func1 = alpha_func_mpb(i_full_points)!*exp(-z_mpb*pot_old/kBT)
                        func2 = alpha_func_anion_mpb(i_full_points)!*exp(z_mpb*pot_old/kBT)
                        freen_est = freen_est - c_mpb*kBT_mpb*&
                            (func1+func2-2d0)*partition_tab_current
                    else
                        freen_est = freen_est -2d0*c_mpb*kBT_mpb*&
                            (alpha_func_mpb(i_full_points)-1d0)*partition_tab_current
                    end if
                end if
                freen_solv_plus_ion = freen_solv_plus_ion + ion_freeen_1

                ! 	      if (alpha_func_mpb(i_full_points).gt.0d0) then
                !full free energy of full system (solute + ions + solvent) according to MPBE free energy expression
                !both f_function and ln_cosh_function contain multiplicative kappa-function, so we only evaluate this
                !expression if kappa function is non-zero
                if (.not.solve_pbe_only) then
                    ion_freeen_2 = &
                        -kBT_mpb/a_mpb**3*ln_cosh_function_current*partition_tab_current
                else !a=0 case
                    ion_freeen_2 = &
                        -kBT_mpb*ln_cosh_function_current*partition_tab_current
                end if
                freeen_sol_plus_solv_plus_ion = freeen_sol_plus_solv_plus_ion + ion_freeen_2
                ! 	      end if

                !difference (ion_freeen_2 - ion_freeen_1) is the solvation energy from the
                !MPBE free energy expression additional to the normal electrostatic energy 1/2 n v
                freeen_MPBE_nonelstat_solv_energy = freeen_MPBE_nonelstat_solv_energy +&
                    (ion_freeen_2 - ion_freeen_1)
                !to obtain purely nonelstat free energy, we substract freeen_MPBE_elstat (elstat ion-potential interaction)
                !which we already added to hartree_delta_energy
                freeen_MPBE_nonelstat_solv_energy = freeen_MPBE_nonelstat_solv_energy - 2d0 * freeen_MPBE_elstat
                !factor 2 here!!!!
            else if (solve_lpbe_only.and..not.use_mpbe_free_energy) then
                !LPBE Nonelstat Free Energy requested
                if (alpha_func_mpb(i_full_points).gt.0d0) then
                    !only for alpha_func_mpb > 0, freeen_LPBE_elstat is > 0
                    if (.not.solve_debye_only) then
                        ion_freeen_2 = &
                            -kBT_mpb/a_mpb**3*ln_cosh_function_current*partition_tab_current
                    else !a=0 case
                        ion_freeen_2 = &
                            -kBT_mpb*ln_cosh_function_current*partition_tab_current
                    end if
                end if
                freeen_LPBE_nonelstat_solv_energy = freeen_LPBE_nonelstat_solv_energy + ion_freeen_2
                !substract electrostatic energy that is already included in hartree_delta_energy
                freeen_LPBE_nonelstat_solv_energy = freeen_LPBE_nonelstat_solv_energy - freeen_LPBE_elstat
            end if !mpbe or lpbe free energy
        end if ! if .not.mpbe_no_ion


    end subroutine evaluate_mpb_energies

    !******
    !---------------------------------------------------------------------
    !****s* lpb_solver_utilities/evaluate_Gnonmf_forces_dens_mat
    !  NAME
    !   evaluate_Gnonmf_forces_dens_mat
    !  SYNOPSIS

    subroutine evaluate_Gnonmf_forces_dens_mat( &
            matrix_shell, hessian_basis_wave, Gnonmf_gradient_deriv, &
            gradient_basis_wave, wave, n_compute, n_points, &
            i_dim, i_spin, partition_tab)

    !  PURPOSE
    !   The subroutine evaluates gga-force term for severel integration points, using the density matrix formalism.
    !   So this is the gga version of evaluate_hamiltonian_shell_PF.

    !  USES

        use dimensions
        !    use runtime_choices
        implicit none

    !  ARGUMENTS

        real*8 :: matrix_shell( n_compute,n_compute )
        real*8, dimension(n_max_compute_dens, 6, n_points):: hessian_basis_wave
        real*8, dimension(3, n_spin, n_points)                    :: Gnonmf_gradient_deriv
        real*8, dimension(n_max_compute_ham*3, n_points)  :: gradient_basis_wave
        real*8, dimension(n_max_compute_ham, n_points)    :: wave
        integer :: n_compute
        integer :: n_points
        integer :: i_dim
        integer :: i_spin
        real*8  :: partition_tab( n_points)
        ! Required for mGGA evaluation of forces
        real*8, dimension(n_spin, n_points)                       :: Gnonmf_tau_deriv
        logical :: meta_gga_forces_on

    !  INPUTS
    !
    ! o hessian_basis_wave -- hessian of the basis waves
    ! o Gnonmf_gradient_deriv -- Gnonmf gradient terms
    ! o gradient_basis_wave -- gradients of the basis functions
    ! o wave -- the basis functions
    ! o n_compute -- the number of basis functions in the grid patch
    ! o n_points -- the number of grid points in the grid patch
    ! o i_dim -- the dimension where forces are calculated
    ! o partition_tab -- values of the partition function
    !
    !  OUTPUT
    ! o matrix_shell -- the results of the force components are added here.
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


        integer :: i_index(3,3)
        integer :: i_coord_1, i_coord_2, i_basis_1, i_basis_2, i_point, i_counter, i_compute
        real*8  :: matrix_tem1(n_compute, n_points)
        real*8  :: matrix_tem2(n_points, n_compute)

        ! mGGA temporary storage arrays
        real*8  :: matrix_tem1_mgga(n_compute, n_points)
        real*8  :: matrix_tem2_mgga(n_points, n_compute)


        i_index = 0
        i_counter = 0

        do i_coord_1 = 1,3

            i_counter = i_counter +1
            i_index(i_coord_1, i_coord_1) = i_counter

            do i_coord_2 = i_coord_1+1,3

                i_counter = i_counter +1
                i_index(i_coord_1, i_coord_2) = i_counter
            end do
        end do


        do i_coord_1 = 1,3
            do i_coord_2 = 1,3

                i_index(i_coord_1, i_coord_2) = max(i_index(i_coord_1, i_coord_2), i_index(i_coord_2, i_coord_1))

            end do
        end do

        do i_basis_2 = 1, n_compute
            do i_point = 1, n_points
                matrix_tem2(i_point, i_basis_2) = &
                    gradient_basis_wave(i_basis_2+(i_dim-1)*n_compute, i_point) &
                    * partition_tab(i_point)
            end do
        end do

        do i_coord_2 = 1,3

            do i_basis_2 = 1, n_compute
                do i_point = 1, n_points

                    matrix_tem1(i_basis_2, i_point) = &
                        gradient_basis_wave(i_basis_2+(i_coord_2-1)*n_compute, i_point) &
                        * Gnonmf_gradient_deriv(i_coord_2, i_spin, i_point)

                end do
            end do

            call dgemm('N', 'N', n_compute, n_compute,  &
                n_points, 1.d0,  &
                matrix_tem1, n_compute, &
                matrix_tem2, n_points, &
                1.d0, matrix_shell, n_compute )

        end do

        do i_basis_2 = 1, n_compute
            do i_point = 1, n_points
                matrix_tem1(i_basis_2, i_point) = &
                    wave(i_basis_2, i_point)* partition_tab(i_point)
            end do
        end do

        do i_coord_2 = 1,3

            do i_basis_2 = 1, n_compute
                do i_point = 1, n_points

                    matrix_tem2(i_point, i_basis_2) = &
                        hessian_basis_wave(i_basis_2, i_index(i_dim,i_coord_2), i_point) &
                        * Gnonmf_gradient_deriv(i_coord_2, i_spin, i_point)

                end do
            end do

            call dgemm('N', 'N', n_compute, n_compute,  &
                n_points, 1.d0,  &
                matrix_tem1, n_compute, matrix_tem2, &
                n_points, 1.0d0, matrix_shell, n_compute )

        end do

    end subroutine evaluate_Gnonmf_forces_dens_mat

    !****s* lpb_solver_utilities/evaluate_orb_hess_dot_rho_grad_Gnonmf
    !  NAME
    !    evaluate_orb_hess_dot_rho_grad_Gnonmf
    !  SYNOPSIS

    subroutine evaluate_orb_hess_dot_rho_grad_Gnonmf(hessian_basis_wave, n_compute, &
            KS_ev_compute, Gnonmf_gradient_deriv, n_points, max_occ_number, i_atom_2, basis_offset, &
            n_local_compute, orb_hess_dot_rho_grad)


    ! PURPOSE
    ! The subroutine evaluates dot product between KS_orbital_gradients and
    ! Gnonmf_gradient_deriv (including gradient rho respectively)
    ! needed for gga-forces.
    !
    ! AJL: We also compute orb_hess_dot_orb_grad in here for the Meta-GGA forces
    ! as all the terms are constructed i.e. this is the most efficient place
    !
    !  USES

        use dimensions
        use runtime_choices
        implicit none

    !  ARGUMENTS

        integer, intent(in) :: n_points
        integer, intent(in) :: n_compute
        integer, intent(in) :: n_local_compute
        integer, intent(in) :: max_occ_number
        integer, intent(in) :: i_atom_2

        real*8,  dimension(n_max_compute_dens, 6, n_points), intent(in) :: hessian_basis_wave
        real*8,  dimension(max_occ_number, n_compute),       intent(in) :: KS_ev_compute
        real*8,  dimension(3, n_points),                     intent(in) :: Gnonmf_gradient_deriv
        integer, dimension(n_atoms+1),                       intent(in) :: basis_offset
        real*8, dimension(n_states, n_max_batch_size, 3),   intent(out) :: orb_hess_dot_rho_grad


    ! INPUTS
    ! o hessian_basis_wave -- ???????????
    ! o n_compute -- number of relevant basis functions
    ! o KS_ev_compute -- relevant Kohn-Sham eigenvectors
    ! o Gnonmf_gradient_deriv -- gradient of Gnonmf energy
    ! o n_points -- number of grid points
    ! o max_occ_number -- number of occupated states
    ! o i_atom_2 -- atom index
    ! o basis_offset -- starting point of the spline
    ! o n_local_compute -- ??????????/
    !
    ! OUTPUT
    ! o orb_hess_dot_rho_grad -- Hessian matrix times gradient of electron density
    ! o orb_hess_dot_orb_grad -- Hessian matrix times gradient of KS orbitals
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

        ! local variables
        real*8, dimension(n_local_compute,n_points) :: aux_wave_hessian
        real*8, dimension(max_occ_number, n_points) :: nuclear_hessian

        ! counter
        integer :: i_state
        integer :: i_coord_1
        integer :: i_coord_2
        integer :: i_point
        integer :: i_counter
        integer :: i_local_compute

        ! Function
        real*8, external :: ddot

        orb_hess_dot_rho_grad = 0.d0

        i_counter = 0
        do i_coord_1 = 1, 3, 1
            do i_coord_2 = i_coord_1, 3, 1
                i_counter = i_counter + 1
                ! condense and reorganize wave_gradient (this is expensive!)
                ! maybe outside the i_atom_2-loop for all i_compute ?
                do i_point = 1, n_points, 1
                    do i_local_compute = 1, n_local_compute, 1
                        aux_wave_hessian(i_local_compute, i_point) = &
                            hessian_basis_wave(basis_offset(i_atom_2) + &
                            i_local_compute-1,i_counter,i_point)
                    end do
                end do

                call dgemm('N','N', max_occ_number, n_points, n_local_compute, 1.0d0, &
                    KS_ev_compute(1, basis_offset(i_atom_2)), max_occ_number, aux_wave_hessian, &
                    n_local_compute, 0.0d0, nuclear_hessian, max_occ_number)

                ! calculate dot product from wave_hessian and rho_gradient
                do i_point = 1, n_points, 1
                    call daxpy(max_occ_number, Gnonmf_gradient_deriv(i_coord_2, i_point), &
                        nuclear_hessian(1,i_point), 1, &
                        orb_hess_dot_rho_grad(1,i_point,i_coord_1), 1)
                end do

                if (i_coord_1 .ne. i_coord_2) then
                    do i_point = 1, n_points, 1
                        call daxpy(max_occ_number, Gnonmf_gradient_deriv(i_coord_1, i_point), &
                            nuclear_hessian(1,i_point), 1, &
                            orb_hess_dot_rho_grad(1,i_point,i_coord_2), 1)
                    end do
                end if
            end do
        end do

    end subroutine evaluate_orb_hess_dot_rho_grad_Gnonmf

    !******


end module lpb_solver_utilities
