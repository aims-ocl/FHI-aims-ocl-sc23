!****h* FHI-aims/runtime_choices
!  NAME
!    runtime_choices
!  SYNOPSIS


      module runtime_choices


!  PURPOSE
!
!
!  Module runtime_choices collects all switches which determine the runtime
!  behaviour of the code. There should be no allocations necessary.
!
!
!  Variable declarations
!
!     charge   : Total charge of the system
!     hydro_cut : specifies whether hydrogenic fns are introduced without or with cutoff potential
!     spin_treatment : 0 - non-spinpolarized calculation, 1 - spin-polarized, collinear spin-density
!     functional theory
!     flag_rel : 0 - non-relativistic, 1 - scalar-relativistic
!     override_relativity : Logical flag that allows to switch off internal checks that would otherwise
!                stop physical bogus settings - like running non-relativistic gold and publishing
!                a PRL about its novel, unexpected properties. Set this to true ONLY if
!                you know what you are doing.
!     zora_threshold: For ZORA scalar relativity, determines the potential threshold below which
!                the relativistic terms are actually calculated
!                (due to the nature of the ZORA scalar relativistic operator, we have
!                 v_threshold = 4 c^2 / (2 + (1/zora_threshold))
!     flag_xc  : 3 - pz-lda, 30 - dft part in xyg3
!     flag_dftpt2_dft_part : flag about the DFT parts in the final DHDFAs
!     flag_atomic_xc  : 8 - pw-lda, 12 - KLI
!     flag_eg       : .true. - even spacing grids for vdwdf calculation
!     flag_kernel   : .true. - kernel data for LL van der Waals is given
!     flag_cell     : .true. - cell_origin is given
!     n_empty  : number of empty KS states per atom
!     flag_KS  : Treatment of KS equations in s.c.-loop:
!               -1 - Lapack treatment, but uncertain whether or not the basis is singular
!                0 - direct diagonalisation by lapack wrapper routine
!                1 - "simultaneous diagonalisation" - lapack, but first transform the
!                    problem to a non-singular one (i.e. use natural orbitals)
!                    This is more stable than "lapack" and should be the default one day.
!                    FIXME - The present implementation is not numerically efficient, it's
!                    just a patch which I cooked up to help me get around the problem,
!                    without any deep optimization of the underlying numerics. Is there
!                    a better/standard implementation?
!                2 - experimental solution of eigenvalue problem by two successive singular
!                    value decomposition steps, first for the overlap matrix, then for the Hamiltonian.
!                    SVD requires extra workspaces because it cannot deal with packed matrices.
!     basis_threshold : Eigenvectors of the overlap matrix with eigenvalues less than
!                    basis_threshold (i.e. the nullspace of the overlap matrix) are removed
!     prodbas_threshold : Eigenvectors of the Coulomb matrix with eigenvalues less than
!                    prodbas_threshold (i.e. the nullspace of the Coulomb matrix) are removed
!                    if flag_KS=1 . This amounts to the use of natural orbitals.
!     prodbas_nb    : Block size for 1D distribution of auxiliary basis functions
!     partition_acc : Accuracy limit for the partition function in partitioned integrals:
!               Where the norm (denominator!)) of the partition function is less than
!               partition_acc, the remaining density is presumably so small that we
!               do not need to integrate there - simply set partition_tab=0.d0 at such points.
!     wave_threshold : Threshold value for radial functions u(r). If u(r) < wave_threshold
!               for all r > outer_rad, integral contributions outside outer_rad are not evaluated.
!     default_initial_moment : default initial moment of each atom in the structure
!               (spin-polarisation must be used)
!     occupation_type : Calculation scheme for occupation numbers ("smearing")
!     occupation_width : Width parameter for smearing
!     n_methfessel_paxton : Definition of Methfessel-Paxton scheme, if needed.
!     occupation_acc : Accuracy to which occupation numbers are required (total charge inaccuracy)
!     occupation_thr : Threshold value below which occupation numbers are simply set to zero
!     entropy_thr    : Numerical threshold relevant for fermi smearing only
!     fermi_acc      : Required accuracy of Fermi level
!     max_zeroin     : Maximum of numerical iterations to find correct Fermi level using Brent's method
!
!     mixer : specify the type of charge density mixer
!             0 - linear lumberjack-style mixing
!             1 - Pulay mixing
!     ini_linear_mixing: In case of Pulay mixing, number of initial iterations for which
!             a linear mixing algorithm is pursued for good preconditioning
!     linear_mix_param: charge mixing parameters for simple linear mixing of electron & spin density
!     charge_mix_param: charge mixing parameters for mixing of electron & spin density
!     max_rho_change
!     sc_accuracy_rho:  maximum change of electron density allowed to define self-consistency
!     sc_accuracy_eev:  maximum change of sum of eigenvalues allowed to define self-consistency
!     sc_accuracy_etot: maximum change of total energy allowed to define self-consistency
!     sc_accuracy_potjump: maximum change of the potential jump allowed to defined self-consistency
!     flag_acc_rho : check whether convergence accuracy of charge density has been set
!     flag_acc_eev : check whether convergence accuracy of sum of eigenvalues has been set
!     flag_acc_etot : check whether convergence accuracy of total energies has been set
!     flag_acc_potjump: check whether convergence accuracy of potential jump has been set
!     sc_iter_limit : Maximum number of self-consistency cycles
!     sc_init_iter : Initial number of iterations after which the scf cycle (mixing) is
!                    completely reinitialized using only the orbitals of the last iteration done.
!                    For technical reasons, this uses exactly the same logic as a
!                    geometry step. The default is initially set to 1000 but should
!                    reduced to approximately 50 (or some other better value) to
!                    help cconverge difficult cases in general.
!     sc_init_factor : Factor by which the density convergence criterion for scf
!                    (sc_accuracy_rho) must not yet have been met to allow the code
!                    to trigger the scf restart at sc_init_iter iterations.
!     scf_restarted : Remembers whether an initial scf cycle restart was already triggered
!                    using the sc_init_iter flag
!     constraint_it_lim : for locally constrained DFT - max number of iterations to find
!                         correct constraint_potentials / Fermi levels.
!     constraint_mix    : for locally constrained DFT - mixing factor(s)  to determine
!                         correct constraint_potentials / Fermi levels
!     constraint_precision : for locally constrained DFT - accuracy to which the requested
!                         constraint must be fulfilled
!     force_potential : Debug flag to specify effective potentials other than the self-consistent KS potential.
!                       0: (keyword sc) self-consistent Kohn-Sham potential
!                       1: (keyword superpos_pot) superposition of atomic effective potentials
!                       2: (keyword superpos_rho) superposition of atomic densities
!                       3: (keyword non-selfconsistent) Only eigenvalues for given potential
!                          will be calculated, no update of potential and no self-consistency loop.
!
!     force_hartree_log_grid: If false then hartree potential is integrated in radial
!                                         grid, if true logarithmic grid is used.
!     ipi_dip : If true passes dipole to i-PI wrapper
!     ipi_hirshfeld: If true passes hirshfeld charges to i-PI wrapper
!     ipi_ftensor: If true passes friction tensor to i-PI wrapper
!     ipi_work: If true passes converged work function to i-PI wrapper
!     out_basis : If true, all basis data (radial fns, ON radial fns, kinetic and potential data)
!                will be written as xmgrace-readable *.dat files.
!     out_onsite : If true, onsite integrands f(r) for <phi_i|H_atom|phi_i> = int dr f(r)
!                will be written as xmgrace-readable *.dat files.
!                This allows to directly check the radial integration grid quality for
!                all basis functions, especially near the nucleus.
!     out_ovlp_spectrum : Output the eigenvalue spectrum of the overlap matrix, to
!                quantify the conditioning of the full basis
!     out_eigenvec : If true, all KS eigenvectors (basis coeffs of eigenstates)
!                will be written
!     out_overlap  : If true, the overlap matrix will be written (periodic only)
!                will be written
!     out_hamiltonian  : If true, the hamiltonian matrix will be written (periodic only)
!     out_matrices : If true, all matrices (Hamiltonian, Overlap)
!                will be written.
!     out_matrices_elsi : If true, the Hamiltonian and overlap matrices will be written (ELSI format).
!     >>> AB: feb 2012
!     out_aitranss : if true, overlap matrix and KS-eigenvectors will written out in a
!                    format compatible with "aitranss" module (Karlsruhe transport code)
!     <<< done with insert: AB: feb 2012

!     out_nconsistent : write the output after each nth consistent step
!                       out_nconsistent=-1 write never
!     out_niteration : write the output after each nth iteration step
!                       out_niteration=-1 write never
!     out_ninitialize : write the output after each nth initialization step
!                       out_ninitialization=-1 write never
!     out_grids : If true, radial grids will be written to .dat files
!     out_dipole_moment : If true, total dipole moment will be written to stdout
!     out_quadrupole_moment: If true, the eigenvalues of the qudrupole tensor will be written out
!     out_ascii : if true, write _some_ output in ascii else binary
!     out_v_eff : If true, effective integration potential will be written to .dat files
!     out_v_eff_new : If true, effective integration potential will be written in new format
!     out_v_hartree : If true, partitioned Hartree potentials will be written to .dat files
!     out_mulliken : If true, Mulliken charge analysis will be written to file 'Mulliken.out'.
!     out_loewdin : If true, Loewdin charge analysis will be written to file 'Mulliken.out'.
!     out_dos : If true, projected density of states will be written to file '*dos.dat'.
!     out_dos_tetrahedron : If true, projected density of states will be written to file '*dos_tetrahedron.dat'.
!     out_polarisability : If true, MP2/RPA dipole polarisability will be
!         written to files "mp2_polarisability.dat" and "rpa_polarisability.dat".
!     safe_minimum : lapack-determined safe minimum real*8 number, so that
!                1.0d0/safe_minimum does not overflow .
!     n_min_points_in_division: Minimum number of points in each angular division
!     angular_div_cutoff: sets the minimum radial cutoff for angular shells to be handled
!     in one single batch in the integrations instead of dividing them
!     multipole_interpolation_style : interpolation style for multipole density
!                                     0 - cubic spline, no use of derivatives
!                                     1 - hermite spline, use derivatives
!     outer_radius_of_multipole : Cutoff radius outside of which the numerically tabulated Hartree
!                                 multipole potential
!               is replaced by an analytical field of fixed constant multipole moment
!
!
!      n_k_points_xyz(3) : number of k-points in k_x, k_y and k_z direction     (periodic systems)
!      k_points_offset(3): ofset of k-point phase in k_x, k_y and k_z direction (periodic systems)
!      read_k_points_from_file : read k-points from the file k_list.in          (periodic systems)
!      real_eigenvalues : is the eigenvectors real or complex ? In clusters real, periodic systems
!                         with k-points can also be complex.
!
!
!
!      far_distance_hartree_multipole_radius_threshold: electron density where the multipole_radius
!                                  is set, (all the charge is inside multipole radius)
!
!      far_distance_hartree_multipole_moment_radius_threshold : electron density threshold for multipole
!                      moments calculations.
!
!        far_distance_hartree_multipole_moment_threshold: if all multipole moment of l, is smaller than
!                  this in the  multipole radius distance, the set if l:s are not used in far distance
!                  hartree calculations.
!
!        distribute_leftover_charge:  distribute leftover charge from multipole expansion
!                                     equally among all atoms to ensure a
!                                     charge-neutral system
!
!      multipole_radius_free_threshold : Threshold that determined the outermost radius of the charge density of
!                                        a free atom. The radius is set to the outermost logarithmic grid shell at
!                                        which the free-atom charge density is still smaller than the chosen threshold.
!
!      use_hartree_non_periodic_ewald: Switch that determines whether to use Ewald's method
!                                      for the Hartree potential in the non-periodic case
!
!      relax_mode : Switch that determines whether or not we are performing geometry relaxation along the way
!      relax_geo : Relaxation scheme that is actually on at a given time
!      relax_accuracy_forces : desired convergence accuracy for the maximum force after geometry relaxation.
!      store_eigenvectors_to_disk_in_relaxation: in order to save memory the eigenvectors are stored to disck.
!      prerelax_accuracy_forces : for a hybrid relaxation scheme (first md_descent, then bfgs), this flag
!           determines when w switch from the first-order to the second-order scheme
!      remove_unitary_force_components : Clean calculated forces from all pure translations / rotations
!                  (unitary transformations)
!      energy_tolerance : value by which the energy may maximally increase in a single relaxation step
!      aggregated_energy_tolerance : value by which the energy may maximally increase across the entire
!                                    structure relaxation, measured beginning from the overall lowest
!                                    energy reached. This is an overall cap to prevent relaxation trajectories
!                                    from increasing over multiple steps without warning.
!      harmonic_length_scale : assume harmonic potential energy surface on length scales smaller than this
!      max_atomic_move  : maximum permissible displacement of a single atom
!      min_line_step : Minimum line step to find decrease in energy before throwing out H
!      line_search_tol : Tolerance for accepting BFGS-suggested line steps when comparing to a-posteriori-determined step
!      line_step_reduce : factor by which the trial line step will be reduced if a BFGS move was rejected
!      line_step_reduce_automatic : use line search instead of simple, stupid reduction. Default: .true.
!      bfgs_mixing_factor : line step mixing in relaxation run to take care of under (over) estimated line steps
!      bfgs_mixing_cap : maximally mixed in line step estimate
!      bfgs_extrapolation_cap : no line step longer than this, deemed to be garbage otherwise
!      bfgs_extrapolation_on : .true. - use mixing algorithm, .false. use usual quadratic update
!      max_relaxation_steps : maximum number of steps in geometry relaxation
!      init_hess : initial Hessian is chosen as diagonal matrix times this factor. Default: 97.173615 eV/A^2 == 1 Ha/a^2
!      hybrid_coeff : the portion of exact-change contribution the hybrid scheme.
!
!      dftpt2_Ex_hf       : the portion of HF-like exchange contribution in double-hybrid density functionals
!      dftpt2_Ec_osPT2    : the portion of osMP2-like correlation contribution in double-hybrid density functionals
!      dftpt2_Ec_ssPT2    : the portion of ssMP2-like correlation contribution in double-hybrid density functionals
!      dftpt2_Ec_oslrcPT2 : the osMP2-like long-range correlation
!      dftpt2_Ec_sslrcPT2 : the ssMP2-like long-range correlation
!      lrc_pt2_started    : Remembers whether a lrc-pt2 evaluation was already triggered
!
!      flag_auxil_basis: PRODBAS_OPT, PRODBAS_FULL, or PRODBAS_SVD,
!                        determines which auxiliary basis will be used
!                        in RI (resolution of identity).
!      hf_version: either HF_DM or HF_EIGEN; defines whether to use an
!                  intermediate density matrix or to use directly the
!                  eigencoefficients to calculate the exchange matrix
!      RI_type: one of RI_SVS, RI_V, and RI_LVL; defines the choice of
!               expansion coefficients of products in the auxiliary
!               basis.
!      use_logsbt: If .true., use the 1D logSBT integration scheme for
!                  2-center auxiliary integegrals instead of the 3D
!                  grid code.
!      sbtgrid_{lnr0,lnk0,lnrange,N}: Parameters for the logarithmic grid
!                  used in the logSBT code.  See logsbt.f90.
!      use_2d_corr: Use efficient 2D distribution for beyond GGA correlation
!                   code where implemented (def: use_scalapack)
!                   WARNING: Just because use_2d_corr is .true. does not
!                            guarantee that a given ovlp_3KS is 2D distributed.
!                            That is in the responsibility of the constructor.
!      walltime_requested : if true, there is a wall time limit on the job (logical)
!      walltime_limit     : if (walltime_requested), this contains the maximal number of seconds in the job; default = 5 years
!
!      precondition_kerker_q0  : damping wave vector for Kerker preconditioning
!      precondition_max_l      : maximum angular momentum value for multipole expansion in Kerker preconditioning


!      packed_matrix_threshold : Only in packed matrix format PM_index.
!                                If a value of the overlap-matrix is absolutelu smaller than this after initialize_integrals
!                                for it is not saved the memory place in the overlap-matrix or the hamiltonian_matrix.
!      legacy_monopole_extrapolation : If .true., use the old monopole extrapolation in update_hartree_potential_p1().  For .false., use an improved one.
!      multiplicity :     multiplicity in the unrestricted Hartree-Fock calculations.
!      dos_n_en_points : DOS calculation, number of points
!      dos_low_energy  : DOS calculation, energy start
!      dos_high_energy : DOS calculation, energy end
!      dos_alpha       : DOS calculation, Gaussian broadening parameter
!      calculate_atom_bsse: If .true., calculate the atomization BSSE of the molecule in question
!      n_frozen_core   : if (n_frozen_core .eq. 0) no frozen core shell
!                        if (n_frozen_core .gt. 0) n_frozen_core represents the number
!                      of the valence shells under which the shells are frozen
!                        if (n_frozen_core .lt. 0) n_frozen_core represents the number
!                      of the frozen core shells starting from the lowest shell (n=0)
!      tau_threshold : cutoff for tau to be calculated in meta-GGA routines
!      ELPA_AEO : keyword for ELPA-AEO specific details, so far development only (integer)
!            Introduce keyword ELPA_AEO int, for now 4 is used for single and 0 as default
!            double precision. For now the keyword is for development purposes only!!
!  USES

  use types, only: dp
  use constants, only: hartree, bohr
  use mpe_constants, only: MPE_CONST, ISC_CONST
  use mbd, only: mbd_input_t

  implicit none

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


      ! allows to only check control.in and geometry.in for consistency but do no work-intensive parts.
      logical :: dry_run

      real*8 charge
      logical hydro_cut
      integer spin_treatment

      integer flag_rel
      logical :: override_relativity
      integer, parameter :: REL_not_set     = -1
      integer, parameter :: REL_none        = 0
      integer, parameter :: REL_zora        = 1
      integer, parameter :: REL_atomic_zora = 2
      integer, parameter :: REL_own = 3
      integer, parameter :: REL_KOLNING_HARMON = 4
      integer, parameter :: REL_zora_spinor = 5
      integer, parameter :: REL_at_zora_spinor = 6
      integer, parameter :: REL_x2c = 7
      integer, parameter :: REL_4c_dks = 8
      logical:: flag_scaled_zora
      logical:: flag_KH_post_correction
      logical:: flag_KH_core_states
      logical:: use_small_component

      logical:: flag_eg
      logical:: flag_kernel
      logical:: flag_cell

      !flags to override some safetey statements in FHI-aims
      logical :: override_warning_nobasis
      logical :: override_warning_negative_nucleus
      logical :: override_warning_libxc

      real*8 zora_threshold
      integer :: flag_xc
      integer :: flag_xc_pre
      real*8  :: flag_old_xc
      integer :: n_steps_xc_pre

! Variables needed for libxc
!      integer :: flag_libxc_x_id
!      integer :: flag_libxc_c_id

! Redundant variable? AJL
      integer flag_xc2
      integer flag_fxc
      integer flag_atomic_xc
      character :: hse_unit
      real*8 :: hse_omega
      real*8 :: hse_omega_hf
      real*8 :: hse_omega_pbe
      ! This lc_dielectric_constant is currently only used in the lc_wpbeh functional as an optional parameter.
      ! It cannot be less than 1 by definition. The asymptotic decay of the exchange is \frac{1}{\epsilon*r}.
      real*8 :: lc_dielectric_constant

      ! The Coulomb operator will be scaled
      ! to 0.5 at cutCb_rcut and
      ! to 0.5*erfc(+/-1.) at cutCb_rcut * cutCb_width**(+/- 1.).
      real*8 :: cutCb_width   ! Cutting width for Cb operator.
      real*8 :: cutCb_rcut    ! Cutting pos for Cb op (8 AA).
      ! The cutCb parameters can either be set explicitly or by a scaling
      ! of default values which depend on other settings.
      real*8 :: cutCb_width_factor ! erfc-width across 8 points.
      real*8 :: cutCb_rcut_factor  ! Inner radius of Wigner-Seitz cell.
      integer, parameter :: CUTCB_FROM_DEFAULT = 0
      integer, parameter :: CUTCB_FROM_EXPLICIT = 1
      integer, parameter :: CUTCB_FROM_FACTOR = 2
      integer :: flag_cutCb_width
      integer :: flag_cutCb_rcut

      ! The Coulomb operator will be attenuated by explicit range separation
      ! ERS, erfc((r_{12}-ers_width)*ers_width) Igor
      real*8 :: ers_width
      real*8 :: ers_omega
      real*8 :: ers_rcut
      real*8 :: erfc_omega
      ! for standard DH-DFAs
      integer flag_dftpt2_dft_part
      integer :: flag_dftpt2_xc
      integer :: flag_post_xc
      real*8 :: dftpt2_Ex_hf
      real*8 :: dftpt2_Ec_osPT2
      real*8 :: dftpt2_Ec_ssPT2
      real*8 :: dftpt2_Ec_oslrcPT2
      real*8 :: dftpt2_Ec_sslrcPT2
      ! for lrc-pt2 in DH-DFAs
      character :: lrc_pt2_unit
      real*8 :: lrc_pt2_omega
      logical :: lrc_pt2_started

      integer :: flag_hartree_partition_type
      integer :: partition_type
      real*8,dimension(2):: hartree_partition_parameters
      real*8,dimension(2):: integral_partition_parameters
      real*8 :: stratmann_a ! default to literature value
      integer :: n_empty    ! -1 now signals to read_geo.f90 that no value has been set in control.in
      logical :: override_illconditioning
      real*8 basis_threshold
      real*8 prodbas_threshold
      integer :: prodbas_nb
      integer, parameter :: FLAG_BASIS_ONLY_PRODS = 0
      integer, parameter :: FLAG_BASIS_INCL_ONEPARTICLE = 1
      integer :: flag_basis
      logical :: use_smooth_prodbas_threshold
      logical :: use_asym_inv_sqrt
      real*8 :: condition_penalty
      real*8 partition_acc
      real*8 wave_threshold
      real*8 default_initial_moment
      logical::  use_aufbau_principle !Adjust atomic occupation numbers to aufbau principle
      logical:: use_density_matrix ! Charge density by dm (or wf)?
      integer, parameter :: UPDATE_ORBITAL = 1 ! switches for density update method
      integer, parameter :: UPDATE_DENS_MATRIX = 2
      integer, parameter :: UPDATE_SPLIT = 3 ! forces will be updated via density matrix and electron density via orbitals
      integer, parameter :: UPDATE_AUTO = 4
      integer :: density_update_method ! default is automatic selection of density update
      logical:: use_density_matrix_hf  ! Used to check for dm mixing.
      logical :: autoselect_density_method
      logical :: split_updates ! if .true., forces will be updated via density matrix and electron density via orbitals
      logical :: lda_update_switch ! if .true., dm is used if no forces are needed, otherwise update via orbitals
      logical :: measure_forces ! measure the time for different force updates and select the faster one?
      logical:: use_ovlp_swap
      logical:: flag_compute_kinetic

      ! The next flag is intended to be a temporary flag which decides how the
      ! average electrostatic potential is calculated for periodic boundary conditions.
      ! There are two ways to do this: Either by simple integration on the
      ! integration grid (always works but may miss some parts of the unit cell volume
      ! for very sparse structures / integration grids), or by an analytic expression
      ! calculated for free atoms on a radial grid.
      logical :: analytic_potential_average

      integer occupation_type
      real*8 occupation_width
      real*8 occupation_acc
      real*8 occupation_thr
      real*8 :: entropy_thr ! VB: No longer configurable in control.in
      real*8 fermi_acc
      integer :: mu_method ! 0: zeroin
                           ! 1: bisection (ELSI)
      integer :: max_zeroin
      integer n_methfessel_paxton

      integer :: atomic_solver ! The radial atomic solver used to generate basis functions
                               ! Choices:
                               ! 0 - sratom (default)
                               ! 1 - atom_sphere
      integer, parameter :: ATOMIC_SOLVER_SRATOM = 0
      integer, parameter :: ATOMIC_SOLVER_ATOM_SPHERE = 1
      ! WPH: Deceptively, aims by default does not use the basis functions coming out
      ! of get_free_atom as the minimal basis set; instead it recalculates them in the
      ! get_species_basis_fn() by integrating free_potential (which was obtained from
      ! get_free_atom) using dftseq
      ! While this works fine for LDAs and GGAs, for basis functions incorporating
      ! exact exchange, this poses difficulties, as aims assumes free_potential has
      ! no l-channel dependence.  While we should implement this in the future (the
      ! official motto of FHI-aims,) for now we bypass the recalculation of the
      ! minimal basis entirely when needed and reuse the wavefunctions and
      ! derivative quantities (eigenvalues, kinetic term, and possibly wavefunction
      ! gradients) that came from the atomic solver
      logical :: use_atomic_solver_fns_as_minimal_basis

      integer mixer
      integer, parameter :: MIX_LINEAR = 0
      integer, parameter :: MIX_PULAY  = 1
      integer, parameter :: MIX_BROYDEN  = 2
      integer mixer_constraint
      integer ini_linear_mixing
      integer ini_linear_mixing_constraint
      real*8, dimension(:), allocatable :: linear_mix_param
      real*8, dimension(:), allocatable :: charge_mix_param
      ! special Pulay mixing parameter while Pulay mixer is on
      real*8, dimension(:), allocatable :: prec_mix_param
      real*8, dimension(:), allocatable :: max_rho_change
      real*8 :: multipole_feed_back_parameter
      logical:: flag_delta_rho_in_multipole
!      logical:: flag_delta_rho_in_multipole =.true.

      ! Flag that determines whether we will try to automatically
      ! adjust our scf mixer settings once the scf cycle has started
      ! and we know something about the system we're studying.
      logical :: adjust_scf
      logical :: adjust_scf_always
      integer :: adjust_scf_iteration
      logical :: charge_mix_param_set
      logical :: occupation_width_set

      real*8 sc_accuracy_rho
      real*8 sc_accuracy_eev
      real*8 sc_accuracy_etot
      real*8 sc_accuracy_forces
      real*8 sc_accuracy_stress
      real*8 sc_accuracy_potjump
      real*8 :: coeff_thr ! no longer configurable in control.in - VB

      logical :: flag_acc_rho
      logical :: flag_acc_eev
      logical :: flag_acc_etot
      logical :: flag_acc_potjump
      logical :: flag_moment        ! If false, ALL initial spin moments are set by Hund's rules. (High spin!)

      integer sc_iter_limit

      integer :: sc_init_iter     ! dummy initialization - actually initialized in read_control.f90
      real*8  :: sc_init_factor   ! dummy initialization - actually initialized in read_control.f90
      logical :: scf_restarted    ! dummy initialization - actually initialized in read_control.f90
      logical :: restarting_scf ! This flag memorizes if we are momentarily restarting the scf cycle
                                ! only to reinitialize the mixing, without any geometry move.

      logical :: check_sc_abandon_etot     ! s.c.f. cycle will be abandoned and program stopped if:
      integer :: sc_abandon_etot_iter      !        for sc_abandon_etot_iter iterations
      real*8  :: sc_abandon_etot_threshold !        the total energy change between cycles has been larger than 1000 eV/atom

      ! Threshold for convergence of real-space supercell summation for Tkatchenko-Scheffler
      ! dispersion correction in periodic systems. If no other value is specified explicitly
      ! in control.in, the value below will be multiplied by n_atoms in read_control.f90 . In read_control.f90, a lower floor of 1.d-6/hartree is also
      ! imposed for the total energy convergence (all atoms) unless an explicit choice is made in control.in .
      real*8 :: vdw_convergence_threshold

      integer constraint_it_lim
      real*8,dimension(2):: constraint_mix
      real*8 constraint_precision

      logical :: flag_frozen_core
      integer :: i_start_mp2

      ! flags in case spin-component scaled (scs) mp2 is used
      logical ::use_scsmp2
      real*8 :: ps_mp2
      real*8 :: pt_mp2

!     Atomic reference symmetry breaking
      character*40 :: ext_safe

      integer :: force_potential
      integer :: force_lebedev
      logical :: force_hartree_log_grid
      logical :: force_smooth_cutoff
      real*8  :: smooth_cutoff_limit

      ! compensate_multipole_errors, if true, normalizes the initial density
      ! and the self-consistent model density to keep the exact charge norm on
      ! the 3D, overlapping atom-centered integration grid.
      !
      ! REMINDER: normalize_initial_density only normalizes the electron density used
      ! for initializing the s.c.f. cycle and should be an exact subset (!) of
      ! compensate_multipole_errors .
      logical :: compensate_multipole_errors
      logical :: normalize_initial_density

      integer :: packed_matrix_format
      integer,parameter :: PM_none    = 0
      integer,parameter :: PM_index   = 2

!     parameters used to output real-space electrostatic potential
      integer :: step_x
      integer :: step_y
      integer :: step_z

      real*8  :: vacuum_length
!     Postprocessing: Output Momentum matrix and linear dielectrix tensor or Kubo-Greenwood

      real*8:: Emin
      real*8:: Emax
      integer:: k_point_mom

!      integer :: n_state_min, n_state_max

      real*8:: widthone
      real*8:: widthtwo
      real*8:: omega_min
      real*8:: omega_max
      integer:: n_omega
      integer:: n_omega_reset
      character*1:: ep1
      character*1:: ep2
      logical :: use_gauss
      logical :: read_q_points
      integer :: k_k_method
!     Global flags associated with the Kubo-Greenwood formalism
      integer:: greenwood_method ! defaults to full transition-element matrix treatment
      integer:: percent_or_ev ! defaults to percent of gap-size above VBM option
      integer:: n_greenenergy
      integer::fermistatesabove
      integer::fermistatesbelow
      real*8::dist_from_vbm
      real*8::fermispacing

      logical :: out_fermisurface

      logical :: out_mommat  ! determines if the momentum matrix will be computed and written

!     ESP charge fitting
      real*8:: esp_Emin
      real*8:: esp_Emax
      integer:: esp_k_point
      real*8:: esp_min_radius
      real*8:: esp_max_radius
      integer:: esp_n_max_radius
      integer:: esp_grid_type
      logical :: flag_esp_restraint
      integer:: esp_restr_reference
      real*8:: esp_restraint_strength

      type(mbd_input_t) :: mbd_input
      real*8 :: mbd_beta

!<<<VVG>>>
      real*8 :: mbd_cfdm_dip_cutoff
      real*8 :: mbd_scs_dip_cutoff
      real*8 :: mbd_supercell_cutoff
      logical,dimension(3) :: mbd_scs_vacuum_axis
!<<<VVG>>>


   !############################################################################
   !#
   !#
   !#    The following block contains output related options.
   !#
   !#
      logical :: out_basis
      logical :: out_onsite
      logical :: out_ovlp_spectrum
      logical :: out_eigenvec
      logical :: out_overlap
      logical :: out_hamiltonian

      logical :: flag_found_unfold
      logical :: flag_hdf5_unfold
      logical :: flag_ASCII_unfold

      ! Elpa timings are written to standard output if the following variable is
      ! set to .true. This is achieved by adding the line "output elpa_timings"
      ! to control.in. Note that this variable is actually not defined here but
      ! in elpa1.f90 because is belongs to the Elpa library.
      ! logical :: elpa_print_times
      logical :: wantDebug

      logical :: out_matrices, out_matrices_format_2005
      logical :: out_matrices_elsi, out_h_elsi, out_s_elsi

      logical :: out_xml
      character(len=100) :: xml_file

!     >>> AB: feb 2012
      logical out_aitranss
!     <<< done with insert: AB, feb 2012

      ! The following variable controls the output of certain matrices in
      ! parallel runs. This variable (and the following) is set by
      ! control.in-option "output matrices_parallel types [format]".
      character(3) :: out_mat_par  ! 'n': No output
                                   ! 'h': Hamiltonian matrix
                                   ! 'o': Overlap matrix
                                   ! 's': System matrix for eigenvalues
                                   ! The last three options can be
                                   ! combined. For example, 'ho' means
                                   ! Hamiltonian and overlap.
      ! Format of the respective output:
      character(3) :: out_mat_par_format ! 'asc': ASCII
                                         ! 'bin': Binary
      character(4) :: coulomb_integral_format ! 'full': all i, j, k, l
                                              ! 'qtmp': quantum pack
                                              ! format
      logical :: out_nuc_matrix
      logical :: out_t_plus_v_matrix
      logical :: out_nuc_nuc_repulsion
      logical out_grids
      logical out_v_eff
      logical out_v_eff_new
      logical out_dielec_func
      logical out_ion_dens
      logical out_delta_v
      logical out_ascii
      integer :: out_nconsistent
      integer :: out_niteration
      integer :: out_ninitialize
      logical out_v_hartree
      logical out_rho_multipole
      logical switch_rho_multipole
      logical out_density
      logical out_band
      logical :: out_band_mulliken   ! output the mulliken decomposition of any k point in band
      logical :: out_band_during_scf ! This option will write only the available band structure
                                     ! information during the s.c.f. cycle
      integer :: exx_band_structure_version ! The version of band structure with exact exchange
                                            ! If not set to >0 in control.in in case out_band is true, the program will stop
      ! CC: Z2 by CC / Adapted from Cmera
      logical out_z2_invariant       ! Z2 invariant.
      logical out_mulliken
      logical out_lowdin
      logical :: out_l_proj_dos
      logical :: out_l_proj_dos_tetrahedron
      logical :: out_dos
      logical :: out_dos_tetrahedron
      logical :: out_atom_dos
      logical :: out_atom_dos_tetrahedron
      logical :: out_species_dos
      logical :: out_species_dos_tetrahedron
      logical out_hirshfeld
      logical out_hirshfeld_iterative
      logical :: out_hirshfeld_always ! default set again in read_control.f90
      logical out_vdwdf
      logical out_polarisability
      logical :: out_hessian
      logical :: out_esp
      logical :: out_esp_full

      logical:: flag_run_mulliken
      logical:: flag_run_lowdin

      logical:: flag_gOpenMol_format

      ! VB: default for k eigenvalue output changed to 1
      integer:: out_k_points_eigenvalues
      logical :: out_soc_eigenvalues

      logical :: out_k_point_list

      ! MR: adding list of flags for more flexible interface with i-PI
      logical :: ipi_dip
      logical :: ipi_hirshfeld
      logical :: ipi_ftensor
      logical :: ipi_work
      logical :: ipi_atomic_stress

      logical out_dipole_moment
      logical out_quadrupole_moment
      logical :: out_sorted_eigenvalue_list
      logical out_cube
      logical out_cube_soc
      logical :: overwrite_existing_cube_files !If true, setup_cube_defaults checks whether cube file exists and changes its filename
      integer*8 :: cube_default_size_safeguard !If true, setup_cube_defaults checks whether cube file exists and changes its filename
      logical :: out_cube_nth_iteration !Allows output for cubes during scf
      integer :: cube_output_nth_iteration
      logical :: use_old_cube_routine ! !Backwards compatibility, disabled by default
      character*10 :: cube_content_unit !Backwards compatibiity, enabled by default. Causes output of cube files in the wrong units (1/A^3 instead of 1/bohr^3)
      real*8 ::average_potential

      ! CC: Output flag for DGrid interface:
      logical :: out_dgrid

   !#
   !#
   !#    End of output related options
   !#
   !#
   !############################################################################

   !#  adding potential well
      logical :: potential_well_requested !not by default
      real*8 :: potential_well_start
      real*8 :: potential_well_end
      real*8 :: potential_well_depth

      real*8 safe_minimum

      integer :: n_min_points_in_division
      real*8 :: angular_div_cutoff
      integer :: grid_partitioning_method
      integer :: grouping_factor
      logical :: use_tetgen
      logical :: use_qhull
      logical :: use_nearest_neighbours
      integer :: n_nearest_neighbours
      integer :: min_batch_size
      integer :: batch_size_hard_limit
      logical :: use_hashed_maxmin
      logical :: use_integer_hash
      logical :: out_batch_statistics

      ! RJ: multipole_interpolation_style is always 0, the hermite style has been removed!
      integer, parameter :: multipole_interpolation_style = 0



   !############################################################################
   !#
   !#    The following block of variables contains (roughly) all settings
   !#    concerning the (matrix based) eigenvalue solution of the Kohn-Sham
   !#    problem.
   !#
   !#    Specifically important variables are:
   !#
   !#    use_scalapack:  Determines whether we use parallel linear algebra, or not.
   !#                    This flag also applies to linear algebra done in many other
   !#                    places (Hartree-Fock et al.)
   !#    use_elpa:       Determines whether we use the custom-developed, rather scalable
   !#                    eigensolver libraries elpa1 and elpa2 developed in the ELPA
   !#                    project ("Eigensolvers for Petaflop Applications"), or an
   !#                    only slightly adjusted version which is completely based on the
   !#                    slower infrastructure provided within scalapack itself.
   !#
      integer flag_KS
      integer, dimension(:), allocatable :: flag_KS_k_points
      logical :: flag_force_complex_eigenvectors !! CC: Force eigenvectors to be complex!

      logical :: use_scalapack
      logical :: use_elpa
      logical :: use_scalapack_mbd
      logical :: use_elpa_mbd

      logical :: use_scalapack_DFPT_phonon

      ! details for the ELPA eigensolver
      ! solver_method decides whether to use 1-stage or 2-stage solver
      ! solver_method = 1: 1 stage solver
      ! solver_method = 2: 2 stage solver
      ! solver_method = 0: decide during runtime which one is faster
      integer :: solver_method
      integer :: ELPA_AEO

      logical :: collect_eigenvectors
      logical :: orthonormalize_evs

      logical :: use_lapack_fast

      logical :: use_cg
      integer, parameter :: prec_diagonal = 1
      integer, parameter :: prec_inv_ovlp = 2
      integer, parameter :: prec_inv_ham_ovlp = 3
      integer :: cg_preconditioner
      integer :: initial_ev_solutions
      integer :: max_cg_iterations
      integer :: lopcg_block_size
      real*8 :: lopcg_tol
      real*8 :: lopcg_start_tol
      logical :: lopcg_adaptive_tolerance
      logical :: lopcg_auto_blocksize
      real*8 :: lopcg_omega
      logical :: lopcg_use_all_evs
      logical :: lopcg_slide_ev_block
      integer :: lopcg_skip_rate
      logical :: use_ritz_scalapack
   !#
   !#
   !#    End of eigensolver-related variable declarations
   !#
   !#
   !############################################################################



      integer,  dimension(3)::  n_k_points_xyz
      integer,  dimension(3) :: n_k_points_xyz_nosym
      integer,  dimension(3)::  n_k_points_dos_xyz
      integer,  dimension(3)::  k_enhancement_factor
      real*8,   dimension(3)::  k_points_offset
      real*8 :: symmetry_thresh
      logical:: use_symmetry_reduced_k_grid
      logical:: use_symmetry_reduced_spg
      logical:: use_spg_mv_mm
      logical:: use_spg_full_Delta
      logical:: reconstruct_proper_only
      logical:: get_full_density_first
      logical:: use_k_phase_exx
      logical:: read_k_points_from_file
      logical:: pert_dos_on
      logical:: postscf_eigenvalues
      logical:: real_eigenvectors
      logical:: numerical_stress_save_scf
      logical:: use_analytical_stress
      logical:: use_inv_symmetry
      logical:: use_full_symmetry
      logical:: flag_k_fine_grid
      ! Number of analytical stress components that are calculated in forces_densmat
      ! 6 means that we calculate only the upper triangle and copy the total result to the lower part in the end.
      integer:: AS_components
      logical:: AS_symmetric_output
      ! Settings for the tag stress_for_relaxation
      integer:: stress_for_relaxation
      integer, parameter :: RELAX_NUMERICAL  = 0
      integer, parameter :: RELAX_ANALYTICAL = 1
      ! Flag for external applied pressure
      logical:: flag_external_pressure

      real*8:: far_distance_hartree_multipole_radius_threshold
      real*8:: far_distance_hartree_multipole_moment_radius_threshold
      real*8:: far_distance_hartree_multipole_moment_threshold ! This value must be set more tightly for the default stratmann partition table later!
      real*8:: far_distance_hartree_radius_threshold
      real*8:: far_distance_adaptive_hartree_radius_threshold
      real*8:: multipole_radius_free_threshold  ! 1.d-20
      real*8:: Ewald_radius             ! VB: Used to be 5.0d0
      logical :: Ewald_radius_automatic ! Determines whether Ewald_radius
                                        ! is fixed or determined adaptively.

      ! VB: Temporary flag to allow for testing and refinement of the new
      !     faster version of the update_hartree_potential_recip routine.
      !     This is now the default. Flag kept for now in case problems
      !     creep up - in that case, just switch back and forth within same
      !     binary to find potential problems.
      logical:: fast_Ewald

      logical :: hartree_fp_function_splines

      logical :: distribute_leftover_charge
      ! VB: this will be the threshold in integrate_hartree_log_grid
      !     by which the outer valid radius of each angular momentum
      !     component of every atom in the structure is determined.
      real*8 :: multipole_threshold
      ! VB: The long-range Hartree potential can produce numerical noise for very high l components,
      !     as observed by Paula long ago - essentially, spuriously small density components at
      !     large radii get weighted up in the calculation of the associated multipole moment,
      !     due to factors r^(l+1) in the calculation. The solution are either denser radial grids,
      !     or a simple cutoff of the problematic components, as done here.
      !     The cutoff can be increased by hand in control.in, but in that case it is hopefully clear that
      !     the user knows what they are doing.
      integer :: l_hartree_far_distance
      logical :: use_hartree_non_periodic_ewald
      integer :: hartree_d_matrix_method
      real*8 :: gridwidth_in_pm

      integer:: n_k_points_group
      real*8:: hartree_potential_fourier_part_threshold
      real*8:: extra_adding_to_hartree_potential_distance
      real*8:: packed_matrix_threshold
      logical:: Adams_Moulton_integrator
      logical :: legacy_monopole_extrapolation
      integer :: relax_mode
      integer :: relax_geo
      integer, parameter :: RELAX_OFF = 0
      integer, parameter :: RELAX_BFGS_TB = 2
      integer, parameter :: RELAX_TRM = 3
      integer, parameter :: RELAX_LTRM = 4
      logical            :: RELAX_LTRM_dev
      real*8 ::  relax_accuracy_forces
      ! FlK: The idea here was to use stress as convergence criterion for very
      ! tight convergence settings during relaxtion. This is, however, currently
      ! not used:
      real*8 ::  relax_accuracy_stress_factor
      logical :: use_stress_to_check_relaxation
      real*8 ::  prerelax_accuracy_forces
      integer :: remove_unitary_force_components
      logical :: store_eigenvectors_to_disk_in_relaxation

      real*8  :: energy_tolerance
      real*8  :: aggregated_energy_tolerance
      real*8  :: aggregated_energy_tolerance_per_atom
      real*8  :: harmonic_length_scale
      real*8  :: max_atomic_move
      real*8  :: min_line_step
      real*8  :: line_step_reduce
      real*8  :: line_search_tol
      logical :: line_step_reduce_automatic
      integer :: max_relaxation_steps
      real*8  :: bfgs_mixing_factor
      real*8  :: bfgs_mixing_cap
      real*8  :: bfgs_extrapolation_cap
      logical :: bfgs_extrapolation_on

      ! Parameters that control the shape of the initial Hessian matrix
      ! used for the structure optimization
      ! Note that most defaults are reset within read_control.f90 unless specified in control.in .
      integer, parameter :: HESS_DIAG = 1
      integer, parameter :: HESS_LINDH = 2
      integer, parameter :: HESS_EXPL = 3
      integer, parameter :: HESS_LATTICE = 4
      ! Type of the initial Hessian (see read_control.f90) :
      integer :: init_hess_type
      ! Diagonal element value for all atomic coordinates in an initially diagonal Hessian:
      real*8 :: init_hess_diag
      ! Diagonal element value for any lattice vector coordinates in an initially diagonal Hessian:
      real*8 :: init_hess_lv_diag
      ! We add a small additional diagonal element between all atomic coordinates to the original Lindh matrix:
      real*8 :: init_hess_lindh_diag
      ! ... and there is also a thresholding parameter to determine the Lindh matrix elements that are actually considered:
      real*8 :: init_hess_lindh_thres

      integer, parameter :: PRODBAS_OPT  = 0
      integer, parameter :: PRODBAS_FULL = 1
      integer, parameter :: PRODBAS_SVD  = 2

      integer :: flag_auxil_basis
      real*8 :: hybrid_coeff
      integer :: hf_multiplicity
      integer, parameter :: HF_DM = 0
      integer, parameter :: HF_EIGEN = 1
      integer :: hf_version
      integer, parameter :: RI_SVS = 1
      integer, parameter :: RI_V = 2
      integer, parameter :: RI_LVL = 3
      integer, parameter :: RI_LVL_2nd = 4
      integer, parameter :: RI_LVL_full = 5
      integer :: RI_type
      logical :: sparse_o3fn ! .eqv. (RI_type == RI_LVL)
      ! The following is in dimensions.f90
      ! logical :: use_lvl = .false.      ! .eqv. (RI_type == RI_LVL .or.
      !                                          RI_type == RI_LVL_full .or.)
      !                                          RI_type == RI_LVL_2nd

      ! Sensible choices for water with 'light' settings
      ! logical :: use_logsbt = .true. - declared in dimensions.f90 !
      real*8 :: sbtgrid_lnr0
      real*8 :: sbtgrid_lnk0
      real*8 :: sbtgrid_lnrange
      integer :: sbtgrid_N ! VB: reduced from 16384
      logical :: use_logsbt_for_radial_hse_integration

      ! Only applicable for (use_lvl)
      logical :: use_logsbt_lvltriples

      logical :: use_2d_corr

      logical :: fixed_spin_moment
      real*8  :: spin_moment
      real*8,dimension(2)  :: fixed_spin_moment_electrons

      logical :: walltime_requested
      integer :: walltime_limit

      ! This is the work space size allowed for the "distributed" hartree potential
      ! (which stores all sorts of things across the grid)
      ! Unit: Megabytes.
      real*8 :: hartree_worksize

      ! If this flag is true, an external embedding potential will act on the
      ! electronic self-consistency loop (this is the default). If explicitly set to
      ! false, only an embedding energy will
      ! be calculated, but the Kohn-Sham potential will not be modified at all.
      logical :: full_embedding

      ! Runtime choices associated with the Kerker preconditioner
      real*8  :: precondition_kerker_q0 ! not too bad default, but could be system-dependent
      integer :: precondition_max_l     ! this seems to be sufficient for all systems where preconditioner matters
      logical :: use_kerker_preconditioner
      logical :: precondition_before_mixer
      logical :: use_preconditioner_turnoff_charge
      real*8  :: preconditioner_turnoff_charge
      real*8  :: preconditioner_turnoff_sum_ev
      real*8  :: preconditioner_turnoff_energy
      logical :: use_kerker_factor
      real*8  :: kerker_coeff

      ! Runtime choices associated with the dielectric preconditioner
      logical :: use_dielectric_preconditioner

      logical :: restart_from_density
      logical :: freeze_restart_density
      logical :: keep_restart_info
      logical :: restart_file_exists
      character*40 :: restart_read_file
      character*40 :: restart_write_file
      character*100 :: restart_3c_file
      character*100 :: restart_exchange_matr_file
      character*100 :: restart_coulomb_matr_for_exx
      logical :: restart_write
      logical :: restart_read
      integer :: restart_save_iterations
      logical :: use_restart_save_iterations
      logical :: restart_exx_write
      logical :: restart_exx_read
      logical :: restart_exx_3c_int
      ! Igor, to restart RPA
      logical :: restart_rpa_file_exists
      character*40 :: restart_rpa_read_file
      character*40 :: restart_rpa_write_file
      logical :: restart_rpa_write
      logical :: restart_rpa_read
      integer :: rpa_freq_start
      real*8  :: current_crpa
      ! Igor, to restart PT2
      logical :: restart_pt2_file_exists
      character*40 :: restart_pt2_read_file
      character*40 :: restart_pt2_write_file
      logical :: restart_pt2_write
      logical :: restart_pt2_read
      integer :: pt2_finish
      real*8  :: current_pt2
      real*8  :: current_pt2_os
      real*8  :: current_pt2_ss
      ! Igor, to skip the scf procedure for postscf calculations
      logical :: skip_scf_for_postscf
      logical :: restart_for_postscf
      integer(kind=4), allocatable, dimension(:) :: fv_orbs
      logical :: fv_filter
      integer :: fv_orbs_n
      character*40 :: fv_orbs_file

      logical :: restart_relaxations
      logical :: hessian_to_restart_geometry
      logical :: write_restart_geometry
      logical :: use_distributed_hessian

      logical :: final_forces_cleaned
      logical :: output_in_original_unit_cell

      logical :: grid_storage_initialized
      logical :: got_n_compute_maxes
      logical :: sc_check(2)

      ! VB - default switched to true.
      logical :: prune_basis_once

      logical :: use_local_index
      logical :: use_load_balancing
      logical :: use_alltoall

      logical :: restrict_kpoint_to_smp_node
      integer :: max_tasks_per_smp_node

      logical :: first_integration

      logical :: parallel_grid_partitioning

      logical :: force_mpi_virtual_topo

      logical :: use_metis_batch_distribution

      logical :: frozen_core_scf
      real*8 :: frozen_core_scf_cutoff
      real*8 :: frozen_core_scf_factor

      logical :: recompute_batches_in_relaxation

      logical:: out_vacuum_potential
      integer:: out_vacuum_potential_x_grid, out_vacuum_potential_y_grid, out_vacuum_potential_z_grid
      real*8::  out_vacuum_potential_z1
      real*8::  out_vacuum_potential_z2
      logical:: out_embedding_potential
      integer:: out_embedding_potential_z_grid
      real*8::  out_embedding_potential_z1
      real*8::  out_embedding_potential_z2
      real*8::  vacuum_z_level
      logical :: flag_set_vacuum_level
      logical::  use_dipole_correction
      logical::  calculate_work_function
      character(LEN=10) :: dipole_correction_method

      integer ::   dos_n_en_points
      real*8  ::   dos_low_energy
      real*8  ::   dos_high_energy
      real*8  ::   dos_alpha
      integer ::   atom_dos_n_en_points
      real*8  ::   atom_dos_low_energy
      real*8  ::   atom_dos_high_energy
      real*8  ::   atom_dos_alpha
      integer ::   l_proj_dos_n_en_points
      real*8  ::   l_proj_dos_low_energy
      real*8  ::   l_proj_dos_high_energy
      real*8  ::   l_proj_dos_alpha

      ! This flag specifies the communication type used for MPI:
      !   BlueGene, or another platform. Currently, this only affects
      !   the distributed_spline_storage flag in the Hartree potential update
      !   but may be extended in the future
      integer :: communication_type
      integer, parameter :: standard_mpi = 1
      integer, parameter :: BlueGene_mpi = 2
      integer, parameter :: mpi_1coeff   = 3
      integer, parameter :: calc_hartree = 4
      integer, parameter :: shmem_comm   = 5
      integer, parameter :: Cray_mpi = 6

      ! temporary flag to force neutralized energy functional also in cluster case
      ! this will be the default at some point.
      logical :: force_new_functional

      ! flags and choices for Born-Oppenheimer MD for use in the module molecular_dynamics
      integer      :: MD_maxsteps
      real*8       :: MD_time
      real*8       :: MD_tstep  ! this is equal to 1 fs ...
      real*8       :: MD_temperature
      real*8       :: MD_init_temperature
      real*8       :: MD_tau_berendsen
      real*8       :: MD_nu_andersen
      real*8       :: MD_tau_BDP
      integer      :: MD_gle_ns
      real*8       :: MD_Q_NP
      integer      :: MD_g_DOF
      character*30 :: MD_ensemble
      real*8       :: NVE_damping_factor
      logical      :: MB_velocity_initialization
      logical      :: MB_clean_rotations
      integer      :: seed
      logical      :: MD_RNG_seeded
      character*40 :: MD_restart_file
      character*20 :: MD_thermostat_units
      logical      :: MD_initialconf_from_restart
      logical      :: MD_time_restart
      logical      :: MD_restart_binary
      logical      :: check_MD_stop
      logical      :: plumed_plugin
      character*40 :: plumed_file
      logical      :: plumed_new_plugin
      character*30, dimension(:)  , allocatable :: MD_schedule_ensemble
      real*8,       dimension(:)  , allocatable :: MD_schedule_temperature
      real*8,       dimension(:)  , allocatable :: MD_schedule_time
      real*8,       dimension(:)  , allocatable :: MD_schedule_tau_berendsen
      real*8,       dimension(:)  , allocatable :: MD_schedule_Q
      real*8,       dimension(:)  , allocatable :: MD_schedule_damping_factor
      real*8,       dimension(:)  , allocatable :: MD_schedule_nu_andersen
      real*8,       dimension(:)  , allocatable :: MD_schedule_tau_BDP
      logical,      dimension(:)  , allocatable :: MD_schedule_random_BDP
      !XZL: added for PIMD
      integer      :: PIMD_maxsteps
      real*8       :: PIMD_time
      real*8       :: PIMD_tstep  ! this is equal to 0.5 fs ...
      real*8       :: PIMD_temperature
      real*8       :: PIMD_init_temperature
      real*8       :: PIMD_Q_NP
      integer      :: PIMD_g_DOF
      character*30 :: PIMD_ensemble
      logical      :: PIMD_RNG_seeded
      character*40 :: PIMD_restart_file
      character*20 :: PIMD_thermostat_units
      logical      :: PIMD_initialconf_from_restart
      logical      :: PIMD_time_restart
      logical      :: PIMD_restart_binary
      logical      :: check_PIMD_stop
      !CC: Additional flags for thermodynamic integration
      real*8,       dimension(:)  , allocatable :: TDI_lambda_start
      real*8,       dimension(:)  , allocatable :: TDI_lambda_end
      real*8,       dimension(:)  , allocatable :: TDI_QHA_E0
      character*100,dimension(:)  , allocatable :: TDI_QHA_file
      logical     , dimension(:)  , allocatable :: TDI_QHA_file_exists
      !CC: Additional variables for MD_QH_init
      integer, dimension(:)       , allocatable :: MD_QH_first_atom
      integer, dimension(:)       , allocatable :: MD_QH_last_atom
      real*8,  dimension(:)       , allocatable :: MD_QH_temperature
      character*100,dimension(:)  , allocatable :: MD_QH_file
      !CC: Flag for skipping the SCF cycles
      logical                                   :: SKIP_SCF
      !CC: Energy density related flags
      logical                                   :: flag_energy_density
      logical                                   :: flag_chetty_martin_energy_density
      logical                                   :: flag_harris_foulkes_energy_density
      logical                                   :: flag_project_energy_density

      integer, parameter :: PP_ANYWAY_NOTHING = -1
      integer, parameter :: PP_ANYWAY_EVERYTHING = 99
      integer :: postprocess_anyway = PP_ANYWAY_NOTHING

      ! for wf-extrapolation
      integer, parameter :: WF_EXTRA_NONE = 1
      integer, parameter :: WF_EXTRA_FUNC = 2
      integer, parameter :: WF_EXTRA_NIKL = 3
      integer :: wf_extra_type
      logical :: use_wf_extrapolation
      logical :: wf_extra_use_densmat

      logical :: out_zero_multipoles
      logical :: out_hartree_multipoles
      ! how much output ... ?
      ! 'normal' = reduced in relaxations and MD
      ! 'full'   = print everything
      ! 'MD_light' = few details per time step in MD
      character*20 :: output_level

      ! basis-dependent cutoff potential default: this does not seem to do anything, test to your own specification.
      real*8 :: basis_dep_cutoff_default

      ! some band plotting data, where and how many bands ...
      logical :: plot_band_allocated
      integer :: band_initialized
      real*8, dimension(:,:),allocatable :: band_begin, band_end
      integer,dimension(:),allocatable   :: n_points_in_band

      ! TZ: some dielectric output data, braodening_type, broadening_parameter,
      ! direction i and direction j
      logical :: plot_dielectric_allocated
      integer :: dielectric_initialized
      integer, dimension(:),allocatable :: broaden_type
      real*8, dimension(:),allocatable  :: broaden_para
      integer, dimension(:),allocatable :: direc_1
      integer, dimension(:),allocatable :: direc_2
      logical :: calc_px_dielectric
      logical :: calc_py_dielectric
      logical :: calc_pz_dielectric
      integer :: die_broaden_type
      real*8  :: die_broaden_para
      character*30 :: die_method
      ! Outputting SOC eigenvectors
      integer, dimension(:), allocatable :: write_soc_eigenvectors

      ! same for band directions requested for output on s.c.f. k-point grid
      real*8, dimension(:,:),allocatable :: band_scf_begin, band_scf_end

      ! DFPT flags
      real*8 :: DFPT_sc_accuracy_dm
      real*8 :: DFPT_mixing
      real*8 :: DFPT_width
      logical :: restart_write_DM1
      logical :: restart_read_DM1
      ! some DFPT_phonon_band plotting data, where and how many bands ...
      real*8, dimension(:,:),allocatable :: DFPT_phonon_band_begin, DFPT_phonon_band_end
      integer,dimension(:),allocatable   :: n_points_in_band_DFPT_phonon

      ! Variables that control the type of a transport calculation
      ! (Paula's version so far)
      logical:: transport_lead_calculation
      logical:: transport_calculation
      logical:: flag_output_effective_potential

      ! flag for atom_bsse
      logical :: calculate_atom_bsse

      ! flags and parameters for time stepping
      real*8  :: TDDFT_time
      real*8  :: TDDFT_dt
      character*30 :: TDDFT_propagator
      integer :: TDDFT_output_every
      integer :: TDDFT_write_info
      integer :: TDDFT_stochastic
      integer :: TDDFT_iseed
      logical :: flag_run_tddft_real_time_propagation

      ! flags and variables for the casida subroutines
      logical :: flag_neutral_excitation
      logical :: neutral_excitation_rpa
      logical :: neutral_excitation_tdhf
      logical :: neutral_excitation_tddft
      logical :: excited_mode_singlet
      logical :: excited_mode_triplet
      integer :: num_excited_states
      integer :: libxc_tddft_x
      integer :: libxc_tddft_c
      integer :: libxc_tddft_xc
      logical :: libxc_tddft_chosen_x
      logical :: libxc_tddft_chosen_c
      logical :: casida_setup_scalapack

      ! flags for bse
      logical :: neutral_excitation_bse
      logical :: read_qpe
      logical :: write_qpe
      logical :: write_ovlp_3ks
      logical :: write_ovlp_3fn
      logical :: write_screened_coulomb
      logical :: bse_singlet
      logical :: bse_triplet
      integer :: bse_lower_limit
      integer :: bse_upper_limit
      logical :: bse_reduce_matrix
      real*8  :: bse_reduce_unocc
      real*8  :: bse_reduce_occ
      logical :: output_ks_coulomb_integral

      ! flags for band mulliken
      integer :: band_mulliken_orbit_num

      ! flags for k point of dipole_trans_mat
      logical :: dipole_trans_mat_out
      integer :: dipole_trans_mat_k ! default gamma point
      integer :: dipole_trans_mat_orbit_num

      !flags and variables for the post-SCF frozen-core
      ! and frozen-virtual calculations
      integer :: n_frozen_shell
      logical :: flag_frozen_core_postSCF
      integer :: n_frozen_virtual_states
      logical :: flag_frozen_virtual_postSCF
      !flags and variables for the energy shifting for post-SCF calculations
      integer :: en_shift_type
      logical :: flag_en_shift
      real*8  :: en_shift_constant
      real*8  :: coupling_pt2_factor
      real*8  :: coupling_pt2_screen
      real*8  :: coupling_pt2_shift
      ! ==============osRPA==========================
      ! Ec[scsRPA] = r1*Ec[osRPA] + r2*Ec[ssRPA] NOTE: might change soon as Ec[ssRPA] is not correct at present
      ! r1 is the rescaled factor for osPRA
      ! r2 is the factor for osRPA
      real*8, dimension(2)  :: scsrpa_coeff
      ! special_radius (sr) are the special radius of non-interacting
      !   linear response matrix in two spin channels for rpa
      !   with omega=0, which are taken as the descriptors
      !   for the strong correlation in a given system
      real*8, dimension(2)  :: special_radius
      ! 1-r1*sr/(r2-sr)*exp(-r3*(eg+r4)**2)
      ! r1 and r2 are the parameters to determine the
      !   the level-shift value and shape w.r.t the strong
      !   correlation (the parameter fo sr)
      ! r3 and r4 are the parameters to determine the
      !   the level-shift value and shape w.r.t. to the
      !   energy gap (eg)
      real*8, dimension(10)  :: renormalized_eg
      ! switch to a renormalization solution of osRPA if special radius > c_osrpa_threshold
      real*8  :: c_osrpa_threshold
      ! for the coupling-constant integration for osRPA
      integer :: n_lambda_osrpa
      ! for osrpa correlation at different orders
      integer:: c_osrpa_order
      real*8,allocatable, dimension(:):: c_osrpa
      ! ==============osRPA==========================
      !flags for the CI methods
      integer :: ci_truncation
      logical :: flag_single_excitation
      logical :: flag_single_spin_channel
      logical :: flag_turn_off_hartree
      real*8  :: ci_accuracy_etot
      integer :: n_scf_ci
      real*8  :: ci_coeff_threshold
      Integer :: n_ci_acc_save
      Integer :: ci_acc_method
      Double precision :: ci_acc_ovlp_thresh
      ! flags for coupled cluster methods
      Integer :: CC_calc_method
      Integer :: CC_max_cycle
      Integer :: CC_acc_method
      Integer :: CC_n_RLE_sv
      Integer :: CC_restart_point
      Integer :: CC_restart_step
      INteger :: CC_mkl_nth
      Double precision :: CC_converg_thresh
      Double precision :: CC_RLE_bas_thresh
      Double precision :: CC_lv_shift
      Integer :: CC_sv_method
      Integer :: CC_n_domain
      Integer :: CC_MPI_max_len
      Integer :: CC_work_tmp_size
      Integer :: CC_DIIS_n_sv
      Integer :: CC_DIIS_step
      Integer :: CC_DIIS_err_method
      Integer :: CC_PT_n_out
      Logical :: flag_cc_general
      Logical :: CC_collective
      Logical :: CC_use_disk
      Logical :: CC_restart_flag
      Logical :: CC_check_flag
      Logical :: CC_gamma_flag
      Logical :: CC_abcd_sv_flag
      Logical :: CC_mem_check
      Logical :: CC_super_system_flag
      Logical :: CC_dynamic_threading
      Character (len = 20) :: CC_solver
      Character (len = 10) :: CC_sv_control
      Character (len = 10) :: CC_method
      ! flags for sampling RPA potentials along AC path
      integer :: rpa_along_ac_path_grid
      real*8, dimension(:), allocatable :: rpa_potentials_along_ac_path
      ! DFT components for printing out
      integer :: n_dft_methods
      integer, dimension(10)             :: dft_print_out

      logical :: add_embedding_grids

      ! parameters for implicit solvent effect
      ! solvation embedding method
      integer, parameter :: SOLVENT_INVALID       = -1
      integer, parameter :: SOLVENT_UNDEF         = 0
      integer, parameter :: SOLVENT_MPE           = 1
      integer, parameter :: SOLVENT_MPB           = 2

      ! implicit solvation
      integer  :: solvent_method
      ! mpb
      real(dp) :: permittivity_switching_width
      ! mpe_reaction_field
      integer  :: mpe_nonel_model
      integer  :: mpe_nonel_model_bulk
      real(dp) :: mpe_nonel_alpha
      real(dp) :: mpe_nonel_beta
      real(dp) :: mpe_nonel_alpha_if
      real(dp) :: mpe_nonel_alpha_bulk
      real(dp) :: mpe_nonel_beta_bulk
      integer  :: mpe_factorization_type
      integer  :: mpe_n_boundary_conditions
      real(dp) :: mpe_solvent_permittivity
      real(dp) :: mpe_bulk_permittivity

      integer  :: mpe_xml_loglevel
      character(len=40) :: mpe_xml_logfile
      ! defaults of the following variables will be set in the module
      !  "mpe_reaction_field"
      ! everything that is still negative by then has not been set by the user
      integer  :: mpe_skip_first_n_scf_steps
      real(dp) :: mpe_degree_of_determination
      integer  :: mpe_lmax_rf
      integer  :: mpe_lmax_ep
      logical  :: mpe_lmax_ep_from_basis
      integer  :: mpe_lmax_ep_additional_order
      integer  :: mpe_n_centers_ep
      logical  :: mpe_lmax_ep_increase_cluster
      real(dp) :: mpe_f_sparsity_threshold

      ! tolerance for bad solutions of mpe equations
      real(dp) :: mpe_tol_adjR2
      logical  :: mpe_tol_adjR2_wait_scf

      ! isc_implicit_solvent_cavity
      integer  :: isc_cavity_type
      logical  :: isc_cavity_restart_read
      character(len=80) :: isc_cavity_restart_read_file
      logical  :: isc_cavity_restart_write
      character(len=80) :: isc_cavity_restart_write_file
      character(len=80) :: isc_record_cavity_creation_file
      integer  :: isc_record_every_n_steps
      logical  :: ifp
      real(dp) :: ifp_dist
      real(dp) :: ifp_normal(3)
      real(dp) :: isc_isodensity_value
      logical  :: isc_calculate_surface_and_volume
      real(dp) :: isc_kill_ratio
      logical  :: isc_surface_curvature_correction
      logical  :: isc_try_restore_convergence
      ! defaults of the following variables will be set in the module
      !  "isc_implicit_solvent_cavity"
      ! everything that remains negative by then has not been set by the user
      ! extend radius of initial Lebedev sphere
      real(dp) :: isc_initial_radius
      integer  :: isc_surface_sampling
      real(dp) :: isc_dt
      real(dp) :: isc_rho_k
      real(dp) :: isc_rho_k_mod_constant
      real(dp) :: isc_g_k
      real(dp) :: isc_rep_k
      real(dp) :: isc_rho_rel_deviation_threshold
      real(dp) :: isc_gradient_threshold
      real(dp) :: isc_dynamics_friction
      integer  :: isc_update_nlist_interval
      integer  :: isc_max_dyn_steps
      real(dp) :: ifp_rmin
      real(dp) :: ifp_rmax
      logical  :: ifp_manual
      real(dp) :: ifp_rlog
      real(dp) :: ifp_area_per_point
      integer  :: ifp_n_shells_log
      integer  :: ifp_n_angular_max
      real(dp) :: ifp_min_dr
      integer  :: ifp_n_shells_lin

      ! Two extra variables that control the accuracy verification of onsite integrals
      ! onsite_accuracy_threshold is the maximum value by which two onsite Hamiltonian integrals
      !   may differ before the calculation is stopped
      ! override_integration_accuracy allows one to force the calculation to continue.
      real*8 :: onsite_accuracy_threshold      ! 30 meV. Could be too tight. Needs testing.
      logical :: override_integration_accuracy

      ! CUDA and CUBLAS flags
      logical :: use_gpu
      logical :: use_gpu_density
      logical :: use_gpu_hamiltonian
      logical :: use_gpu_forces

      ! Identification strings set by user to identify the current calculation
      ! as belonging to a set and/or subset of calculations
      character(len=:), allocatable :: calculation_set
      character(len=:), allocatable :: calculation_subset

      ! CPU consistency checks
      logical :: flag_cpu_consistency
      real*8 :: cpu_consistency_threshold

      ! Whether myid.eq.0 should output log files in json format
      logical :: out_aims_json_log

      ! ScaLAPACK block size
      ! Negative value: not defined in control.in file
      integer :: scalapack_block_size

      logical :: output_matrix_io_timing

      ! use ELSI (ELectronic Structure Infrastructure)
      logical :: first_elsi_call
      logical :: use_elsi
      logical :: use_elsi_dm
      logical :: use_elsi_ev

      ! 0: auto
      ! 1: ELPA
      ! 2: libOMM
      ! 3: PEXSI
      ! 4: Not available
      ! 5: SLEPc-SIPs
      ! 6: NTPoly
      integer :: elsi_solver

      ! ELSI output level
      integer :: elsi_out_level
      integer :: elsi_out_json
      ! ELSI restart
      logical :: elsi_read_dm
      logical :: elsi_read_dm_done
      logical :: elsi_write_dm
      logical :: elsi_write_dm_this
      integer :: elsi_write_dm_freq
      logical :: do_elsi_rw
      logical :: elsi_extrap_dm
      real*8 :: elsi_rw_filter
      real*8 :: elsi_filter
      ! ELPA solver
      integer :: elsi_elpa_solver
      ! Number of single precision ELPA steps
      integer :: elsi_elpa_n_single
      ! Use GPU in ELPA
      integer :: elsi_elpa_gpu
      ! Number of ELPA steps before OMM
      integer :: elsi_omm_n_elpa
      ! OMM flavor
      integer :: elsi_omm_flavor
      ! OMM minimization tolerance
      real*8 :: elsi_omm_tol
      ! Number of processors for each PEXSI pole
      integer :: elsi_pexsi_np_per_pole
      ! Number of processors for symbolic factorization
      integer :: elsi_pexsi_np_symbo
      ! EigenExa method
      integer :: elsi_eigenexa_method
      ! Number of slices in SIPs
      integer :: elsi_sips_n_slice
      ! Slicing type in SIPs
      integer :: elsi_sips_slice_type
      ! Buffer to expand the slices
      real*8 :: elsi_sips_buffer
      ! Number of ELPA steps before SIPs
      integer :: elsi_sips_n_elpa
      ! Density matrix purification method
      integer :: elsi_ntpoly_method
      ! Density matrix purification accuracy
      real*8 :: elsi_ntpoly_tol
      ! Density matrix purification accuracy
      real*8 :: elsi_ntpoly_filter
      ! MAGMA solver
      integer :: elsi_magma_solver

      ! related to PIMD wrapper, initialization of electronic density
      logical :: need_dens_superpos
      logical :: recalc_dens

      ! related to constraining potential
      character*40 :: potconstype
      real*8 :: potconspar

      ! Magnetic response flags
      logical :: magnetic_response
      ! This is the line starting with 'magnetic_respone' in
      ! control.in. It is processed in check_consistency_of_keywords.
      character(132), public :: magnetic_response_pars
      logical :: mr_experimental
      integer :: dfpt_iter_limit
      real*8 :: dfpt_accuracy_n1
      real*8 :: dfpt_linear_mix_param
      integer :: dfpt_pulay_steps
      logical :: use_dfpt_pulay
      logical :: output_sxml
      real*8 :: mr_gauge_origin(3) ! Default is the center of mass!
      character(100), public :: sxml_filename

      ! flags for fo_dft
      character*10 :: fo_dft_option, fo_info_file
      character*6  :: fo_comb_m
      character(:), allocatable :: fo_folder1, fo_folder2
      character*4 :: fo_type
      character*5 :: fo_fragment_number
      logical :: fo_dft, fo_fragment, fo_density
      logical :: fo_finalise, fo_orbitals, use_fo_potential
      integer :: fo_orb1, fo_orb2, fo_range1, fo_range2, fo_verbosity
      logical :: occ_fo, fo_potential_file, fo_file_exists, fo_deltaplus

      logical :: force_single_restartfile
      logical :: restart_eigenvectors_periodic

      ! cutoff for tau in meta-GGAs. Default is zero i.e. not applied.
      ! Different meta-GGAs have different accuracy when introducing this
      ! parameter, hence it is not applied as standard, but can be set manually
      real*8 :: tau_threshold

      ! Flags for symmetry
      logical :: use_rlsy_symmetry
      logical :: skip_rlsy_symmetry

      contains

      subroutine set_runtime_choices_defaults ( )
         implicit none

         dry_run = .false.

         override_relativity = .false.
         flag_scaled_zora = .false.
         flag_KH_post_correction = .false.
         flag_KH_core_states = .false.
         use_small_component = .false.

         flag_eg       = .false.
         flag_kernel   = .false.
         flag_cell     = .false.

         override_warning_nobasis = .false.
         override_warning_negative_nucleus = .false.
         override_warning_libxc = .false.

         hse_unit = 'A'
         hse_omega = 0.11d0
         hse_omega_hf = 0.0d0
         hse_omega_pbe = 0.0d0
         lc_dielectric_constant = 1.0d0

         cutCb_rcut = 15.117809d0
         cutCb_width_factor = 8.0d0
         cutCb_rcut_factor = 1.d0
         flag_cutCb_width = CUTCB_FROM_DEFAULT
         flag_cutCb_rcut  = CUTCB_FROM_DEFAULT

         ers_width = 1.0d0
         ers_omega = 1.0d0
         ers_rcut = 1.0d0*bohr
         erfc_omega = 1.0d0
         flag_dftpt2_xc = 0
         flag_post_xc = 0
         dftpt2_Ex_hf       = 0.0d-1
         dftpt2_Ec_osPT2    = 0.0d-1
         dftpt2_Ec_ssPT2    = 0.0d-1
         dftpt2_Ec_oslrcPT2 = 0.0d-1
         dftpt2_Ec_sslrcPT2 = 0.0d-1
         lrc_pt2_unit = 'B'
         lrc_pt2_omega = 0.20d0
         lrc_pt2_started = .false.

         partition_type = 8
         stratmann_a = 0.64d0
         n_empty = -1
         override_illconditioning = .false.
         prodbas_nb = 0
         flag_basis = FLAG_BASIS_ONLY_PRODS
         use_smooth_prodbas_threshold = .false.
         use_asym_inv_sqrt = .false.
         condition_penalty = 0.d0
         use_aufbau_principle = .false.
         use_density_matrix = .false.
         density_update_method = UPDATE_AUTO
         use_density_matrix_hf = .false.
         autoselect_density_method = .false.
         split_updates = .false.
         lda_update_switch = .false.
         measure_forces = .false.
         use_ovlp_swap = .false.
         flag_compute_kinetic = .false.

         analytic_potential_average = .false.

         entropy_thr = 1e-15

         atomic_solver = 0
         use_atomic_solver_fns_as_minimal_basis = .false.

         multipole_feed_back_parameter = 1.0d0
         flag_delta_rho_in_multipole =.false.

         adjust_scf = .true.
         adjust_scf_always = .false.
         adjust_scf_iteration = 2
         charge_mix_param_set = .false.
         occupation_width_set = .false.

         coeff_thr = 1e-10

         sc_init_iter =   1001
         sc_init_factor = 1.d0
         scf_restarted = .false.
         restarting_scf = .false.

         check_sc_abandon_etot = .true.
         sc_abandon_etot_iter = 5
         sc_abandon_etot_threshold = 1000

         vdw_convergence_threshold = 1.d-8/hartree

         use_scsmp2 = .false.
         ps_mp2 = 1.0
         pt_mp2 = 1.0

         force_lebedev = 0
         force_hartree_log_grid = .true.

         ! SK: Changed Aug 4, 19. New default. Switch to .false. only if DFPT is requested
         compensate_multipole_errors = .true.
         normalize_initial_density = .true.

         packed_matrix_format = 0

         step_x = 0
         step_y = 0
         step_z = 0

         vacuum_length = 0

         Emin = 0
         Emax = 0
         k_point_mom = 1

         widthone = 0
         widthtwo = 0
         omega_min = 0
         omega_max = 0
         n_omega = 0
         n_omega_reset = 0
         ep1 = 'x'
         ep2 = 'x'
         use_gauss = .False.
         read_q_points = .False.
         k_k_method = 2
         greenwood_method = 1
         percent_or_ev = 2
         n_greenenergy = 1
         fermistatesabove=0
         fermistatesbelow=0
         dist_from_vbm=0.0
         fermispacing=0.0

         out_fermisurface = .false.

         out_mommat = .false.

         esp_Emin = 0
         esp_Emax = 0
         esp_k_point = 1
         esp_min_radius = 3.0
         esp_max_radius = 8.0
         esp_n_max_radius = 5
         esp_grid_type = 1
         flag_esp_restraint = .false.
         esp_restr_reference = 1
         esp_restraint_strength = 1.0d0

         mbd_beta = 0d0

         out_onsite = .false.

         flag_found_unfold = .false.
         flag_hdf5_unfold  = .false.
         flag_ASCII_unfold = .true.

         wantDebug = .false.

         out_xml = .false.

         out_mat_par = 'n'
         out_mat_par_format = 'asc'
         coulomb_integral_format = 'full'
         out_nuc_matrix = .false.
         out_t_plus_v_matrix = .false.
         out_nuc_nuc_repulsion = .false.
         exx_band_structure_version = 0


         out_l_proj_dos  = .false.
         out_l_proj_dos_tetrahedron  = .false.
         out_dos         = .false.
         out_dos_tetrahedron  = .false.
         out_atom_dos    = .false.
         out_atom_dos_tetrahedron    = .false.
         out_species_dos = .false.
         out_species_dos_tetrahedron = .false.
         out_hirshfeld_always = .false.
         out_hessian = .false.
         out_esp = .false.
         out_esp_full = .false.

         flag_run_mulliken = .false.
         flag_run_lowdin = .false.

         flag_gOpenMol_format = .false.

         out_k_points_eigenvalues = 1
         out_soc_eigenvalues = .false.

         out_k_point_list = .false.

         out_sorted_eigenvalue_list = .false.
         overwrite_existing_cube_files = .false.
         cube_default_size_safeguard = 5e7
         out_cube_nth_iteration = .false.
         use_old_cube_routine =.false.
         cube_content_unit = 'legacy'

         out_dgrid = .FALSE.

         potential_well_requested = .false.
         potential_well_start = 0
         potential_well_end = 0
         potential_well_depth = 0

         n_min_points_in_division = 1
         grid_partitioning_method = 1
         grouping_factor = 2
         use_tetgen = .false.
         use_qhull = .false.
         use_nearest_neighbours = .false.
         n_nearest_neighbours = 5
         min_batch_size = 1
         batch_size_hard_limit = 200
         use_hashed_maxmin = .false.
         use_integer_hash = .false.
         out_batch_statistics = .false.

         flag_force_complex_eigenvectors = .false.

         use_scalapack = .false.
         use_elpa = .false.
         use_scalapack_mbd = .false.
         use_elpa_mbd = .false.

         use_scalapack_DFPT_phonon = .false.
         solver_method = 0
         ELPA_AEO = 0

         collect_eigenvectors = .true.
         orthonormalize_evs = .true.

         use_lapack_fast = .false.

         use_cg = .false.
         cg_preconditioner = prec_inv_ovlp
         initial_ev_solutions = 5
         max_cg_iterations = 100
         lopcg_block_size = 1
         lopcg_tol = 1.0e-6
         lopcg_start_tol = 0.0d0
         lopcg_adaptive_tolerance = .false.
         lopcg_auto_blocksize = .false.
         lopcg_omega = 1.2d0
         lopcg_use_all_evs = .false.
         lopcg_slide_ev_block = .false.
         lopcg_skip_rate = 0

         n_k_points_xyz  = 1
         n_k_points_xyz_nosym  = 1
         n_k_points_dos_xyz  = 1
         k_enhancement_factor  = 2
         k_points_offset = 0.0
         symmetry_thresh                 = 1d-4
         use_symmetry_reduced_k_grid     = .true.
         use_symmetry_reduced_spg        = .false.
         use_spg_mv_mm                   = .false.
         use_spg_full_Delta              = .true.
         reconstruct_proper_only         = .true.
         get_full_density_first          = .true.
         use_k_phase_exx                 = .true.
         read_k_points_from_file         = .false.
         pert_dos_on                     = .false.
         postscf_eigenvalues             = .false.
         real_eigenvectors               = .true.
         numerical_stress_save_scf       = .true.
         use_analytical_stress           = .false.
         use_inv_symmetry                = .true.
         use_full_symmetry               = .false.
         flag_k_fine_grid                = .true.
         AS_components       = 6
         AS_symmetric_output = .false.
         stress_for_relaxation = 0
         flag_external_pressure = .false.

         far_distance_hartree_multipole_radius_threshold         = 1e-12
         far_distance_hartree_multipole_moment_radius_threshold  = 1e-10
         far_distance_hartree_multipole_moment_threshold         = 1e-10
         far_distance_hartree_radius_threshold                   = 1e-10
         far_distance_adaptive_hartree_radius_threshold          = 1e-8
         multipole_radius_free_threshold                         = 0.d0
         Ewald_radius                                            = 3.0d0
         Ewald_radius_automatic                                = .true.
         fast_Ewald = .true.

         hartree_fp_function_splines                           = .true.

         distribute_leftover_charge     = .false.
         l_hartree_far_distance = 10
         use_hartree_non_periodic_ewald = .false.
         hartree_d_matrix_method = 1 ! 1:submatrix 2:pseudo inverse
         
         gridwidth_in_pm = 60

         n_k_points_group = 1
         hartree_potential_fourier_part_threshold = 5e-7
         extra_adding_to_hartree_potential_distance = 2.0d0
         packed_matrix_threshold = 1e-13
         Adams_Moulton_integrator = .true.
         legacy_monopole_extrapolation = .false.
         RELAX_LTRM_dev = .false.
         relax_accuracy_stress_factor = 10.0d0
         use_stress_to_check_relaxation = .false.
         remove_unitary_force_components = 0
         store_eigenvectors_to_disk_in_relaxation =.false.

         energy_tolerance           = 0.001/hartree
         aggregated_energy_tolerance= 0.003/hartree
         aggregated_energy_tolerance_per_atom = 0.0001/hartree
         harmonic_length_scale      = 0.025/bohr
         max_atomic_move            = 0.2/bohr
         min_line_step              = 0.1
         line_step_reduce           = 0.25
         line_search_tol            = 3d0
         line_step_reduce_automatic = .true.
         max_relaxation_steps       = 1000
         bfgs_mixing_factor         = 0.2d0
         bfgs_mixing_cap            = 2d0
         bfgs_extrapolation_cap     = 4d0
         bfgs_extrapolation_on      = .false.

         init_hess_type             = HESS_DIAG
         init_hess_diag              = 25.d0 / (hartree/bohr**2)
         init_hess_lv_diag              = 25.d0 / (hartree/bohr**2)
         init_hess_lindh_diag        = 2.0d0 / (hartree/bohr**2)
         init_hess_lindh_thres       = 15.d0

         flag_auxil_basis = PRODBAS_FULL
         hybrid_coeff = 0.0d-1
         hf_multiplicity = 1
         hf_version = HF_DM
         RI_type = RI_V
         sparse_o3fn = .false.

         sbtgrid_lnr0 = -38.d0
         sbtgrid_lnk0 = -25.d0
         sbtgrid_lnrange = 45.d0
         sbtgrid_N = 4096
         use_logsbt_for_radial_hse_integration = .false.

         use_logsbt_lvltriples = .true.

         fixed_spin_moment = .false.

         walltime_requested = .false.
         walltime_limit = 157784630

         precondition_kerker_q0 = 2.0d0
         precondition_max_l = 0
         use_kerker_preconditioner = .false.
         precondition_before_mixer = .false.
         use_preconditioner_turnoff_charge = .false.
         preconditioner_turnoff_charge = 1d-3
         preconditioner_turnoff_sum_ev = 1d-1
         preconditioner_turnoff_energy = 1d-3
         use_kerker_factor = .false.
         kerker_coeff = 1.0d0

         use_dielectric_preconditioner = .false.

         restart_from_density = .false.
         freeze_restart_density = .false.
         keep_restart_info = .false.
         restart_file_exists = .false.
         restart_read_file = 'aims_restart_information.dat'
         restart_write_file='aims_restart_information.dat'
         restart_3c_file='aims_3c_int_'
         restart_exchange_matr_file='aims_exx_matr_'
         restart_coulomb_matr_for_exx='aims_cm_'
         restart_write = .false.
         restart_read  = .false.
         restart_save_iterations = 0
         use_restart_save_iterations = .false.
         restart_exx_write = .false.
         restart_exx_read = .false.
         restart_exx_3c_int = .false.
         restart_rpa_file_exists = .false.
         restart_rpa_read_file = 'aims_restart_rpa_information.dat'
         restart_rpa_write_file='aims_restart_rpa_information.dat'
         restart_rpa_write = .false.
         restart_rpa_read  = .false.
         rpa_freq_start = 0
         current_crpa = 0.0d0
         restart_pt2_file_exists = .false.
         restart_pt2_read_file = 'aims_restart_pt2_information.dat'
         restart_pt2_write_file='aims_restart_pt2_information.dat'
         restart_pt2_write = .false.
         restart_pt2_read  = .false.
         pt2_finish = 0
         current_pt2 = 0.0d0
         current_pt2_os = 0.0d0
         current_pt2_ss = 0.0d0
         skip_scf_for_postscf = .false.
         restart_for_postscf  = .false.
         fv_filter = .false.
         fv_orbs_file='aims_fv_orbs.dat'
         ! TODO: deallocate arrays in this module at end of calculation, rather
         !       than beginning of next calculation
         if (allocated(fv_orbs)) deallocate(fv_orbs)
         if (allocated(band_begin)) deallocate(band_begin)
         if (allocated(band_end)) deallocate(band_end)
         if (allocated(n_points_in_band)) deallocate(n_points_in_band)
         if (allocated(broaden_type)) deallocate(broaden_type)
         if (allocated(broaden_para)) deallocate(broaden_para)
         if (allocated(direc_1)) deallocate(direc_1)
         if (allocated(direc_2)) deallocate(direc_2)

         restart_relaxations  = .false.
         hessian_to_restart_geometry = .false.
         write_restart_geometry = .false.
         use_distributed_hessian = .false.

         final_forces_cleaned = .true.
         output_in_original_unit_cell = .true.

         grid_storage_initialized = .false.
         got_n_compute_maxes = .false.

         prune_basis_once = .true.

         use_local_index = .false.
         use_load_balancing = .false.
         use_alltoall = .false.

         restrict_kpoint_to_smp_node = .false.

         first_integration = .true.

         parallel_grid_partitioning = .false.

         force_mpi_virtual_topo = .false.

         use_metis_batch_distribution = .false.

         frozen_core_scf = .false.
         frozen_core_scf_cutoff = -5.d2
         frozen_core_scf_factor = -1.d0

         recompute_batches_in_relaxation = .true.

         out_vacuum_potential = .false.
         out_embedding_potential = .false.
         flag_set_vacuum_level = .false.
         use_dipole_correction = .false.
         calculate_work_function = .false.
         dipole_correction_method='potential'

         dos_n_en_points = 500
         dos_low_energy  = -20d0
         dos_high_energy = 20d0
         dos_alpha       = 3d-1
         atom_dos_n_en_points   = 500
         atom_dos_low_energy    = -20d0
         atom_dos_high_energy   = 20d0
         atom_dos_alpha         = 3d-1
         l_proj_dos_n_en_points = 500
         l_proj_dos_low_energy  = -20d0
         l_proj_dos_high_energy = 20d0
         l_proj_dos_alpha       = 3d-1

         MD_maxsteps         = -1
         MD_time             = 0d0
         MD_tstep            = 1d-3
         MD_temperature      = 0d0
         MD_init_temperature = 0d0
         MD_tau_berendsen    = 1d0
         MD_nu_andersen      = 1d0
         MD_tau_BDP          = 1d0
         MD_gle_ns           = 0
         MD_Q_NP             = 0d0
         MD_ensemble         = 'NVE'
         NVE_damping_factor  = 1d0
         MB_velocity_initialization = .false.
         MD_RNG_seeded       = .false.
         MD_restart_file     = 'aims_MD_restart.dat'
         MD_thermostat_units = 'amu*bohr^2'
         MD_initialconf_from_restart = .false.
         MD_time_restart     = .true.
         MD_restart_binary   = .false.
         check_MD_stop       = .true.
         plumed_plugin       = .false.
         plumed_file         = 'plumed.dat'
         plumed_new_plugin   = .false.
         PIMD_maxsteps         = -1
         PIMD_time             = 0d0
         PIMD_tstep            = 5d-4
         PIMD_temperature      = 0d0
         PIMD_init_temperature = 0d0
         PIMD_Q_NP             = 0d0
         PIMD_ensemble         = 'NVE'
         PIMD_RNG_seeded       = .false.
         PIMD_restart_file     = 'aims_PIMD_restart.dat'
         PIMD_thermostat_units = 'amu*bohr^2'
         PIMD_initialconf_from_restart = .false.
         PIMD_time_restart     = .true.
         PIMD_restart_binary   = .true.
         check_PIMD_stop       = .true.
         flag_energy_density                = .false.
         flag_chetty_martin_energy_density  = .false.
         flag_harris_foulkes_energy_density = .false.
         flag_project_energy_density        = .false.

         postprocess_anyway = PP_ANYWAY_NOTHING

         wf_extra_type = WF_EXTRA_NONE
         wf_extra_use_densmat = .false.

         out_zero_multipoles = .false.
         out_hartree_multipoles = .false.
         output_level = 'normal'

         basis_dep_cutoff_default    = 1d-4

         plot_band_allocated = .false.
         band_initialized = 0

         plot_dielectric_allocated = .false.
         dielectric_initialized = 0
         calc_px_dielectric = .false.
         calc_py_dielectric = .false.
         calc_pz_dielectric = .false.
         die_broaden_type  = 2
         die_broaden_para  = 0.1

         DFPT_sc_accuracy_dm  = 1.0d-3
         DFPT_mixing          = 0.2d0
         DFPT_width          = 0.0
         restart_write_DM1 = .false.
         restart_read_DM1  = .false.

         transport_lead_calculation = .false.
         transport_calculation  = .false.
         flag_output_effective_potential = .false.

         calculate_atom_bsse = .false.

         TDDFT_time                      = 0d0
         TDDFT_dt                        = 0.1d0
         TDDFT_propagator           = 'crank_nicholson'
         TDDFT_output_every              = 100
         TDDFT_write_info                = 1
         TDDFT_stochastic                = 0
         TDDFT_iseed                     = 12345
         flag_run_tddft_real_time_propagation = .false.

         flag_neutral_excitation = .false.
         neutral_excitation_rpa = .false.
         neutral_excitation_tdhf = .false.
         neutral_excitation_tddft = .false.
         excited_mode_singlet = .true.
         excited_mode_triplet = .true.
         num_excited_states = -1
         libxc_tddft_x = 0
         libxc_tddft_c = 0
         libxc_tddft_xc = 0
         libxc_tddft_chosen_x = .false.
         libxc_tddft_chosen_c = .false.
         casida_setup_scalapack = .false.

         neutral_excitation_bse = .false.
         read_qpe = .false.
         write_qpe = .false.
         write_ovlp_3ks = .false.
         write_ovlp_3fn = .false.
         write_screened_coulomb = .false.
         bse_singlet = .false.
         bse_triplet = .false.
         bse_lower_limit = 1
         bse_upper_limit = 25000
         bse_reduce_matrix = .false.
         bse_reduce_unocc = 3
         bse_reduce_occ = -3
         output_ks_coulomb_integral = .false.

         band_mulliken_orbit_num = 100

         dipole_trans_mat_out = .false.
         dipole_trans_mat_k = 1
         dipole_trans_mat_orbit_num = 20

         n_frozen_shell = 0
         flag_frozen_core_postSCF = .false.
         n_frozen_virtual_states = 0
         flag_frozen_virtual_postSCF = .false.
         en_shift_type       = 0
         flag_en_shift       = .false.
         en_shift_constant   = 0.0d0
         coupling_pt2_factor = 0.0d0
         coupling_pt2_screen = 0.0d0
         coupling_pt2_shift  = 0.0d0
         scsrpa_coeff = 1.0d0
         special_radius = 0.0d0
         renormalized_eg = 0.0d0
         c_osrpa_threshold   = 0.9d0
         n_lambda_osrpa      = 10
         ci_truncation             = 0
         flag_single_excitation    = .true.
         flag_single_spin_channel  = .false.
         flag_turn_off_hartree     = .false.
         ci_accuracy_etot          = 1.0d-6
         n_scf_ci                  = 100
         ci_coeff_threshold        = 1.0d-3
         n_ci_acc_save             = 10
         ci_acc_method             = 2
         ci_acc_ovlp_thresh    =   1.0E-5
         CC_calc_method             = 1
         CC_max_cycle               = 50
         CC_acc_method              = 2
         CC_n_RLE_sv                = 10
         CC_restart_point           = 0
         CC_restart_step            = 10
         CC_mkl_nth                 = 1
         CC_converg_thresh = 1.0D-6
         CC_RLE_bas_thresh = 1.0D-6
         CC_lv_shift       = 0.0D0
         CC_sv_method               = 1
         CC_n_domain                = 1
         CC_MPI_max_len             = 100000000
         CC_work_tmp_size           = 100000000
         CC_DIIS_n_sv               = 10
         CC_DIIS_step               = 1
         CC_DIIS_err_method         = 0
         CC_PT_n_out                = 5
         flag_cc_general            = .false.
         CC_collective              = .false.
         CC_use_disk                = .false.
         CC_restart_flag            = .false.
         CC_check_flag              = .false.
         CC_gamma_flag              = .false.
         CC_abcd_sv_flag            = .true.
         CC_mem_check               = .false.
         CC_super_system_flag       = .false.
         CC_dynamic_threading       = .true.
         CC_solver     = 'cluster'
         CC_sv_control = 'full'
         CC_method     = 'CCSD'
         rpa_along_ac_path_grid    = 0
         n_dft_methods = 0

         add_embedding_grids = .true.

         solvent_method                  = SOLVENT_UNDEF
         permittivity_switching_width    = 0.e0_dp
         mpe_nonel_model                 = MPE_CONST % NONEL_UNDEF
         mpe_nonel_model_bulk            = MPE_CONST % NONEL_UNDEF
         mpe_nonel_alpha                 = 0.e0_dp
         mpe_nonel_beta                  = 0.e0_dp
         mpe_nonel_alpha_if              = -1.e0_dp
         mpe_nonel_alpha_bulk            = 0.e0_dp
         mpe_nonel_beta_bulk             = 0.e0_dp
         mpe_factorization_type          = MPE_CONST % FACTZN_UNDEF
         mpe_n_boundary_conditions       = 2
         mpe_solvent_permittivity        = 1.0e0_dp
         mpe_bulk_permittivity           = 1.0e0_dp

         mpe_xml_loglevel                = MPE_CONST % XML_NOLOG
         mpe_xml_logfile        = 'mpe_interface.xml'
         mpe_skip_first_n_scf_steps      = -1
         mpe_degree_of_determination     = -1.e0_dp
         mpe_lmax_rf                     = -1
         mpe_lmax_ep                     = -1
         mpe_lmax_ep_from_basis          = .false.
         mpe_lmax_ep_additional_order    = 0
         mpe_n_centers_ep                = -1
         mpe_lmax_ep_increase_cluster    = .false.
         mpe_f_sparsity_threshold        = -1.e0_dp

         mpe_tol_adjR2                   = 0.075
         mpe_tol_adjR2_wait_scf          = .false.

         isc_cavity_type                 = ISC_CONST % CAVITY_UNDEF
         isc_cavity_restart_read         = .false.
         isc_cavity_restart_read_file  = 'cavity_restart.dat'
         isc_cavity_restart_write        = .false.
         isc_cavity_restart_write_file = 'cavity_restart.dat'
         isc_record_cavity_creation_file = 'cavity_creation.xyz'
         isc_record_every_n_steps        = 0
         ifp = .false.
         ifp_dist = 0.e0_dp
         ifp_normal = (/ 0.e0_dp,0.e0_dp,1.e0_dp /)
         isc_isodensity_value            = 1.e-3_dp * (bohr**3)
         isc_calculate_surface_and_volume = .true.
         isc_kill_ratio                  = 0.e0_dp
         isc_surface_curvature_correction = .true.
         isc_try_restore_convergence     = .false.
         isc_initial_radius                = -1.e0_dp
         isc_surface_sampling              = ISC_CONST%CAV_SAMP_Undef
         isc_dt                            = -1.e0_dp
         isc_rho_k                         = -1.e0_dp
         isc_rho_k_mod_constant            = -1.e0_dp
         isc_g_k                           = -1.e0_dp
         isc_rep_k                         = -1.e0_dp
         isc_rho_rel_deviation_threshold   = -1.e0_dp
         isc_gradient_threshold            = -1.e0_dp
         isc_dynamics_friction             = -1.e0_dp
         isc_update_nlist_interval         = -1
         isc_max_dyn_steps                 = -1
         ifp_rmin                          = -1.e0_dp
         ifp_rmax                          = -1.e0_dp
         ifp_manual                        = .false.
         ifp_rlog                          = -1.e0_dp
         ifp_area_per_point                = -1.e0_dp
         ifp_n_shells_log                  = -1
         ifp_n_angular_max                 = -1
         ifp_min_dr                        = -1.e0_dp
         ifp_n_shells_lin                  = -1

         onsite_accuracy_threshold = 0.030/hartree
         override_integration_accuracy = .true.

         use_gpu = .false.
         use_gpu_density = .false.
         use_gpu_hamiltonian = .false.
         use_gpu_forces = .false.

         flag_cpu_consistency = .true.
         cpu_consistency_threshold = 1d-11

         scalapack_block_size = -1

         out_matrices_elsi = .false.
         out_h_elsi = .false.
         out_s_elsi = .false.
         out_aims_json_log = .false.
         output_matrix_io_timing = .false.

         first_elsi_call = .true.
         use_elsi = .false.
         use_elsi_dm = .false.
         use_elsi_ev = .false.

         elsi_solver = 1

         elsi_out_level = 1
         elsi_out_json = 0
         elsi_read_dm = .false.
         elsi_read_dm_done = .false.
         elsi_write_dm = .false.
         elsi_write_dm_this = .false.
         elsi_write_dm_freq = 5
         do_elsi_rw = .false.
         elsi_extrap_dm = .false.
         elsi_rw_filter = 1.d-15
         elsi_filter = 1.d-12
         elsi_elpa_solver = 2
         elsi_elpa_n_single = 0
         elsi_elpa_gpu = 0
         elsi_omm_n_elpa = 6
         elsi_omm_flavor = 0
         elsi_omm_tol = 1d-12
         elsi_pexsi_np_per_pole = -1
         elsi_pexsi_np_symbo = 1
         elsi_eigenexa_method = 2
         elsi_sips_n_slice = -1
         elsi_sips_slice_type = 3
         elsi_sips_buffer = 0.02d0
         elsi_sips_n_elpa = 0
         elsi_ntpoly_method = 2
         elsi_ntpoly_tol = 1d-4
         elsi_ntpoly_filter = 1d-8
         elsi_magma_solver = 1

         need_dens_superpos = .false.
         recalc_dens = .false.

         potconspar = 0.0

         magnetic_response = .false.
         mr_experimental = .false.
         dfpt_iter_limit = 40
         dfpt_accuracy_n1 = 1d-9
         dfpt_linear_mix_param = 1d0
         dfpt_pulay_steps = 8
         use_dfpt_pulay = .true.
         output_sxml = .false.
         mr_gauge_origin = 1d300
         sxml_filename = 'aims.sxml'

         fo_fragment = .false.
         fo_density = .false.

         force_single_restartfile = .false.
         restart_eigenvectors_periodic = .false.

         tau_threshold = 0.d0

         use_rlsy_symmetry = .false.
         skip_rlsy_symmetry = .false.

      end subroutine set_runtime_choices_defaults

      end module runtime_choices
!******
