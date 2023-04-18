!****** FHI-aims/dimensions
!  NAME
!    dimensions
!  SYNOPSIS

module dimensions

!  PURPOSE
!
!  Module dimensions obtains all directly accessible array dimensions
!  for dynamic allocation by parsing the input files:
!  * control.in
!  * geometry.in
!
!  SINCE OTHER MODULES DEPEND ON THIS MODULE, THIS MODULE MUST NOT USE ANY OTHER MODULE!
!
!  Subroutines and functions:
!  * parse_control
!  * parse_species
!  * parse_geometry
!
!  FIXME: This module is rather rudimentary at present - it will grow as
!  more parts of the program become allocatable. At the same time,
!  all static array dimensions in param.f should vanish.
!
!  global variable declarations
!
!     n_species:  Total number of different atom types (species) in control.in
!     n_atoms:    Total number of atoms in geometry.in
!     n_centers:  Total number of centers in system (meaning in periodic systems)
!     n_k_points: Total number k points   in system (meaning in periodic systems)
!
!     Nadia
!     n_multipoles: Total number of multipole coordinates for a possible embedding potential, in geometry.in
!
!     n_periodic: Number of dimensions (0,1,2,3) in which we have periodic boundary conditions
!     n_spin:     Number of spin channels in the program: 1 for unpolarized, 2 for collinear spin-density functional theory
!     spin_degeneracy:  Maximum degeneracy of a spin channel.  2.0d0 for unpolarized (spin up and down degenerate), 1.0d0 for
!                       collinearly-polarized and non-spin-collinear calculations (example of later being spin-orbit coupling)
!     l_ext_max: Maximum angular momentum quantum number for extra basis functions for auxliriary basis in RI
!     l_wave_max: Maximum angular momentum quantum number for basis functions
!     n_wave_max: Maximum radial quantum number for occupied atomic / ionic
!                 states
!     n_max_ind_fns: Maximum number of distinct radial functions per species and type
!                 FIXME: This is a crudge and could be replaced entirely by using
!                 n_max_basis_fns and index arrays.
!     n_max_aux_fns: Maximum number of distinct extra radial functions for auxiliary basis
!                    per species and type
!     n_basis_types: Number of actual basis function types (atomic, ionic, confined, ...)
!                 used in the calculation. Needed to avoid a fixed dimension
!                 in shrink_*_basis, could be used to unify the various
!                 species basis fn arrays confined_*, hydro_*, ionic_*, etc.
!     n_max_contracted: [for Gaussian basis functions only] Number of elementary wave functions
!     n_max_aux_contracted: [for auxiliary Gaussian basis functions only]
!                  Number of elementary wave functions in the contracted functions.
!                 in a single contracted basis function
!     n_max_grid: Maximum number of points in the logarithmic grid
!                 of any one species (used for saving splines etc.)
!     n_max_radial: Maximum number of radial integration shells
!                 for any one species
!     n_max_angular: Maximum number of angular integration points
!                 in any one radial integration shell of any one species
!     n_max_angular_division: Maximum number of divisions for angular integration shell
!     l_pot_max : Maximum angular momentum quantum number needed for partitioned Hartree potential
!     n_max_spline : Maximum number of spline parameters at a given spline point.
!                 Currently only cubic splines are implemented, hence n_max_spline = 4 always.
!     n_channels : Auxiliary dimension variable to remove static dimensions from
!                 atomic solver (ancient fhipp part of program).
!                 Without spin, there is only one effective potential, hence one "channel."
!                 With spin, there would be two "channels." If the atomic effective potential
!                 were orbital-dependent (e.g. Hartree-Fock), we would need as many channels
!                 as there are orbitals for a given element.

!     n_max_basis_fns: Maximum total number of distinct radial basis functions
!                 according to control.in, before redundancy check.
!     n_max_basis: Maximum number of atomic basis states u_{at,n,l,m}(r)*Y_{l,m}(Omega)
!                 according to control.in, before redundancy check. Determined in read_geo() .
!     n_states : Total number of Kohn-Sham eigenstates to be considered
!                This variable might be reduced during an FHI-aims run if its initial value
!                is larger than the number of basis functions after reduction by the check
!                of ill-conditioning
!     n_states_init : Total number of Kohn-Sham eigenstates specified by the user
!                     This variable won't be modified during an FHI-aims run
!     n_states_k  : Total number of k-dependent Kohn-Sham eigenstates to be considered in periodic
!                   post-DFT calculations

!     n_basis_fns : Number of distinct radial basis functions after orthonormalisation and
!                 redundancy check
!     n_basis    : total number of basis functions in Hamiltonian matrix
!                  (for each radial fn, need 2l+1 m components) after orthonormalisation
!                 and redundancy check
!     n_max_pulay : Maximum number of stored charge density iterations for Pulay mixer
!     n_max_broyden : Maximum number of stored charge density iterations for Broyden mixer
!     relative_fp_charge_mix: Relative under-relaxation of fixed point part of
!                  quasi-Newton method.
!     n_full_points: Number of integration points in the grid known to each thread.
!     n_int_points: Actual number of integration points which are evaluated. We skip all
!                  points on which the integral partition function is zero.
!     n_spline_atoms : Number of atoms that each thread has the spline informaion on
!     n_cube      : Graphical output option - number of distinct requested .cube files
!     n_hartree_grid: Specifies the 1D radial grid used to calculate the Hartree potential
!                   either the regular radial integration grid, n_radial, or the free-atom
!                   logarithmic grid, n_grid
!                   FIXME: n_hartree_grid should be removed at some point in the future,
!                   leaving only the better option (n_radial or n_grid) behind.
!     n_max_coeff_3fn: maximal unit cells for the array coeff_3fn used in periodic Hartree-Fock calculations.

!     n_region     : number of different spatial regions for the locally-constrained DFT
!     n_ini_type : number of different initialization types for initial_rho due to different
!                      moments and charges

!     n_pp_species:  Total number of different pseudoized atom types (pp_species) in control.in

!     FORCED OCCUPATION

!     n_force_occ : number of states to follow by projection and force their occupation

!     BSSE
!     n_empty_atoms : number of empty atoms
!     n_occ_atoms: number of occupied atoms
!     n_occ_centers : number of occupied centers in PB
!     n_occ_centers_basis_integr : number of occupied basis integral centers in PB

!     External perturbation
!     ext_cycle: number of cycles where perturbation (Electric field) is applied
!
!
!  USES

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

      integer :: n_species
      integer :: n_atoms
      integer :: n_centers
      integer :: n_beads   ! XZL: added for PIMD

      integer :: n_empty_atoms
      integer :: n_occ_atoms
      integer :: n_occ_centers
      integer :: n_occ_centers_basis_integr
      integer :: n_occ_centers_integr

      ! See pbc_lists.f90.
      integer :: n_centers_hartree_potential
      integer :: n_centers_hartree_multipole
      integer :: n_centers_ele_summation
      integer :: n_centers_basis_integrals
      integer :: n_centers_integrals
      integer :: n_centers_in_hamiltonian  ! can be used in 1D
      integer :: n_centers_in_sc_DFPT      ! can be used in 1D ~ 3D

      ! JW: There (supposedly) used to be two different lists of basis
      ! functions touching the 0-cell: All of them and only (about) half of
      ! them where the rest was equivalent by inversion symmetry.  Now, there
      ! is only the list with all of them.  I guess that parts of the code
      ! might signal something by using one or the other, so I did not want to
      ! replace them by a unified version.

      ! Number of basis functions in super cell:
      ! (if in doubt, use this):
      integer :: n_Cbasis
      ! Number of basis functions in super cell
      ! (used to explicitly request ALL of them):
      integer :: n_centers_basis_T
      ! Number of basis functions in super cell
      ! (used to explicitly state that they are only needed modulo inversion,
      ! but do not rely on this in old code!):
      integer :: n_centers_basis_I

      integer :: n_k_points
      integer :: n_k_points_nosym
      integer :: n_k_points_task
      integer :: n_k_points_original
      integer :: n_ks_points_task
      integer :: n_irk_points ! # of irreducible k points
      integer :: n_irk_points_task
      integer :: n_max_coeff_3fn

      integer, allocatable :: ik2irred_map(:,:,:) ! mapping from nosym kpoints to reduced kpoints

      integer :: n_irkq_points_task !local # of irreducible k points for polarization matrix
      integer :: n_kq_points_task  !local # of k points for coulomb matrix
      integer :: lbb_row, ubb_row, lbb_col, ubb_col
      integer :: n_bb_row, n_bb_col, bb_bl_row, bb_bl_col

      integer :: n_band_kpoints  ! # of k points in band plotting
      integer :: n_band_kpoints_task

! SVL For periodic Hartree-Fock, one may need to have a different grid
      integer :: n_q_points
      integer :: n_q_points_task
!     Nadia
      integer :: n_multipoles

      integer :: vdw_pairs
      integer :: n_vdw_pairs_ignore

      integer :: cell_edge_steps(3)
      real*8 :: cell_edge_units(3)

      integer :: n_periodic

      integer :: n_spin
      real*8  :: spin_degeneracy

      integer :: n_region

      integer :: l_ext_max
      integer :: l_wave_max
      integer :: n_wave_max
      integer :: n_max_ind_fns
      integer :: n_max_aux_fns
      integer :: n_basis_types
      integer :: n_max_contracted
      integer :: n_max_aux_contracted
      integer :: n_max_grid
      integer :: n_max_radial
      integer :: n_max_angular
      integer :: n_max_angular_division

      integer :: l_pot_max
      integer :: n_max_spline
      integer :: n_channels

      integer :: n_max_basis_fns
      integer :: n_max_basis
      integer :: n_states
      integer :: n_states_init
      integer :: n_core_states
      integer :: n_valence_basis
      real*8 :: min_energy_include_in_soc
      logical :: min_energy_include_in_soc_set
      real*8 :: min_energy_save_in_soc
      logical :: min_energy_save_in_soc_set
      real*8 :: max_energy_include_in_soc
      logical :: max_energy_include_in_soc_set
      real*8 :: max_energy_save_in_soc
      logical :: max_energy_save_in_soc_set
      real*8 :: gap_for_min_energy_in_soc
      logical :: gap_for_min_energy_in_soc_set
      real*8 :: gap_for_saved_min_energy_in_soc
      logical :: gap_for_saved_min_energy_in_soc_set
!     >>> inserted by AB: march 2012
      integer :: n_states_save
!     <<< done with insert, AB: march 2012
      integer :: n_states_sc_DFPT

      integer :: n_basis_fns
      integer :: n_ext_fns
      integer :: n_basis
      integer :: n_ext
      integer :: n_hamiltonian_matrix_size
      integer :: n_hamiltonian_matrix_size_OLD
      integer :: n_hamiltonian_matrix_size_no_symmetry
      integer :: n_basis_supercell  ! can be used in 1D
      integer :: n_basis_sc_DFPT    ! can be used in 1D ~ 3D

      integer :: n_max_pulay
      integer :: n_max_pulay_constraint
      integer :: n_max_broyden
      real*8 :: relative_fp_charge_mix
      integer :: n_int_points
      integer :: n_int_points_total
      integer :: n_full_points
      integer :: n_full_points_total
      integer :: n_full_points_hamiltonian_integrals

      integer :: n_spline_atoms

      integer :: n_cube
      integer :: n_esp
      integer :: n_plot_band
      integer :: n_plot_band_scf
      integer :: n_plot_band_DFPT_phonon

      integer :: n_plot_dielectric
      integer :: n_plot_absorption

      integer :: n_write_soc_eigenvectors
      ! CC: Z2 by CC / Adapted from Cmera:
      integer      :: Z2_n_plot     ! Number of points in the Z2 invariant calculation
      integer      :: Z2_n_k_points ! k-points for each Wannier center
      integer      :: Z2_n_planes   ! Wilson loop index

      integer :: n_hartree_grid

      integer :: ext_cycle

!  Infrastructure to reduce dimensions of temporary arrays
!  during integration, density update etc on a grid to the minimum:

      integer :: n_max_points_per_div

      integer :: n_max_compute_atoms
      integer :: n_max_compute_ham
      integer :: n_max_compute_fns_ham
      integer :: n_max_compute_dens
      integer :: n_max_compute_fns_dens
      integer :: n_max_compute_ang
      real*8 :: n_avg_compute_dens
      integer :: n_max_compute_missing_dens
      integer :: n_max_compute_fns_missing_dens
      integer :: n_max_compute_missing_atoms
      real*8 ::  n_avg_compute_missing_dens

!  Dimension variables for integration over arbitrary batches
      integer :: n_grid_batches
      integer :: n_points_in_batch
      integer :: n_max_batch_size
      integer :: n_my_batches

      integer :: n_ini_type

!  Infrastructure for new external grids

      integer :: n_max_ang_shells

!  Maximal number of (n,l) shells per species for which LDA(GGA)+U is requested

      integer :: n_max_shells_plus_u

      integer :: n_hess_blocks
      integer :: n_hess_blocks_lv
      integer :: n_hess_blocks_lv_atom
      logical       :: hess_in_file

!     FORCED OCCUPATION
      integer :: n_force_occ

!     Boys localization
      integer :: n_sub_boys

      logical :: flag_out_ev_scalapack_hdf5

!     output the atom electrostati potential
      logical :: flag_out_elec_real
      logical :: flag_out_locpot_atom
!     Momentum matrix and linear dielctric tensor and localized orbitals
      logical :: flag_out_dielectric
      logical :: flag_dielectric_test
      logical :: flag_broadening
      logical :: flag_out_absorption
      logical :: flag_out_greenwood
      logical :: flag_greenwood_method
      logical :: flag_out_dipmat
      logical :: flag_out_dipmat_k_k
      logical :: flag_out_coulmat_ovl
      logical :: flag_out_coulmat_lvl
      logical :: flag_out_boys
!     Variables associated with the Kubo-Greenwood formalism for
!     electronic transport
      logical :: flag_out_dclimit
      logical :: flag_flex_fermi
      logical :: flag_explicit_fermi

      logical :: set_blacsdim, use_threadsafe_gwinit
      logical :: set_blacsdim_cpt2, set_pair_block_cpt2

      logical :: gamma_cut_coulomb
      logical :: use_gw_gamma_corr
      logical :: flag_adjust_qp_pre
      logical :: flag_adjust_qp_post

      !  The following flags are set to avoid unnecessary allocations and also unnecessary
!  code executions in case a run does not need them:
!    use_min_basis determines whether we need to deal with the minimal basis at all
!    use_confined determines whether we need confined basis function arrays
!    use_ionic determines whether we need ionic basis function arrays
!    use_hydro determines whether we need hydrogenic basis function arrays
!    use_gaussian determines whether we need gaussian basis function arrays
!    use_aux_gaussian determines whether we need auxiliary gaussian basis function
!                      arrays for RI.
!    use_ext_basis determines whether we need extra basis function
!                      arrays for RI.
!    use_basis_gradients determines whether we need the x,y,z gradients of
!                        basis functions, and the associated maths
!    use_density_gradient determines whether we need the gradient of the
!                        charge density, and the associated maths
!    use_hartree_log_grid     specifies whether the 1D Hartree multipole components are
!                        given on the logarithmic grid or on the radial integration grid
!    use_cube_output     specifies whether cube-style output is requested
!    use_out_eigenvec_ovlp    specifies whether unformatted files containing eigenvector and overlap matrix
!                             used for MODOS calculations are requested
!    use_bare_ci         specifies whether molecular bare coulomb Integral output is requested
!    use_distributed_spline_storage     specifies if splines arrays are to be distributed over MPI threads
!    use_embedding_potential  initially, intended for embedding of fixed charges into the
!                             external potential
!    use_vdw_correction       Empirical vdW C6 correction (fitted
!                             separately for each functional)
!    use_ll_vdwdf             calculating non_local correlation energy from
!                             Langreth-Lundqvist van der Waals density functional
!    use_nlcorr_post          Calculating nonlocal correlation energy from LL df (2010).
!    use_vdw_post             Post-process with either pbe_vdw or revpbe_vdw.
!    use_nlcorr_in_xc         flag to turn off nlcorr in xc functional.
!    use_qmmm                 QM/MM switch for calculating energy/forces
!                             on external multipoles
!    use_meta_gga             flag for meta-gga
!    use_meta_gga_post        flag for meta-gga post processing
!    use_meta_gga_printout    flag for printing meta-gga components
!    use_hybrid_meta          flag for post processing hybrid meta-gga's
!    force_occupation_projector : whether do force occupation or not
!    force_occupation_projector_redsub : whether we reduce the subspace for MOM during SCF or not
!    force_occupation_projector_auto : whether we automatically detect oscillations and adapt the subspace
!    start_force_occ  : needed to avoid occupation forcing during first scf setup before reading restart
!    calculate_perturbative_soc : logical variable
!     use_periodic_hf : periodic version of Hartree-Fock exchange; can be used also for cluster calculations
!     use_hf_realspace : real-space periodic Hartree-Fock implementation
!     use_hf_kspace :  k-space periodic Hartree-Fock implementation
!     use_hf_kspace_with_rpa : use k-space HF to calculate exchange matrix and energy in EX+RPA
!    use_el_ph_coupling : calculate electron phonon coupling matrix element for
!                         specified electron initial state <i1,k1>, final state
!                         <i2,k2> and accompanying phonon state <j,q>
!     use_vdw_method      using LL vdw DF (2010). SAG
!    use_embedding_pp         switch to include pseudopotentials

      logical :: use_min_basis
      logical :: use_confined
      logical :: use_ionic
      logical :: use_hydro
      logical :: use_gaussian
      logical :: use_sto = .false.
      logical :: use_aux_gaussian
      logical :: use_ext_basis
      logical :: use_sph_gaussian
      logical :: use_basis_gradients
      logical :: use_relativistic_basis
      logical :: use_density_gradient
      logical :: use_gga
      logical :: use_hse
      logical :: use_lc_wpbeh
      logical :: use_hartree_log_grid
      logical :: use_cube_output
      logical :: use_out_eigenvec_ovlp
      logical :: use_constraint
      logical :: use_hf_multiplicity
      logical :: use_angular_division
      logical :: use_forces
      logical :: use_numerical_stress
      logical :: use_AS_Jac_in_pulay
      logical :: compute_heat_flux

      logical :: use_distributed_spline_storage

      logical :: use_initial_rho
      logical :: use_initial_moment

      logical :: use_specified_grid

      logical :: use_mixer_threshold

      logical :: use_el_ph_coupling

!     Nadia
      logical :: use_embedding_potential

!     Empirical vdW correction
      logical :: use_vdw_correction
      logical :: use_vdw_correction_hirshfeld
      logical :: use_vdw_correction_hirshfeld_sc

!     Many-body-dispersion(MBD)
      logical :: use_mbd_old
      logical :: use_mbd_dev
      logical :: use_mbd_std
      logical :: use_libmbd
!     LL vdW density functional
      logical :: use_ll_vdwdf

!     LL vdw density functional, self-consistent, (2010).  SAG
      logical :: use_vdw_method
      logical :: use_nlcorr_post
      logical :: use_vdw_post
      logical :: use_nlcorr_in_xc

!     QM/MM switch
      logical :: use_qmmm

!     Meta-gga
      logical :: use_meta_gga
      logical :: use_meta_gga_post
      logical :: use_meta_gga_printout
      logical :: use_hybrid_meta

!     DFPT
      logical :: use_DFPT, use_partition_deriv, use_ASR,  &
                      use_DFPT_reduce_memory,  &
                      use_DFPT_polarizability, &
                      use_DFPT_polar_reduce_memory, &
                      use_DFPT_dielectric, &
                      use_DFPT_phonon_gamma,   &
                      use_DFPT_phonon,  &
                      use_DFPT_phonon_reduce_memory

!     friction
      logical :: use_friction

!     Hartree-Fock and GW (see prodbas.f90)
      integer :: n_basbas
      integer :: n_basbas_supercell
      integer :: n_basbas_fns
      integer :: n_loc_prodbas
      integer :: n_loc_prodbas_supercell
      integer :: n_max_loc_prodbas
      integer :: n_max_loc_prodbas_supercell
      integer, dimension(:), allocatable :: n_prodbas_per_proc

      logical :: use_hartree_fock
      logical :: use_periodic_hf
      logical :: use_lvl_fast
      logical :: use_hf_realspace
      logical :: use_hf_kspace
      logical :: use_hf_kspace_with_rpa
      logical :: flag_KS_eigenfunc_conjg
      logical :: use_pbe0
      logical :: use_b3lyp
      logical :: use_b1lyp
      logical :: use_pbesol0
      logical :: use_minnesota
      logical :: use_gw
      logical :: use_periodic_gw
      logical :: use_gw_expt
      logical :: use_ev_scgw
      logical :: use_ev_scgw0
      logical :: use_gw2ox
      logical :: use_gws2ox
      logical :: use_gwsoxw
      logical :: use_sosex_selfe
      logical :: use_sosex2w_selfe
      logical :: use_scgw
      logical :: use_scgw0
      logical :: use_gw_and_hse
      logical :: use_split_xc_gw
      logical :: flag_print_self_energy
      logical :: flag_calc_spectral_func
      logical :: specify_mu
      logical :: plot_self_energy
      integer :: n_loc_grid
      logical :: use_screx
      logical :: use_cohsex
      logical :: use_dmft_pbe0
      logical :: use_dmft_gw
      logical :: use_contour_def_gw

      logical :: use_qpe
      logical :: use_corr

      logical :: use_cutCb
      logical :: use_ers
      logical :: use_erfc

      logical :: use_logsbt
      logical :: use_lvl      ! see runtime_choices.f90 for RI_type

      logical, private :: flag_RI       ! private flag only - carries information about whether
                                              ! an explicit choice for RI_method was made between sub-
                                              ! routines parse_control and parse_geometry
      logical :: flag_default_lvl       ! if true and RI_method not explicitly set, will use LVL
      logical, private :: flag_neutral_excitation ! set true if neutral excitation calculations are required.
      logical :: flag_max_coeff_3fn ! if true and RI_method not explicitly set, will use LVL

      logical :: use_bse
      logical :: use_coulomb_integral_ks
!      MP2 et al.
      ! use_full_spectrum is a short flag that indicates whether the full spectrum of
      !     eigenvalues and -vectors should be returned. In principle, it'll be a simple grouping
      !     of the remaining ones.
      logical :: use_full_spectrum
      logical :: use_hf_post
      logical :: use_gw_energy
      logical :: use_mp2
      logical :: use_mp2_blacs
      logical :: use_os_mp2
      logical :: use_os_mp2_qpe
      logical :: use_os_rpa_qpe
      logical :: use_sic_rpa_qpe
      logical :: use_mp2sf
      logical :: use_rpa_ene
      logical :: use_rpa_plus_2ox
      logical :: use_rpa_plus_sosex
      logical :: use_C6_coef
      logical :: use_prodbas
      logical :: use_dftpt2
      logical :: use_dftpt2_qpe
      logical :: use_dftpt2_and_hse
      logical :: use_dftpt2_and_lrc

! periodic post-DFT (RPA, GW, etc.)
      integer, dimension (:), allocatable :: n_states_k

      ! FCI relavent flags
      logical :: use_ci

      integer, dimension (:), allocatable :: n_fc_shell

      ! coupled cluster flags
      logical :: use_cc

      logical :: use_geo_relaxation
      logical :: use_relaxation_constraints
      integer :: relax_unit_cell     ! Enumerate options for cell relaxation:
                                     !    0 ->  no cell relaxation
                                     !    1 -> 'full'
                                     !    2 -> 'shape'

!======== RRS-PBC scheme, igor ===============
! rrs_pbc_center_atom          :: an index list of species, which are selected as the
!                                 center unit cell.
!
! rrs_pbc_center_atom(<1,2,3>,i_center) contains all info. of i_center
!             here "i_center" is the index of the atom in the center unit cell,
!                   1 for the index of i_center
!                   2 for the start basis index of i_center
!                   3 for the basis number of i_center
!
! rrs_pbc_equal_atom           :: an index matrix of species, which are the equal
!                                 atoms for the center unit cell.
!
! rrs_pbc_equal_atom(<1,2>,i_equal,i_center) contains all equal atoms of i_center
!            here "i_center" is the index of the atom in the center unit cell,
!                 "i_equal" is the equal atom associated with i_center
!                  1 for the index of i_equal
!                  2 for the start basis index of i_equal
!
! rrs_pbc_n_center_atom        :: the atom number in the center unit cell
! rrs_pbc_n_equal              :: the max number of equal points in RRS-PBC projection
! rrs_pbc_n_center_basis       :: the basis number about the center unit cell
! rrs_pbc_n_electron           :: the electron number in the unit cell
! rrs_pbc_n_hamiltonian_k      :: the length of the unit cell matrix,
!                                 rrs_pbc_hamiltonian_k
! rrs_pbc_n_k_points           :: total number of k points
! rrs_pbc_n_k_points_xyz       :: the k point number along three directions
! rrs_pbc_k_point_list         :: the k point list in the representation of
!                                 recip_lattice_vector
! rrs_pbc_n_plot_band          :: the band number we want to plot
!
! rrs_pbc_k_point_list(i_k_point,<1,2,3>)
!
      logical :: use_rrs_pbc
      integer, dimension (:,:), allocatable   :: rrs_pbc_center_atom
      integer, dimension (:,:,:), allocatable :: rrs_pbc_equal_atom
      integer                                 :: rrs_pbc_n_center_atom
      integer                                 :: rrs_pbc_n_center_basis
      integer                                 :: rrs_pbc_n_equal
      integer                                 :: rrs_pbc_n_hamiltonian_k
      integer                                 :: rrs_pbc_n_lattice_vector
      real*8, dimension(2)                         :: rrs_pbc_n_electron
      integer, dimension(2)                        :: rrs_pbc_n_electron_int
      integer                                 :: rrs_pbc_n_k_points
      integer,dimension(3)                         :: rrs_pbc_n_k_points_xyz
      real*8,dimension(:,:), allocatable      :: rrs_pbc_k_point_list
      ! some band plotting data, where and how many bands ...
      integer                                 :: rrs_pbc_n_plot_band
      real*8, dimension(:,:), allocatable     :: rrs_pbc_band_begin
      real*8, dimension(:,:), allocatable     :: rrs_pbc_band_end
      integer, dimension(:), allocatable      :: rrs_pbc_n_points_in_band
      real*8,dimension(:,:),allocatable       :: rrs_pbc_occ_num

!  MD and its features
      logical :: use_molecular_dynamics
      ! flags for the MD schedule
      integer      :: MD_segments
      logical      :: MD_use_schedule
      logical      :: MD_velocity_in_input

!  XZL: added for PIMD
      logical :: use_pathint
      logical :: use_nm
      logical :: use_pimd
      logical :: use_cmd
      logical :: use_rpmd
      ! flags for the PIMD schedule
      integer      :: PIMD_segments
      logical      :: PIMD_use_schedule
      logical      :: PIMD_velocity_in_input
! MR: PIMD wrapper
       logical :: use_pimd_wrapper
       character*40  :: pimdwrapper_server
       character*40  :: pimdwrapper_port

       character*40  :: calc_dens_superpos
! Add a constraining potential (not a bias potential!)
      logical :: use_potential_constrain

! TDI and its features
      logical :: use_thermodynamic_integration
      ! flags for the TDI in the MD_schedule environment
      integer      :: TDI_segments

! Reference free Adiabatic Switching (AS)
      logical :: use_reffree_AS
      ! in the present implementation it is a variation of TDI
      integer       :: AS_segments

! Do MD in harmonic potential only
      logical :: use_harmonic_pot_only

! MD_QH_init:
      logical :: use_MD_QH_init
      integer       :: MD_QH_init_segments

!      Atomic Reference
      logical :: use_atom_ref

!     LDA(GGA)+U
      logical :: use_plus_u
      logical :: plus_u_mulliken_charges
      logical :: plus_u_full

!     molecular bare Coulomb integrals
      logical :: use_bare_ci

!  private variables, for intramodule consistency

!     log_r_min, log_r_max, log_inc_min are only needed if we have external integration grids;
!     in that case, we need a waterproof worst-case dimension for the log grid
!     which is potentially too large

      real*8, private :: log_r_min
      real*8, private :: log_r_max
      real*8, private :: log_inc_min
      logical, private :: first_species

      logical :: force_occupation_projector
      logical :: force_occupation_projector_redsub
      logical :: force_occupation_projector_auto
      logical :: force_occupation_basis
      logical :: force_n_electrons
      real*8  :: forced_n_electrons
      logical :: start_force_occ

      logical :: force_occupation_smearing

      logical :: apply_boys_flag

      ! perturbative spin-orbit coupling stuff
      ! By assumption, if only calculate_perturbative_soc is set, then
      ! it is the non-self-consistent version that is requested
      ! These variables really belong in runtime_choices...
      logical :: calculate_perturbative_soc
      logical :: save_soc_perturbed_eigenvectors
      logical :: include_pw_lda_in_v_soc ! Include the PW-LDA XC potential evaluated at the SCF electron density in the
                                         ! effective potential used for v_soc.  By default, the XC potential is omitted
                                         ! from the effective potential.  This flag can be used regardless of the
                                         ! XC functional chosen by the user for the SCF cycle, but it is only supported
                                         ! for non-spin-polarized calculations.  This may become a default once more
                                         ! testing is done, but right now it is an undocumented feature.  Early tests
                                         ! indicate that this changes spin-orbit splittings by 1%.
      logical :: calculate_all_eigenstates ! This is similar to use_full_spectrum (sets n_states = n_basis) but does not
                                           ! assume that the user is using a correlated method.

      ! ESP charges
      logical :: use_esp
      integer :: esp_constraint
!     switch to include pseudopotentials
      logical :: use_embedding_pp                              !global and very important flag for wrappers

      character*40, dimension(:), allocatable :: pp_path_string            !pp_path read in parse_species

      integer :: n_pp_species
      integer :: n_pp_basis_fns
      integer :: n_max_pp_fn                              ! maximum number of different radial parts per pseudoized species
      integer :: n_max_points_pp_fn                 ! maximum number of grid points for pseudoized species data

      logical :: use_nonlinear_core                          !use nonlinear core-valence interaction
      logical,dimension(:), allocatable :: species_nonlinear_core         !array which saves whether species is a pseudospaces



      integer :: n_pp_atoms
      logical :: pp_path_exist


      integer :: n_real_atoms
      integer :: n_pp_in_qm
      integer :: n_mp_in_qm

      logical :: use_libxc
      logical :: use_libxc_tddft

      logical :: casida_reduce_matrix
      real*8 :: casida_occ_limit
      real*8 :: casida_unocc_limit

      ! SPGlib variables
      real*8  :: sym_precision
      logical :: use_spglib
      logical :: set_symmetry_analysis
      logical :: use_symmetry_analysis
      logical :: use_symmetric_forces

      logical :: flag_force_single_restartfile
      logical :: flag_restart_eigenvectors_periodic


      ! MOL: symmetry-constrained relaxation
      logical :: use_symm_const_geo_relaxation
      character*10, dimension(:), allocatable    :: SCR_params
      character*100, dimension(:,:), allocatable  :: SCR_lv_str
      character*100, dimension(:,:), allocatable  :: SCR_coord_str
      integer  :: SCR_n_params
      integer  :: SCR_n_params_coords
      integer  :: SCR_n_params_lv

      logical :: flag_split_atoms
      logical :: flag_split_min_val
      logical :: flag_split_max_val
      logical :: flag_split_batch
      integer :: split_min_val
      integer :: split_max_val
      real*8  :: split_batch

   contains

      subroutine set_dimensions_defaults ()

         n_k_points = 1
         n_k_points_nosym = 1
         n_k_points_task = 1
         n_k_points_original = 1
         n_ks_points_task = 1
         n_irk_points = 1
         n_irk_points_task = 1
         n_max_coeff_3fn =16

         n_irkq_points_task = 1
         n_kq_points_task = 1

         n_band_kpoints = 1
         n_band_kpoints_task = 1

         n_q_points = 1
         n_q_points_task = 1

         n_spin = 1
         spin_degeneracy = 2.0d0

         n_max_spline = 4
         n_channels = 1

         min_energy_include_in_soc_set = .false.
         min_energy_save_in_soc_set = .false.
         max_energy_include_in_soc_set = .false.
         max_energy_save_in_soc_set = .false.
         gap_for_min_energy_in_soc_set = .false.
         gap_for_saved_min_energy_in_soc_set = .false.

         hess_in_file = .false.
         flag_out_elec_real = .false.
         flag_out_locpot_atom = .false.
         flag_dielectric_test = .false.
         flag_broadening = .false.

         gamma_cut_coulomb = .false.
         use_gw_gamma_corr = .false.

         use_out_eigenvec_ovlp = .false.
         use_AS_Jac_in_pulay = .true.

         use_distributed_spline_storage = .false.

         use_mbd_dev = .false.
         use_mbd_std = .false.
         use_libmbd = .false.

         flag_print_self_energy = .false.
         flag_calc_spectral_func = .false.

         use_logsbt = .true.
         use_lvl = .false.

         flag_RI = .false.
         flag_default_lvl = .false.
         flag_neutral_excitation = .false.
         flag_max_coeff_3fn = .false.

         relax_unit_cell = 0

         use_rrs_pbc = .false.
         rrs_pbc_n_plot_band = 0

         MD_segments = 0
         MD_use_schedule = .false.
         MD_velocity_in_input = .false.

         PIMD_segments = 0
         PIMD_use_schedule = .false.
         PIMD_velocity_in_input = .false.

         use_potential_constrain = .false.

         TDI_segments = 0

         use_reffree_AS = .false.
         AS_segments = 0

         use_harmonic_pot_only = .false.

         use_MD_QH_init      = .false.
         MD_QH_init_segments = 0

         first_species = .true.

         force_occupation_projector = .false.
         force_occupation_projector_redsub = .false.
         force_occupation_projector_auto = .false.
         force_occupation_basis = .false.
         force_n_electrons = .false.
         start_force_occ = .false.

         force_occupation_smearing = .false.

         apply_boys_flag = .false.

         calculate_perturbative_soc = .false.
         save_soc_perturbed_eigenvectors = .false.
         include_pw_lda_in_v_soc = .false.
         calculate_all_eigenstates = .false.

         pp_path_exist = .false.

         use_libxc = .false.
         use_libxc_tddft = .false.

         casida_reduce_matrix = .false.
         casida_occ_limit = 0d0
         casida_unocc_limit = 0d0

         sym_precision         = 1d-5
         use_spglib            = .false.
         set_symmetry_analysis = .false.
         use_symmetry_analysis = .true.
         use_symmetric_forces  = .false.

         flag_force_single_restartfile = .false.
         flag_restart_eigenvectors_periodic = .false.

         use_symm_const_geo_relaxation = .false.

         flag_split_atoms = .true.
         flag_split_min_val = .false.
         flag_split_max_val = .false.
         flag_split_batch = .false.
         split_min_val = 1

      end subroutine

!******
!---------------------------------------------------------------------
!****s* dimensions/parse_control
!  NAME
!   parse_control
!  SYNOPSIS

      subroutine parse_control ( )

!  PURPOSE
!  Subroutine parse_control parses the input file control.in to
!  determine the sizes of allocatable arrays, and does a few
!  integrity checks.
!
!
!  FIXME: Must make sure that species defaults, as given in read_species.f,
!  actually work, even if a given quantity is not defined in
!  control.in. Right now this does not work.
!

!  USES
  use localorb_io
  use mpi_tasks, only: aims_stop, aims_stop_coll, myid, use_mpi_in_place
  use xc_f03_lib_m

!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE





         implicit none



!  imported variables

!  local variables

         integer :: i_code
         character*40 :: desc_str
         logical :: inputlogical
         character*250 :: inputline

         character*100 :: info_str

         logical :: pseudoized = .false.

         logical :: verbatim_writeout            ! Governs if the control file should be written to output or not
         logical :: coming_from_parse_species    ! Addition to verbatim_writeout - prevents the double writing of
                                                 ! a line after a species description is finished

         ! So we can work out what type of XC functional we have with libxc. AJL/Aug 2016
         ! It is also reused for the same purpose by the dfauto infrastructure.
         ! We can possibly rename this to xc_family, etc. JH/Sep 2016
         character*40 :: libxc_x_str, libxc_c_str
         integer :: libxc_index_for_split = 0
         integer :: libxc_family_x = 0
         integer :: libxc_family_c = 0
         integer :: libxc_family = 0
         integer :: libxc_family_post = 0
!  counters

         integer :: i_species = 0
         integer :: i_pp_species = 0

         character(*), parameter :: func = 'parse_control'

!  begin work

         call localorb_info(" ", &
               use_unit,'(2X,A)' )

         call localorb_info("-----------------------------------------------------------------------", &
               use_unit,'(2X,A)' )

         call localorb_info("Parsing control.in (first pass over file, find array dimensions only).", &
               use_unit,'(2X,A)' )

         call localorb_info("The contents of control.in will be repeated verbatim below", &
               use_unit,'(2X,A)' )

         call localorb_info("unless switched off by setting 'verbatim_writeout .false.' .", &
               use_unit,'(2X,A)' )

         call localorb_info("in the first line of control.in .", &
               use_unit,'(2X,A)' )

         call localorb_info("-----------------------------------------------------------------------", &
               use_unit,'(2X,A)' )

         call localorb_info(" ", &
               use_unit,'(2X,A)' )




         i_species = 0
         n_max_grid = 0
         n_max_radial = 0
         n_max_angular = 0

         n_wave_max = 0
         l_wave_max = 0
         l_ext_max = 0

         l_pot_max = 0

         n_max_ind_fns = 0
         n_max_aux_fns = 0
         n_basis_types = 0
         n_max_contracted = 0
         n_max_aux_contracted = 0

         n_max_pulay = 8
         n_max_pulay_constraint = 0

         n_max_broyden = 0
         relative_fp_charge_mix = 0.d0

         n_cube = 0
         n_esp = 0
         n_plot_band = 0
         n_plot_band_scf = 0
         n_plot_band_DFPT_phonon = 0

         n_plot_dielectric = 0
         n_plot_absorption = 0

         n_region = 0

         n_force_occ = 0
         n_sub_boys = 0

         use_min_basis = .false.
         use_ionic = .false.
         use_confined = .false.
         use_hydro = .false.
         use_gaussian = .false.
         use_aux_gaussian = .false.
         use_ext_basis = .false.
         use_sph_gaussian = .false.

         use_basis_gradients = .false.
         use_relativistic_basis = .false.
         use_density_gradient = .false.
         use_gga = .false.
         use_forces = .false.
         use_numerical_stress = .false.
         compute_heat_flux     = .false.
         use_hartree_log_grid = .true.


         use_cube_output = .false.

         use_esp = .false.
         esp_constraint  = 0

         use_constraint = .false.

         use_hf_multiplicity = .false.

         use_angular_division = .true.
         n_max_angular_division = 8

         use_geo_relaxation = .false.
         use_molecular_dynamics = .false.

         use_pathint = .false.  ! XZL: added for PIMD
         use_nm = .false.       ! XZL: added for PIMD
         use_pimd = .false.     ! XZL: added for PIMD
         use_cmd = .false.      ! XZL: added for PIMD
         use_rpmd = .false.     ! XZL: added for PIMD

         use_DFPT=.false.
         use_partition_deriv= .false.
         use_ASR=.false.
         use_DFPT_reduce_memory= .false.
         use_DFPT_polarizability=.false.
         use_DFPT_polar_reduce_memory=.false.
         use_DFPT_dielectric=.false.
         use_DFPT_phonon_gamma=.false.
         use_DFPT_phonon=.false.
         use_DFPT_phonon_reduce_memory=.false.

         use_friction = .false.

         use_hartree_fock = .false.
         use_periodic_hf = .false.
         use_lvl_fast = .false.
         use_hf_realspace = .false.
         use_hf_kspace = .false.
         use_hf_kspace_with_rpa = .false.
         use_pbe0 = .false.
         use_pbesol0 = .false.
         use_b3lyp = .false.
         use_b1lyp = .false.
         use_minnesota = .false.
         use_gw = .false.
         use_gw_expt = .false.
         use_periodic_gw = .false.
         use_gw_and_hse = .false.
         use_ev_scgw = .false.
         use_ev_scgw0 = .false.
         use_gw2ox = .false.
         use_gws2ox = .false.
         use_gwsoxw = .false.
         use_contour_def_gw = .false.
! add bse
         use_bse = .false.
         use_coulomb_integral_ks = .false.
         use_sosex_selfe = .false.
         use_sosex2w_selfe = .false.
         use_split_xc_gw = .false.
         use_scgw = .false.
         use_scgw0 = .false.
         specify_mu = .false.
         use_screx = .false.
         use_cohsex = .false.
         use_mp2sf = .false.

         use_qpe = .false.
         use_corr = .false.
         use_dmft_pbe0 = .false.
         use_dmft_gw = .false.
         use_cutCb = .false.
         use_ers = .false.
         use_erfc = .false.

         use_full_spectrum = .false.
         use_hf_post = .false.
         use_hse = .false.
         use_mp2 = .false.
         use_mp2_blacs = .false.
         use_os_mp2 = .false.
         use_os_mp2_qpe = .false.
         use_os_rpa_qpe = .false.
         use_sic_rpa_qpe = .false.
         use_ci  = .false.
         use_dftpt2 = .false.
         use_dftpt2_qpe = .false.
         use_dftpt2_and_lrc = .false.
         use_dftpt2_and_hse = .false.
         use_rpa_ene = .false.
         use_gw_energy = .false.
         use_rpa_plus_2ox = .false.
         use_rpa_plus_sosex = .false.
         use_C6_coef = .false.
         use_prodbas = .false.

         use_atom_ref= .false.
         ext_cycle=0

         use_initial_rho = .false.
         use_initial_moment = .false.

         use_bare_ci = .false.

         use_specified_grid = .false.
         n_max_ang_shells = 0

         use_mixer_threshold = .false.

         use_el_ph_coupling = .false.

         flag_adjust_qp_pre = .false.
         flag_adjust_qp_post = .false.

         ! pseudopot infrastructure. by default, we set use_embedding_pp FALSE.
         ! However, the actual and final decision is made in read_geo.f90.
         use_embedding_pp = .false.
         i_pp_species = 0
         use_nonlinear_core = .false.
         n_max_pp_fn = 0
         n_max_points_pp_fn = 0
         n_pp_basis_fns = 0

         ! initialize
         use_vdw_correction = .false.
         use_vdw_correction_hirshfeld = .false.
         use_vdw_correction_hirshfeld_sc = .false.
         n_vdw_pairs_ignore = 0
         use_ll_vdwdf = .false.
         use_nlcorr_post = .false.
         use_nlcorr_in_xc = .true.
         use_vdw_post = .false.
         use_qmmm = .false.
         use_meta_gga = .false.
         use_meta_gga_post = .false.
         use_meta_gga_printout = .false.
         use_hybrid_meta=.false.
         use_plus_u = .false.
         plus_u_mulliken_charges = .false.
         plus_u_full = .false.
         n_max_shells_plus_u = 0
         !PIMD wrapper variables
         use_pimd_wrapper=.false.
         calc_dens_superpos='none'
         pimdwrapper_server='none'
         pimdwrapper_port='none'
         verbatim_writeout = .true.
         coming_from_parse_species = .false.

         use_mpi_in_place = .true.

         !Momentum matrix and linear dielectric tensor and localized orbitals
         flag_out_dielectric = .false.
         flag_out_absorption = .false.
         flag_out_dipmat = .false.
         flag_out_dipmat_k_k = .false.
         flag_out_coulmat_ovl = .false.
         flag_out_coulmat_lvl = .false.
         flag_out_boys = .false.
         ! Variables for Kubo-Greenwood
         flag_out_dclimit= .false.
         flag_flex_fermi= .false.
         flag_explicit_fermi= .false.

         set_blacsdim = .false.
         use_threadsafe_gwinit = .false.

         set_blacsdim_cpt2 = .false.
         set_pair_block_cpt2 = .false.

         open(777,file="control.in",status='OLD',iostat=i_code)

         if (i_code/=0) then
            call aims_stop("* Input file control.in not found.", func)
         endif

         lineloop: do
            read(777,'(A)',iostat=i_code) inputline

            if(i_code<0) exit lineloop      ! end of file reached
            if(i_code>0) then
               call aims_stop("Unknown error reading file 'control.in'...", '')
            endif

            if (verbatim_writeout.and.(.not.coming_from_parse_species) ) then
               ! CC: Commands that span multiple lines have to be
               !      printed out explicitly in the checks below, e.g.
               !      MD_schedule
               call localorb_info(inputline,use_unit,'(2X,A)') ! Write the control file to the output line by line
                                                        ! ... but only if we didn't already write that same line while parsing a species description
            end if
            coming_from_parse_species = .false.

            read(inputline,*,iostat=i_code) desc_str
            if (i_code /= 0) then
               cycle ! empty line
            elseif (desc_str(1:1).eq.'#') then
               cycle ! comment
            elseif (desc_str.eq."verbatim_writeout") then
                  read(inputline,*,end=88,err=99) desc_str, verbatim_writeout
            elseif (desc_str.eq."n_max_pulay") then
               read(inputline,*,end=88,err=99) desc_str, n_max_pulay
            elseif (desc_str.eq."n_max_broyden") then
               read(inputline,*,end=88,err=99) desc_str, n_max_broyden
            elseif (desc_str.eq."relative_fp_charge_mix") then
               read(inputline,*,end=88,err=99) desc_str, relative_fp_charge_mix
            elseif (desc_str.eq."mixer_threshold") then
               use_mixer_threshold = .true.
            elseif (desc_str.eq."n_max_pulay_constraint") then
               read(inputline,*,end=88,err=99) desc_str, n_max_pulay_constraint
            elseif (desc_str.eq."spin") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if (desc_str.eq."collinear") then
                  n_spin = 2
                  spin_degeneracy = 1.0d0
                  use_initial_rho = .true.
               endif

            elseif (desc_str.eq."compute_forces") then
               read(inputline,*,end=88,err=99) desc_str, inputlogical
               if (inputlogical) then
                  use_basis_gradients = .true.
                  use_forces = .true.
               endif
            elseif (desc_str.eq."compute_numerical_stress") then
               read(inputline,*,end=88,err=99) desc_str, inputlogical
               if (inputlogical) then
                  use_numerical_stress = .true.
               endif
            elseif (desc_str.eq."compute_heat_flux") then
               read(inputline,*,end=88,err=99) desc_str, inputlogical
               if (inputlogical) then
                  compute_heat_flux     = .true.
                  use_basis_gradients   = .true.
               endif
            elseif (desc_str.eq."relativistic") then
               read(inputline,*) desc_str, desc_str
               if (desc_str.eq."zora") then
!             require wave function gradients in this case
                  use_basis_gradients = .true.
!             check do we use relativistic basis functions
                  read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str
                  if (desc_str.eq.'spinor') then
                     use_relativistic_basis = .true.
                  endif
               elseif (desc_str.eq."atomic_zora") then
!             require wave function gradients in this case
                  use_basis_gradients = .true.
!             check do we use relativistic basis functions
                  read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str
                  if (desc_str.eq.'spinor') then
                     use_relativistic_basis = .true.
                  endif
               elseif (desc_str.eq."kolning_harmon") then
!             require wave function gradients in this case
                  use_basis_gradients = .true.
                  ! use_second_density_gradient =.true.
               elseif (desc_str.eq."own") then
!             require wave function gradients in this case
                  use_basis_gradients = .true.
               elseif (desc_str.eq."x2c".or.desc_str.eq."4c_dks") then
!             require wave function gradients in this case
!             require relativistic basis functions
                  use_basis_gradients = .true.
                  use_relativistic_basis = .true.
               endif
            elseif (desc_str.eq."xc") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if (desc_str.eq."libxc") then
                   use_libxc = .true.
                   ! We need to work out the class of functional for initialisations
                   read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str
                   libxc_index_for_split = SCAN(desc_str,'+')
                   if (libxc_index_for_split.gt.0) then
                     libxc_x_str = desc_str(1:libxc_index_for_split-1)
                     libxc_c_str = desc_str(libxc_index_for_split+1:)
                   else
                     libxc_x_str = desc_str
                     libxc_c_str = ""
                   endif
                   libxc_family_x = xc_f03_family_from_id(xc_f03_functional_get_number(libxc_x_str))
                   libxc_family_c = xc_f03_family_from_id(xc_f03_functional_get_number(libxc_c_str))
                   ! I'm not going to worry as to whether this is valid - we'll do this in
                   ! read_control.f90 - so for now we'll initialise assuming its all good.
                   if (libxc_family_x.gt.libxc_family_c) then
                     libxc_family = libxc_family_x
                   else
                     libxc_family = libxc_family_c
                   endif
               else if (desc_str == "dfauto") then
                   read (inputline, *, end=88, err=99) desc_str, desc_str, desc_str
                   select case (str_lower(desc_str))
                   case ('pbe')
                       libxc_family = XC_FAMILY_GGA
                   case ('tpss', 'scan')
                       libxc_family = XC_FAMILY_MGGA
                   case ('pbe0')
                       libxc_family = XC_FAMILY_HYB_GGA
                   case ('scan0')
                       libxc_family = XC_FAMILY_HYB_MGGA
                   end select
               end if
               if (desc_str.eq."pbe" .or. &
                   desc_str.eq."r48pbe" .or. &
                   desc_str.eq."rpbe" .or. &
                   desc_str.eq."revpbe" .or. &
                   desc_str.eq."revpbe_vdw" .or. &
                   desc_str.eq."pbe_vdw" .or. &
                   desc_str.eq."pw91_gga" .or. &
                   desc_str.eq."blyp" .or. &
                   desc_str.eq."am05" .or. &
                   desc_str.eq."pbesol" .or. desc_str.eq."PBESOL" .or. &
                   desc_str.eq."pbeint" .or. &
                   desc_str.eq."xpbe" .or. &
                   libxc_family.eq.XC_FAMILY_GGA) then
                  ! require the density gradient
                  use_density_gradient = .true.
                  ! require gradient terms of the XC potential
                  use_gga = .true.
                  use_basis_gradients = .true.
               elseif (desc_str.eq."m06-l".or.desc_str.eq."M06-L".or. &
                       desc_str.eq."m11-l".or.desc_str.eq."M11-L".or. &
                       desc_str.eq."tpss".or.desc_str.eq."TPSS".or. &
                       desc_str.eq."revtpss".or.desc_str.eq."REVTPSS".or. &
                       desc_str.eq."tpssloc".or.desc_str.eq."TPSSLOC".or. &
                       desc_str.eq."scan".or.desc_str.eq."SCAN".or. &
                       libxc_family.eq.XC_FAMILY_MGGA) then
                  use_density_gradient = .true.
                  use_gga = .true.
                  use_basis_gradients = .true.
               ! Include Meta-GGA contributions i.e. kinetic energy density
                  use_meta_gga = .true.
               elseif (desc_str.eq."m06" .or. desc_str.eq."M06" .or. &
                       desc_str.eq."m06-2x" .or. desc_str.eq."M06-2X" .or. &
                       desc_str.eq."m06-hf" .or. desc_str.eq."M06-HF" .or. &
                       desc_str.eq."m08-hx" .or. desc_str.eq."M08-HX" .or. &
                       desc_str.eq."m08-so" .or. desc_str.eq."M08-SO" .or. &
                       desc_str.eq."m11" .or. desc_str.eq."M11" .or. &
                       libxc_family.eq.XC_FAMILY_HYB_MGGA) then
               ! Need this to stop hybrid_coeff being altered
                  use_minnesota = .true.
                  use_density_gradient = .true.
                  use_gga = .true.
                  use_basis_gradients = .true.
               ! Include Meta-GGA contributions i.e. kinetic energy density
                  use_meta_gga = .true.
               ! require exact exchange
                  use_hartree_fock = .true.
                  if (desc_str.eq."m11" .or. desc_str.eq."M11") then
                     use_hse=.true.
                     use_lc_wpbeh = .true.
                  else
                     flag_default_lvl = .true.
                  endif
               ! I've consolidated all the hybrids in to this one command so
               ! we can also do the hybrids from LibXC (hopefully) AJL/August 2016
               elseif (desc_str.eq."hse03" .or. desc_str.eq."HSE03" .or. &
                       desc_str.eq."hse06" .or. desc_str.eq."HSE06" .or. &
                       desc_str.eq."lc_wpbeh".or.desc_str.eq."LC_wPBEh" .or. &
                       desc_str.eq."pbe0" .or. desc_str.eq.'PBE0' .or. &
                       desc_str.eq."pbesol0" .or. desc_str.eq.'PBESOL0' .or. &
                       desc_str.eq."b3lyp" .or. desc_str.eq."B3LYP" .or. &
                       desc_str.eq."b1lyp" .or. &
                       libxc_family.eq.XC_FAMILY_HYB_GGA) then

                  ! require the density gradient
                  use_density_gradient = .true.
                  ! require gradient terms of the XC potential
                  use_gga = .true.
                  use_basis_gradients = .true.
                  ! require exact exchange
                  use_hartree_fock = .true.
                  ! Now for more specific definitions
                  if (desc_str.eq."lc_wpbeh".or.desc_str.eq."LC_wPBEh") then
                    use_lc_wpbeh = .true.
                  else
                    flag_default_lvl = .true.
                  endif

                  if (desc_str.eq."lc_wpbeh".or.desc_str.eq."LC_wPBEh".or. &
                      desc_str.eq."hse03" .or. desc_str.eq."HSE03" .or. &
                      desc_str.eq."hse06" .or. desc_str.eq."HSE06") then
                    use_hse = .true.
                  elseif (desc_str.eq."pbe0" .or. desc_str.eq.'PBE0') then
                    use_pbe0 = .true.
                  elseif (desc_str.eq."pbesol0" .or. desc_str.eq.'PBESOL0') then
                    use_pbesol0 = .true.
                  elseif (desc_str.eq."b3lyp" .or. desc_str.eq."B3LYP") then
                    use_b3lyp = .true.
                  elseif (desc_str.eq."b1lyp") then
                    use_b1lyp = .true.
                  endif

!               elseif (desc_str.eq."hse03" .or. &
!                       desc_str.eq."HSE03" .or. &
!                       desc_str.eq."hse06" .or. &
!                       desc_str.eq."HSE06" .or. ) then
!                  use_hse = .true.
!               ! require the density gradient
!                  use_density_gradient = .true.
!               ! require gradient terms of the XC potential
!                  use_gga = .true.
!                  use_basis_gradients = .true.
!               ! require exact exchange
!                  use_hartree_fock = .true.
!                  flag_default_lvl = .true.
!               elseif (desc_str.eq."lc_wpbeh".or.&
!                        desc_str.eq."LC_wPBEh") then
!                  use_hse = .true.
!                  use_lc_wpbeh = .true.
!                 ! require the density gradient
!                  use_density_gradient = .true.
!                 ! require gradient terms of the XC potential
!                  use_gga = .true.
!                  use_basis_gradients = .true.
!                 ! require exact exchange
!                  use_hartree_fock = .true.
!               elseif (desc_str.eq."pbe0" .or. desc_str.eq.'PBE0') then
!                  use_pbe0 = .true.
!                  ! require the density gradient
!                  use_density_gradient = .true.
!                  ! require gradient terms of the XC potential
!                  use_gga = .true.
!                  use_basis_gradients = .true.
!                  ! require exact exchange
!                  use_hartree_fock = .true.
!                  flag_default_lvl = .true.
!               elseif (desc_str.eq."b3lyp") then
!                  use_b3lyp = .true.
!                  ! require the density gradient
!                  use_density_gradient = .true.
!                  ! require gradient terms of the XC potential
!                  use_gga = .true.
!                  use_basis_gradients = .true.
!                  ! require exact exchange
!                  use_hartree_fock = .true.
!                  flag_default_lvl = .true.
               elseif (desc_str.eq."lrc-xyg3" .or. desc_str.eq."lrc-XYG3" .or. &
                       desc_str.eq."lrc-xyg6" .or. desc_str.eq."lrc-XYG6") then
                  use_b3lyp = .true.
                  ! require the density gradient
                  use_density_gradient = .true.
                  ! require gradient terms of the XC potential
                  use_gga = .true.
                  use_basis_gradients = .true.
                  ! require exact exchange
                  use_hartree_fock = .true.
                  flag_default_lvl = .true.
                  ! require post-SCF correlation
                  use_mp2    = .true.
                  use_dftpt2 = .true.
                  use_dftpt2_and_lrc = .true.
                  use_full_spectrum = .true.
               elseif (desc_str.eq."xyg3" .or. desc_str.eq."XYG3") then
                  use_b3lyp = .true.
                  ! require the density gradient
                  use_density_gradient = .true.
                  ! require gradient terms of the XC potential
                  use_gga = .true.
                  use_basis_gradients = .true.
                  ! require exact exchange
                  use_hartree_fock = .true.
                  flag_default_lvl = .true.
                  ! require post-SCF correlation
                  use_mp2    = .true.
                  use_dftpt2 = .true.
                  use_full_spectrum = .true.
               elseif (desc_str.eq."xygjos" .or. desc_str.eq."XYGJOS") then
                  use_pbe0 = .true.
                  ! require the density gradient
                  use_density_gradient = .true.
                  ! require gradient terms of the XC potential
                  use_gga = .true.
                  use_basis_gradients = .true.
                  ! require exact exchange
                  use_hartree_fock = .true.
                  flag_default_lvl = .true.
                  ! require post-SCF correlation
                  use_os_mp2    = .true.
                  use_dftpt2 = .true.
                  use_full_spectrum = .true.
               elseif (desc_str.eq."qpe-xygjos" .or. desc_str.eq."QPE-XYGJOS") then
                  use_pbe0 = .true.
                  ! require the density gradient
                  use_density_gradient = .true.
                  ! require gradient terms of the XC potential
                  use_gga = .true.
                  use_basis_gradients = .true.
                  ! require exact exchange
                  use_hartree_fock = .true.
                  flag_default_lvl = .true.
                  ! require post-SCF correlation
                  use_os_mp2_qpe    = .true.
                  use_dftpt2_qpe = .true.
                  use_full_spectrum = .true.
               elseif (desc_str.eq."xygjros" .or. desc_str.eq."XYGJROS") then
                  use_pbe0 = .true.
                  ! require the density gradient
                  use_density_gradient = .true.
                  ! require gradient terms of the XC potential
                  use_gga = .true.
                  use_basis_gradients = .true.
                  ! require exact exchange
                  use_hartree_fock = .true.
                  flag_default_lvl = .true.
                  ! require post-SCF correlation
                  use_os_rpa_qpe    = .true.
                  use_dftpt2_qpe = .true.
                  use_full_spectrum = .true.
               elseif (desc_str.eq."xdh-pbe0" .or. desc_str.eq."xDH-PBE0") then
                  use_pbe0 = .true.
                  ! require the density gradient
                  use_density_gradient = .true.
                  ! require gradient terms of the XC potential
                  use_gga = .true.
                  use_basis_gradients = .true.
                  ! require exact exchange
                  use_hartree_fock = .true.
                  flag_default_lvl = .true.
                  ! require post-SCF correlation
                  use_os_mp2    = .true.
                  use_dftpt2 = .true.
                  use_full_spectrum = .true.
               elseif (desc_str.eq."zrps" .or. desc_str.eq."ZRPS") then
                  use_pbe0 = .true.
                  ! require the density gradient
                  use_density_gradient = .true.
                  ! require gradient terms of the XC potential
                  use_gga = .true.
                  use_basis_gradients = .true.
                  ! require exact exchange
                  use_hartree_fock = .true.
                  flag_default_lvl = .true.
                  ! require post-SCF correlation
                  use_os_mp2    = .true.
                  use_dftpt2 = .true.
                  use_full_spectrum = .true.
!               elseif (desc_str.eq."pbesol0") then
!                  use_pbesol0 = .true.
!                  ! require the density gradient
!                  use_density_gradient = .true.
!                  ! require gradient terms of the XC potential
!                  use_gga = .true.
!                  use_basis_gradients = .true.
!                  ! require exact exchange
!                  use_hartree_fock = .true.
!                  flag_default_lvl = .true.
               elseif (desc_str.eq."hf".or.desc_str.eq."HF") then
                  use_hartree_fock = .true.
                  flag_default_lvl = .true.
               elseif (desc_str.eq."hf_kspace".or.desc_str.eq."HF_kspace") then
                  use_hf_kspace = .true.
                  use_hartree_fock = .true.
               elseif (desc_str.eq.'mp2'.or.desc_str.eq.'MP2') then
                  use_hartree_fock = .true.
                  use_mp2 = .true.
                  use_full_spectrum = .true.
               elseif (desc_str.eq.'screx') then
                  use_hartree_fock = .true.
                  use_screx = .true.
               elseif (desc_str.eq.'cohsex') then
                  use_hartree_fock = .true.
                  use_screx = .true.
                  use_cohsex = .true.
                  use_full_spectrum = .true.
               endif
            elseif (desc_str.eq.'neutral_excitation') then
                  use_prodbas = .true.
!                  use_density_gradient = .true.
                  flag_neutral_excitation = .true.
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if(desc_str.eq.'bse') then
                  use_bse = .true.
               endif
            elseif(desc_str.eq.'tddft_kernel') then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if(desc_str.eq.'libxc') then
                  use_libxc_tddft = .true.
               endif
            elseif (desc_str.eq."qpe_calc") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if(desc_str.eq.'gw'.or.desc_str.eq.'GW') then
                  use_gw = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'gw_expt'.or.desc_str.eq.'GW_EXPT') then
                  use_gw_expt = .true.
                  use_gw = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'ev_scgw') then
                  use_gw = .true.
                  use_ev_scgw = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'ev_scgw0') then
                  use_gw = .true.
                  use_ev_scgw0 = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'gw+2ox'.or.desc_str.eq.'GW+2OX') then
                  use_gw2ox = .true.
                  use_gw = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'gw+soxw'.or.desc_str.eq.'GW+SOXW') then
                  use_gwsoxw = .true.
                  use_gw = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'gw+s2ox'.or.desc_str.eq.'GW+S2OX') then
                  use_gws2ox = .true.
                  use_sosex_selfe = .true.
                  use_gw = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'gw+sosex2w'.or.desc_str.eq.'GW+SOSEX2W') then
                  use_gws2ox = .true.
                  use_sosex2w_selfe = .true.
                  use_gw = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'mp2'.or.desc_str.eq.'MP2') then
                  use_mp2sf = .true.
                  use_full_spectrum = .true.
               endif
            elseif (desc_str.eq."contour_def_gw".or.desc_str.eq."CONTOUR_DEF_GW") then
                  use_contour_def_gw = .true.
            elseif (desc_str.eq."use_pimd_wrapper") then
                use_pimd_wrapper = .true.
                use_forces = .true.
            elseif (desc_str.eq."path_integral") then     ! XZL: added for PIMD
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if(desc_str.eq.'pimd'.or.desc_str.eq.'PIMD') then
                  use_pathint = .true.
                  use_pimd = .true.
               elseif(desc_str.eq.'cmd'.or.desc_str.eq.'CMD') then
                  use_pathint = .true.
                  use_cmd = .true.
               elseif(desc_str.eq.'rpmd'.or.desc_str.eq.'RPMD') then
                  use_pathint = .true.
                  use_rpmd = .true.
               endif
            elseif (desc_str.eq."normal_mode") then      ! XZL: added for PIMD
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if(desc_str.eq.'nm'.or.desc_str.eq.'NM') then
                 use_nm = .true.
               endif
            elseif (desc_str.eq."sc_self_energy") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if(desc_str.eq.'scgw'.or.desc_str.eq.'SCGW') then
                  use_scgw = .true.
                  use_rpa_ene = .false.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'scgw0'.or.desc_str.eq.'SCGW0') then
                  use_scgw0 = .true.
                  use_rpa_ene = .false.
                  use_full_spectrum = .true.
               elseif (desc_str.eq."dmft_gw".or.desc_str.eq."DMFT_GW"&
                       .or.desc_str.eq."DMFT_gw".or.desc_str.eq."dmft_GW") then
                  print *, " ----- > test"
                  use_dmft_gw = .true.
               elseif (desc_str.eq."dmft_pbe0".or.desc_str.eq."DMFT_PBE0"&
                       .or.desc_str.eq."DMFT_pbe0".or.desc_str.eq."dmft_PBE0") then
                  print *, " ----- > test"
                  use_dmft_pbe0 = .true.
               endif
            elseif (desc_str.eq."DFPT") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if(desc_str.eq.'vibration') then
                  use_DFPT=.true.
                  use_full_spectrum = .true.
                  use_partition_deriv = .false.
                  use_ASR = .true.
                  use_basis_gradients = .true.
               elseif(desc_str.eq.'vibration_reduce_memory') then
                  use_DFPT_reduce_memory=.true.
                  use_full_spectrum = .true.
                  use_partition_deriv = .false.
                  use_ASR = .true.
                  use_basis_gradients = .true.
               elseif(desc_str.eq.'vibration_with_moving_grid_effect') then
                  use_DFPT=.true.
                  use_full_spectrum = .true.
                  use_partition_deriv = .true.
                  use_ASR = .false.
                  use_basis_gradients = .true.
               elseif(desc_str.eq.'vibration_without_moving_grid_effect') then
                  use_DFPT=.true.
                  use_full_spectrum = .true.
                  use_partition_deriv = .false.
                  use_ASR = .false.
                  use_basis_gradients = .true.
               elseif(desc_str.eq.'polarizability') then
                  use_DFPT_polarizability=.true.
                  use_full_spectrum = .true.
                  use_partition_deriv = .false.
               elseif(desc_str.eq.'polar_reduce_memory') then
                  use_DFPT_polar_reduce_memory=.true.
                  use_full_spectrum = .true.
                  use_partition_deriv = .false.
               elseif(desc_str.eq.'dielectric') then
                  use_DFPT_dielectric=.true.
                  use_full_spectrum = .true.
                  use_partition_deriv = .false.
               elseif(desc_str.eq.'phonon_gamma') then
                  use_DFPT_phonon_gamma=.true.
                  use_full_spectrum = .true.
                  use_partition_deriv = .false.
                  use_basis_gradients = .true.
               elseif(desc_str.eq.'phonon') then
                  use_DFPT_phonon=.true.
                  use_full_spectrum = .true.
                  use_partition_deriv = .false.
                  use_basis_gradients = .true.
               elseif(desc_str.eq.'phonon_reduce_memory') then
                  use_DFPT_phonon_reduce_memory=.true.
                  use_full_spectrum = .true.
                  use_partition_deriv = .false.
                  use_basis_gradients = .true.
               endif
            elseif (desc_str.eq."calculate_friction") then
               use_friction = .true.
               use_full_spectrum = .true.
               use_partition_deriv = .false.
               use_basis_gradients = .true.
            elseif (desc_str.eq."friction_numeric_disp") then
            elseif (desc_str.eq."friction_iter_limit") then
            elseif (desc_str.eq."friction_delta_type") then
            elseif (desc_str.eq."friction_window_size") then
            elseif (desc_str.eq."friction_broadening_width") then
            elseif (desc_str.eq."friction_discretization_length") then
            elseif (desc_str.eq."friction_max_energy") then
            elseif (desc_str.eq."friction_methfessel_n") then
            elseif (desc_str.eq."friction_temperature") then
            elseif (desc_str.eq."friction_q_integration") then
            elseif (desc_str.eq."friction_output_spectrum") then
            elseif (desc_str.eq."friction_read_matrices") then
            elseif (desc_str.eq."friction_accuracy_etot") then
            elseif (desc_str.eq."friction_accuracy_rho") then
            elseif (desc_str.eq."friction_accuracy_eev") then
            elseif (desc_str.eq."friction_accuracy_potjump") then
            elseif (desc_str.eq."friction_coupling_matrix_mode") then

            elseif (desc_str.eq."total_energy_method") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if (desc_str == "dfauto") then
                   read (inputline, *, end=88, err=99) desc_str, desc_str, desc_str
                   select case (str_lower(desc_str))
                   case ('pbe')
                       libxc_family_post = XC_FAMILY_GGA
                   case ('tpss', 'scan')
                       libxc_family_post = XC_FAMILY_MGGA
                   case ('pbe0')
                       libxc_family_post = XC_FAMILY_HYB_GGA
                   case ('scan0')
                       libxc_family_post = XC_FAMILY_HYB_MGGA
                   end select
               end if
               if(desc_str.eq.'HF'.or.desc_str.eq.'hf') then
                  use_hf_post = .true.
               elseif(desc_str.eq.'RPA'.or.desc_str.eq.'rpa' .or.desc_str.eq.'klein') then
                  use_rpa_ene = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'RPA+2OX'.or.desc_str.eq.'rpa+2ox') then
                  use_rpa_ene = .true.
                  use_rpa_plus_2ox = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'RPA+SOSEX'.or.desc_str.eq.'rpa+sosex') then
                  use_rpa_ene = .true.
                  use_rpa_plus_2ox = .true.
                  use_rpa_plus_sosex = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'rPT2'.or.desc_str.eq.'rpt2') then
                  use_rpa_ene = .true.
                  use_rpa_plus_2ox = .true.
                  use_rpa_plus_sosex = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'gw'.or.desc_str.eq.'GW') then
                  use_gw = .true.
                  use_gw_energy = .true.
               elseif(desc_str.eq.'mp2'.or.desc_str.eq.'MP2') then
                  use_mp2 = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'mp2_blacs'.or.desc_str.eq.'MP2_BLACS') then
                  use_mp2 = .true.
                  use_mp2_blacs = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'os-pt2'.or.desc_str.eq.'OS-PT2') then
                  use_mp2 = .true.
                  use_os_mp2 = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'ci'.or.desc_str.eq.'CI') then
                  use_ci = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'cc'.or.desc_str.eq.'CC') then
                  use_cc = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'lrc-xyg3'.or.desc_str.eq.'lrc-XYG3' .or. &
                      desc_str.eq.'lrc-xyg6'.or.desc_str.eq.'lrc-XYG6' ) then
                  use_mp2    = .true.
                  use_dftpt2 = .true.
               elseif(desc_str.eq.'xyg3'.or.desc_str.eq.'XYG3') then
                  use_mp2    = .true.
                  use_dftpt2 = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'osRPA'.or.desc_str.eq.'osrpa') then
                  use_os_rpa_qpe = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'osPT2'.or.desc_str.eq.'ospt2') then
                  use_os_mp2_qpe = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'SIC-RPA'.or.desc_str.eq.'sic-rpa') then
                  use_sic_rpa_qpe = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'xdh-pbe0'.or.desc_str.eq.'xDH-PBE0') then
                  use_os_mp2 = .true.
                  use_dftpt2 = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'zrps'.or.desc_str.eq.'ZRPS') then
                  use_os_mp2 = .true.
                  use_dftpt2 = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'C6_coef') then
                  use_C6_coef = .true.
                  use_full_spectrum = .true.
               elseif(desc_str.eq.'ll_vdwdf'.or.desc_str.eq.'LL_vdwdf') then
                  use_ll_vdwdf = .true.
               elseif(desc_str.eq.'M06-L'.or.desc_str.eq.'m06-l' .or. &
                      desc_str.eq.'M11-L'.or.desc_str.eq.'m11-l' .or. &
                      desc_str.eq.'TPSS'.or.desc_str.eq.'tpss' .or. &
                      desc_str.eq.'revTPSS'.or.desc_str.eq.'revtpss' .or. &
                      desc_str == 'scan!' .or. desc_str == 'SCAN!' .or. &
                      desc_str.eq.'TPSSloc'.or.desc_str.eq.'tpssloc' &
                      .or. libxc_family_post == XC_FAMILY_GGA &
                      .or. libxc_family_post == XC_FAMILY_MGGA &
               ) then
                  use_meta_gga_post = .true.
                  use_basis_gradients = .true.
                  use_density_gradient = .true.
               elseif(desc_str.eq.'M06'.or.desc_str.eq.'m06' .or. &
                      desc_str.eq.'M06-2X'.or.desc_str.eq.'m06-2x' .or. &
                      desc_str.eq.'M06-HF'.or.desc_str.eq.'m06-hf' .or. &
                      desc_str.eq.'M08-HX'.or.desc_str.eq.'m08-hx' .or. &
                      desc_str.eq.'M08-SO'.or.desc_str.eq.'m08-so' &
                      .or. libxc_family_post == XC_FAMILY_HYB_GGA &
                      .or. libxc_family_post == XC_FAMILY_HYB_MGGA &
               ) then
                  ! Need to think of some new flags for M11? On the to-do list. AJL
                  use_meta_gga_post = .true.
                  use_hybrid_meta = .true.
                  use_basis_gradients = .true.
                  use_density_gradient = .true.
               elseif(desc_str.eq.'nlcorr'.or.desc_str.eq.'NLCORR') then  !SAG
                  use_nlcorr_post = .true.
                  use_nlcorr_in_xc = .false.
               elseif( desc_str.eq.'PBE_VDW'.or.desc_str.eq.'pbe_vdw'.or.&
                    desc_str.eq.'REVPBE_VDW'.or.desc_str.eq.'revpbe_vdw') then !SAG
                  use_vdw_post = .true.
             endif
            elseif (desc_str.eq."use_hf_kspace") then
               read(inputline,*,end=88,err=99) desc_str, use_hf_kspace
            elseif (desc_str.eq."use_hf_kspace_with_rpa") then
               read(inputline,*,end=88,err=99) desc_str, use_hf_kspace_with_rpa
            elseif (desc_str.eq."use_angular_division") then
               read(inputline,*,end=88,err=99) desc_str, use_angular_division
            elseif (desc_str.eq."esp_constraint") then
               read(inputline,*,end=88,err=99) desc_str, esp_constraint
            elseif (desc_str.eq. "output_vacuum_length") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
           !    n_cube = n_cube+1
           !    use_cube_output = .true.
            elseif (desc_str.eq."output") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if (desc_str.eq."cube") then
                  n_cube = n_cube+1
                  use_cube_output = .true.
               elseif (desc_str.eq. "realspace_esp") then
                  read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str,desc_str,desc_str
                  n_cube = n_cube+1
                  use_cube_output = .true.
               elseif (desc_str.eq. "on-site_esp") then
               elseif (desc_str.eq."esp") then
                  n_esp = n_esp+1
                 use_esp = .true.
               elseif (desc_str.eq."bare_coulomb_integrals") then
                  use_bare_ci = .true.
                  use_prodbas = .true.
               elseif (desc_str.eq."band") then
                  n_plot_band = n_plot_band +1
               elseif (desc_str.eq."band_mulliken") then
                  n_plot_band = n_plot_band +1
               !TZ
               elseif (desc_str .eq. "dielectric") then
                  n_plot_dielectric = n_plot_dielectric + 1
               elseif (desc_str .eq. "absorption") then
                  n_plot_absorption = n_plot_absorption + 1
               !TZ
               elseif (desc_str.eq."soc_eigenvectors") then
                  n_write_soc_eigenvectors = n_write_soc_eigenvectors + 1
               elseif (desc_str.eq."band_during_scf") then
                  n_plot_band_scf = n_plot_band_scf +1
               elseif (desc_str.eq."DFPT_phonon_band") then
                  n_plot_band_DFPT_phonon = n_plot_band_DFPT_phonon + 1
               elseif (desc_str.eq."rrs_pbc_band") then
                  rrs_pbc_n_plot_band = rrs_pbc_n_plot_band +1
               elseif (desc_str.eq."friction_matrices") then
               elseif (desc_str.eq."friction_eigenvectors") then
               elseif (desc_str.eq."ks_coulomb_integral") then
                  use_coulomb_integral_ks = .true.
               endif
            elseif (desc_str.eq."constraint_electrons") then
               use_constraint = .true.
               n_region=n_region+1
            elseif (desc_str.eq."multiplicity") then
               use_hf_multiplicity = .true.
            elseif (desc_str.eq."relax_geometry") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if (desc_str.ne.'none') then
               ! require the gradient of the charge density later on
                  use_geo_relaxation = .true.
                  use_forces = .true.
                  use_basis_gradients = .true.
               endif
            elseif (desc_str.eq."sc_accuracy_forces") then
               ! if sc_accuracy_forces was defined, assume reasonably that the
               ! user wishes to compute forces and just forgot the compute_forces flag ...
               use_basis_gradients = .true.
               use_forces = .true.
            else if (desc_str .eq. 'relax_unit_cell') then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if (desc_str.eq.'none') then
                   relax_unit_cell = 0
               ! Here come other types of unit cell relaxations
               else if (desc_str.eq.'full') then
                   relax_unit_cell = 1
               else if ( (desc_str.eq.'fixed_angles') .or. (desc_str.eq.'shape') ) then
                   relax_unit_cell = 2
               end if
            elseif (desc_str.eq."species") then
               i_species = i_species+1
               call parse_species ( i_species, pseudoized, verbatim_writeout )
!              check whether species is actually a pseudoized species
               if(pseudoized) i_pp_species = i_pp_species + 1
               pseudoized = .false.
               coming_from_parse_species = .true. ! only true to prevent double writing of a line from control.in
            elseif (desc_str.eq."switch_external_pert") then
               use_embedding_potential = .true.
               use_atom_ref=.true.
            elseif (desc_str.eq."vdw_correction") then
               ! vdW correction
               use_vdw_correction = .true.
               call parse_vdw ()
            elseif (desc_str.eq."vdw_correction_hirshfeld") then
               read(inputline,*, iostat=i_code) desc_str, desc_str
               use_vdw_correction_hirshfeld = .true.
               if (i_code == 0 .and. desc_str /= '') then
                  read(inputline,*,end=88,err=99) &
                       & desc_str, use_vdw_correction_hirshfeld
               end if
            elseif (desc_str.eq."vdw_correction_hirshfeld_sc") then
               read(inputline,*, iostat=i_code) desc_str, desc_str
               use_vdw_correction_hirshfeld_sc = .true.
               if (i_code == 0 .and. desc_str /= '') then
                  read(inputline,*,end=88,err=99) &
                       & desc_str, use_vdw_correction_hirshfeld_sc
               end if
            elseif (desc_str.eq.'vdw_pair_ignore') then
               n_vdw_pairs_ignore = n_vdw_pairs_ignore + 1
            elseif (desc_str.eq."qmmm") then
               use_qmmm = .true.
            else if (desc_str.eq."TDDFT_run") then
            else if (desc_str.eq."TDDFT_propagator") then
            else if (desc_str.eq."MD_run") then
               use_molecular_dynamics = .true.
               use_forces             = .true.
               use_basis_gradients = .true.
            else if (desc_str.eq."MD_QH_init") then
               use_MD_QH_init = .true.
               MD_QH_init_segments = MD_QH_init_segments+1
            else if (desc_str == "many_body_dispersion_alt") then
                use_mbd_old = .true.
            else if (desc_str == "many_body_dispersion_dev") then
                use_mbd_dev = .true.
            else if ( &
                desc_str == "many_body_dispersion" &
                .or. desc_str == "many_body_dispersion_nl" &
            ) then
                use_libmbd = .true.
            else if (desc_str == "many_body_dispersion_pre2019") then
                use_mbd_std = .true.
            else if (desc_str == "many_body_dispersion_rsscs") then
                write(info_str,'(A)') &
                     "* Please note: The keyword many_body_dispersion_rsscs has changed"
                call localorb_info(info_str)
                write(info_str,'(A)') &
                     "* to many_body_dispersion"
                call localorb_info(info_str)
                write(info_str,'(A)') &
                     "* Stopping the code for now."
                call localorb_info(info_str)
                call aims_stop_coll('', '') 
            else if (desc_str.eq."MD_schedule") then
               MD_segments     = 0
               TDI_segments    = 0
               MD_segment_loop: do
                  read(777,'(A)',iostat=i_code) inputline
                  if(i_code/=0) then
                     ! handle error at main level
                     backspace(777)
                     exit MD_segment_loop
                  endif
                  read(inputline,*,iostat=i_code) desc_str
                  if (i_code/=0)cycle MD_segment_loop ! skip empty line
                  if (desc_str(1:1)=='#')cycle MD_segment_loop ! skip comment
                  if (desc_str.eq.'MD_segment') then
                     MD_segments = MD_segments + 1
                     ! CC: Ouput line to control.in
                     call localorb_info(inputline,use_unit,'(2X,A)') ! Write the control file to the output line by line
                  elseif ( (desc_str.eq."thermodynamic_integration") .or. (desc_str.eq."harmonic_potential_only") ) then
                     TDI_segments = TDI_segments + 1
                     ! CC: Ouput line to control.in
                     call localorb_info(inputline,use_unit,'(2X,A)') ! Write the control file to the output line by line
                     if ( desc_str.eq."harmonic_potential_only") then
                       use_harmonic_pot_only = .true.
                     end if
                  elseif (desc_str.eq."adiabatic_switching") then
                     AS_segments = AS_segments + 1
                     ! CC: Ouput line to control.in
                     call localorb_info(inputline,use_unit,'(2X,A)') ! Write the control file to the output line by line
                  else
                     backspace(777)
                     exit MD_segment_loop
                  endif
               enddo MD_segment_loop
               if (MD_segments.gt.0) then
                  MD_use_schedule = .true.
                  use_molecular_dynamics = .true.
                  use_forces             = .true.
                  use_basis_gradients    = .true.
               endif
               if (TDI_segments.gt.0) then
                  use_thermodynamic_integration = .true.
               else
                  use_thermodynamic_integration = .false.
               endif
               if (AS_segments.gt.0) then
                     use_reffree_AS = .true.
               endif
            elseif (desc_str.eq."force_occupation_projector") then
               n_force_occ = n_force_occ + 1
               force_occupation_projector = .true.
            elseif (desc_str .eq. 'force_occupation_basis') then
               n_force_occ = n_force_occ + 1
               force_occupation_basis = .true.
            elseif (desc_str .eq. 'force_occupation_projector_redsub') then
               force_occupation_projector_redsub = .true.
            elseif (desc_str .eq. 'force_occupation_projector_auto') then
               force_occupation_projector_auto = .true.
            elseif (desc_str .eq. 'apply_boys') then
               n_sub_boys = n_sub_boys + 1
               apply_boys_flag = .true.
            elseif (desc_str.eq."compute_dielectric") then
               read(inputline,*,end=88,err=99) desc_str, desc_str,&
                                               desc_str
               flag_out_dielectric = .true.
               use_basis_gradients = .true.
               !use_scalapack = .true.
            elseif (desc_str .eq. "dielectric_broadening") then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str
               flag_broadening = .true.
            elseif (desc_str.eq."compute_kubo_greenwood") then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, &
                  desc_str, desc_str, desc_str, desc_str, desc_str, desc_str, desc_str
               flag_out_greenwood = .true.
            elseif (desc_str.eq."greenwood_method") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               flag_greenwood_method = .true.
            elseif (desc_str.eq."flexible_fermi_level") then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str,&
                                               desc_str
               flag_flex_fermi = .true.
            elseif (desc_str.eq."explicit_fermi") then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str
               flag_explicit_fermi = .true.
            elseif (desc_str.eq."compute_momentummatrix") then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str,&
                                               desc_str
               use_basis_gradients = .true.
            elseif (desc_str.eq."compute_dipolematrix") then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str,&
                                                                   desc_str
               flag_out_dipmat = .true.
            elseif (desc_str.eq."compute_dipolematrix_k_k") then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str
               flag_out_dipmat_k_k = .true.
            elseif (desc_str.eq."compute_absorption") then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, &
                      desc_str, desc_str, desc_str, desc_str, desc_str, desc_str
               flag_out_absorption = .true.
               use_basis_gradients = .true.
            elseif (desc_str.eq."compute_coulomb_matrix_ovl") then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str
                  flag_out_coulmat_ovl = .true.
            elseif (desc_str.eq."compute_coulomb_matrix_lvl") then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str
                  flag_out_coulmat_lvl = .true.
                  use_lvl = .true.
                  use_prodbas=.true.
!                  use_hartree_fock=.true.
                  use_periodic_hf=.true.
                  use_hf_kspace=.true.
                  use_cutCb=.true.
                  use_logsbt=.true.
            elseif (desc_str.eq."compute_greenwood_dc_transport") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               flag_out_dclimit = .true.
            elseif (desc_str .eq. 'el_ph_coupling') then
               use_el_ph_coupling = .true.
            elseif ( desc_str.eq."use_logsbt" ) then
               read(inputline,*,end=88,err=99) desc_str, use_logsbt
            else if (desc_str.eq."solvent") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if (desc_str.eq."mpb") then
                  call parse_mpb(verbatim_writeout)
               end if
            elseif ( desc_str.eq."RI_method" ) then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if ( (desc_str.eq.'LVL_full') .or. (desc_str.eq.'lvl_full') &
                    .or. (desc_str.eq.'LVL_2nd') .or. (desc_str.eq.'lvl_2nd') ) then
                 use_lvl = .true.
               else if ( (desc_str.eq.'LVL') .or. (desc_str.eq.'lvl') &
                    .or. (desc_str.eq.'LVL_fast') .or. (desc_str.eq.'lvl_fast') ) then
                 ! This case must be handled analogous to read_control.f90.
                 ! Check there for the current defaults associated with each keyword.
                 use_lvl = .true.
                 use_lvl_fast = .true.
               else if ( (desc_str.eq.'V') .or. (desc_str.eq.'v') &
                         .or. (desc_str.eq.'SVS') .or. (desc_str.eq.'svs') ) then
                 use_lvl = .false.
               else
                 call aims_stop_coll("* dimensions.f90: Unknown RI_method encountered in control.in", func)
               end if
               flag_RI = .true.
            ! RRS-PBC scheme, igor
            elseif (desc_str .eq. 'rrs_pbc') then
               use_rrs_pbc = .true.
            elseif (desc_str .eq. 'rrs_pbc_unit_cell') then
               read(inputline,*,end=88,err=99) desc_str, &
                   rrs_pbc_n_center_atom, rrs_pbc_n_equal
               call read_rrs_pbc_unit_cell(777,1)
            elseif (desc_str .eq. 'rrs_pbc_k_grid') then
               continue
            elseif (desc_str .eq. "memtrace_use") then
               call aims_stop_coll("* memtrace is no longer supported by FHI-aims.  Please remove &
                                   &memtrace-related flags from your control.in file.", func)
            elseif (desc_str .eq. "memtrace_summary") then
               call aims_stop_coll("* memtrace is no longer supported by FHI-aims.  Please remove &
                                   &memtrace-related flags from your control.in file.", func)
            elseif (desc_str .eq. "memtrace_targetdir") then
               call aims_stop_coll("* memtrace is no longer supported by FHI-aims.  Please remove &
                                   &memtrace-related flags from your control.in file.", func)
            elseif (desc_str .eq. 'KS_eigenfunc_conjg') then
               read(inputline,*,end=88,err=99) desc_str, flag_KS_eigenfunc_conjg
            elseif (desc_str .eq. 'force_single_restartfile') then
               read(inputline,*,end=88,err=99) desc_str, flag_force_single_restartfile
            elseif (desc_str .eq. 'restart_eigenvectors_periodic') then
               read(inputline,*,end=88,err=99) desc_str, flag_restart_eigenvectors_periodic
            else if (desc_str == 'magnetic_response') then
               use_basis_gradients = .true.
            ! keyword for evaluating xc values for various DFTs
            else if (desc_str.eq."printout_dft_components") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if (desc_str.eq.'SCAN'.or.desc_str.eq.'scan' .or. &
                   desc_str.eq.'TPSS'.or.desc_str.eq.'tpss' .or. &
                   desc_str.eq.'revTPSS'.or.desc_str.eq.'revtpss' .or. &
                   desc_str.eq.'TPSSloc'.or.desc_str.eq.'tpssloc') then
                  use_meta_gga_printout = .true.
                  use_basis_gradients = .true.
                  use_density_gradient = .true.
               end if
            elseif (desc_str.eq.'split_atoms') then
               read(inputline,*,end=88,err=99) desc_str, flag_split_atoms
            elseif (desc_str.eq.'split_min_val') then
               flag_split_min_val = .true.
            elseif (desc_str.eq.'split_max_val') then
               flag_split_max_val = .true.
            elseif (desc_str.eq.'split_batch') then
               flag_split_batch = .true.
            elseif ( (desc_str.eq."empty_states") .or. &
                     (desc_str.eq."use_full_spectrum") .or. &
                     (desc_str.eq."charge") .or. &
                     (desc_str.eq."KS_method") .or. &
                     (desc_str.eq."elpa_settings") .or. &
                     (desc_str.eq."ELPA_AEO") .or. &
                     (desc_str.eq."default_initial_moment") .or. &
                     (desc_str.eq."basis_threshold") .or. &
                     (desc_str.eq."prodbas_nb") .or. &
                     (desc_str.eq."prodbas_threshold") .or. &
                     (desc_str.eq."use_smooth_prodbas_threshold") .or. &
                     (desc_str.eq."default_prodbas_acc") .or. &
                     (desc_str.eq."default_max_l_prodbas") .or. &
                     (desc_str.eq."default_max_n_prodbas") .or. &
                     (desc_str.eq."use_asym_inv_sqrt") .or. &
                     (desc_str.eq."hybrid_xc_coeff") .or. &
                     (desc_str.eq."partition_acc") .or. &
                     (desc_str.eq."wave_threshold") .or. &
                     (desc_str.eq."packed_matrix_threshold") .or. &
                     (desc_str.eq."legacy_monopole_extrapolation") .or. &
                     (desc_str.eq."mixer") .or. &
                     (desc_str.eq."mixer_swap_boundary") .or. &
                     (desc_str.eq."nbeads") .or. &             ! XZL: added for PIMD
                     (desc_str.eq."normal_mode") .or. &        ! XZL: added for PIMD
                     (desc_str.eq."PIMD_run") .or. &           ! XZL: added for PIMD
                     (desc_str.eq."PIMD_time_step") .or. &     ! XZL: added for PIMD
                     (desc_str.eq."communicate_pimd_wrapper") .or. & !MR: added for comm with ipi
                     (desc_str.eq."ini_linear_mixing") .or. &
                     (desc_str.eq."ini_linear_mix_param") .or. &
                     (desc_str.eq."ini_spin_mix_param") .or. &
                     (desc_str.eq."charge_mix_param") .or. &
                     (desc_str.eq."prec_mix_param") .or. &
                     (desc_str.eq."spin_mix_param") .or. &
                     (desc_str.eq."adjust_scf") .or. &
                     (desc_str.eq."frozen_core") .or. &
                     (desc_str.eq."occupation_type") .or. &
                     (desc_str.eq."set_blacsdim") .or. &
                     (desc_str.eq."set_blacsdim_cpt2") .or. &
                     (desc_str.eq."set_pair_block_cpt2") .or. &
                     (desc_str.eq."use_threadsafe_gwinit") .or. &
                     (desc_str.eq."occupation_acc") .or. &
                     (desc_str.eq."occupation_thr") .or. &
                     (desc_str.eq."mu_determination_method") .or. &
                     (desc_str.eq."fermi_acc") .or. &
                     (desc_str.eq."max_zeroin") .or. &
                     (desc_str.eq."sc_accuracy_rho") .or. &
                     (desc_str.eq."sc_accuracy_eev") .or. &
                     (desc_str.eq."sc_accuracy_etot") .or. &
                     (desc_str.eq."sc_accuracy_potjump") .or. &
                     (desc_str.eq."sc_iter_limit") .or. &
                     (desc_str.eq."sc_init_iter") .or. &
                     (desc_str.eq."constraint_it_lim") .or. &
                     (desc_str.eq."k_grid") .or. &
                     (desc_str.eq."k_offset") .or. &
                     (desc_str.eq."symmetry_reduced_k_grid") .or. &
                     (desc_str.eq."symmetry_reduced_k_grid_spg") .or. &
                     (desc_str.eq."symmetry_spg") .or. &
                     (desc_str.eq."use_spg_full_Delta") .or. &
                     (desc_str.eq."reconstruct_proper_only") .or. &
                     (desc_str.eq."use_spg_mv_mm") .or. &
                     (desc_str.eq."get_full_density_first") .or. &
                     (desc_str.eq."use_k_phase_exx") .or. &
                     (desc_str.eq."symmetry_thresh") .or. &
                     (desc_str.eq."use_full_symmetry") .or. &
                     (desc_str.eq."k_points_external") .or. &
                     (desc_str.eq."frozen_core_scf") .or. &
                     (desc_str.eq."frozen_core_scf_cutoff") .or. &
                     (desc_str.eq."frozen_core_scf_factor") .or. &
                     (desc_str.eq."constraint_precision") .or. &
                     (desc_str.eq."constraint_mix") .or. &
                     (desc_str.eq."constraint_debug") .or. &
                     (desc_str.eq."force_potential") .or. &
                     (desc_str.eq."force_lebedev") .or. &
                     (desc_str.eq."force_new_functional") .or. &
                     (desc_str.eq."force_occupation_smearing") .or. &
                     (desc_str.eq."hartree_fourier_part_th") .or. &
                     (desc_str.eq."hartree_partition_type") .or. &
                     (desc_str.eq."partition_type") .or. &
                     (desc_str.eq."stratmann_a") .or. &
                     (desc_str.eq."hartree_worksize") .or. &
                     (desc_str.eq."force_smooth_cutoff") .or. &
                     (desc_str.eq."sorted_eigenvalue_list") .or. &
                     (desc_str.eq."cube") .or. &
                     (desc_str.eq."fermisurface") .or. &
                     (desc_str.eq."cube_content_unit") .or. &
                     (desc_str.eq."output_cube_nth_iteration") .or. &
                     (desc_str.eq."overwrite_existing_cube_files") .or. &
                     (desc_str.eq."cube_default_size_safeguard") .or. &
                     (desc_str.eq."hydro_cut") .or. &
                     (desc_str.eq."packed_matrix_format") .or. &
                     (desc_str.eq."ini_linear_mixing_constraint") .or. &
                     (desc_str.eq."mixer_constraint") .or. &
                     (desc_str.eq."task_distribution") .or. &
                     (desc_str.eq."compute_kinetic") .or. &
                     (desc_str.eq."use_small_system_solver") .or. &
                     (desc_str.eq."use_density_matrix") .or. & ! VB: Legacy only
                     (desc_str.eq."use_density_matrix_hf") .or. &
                     (desc_str.eq."use_stress_tensor") .or. &
                     (desc_str.eq."numerical_stress_save_scf") .or. &
                     (desc_str.eq."collect_eigenvectors") .or. &
                     (desc_str.eq."hartree_convergence_parameter") .or. &
                     (desc_str.eq."ewald_radius") .or. &
                     (desc_str.eq."Ewald_radius") .or. &
                     (desc_str.eq."hartree_fp_function_splines") .or. &
                     (desc_str.eq."multip_radius_threshold") .or. &
                     (desc_str.eq."periodic_b_cut_threshold").or. &
                     (desc_str.eq."hartree_radius_threshold") .or. &
                     (desc_str.eq."adaptive_hartree_radius_th") .or. &
                     (desc_str.eq."multip_moments_rad_threshold") .or. &
                     (desc_str.eq."use_hartree_non_periodic_ewald") .or. &
                     (desc_str.eq."hartree_d_matrix_method") .or. &
                     (desc_str.eq."aggregated_energy_tolerance") .or. &
                     (desc_str.eq."energy_tolerance") .or. &
                     (desc_str.eq."harmonic_length_scale") .or. &
                     (desc_str.eq."max_atomic_move") .or. &
                     (desc_str.eq."min_line_step") .or. &
                     (desc_str.eq."line_step_reduce") .or. &
                     (desc_str.eq."max_relaxation_steps") .or. &
                     (desc_str.eq."init_hess") .or. &
                     (desc_str.eq."clean_forces").or. &
                     (desc_str.eq."multip_moments_threshold") .or. &
                     (desc_str.eq."multip_radius_free_threshold") .or. &
                     (desc_str.eq."multipole_threshold") .or. &
                     (desc_str.eq."distributed_spline_storage") .or. &
                     (desc_str.eq."communication_type") .or. &
                     (desc_str.eq."points_in_batch") .or. &
                     (desc_str.eq."grid_partitioning_method") .or. &
                     (desc_str.eq."walltime") .or. &
                     (desc_str.eq."KH_post_correction") .or. &
                     (desc_str.eq."preconditioner") .or. &
                     (desc_str.eq."precondition_before_mixer") .or. &
                     (desc_str.eq."precondition_max_l") .or. &
                     (desc_str.eq."density_update_method") .or. &
                     (desc_str.eq."density_restart").or. &
                     (desc_str.eq."keep_restart_density_frozen").or. &
                     (desc_str.eq."restart").or. &
                     (desc_str.eq."restart_rpa").or. &
                     (desc_str.eq."restart_pt2").or. &
                     (desc_str.eq."restart_read_only").or. &
                     (desc_str.eq."restart_write_only").or. &
                     (desc_str.eq."restart_save_iterations").or. &
                     (desc_str.eq."restart_relaxations").or. &
                     (desc_str.eq."restart_exx_write").or. &
                     (desc_str.eq."restart_exx_read").or. &
                     (desc_str.eq."restart_exx_3c_int").or. &
                     (desc_str.eq."write_restart_geometry").or. &
                     (desc_str.eq."hessian_to_restart_geometry").or. &
                     (desc_str.eq."distributed_hessian").or. &
                     (desc_str.eq."restart_eigenvectors_periodic").or. &
                     (desc_str.eq."orthonormalize_eigenvectors") .or. &
                     (desc_str.eq."auxil_basis").or. &
                     (desc_str.eq."final_forces_cleaned") .or. &
                     (desc_str.eq."grouping_factor") .or. &
                     (desc_str.eq."min_batch_size") .or. &
                     (desc_str.eq."batch_size_limit").or. &
                     (desc_str.eq."store_EV_to_disk_in_relaxation").or. &
                     (desc_str.eq."Adams_Moulton_integrator").or. &
                     (desc_str.eq."prune_basis_once") .or. &
                     (desc_str.eq."full_embedding") .or. &
                     (desc_str.eq."force_mpi_virtual_topo") .or. &
                     (desc_str.eq."batch_distribution_method") .or. &
                     (desc_str.eq."recompute_batches_in_relaxation") .or. &
                     (desc_str.eq."initial_ev_solutions") .or. &
                     (desc_str.eq."max_lopcg_iterations") .or. &
                     (desc_str.eq."lopcg_tolerance") .or. &
                     (desc_str.eq."lopcg_start_tolerance") .or. &
                     (desc_str.eq."lopcg_preconditioner") .or. &
                     (desc_str.eq."lopcg_omega") .or. &
                     (desc_str.eq."lopcg_skip_rate") .or. &
                     (desc_str.eq."lopcg_block_size") .or. &
                     (desc_str.eq."lopcg_use_all_evs") .or. &
                     (desc_str.eq."lopcg_slide_ev_block") .or. &
                     (desc_str.eq."lopcg_ritz_scalapack") .or. &
                     (desc_str.eq."lopcg_adaptive_tolerance") .or. &
                     (desc_str.eq."lopcg_auto_blocksize") .or. &
                     (desc_str.eq."squeeze_memory") .or. &
                     (desc_str.eq."use_local_index") .or. &
                     (desc_str.eq."load_balancing") .or. &
                     (desc_str.eq."heterostructure") .or. &
                     (desc_str.eq."use_mpi_in_place") .or. &
                     (desc_str.eq."use_alltoall") .or. &
                     (desc_str.eq."output_in_original_unit_cell") .or. &
                     (desc_str.eq."plus_u_petukhov_mixing") .or. &
                     (desc_str.eq."plus_u_out_eigenvalues") .or. &
                     (desc_str.eq."plus_u_matrix_control") .or. &
                     (desc_str.eq."plus_u_matrix_release") .or. &
                     (desc_str.eq."plus_u_use_hydros") .or. &
                     (desc_str.eq."plus_u_use_mulliken").or. &
                     (desc_str.eq."plus_u_use_full").or. &
                     (desc_str.eq."plus_u_ramping").or. &
                     (desc_str.eq."plus_u_ramping_accuracy").or. &
                     (desc_str.eq."hubbard_coefficient") .or. &
                     (desc_str.eq."plus_u_occ_mat_error") .or. &
                     (desc_str.eq."wf_extra_use_densmat") .or. &
                     (desc_str.eq."wf_extrapolation") .or. &
                     (desc_str.eq."wf_func") .or. &
                     (desc_str.eq."hf_version") .or. &
                     (desc_str.eq."hse_unit") .or. &
                     (desc_str.eq."sbtgrid_lnr0") .or. &
                     (desc_str.eq."sbtgrid_lnk0") .or. &
                     (desc_str.eq."sbtgrid_lnrange") .or. &
                     (desc_str.eq."sbtgrid_N") .or. &
                     (desc_str.eq."cutCb_width") .or. &
                     (desc_str.eq."cutCb_width_factor") .or. &
                     (desc_str.eq."cutCb_rcut") .or. &
                     (desc_str.eq."cutCb_rcut_factor") .or. &
                     (desc_str.eq."use_logsbt_for_radial_hse_integration").or.&
                     (desc_str.eq."logsbt_mp_far") .or. &
                     (desc_str.eq."use_2d_corr") .or. &
                     (desc_str.eq."line_search_tol").or. &
                     (desc_str.eq."bfgs_step") .or. &
                     (desc_str.eq."bfgs_extrapolation_cap") .or. &
                     (desc_str.eq."phonon") .or. &
                     (desc_str.eq."unfold") .or. &
                     (desc_str.eq."vibrations") .or. &
                     (desc_str.eq."use_ovlp_swap").or. &
                     (desc_str.eq."MD_time_step").or. &
                     (desc_str.eq."MD_MB_init").or. &
                     (desc_str.eq."MD_maxsteps").or. &
                     (desc_str.eq."check_MD_stop").or. &
                     (desc_str.eq."output_level").or. &
                     (desc_str.eq."MD_restart").or. &
                     (desc_str.eq."MD_restart_binary").or. &
                     (desc_str.eq."plumed").or. &
                     (desc_str.eq."plumed_new").or. &
                     (desc_str.eq."MD_clean_rotations").or. &
                     (desc_str.eq."MD_thermostat_units").or. &
                     (desc_str.eq."MD_RNG_seed").or. &
                     (desc_str.eq."scs_mp2_parameters").or. &
                     (desc_str.eq."n_max_coeff_3fn").or. &
                     (desc_str.eq."frequency_points").or. &
                     (desc_str.eq."plot_self_energy").or. &
                     (desc_str.eq."time_points").or. &
                     (desc_str.eq."maximum_frequency").or. &
                     (desc_str.eq."maximum_time").or. &
                     (desc_str.eq."freq_grid_type").or. &
                     (desc_str.eq."scgw_no_core").or. &
                     ! frozen_core and frozen virtual algorithm for post SCF methods Igor
                     (desc_str.eq."frozen_core_postSCF").or. &
                     (desc_str.eq."frozen_core_postscf").or. &
                     (desc_str.eq."frozen_virtual_postSCF").or. &
                     ! keywords related to the screening and coupling MP2 method
                     (desc_str.eq."en_shift").or. &
                     (desc_str.eq."en_shift_constant").or. &
                     (desc_str.eq."coupling_pt2_factor").or. &
                     (desc_str.eq."coupling_pt2_screen").or. &
                     (desc_str.eq."coupling_pt2_shift").or. &
                     ! keywords for osrpa and scsrpa
                     (desc_str.eq."scsrpa_coeff").or. &
                     (desc_str.eq."renormalized_eg").or. &
                     (desc_str.eq."n_lambda_osrpa").or. &
                     (desc_str.eq."osrpa_threshold").or. &
                     (desc_str.eq."osrpa_orders").or. &
                     ! keyword for the truncated CI method
                     (desc_str.eq."ci_truncation").or. &
                     (desc_str.eq."ci_single_excitations").or. &
                     (desc_str.eq."ci_single_spin_channel").or. &
                     (desc_str.eq."ci_accuracy_etot").or. &
                     (desc_str.eq."n_scf_ci").or. &
                     (desc_str.eq."ci_coeff_threshold").or. &
                     (desc_str.eq."ci_acc_method").or. &
                     (desc_str.eq."n_ci_acc_save").or. &
                     ! keyword fo coupled cluster calculation
                     (desc_str.eq."cc_method").or. &
                     (desc_str.eq."cc_scf_max_cycle").or. &
                     (desc_str.eq."cc_acc_rle").or. &
                     (desc_str.eq."cc_acc_diis").or. &
                     (desc_str.eq."cc_rle_vector_save").or. &
                     (desc_str.eq."cc_rle_thresh").or. &
                     (desc_str.eq."cc_converg_thresh").or. &
                     (desc_str.eq."cc_general").or. &
                     (desc_str.eq."cc_collective").or. &
                     (desc_str.eq."cc_use_disk").or. &
                     (desc_str.eq."cc_n_domain").or. &
                     (desc_str.eq."cc_solver").or. &
                     (desc_str.eq."cc_mpi_max_len").or. &
                     (desc_str.eq."cc_work_tmp_size").or. &
                     (desc_str.eq."cc_sv_control").or. &
                     (desc_str.eq."cc_restart").or. &
                     (desc_str.eq."cc_diis_step").or. &
                     (desc_str.eq."cc_diis_n_sv").or. &
                     (desc_str.eq."cc_init_level_shift").or. &
                     (desc_str.eq."cc_check").or. &
                     (desc_str.eq."cc_abcd_not_save").or. &
                     (desc_str.eq."cc_gamma").or. &
                     (desc_str.eq."cc_mem_check").or. &
                     (desc_str.eq."cc_super_system").or. &
                     (desc_str.eq."cc_pt_n_out").or. &
                     ! keyword for sampling RPA potentials along AC path
                     (desc_str.eq."rpa_along_ac_path").or. &
                     ! keywords for explict range separation scheme
                     (desc_str.eq."ers_width") .or. &
                     (desc_str.eq."ers_rcut") .or. &
                     (desc_str.eq."ers_omega") .or. &
                     (desc_str.eq."erfc_omega_postSCF") .or. &
                     ! keyword for skipping scf for post-scf calculation
                     (desc_str.eq."skip_scf_for_postscf") .or. &
                     (desc_str.eq."restart_for_postscf") .or. &
                     ! keyword for frozen virtual orbitals
                     (desc_str.eq."frozen_virtual_file") .or. &
                     !
                     (desc_str.eq."state_lower_limit").or. &
                     (desc_str.eq."state_upper_limit").or. &
                     (desc_str.eq."n_anacon_par").or. &
                     (desc_str.eq."anacon_type") .or.&
                     (desc_str.eq."set_vacuum_level") .or.  &
                     (desc_str.eq."use_dipole_correction") .or.&
                     (desc_str.eq."evaluate_dipole_correction_method") .or.&
                     (desc_str.eq."evaluate_work_function") .or. &
                     (desc_str.eq."use_small_component") .or. &
                     (desc_str.eq."mc_int") .or. &
                     (desc_str.eq."vdwdf") .or. &
                     (desc_str.eq."transport") .or.&
                     (desc_str.eq."vdw_cells") .or. &
                     (desc_str.eq."vdw_convergence_threshold") .or. &
                     (desc_str.eq."force_n_electrons") .or. &
                     (desc_str.eq."fixed_spin_moment") .or. &
                     (desc_str.eq."dos_kgrid_factors") .or. &
                     (desc_str.eq."override_relativity") .or. &
                     (desc_str.eq."override_warning_nobasis") .or. &
                     (desc_str.eq."override_warning_libxc") .or. &
                     (desc_str.eq."override_warning_negative_nucleus") .or. &
                     (desc_str.eq."override_illconditioning") .or.  &
                     (desc_str.eq."condition_penalty") .or.  &
                     (desc_str.eq."atomic_solver_xc")  .or. &
                     (desc_str.eq."calculation_set") .or. &
                     (desc_str.eq."calculation_subset") .or. &
                     !(desc_str.eq."output_locpot_atom") .or. &
                     (desc_str.eq."calculate_perturbative_soc") .or. &
                     (desc_str.eq."calculate_perturbative_soc_old_cluster") .or. &
                     (desc_str.eq."include_spin_orbit") .or. &
                     (desc_str.eq."include_spin_orbit_sc") .or. &
                     (desc_str.eq."n_core_states_omit_from_soc") .or. &
                     (desc_str.eq."min_energy_include_in_soc") .or. &
                     (desc_str.eq."min_energy_save_in_soc") .or. &
                     (desc_str.eq."gap_for_min_energy_in_soc") .or. &
                     (desc_str.eq."gap_for_saved_min_energy_in_soc") .or. &
                     (desc_str.eq."n_high_states_omit_from_soc") .or. &
                     (desc_str.eq."max_energy_include_in_soc") .or. &
                     (desc_str.eq."max_energy_save_in_soc") .or. &
                     (desc_str.eq."soc_subspace") .or. &
                     (desc_str.eq."soc_perturb_order") .or. &
                     (desc_str.eq."include_pw_lda_in_v_soc") .or. &
                     (desc_str.eq."calculate_all_eigenstates") .or. &
                     (desc_str.eq."atomic_solver") .or. &
                     (desc_str.eq."finite_nuclear_radius") .or. &
                     (desc_str.eq."magnetic_response") .or. &
                     (desc_str.eq."mr_experimental") .or. &
                     (desc_str.eq."dfpt_iter_limit") .or. &
                     (desc_str.eq."dfpt_linear_mix_param") .or. &
                     (desc_str.eq."dfpt_pulay_steps") .or. &
                     (desc_str.eq."output_sxml") .or. &
                     (desc_str.eq."mr_gauge_origin") .or. &
                     (desc_str.eq."mr_restart") .or. &
                     (desc_str.eq."dfpt_accuracy_n1") .or. &
                     (desc_str.eq."skip_SCF") .or.  &
                     (desc_str.eq."dry_run") .or.  &
                     (desc_str.eq."delta_numerical_stress") .or. &
                     (desc_str.eq."calculate_atom_bsse") .or.  &
                     (desc_str.eq."postprocess_anyway") .or. &
                     (desc_str.eq."nlcorr_nrad") .or. &
                     (desc_str.eq."nlcorr_i_leb") .or. &
                     (desc_str.eq."vdw_method") .or. &
                     (desc_str.eq."sc_self_energy").or.&
                     (desc_str.eq."n_poles").or. &
                     (desc_str.eq."pole_min").or. &
                     (desc_str.eq."pole_max").or. &
                     (desc_str.eq."dist_pole_type").or. &
                     (desc_str.eq."scgw_it_limit").or. &
                     (desc_str.eq."scgw_scf_accuracy").or. &
                     (desc_str.eq."scgw_print_all_spectrum").or. &
                     (desc_str.eq."scgw_mix_param").or. &
                     (desc_str.eq."use_aufbau_principle") .or.  &
                     (desc_str.eq."verbatim_writeout") .or. &
                     (desc_str.eq."distribute_leftover_charge").or. &
                     (desc_str.eq."compensate_multipole_errors").or. &
                     (desc_str.eq."normalize_initial_density").or. &
                     (desc_str.eq."use_old_cube_routine").or. &
                     (desc_str.eq."nomad_comment").or. &
                     (desc_str.eq."cube_units").or. &
                     (desc_str.eq."scgw_chem_pot").or. &
                     (desc_str.eq."use_aufbau_principle").or. &
                     (desc_str.eq."solvent").or. &
                     (desc_str.eq."mpe_nonelectrostatic_model").or. &
                     (desc_str.eq."mpe_xml_logging").or. &
                     (desc_str.eq."mpe_factorization_type").or. &
                     (desc_str.eq."mpe_solvent_permittivity").or. &
                     (desc_str.eq."mpe_skip_first_n_scf_steps") .or. &
                     (desc_str.eq."mpe_degree_of_determination").or. &
                     (desc_str.eq."mpe_f_sparsity_threshold") .or. &
                     (desc_str.eq."mpe_tol_adjR2") .or. &
                     (desc_str.eq."mpe_tol_adjR2_wait_scf").or. &
                     (desc_str.eq."mpe_lmax_rf") .or. &
                     (desc_str.eq."mpe_lmax_ep") .or. &
                     (desc_str.eq."mpe_n_centers_ep") .or. &
                     (desc_str.eq."mpe_n_boundary_conditions") .or. &
                     (desc_str.eq."isc_cavity_type").or. &
                     (desc_str.eq."isc_cavity_restart").or. &
                     (desc_str.eq."isc_cavity_restart_read").or. &
                     (desc_str.eq."isc_cavity_restart_write").or. &
                     (desc_str.eq."isc_record_cavity_creation") .or. &
                     (desc_str.eq."isc_rho_rel_deviation_threshold") .or. &
                     (desc_str.eq."isc_dt") .or. &
                     (desc_str.eq."isc_rho_k") .or. &
                     (desc_str.eq."isc_g_k") .or. &
                     (desc_str.eq."isc_rep_k") .or. &
                     (desc_str.eq."isc_gradient_threshold") .or. &
                     (desc_str.eq."isc_dynamics_friction") .or. &
                     (desc_str.eq."isc_update_nlist_interval") .or. &
                     (desc_str.eq."isc_max_dyn_steps") .or. &
                     (desc_str.eq."isc_calculate_surface_and_volume") .or. &
                     (desc_str.eq."isc_surface_curvature_correction") .or. &
                     (desc_str.eq."isc_kill_ratio") .or. &
                     (desc_str.eq."isc_try_restore_convergence") .or. &
                     (desc_str.eq."periodic_hf").or. &
                     (desc_str.eq."use_split_xc_gw").or. &
                     (desc_str.eq."neutral_excitation").or. &
                     (desc_str.eq."excited_states").or. &
                     (desc_str.eq."output_tddft_spectrum") .or. &
                     (desc_str.eq."tddft_spectrum_broaden") .or. &
                     (desc_str.eq."excited_mode") .or. &
                     (desc_str.eq."tddft_x") .or. &
                     (desc_str.eq."tddft_c") .or. &
                     (desc_str.eq."tddft_xc") .or. &
                     (desc_str == "xc_pre") .or. &
                     (desc_str.eq."casida_reduce_matrix") .or. &
                     (desc_str.eq."casida_reduce_occ") .or. &
                     (desc_str.eq."casida_reduce_unocc") .or. &
                     (desc_str.eq."casida_setup_scalapack") .or. &
                     (desc_str.eq."print_self_energy").or. &
                     (desc_str.eq."restrict_kpoint_to_smp_node").or. &
                     (desc_str.eq."max_tasks_per_smp_node") .or. &
                     (desc_str.eq."exx_band_structure_version") .or.&
                     (desc_str.eq."use_dipole_correction").or. &
                     (desc_str.eq."energy_density") .or. &
                     (desc_str.eq."exx_band_structure_version") .or. &
                     (desc_str.eq."compute_analytical_stress") .or. &
                     (desc_str.eq."sc_accuracy_stress") .or. &
                     (desc_str.eq."calc_analytical_stress_symmetrized") .or. &
                     (desc_str.eq."output_analytical_stress_symmetrized") .or. &
                     (desc_str.eq."external_pressure") .or. &
                     (desc_str.eq."MD_gle_A") .or. &
                     (desc_str.eq."MD_gle_C") .or. &
                     (desc_str.eq."stress_for_relaxation") .or. &
                     (desc_str.eq."analytic_potential_average") .or. &
                     (desc_str.eq."mbd_cfdm_dip_cutoff") .or. &
                     (desc_str.eq."mbd_scs_dip_cutoff").or.  &
                     (desc_str.eq."mbd_supercell_cutoff") .or. &
                     (desc_str.eq."mbd_scs_vacuum_axis") .or. &
                     (desc_str.eq."mbd_eigensolver") .or. &
                     (desc_str.eq."add_embedding_grids") .or. &
                     (desc_str.eq."output_KS_eigenvector_scalapack_hdf5") .or. &
                     (desc_str.eq."esp") .or. &
                     (desc_str.eq."compute_esp_charges") .or. &
                     (desc_str.eq."harmonic_potential_only") .or. &
                     (desc_str.eq."sc_abandon_etot") .or. &
                     (desc_str.eq."kerker_factor") .or. &
                     (desc_str.eq."use_integer_hash") .or. &
                     (desc_str.eq."calc_dens_superpos") .or. &
                     (desc_str.eq."use_potential_constrain") .or. &
                     (desc_str.eq."output_boys_centers") .or. &
                     (desc_str.eq."fast_Ewald") .or. &
                     ! For GPU acceleration
                     (desc_str.eq."use_gpu") .or. &
                     (desc_str.eq."gpu_density") .or. &
                     (desc_str.eq."gpu_hamiltonian") .or. &
                     (desc_str.eq."gpu_forces") .or. &
                     ! ELSI settings
                     (desc_str.eq."elsi_method") .or. &
                     (desc_str.eq."elsi_solver") .or. &
                     (desc_str.eq."elsi_output") .or. &
                     (desc_str.eq."elsi_restart") .or. &
                     (desc_str.eq."elsi_restart_use_overlap") .or. &
                     (desc_str.eq."elsi_rw_filter") .or. &
                     (desc_str.eq."elsi_filter") .or. &
                     (desc_str.eq."elsi_elpa_solver") .or. &
                     (desc_str.eq."elsi_elpa_n_single") .or. &
                     (desc_str.eq."elsi_elpa_gpu") .or. &
                     (desc_str.eq."elsi_omm_n_elpa") .or. &
                     (desc_str.eq."elsi_omm_flavor") .or. &
                     (desc_str.eq."elsi_omm_tol") .or. &
                     (desc_str.eq."elsi_pexsi_np_per_pole") .or. &
                     (desc_str.eq."elsi_pexsi_np_symbo") .or. &
                     (desc_str.eq."elsi_eigenexa_method") .or. &
                     (desc_str.eq."elsi_sips_n_slice") .or. &
                     (desc_str.eq."elsi_sips_slice_type") .or. &
                     (desc_str.eq."elsi_sips_n_elpa") .or. &
                     (desc_str.eq."elsi_ntpoly_method") .or. &
                     (desc_str.eq."elsi_ntpoly_tol") .or. &
                     (desc_str.eq."elsi_ntpoly_filter") .or. &
                     (desc_str.eq."elsi_magma_solver") .or. &
                     ! Block size
                     (desc_str.eq."scalapack_block_size") .or. &
                     (desc_str.eq."periodic_soc") .or. &
                     (desc_str.eq."check_cpu_consistency") .or. &
                     (desc_str.eq."cpu_consistency_threshold") .or. &
                     (desc_str.eq."force_complex_eigenvectors") .or. &
                     (desc_str.eq."gamma_cut_coulomb") .or. &
                     (desc_str.eq."use_gw_gamma_corr") .or. &
                     (desc_str.eq."qpe_multi_solu") .or. &
                     (desc_str.eq."gw_zshot") .or. &
                     (desc_str.eq."contour_eta") .or. &
                     (desc_str.eq."contour_zshot_offset") .or. &
                     (desc_str.eq."contour_spin_channel") .or. &
                     (desc_str.eq."contour_restart") .or. &
                     (desc_str.eq."pre_adjust_qp_relativistic") .or. &
                     (desc_str.eq."post_adjust_qp_relativistic") .or. &
                     (desc_str.eq."calc_spectral_func") .or. &
                     (desc_str.eq."spectral_func_state") .or. &
                     (desc_str.eq."full_cmplx_sigma") .or. &
                     (desc_str.eq."iterations_sc_cd") .or. &
                     (desc_str.eq."nocc_sc_cd") .or. &
                     (desc_str.eq."nvirt_sc_cd") .or. &
                     (desc_str.eq."sc_reiterate") .or. &
                     (desc_str.eq."try_zshot") .or. &
                     (desc_str.eq."gw_hedin_shift") .or. &
                     (desc_str.eq."debug_module") .or. &
                     (desc_str.eq."global_memory_tracking") .or. &
                     (desc_str.eq."switch_rho_multipole") .or. &
                     (desc_str.eq."use_symmetry_analysis") .or. &
                     (desc_str.eq."use_symmetric_forces") .or. &
                     (desc_str.eq."sym_precision") .or. &
                     (desc_str.eq."switch_rho_multipole") .or. &
                     ! For all DFPT
                     (desc_str .eq. "DFPT_mixing") .or. &
                     (desc_str .eq. "DFPT_sc_accuracy_dm") .or. &
                     (desc_str .eq. "DFPT_restart_read_DM1") .or. &
                     (desc_str .eq. "DFPT_restart_write_DM1") .or. &
                     (desc_str .eq. "DFPT_width") .or. &
                     ! For FODFT
                     (desc_str.eq."fo_dft") .or. &
                     (desc_str.eq."fo_orbitals") .or. &
                     (desc_str.eq."fo_options") .or. &
                     (desc_str.eq."fo_flavour") .or. &
                     (desc_str.eq."fo_folders") .or. &
                     (desc_str.eq."fo_verbosity") .or. &
                     (desc_str.eq."fo_deltaplus") .or. &
                     (desc_str.eq."force_single_restartfile") .or. &
                     (desc_str.eq."lc_dielectric_constant") .or. &
                     (desc_str == "xml_file") &
                     .or. (desc_str == 'python_hook') .or. &
                     ! Managing meta-GGA SCF convergence
                     (desc_str.eq."tau_threshold") .or. &
                     ! For write integrals, ovlp_3ks, ovlp_3fn, screened_coulomb
                     (desc_str.eq."write_integral") .or. &
                     ! parameters for bse
                     (desc_str.eq."read_write_qpe") .or. &
                     (desc_str.eq."bse_s_t") .or. &
                     (desc_str.eq."bse_lower_limit") .or. &
                     (desc_str.eq."bse_upper_limit") .or. &
                     (desc_str.eq."bse_reduce_matrix") .or. &
                     (desc_str.eq."bse_reduce_occ") .or. &
                     (desc_str.eq."bse_reduce_unocc") .or. &
                     ! parameters for band_mulliken
                     (desc_str.eq."band_mulliken_orbit_num") .or. &
                     ! parameters for dipole transition matrix, used with
                     ! compute_dielectric
                     (desc_str.eq."dipole_trans_mat") .or. &
                     (desc_str.eq."potential_well") .or.&
                     (desc_str.eq."rlsy_symmetry") &
                     ) then
   !           all other keywords do not matter
   !           BUT list them anyway - if you introduce an unknown keyword and the
   !           prg does not warn you, it will instead quit parsing the control.in
   !           file and give you the most bizarre segfaults - you really do not want that.
               continue
            else
               write(info_str,'(A)') &
                     "***************************** WARNING ******************************"
               call localorb_info(info_str)
               call localorb_info('* ')
               write(info_str,'(A,A,A)') &
                     "* Unknown keyword ", desc_str, " in input file control.in."
               call localorb_info(info_str)
               write(info_str,'(A)') &
                     "* Please check the syntax of this keyword - is it a typo?"
               call localorb_info(info_str)
               write(info_str,'(A)') &
                     "* If this is a valid keyword that was just introduced to the code, "
               call localorb_info(info_str)
               write(info_str,'(A)') &
                     "* please also make sure that the keyword is listed in subroutine parse_control()"
               call localorb_info(info_str)
               write(info_str,'(A)') &
                     "* in the dimensions.f90 module."
               call localorb_info(info_str)
               write(info_str,'(A)') &
                     "* Stopping the code for now."
               call localorb_info(info_str)
               call localorb_info('* ')
               write(info_str,'(A)') &
                     "********************************************************************"
               call localorb_info('')
               call aims_stop_coll('', '')
            endif
         enddo lineloop

         close(777)

         call localorb_info(" ", &
               use_unit,'(2X,A)' )
         call localorb_info("-----------------------------------------------------------------------", &
               use_unit,'(2X,A)' )
         call localorb_info("Completed first pass over input file control.in .", &
               use_unit,'(2X,A)' )
         call localorb_info("-----------------------------------------------------------------------", &
               use_unit,'(2X,A)' )
         call localorb_info(" ", &
               use_unit,'(2X,A)' )

!       verify number of listed species
         if (i_species.gt.0) then
            n_species = i_species
         else
            call aims_stop_coll("! No species in control.in.", func)
         endif

!       count number of pseudoized species
         if(i_pp_species .gt. 0) then
            use_embedding_pp = .true.
            n_pp_species = i_pp_species
         endif

!       if minimal basis was not switched off, increase number of basis function types
         if (use_min_basis) then
            n_basis_types = n_basis_types+1
         endif

!       Set the correct grid dimension for the Hartree multipole components
         if (use_hartree_log_grid) then
            n_hartree_grid = n_max_grid
         else
            n_hartree_grid = n_max_radial
         endif

         ! multiplicity only makes sense together with spin, check!
         if (use_hf_multiplicity  ) then
            if (n_spin.eq.2) then
               use_constraint = .true.
               n_region=n_region+1
            else
               use_hf_multiplicity = .false.
            endif
         endif

!       Verify that the angular division is sensible
!       and if it is not forced then default to one for allocations

         if ((use_angular_division).and. &
               (n_max_angular_division.le.0)) then
            write (info_str,'(1X,A,I4)') &
                  "* You have requested maximum of", &
                  n_max_angular_division
            call localorb_info(info_str)
            write (info_str,'(1X,A,/,2X,A,/,2X,A)') &
                  "* divisions of the angular shells. ", &
                  "This does not make sense. ", &
                  "Please, adjust control.in."
            call aims_stop_coll(info_str, func)
         elseif (.not.use_angular_division) then
            n_max_angular_division = 1
         endif

         if( use_hartree_fock .or. &
            use_hf_post .or. &
            use_gw .or. &
            use_gw_expt .or. &
            use_scgw .or. &
            use_scgw0 .or. &
            use_mp2.or. &
            use_os_mp2_qpe.or. &
            use_os_rpa_qpe.or. &
            use_sic_rpa_qpe.or. &
            use_dftpt2.or. &
            use_ci.or. &
            use_cc.or. &
            use_mp2sf.or. &
            use_rpa_ene.or. &
            use_C6_coef.or. &
            use_bare_ci .or. &
            use_hybrid_meta .or. &
            use_dmft_gw .or. &
            use_dmft_pbe0 ) then

            use_prodbas = .true.

         endif
         use_qpe = (use_gw .or. use_gw_expt .or. use_mp2sf .or. use_rpa_ene .or. use_C6_coef)
         use_corr = (use_qpe .or. use_scgw .or. use_scgw0 .or. use_screx .or. use_cohsex .or. &
         &           use_mp2 .or. use_bare_ci .or. use_dftpt2 .or. use_ci .or. use_os_mp2 .or. &
         &           use_os_mp2_qpe .or. use_os_rpa_qpe .or. use_sic_rpa_qpe .or. &
         &           use_cc .or. use_bse .or. use_coulomb_integral_ks)

!XR: I hope I didn't miss anything, otherwisse please add new entries here.
         if(use_corr .or. flag_neutral_excitation) flag_default_lvl = .false.
         if(use_prodbas .and. (.not. flag_RI)) then
             if( flag_default_lvl .or. n_periodic .gt. 0) then
                 use_lvl=.true.
                 use_lvl_fast=.true.
             else
                 use_lvl=.false.
                 use_lvl_fast=.false.
             endif
         endif

         return

      88 continue
         if (myid == 0) then
            write(use_unit,*) "Syntax error reading 'control.in'."
            write(use_unit,*) "The keyword is valid but it is missing one or more mandatory arguments."
            write(use_unit,*) "Check the input file, the manual or the source code in read_control.f90 if needed."
            write(use_unit,*) "line: '"//trim(inputline)//"'"
         end if
         call aims_stop_coll('Missing keyword argument in control.in', func)

      99 continue
         if (myid == 0) then
            write(use_unit,*) "Syntax error reading 'control.in'"
            write(use_unit,*) "The keyword is valid but one or more of its mandatory arguments could not be read."
            write(use_unit,*) "Check the input file, the manual or the source code in read_control.f90 if needed."
            write(use_unit,*) "line: '"//trim(inputline)//"'"
         end if
         call aims_stop_coll('Error reading keyword argument in control.in', func)

      end subroutine parse_control

!******
!---------------------------------------------------------------------
!****s* dimensions/parse_mpb
!  NAME
!    parse_mpb
!  SYNOPSIS

      subroutine parse_mpb  ( verbatim )


         use mpi_tasks
         use localorb_io

         implicit none

         character(*), parameter :: func = 'parse_mpb'
         character*40 :: desc_str
         character*132 :: inputline
	 integer :: i_code
	 logical :: verbatim

	 !we need the gradient of the density for MBPE calculations
	 use_density_gradient = .true.

         lineloop: do

            read(777,'(A)',iostat=i_code) inputline

            if(i_code<0) then ! EOF
               backspace(777)
               exit lineloop
            elseif(i_code>0) then
               ! handle error at main level
               backspace(777)
               return
            endif

	    if (verbatim) call localorb_info(inputline,use_unit,'(2X,A)')

!             if (verbatim) call localorb_info(inputline,use_unit,'(2X,A)')

            read(inputline,*,iostat=i_code) desc_str
            if (i_code.ne.0) then
               cycle  ! empty line
            elseif (desc_str(1:1).eq."#") then
               cycle  ! comment
            elseif (desc_str.eq."dielec_func".or.desc_str.eq."ions_temp".or.&
		desc_str.eq."ions_conc".or.desc_str.eq."ions_charge".or.&
		desc_str.eq."ions_kind".or.desc_str.eq."ions_size".or.&
		desc_str.eq."SPE_lmax".or.&
		desc_str.eq."lpbconstcoeff".or.desc_str.eq."not_converge_rho_mpb".or.&
		desc_str.eq."SPE_conv".or.desc_str.eq."dynamic_cavity_off".or.&
		desc_str.eq."start_mpb_solver".or.&
		desc_str.eq."delta_rho_in_sclpb".or.desc_str.eq."surface_and_volume_calc".or.&
		desc_str.eq."SPE_cut_and_lmax_ff".or.desc_str.eq."ions_mod_alpha".or.&
		desc_str.eq."solve_lpbe_only".or.desc_str.eq."reg_method".or.&
		desc_str.eq."MERM_in_SPE_solver".or.desc_str.eq."correct_dielectric_decrement".or.&
		desc_str.eq."set_nonelstat_params".or.desc_str.eq."dynamic_ions_off".or.&
		desc_str.eq."ions_mod_alpha_anion" .or. desc_str.eq."forces_mpb_force_off".or.&
		desc_str.eq."nonsc_Gnonmf".or.desc_str.eq."Gnonmf_FD_delta".or.&
		desc_str.eq."MERM_atom_wise") then
               continue
            else
!           must have reached end of mpb definitions
               backspace(777)
               exit lineloop
            endif

         enddo lineloop


	 return

      88 continue
         if (myid == 0) then
            write(use_unit,*) "Syntax error reading 'control.in' (missing arguments)"
            write(use_unit,*) "line: '"//trim(inputline)//"'"
         end if
         call aims_stop_coll('', func)

      99 continue
         if (myid == 0) then
            write(use_unit,*) "Syntax error reading 'control.in'"
            write(use_unit,*) "line: '"//trim(inputline)//"'"
         end if
         call aims_stop_coll('', func)

     end subroutine parse_mpb

!******
!---------------------------------------------------------------------
!****s* dimensions/parse_species
!  NAME
!    parse_species
!  SYNOPSIS

      subroutine parse_species  ( i_species, flag_pseudoized, verbatim )

!  PURPOSE

!  Subroutine parse_species parses only one species in the input file
!  control.in to determine the sizes of allocatable arrays, and does a few
!  integrity checks.
!
!  USES

         use constants, only: bohr
         use localorb_io
         use mpi_tasks

!  ARGUMENTS
!  INPUTS
!    o verbatim:   inherited from parse_control, governs if each read line should be written out
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

         implicit none


!  imported variables

         integer :: i_species
         logical :: verbatim

!  local variables

         integer :: i_code
         character*40 :: desc_str
         character*132 :: inputline
         real*8 :: nuclear_charge
         integer :: n_radial_sp
         real*8 :: scale_radial_sp
         integer :: radial_multiplier
         ! n_angular_sp will be the maximum number of angular grid points
         ! for the species in question
         integer :: n_angular_sp
         real*8 :: r_grid_min_sp
         real*8 :: r_grid_max_sp
         real*8 :: r_grid_inc_sp

         integer :: n_shell
         character :: l_shell_str
         integer :: l_shell

         integer :: n_contracted
         integer :: cartesian_l

         integer :: l_pot_sp

         logical :: use_min_basis_sp

         integer :: n_atomic_sp
         integer :: n_ionic_sp
         integer :: n_hydro_sp
!         integer :: n_hydro_sp_aux
         integer :: n_conf_sp
         integer :: n_gaussian_sp
         integer :: n_aux_gaussian_sp
         integer :: n_ext_basis_sp
         integer :: n_sto_sp

         real*8 :: r_radial_min

         integer :: n_grid_sp
         integer :: n_max_ang_shells_sp

!       auxiliary
         integer :: n_angular_new
         logical :: flag_verify
         integer :: n_ang_points


!   pseudopot infrastrcuture
         logical :: flag_pseudoized
         logical :: flag_nonlocal_core
         integer :: n_pseudofn
         integer :: n_points_pseudofn

!  counters

         integer :: i_contracted
         integer :: i_shell_plus_u

!  functions

         integer str_to_l

         character(*), parameter :: func = 'parse_species'

!  begin work

         nuclear_charge = 0.d0
         n_radial_sp = 0
         scale_radial_sp = 0.d0
         radial_multiplier = 1
         n_angular_sp = 0

         ! Set the same defaults for logarithmic grid here as in read_species_data.f90 .
         r_grid_min_sp = 0.0001
         r_grid_max_sp = 100.
         r_grid_inc_sp = 1.0123

         use_min_basis_sp = .true.

         n_atomic_sp = 0
         n_ionic_sp = 0
         n_conf_sp = 0
         n_hydro_sp = 0
         n_gaussian_sp = 0
         n_aux_gaussian_sp = 0
         n_ext_basis_sp = 0
         n_sto_sp = 0

         n_max_ang_shells_sp = 0

         i_shell_plus_u = 0

         flag_pseudoized = .false.
         flag_nonlocal_core = .false.

         n_pseudofn = 0
         n_points_pseudofn = 0

         lineloop: do

            read(777,'(A)',iostat=i_code) inputline

            if(i_code<0) then ! EOF
               backspace(777)
               exit lineloop
            elseif(i_code>0) then
               ! handle error at main level
               backspace(777)
               return
            endif

            if (verbatim) call localorb_info(inputline,use_unit,'(2X,A)')

            read(inputline,*,iostat=i_code) desc_str
            if (i_code.ne.0) then
               cycle  ! empty line
            elseif (desc_str(1:1).eq."#") then
               cycle  ! comment
            elseif (desc_str.eq."pseudo") then
               flag_pseudoized = .true.
               read(inputline,*,end=88,err=99) desc_str, desc_str
               desc_str = trim(desc_str)
               if(len_trim(desc_str).eq.0) then
                 if (myid.eq.0) then
                    write(use_unit,*) "AIMS is missing a link to the pseudopotential file (*.cpi) of species", &
                      i_species,"."
                    write(use_unit,*) "! Please specify a link in your control.in file before continuing."
                    call aims_stop_coll('', func)
                 endif
               endif
               call parse_pseudopot_dimensions( desc_str, n_pseudofn, &
                                            n_points_pseudofn)

               if(n_pseudofn .gt. n_max_pp_fn) then
                  n_max_pp_fn = n_pseudofn
               endif

               n_pp_basis_fns = n_pp_basis_fns + n_pseudofn

               if(n_points_pseudofn .gt. n_max_points_pp_fn) then
                  n_max_points_pp_fn = n_points_pseudofn
               endif

            elseif (desc_str.eq."nonlinear_core") then
               read(inputline,*,end=88,err=99) desc_str, flag_nonlocal_core
               if(flag_nonlocal_core) use_nonlinear_core = .true.

            elseif (desc_str.eq."nucleus") then
               read(inputline,*,end=88,err=99) desc_str, nuclear_charge
            elseif (desc_str.eq."radial_base") then
               read(inputline,*,end=88,err=99) desc_str, n_radial_sp, scale_radial_sp
               scale_radial_sp = scale_radial_sp/bohr
!           convert scale_radial from maximal radius to scaling factor ...
               scale_radial_sp = - scale_radial_sp / &
               log( 1.d0 - &
                  (n_radial_sp/(1.d0+n_radial_sp))**2.d0 )
            elseif (desc_str.eq."radial_multiplier") then
               read(inputline,*,end=88,err=99) desc_str, radial_multiplier
            elseif (desc_str.eq."angular") then
               read(inputline,*,end=88,err=99) desc_str, n_ang_points
               if (n_ang_points.gt.n_angular_sp) then
                  n_angular_sp = n_ang_points
               endif
            elseif (desc_str.eq."logarithmic") then
               read(inputline,*,end=88,err=99) desc_str, r_grid_min_sp, r_grid_max_sp, r_grid_inc_sp
            elseif (desc_str.eq."l_hartree") then
               read(inputline,*,end=88,err=99) desc_str, l_pot_sp
               l_pot_max = max(l_pot_max, l_pot_sp)
            elseif (desc_str.eq."core_states") then
            elseif (desc_str.eq."KH_core_states") then
               use_density_gradient = .true.
               use_basis_gradients = .true.
            elseif (desc_str.eq."include_min_basis") then
               read(inputline,*,end=88,err=99) desc_str, use_min_basis_sp
            elseif (desc_str.eq."valence") then
               read(inputline,*,end=88,err=99) desc_str,n_shell,l_shell_str
               l_shell = str_to_l(l_shell_str)
               if (n_wave_max.lt.n_shell) then
                  n_wave_max = n_shell
               endif
               if (l_wave_max.lt.l_shell) then
                  l_wave_max = l_shell
               endif
               if(use_relativistic_basis) then
              ! For relativistic cases, more basis functions are needed because
              ! of the splitting of orbitals (p -> p_1/2 and p_3/2; d -> d_3/2
              ! and d_5/2; f -> f_5/2 and f7/2).
                  if (l_shell.gt.0) then
                     n_atomic_sp = n_atomic_sp + 2* (n_shell-l_shell)
                  else
                     n_atomic_sp = n_atomic_sp + (n_shell-l_shell)
                  endif
               else
                  n_atomic_sp = n_atomic_sp + (n_shell-l_shell)
               endif
            elseif (desc_str.eq."ion_occ") then
               read(inputline,*,end=88,err=99) desc_str,n_shell,l_shell_str
               l_shell = str_to_l(l_shell_str)
               if (n_wave_max.lt.n_shell) then
                  n_wave_max = n_shell
               endif
               if (l_wave_max.lt.l_shell) then
                  l_wave_max = l_shell
               endif
            elseif (desc_str.eq."ionic") then
               read(inputline,*,end=88,err=99) desc_str,n_shell,l_shell_str
               if (.not.use_ionic) then
                  use_ionic = .true.
                  n_basis_types = n_basis_types+1
               endif
               l_shell = str_to_l(l_shell_str)
               if (l_wave_max.lt.l_shell) then
                  l_wave_max = l_shell
               endif
               if(use_relativistic_basis) then
                  if (l_shell.gt.0) then
                     n_ionic_sp = n_ionic_sp + 2
                  else
                     n_ionic_sp = n_ionic_sp + 1
                  endif
               else
                  n_ionic_sp = n_ionic_sp+1
               endif
            elseif (desc_str.eq."confined") then
               read(inputline,*,end=88,err=99) desc_str,n_shell,l_shell_str
               if (.not.use_confined) then
                  use_confined = .true.
                  n_basis_types = n_basis_types+1
               endif
               l_shell = str_to_l(l_shell_str)
               if (l_wave_max.lt.l_shell) then
                  l_wave_max = l_shell
               endif
               if(use_relativistic_basis) then
                  if (l_shell.gt.0) then
                     n_conf_sp = n_conf_sp + 2
                  else
                     n_conf_sp = n_conf_sp + 1
                  endif
               else
                  n_conf_sp = n_conf_sp+1
               endif
            elseif (desc_str.eq."hydro") then
               read(inputline,*,end=88,err=99) desc_str,n_shell,l_shell_str
               if (.not.use_hydro) then
                  use_hydro = .true.
                  n_basis_types = n_basis_types+1
               endif
               l_shell = str_to_l(l_shell_str)
               if (l_wave_max.lt.l_shell) then
                  l_wave_max = l_shell
               endif
               if(use_relativistic_basis) then
                  if (l_shell.gt.0) then
                     n_hydro_sp = n_hydro_sp + 2
                  else
                     n_hydro_sp = n_hydro_sp + 1
                  endif
               else
                  n_hydro_sp = n_hydro_sp+1
               endif
            elseif (desc_str.eq."pure_gauss") then
               read(inputline,*,end=88,err=99) desc_str, use_sph_gaussian
            elseif (desc_str.eq."gaussian") then
               read(inputline,*,end=88,err=99) desc_str,cartesian_l,n_contracted
               if (.not.use_gaussian) then
                  use_gaussian = .true.
                  n_basis_types = n_basis_types+1
               endif
               if (n_max_contracted.lt.n_contracted) then
                  n_max_contracted = n_contracted
               endif
               if (l_wave_max.lt.cartesian_l) then
                  l_wave_max = cartesian_l
               endif
               n_gaussian_sp = n_gaussian_sp &
                           + floor(real(cartesian_l)/2.) + 1
               if (n_contracted.gt.1) then
                  do i_contracted = 1, n_contracted, 1
                     read(777,'(A)',iostat=i_code) inputline
                     if(i_code<0)then
                        write(use_unit,*)"ERROR reading 'config.in': EOF within list of Gaussian contractions"
                        stop
                     endif
                     if(i_code>0)then
                        ! handle error at main level
                        backspace(777)
                        exit lineloop
                     endif
                  call localorb_info(inputline) !write out the contracted factors
                  enddo
               endif

            elseif (desc_str.eq."aux_gaussian") then
               read(inputline,*,end=88,err=99) &
               desc_str,cartesian_l,n_contracted
               if (.not.use_aux_gaussian) then
                  use_aux_gaussian = .true.
               endif
               if (n_max_aux_contracted.lt.n_contracted) then
                  n_max_aux_contracted = n_contracted
               endif
!            if (l_wave_max.lt.cartesian_l) then
!              l_wave_max = cartesian_l
!            endif
               n_aux_gaussian_sp = n_aux_gaussian_sp &
                           + floor(real(cartesian_l)/2.) + 1
               if (n_contracted.gt.1) then
                  do i_contracted = 1, n_contracted, 1
                     read(777,'(A)',iostat=i_code) inputline
                     if(i_code<0)then
                        write(use_unit,*)"ERROR reading 'config.in': EOF within list of Gaussian contractions"
                        stop
                     endif
                     if(i_code>0)then
                        ! handle error at main level
                        backspace(777)
                        exit lineloop
                     endif
                  enddo
               endif

            elseif (desc_str.eq."for_aux") then
! Read information for extra basis basis functions for auxiliary RI basis
               read(inputline,*,end=88, err=99) desc_str, desc_str, n_shell, l_shell_str
!               write(use_unit,*) "parse ext basis:",  desc_str, n_shell, l_shell_str
               if (.not.use_ext_basis) then
                  use_ext_basis = .true.
               endif

               if (desc_str.eq."ionic") then
                    if (.not.use_ionic) then
                        use_ionic = .true.
                        n_basis_types = n_basis_types+1
                    endif
                    l_shell = str_to_l(l_shell_str)
                    if (l_wave_max.lt.l_shell) then
                        l_wave_max = l_shell
                    endif
                    n_ionic_sp = n_ionic_sp+1
               elseif (desc_str.eq."hydro") then
                    if (.not.use_hydro) then
                        use_hydro = .true.
                        n_basis_types = n_basis_types+1
                    endif
                    l_shell = str_to_l(l_shell_str)
                    if (l_wave_max.lt.l_shell) then
                        l_wave_max = l_shell
                    endif
                    n_hydro_sp = n_hydro_sp+1
               end if
!---------------------

            else if (desc_str == 'sto') then
               read(inputline, *, end=88, err=99) desc_str, n_shell, l_shell
               if (.not. use_sto) then
                  use_sto = .true.
                  n_basis_types = n_basis_types + 1
               endif
               n_wave_max = max(n_wave_max, n_shell)
               l_wave_max = max(l_wave_max, l_shell)
               if(use_relativistic_basis) then
                  if (l_shell.gt.0) then
                     n_sto_sp = n_sto_sp + 2
                  else
                     n_sto_sp = n_sto_sp + 1
                  endif
               else
                  n_sto_sp = n_sto_sp + 1
               endif

            elseif (desc_str.eq."angular_grids") then
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if (desc_str.eq."specified") then
                  use_specified_grid = .true.
               endif

            elseif (desc_str.eq."division") then
               if (use_specified_grid) then

                  read(inputline,*,end=88,err=99) desc_str, desc_str, n_ang_points
                  if (n_ang_points.gt.n_angular_sp) then
                     n_angular_sp = n_ang_points
                  endif
                  n_max_ang_shells_sp = n_max_ang_shells_sp+1

               endif
            elseif (desc_str.eq."outer_grid") then
               if (use_specified_grid) then

                  read(inputline,*,end=88,err=99) desc_str, n_ang_points
                  if (n_ang_points.gt.n_angular_sp) then
                     n_angular_sp = n_ang_points
                  endif

               endif

            elseif ((desc_str .eq. "plus_u") .or. (desc_str .eq. "plus_U")) then
               use_plus_u = .true.
               i_shell_plus_u = i_shell_plus_u + 1

            elseif &
               ( (desc_str.eq."basis_acc") .or. &
               (desc_str.eq."prodbas_acc") .or. &
               (desc_str.eq."angular_acc") .or. &
               (desc_str.eq."angular_min") .or. &
               (desc_str.eq."fixed_grid") .or. &
               (desc_str.eq."innermost_max") .or. &
               (desc_str.eq."multipole_radius") .or. &
               (desc_str.eq."core") .or. &
               (desc_str.eq."cut_pot") .or. &
               (desc_str.eq."cutoff_type") .or. &
               (desc_str.eq."outer_grid") .or. &
               (desc_str.eq."cut_free_atom") .or. &
               (desc_str.eq."cut_core") .or. &
               (desc_str.eq."basis_dep_cutoff") .or. &
               (desc_str.eq."cut_atomic_basis") .or. &
               (desc_str.eq."max_n_prodbas") .or. &
               (desc_str.eq."max_l_prodbas") .or. &
               (desc_str.eq."mass") .or. &
               (desc_str.eq."pp_charge") .or. &
               (desc_str.eq."pp_local_component") .or. &
               (desc_str.eq."cite_reference") .or. &
               (desc_str.eq."hirshfeld_param") &
               .or. (desc_str == 'element') &
               ) then
!           all other species keywords do not matter
               continue
            else
!           must have reached end of current species
               backspace(777)
               exit lineloop
            endif

         enddo lineloop

!       now verify input data

!       radial integration grid
         if ( (n_radial_sp.gt.0) .and. &
               (scale_radial_sp.gt.0.0) .and. &
               (radial_multiplier.ge.1) ) then
!     increment number of radial grid points according to radial_multiplier
	    n_radial_sp = radial_multiplier * (n_radial_sp+1) - 1
	    if (n_max_radial.lt.n_radial_sp) then
	      n_max_radial = n_radial_sp
	    endif
         else
            if (myid.eq.0) then
               write(use_unit,*) "! Ill-defined radial grid in species ", &
                     i_species,"."
               write(use_unit,*) "! Alternatively, this could be a typo in the last found keyword: "
               write(use_unit,*) "! ", desc_str
               write(use_unit,*) "! Please check your control.in file before continuing."
            endif
            call aims_stop_coll('', func)
         endif

!       angular integration grid
!       make sure angular integration grid matches the requirements of
!       Hartree potential
         call verify_angular_grid &
         ( n_angular_sp, l_pot_sp, n_angular_new, flag_verify )
         if (flag_verify) then
            n_angular_sp = n_angular_new
         endif

         if (n_angular_sp.gt.0) then
            if (n_max_angular.lt.n_angular_sp) then
               n_max_angular = n_angular_sp
            endif
         else
            if (myid.eq.0) then
               write(use_unit,*) "! Ill-defined angular grid in species ", &
                     i_species,"."
               write(use_unit,*) "! Alternatively, this could be a typo in the last found keyword: "
               write(use_unit,*) "! ", desc_str
               write(use_unit,*) "! Please check your control.in file before continuing."
            endif
            call aims_stop_coll('', func)
         endif

         if (use_specified_grid) then
            if ( (n_max_ang_shells_sp+1) .gt. n_max_ang_shells ) then
               n_max_ang_shells = n_max_ang_shells_sp+1
            endif
         endif

!       logarithmic integration grid
!       This is complicated by the fact that r_grid_min is adjusted
!       to the elements, and the actual innermost point of the
!       radial integration grid. We now have an unfortunate implicit
!       interdependence with module grids.

! DB 101512: the nuclear_charge check is done in read_species_data once again.
!            We remove this check here, since we need to allow nuclear_charge e.q.
!            for the species "Emptium"

!         if ((nuclear_charge .le. 0.0) .and. (.not.(flag_pseudoized))) then
!            if (myid.eq.0) then
!               write(use_unit,*) "! Ill-defined nuclear charge in species ", &
!                     i_species,"."
!               write(use_unit,*) "! Alternatively, this could be a typo in the last found keyword: "
!               write(use_unit,*) "! ", desc_str
!               write(use_unit,*) "! Please check your control.in file before continuing."
!            endif
!            call aims_stop_coll('', func)
!         endif



         if ( (r_grid_min_sp.gt.0.0) .and. &
               (r_grid_max_sp.gt.r_grid_min_sp) .and. &
               (r_grid_inc_sp.gt.1.0) ) then

!         reduction due to nuclear charge
            r_grid_min_sp = r_grid_min_sp/nuclear_charge

!         verify that innermost logarithmic grid point is inside the
!         innermost radial integration shell, exactly as in module grids.
            r_radial_min = dble(1/(n_radial_sp+1.0d0))
            r_radial_min = -log(1.d0-r_radial_min**2.0d0)
            r_radial_min = scale_radial_sp*r_radial_min

            if (r_grid_min_sp.gt.r_radial_min) then
!           Grid will be modified just like this in grids.f - unfortunately
!           this is now hardwired in two different places instead of being
!           accessible as one subroutine! Be very careful!

               r_grid_min_sp = 0.5d0*r_radial_min

            endif

            n_grid_sp = int ( log(r_grid_max_sp/r_grid_min_sp) / &
                              log(r_grid_inc_sp) ) + 1

            if (n_max_grid.lt.n_grid_sp) then
               n_max_grid = n_grid_sp
            endif

            if (first_species) then
               log_r_min = r_grid_min_sp
               log_r_max = r_grid_max_sp
               log_inc_min = r_grid_inc_sp
               first_species = .false.
            else
               log_r_min = min(log_r_min, r_grid_min_sp)
               log_r_max = max(log_r_max, r_grid_max_sp)
               log_inc_min = min(log_inc_min, r_grid_inc_sp)
            endif

         else
            if (myid.eq.0) then
               write(use_unit,*) "! Ill-defined logarithmic grid in species ", &
                     i_species,"."
               write(use_unit,*) "! Alternatively, this could be a typo in the last found keyword: "
               write(use_unit,*) "! ", desc_str
               write(use_unit,*) "! Please check your control.in file before continuing."
            endif
            call aims_stop_coll('', func)
         endif

!  Allocation decision for minimal ("atomic") basis fn arrays

         use_min_basis = use_min_basis.or.use_min_basis_sp

!  Array dimension for individual basis function types

         if (n_atomic_sp.gt.n_max_ind_fns) then
            n_max_ind_fns = n_atomic_sp
         endif
         if (n_ionic_sp.gt.n_max_ind_fns) then
            n_max_ind_fns = n_ionic_sp
         endif
         if (n_conf_sp.gt.n_max_ind_fns) then
            n_max_ind_fns = n_conf_sp
         endif
         if (n_hydro_sp.gt.n_max_ind_fns) then
            n_max_ind_fns = n_hydro_sp
         endif
         if (n_gaussian_sp.gt.n_max_ind_fns) then
            n_max_ind_fns = n_gaussian_sp
         endif
         n_max_ind_fns = max(n_max_ind_fns, n_sto_sp)

         if (n_aux_gaussian_sp.gt.n_max_aux_fns) then
            n_max_aux_fns = n_aux_gaussian_sp
         endif


!  Maximal number of (n,l) shells per species which are treated with LDA(GGA)+U

         if (i_shell_plus_u .gt. n_max_shells_plus_u) then
            n_max_shells_plus_u = i_shell_plus_u
         endif

         return

      88 continue
         if (myid == 0) then
            write(use_unit,*) "Syntax error reading 'control.in' (missing arguments)"
            write(use_unit,*) "line: '"//trim(inputline)//"'"
         end if
         call aims_stop_coll('', func)

      99 continue
         if (myid == 0) then
            write(use_unit,*) "Syntax error reading 'control.in'"
            write(use_unit,*) "line: '"//trim(inputline)//"'"
         end if
         call aims_stop_coll('', func)

      end subroutine parse_species
!******
!---------------------------------------------------------------------
!****s* dimensions/parse_geometry
!  NAME
!   parse_geometry
!  SYNOPSIS

      subroutine parse_geometry

!  PURPOSE

!  Subroutine parse_geometry parses the input file geometry.in to
!  determine the sizes of allocatable arrays, and does a few
!  integrity checks.
!
!  USES

         use localorb_io
         use mpi_tasks, only: aims_stop, aims_stop_coll, myid

!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

         implicit none

   !  imported variables

   !  local variables

         integer :: i_code
         character*40 :: desc_str
         character*132 :: inputline

         logical :: constrain_relaxation

         real*8 :: initial_charge
         logical :: verbatim_writeout = .true.

   !  counters
   !  Reset below - this seems to work better in the library setup
         integer :: i_atom != 0
         integer :: i_empty_atoms !=0
         integer :: i_periodic != 0
   !       Nadia
         integer :: i_multipole != 0
   !       PSPs
         integer :: i_pp_atom != 0

         character(*), parameter :: func = 'parse_geometry'
   !  begin work

         call localorb_info(" ", &
               use_unit,'(2X,A)' )

         call localorb_info("-----------------------------------------------------------------------", &
               use_unit,'(2X,A)' )

         call localorb_info("Parsing geometry.in (first pass over file, find array dimensions only).", &
               use_unit,'(2X,A)' )

         call localorb_info("The contents of geometry.in will be repeated verbatim below", &
               use_unit,'(2X,A)' )

         call localorb_info("unless switched off by setting 'verbatim_writeout .false.' .", &
               use_unit,'(2X,A)' )

         call localorb_info("in the first line of geometry.in .", &
               use_unit,'(2X,A)' )

         call localorb_info("-----------------------------------------------------------------------", &
               use_unit,'(2X,A)' )

         call localorb_info(" ", &
               use_unit,'(2X,A)' )


! Reset all the counters
         i_atom = 0
         i_empty_atoms = 0
         i_periodic = 0
         i_multipole = 0
         i_pp_atom = 0

   !       Nadia
   !       Default is .false.
         if (.not.use_atom_ref) then
            use_embedding_potential = .false.
         endif

         use_relaxation_constraints = .false.
         constrain_relaxation = .false.
         n_hess_blocks = 0
         n_hess_blocks_lv       = 0
         n_hess_blocks_lv_atom  = 0

         open(888,file="geometry.in",status='OLD',iostat=i_code)

         if (i_code.ne.0) then
            call aims_stop("* Input file geometry.in not found.")
         endif

         ! counts number of entries (atoms, lattice_vectors ...)
         lineloop: do
            read(888,'(A)',iostat=i_code) inputline
            if(i_code<0) exit lineloop      ! end of file reached
            if(i_code>0) then
               call aims_stop_coll("Unknown error reading file 'geometry.in'...", func)
            endif

            if (verbatim_writeout) call localorb_info(inputline,use_unit,'(2X,A)') !Write input to output

            read(inputline,*,iostat=i_code) desc_str
            if (i_code.ne.0) cycle ! skip empty line
            if (desc_str(1:1)=='#') cycle ! skip comment

            if ((desc_str.eq."atom") .or. (desc_str.eq."atom_frac")) then
               i_atom = i_atom+1
            elseif(desc_str.eq.'verbatim_writeout') then
                  read(inputline,*,end=88,err=99) desc_str, verbatim_writeout

            elseif (desc_str.eq."pseudocore") then
               i_pp_atom = i_pp_atom+1
               use_embedding_pp = .true.
            elseif (desc_str.eq."empty") then
               i_empty_atoms = i_empty_atoms+1
            elseif (desc_str.eq."velocity") then
               MD_velocity_in_input = .true.    ! just remember that a nonzero number of v's have been specified; need that in read_control
            elseif (desc_str.eq."lattice_vector") then
               i_periodic = i_periodic+1
            elseif (desc_str.eq."initial_moment") then
               if (n_spin.eq.2) then
               use_initial_rho = .true.
               use_initial_moment = .true.
               endif
            elseif (desc_str.eq."initial_charge") then
               read(inputline,*,end=88,err=99) desc_str, initial_charge
               if (initial_charge.ne.0.d0) then
                  use_initial_rho = .true.
               endif
            elseif (desc_str.eq."homogeneous_field") then
                !TODO: Remove testwise change
                if (n_periodic.eq.0) use_embedding_potential = .true.
   !            use_embedding_potential = .true.

            elseif(desc_str.eq."image_potential") then !testwise
                    use_embedding_potential = .true.
   !         Nadia
            elseif (desc_str.eq."multipole") then
               use_embedding_potential = .true.
               i_multipole = i_multipole+1
            elseif (desc_str.eq."constrain_relaxation") then
               if (use_geo_relaxation.or.use_molecular_dynamics) then
                  read(inputline,*,end=88,err=99) desc_str, desc_str
                  if (desc_str.ne.".false.") then
                     use_relaxation_constraints = .true.
                  endif
               endif
            elseif (desc_str.eq."calculate_friction") then
            elseif (desc_str.eq."hessian_block") then
               n_hess_blocks = n_hess_blocks + 1
            else if (desc_str.eq."hessian_block_lv") then
               n_hess_blocks_lv = n_hess_blocks_lv + 1
            else if (desc_str.eq."hessian_block_lv_atom") then
               n_hess_blocks_lv_atom = n_hess_blocks_lv_atom + 1
            else if (desc_str.eq."hessian_file") then
               hess_in_file = .true.
            ! read in number of parameters for symmetry-constrained relaxation
            else if (desc_str.eq."symmetry_n_params") then
               backspace(888)
               read(888,*,err=99) desc_str, SCR_n_params, SCR_n_params_lv, SCR_n_params_coords
               use_symm_const_geo_relaxation = .true.
               if ( SCR_n_params_lv .eq. 0 ) then
                 relax_unit_cell = 0
               end if
            else
               continue  ! silently ignore anything else
            endif
         enddo lineloop

         close(888)

         call localorb_info(" ", &
               use_unit,'(2X,A)' )
         call localorb_info("-----------------------------------------------------------------------", &
               use_unit,'(2X,A)' )
         call localorb_info("Completed first pass over input file geometry.in .", &
               use_unit,'(2X,A)' )
         call localorb_info("-----------------------------------------------------------------------", &
               use_unit,'(2X,A)' )
         call localorb_info(" ", &
               use_unit,'(2X,A)' )

   !       verify number of atoms
         if (i_atom.gt.0) then
            n_atoms = i_atom
         else
            call localorb_info("* No atoms in geometry.in.")
            call localorb_info("* Please check your input file.")
            call aims_stop_coll('', func)
         endif

!       BSSE
         n_occ_atoms=i_atom
         n_empty_atoms=i_empty_atoms
         n_atoms=i_atom+i_empty_atoms

!       verify number of lattice vector
         if(i_periodic >  0)then
            if (i_periodic<3) then
               call localorb_info("* Less than three lattice vectors in geometry.in.")
               call localorb_info("* Please check your input file.")
               call aims_stop_coll('', func)
            elseif(i_periodic>3)then
               call localorb_info("* More than three lattice vectors in geometry.in.")
               call localorb_info("* Please check your input file.")
               call aims_stop_coll('', func)
            else
               n_periodic = 3
            endif
         endif

         ! embedding with periodic boundary conditions must be forbidden.
          if ( (n_periodic.gt.0).and.i_multipole.gt.0) then !Changed testwise to allow fields w/ PBCs
            call localorb_info("* Attention. An external embedding field (multipoles or homogeneous)")
            call localorb_info("* was specified for a periodic system.")
            call localorb_info("* At present, this combination is simply not supported.")
            call localorb_info("* Implementing the missing part (Ewald summation) can be done,")
            call localorb_info("* but please talk to us beforehand.")
            call localorb_info("* Aborting the present run to alert you to the problem.")
            call aims_stop_coll('', func)
         endif

          if ( (n_periodic.gt.0).and.(i_pp_atom.gt.0)) then
            call localorb_info("* Attention. An pseudocore was specified for a periodic system.")
            call localorb_info("* At present, this combination is simply not supported.")
            call localorb_info("* Implementing the missing part (Ewald summation) can be done,")
            call localorb_info("* but please talk to us beforehand.")
            call localorb_info("* Aborting the present run to alert you to the problem.")
            call aims_stop_coll('', func)
         endif

         if ( (n_periodic.gt.0).and.use_embedding_potential) then
            call localorb_info("* Attention. An external embedding field")
            call localorb_info("* was specified for a periodic system.")
            call localorb_info("* At present, this combination is experimental,")
            call localorb_info("* and incorrect results could be obtained.")
            call localorb_info("* Please proceed with caution.")
         endif

         n_real_atoms = n_atoms

         n_multipoles = 0
         if (use_embedding_potential) then
            n_multipoles = i_multipole
         endif

         n_pp_atoms = i_pp_atom

         if(n_periodic > 0) then
            use_periodic_hf = use_hartree_fock
            use_hf_realspace = use_periodic_hf .and. (.not. use_hf_kspace)
         else
            use_periodic_hf = use_hartree_fock .and. use_lvl_fast
            if(use_periodic_hf) then
              if(use_hf_kspace) call aims_stop_coll('use_hf_kspace is only implemented for periodic systems', func)
              use_hf_realspace = .true.
            endif
         endif

         use_cutCb = use_periodic_hf
         use_cutCb = use_cutCb .or. (n_periodic > 0 .and. (use_rpa_ene .or. use_gw .or. use_mp2)) .or. flag_out_coulmat_lvl

         ! for periodic Hartree-Fock, we use LVL unless a specific other choice was made.
         if (.not. flag_RI) then
            if (use_periodic_hf .or. (n_periodic>0 .and. use_corr)) then
               use_lvl = .true.
            ! ... but no else. "Else" we leave it set the way it came in.
            end if
         end if

         ! MOL: needs allocation before reading geometry file
         if (use_symm_const_geo_relaxation) then
            allocate(SCR_params(SCR_n_params))
            allocate(SCR_lv_str(3,n_periodic))
            allocate(SCR_coord_str(3,n_atoms))
         end if

     return

      88 continue
         if (myid == 0) then
            write(use_unit,*) "Syntax error reading 'geometry.in' (missing arguments)"
            write(use_unit,*) "line: '"//trim(inputline)//"'"
         end if
         call aims_stop_coll('', func)

      99 continue
         if (myid == 0) then
            write(use_unit,*) "Syntax error reading 'geometry.in'"
            write(use_unit,*) "line: '"//trim(inputline)//"'"
         end if
         call aims_stop_coll('', func)



      end subroutine parse_geometry
!******
!---------------------------------------------------------------------
!****s* dimensions/parse_vdw
!  NAME
!   parse_vdw
!  SYNOPSIS

      subroutine parse_vdw

!  PURPOSE
!  Subroutine parse_vdw parses the vdW correction definitions in the input
!  file control.in to determine the sizes of allocatable arrays, and does a
!  few integrity checks.
!
!  USES

         use localorb_io
         use mpi_tasks

!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

         implicit none

! local variables
         integer :: i_code
         character*40 :: desc_str
         character*132 :: inputline
         real*8 :: C6
         real*8 :: R0
         real*8 :: d
         character*2  :: atom1
         character*2  :: atom2

   !  counters
         integer :: i_pair

         character(*), parameter :: func = 'parse_vdw'

         read(777,'(A)',iostat=i_code) inputline
         if(i_code/=0) then
            ! handle error at main level
            backspace(777)
            return
         endif

         read(inputline,*,end=88,err=99) desc_str, vdw_pairs
         if (desc_str/="vdw_pairs") then
            call aims_stop_coll("* Error in control.in: expected vdw_pairs line", func)
         endif

         do i_pair = 1, vdw_pairs, 1
            read(777,'(A)',iostat=i_code) inputline
            if(i_code<0)then
               call aims_stop_coll("* Incorrect number of vdw_coeff lines (found EOF).", func)
            endif
            if(i_code>0)then
               ! error handled at top level
               backspace(777)
               return
            endif

            read(inputline,*,iostat=i_code) desc_str,atom1,atom2,C6,R0,d
            if (i_code/=0 .or. desc_str /= "vdw_coeff") then
               call aims_stop_coll("* Incorrect number of vdw_coeff lines.", func)
            endif
         enddo

         return

      88 continue
         if (myid == 0) then
            write(use_unit,*) "Syntax error reading 'control.in' (missing arguments)"
            write(use_unit,*) "line: '"//trim(inputline)//"'"
         end if
         call aims_stop_coll('', func)

      99 continue
         if (myid == 0) then
            write(use_unit,*) "Syntax error reading 'control.in'"
            write(use_unit,*) "line: '"//trim(inputline)//"'"
         end if
         call aims_stop_coll('', func)

      end subroutine parse_vdw
!******
!---------------------------------------------------------------------


      subroutine parse_pseudopot_dimensions( pp_path, n_func, &
                                            n_func_points)

      !Purpose: parses number of gridpoints the pseudopotentials are defined on and
      !         the number of pseudowavefunction of species pp_species
      !         needed for proper allocation of pseudopotential array


      use mpi_tasks


      implicit none

      character*40 :: pp_path

      !local variables, counters
      integer :: n_func           !number of pseudowavefunction of certain species
      integer :: n_func_points    !number of gridpoints the pseudopotential are defined on
      integer :: i, i_code
      real*8 :: dummy
      character(*), parameter :: func = 'parse_pseudopot_dimensions'

      open(146,file=trim(pp_path),status='OLD',iostat=i_code)

      if (i_code/=0) then
         call aims_stop("* AIMS is missing a *.cpi file containing pseudopotential data.", func)
      endif

      read(146,*) dummy, n_func

      do i = 1,10
         read(146,*) dummy
      end do

      read(146,*)  n_func_points, dummy


      close(146)

      end subroutine parse_pseudopot_dimensions


!******
!-----------------------------------------------------------------------------------

end module dimensions
