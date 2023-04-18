!****h* FHI-aims/species_data
!  NAME
!    species_data
!  SYNOPSIS
 module species_data
!  PURPOSE
!     This module contains all data structures & (de)allocation routines relevant for the different 
!     species
!  USES
      implicit none
!
!
!     Subroutines:

!     * allocate_species_data
!     * allocate_species_waves
!     * cleanup_species_data
!     * cleanup_species_waves
!
!     Atomic variables:

!     * species_name : List of all element names in the project.
!     * species_element : List of actual chemical identities of the species.
!     * species_z : atomic numbers for each species in the game.
!     * atoms_in_structure : The number of occurrences of a given species i_species in the
!                 structure
! SVL * atoms_in_supercell : The number of occurrences of a given species i_species in the
!                 supercell defined by centers_ele_summation (all atoms that touch
!                 original unit cell according to atom_radius)
!     * basis_acc : a basis state |a> is included if
!                       [1.-sum(<b_i|a>)] > basis_acc
!                 for all previous basis states |b_i> of the same atom and angular momentum
!     * prodbas_acc : Similar to basis_acc, but for auxiliary basis
!     * innermost_max : Program stops if the innermost maximum of a radial
!                 function lies within the inner_max'th radial integration
!                 grid shell
!     * r_cutoff : for each species, onset of cutoff potential in basis defining
!                potential
!     * w_cutoff : for each species, width of cutoff potential in basis defining
!                potential

!     * scale_cutoff : Factor to scale the cutoff potential up and down, to allow more or less
!                rapid decay of basis functions. In fact, we scale (r-r_cut)*pot_scale, so
!                scale_cutoff can be directly related to the presumed width of the wave function tail.
!     * cut_free_atom : Flag that determines whether or not cutoff potential will be applied to the free-atom density.
!     * free_r_cut : If applicable, onset radius of cutoff potential for the free atom.

!     * l_shell_max   : maximum angular momentum shell for atomic-like basis functions for each species.
!     * l_hartree : maximum angular momentum component used for partitioned Hartree potential around species

!     * multipole_radius : Cutoff radius outside of which the numerically tabulated Hartree multipole potential
!               is replaced by an analytical field of fixed constant multipole moment
!               this radius is defined by the extent of charge partitioned to that atom.
!     * multipole_radius_free : Cutoff radius outside of which the free-atom charge is zero.

!     * valence_n_max : maximum occupied quantum number for given element and
!                     angular momentum
!     * valence_occ   : occupation of (atomic!) shells for free atom of each species
!     * ion_n_max : maximum shell number for given angular momentum in ion
!     * ion_occ  : Occupation of ionic shell (l,n). This specifies the charge state, and also
!                arbitrary excited states could be included if desired.
!     * include_min_basis : Flag which determines whether the "minimal" basis functions
!                 (i.e. the numeric free-atom basis functions) are part of the full basis set used.
!                 This flag is mainly a test flag to compare with Gaussian basis function
!                 implementations but may have other uses.
!     * n_atomic : Number of occupied free-atom basis functions per species
!     * atomic_n : radial quantum numbers of occupied free-atom basis functions
!     * atomic_l : angular momentum quantum numbers of occupied free-atom basis functions
!     * atomic_k : kappa quantum number for relativistic basis set
!     * n_ionic  : Number of ionic radial basis functions per species
!     * ionic_n  : radial quantum number of ionic function i_species, i_ionic
!     * ionic_l  : angular momentum quantum number of ionic function i_species, i_ionic
!     * ionic_rad: confinement radius of ionic function i_species, i_ionic
!     * n_conf   : number of empty/confined radial basis functions per species
!     * conf_n   : radial quantum number of confined function i_species, i_conf
!     * conf_l   : angular momentum quantum number of confined function i_species, i_conf
!     * conf_rad : confinement radius of confined function i_species, i_conf
!     * n_hydro  : Number of hydrogenic radial functions per species
!     * hydro_n  : radial quantum number of hydrogenic function i_species, i_hydro
!     * hydro_l  : angular momentum quantum number of hydrogenic function i_species, i_hydro
!     * hydro_scale : Scaling parameter (effective z) of hydrogenic function i_species, i_hydro
!     * n_gaussian  : Number of Gaussian radial functions per species
!     * gaussian_l  : cartesian quantum number of Gaussian radial function i_species, i_gaussian
!     * gaussian_n  : Number of elementary Gaussians which form Gaussian radial fn. # i_gaussian
!     * gaussian_alpha : exponent of Gaussian function i_species, i_gaussian
!     * gaussian_coeff : weight factor for elementary Gaussian function with exponent alpha
!     * max_n_prodbas  : The maximal radial quantum number of the single basis fns used to
!                       construct the product basis
!     * max_l_prodbas  : The maximal angular quantum number of the product basis fns.
!     * n_aux_gaussian  : Number of auxiliary Gaussian radial functions per species
!     * aux_gaussian_l  : cartesian quantum number of auxiliary Gaussian radial function i_species, i_gaussian
!     * aux_gaussian_n  : Number of elementary Gaussians which form auxiliary Gaussian radial fn. # i_gaussian. Here we use pure Gaussians so that aux_g
!     * aux_gaussian_alpha : exponent of auxiliary Gaussian function i_species, i_gaussian
!     * aux_gaussian_coeff : weight factor for elementary Gaussian function with exponent alpha. This is not needed for the moment.
!
!     * plus_u_flag : true if LDA(GGA)+U treatment requested for given species  
!     * plus_u_n : n shell to which Hubbard U is applied if LDA(GGA)+U is requested
!     * plus_u_l_str : string expression for angular momentum shell l 
!     * plus_u_l : l shell to which Hubbard U is applied if LDA(GGA)+U is requested
!     * plus_u_value : value of Hubbard U in eV if LDA(GGA)+U is requested
!     * plus_u_hubc : linear coeffietens for hubbard projectors 
!
!     Basis function related data:

!     * basis_potential : For each species, the basis-defining potential for
!                    atomic-like basis states
!     * atomic_wave  : tabulated minimal basis wave functions (with cutoff potential)
!                    for all species
!     * atomic_wave_deriv : first derivative of minimal basis wave function
!     * atm_wave_large : tabulated minimal large-component basis wave functions (with cutoff potential)
!                    for all species for relativistic cases (4c DKS and X2C)
!     * atm_wave_small : tabulated minimal small-component basis wave functions (with cutoff potential)
!                    for all species for relativistic cases (4c DKS and X2C)
!     * atomic_kinetic : non-rel. kinetic energy terms associated with atomic waves
!     * atomic_eigenval : tabulated eigenvalues (with cutoff potential) of all minimal
!                    basis wave functions
!     * ionic_pot : for each species and if requested, self-consistent potential of free ion
!     * ionic_wave : self-consistent wave functions of free ion(s)
!     * ionic_wave_deriv : first derivative of ionic basis radial fn.
!     * ionic_kinetic : non-rel. kinetic energy terms associated with ionic waves
!     * ionic_eigenval : self-consistent eigenvalues of free ion(s)
!     * confined_pot : for each species and if requested, confined atomic potential
!     * confined_wave : confined wave functions of free atom(s)
!     * confined_kinetic : non-rel. kinetic energy terms associated with confined waves
!     * confined_eigenval : confined eigenvalues of free atom(s)
!     * confined_wave_deriv : first derivative of confined basis radial fn.
!     * hydro_wave : hydrogenic radial functions
!     * hydro_wave_deriv : first derivative of hydrogenic radial functions
!     * hydro_kinetic : kinetic energy terms associated with hydrogenic waves
!     * gaussian_wave : cartesian Gaussian style radial functions
!           Format : u(r) = r^(L+1) * exp(-alpha*r^2)
!           where L is the cartesian Gaussian value L.
!           (The definition of cartesian Gaussians can be found e.g. in
!           Koch, Holthausen: "A Chemist's Guide to Density Functional Theory" p. 98.
!           We use here spherical harmonic basis functions 1/r u_l(r) Y_lm(theta,phi)
!           instead, i.e. each cartesian Gaussian fn L gives rise to potentially
!           multiple angular momenum basis functions for l = L, L-2, L-4, ... )
!     * gaussian_wave_deriv : first derivative of cartesian Gaussian style radial functions
!     * gaussian_kinetic : kinetic energy terms associated with Gaussian radial functions
!           Note that gaussian_kinetic is not the equivalent of the other kinetic energy terms -
!           it is only the second derivative of the radial function u(r), but does not include
!           the angular momentum barrier l(l+1)/r^2 - that is added only in shrink_fixed_basis.
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
      character*20, dimension(:), allocatable :: species_name
      character(len=2), dimension(:), allocatable :: species_element
      real*8, dimension(:),     allocatable :: species_z
      real*8,  dimension(:),     allocatable :: species_m
      integer, dimension(:),     allocatable :: atoms_in_structure
      integer, dimension(:),     allocatable :: atoms_in_supercell
      real*8,  dimension(:),     allocatable :: basis_acc
      real*8,  dimension(:),     allocatable :: prodbas_acc
      integer, dimension(:),     allocatable :: innermost_max
      real*8,  dimension(:),     allocatable :: r_cutoff
      real*8,  dimension(:),     allocatable :: w_cutoff
      real*8,  dimension(:),     allocatable :: scale_cutoff
      integer, dimension(:),     allocatable :: cutoff_type
      logical, dimension(:),     allocatable :: cut_free_atom
      real*8,  dimension(:),     allocatable :: free_r_cut
      real*8,  dimension(:),     allocatable :: basis_dep_cutoff_thresh
      logical, dimension(:),     allocatable :: cut_atomic_basis_functions
      logical, dimension(:),     allocatable :: cut_core
      real*8,  dimension(:),     allocatable :: core_r_cut
      integer, dimension(:),     allocatable :: l_shell_max
      integer, dimension(:),     allocatable :: ext_l_shell_max
      integer, dimension(:),     allocatable :: l_hartree
      real*8,  dimension(:),     allocatable :: multipole_radius_free

      real*8,  dimension(:),     allocatable :: multipole_radius_free_sq

      ! outer_partition_radius :
      ! This marks the farthest distance from a given atom where it should still
      ! influence the partitioning function ("partition_tab") used to partition
      ! integrands and the elecrton density.
      ! At the time of this writing, only needed in the construction of the 
      ! stratmann_smooth and stratmann_smoother partition tables.
      ! It may however be that this value will also be needed as a rigorous bound
      ! of the extent of local Hartree potential components elsewhere, later.
      ! Now that is a tricky question .... the present bound is multipole_radius_free,
      ! which may be slightly too small.
      real*8, dimension(:), allocatable :: outer_partition_radius

      integer, dimension(:,:),   allocatable :: valence_n_max
      real*8,  dimension(:,:,:), allocatable :: valence_occ
      integer, dimension(:,:),   allocatable :: core_n_max
      logical, dimension(:),     allocatable :: include_min_basis
      integer, dimension(:),     allocatable :: n_atomic
      integer, dimension(:,:),   allocatable :: atomic_n
      integer, dimension(:,:),   allocatable :: atomic_l
      integer, dimension(:,:),   allocatable :: atomic_k
      logical, dimension(:,:),   allocatable :: core_fn
      real*8,  dimension(:,:),   allocatable :: atomic_cutoff
      integer, dimension(:,:),   allocatable :: ion_n_max
      real*8,  dimension(:,:,:), allocatable :: ion_occ
      integer, dimension(:),     allocatable :: n_ionic
      integer, dimension(:,:),   allocatable :: ionic_n
      integer, dimension(:,:),   allocatable :: ionic_l
      integer, dimension(:,:),   allocatable :: ionic_kappa
      real*8,  dimension(:,:),   allocatable :: ionic_rad
      real*8,  dimension(:,:),   allocatable :: ionic_cutoff      
      logical, dimension(:,:),   allocatable :: ionic_in_large_basis
      integer, dimension(:),     allocatable :: n_conf
      integer, dimension(:,:),   allocatable :: conf_n
      integer, dimension(:,:),   allocatable :: conf_l
      integer, dimension(:,:),   allocatable :: conf_kappa
      real*8,  dimension(:,:),   allocatable :: conf_rad
      real*8,  dimension(:,:),   allocatable :: conf_cutoff
      integer, dimension(:),     allocatable :: n_hydro    
      integer, dimension(:,:),   allocatable :: hydro_n
      integer, dimension(:,:),   allocatable :: hydro_l
      integer, dimension(:,:),   allocatable :: hydro_kappa
      real*8,  dimension(:,:),   allocatable :: hydro_scale
      real*8,  dimension(:,:),   allocatable :: hydro_cutoff
      logical, dimension(:,:),   allocatable :: hydro_in_large_basis
      integer, dimension(:),     allocatable :: n_gaussian
      integer, dimension(:,:),   allocatable :: gaussian_n
      integer, dimension(:,:),   allocatable :: gaussian_l
      integer, dimension(:,:),   allocatable :: gaussian_n_contr
      real*8,  dimension(:,:,:), allocatable :: gaussian_alpha
      real*8,  dimension(:,:,:), allocatable :: gaussian_coeff
      logical, dimension(:),     allocatable :: pure_gaussian
      integer, allocatable :: n_sto(:)
      integer, allocatable :: sto_n(:,:)
      integer, allocatable :: sto_l(:,:)
      integer, allocatable :: sto_k(:,:)
      real*8, allocatable :: sto_zeta(:,:)
      integer, dimension(:),     allocatable :: max_n_prodbas
      integer, dimension(:),     allocatable :: max_l_prodbas
      logical, dimension(:),     allocatable :: specified_grid
      integer, dimension(:),     allocatable :: n_ang_shells
      real*8,  dimension(:,:),   allocatable :: r_ang_shell
      integer, dimension(:,:),   allocatable :: n_ang_points
      integer, dimension(:),     allocatable :: n_aux_gaussian
      integer, dimension(:,:),   allocatable :: aux_gaussian_n
      integer, dimension(:,:),   allocatable :: aux_gaussian_l
      integer, dimension(:,:),   allocatable :: aux_gaussian_n_contr
      real*8,  dimension(:,:,:), allocatable :: aux_gaussian_alpha
      real*8,  dimension(:,:,:), allocatable :: aux_gaussian_coeff
      real*8,  dimension(:,:),   allocatable :: basis_potential
      real*8,  dimension(:,:,:), allocatable :: atomic_wave
      real*8,  dimension(:,:,:), allocatable :: atomic_wave_deriv
      real*8,  dimension(:,:,:), allocatable :: atomic_large_deriv
      real*8,  dimension(:,:,:), allocatable :: atomic_small_deriv
      real*8,  dimension(:,:,:), allocatable :: atm_wave_large
      real*8,  dimension(:,:,:), allocatable :: atm_wave_small
      real*8,  dimension(:,:,:), allocatable :: atomic_kinetic
      real*8,  dimension(:,:,:), allocatable :: atomic_kinetic_small
      real*8,  dimension(:,:),   allocatable :: atomic_eigenval
      real*8,  dimension(:,:),   allocatable :: atomic_outer_radius
      real*8,  dimension(:,:,:), allocatable :: ionic_pot
      real*8,  dimension(:,:,:), allocatable :: ionic_wave
      real*8,  dimension(:,:,:), allocatable :: ionic_wave_large
      real*8,  dimension(:,:,:), allocatable :: ionic_wave_small
      real*8,  dimension(:,:,:), allocatable :: ionic_wave_deriv
      real*8,  dimension(:,:,:), allocatable :: ionic_large_deriv
      real*8,  dimension(:,:,:), allocatable :: ionic_small_deriv
      real*8,  dimension(:,:,:), allocatable :: ionic_kinetic
      real*8,  dimension(:,:,:), allocatable :: ionic_kinetic_small
      real*8,  dimension(:,:),   allocatable :: ionic_eigenval
      real*8,  dimension(:,:),   allocatable :: ionic_outer_radius
      real*8,  dimension(:,:,:), allocatable :: confined_pot
      real*8,  dimension(:,:,:), allocatable :: confined_wave
      real*8,  dimension(:,:,:), allocatable :: confined_wave_large
      real*8,  dimension(:,:,:), allocatable :: confined_wave_small
      real*8,  dimension(:,:,:), allocatable :: confined_wave_deriv
      real*8,  dimension(:,:,:), allocatable :: confined_kinetic
      real*8,  dimension(:,:),   allocatable :: confined_eigenval
      real*8,  dimension(:,:),   allocatable :: confined_outer_radius
      real*8,  dimension(:,:,:), allocatable :: hydro_wave
      real*8,  dimension(:,:,:), allocatable :: hydro_wave_deriv
      real*8,  dimension(:,:,:), allocatable :: hydro_large_deriv
      real*8,  dimension(:,:,:), allocatable :: hydro_small_deriv
      real*8,  dimension(:,:,:), allocatable :: hydro_kinetic
      real*8,  dimension(:,:,:), allocatable :: hydro_kinetic_small
      real*8,  dimension(:,:),   allocatable :: hydro_outer_radius
      real*8,  dimension(:,:,:), allocatable :: hydro_wave_large
      real*8,  dimension(:,:,:), allocatable :: hydro_wave_small
      real*8,  dimension(:,:,:), allocatable :: gaussian_wave
      real*8,  dimension(:,:,:), allocatable :: gaussian_wave_deriv
      real*8,  dimension(:,:,:), allocatable :: gaussian_kinetic
      real*8,  dimension(:,:),   allocatable :: gaussian_outer_radius
      real*8, allocatable :: sto_wave(:,:,:)
      real*8, allocatable :: sto_wave_small(:,:,:)
      real*8, allocatable :: sto_wave_deriv(:,:,:)
      real*8, allocatable :: sto_kinetic(:,:,:)
      real*8, allocatable :: sto_kinetic_small(:,:,:)
      real*8,  dimension(:,:,:), allocatable :: aux_gaussian_wave
      integer, dimension(:),     allocatable :: n_core_states_species
      integer, dimension(:),     allocatable :: n_core_basis_species

      integer, dimension(:), allocatable :: n_KH_core_states_species
      integer, dimension(:), allocatable :: n_KH_core_basis_species

      logical, dimension(:),     allocatable :: plus_u_flag
      integer, dimension(:,:),   allocatable :: plus_u_n
      character, dimension(:,:), allocatable :: plus_u_l_str
      integer, dimension(:,:),   allocatable :: plus_u_l
      real*8,  dimension(:,:),   allocatable :: plus_u_value
      real*8,  dimension(:,:),   allocatable :: plus_u_hubc
      real*8,  dimension(:),     allocatable :: plus_u_ramping_increment

      logical,dimension(:), allocatable :: species_pseudoized         !array which saves whether species is a pseudospecies
      logical, dimension(:), allocatable :: no_basis   ! is set true if species carries neither min basis nor any basis function

      ! convenice function to convert between nucleus charges and element names
      interface periodic_table
          module procedure get_nucleus
          module procedure get_element
      end interface

      contains
!******

!****s* species_data/allocate_species_data
!  NAME
!    allocate_species_data
!  SYNOPSIS
        subroutine allocate_species_data( )
!  PURPOSE
!    allocation of species data, after initial determination of the necessary dimensions
!
!    AJL: If you add a new allocation *PLEASE PLEASE* include a
!    deallocation in the cleanup routine, otherwise problems are
!    encountered when treating aims as a library
!  USES
        use dimensions, only : l_wave_max, n_max_ang_shells, n_max_aux_contracted, &
                               n_max_aux_fns, n_max_contracted, n_max_ind_fns, &
                               n_max_shells_plus_u, n_species, n_wave_max, &
                               use_aux_gaussian, use_confined, &
                               use_gaussian, use_hydro, use_ionic, use_plus_u, use_prodbas, &
                               use_relativistic_basis, use_specified_grid, use_sto

        implicit none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!  local variables

        integer :: i_shell_plus_u, i_species 
!  begin work


        allocate (species_name(n_species))
        allocate (species_element(n_species))
        allocate (species_z(n_species))
        allocate (species_m(n_species))
        allocate (atoms_in_structure(n_species))
        allocate (atoms_in_supercell(n_species)) !SVL
        allocate (basis_acc(n_species))
        allocate (prodbas_acc(n_species))
        allocate (innermost_max(n_species))


        allocate (r_cutoff(n_species))
        allocate (w_cutoff(n_species))
        allocate (scale_cutoff(n_species))
        allocate (cutoff_type(n_species))


        allocate(cut_free_atom(n_species))
        allocate(free_r_cut(n_species))
        allocate(basis_dep_cutoff_thresh(n_species))
        allocate(cut_atomic_basis_functions(n_species))


        allocate(cut_core(n_species))
        allocate(core_r_cut(n_species))

        allocate (l_shell_max(n_species))


        allocate (ext_l_shell_max(n_species))
        allocate (l_hartree(n_species))
!        allocate (multipole_radius(n_species))
        allocate (multipole_radius_free(n_species))
        allocate (multipole_radius_free_sq(n_species))


        allocate (outer_partition_radius(n_species))
        outer_partition_radius = 0.d0  ! initialize; zero is meaningless in practice

        allocate (valence_n_max(0:l_wave_max,n_species))
        allocate (valence_occ(n_wave_max,0:l_wave_max,n_species))

        allocate (core_n_max(0:l_wave_max,n_species))


        allocate (include_min_basis(n_species))
        allocate (n_atomic(n_species))
        allocate (atomic_n(n_species,n_max_ind_fns))
        allocate (atomic_l(n_species,n_max_ind_fns))
        if (use_relativistic_basis) then
            allocate (atomic_k(n_species,n_max_ind_fns))
        end if
        allocate (core_fn(n_species,n_max_ind_fns))
        allocate (atomic_cutoff(n_species,n_max_ind_fns))


        allocate(n_ionic(n_species))
        if (use_ionic) then
          allocate(ionic_in_large_basis(n_species, n_max_ind_fns))
          allocate(ion_n_max(0:l_wave_max,n_species))
          allocate(ion_occ(n_wave_max,0:l_wave_max,n_species))
          allocate(ionic_n(n_species,n_max_ind_fns))
          allocate(ionic_l(n_species,n_max_ind_fns))
          allocate(ionic_kappa(n_species,n_max_ind_fns))
          allocate(ionic_rad(n_species,n_max_ind_fns))
          allocate(ionic_cutoff(n_species,n_max_ind_fns))
        end if


        allocate(n_conf(n_species))
        if (use_confined) then
          allocate(conf_n(n_species,n_max_ind_fns))
          allocate(conf_l(n_species,n_max_ind_fns))
          allocate(conf_kappa(n_species,n_max_ind_fns))
          allocate(conf_rad(n_species,n_max_ind_fns))
          allocate(conf_cutoff(n_species,n_max_ind_fns))
        end if


        allocate(n_hydro(n_species))
        if (use_hydro) then
          allocate(hydro_in_large_basis(n_species, n_max_ind_fns))
          allocate(hydro_n(n_species,n_max_ind_fns))
          allocate(hydro_l(n_species,n_max_ind_fns))
          allocate(hydro_kappa(n_species,n_max_ind_fns))
          allocate(hydro_scale(n_species,n_max_ind_fns))
          allocate(hydro_cutoff(n_species,n_max_ind_fns))
        end if


        allocate(n_gaussian(n_species))
        if (use_gaussian) then
          allocate(gaussian_n(n_species,n_max_ind_fns))
          allocate(gaussian_l(n_species,n_max_ind_fns))
          allocate(gaussian_n_contr(n_species,n_max_ind_fns))
          allocate( &
            gaussian_alpha( n_species,n_max_ind_fns, n_max_contracted) )
          allocate( &
            gaussian_coeff( n_species,n_max_ind_fns, n_max_contracted) )
          allocate (pure_gaussian(n_species))
          pure_gaussian(:) = .true.
        end if


        allocate(n_aux_gaussian(n_species))
        if (use_aux_gaussian) then
          allocate(aux_gaussian_n(n_species,n_max_aux_fns))
          allocate(aux_gaussian_l(n_species,n_max_aux_fns))
          allocate(aux_gaussian_n_contr(n_species,n_max_aux_fns))
          allocate(aux_gaussian_alpha &
                   ( n_species,n_max_aux_fns, n_max_aux_contracted) )
          allocate(aux_gaussian_coeff &
                   ( n_species,n_max_aux_fns, n_max_aux_contracted) )
        end if

        allocate(n_sto(n_species))
        n_sto = 0
        if (use_sto) then
           allocate(sto_n(n_species,n_max_ind_fns))
           allocate(sto_l(n_species,n_max_ind_fns))
           allocate(sto_k(n_species,n_max_ind_fns))
           allocate(sto_zeta(n_species,n_max_ind_fns))
        end if

        if (use_prodbas) then
          allocate(max_n_prodbas(n_species))
          allocate(max_l_prodbas(n_species))
        end if


        allocate (specified_grid(n_species))
        specified_grid = .false.

        if (use_specified_grid) then
          allocate(n_ang_shells(n_species))
          allocate(r_ang_shell(n_max_ang_shells,n_species))
          allocate(n_ang_points(n_max_ang_shells,n_species))
        end if


        allocate(n_core_states_species(n_species))
        n_core_states_species = 0
        allocate(n_core_basis_species(n_species))

        allocate(n_KH_core_states_species(n_species))
        n_KH_core_states_species = 0
        allocate(n_KH_core_basis_species(n_species))


        if(use_plus_u) then
           allocate(plus_u_flag(n_species))
           allocate(plus_u_n(n_species, n_max_shells_plus_u))
           allocate(plus_u_l_str(n_species, n_max_shells_plus_u))
           allocate(plus_u_l(n_species, n_max_shells_plus_u))
           allocate(plus_u_value(n_species, n_max_shells_plus_u))
           allocate(plus_u_ramping_increment(n_species))
           ! at the moment only linear combination of 4 basis functions are
           ! allowed
           allocate(plus_u_hubc(n_species,4))
 
          ! Initialize all array elements to 0
           do i_species = 1, n_species, 1
              plus_u_flag(i_species) = .false.
              do i_shell_plus_u = 1, n_max_shells_plus_u, 1
                 plus_u_n(i_species, i_shell_plus_u) = 0
                 plus_u_l_str(i_species, i_shell_plus_u) = ""
                 plus_u_l(i_species, i_shell_plus_u) = 0
                 plus_u_value(i_species, i_shell_plus_u) = 0.0d0
                 plus_u_hubc(i_species,:) = 1.0d0
                 plus_u_ramping_increment(i_species) =0.0d0
              end do
           end do
        end if


        allocate(species_pseudoized(n_species))
        species_pseudoized = .false.

        allocate(no_basis(n_species))  
        no_basis = .false.


        end subroutine allocate_species_data
!******

!****s* species_data/allocate_species_waves
!  NAME
!    allocate_species_waves
!  SYNOPSIS
        subroutine allocate_species_waves( )
!  PURPOSE
!    allocation of species wave function data
!  USES
        use dimensions, only : use_min_basis, use_basis_gradients, &
                               use_ionic, use_confined, use_hydro, &
                               use_gaussian, &
                               use_aux_gaussian, n_max_aux_fns, &
                               n_max_grid, n_max_ind_fns, n_species, use_sto
        use runtime_choices, only: flag_rel, REL_x2c, REL_4c_dks
        implicit none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
        allocate(basis_potential(n_max_grid, n_species))
        if (use_min_basis) then

          if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
            allocate(atm_wave_large(n_max_grid,n_species,n_max_ind_fns))
            allocate(atm_wave_small(n_max_grid,n_species,n_max_ind_fns))
            allocate(atomic_large_deriv(n_max_grid,n_species,n_max_ind_fns))
            allocate(atomic_small_deriv(n_max_grid,n_species,n_max_ind_fns))
           ! For fully-relativistic case, atomic_kinetic saves the kinetic
           ! basis for large comp. This is for small comp.:
            allocate(atomic_kinetic_small(n_max_grid,n_species,n_max_ind_fns))
          else
            allocate(atomic_wave(n_max_grid,n_species,n_max_ind_fns))
          endif
          allocate(atomic_kinetic(n_max_grid,n_species,n_max_ind_fns))
          if (use_basis_gradients) then
             allocate(atomic_wave_deriv(n_max_grid,n_species,n_max_ind_fns))
          end if
          allocate(atomic_eigenval(n_species,n_max_ind_fns))

        end if
        allocate(atomic_outer_radius(n_species,n_max_ind_fns))

        if (use_ionic) then
          if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
            allocate(ionic_wave_large(n_max_grid,n_species,n_max_ind_fns))
            allocate(ionic_wave_small(n_max_grid,n_species,n_max_ind_fns))
            allocate(ionic_large_deriv(n_max_grid,n_species,n_max_ind_fns))
            allocate(ionic_small_deriv(n_max_grid,n_species,n_max_ind_fns))
            allocate(ionic_kinetic_small(n_max_grid,n_species,n_max_ind_fns))
          else
            allocate(ionic_wave(n_max_grid,n_species,n_max_ind_fns))
          endif
          allocate(ionic_pot(n_max_grid,n_species,n_max_ind_fns))
          allocate(ionic_kinetic(n_max_grid,n_species,n_max_ind_fns))
          allocate(ionic_outer_radius(n_species,n_max_ind_fns))

          if (use_basis_gradients) then
             allocate(ionic_wave_deriv(n_max_grid,n_species,n_max_ind_fns))
          end if
          allocate(ionic_eigenval(n_species,n_max_ind_fns))
        end if

        if (use_confined) then
          if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
            allocate(confined_wave_large(n_max_grid,n_species,n_max_ind_fns))
            allocate(confined_wave_small(n_max_grid,n_species,n_max_ind_fns))
          else
            allocate(confined_wave(n_max_grid,n_species,n_max_ind_fns))
          endif
          allocate(confined_pot(n_max_grid,n_species,n_max_ind_fns))
          allocate(confined_kinetic(n_max_grid,n_species,n_max_ind_fns))
          allocate(confined_outer_radius(n_species,n_max_ind_fns))

          if (use_basis_gradients) then
            allocate &
            ( confined_wave_deriv (n_max_grid,n_species,n_max_ind_fns) )
          end if
          allocate(confined_eigenval(n_species,n_max_ind_fns))
        end if

        if (use_hydro) then
          if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
            allocate(hydro_wave_large(n_max_grid,n_species,n_max_ind_fns))
            allocate(hydro_wave_small(n_max_grid,n_species,n_max_ind_fns))
            allocate(hydro_large_deriv(n_max_grid,n_species,n_max_ind_fns))
            allocate(hydro_small_deriv(n_max_grid,n_species,n_max_ind_fns))
           ! For fully-relativistic case, hydro_kinetic saves the kinetic
           ! basis for large comp. This is for small comp.:
            allocate(hydro_kinetic_small(n_max_grid,n_species,n_max_ind_fns))
          else
            allocate(hydro_wave(n_max_grid,n_species,n_max_ind_fns))
          endif
          if (use_basis_gradients) then
             allocate(hydro_wave_deriv(n_max_grid,n_species,n_max_ind_fns))
          end if
          
          allocate(hydro_kinetic(n_max_grid,n_species,n_max_ind_fns))
          allocate(hydro_outer_radius(n_species,n_max_ind_fns))
        end if


        if (use_gaussian) then
          allocate(gaussian_wave(n_max_grid,n_species,n_max_ind_fns))
          if (use_basis_gradients) then
            allocate &
            ( gaussian_wave_deriv (n_max_grid,n_species,n_max_ind_fns) )
          end if
          
          allocate(gaussian_kinetic(n_max_grid,n_species,n_max_ind_fns))
          allocate(gaussian_outer_radius(n_species,n_max_ind_fns))
        end if

        if (use_aux_gaussian) then
          allocate(aux_gaussian_wave &
                   (n_max_grid,n_species,n_max_aux_fns))
        end if

        if (use_sto) then
          allocate(sto_wave(n_max_grid,n_species,n_max_ind_fns))
          if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
            allocate(sto_wave_small(n_max_grid,n_species,n_max_ind_fns))
            allocate(sto_kinetic_small(n_max_grid,n_species,n_max_ind_fns))
          endif
          if (use_basis_gradients) &
                & allocate(sto_wave_deriv(n_max_grid,n_species,n_max_ind_fns))
          allocate(sto_kinetic(n_max_grid,n_species,n_max_ind_fns))
        end if

        end subroutine allocate_species_waves
!******

!****s* species_data/cleanup_species_data
!  NAME
!    cleanup_species_data
!  SYNOPSIS
        subroutine cleanup_species_data ( )
!  PURPOSE
!    deallocation of species data
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        if (allocated(species_name)) then
          deallocate (species_name)
        end if
        if (allocated(species_element)) then
          deallocate (species_element)
        end if
        if (allocated(species_z)) then
          deallocate (species_z)
        end if
        if (allocated(species_m)) then
          deallocate (species_m)
        end if
        if (allocated(atoms_in_structure)) then
          deallocate (atoms_in_structure)
        end if
! SVL
        if (allocated(atoms_in_supercell)) then
           deallocate (atoms_in_supercell)
        end if
        if (allocated(basis_acc)) then
          deallocate (basis_acc)
        end if
        if (allocated(prodbas_acc)) then
          deallocate (prodbas_acc)
        end if
        if (allocated(innermost_max)) then
          deallocate (innermost_max)
        end if
        if (allocated(r_cutoff)) then
          deallocate (r_cutoff)
        end if
        if (allocated(w_cutoff)) then
          deallocate (w_cutoff)
        end if
        if (allocated(scale_cutoff)) then
          deallocate (scale_cutoff)
        end if
        if (allocated(cutoff_type)) then
          deallocate (cutoff_type)
        end if

        if (allocated(cut_free_atom)) then
          deallocate(cut_free_atom)
        end if
        if (allocated(free_r_cut)) then
          deallocate(free_r_cut)
        end if

        if (allocated(basis_dep_cutoff_thresh)) then
           deallocate(basis_dep_cutoff_thresh)
        end if
        if (allocated(cut_atomic_basis_functions)) then
           deallocate(cut_atomic_basis_functions)
        end if

        if (allocated(cut_core)) then
          deallocate(cut_core)
        end if
        if (allocated(core_r_cut)) then
          deallocate(core_r_cut)
        end if

        if (allocated(l_shell_max)) then
          deallocate (l_shell_max)
        end if
        if (allocated(ext_l_shell_max)) then
          deallocate (ext_l_shell_max)
        end if
        
        if (allocated(l_hartree)) then
          deallocate (l_hartree)
        end if
        if (allocated(multipole_radius_free)) then
          deallocate (multipole_radius_free)
        end if
        if (allocated(multipole_radius_free_sq)) then
          deallocate (multipole_radius_free_sq)
        end if
        if (allocated(outer_partition_radius)) then
          deallocate (outer_partition_radius)
        end if

        if (allocated(valence_n_max)) then
          deallocate (valence_n_max)
        end if
        if (allocated(valence_occ)) then
          deallocate (valence_occ)
        end if

        if (allocated(core_n_max)) then
          deallocate (core_n_max)
        end if

        if (allocated(include_min_basis)) then
          deallocate (include_min_basis)
        end if
        if (allocated(n_atomic)) then
          deallocate (n_atomic)
        end if
        if (allocated(atomic_n)) then
          deallocate (atomic_n)
        end if
        if (allocated(atomic_l)) then
          deallocate (atomic_l)
        end if
        if (allocated(atomic_k)) then
          deallocate (atomic_k)
        end if
        if (allocated(core_fn)) then
          deallocate (core_fn)
        end if
        if (allocated(atomic_cutoff)) then
           deallocate(atomic_cutoff)
        end if

        if (allocated(ionic_in_large_basis)) then
            deallocate(ionic_in_large_basis)
        endif
        if (allocated(ion_n_max)) then
            deallocate(ion_n_max)
        end if
        if (allocated(ion_occ)) then
            deallocate(ion_occ)
        end if
        if (allocated(n_ionic)) then
            deallocate(n_ionic)
        end if
        if (allocated(ionic_n)) then
            deallocate(ionic_n)
        end if
        if (allocated(ionic_l)) then
            deallocate(ionic_l)
        end if
        if (allocated(ionic_kappa)) then
            deallocate(ionic_kappa)
        end if
        if (allocated(ionic_rad)) then
            deallocate(ionic_rad)
        end if
        if (allocated(ionic_cutoff)) then
           deallocate(ionic_cutoff)
        end if

        if (allocated(n_conf)) then
            deallocate(n_conf)
        end if
        if (allocated(conf_n)) then
            deallocate(conf_n)
        end if
        if (allocated(conf_l)) then
            deallocate(conf_l)
        end if
        if (allocated(conf_kappa)) then
            deallocate(conf_kappa)
        end if
        if (allocated(conf_rad)) then
            deallocate(conf_rad)
        end if
        if (allocated(conf_cutoff)) then
           deallocate(conf_cutoff)
        end if

        if (allocated(n_hydro)) then
            deallocate(n_hydro)
        end if
        if (allocated(hydro_in_large_basis)) then
            deallocate(hydro_in_large_basis)
        end if
        if (allocated(hydro_n)) then
            deallocate(hydro_n)
        end if
        if (allocated(hydro_l)) then
            deallocate(hydro_l)
        end if
        if (allocated(hydro_kappa)) then
            deallocate(hydro_kappa)
        end if
        if (allocated(hydro_scale)) then
            deallocate(hydro_scale)
        end if
        if (allocated(hydro_cutoff)) then
           deallocate(hydro_cutoff)
        end if

        if (allocated(n_gaussian)) then
            deallocate(n_gaussian)
        end if
        if (allocated(gaussian_n)) then
            deallocate(gaussian_n)
        end if
        if (allocated(gaussian_l)) then
            deallocate(gaussian_l)
        end if
        if (allocated(gaussian_n_contr)) then
            deallocate(gaussian_n_contr)
        end if
        if (allocated(gaussian_alpha)) then
            deallocate(gaussian_alpha)
        end if
        if (allocated(gaussian_coeff)) then
            deallocate(gaussian_coeff)
        end if
        if (allocated(pure_gaussian)) then
          deallocate (pure_gaussian)
        end if

        if (allocated(n_aux_gaussian)) then
            deallocate(n_aux_gaussian)
        end if
        if (allocated(aux_gaussian_n)) then
            deallocate(aux_gaussian_n)
        end if
        if (allocated(aux_gaussian_l)) then
            deallocate(aux_gaussian_l)
        end if
        if (allocated(aux_gaussian_n_contr)) then
            deallocate(aux_gaussian_n_contr)
        end if
        if (allocated(aux_gaussian_alpha)) then
            deallocate(aux_gaussian_alpha)
        end if
        if (allocated(aux_gaussian_coeff)) then
            deallocate(aux_gaussian_coeff)
        end if

        if (allocated(n_sto)) deallocate(n_sto)
        if (allocated(sto_n)) deallocate(sto_n)
        if (allocated(sto_l)) deallocate(sto_l)
        if (allocated(sto_k)) deallocate(sto_k)
        if (allocated(sto_zeta)) deallocate(sto_zeta)

        if (allocated(max_n_prodbas)) then
          deallocate(max_n_prodbas)
        end if
        if (allocated(max_l_prodbas)) then
          deallocate(max_l_prodbas)
        end if

        if(allocated (specified_grid))then
           deallocate (specified_grid)
        end if
        if (allocated(n_ang_shells)) then
          deallocate(n_ang_shells)
        end if
        if (allocated(r_ang_shell)) then
          deallocate(r_ang_shell)
        end if
        if (allocated(n_ang_points)) then
          deallocate(n_ang_points)
        end if


        if(allocated(n_core_states_species))then
           deallocate(n_core_states_species)
        end if
        if(allocated(n_core_basis_species))then
           deallocate(n_core_basis_species)
        end if

        if(allocated(n_KH_core_states_species))then
           deallocate(n_KH_core_states_species)
        end if
        if(allocated(n_KH_core_basis_species))then
           deallocate(n_KH_core_basis_species)
        end if

        if(allocated(plus_u_flag)) then
           deallocate(plus_u_flag)
        end if
        if(allocated(plus_u_n)) then
           deallocate(plus_u_n)
        end if
        if(allocated(plus_u_l_str)) then
           deallocate(plus_u_l_str)
        end if
        if(allocated(plus_u_l)) then
           deallocate(plus_u_l)
        end if
        if(allocated(plus_u_value)) then
           deallocate(plus_u_value)
        end if
        if(allocated(plus_u_ramping_increment)) then
           deallocate(plus_u_ramping_increment)
        end if
        if(allocated(plus_u_hubc)) then
           deallocate(plus_u_hubc)
        endif
        if (allocated(species_pseudoized)) then
           deallocate(species_pseudoized)
        end if
        if (allocated(no_basis)) then
           deallocate(no_basis)
        end if
        end subroutine cleanup_species_data
!******

!****s* species_data/cleanup_species_waves
!  NAME
!    cleanup_species_waves
!  SYNOPSIS
        subroutine cleanup_species_waves( )
!  PURPOSE
!    deallocation of wave function data that are species-dependent
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        if (allocated(basis_potential)) then
          deallocate(basis_potential)
        end if

        if (allocated(atomic_wave)) then
          deallocate(atomic_wave)
        end if

        if (allocated(atm_wave_large)) then
          deallocate(atm_wave_large)
        end if

        if (allocated(atm_wave_small)) then
          deallocate(atm_wave_small)
        end if

        if (allocated(atomic_wave_deriv)) then
          deallocate(atomic_wave_deriv)
        end if
        if (allocated(atomic_large_deriv)) then
          deallocate(atomic_large_deriv)
        end if
        if (allocated(atomic_small_deriv)) then
          deallocate(atomic_small_deriv)
        end if
        if (allocated(atomic_kinetic)) then
          deallocate(atomic_kinetic)
        end if
        if (allocated(atomic_kinetic_small)) then
          deallocate(atomic_kinetic_small)
        end if
        if (allocated(atomic_eigenval)) then
          deallocate(atomic_eigenval)
        end if
        if (allocated(atomic_outer_radius)) then
          deallocate(atomic_outer_radius)
        end if

        if (allocated(ionic_pot)) then
          deallocate(ionic_pot)
        end if
        if (allocated(ionic_wave)) then
          deallocate(ionic_wave)
        end if
        if (allocated(ionic_wave_large)) then
          deallocate(ionic_wave_large)
        end if
        if (allocated(ionic_wave_small)) then
          deallocate(ionic_wave_small)
        end if
        if (allocated(ionic_kinetic)) then
          deallocate(ionic_kinetic)
        end if
        if (allocated(ionic_kinetic_small)) then
          deallocate(ionic_kinetic_small)
        end if
        if (allocated(ionic_outer_radius)) then
          deallocate(ionic_outer_radius)
        end if

        if (allocated(ionic_wave_deriv)) then
          deallocate(ionic_wave_deriv)
        end if
        if (allocated(ionic_large_deriv)) then
          deallocate(ionic_large_deriv)
        end if
        if (allocated(ionic_small_deriv)) then
          deallocate(ionic_small_deriv)
        end if
        if (allocated(ionic_eigenval)) then
          deallocate(ionic_eigenval)
        end if

        if (allocated(confined_pot)) then
          deallocate(confined_pot)
        end if
        if (allocated(confined_wave)) then
         deallocate(confined_wave)
        end if
        if (allocated(confined_wave_large)) then
         deallocate(confined_wave_large)
        end if
        if (allocated(confined_wave_small)) then
         deallocate(confined_wave_small)
        end if
        if (allocated(confined_kinetic)) then
         deallocate(confined_kinetic)
        end if
        if (allocated(confined_wave_deriv)) then
         deallocate(confined_wave_deriv)
        end if
        if (allocated(confined_eigenval)) then
          deallocate(confined_eigenval)
        end if
        if (allocated(confined_outer_radius)) then
          deallocate(confined_outer_radius)
        end if

        if (allocated(hydro_wave)) then
          deallocate(hydro_wave)
        end if
        if (allocated(hydro_wave_deriv)) then
          deallocate(hydro_wave_deriv)
        end if
        if (allocated(hydro_large_deriv)) then
          deallocate(hydro_large_deriv)
        end if
        if (allocated(hydro_small_deriv)) then
          deallocate(hydro_small_deriv)
        end if
        if (allocated(hydro_kinetic)) then
          deallocate(hydro_kinetic)
        end if
        if (allocated(hydro_kinetic_small)) then
          deallocate(hydro_kinetic_small)
        end if
        if (allocated(hydro_outer_radius)) then
          deallocate(hydro_outer_radius)
        end if
        if (allocated(hydro_wave_large)) then
          deallocate(hydro_wave_large)
        endif
        if (allocated(hydro_wave_small)) then
          deallocate(hydro_wave_small)
        endif

        if (allocated(gaussian_wave)) then
          deallocate(gaussian_wave)
        end if
        if (allocated(gaussian_wave_deriv)) then
          deallocate(gaussian_wave_deriv)
        end if
        if (allocated(gaussian_kinetic)) then
          deallocate(gaussian_kinetic)
        end if
        if (allocated(gaussian_outer_radius)) then
          deallocate(gaussian_outer_radius)
        end if

        if (allocated(aux_gaussian_wave)) then
          deallocate(aux_gaussian_wave)
        end if

        if (allocated(sto_wave)) then
          deallocate(sto_wave)
        end if
        if (allocated(sto_wave_small)) then
          deallocate(sto_wave_small)
        end if
        if (allocated(sto_wave_deriv)) then
          deallocate(sto_wave_deriv)
        end if
        if (allocated(sto_kinetic)) then
          deallocate(sto_kinetic)
        end if
        if (allocated(sto_kinetic_small)) then
          deallocate(sto_kinetic_small)
        end if

        end subroutine cleanup_species_waves
!******

 integer elemental function get_nucleus(element) result(nucleus)
    character(len=*), intent(in) :: element

    select case(lower(trim(element)))
    case ('h');  nucleus = 1;   case ('he'); nucleus = 2;   case ('li'); nucleus = 3
    case ('be'); nucleus = 4;   case ('b');  nucleus = 5;   case ('c');  nucleus = 6
    case ('n');  nucleus = 7;   case ('o');  nucleus = 8;   case ('f');  nucleus = 9
    case ('ne'); nucleus = 10;  case ('na'); nucleus = 11;  case ('mg'); nucleus = 12
    case ('al'); nucleus = 13;  case ('si'); nucleus = 14;  case ('p');  nucleus = 15
    case ('s');  nucleus = 16;  case ('cl'); nucleus = 17;  case ('ar'); nucleus = 18
    case ('k');  nucleus = 19;  case ('ca'); nucleus = 20;  case ('sc'); nucleus = 21
    case ('ti'); nucleus = 22;  case ('v');  nucleus = 23;  case ('cr'); nucleus = 24
    case ('mn'); nucleus = 25;  case ('fe'); nucleus = 26;  case ('co'); nucleus = 27
    case ('ni'); nucleus = 28;  case ('cu'); nucleus = 29;  case ('zn'); nucleus = 30
    case ('ga'); nucleus = 31;  case ('ge'); nucleus = 32;  case ('as'); nucleus = 33
    case ('se'); nucleus = 34;  case ('br'); nucleus = 35;  case ('kr'); nucleus = 36
    case ('rb'); nucleus = 37;  case ('sr'); nucleus = 38;  case ('y');  nucleus = 39
    case ('zr'); nucleus = 40;  case ('nb'); nucleus = 41;  case ('mo'); nucleus = 42
    case ('tc'); nucleus = 43;  case ('ru'); nucleus = 44;  case ('rh'); nucleus = 45
    case ('pd'); nucleus = 46;  case ('ag'); nucleus = 47;  case ('cd'); nucleus = 48
    case ('in'); nucleus = 49;  case ('sn'); nucleus = 50;  case ('sb'); nucleus = 51
    case ('te'); nucleus = 52;  case ('i');  nucleus = 53;  case ('xe'); nucleus = 54
    case ('cs'); nucleus = 55;  case ('ba'); nucleus = 56;  case ('la'); nucleus = 57
    case ('ce'); nucleus = 58;  case ('pr'); nucleus = 59;  case ('nd'); nucleus = 60
    case ('pm'); nucleus = 61;  case ('sm'); nucleus = 62;  case ('eu'); nucleus = 63
    case ('gd'); nucleus = 64;  case ('tb'); nucleus = 65;  case ('dy'); nucleus = 66
    case ('ho'); nucleus = 67;  case ('er'); nucleus = 68;  case ('tm'); nucleus = 69
    case ('yb'); nucleus = 70;  case ('lu'); nucleus = 71;  case ('hf'); nucleus = 72
    case ('ta'); nucleus = 73;  case ('w');  nucleus = 74;  case ('re'); nucleus = 75
    case ('os'); nucleus = 76;  case ('ir'); nucleus = 77;  case ('pt'); nucleus = 78
    case ('au'); nucleus = 79;  case ('hg'); nucleus = 80;  case ('tl'); nucleus = 81
    case ('pb'); nucleus = 82;  case ('bi'); nucleus = 83;  case ('po'); nucleus = 84
    case ('at'); nucleus = 85;  case ('rn'); nucleus = 86;  case ('fr'); nucleus = 87
    case ('ra'); nucleus = 88;  case ('ac'); nucleus = 89;  case ('th'); nucleus = 90
    case ('pa'); nucleus = 91;  case ('u');  nucleus = 92;  case ('np'); nucleus = 93
    case ('pu'); nucleus = 94;  case ('am'); nucleus = 95;  case ('cm'); nucleus = 96
    case ('bk'); nucleus = 97;  case ('cf'); nucleus = 98;  case ('es'); nucleus = 99
    case ('fm'); nucleus = 100; case ('md'); nucleus = 101; case ('no'); nucleus = 102
    case ('lr'); nucleus = 103; case ('rf'); nucleus = 104; case ('db'); nucleus = 105
    case ('sg'); nucleus = 106; case ('bh'); nucleus = 107; case ('hs'); nucleus = 108
    case ('mt'); nucleus = 109; case ('ds'); nucleus = 110; case ('rg'); nucleus = 111
    case ('cn'); nucleus = 112
    case default
        nucleus = -1
    end select

    contains

    elemental function lower(str)
        character(len=*), intent(in) :: str
        character(len=len(str)) :: lower

        integer :: i

        lower = '  '
        do i = 1, len(str)
            select case (str(i:i))
                case ('A':'Z')
                    lower(i:i) = achar(iachar(str(i:i))+32)
                case default
                    lower(i:i) = str(i:i)
            end select
        end do
    end function
 end function get_nucleus

 character(len=2) elemental function get_element(nucleus) result(elem)
    integer, intent(in) :: nucleus

    select case(nucleus)
    case (1);   elem = 'H';  case (2);   elem = 'He'; case (3);   elem = 'Li'
    case (4);   elem = 'Be'; case (5);   elem = 'B';  case (6);   elem = 'C'
    case (7);   elem = 'N';  case (8);   elem = 'O';  case (9);   elem = 'F'
    case (10);  elem = 'Ne'; case (11);  elem = 'Na'; case (12);  elem = 'Mg'
    case (13);  elem = 'Al'; case (14);  elem = 'Si'; case (15);  elem = 'P'
    case (16);  elem = 'S';  case (17);  elem = 'Cl'; case (18);  elem = 'Ar'
    case (19);  elem = 'K';  case (20);  elem = 'Ca'; case (21);  elem = 'Sc'
    case (22);  elem = 'Ti'; case (23);  elem = 'V';  case (24);  elem = 'Cr'
    case (25);  elem = 'Mn'; case (26);  elem = 'Fe'; case (27);  elem = 'Co'
    case (28);  elem = 'Ni'; case (29);  elem = 'Cu'; case (30);  elem = 'Zn'
    case (31);  elem = 'Ga'; case (32);  elem = 'Ge'; case (33);  elem = 'As'
    case (34);  elem = 'Se'; case (35);  elem = 'Br'; case (36);  elem = 'Kr'
    case (37);  elem = 'Rb'; case (38);  elem = 'Sr'; case (39);  elem = 'Y'
    case (40);  elem = 'Zr'; case (41);  elem = 'Nb'; case (42);  elem = 'Mo'
    case (43);  elem = 'Tc'; case (44);  elem = 'Ru'; case (45);  elem = 'Rh'
    case (46);  elem = 'Pd'; case (47);  elem = 'Ag'; case (48);  elem = 'Cd'
    case (49);  elem = 'In'; case (50);  elem = 'Sn'; case (51);  elem = 'Sb'
    case (52);  elem = 'Te'; case (53);  elem = 'I';  case (54);  elem = 'Xe'
    case (55);  elem = 'Cs'; case (56);  elem = 'Ba'; case (57);  elem = 'La'
    case (58);  elem = 'Ce'; case (59);  elem = 'Pr'; case (60);  elem = 'Nd'
    case (61);  elem = 'Pm'; case (62);  elem = 'Sm'; case (63);  elem = 'Eu'
    case (64);  elem = 'Gd'; case (65);  elem = 'Tb'; case (66);  elem = 'Dy'
    case (67);  elem = 'Ho'; case (68);  elem = 'Er'; case (69);  elem = 'Tm'
    case (70);  elem = 'Yb'; case (71);  elem = 'Lu'; case (72);  elem = 'Hf'
    case (73);  elem = 'Ta'; case (74);  elem = 'W';  case (75);  elem = 'Re'
    case (76);  elem = 'Os'; case (77);  elem = 'Ir'; case (78);  elem = 'Pt'
    case (79);  elem = 'Au'; case (80);  elem = 'Hg'; case (81);  elem = 'Tl'
    case (82);  elem = 'Pb'; case (83);  elem = 'Bi'; case (84);  elem = 'Po'
    case (85);  elem = 'At'; case (86);  elem = 'Rn'; case (87);  elem = 'Fr'
    case (88);  elem = 'Ra'; case (89);  elem = 'Ac'; case (90);  elem = 'Th'
    case (91);  elem = 'Pa'; case (92);  elem = 'U';  case (93);  elem = 'Np'
    case (94);  elem = 'Pu'; case (95);  elem = 'Am'; case (96);  elem = 'Cm'
    case (97);  elem = 'Bk'; case (98);  elem = 'Cf'; case (99);  elem = 'Es'
    case (100); elem = 'Fm'; case (101); elem = 'Md'; case (102); elem = 'No'
    case (103); elem = 'Lr'; case (104); elem = 'Rf'; case (105); elem = 'Db'
    case (106); elem = 'Sg'; case (107); elem = 'Bh'; case (108); elem = 'Hs'
    case (109); elem = 'Mt'; case (110); elem = 'Ds'; case (111); elem = 'Rg'
    case (112); elem = 'Cn'
    case default
        elem = 'X'
    end select
 end function get_element

 end module species_data
