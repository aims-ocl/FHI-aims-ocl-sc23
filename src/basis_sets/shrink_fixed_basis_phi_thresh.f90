!  FIXME: For optimum integrations, the present subroutine should also interpolate all
!         wave functions as splines onto the radial integration grid, and _then_ perform
!         all normalisation integrals on this grid. THEN, we would have a numerically
!         well-defined Hamiltonian matrix later on.
!         The present case will produce exact orthonormality on the log grid, but only
!         approximate orthonormality on the later radial integration grid.
!****s* FHI-aims/shrink_fixed_basis_phi_thresh
!  NAME
!   shrink_fixed_basis_phi_thresh
!  SYNOPSIS
subroutine shrink_fixed_basis_phi_thresh ( )
!  PURPOSE
!   Subroutine shrink_fixed_basis:
!   * orthogonalizes extra basis functions to atomic ones, for each species.
!   * throws out all basis functions which are substantially identical to another one
!   * gathers basis functions in a compressed array for numerics
!   * implicitly sort all basis functions into blocks of identical l per species
!   * sorts those blocks according to basis-dependent cutoff potential width 
!   * labels all basis functions to remember their origin
!  USES
  use runtime_choices, only : atomic_solver, ATOMIC_SOLVER_SRATOM, calculate_atom_bsse, &
                              flag_rel, REL_atomic_zora, REL_x2c, REL_4c_dks, &
                              wave_threshold, force_smooth_cutoff, n_empty, &
                              out_aitranss, out_basis, REL_own, smooth_cutoff_limit
  use dimensions,      only : n_species, calculate_all_eigenstates, n_max_basis_fns, &
                              n_atoms, n_basis, n_empty_atoms, n_max_grid, &
                              n_max_spline, n_pp_atoms, n_states, n_states_save, &
                              use_basis_gradients, use_ext_basis, &
                              use_hydro, use_ionic, use_plus_u, l_wave_max, n_basis_fns
  use basis,           only : basis_wave_spl, basis_wave_s_spl, basis_kinetic_spl, &
                              basis_kinetic_s_spl, basis_deriv_spl, basis_deriv_s_spl, &
                              basis_kinetic_scaled_zora_spl, basis_mapping, basis_fn_atom, &
                              basis_fn_start_spl, basis_wave_ordered, basis_wave_s_ordered, &
                              basis_kinetic_ordered, basis_kinetic_s_ordered, &
                              basis_kinetic_scaled_zora_ordered, basis_deriv_ordered, &
                              basis_deriv_s_ordered, outer_radius, outer_radius_sq, &
                              atom_radius_sq, basis_fn_atom, spline_offset, atom2basis_off, &
                              atom_radius, basis_fn, basis_small_fn, basis_l, basis_k, basis_small_l, &
                              basis_small_k, basis_m, basis_small_m, basis_atom, basis_small_atom, &
                              basisfn_l, n_basis_fn_species, Lsp2basis_fn, Lsp2basis_sp, &
                              Lsp2n_basis_fnLsp, max_basis_L, max_n_basis_fnLsp, &
                              max_n_basis_sp, sp2n_basis_sp, basis_fn_atom, basisfn_n, basisfn_k, &
                              basisfn_species, basisfn_type, perm_basis_fns_spl, &
                              i_radial_fn, allocate_basis
  use species_data,    only : n_ionic, n_hydro, l_shell_max, &
                              ionic_in_large_basis, hydro_in_large_basis, &
                              ext_l_shell_max, &
                              multipole_radius_free, multipole_radius_free_sq, &
                              no_basis, innermost_max, r_cutoff, include_min_basis, &
                              atomic_l, atomic_n, atomic_k, atomic_outer_radius, core_n_max, &
                              atomic_eigenval, atomic_wave, atomic_kinetic, atomic_kinetic_small, &
                              atomic_wave_deriv, atm_wave_large, atm_wave_small, &
                              atomic_large_deriv, atomic_small_deriv, species_pseudoized, species_name, &
                              n_atomic, n_conf, conf_l, conf_n, conf_kappa, confined_outer_radius, &
                              confined_eigenval, confined_wave, confined_kinetic, &
                              confined_wave_deriv, confined_wave_large, &
                              confined_wave_small, n_conf, ionic_l, ionic_n, ionic_kappa, &
                              ionic_outer_radius, ionic_eigenval, ionic_wave, &
                              ionic_kinetic, ionic_kinetic_small, ionic_wave_deriv, ionic_wave_large, &
                              ionic_wave_small, ionic_large_deriv, ionic_small_deriv, &
                              hydro_l, hydro_n, hydro_kappa, hydro_outer_radius, hydro_wave, &
                              hydro_kinetic, hydro_kinetic_small, hydro_wave_deriv, &
                              hydro_large_deriv, hydro_small_deriv, hydro_wave_large, hydro_wave_small, &
                              n_gaussian, gaussian_l, gaussian_n, gaussian_outer_radius, &
                              gaussian_wave, gaussian_kinetic, gaussian_wave_deriv, &
                              n_gaussian, basis_acc, atoms_in_structure, &
                              basis_dep_cutoff_thresh, n_sto, sto_n, sto_l, sto_k, &
                              sto_wave, sto_wave_small, sto_kinetic, sto_kinetic_small, sto_wave_deriv
  use aims_memory_tracking, only: aims_allocate
  use constants,       only : light_speed_sq, sqrt_pi
  use timing,          only : current_atom_num_bsse
  use grids,           only : r_grid, n_grid, r_radial, r_grid_inc
  use geometry,        only : species
  use spline,          only : val_spline_deriv, cubic_spline
  use mpi_tasks,       only : myid
  use localorb_io,     only : use_unit, OL_norm, localorb_info
  use free_atoms,      only : free_potential_spl
  use psi_at_nucleus_mod, only: psi_at_nucleus_atomic, psi_at_nucleus_hydro, &
       & psi_at_nucleus_gauss, psi_at_nucleus_ionic, psi_at_nucleus
  use rel_x2c_mod,     only : n_basis_small

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
  !  imported variables
  !     input
  !     most "input" is found in module species_data
  !     output
  !     all output is found in module basis
  !  local variables
  !     n_max : aux variables to determine loop counter maxima
  !     i_first_fn : First basis function for current atom and current angular momentum shell
  !     i_reject_wave : If the innermost maximum of a given wave function lies inside of
  !                  r_radial(i_reject_wave), then that wave function will not be used.  
  !     basis_l_max : this will save the maximum l number for each species in the normal basis
  !                   to replance l_shell_max values after this routine. l_shell_max contains 
  !                   l number also from possible auxiliary basis set.
  integer   :: n_max
  integer   :: i_first_fn, i_last_fn
  character :: l_shell_str
  integer   :: i_reject_wave
  real*8    :: r_inner_max(n_max_basis_fns)
  integer   :: function_index( n_species, n_max_basis_fns)
  integer   :: function_index_last
  integer   :: basis_l_max(n_species)
  ! Arrays for intermediate storage/handling of radial functions before spline
  ! basis_wave is compacted array of basis functions
  ! For relativistic cases, basis_wave saves the large component basis, while basis_wave_s saves the small component part.
  ! basis_deriv is compacted array of basis function derivatives (if needed)
  ! basis_kinetic : directly calculated kinetic energy term for all basis functions
  ! which are not free-atom like basis functions. [Treat free-atom like basis functions
  ! separately so we can subtract out the free-atom basis-defining potential later.]
  real*8, dimension(:,:), allocatable :: basis_wave
  real*8, dimension(:,:), allocatable :: basis_wave_s
  real*8, dimension(:,:), allocatable :: basis_deriv
  real*8, dimension(:,:), allocatable :: basis_deriv_s
  real*8, dimension(:,:), allocatable :: basis_kinetic
  real*8, dimension(:,:), allocatable :: basis_kinetic_s
  real*8, dimension(:,:), allocatable :: basis_kinetic_scaled_zora
  !     Arrays for a possible reindexing of all spline arrays
  integer :: i_species_index(n_max_basis_fns)
  integer :: perm_basis_fns(n_max_basis_fns)
  integer :: perm_basis_fns_inv(n_max_basis_fns)
  real*8  :: temp_radius(n_max_basis_fns)
  integer :: n_function, i_offset_spl, i_basis_2
  integer :: i_offset_spl_density
  integer :: i_offset
  integer :: species_offset
  integer :: outer_point
  !  counters
  integer :: i_species, i_basis, i_basis_small, i_function, i_atom, i_shell
  integer :: i_m, i_l, i_grid, i_type, i_function2
  integer :: i_ionic, i_conf, i_hydro, i_gaussian, i_sto
  integer :: i_atomic
  integer :: i_spline
  !  functions
  real*8    :: int_log_mesh
  character :: l_to_str
  real*8    :: get_inner_max
  integer,     dimension(n_species)       :: species_first_function, species_last_function   
  integer,     dimension(n_max_basis_fns) :: fn_species
  real*8,      dimension(n_max_basis_fns) :: fn_eigenval
  real*8,      dimension(n_max_basis_fns) :: fn_cutoff_radius
  character*8, dimension(n_max_basis_fns) :: fn_type
  integer,     dimension(n_max_basis_fns) :: fn_l
  integer,     dimension(n_max_basis_fns) :: fn_n
  integer,     dimension(n_max_basis_fns) :: fn_k
  character*16                            :: output_name
  real*8,      dimension(n_max_basis_fns) :: scalar_prod
  integer                                 :: i_prev
  real*8 :: write_radius, write_function, write_i_r
  ! Paula: these are for ordering the basis functions, so that the memory use in
  ! packing would be smaller.
  integer,    allocatable, dimension(:)   :: perm_basis, perm_basis_inv
  integer,    allocatable, dimension(:)   :: work, work_k
  character*8,allocatable, dimension(:)   :: work_c
  real*8,     allocatable, dimension(:,:) :: work_r2D
  real*8,     allocatable, dimension(:)   :: work_r
  real*8,     allocatable                 :: work_nuc0(:)
  real*8,     allocatable, dimension(:)   :: radius
  real*8 :: V_radial_deriv, basis_radial_deriv
  character*100 :: info_str
  logical :: need_section
  
  ! for atom bsse: may need cleanup
  integer :: j_basis, basis_counter, empty_basis  
                                        real*8 :: scalar


  call localorb_info("Assembling full basis from fixed parts.",use_unit,'(2X,A)')
  !     allocations - local variables for basis function storage
  !     basis_deriv is strictly _not_ needed unless use_basis_gradients is true.
  !                 Allocate it here anyway (perhaps unnecessarily) so that
  !                 orthonormalize_basis_fn can be called properly
  allocate ( basis_wave    (n_max_grid, n_max_basis_fns) )
  allocate ( basis_deriv   (n_max_grid, n_max_basis_fns) )
  allocate ( basis_kinetic (n_max_grid, n_max_basis_fns) )
  if (.not. allocated(psi_at_nucleus)) call aims_allocate(psi_at_nucleus, n_max_basis_fns, 'shrink_fixed_basis_phi_thresh::psi_at_nucleus')
  psi_at_nucleus = 0d0
  if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
    allocate(basis_wave_s (n_max_grid, n_max_basis_fns))
    allocate(basis_deriv_s(n_max_grid, n_max_basis_fns))
    allocate(basis_kinetic_s(n_max_grid, n_max_basis_fns))
    basis_wave_s=0.d0; basis_deriv_s=0.d0; basis_kinetic_s=0.d0
  endif
  if (flag_rel .eq. REL_atomic_zora) then
     allocate( basis_kinetic_scaled_zora(n_max_grid,n_max_basis_fns))
  end if

  !     determine basis _functions_ first
  !     i_function counts distinct radial functions (at, type, n, l)
  i_function             = 0
  i_basis                = 0
  i_basis_small          = 0
  n_basis                = 0
  n_basis_small          = 0
  fn_eigenval            = 0.d0
  function_index         = 0
  species_first_function = 0
  species_last_function  = 0
  r_inner_max(:)         = 0d0
  basis_l_max(:)         = 0
  do i_species = 1, n_species, 1
  if(.not.species_pseudoized(i_species).and.(.not.(no_basis(i_species)))) then
     ! species dependent counters & such
     species_first_function(i_species) = i_function + 1           ! remember where in the basis function list a certain species starts.
     i_reject_wave                     = innermost_max(i_species) ! index for wave function rejection
     i_basis                           = 0                        ! radial function counter for species
     i_basis_small                     = 0
     fn_cutoff_radius                  = r_cutoff(i_species)      ! initialize cutoff radius for potential sorting later: this is the biggest value
     ! reat each l shell as a block
     do i_l = 0, l_shell_max(i_species), 1
        i_type = 0
        i_first_fn = i_function+1    ! first function in a given l-shell of a given species, need to remember that throughout this loop!

        if (include_min_basis(i_species)) then
           ! treat atomic-like wave functions first
           i_type = i_type + 1
           do i_atomic = 1, n_atomic(i_species), 1
              if (atomic_l(i_species,i_atomic).eq.i_l) then
                 i_function                   = i_function+1
                 fn_species(i_function)       = i_species
                 fn_type(i_function)          = "atomic"
                 fn_l(i_function)             = i_l
                 fn_n(i_function)             = atomic_n(i_species,i_atomic)
                 if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) fn_k(i_function) = atomic_k(i_species,i_atomic)
                 fn_cutoff_radius(i_function) = atomic_outer_radius(i_species,i_atomic)
                 if (basis_l_max(i_species) .le. i_l) then
                 !  save species maximun l value
                    basis_l_max(i_species) = i_l
                 end if
                 if (core_n_max(i_l,i_species) .ge. atomic_n(i_species,i_atomic)) then
                    ! only for our own relativity, store the correct eigenvalue -
                    ! for all but core functions, this has to be zero.
                    fn_eigenval(i_function) = atomic_eigenval(i_species,i_atomic)
                 end if
                 ! store current basis function in basis_wave
                 ! also store kinetic energy term of basis function
                 do i_grid = 1, n_grid(i_species),1
                   if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
                     ! For relativistic case, basis_wave saves the large component basis, while basis_wave_s saves the small component part.
                     basis_wave(i_grid, i_function) = atm_wave_large(i_grid, i_species, i_atomic)
                     basis_wave_s(i_grid, i_function) = atm_wave_small(i_grid, i_species, i_atomic)
                     basis_kinetic(i_grid, i_function) = atomic_kinetic(i_grid, i_species, i_atomic)
                     basis_kinetic_s(i_grid, i_function) = atomic_kinetic_small(i_grid, i_species, i_atomic)
                   else
                     basis_wave(i_grid, i_function) = atomic_wave(i_grid, i_species, i_atomic)
                     basis_kinetic(i_grid, i_function) = atomic_kinetic(i_grid, i_species, i_atomic)
                   endif
                 enddo
                 psi_at_nucleus(i_function) = psi_at_nucleus_atomic(i_atomic,i_species)
                 ! psi_at_nucleus_atomic is undefined with sratom. This is an inaccurate workaround.
                 if (atomic_solver == ATOMIC_SOLVER_SRATOM) psi_at_nucleus(i_function) = basis_wave(1,i_function)/r_grid(1,i_species)
                 if (use_basis_gradients) then
                    do i_grid = 1, n_grid(i_species),1
                       basis_deriv(i_grid, i_function) = atomic_wave_deriv(i_grid, i_species, i_atomic)
                    enddo
                 end if
                 if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
                    do i_grid = 1, n_grid(i_species),1
                       basis_deriv(i_grid,i_function)   = atomic_large_deriv(i_grid, i_species, i_atomic)
                       basis_deriv_s(i_grid,i_function) = atomic_small_deriv(i_grid, i_species, i_atomic)
                    enddo
                 endif
                !if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
                !  r_inner_max(i_function) = min ( get_inner_max ( basis_wave(1,i_function), r_grid(1,i_species), n_grid(i_species) ), &
                !                            get_inner_max ( basis_wave_s(1,i_function), r_grid(1,i_species), n_grid(i_species) ) )
                !else
                   r_inner_max(i_function) = get_inner_max ( basis_wave(1,i_function), r_grid(1,i_species), n_grid(i_species) )
                !endif

                 ! if a minimal basis function does not meet the criterion, stop the calculation
                 ! entirely - the grids must be chosen differently.
                 if (r_inner_max(i_function).lt.r_radial(i_reject_wave,i_species)) then
                    l_shell_str = l_to_str(i_l)
                    if (myid.eq.0) then
                       write(use_unit,'(1X,A,A)') "* Warning - the innermost maximum ", "of the minimal basis function"
                       write(use_unit,'(1X,A,I2,3A)') "* ", atomic_n(i_species,i_atomic), l_shell_str, " of species ",trim(species_name(i_species))
                       write(use_unit,'(1X,A,I1,A)') "* has its innermost maximum inside the ", i_reject_wave, "th radial integration shell."
                       write(use_unit,'(1X,A,A)') "* Adjust your integration grid to be ", "accurate enough before continuing."
                    end if
                    stop
                 end if
              end if
           enddo
        end if

        !     treat confined basis functions next
        if (n_conf(i_species).gt.0) then
           i_type = i_type+1
           do i_conf = 1, n_conf(i_species), 1
              if (conf_l(i_species,i_conf).eq.i_l) then
                 i_function                   = i_function+1
                 i_shell                      = conf_n(i_species, i_conf)
                 fn_species(i_function)       = i_species
                 fn_type(i_function)          = "confined"
                 fn_l(i_function)             = i_l
                 fn_n(i_function)             = i_shell
                 if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) fn_k(i_function) = conf_kappa(i_species,i_conf)
                 fn_cutoff_radius(i_function) = confined_outer_radius(i_species, i_conf)
                 if (basis_l_max(i_species) .lt. i_l) then
                 !  save species maximun l value
                    basis_l_max(i_species) = i_l
                 end if
                 if (core_n_max(i_l,i_species) .ge. i_shell) then
                    ! only for our own relativity, store the correct eigenvalue -
                    ! for all but core functions, this has to be zero.
                    fn_eigenval(i_function) = confined_eigenval(i_species,i_conf)
                 end if
                 ! copy current basis function to basis_wave
                 do i_grid = 1, n_grid(i_species), 1
                   if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
                     ! For relativistic case, basis_wave saves the large component basis, while basis_wave_s saves the small component part.
                     basis_wave(i_grid, i_function) = confined_wave_large(i_grid, i_species, i_conf)
                     basis_wave_s(i_grid, i_function) = confined_wave_small(i_grid, i_species, i_conf)
                     basis_kinetic(i_grid, i_function) = confined_kinetic(i_grid, i_species, i_conf)
                   else
                     basis_wave(i_grid, i_function)    = confined_wave(i_grid, i_species, i_conf)
                     basis_kinetic(i_grid, i_function) = confined_kinetic(i_grid, i_species, i_conf)
                   endif
                 enddo
                 if (use_basis_gradients) then
                    do i_grid = 1, n_grid(i_species), 1
                       basis_deriv(i_grid, i_function) = confined_wave_deriv(i_grid, i_species, i_conf)
                    enddo
                 end if
                 r_inner_max(i_function) = get_inner_max ( basis_wave(1,i_function), r_grid(1,i_species), n_grid(i_species) )
              end if
           enddo
        end if           ! end confined function

        ! treat ionic basis functions next
        if (use_ionic .and. n_ionic(i_species).gt.0) then
           ! JW: Avoid test for ionic_in_large_basis if .not. use_ionic.
           need_section = (.not. use_ext_basis) .or. any(.not. ionic_in_large_basis(i_species, 1:n_ionic(i_species)))
        else
           need_section = .false.
        end if

        if (need_section) then
           i_type = i_type+1
           do i_ionic = 1, n_ionic(i_species), 1
              if (ionic_l(i_species, i_ionic).eq.i_l .and. (.not. ionic_in_large_basis(i_species,i_ionic)) ) then
                 i_function                   = i_function+1
                 i_shell                      = ionic_n(i_species, i_ionic)
                 fn_species(i_function)       = i_species
                 fn_type(i_function)          = "ionic"
                 fn_l(i_function)             = i_l
                 fn_n(i_function)             = i_shell
                 if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) fn_k(i_function) = ionic_kappa(i_species,i_ionic)
                 fn_cutoff_radius(i_function) = ionic_outer_radius(i_species, i_ionic)
                 if (basis_l_max(i_species) .lt. i_l) then
                 !  save species maximun l value
                    basis_l_max(i_species) = i_l
                 end if
                 if (core_n_max(i_l,i_species) .ge. ionic_n(i_species,i_ionic)) then
                    ! only for our own relativity, store the correct eigenvalue -
                    ! for all but core functions, this has to be zero.
                    fn_eigenval(i_function) = ionic_eigenval(i_species,i_ionic)
                 end if
                 ! copy current basis function to basis_wave
                 do i_grid = 1, n_grid(i_species), 1
                   if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
                     ! For relativistic case, basis_wave saves the large component basis, while basis_wave_s saves the small component part.
                     basis_wave(i_grid, i_function) = ionic_wave_large(i_grid, i_species, i_ionic)
                     basis_wave_s(i_grid, i_function) = ionic_wave_small(i_grid, i_species, i_ionic)
                     basis_kinetic(i_grid, i_function) = ionic_kinetic(i_grid, i_species, i_ionic)
                     basis_kinetic_s(i_grid, i_function) = ionic_kinetic_small(i_grid, i_species, i_ionic)
                   else
                     basis_wave(i_grid, i_function) = ionic_wave(i_grid, i_species, i_ionic)
                     basis_kinetic(i_grid, i_function) = ionic_kinetic(i_grid, i_species, i_ionic)
                   endif
                 enddo
                 psi_at_nucleus(i_function) = psi_at_nucleus_ionic(i_ionic,i_species)
                 if (use_basis_gradients) then
                    do i_grid = 1, n_grid(i_species), 1
                       basis_deriv(i_grid, i_function) = ionic_wave_deriv(i_grid, i_species, i_ionic)
                    enddo
                 end if
                 if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
                    do i_grid = 1, n_grid(i_species), 1
                       basis_deriv(i_grid,i_function)   = ionic_large_deriv(i_grid, i_species, i_ionic)
                       basis_deriv_s(i_grid,i_function) = ionic_small_deriv(i_grid, i_species, i_ionic)
                    enddo
                 endif
                 r_inner_max(i_function) = get_inner_max ( basis_wave(1,i_function), r_grid(1,i_species), n_grid(i_species) )
              end if
           enddo
           !     end ionic function
        end if
        
        !     treat hydrogenic basis functions next
        !     check also that does all hydrogenic function belong to the large or small basis
        if (use_hydro .and. n_hydro(i_species).gt.0) then 
           ! JW: Avoid test for hydro_in_large_basis if .not. use_hydro.
           need_section = (.not. use_ext_basis) .or. any(.not. hydro_in_large_basis(i_species, 1:n_hydro(i_species)))
        else
           need_section = .false.
        end if
        if (need_section) then
           i_type = i_type+1
           do i_hydro = 1, n_hydro(i_species), 1
              if (hydro_l(i_species, i_hydro).eq.i_l .and. .not. hydro_in_large_basis(i_species,i_hydro)) then
                 i_function                   = i_function+1
                 i_shell                      = hydro_n(i_species, i_hydro)
                 fn_species(i_function)       = i_species
                 fn_type(i_function)          = "hydro"
                 fn_l(i_function)             = i_l
                 fn_n(i_function)             = i_shell
                 if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) fn_k(i_function) = hydro_kappa(i_species,i_hydro)
                 fn_cutoff_radius(i_function) = hydro_outer_radius(i_species,i_hydro)
                 if (basis_l_max(i_species) .lt. i_l) then
                 !  save species maximun l value
                    basis_l_max(i_species) = i_l
                 end if
                 !     copy current basis function to basis_wave
                 do i_grid = 1, n_grid(i_species), 1
                   if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
                     ! For relativistic case, basis_wave saves the large component basis, while basis_wave_s saves the small component part.
                     basis_wave(i_grid, i_function) = hydro_wave_large(i_grid, i_species, i_hydro)
                     basis_wave_s(i_grid, i_function) = hydro_wave_small(i_grid, i_species, i_hydro)
                     basis_kinetic(i_grid, i_function) = hydro_kinetic(i_grid, i_species, i_hydro)
                     basis_kinetic_s(i_grid, i_function) = hydro_kinetic_small(i_grid, i_species, i_hydro)
                   else
                     basis_wave(i_grid, i_function) = hydro_wave(i_grid, i_species, i_hydro)
                     basis_kinetic(i_grid, i_function) = hydro_kinetic(i_grid, i_species, i_hydro)
                   endif
                 enddo
                 psi_at_nucleus(i_function) = psi_at_nucleus_hydro(i_hydro,i_species)
                 if (use_basis_gradients) then
                    do i_grid = 1, n_grid(i_species), 1
                       basis_deriv(i_grid, i_function) = hydro_wave_deriv(i_grid, i_species, i_hydro)
                    enddo
                 end if
                 if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
                    do i_grid = 1, n_grid(i_species), 1
                       basis_deriv(i_grid,i_function)   = hydro_large_deriv(i_grid, i_species, i_hydro)
                       basis_deriv_s(i_grid,i_function) = hydro_small_deriv(i_grid, i_species, i_hydro)
                    enddo
                 endif
                 r_inner_max(i_function) = get_inner_max ( basis_wave(1,i_function), r_grid(1,i_species), n_grid(i_species) )
              end if
           enddo
           !         end hydrogenic function
        end if

        !         treat Gaussian basis functions next
        if (n_gaussian(i_species).gt.0) then
           i_type = i_type+1
           do i_gaussian = 1, n_gaussian(i_species), 1
              if (gaussian_l(i_species,i_gaussian).eq.i_l) then
                 i_function             = i_function + 1
                 i_shell                = gaussian_n(i_species, i_gaussian)
                 fn_species(i_function) = i_species
                 fn_type(i_function)    = "gaussian"
                 fn_l(i_function)       = i_l
                 fn_n(i_function)       = i_shell
                !if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) fn_k(i_function) = gaussian_kappa(i_species,i_gaussian)
                 fn_cutoff_radius(i_function) = gaussian_outer_radius(i_species,i_gaussian)
                 if (basis_l_max(i_species) .lt. i_l) then
                 !  save species maximun l value
                    basis_l_max(i_species) = i_l
                 end if
                 !               copy current basis function to basis_wave
                 do i_grid = 1, n_grid(i_species), 1
                    basis_wave(i_grid, i_function) = gaussian_wave(i_grid, i_species, i_gaussian)
                    basis_kinetic(i_grid, i_function) = gaussian_kinetic(i_grid, i_species, i_gaussian)
                 enddo
                 psi_at_nucleus(i_function) = psi_at_nucleus_gauss(i_gaussian,i_species)
                 if (use_basis_gradients) then
                    do i_grid = 1, n_grid(i_species), 1
                       basis_deriv(i_grid, i_function) = gaussian_wave_deriv(i_grid, i_species, i_gaussian)
                    enddo
                 end if
                 r_inner_max(i_function) = get_inner_max ( basis_wave(1,i_function), r_grid(1,i_species), n_grid(i_species) )
              end if
           enddo
           !         end Gaussian function
        end if

        if (n_sto(i_species) > 0) then
           i_type = i_type + 1
           do i_sto = 1, n_sto(i_species)
              if (sto_l(i_species,i_sto) == i_l) then
                 i_function             = i_function + 1
                 i_shell                = sto_n(i_species, i_sto)
                 fn_species(i_function) = i_species
                 fn_type(i_function)    = 'sto'
                 fn_l(i_function)       = i_l
                 fn_n(i_function)       = i_shell
                 if (flag_rel == REL_x2c .or. flag_rel == REL_4c_dks) &
                      & fn_k(i_function) = sto_k(i_species,i_sto)
                 if (basis_l_max(i_species) < i_l) basis_l_max(i_species) = i_l
                 if(flag_rel == REL_x2c .or. flag_rel == REL_4c_dks)then
                    ! For relativistic case, basis_wave saves the
                    ! large component basis, while basis_wave_s saves
                    ! the small component part.
                    basis_wave(:n_grid(i_species), i_function) = &
                         & sto_wave(:n_grid(i_species), i_species, i_sto)
                    basis_wave_s(:n_grid(i_species), i_function) = &
                         & sto_wave_small(:n_grid(i_species), i_species, i_sto)
                    basis_kinetic(:n_grid(i_species), i_function) = &
                         & sto_kinetic(:n_grid(i_species), i_species, i_sto)
                    basis_kinetic_s(:n_grid(i_species), i_function) = &
                         & sto_kinetic_small(:n_grid(i_species), i_species, i_sto)
                 else
                    basis_wave(:n_grid(i_species), i_function) = &
                         & sto_wave(:n_grid(i_species), i_species, i_sto)
                    basis_kinetic(:n_grid(i_species), i_function) = &
                         & sto_kinetic(:n_grid(i_species), i_species, i_sto)
                 endif
                 ! psi_at_nucleus(i_function) = psi_at_nucleus_sto(i_sto,i_species)
                 if (use_basis_gradients) &
                      & basis_deriv(:n_grid(i_species), i_function) = &
                      & sto_wave_deriv(:n_grid(i_species), i_species, i_sto)
                 if (flag_rel == REL_x2c .or. flag_rel == REL_4c_dks) then
                    r_inner_max(i_function) = &
                         & min(get_inner_max(basis_wave(1,i_function), &
                         & r_grid(1,i_species), n_grid(i_species)), &
                         & get_inner_max(basis_wave_s(1,i_function), &
                         & r_grid(1,i_species), n_grid(i_species)))
                 else
                    r_inner_max(i_function) = &
                         & get_inner_max(basis_wave(1,i_function), &
                         & r_grid(1,i_species), n_grid(i_species))
                 end if
              end if
           end do
        end if

        ! now have the basis functions and their associated information sorted into different arrays. 
        ! sort the cutoff potentials for each l-shell and then normalize according to cutoff!
        ! i_first_fn indexes the first function in a certain l-shell
        ! i_function indexes the last function in a given l-shell 
        i_last_fn  = i_function
        i_function = i_first_fn

        if (use_plus_u) then
           if (basis_dep_cutoff_thresh(i_species) .ne. 0d0) then
              basis_dep_cutoff_thresh(i_species) = 0d0
              write(info_str,'(1X,A)') "* Warning - Since '+U' treatment was requested,"
              call localorb_info(info_str)
              write(info_str,'(1X,A)') "* the basis-dependent cutoff is turned off automatically."
              call localorb_info(info_str)
           end if
        end if

        if (basis_dep_cutoff_thresh(i_species).gt.0d0) then
           ! allocate temporary storage arrays and switches if not already done so
           if (.not.allocated(perm_basis))     allocate(perm_basis    (n_max_basis_fns))
           if (.not.allocated(perm_basis_inv)) allocate(perm_basis_inv(n_max_basis_fns))
           if (.not.allocated(work))           allocate(work          (n_max_basis_fns))             ! integer
           if (.not.allocated(work_r))         allocate(work_r        (n_max_basis_fns))             ! real
           if (.not.allocated(work_c))         allocate(work_c        (n_max_basis_fns))             ! character
           if (.not.allocated(work_r2D))       allocate(work_r2D      (n_max_grid, n_max_basis_fns)) ! basis function
           ! wavefunction values at nucleus
           if (.not.allocated(work_nuc0))      allocate(work_nuc0     (n_max_basis_fns))
           if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
              if (.not.allocated(work_k))      allocate(work_k        (n_max_basis_fns))             ! kappa value for relativistic basis
           endif
           ! sort single l-shell according to array fn_cutoff_radius: from i_function to i_last_fn         
           call insertionsort( fn_cutoff_radius(i_function:i_last_fn), i_last_fn-i_function+1, &
                perm_basis(i_function:i_last_fn), perm_basis_inv(i_function:i_last_fn) )
           ! all entries in perm_basis start counting @ 1 - but I would really like them to start at i_function. 
           perm_basis_inv(i_function:i_last_fn) = perm_basis_inv(i_function:i_last_fn) + i_function - 1 
           ! permute all other necessary data according to sorting:
           do i_function2 = i_function, i_last_fn
              work      (i_function2) = fn_n        (perm_basis_inv(i_function2))
              work_r    (i_function2) = r_inner_max (perm_basis_inv(i_function2))
              work_c    (i_function2) = fn_type     (perm_basis_inv(i_function2))
              work_nuc0 (i_function2) = psi_at_nucleus(perm_basis_inv(i_function2))
              work_r2D(:,i_function2) = basis_wave(:,perm_basis_inv(i_function2)) 
              if((flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks))then
                  work_k(i_function2) = fn_k(perm_basis_inv(i_function2))
              endif
           end do
           fn_n          (i_function:i_last_fn)   = work      (i_function:i_last_fn)
           fn_type       (i_function:i_last_fn)   = work_c    (i_function:i_last_fn)
           r_inner_max   (i_function:i_last_fn)   = work_r    (i_function:i_last_fn)
           psi_at_nucleus(i_function:i_last_fn)   = work_nuc0 (i_function:i_last_fn)
           basis_wave    (:,i_function:i_last_fn) = work_r2D(:,i_function:i_last_fn)
           if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
             fn_k(i_function:i_last_fn) = work_k(i_function:i_last_fn)
             do i_function2 = i_function, i_last_fn
                work_r2D(:,i_function2) = basis_wave_s(:,perm_basis_inv(i_function2)) 
             end do
             basis_wave_s(:,i_function:i_last_fn) = work_r2D(:,i_function:i_last_fn)
           endif
            ! for basis derivative:
           do i_function2 = i_function, i_last_fn
              work_r2D(:,i_function2) = basis_deriv(:,perm_basis_inv(i_function2))
           end do
           basis_deriv(:,i_function:i_last_fn) = work_r2D(:,i_function:i_last_fn)
           if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
             do i_function2 = i_function, i_last_fn
                work_r2D(:,i_function2) = basis_deriv_s(:,perm_basis_inv(i_function2)) 
             end do
             basis_deriv_s(:,i_function:i_last_fn) = work_r2D(:,i_function:i_last_fn)
           endif
            ! for kinetic basis
           do i_function2 = i_function, i_last_fn
              work_r2D(:,i_function2) = basis_kinetic(:,perm_basis_inv(i_function2))
           end do
           basis_kinetic(:,i_function:i_last_fn) = work_r2D(:,i_function:i_last_fn)
           if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
             do i_function2 = i_function, i_last_fn
                work_r2D(:,i_function2) = basis_kinetic_s(:,perm_basis_inv(i_function2)) 
             end do
             basis_kinetic_s(:,i_function:i_last_fn) = work_r2D(:,i_function:i_last_fn)
           endif
        end if

        ! keep for output in this particular l-shell
        l_shell_str = l_to_str(i_l)
        
        ! loop through all basis functions in the current l-shell, to orthonormalize the bunch according to current sorting 
        ! take special attention to the fact that their number might be reduced as we reject one here or there. 
        i_function = i_first_fn
        do while (i_function.le.i_last_fn)

           ! check inner radius
           if (r_inner_max(i_function).ge.r_radial(i_reject_wave,i_species)) then
              ! if alright, orthonormalize with the rest
              call orthonormalize_basis_fn &
                 ( i_function, i_first_fn, basis_wave, n_grid(i_species), r_grid(1,i_species), basis_acc(i_species), &
                 basis_kinetic, basis_deriv, function_index(i_species, i_function), .true., psi_at_nucleus)
           else
              function_index(i_species,i_function) = -1   ! if not, reject function
           end if           
           ! evaluate function according to the outcome of the orthornomalization procedure
           if (function_index(i_species,i_function).gt.0) then
              ! accepted, that's good I suppose
              write(info_str,'(2X,3A,A8,A,I3,3A)') &
                   "| Species ", trim(species_name(i_species)), " : ", trim(fn_type(i_function)), &
                   " orbital ", fn_n(i_function)," ",l_shell_str," accepted."
              call localorb_info(info_str,use_unit,'(A)')                   
              i_basis = i_basis + (2*i_l+1)
              if((flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks))then
               ! For relativistic cases, the n_basis here denotes the scalar large comp. basis, 
               ! and the n_basis_small denotes its small comp. counterpart.
               ! For detailed interpretation, see the annotations in 
               ! /src/relativity/rel_integrations.f90/prune_radial_basis_rel
                 if(i_l.ne.0)then
                    if(fn_k(i_function).eq.i_l)    i_basis_small = i_basis_small + 2*(i_l-1)+1  ! spin down
                    if(fn_k(i_function).eq.-i_l-1) i_basis_small = i_basis_small + 2*(i_l+1)+1  ! spin up
                 else
                    ! For s orbital, the number of large comp. does not increase, 
                    ! but the number of small comp. increases by 3 (corresponding to py, pz, px).
                    i_basis_small = i_basis_small + 3
                 endif
              endif
              i_function = i_function + 1
           else 
              if ( function_index(i_species,i_function).eq.0 ) then
                 write(info_str,'(2X,3A,A8,A,I3,1X,A,A)') &
                      "| Species ", trim(species_name(i_species)), " : ", trim(fn_type(i_function)), &
                      " orbital ", fn_n(i_function), l_shell_str, " rejected: Linear dependence."
                 call localorb_info(info_str,use_unit,'(A)')
              else
                 write(info_str,'(2X,3A,A8,A,I3,1X,A,A)') &
                      "| Species ", trim(species_name(i_species)), " : ",trim(fn_type(i_function)),&
                      " orbital ", fn_n(i_function), l_shell_str, " rejected: Insufficient integration grid."
                 call localorb_info(info_str,use_unit,'(A)')
              end if
              ! Function got rejected outright. 
              ! if this was an atomic function, there is something seriously wrong. Notify user and quit.
              if (fn_type(i_function).eq.'atomic') then
                 call localorb_info("* WARNING: Attempting to reject atomic basis function.",use_unit,'(A)')
                 call localorb_info("*          You may be running with a very large basis set, ",use_unit,'(A)')
                 call localorb_info("*          but non-zero basis_dep_cutoff in control.in.",use_unit,'(A)')
                 call localorb_info("*          If so, set basis_dep_cutoff to zero in control.in and retry.",use_unit,'(A)')
                 call localorb_info("*          Please also consult the manual on basis_dep_cutoff .",use_unit,'(A)')
                 call localorb_info("*          Check settings and restart.",use_unit,'(A)')
                 stop
              end if
              ! Before continuing, sort every other function in this
              ! l-shell back one step. In particular, this requires moving ...
              !    basis_wave(n_max_grid, n_max_basis_fns)
              !    basis_kinetic(n_max_grid, n_max_basis_fns)
              !    basis_deriv(n_max_grid, n_max_basis_fns)
              do i_function2 = i_function + 1, i_last_fn
                 basis_wave   (:,i_function2-1) = basis_wave   (:,i_function2)
                 basis_kinetic(:,i_function2-1) = basis_kinetic(:,i_function2)
                 basis_deriv  (:,i_function2-1) = basis_deriv  (:,i_function2)
                 r_inner_max  (  i_function2-1) = r_inner_max  (  i_function2)
                 fn_n         (  i_function2-1) = fn_n         (  i_function2)
                 fn_type      (  i_function2-1) = fn_type      (  i_function2)
                 fn_l         (  i_function2-1) = fn_l         (  i_function2)
                 fn_species   (  i_function2-1) = fn_species   (  i_function2)
                 psi_at_nucleus(i_function2-1)  = psi_at_nucleus(i_function2)
                 if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
                   basis_wave_s   (:,i_function2-1) = basis_wave_s   (:,i_function2)
                   basis_deriv_s  (:,i_function2-1) = basis_deriv_s  (:,i_function2)
                   basis_kinetic_s(:,i_function2-1) = basis_kinetic_s(:,i_function2)
                   fn_k             (i_function2-1) = fn_k             (i_function2)
                 endif
              end do
              ! finally, the last function must be decreased by 1 as we took one out in the middle. 
              i_last_fn = i_last_fn - 1
           end if
        end do     ! loop over basis functions within a given l-shell
        i_function = i_function - 1  ! This points to the last known accepted function. 

     enddo      ! end loop over l-shells
     species_last_function(i_species) = i_function    ! remember where a species stops, required for labeling later.
     n_basis = n_basis + i_basis * atoms_in_structure(i_species)
     if (flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) n_basis_small = n_basis_small + i_basis_small * atoms_in_structure(i_species)
  end if   ! .not.species_pseudoized
  enddo    ! first loop over species 


  ! Normalize wavefunction value at the nucleus with the first spherical harmonic.
  psi_at_nucleus = psi_at_nucleus/2/sqrt_pi

  ! deallocate work array from function sorting within a single l-shell
  if (allocated(perm_basis))     deallocate(perm_basis)
  if (allocated(perm_basis_inv)) deallocate(perm_basis_inv)
  if (allocated(work))           deallocate(work)
  if (allocated(work_k))         deallocate(work_k)
  if (allocated(work_r))         deallocate(work_r)
  if (allocated(work_c))         deallocate(work_c)
  if (allocated(work_r2D))       deallocate(work_r2D)

  !     now we know the actual basis size, and hence all necessary array dimensions
  n_basis_fns = i_function

! >>> AB: feb 2012
!
!     in case, the overlap matrix and KS-eigenvectors are requested to be
!     written out to perform transport calculations as done by "aitranss"
!     module (Karlsruhe transport code), keyword 'output aitranss', all
!     KS-states must be calculated: i.e. n_states = n_basis independent on
!     the user's settings
  if (out_aitranss) then
!    put warning/info message ...
     if (myid.eq.0) then
        if (n_empty /= -1 .or. calculate_all_eigenstates) then  ! explicitly requested
           write(use_unit,*) ' '
           write(use_unit,'(2x,a,i8,a)') '* You have requested ', n_states, ' Kohn-Sham states in the calculation, however,'
           write(use_unit,'(2x,a)') '* for the output of the overlap matrix and KS-eigenvectors in a format compatible'
           write(use_unit,'(2x,a,i8,a)') '* with "aitranss" module, all ', n_basis, ' available KS-states will be calculated'
        else
           write(use_unit,*) ' '
           write(use_unit,'(2x,a)') '* For the output of the overlap matrix and KS-eigenvectors in a format compatible'
           write(use_unit,'(2x,a,i8,a)') '* with "aitranss" module, all ', n_basis, ' available KS-states will be computed'
        end if
     end if
!    >>> inserted by AB: march 2012
     n_states_save = n_states
     n_states = n_basis
!    <<< done with insert: AB, march 2012

! otherwise, verify n_states against n_basis
  else if (n_states.gt.n_basis) then
     if (myid.eq.0) then
        if (n_empty /= -1 .or. calculate_all_eigenstates) then    ! explicitly requested
           write(use_unit,"()")
           write(use_unit,"(2X,A,I8,A)") "* You requested ", n_states, " Kohn-Sham states in the calculation,"
           write(use_unit,"(2X,A,I8,A)") "* but there are only ", n_basis, " available basis states."
           write(use_unit,"(2X,A,A,I8,A)") "* Reducing total number of ", " Kohn-Sham states to ", n_basis, "."
           write(use_unit,*)
        else
           write(use_unit,"(2X,A,A,I8,A)") "Reducing total number of ", " Kohn-Sham states to ", n_basis, "."
        end if
     end if
     n_states = n_basis
     if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) n_states = (n_basis + n_basis_small)/2 ! = dim_matrix_rel*2
  end if
! DB: hack for a system only consisting of pseudocore+empty
  if((n_pp_atoms.gt.0).and.(n_empty_atoms.gt.0) .and.(n_atoms.eq.0))then
    n_states = n_basis
  end if
! <<< done with modification, AB: feb 2012

  ! VB: Deal with specialities for kinetic energy here. If atomic_zora needed, must amend
  !     basis_kinetic by various extra terms.
  !
  if (flag_rel .eq. REL_atomic_zora)then
     do i_function = 1, n_basis_fns, 1
        do i_grid = 1, n_grid(fn_species(i_function)),1
           ! first, obtain radial derivative of atomic potential
           V_radial_deriv =  val_spline_deriv( dble(i_grid), free_potential_spl(1,1, fn_species(i_function)), n_grid(fn_species(i_function)) ) &
                / (log(r_grid_inc( fn_species(i_function))) * r_grid(i_grid, fn_species(i_function)))
           ! ... and radial derivative of this radial function
           basis_radial_deriv = basis_deriv(i_grid, i_function)
           ! first, use unmodified kinetic energy part to create renormalization
           ! expression for "scaled ZORA"
           basis_kinetic_scaled_zora(i_grid,i_function) = &
             basis_kinetic(i_grid,i_function) * 2 * light_speed_sq / (  2 * light_speed_sq - free_potential_spl(1,i_grid,fn_species(i_function)))**2 &
             !
             -  2*light_speed_sq &
             !
             * V_radial_deriv * (   basis_radial_deriv - basis_wave(i_grid,i_function) / r_grid(i_grid, fn_species(i_function))) &
             !
             / (  2*light_speed_sq - free_potential_spl(1,i_grid,fn_species(i_function)))**3
           
           ! now, modify the kinetic energy part itself for on-site ("atomic") ZORA corrections
           basis_kinetic(i_grid,i_function) = &
             basis_kinetic(i_grid,i_function) * 2 * light_speed_sq / (  2 * light_speed_sq - free_potential_spl(1,i_grid,fn_species(i_function))) &
             !
             - light_speed_sq / (  2*light_speed_sq - free_potential_spl(1,i_grid,fn_species(i_function)))**2 &
             !
             * V_radial_deriv * (   basis_radial_deriv - basis_wave(i_grid,i_function) / r_grid(i_grid, fn_species(i_function)))
        end do
     end do
     
  end if ! REL_atomic_zora
  
  if (flag_rel .eq. REL_own) then
     do i_function = 1, n_basis_fns, 1
        ! amend kinetic energy by core eigenvalue (all other eigenvales are set to zero above)
        
        ! Notice that this treatment is still outrightly wrong. Really, we should be
        ! setting the correct kinetic energy term already when the basis functions are
        ! generated. Then, we could simply use the orthonormalizaton above to project out all core
        ! contributions correctly. Now, confined and ionic core radial functions will
        ! get the wrong kinetic energy ... as will any accidentally overlapping core contributions from
        ! the original non-orthonormalized basis functions.
        
        do i_grid = 1, n_grid(fn_species(i_function)),1
           ! first, obtain radial derivative of atomic potential
           V_radial_deriv =  val_spline_deriv( dble(i_grid), free_potential_spl(1,1, fn_species(i_function)), n_grid(fn_species(i_function)) ) &
                / (log(r_grid_inc( fn_species(i_function))) * r_grid(i_grid, fn_species(i_function)))
           ! ... and radial derivative of this radial function
           basis_radial_deriv = basis_deriv(i_grid, i_function)
           ! now, modify the kinetic energy part itself for on-site ("atomic") ZORA corrections
           basis_kinetic(i_grid,i_function) = &
             basis_kinetic(i_grid,i_function) * 2 * light_speed_sq / ( 2 * light_speed_sq + fn_eigenval(i_function) - free_potential_spl(1,i_grid,fn_species(i_function))) &
             !
             - light_speed_sq / (  2*light_speed_sq + fn_eigenval(i_function) - free_potential_spl(1,i_grid,fn_species(i_function)))**2 &
             !
             * V_radial_deriv * (   basis_radial_deriv - basis_wave(i_grid,i_function) / r_grid(i_grid, fn_species(i_function)))
        end do
     end do
  end if ! REL_own
  
  if (out_basis) then
     !       Output accepted basis functions for consistency
     do i_function = 1, n_basis_fns, 1
        l_shell_str = l_to_str(fn_l(i_function))
        if (i_function.lt.10) then
           if (fn_n(i_function).lt.10) then
              write(output_name,'(A2,A1,I1,A1,I1,A1,I1,A1,A1,A4)') &
                   fn_type(i_function),"_",fn_species(i_function),"_", i_function, "_", fn_n(i_function),"_",l_shell_str,".dat"
           else if  (fn_n(i_function).lt.100) then
              write(output_name,'(A2,A1,I1,A1,I1,A1,I2,A1,A1,A4)') &
                   fn_type(i_function),"_",fn_species(i_function),"_", i_function, "_", fn_n(i_function),"_",l_shell_str,".dat"
           end if
        else if (i_function.lt.100) then
           if (fn_n(i_function).lt.10) then
              write(output_name,'(A2,A1,I1,A1,I2,A1,I1,A1,A1,A4)') &
                   fn_type(i_function),"_",fn_species(i_function),"_", i_function, "_", fn_n(i_function),"_",l_shell_str,".dat"
           else if  (fn_n(i_function).lt.100) then
              write(output_name,'(A2,A1,I1,A1,I2,A1,I2,A1,A1,A4)') &
                   fn_type(i_function),"_",fn_species(i_function),"_", i_function, "_", fn_n(i_function),"_",l_shell_str,".dat"
           end if
        end if
        
        if (myid.eq.0) then
           open(50, file=output_name)
           write(50,*) "# ",n_grid(fn_species(i_function))
           do i_grid = 1, n_grid(fn_species(i_function)), 1
              write(50,*) r_grid(i_grid, fn_species(i_function)), basis_wave(i_grid, i_function)
           enddo
           close(50)
        end if
        !         output kinetic energy term also
        write(output_name,'(A4,A2,A1,I1,A1,I1,A1,A1,A4)') &
             "kin_", fn_type(i_function),"_",fn_species(i_function),"_", fn_n(i_function),"_",l_shell_str,".dat"
        open(50, file=output_name)
        write(50,*) "# ",n_grid(fn_species(i_function))
        do i_grid = 1, n_grid(fn_species(i_function)), 1
           write(50,*) r_grid(i_grid, fn_species(i_function)), basis_kinetic(i_grid, i_function)
        enddo
        close(50)
     enddo
  end if

  
  ! Well, n_basis is already set.  But resetting doesn't hurt.
  call get_bas_dimensions(n_basis_fns, fn_l(1:n_basis_fns), fn_species(1:n_basis_fns), max_basis_L, max_n_basis_fnLsp, n_basis)

  ! VB: BEFORE this point, anything to do with basis function manipulation is handled.
  !     AFTER this point, only splining, and manipulation and reorganizing of splines.
    
  !     can now allocate the actual basis storage arrays from module basis .
  call allocate_basis ( )

  !     and can store the splined versions of each basis function.
  !
  !     Note. We also "doctor" the splines, so that they are strictly zero outside
  !     a given radius, and so that they provide a smooth non-oscillating transition to zero.
  !     This is not identical with outer_radius further below, which is detemined by a threshold
  !     parameter, rather than a "strictly zero" criterion.
  !


  do i_function = 1, n_basis_fns, 1
     ! determine first "strictly zero" point on logarithmic grid
     outer_point = n_grid(fn_species(i_function))
     if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
       do while ( ( outer_point.ge.1 ) .and.(basis_wave(outer_point,i_function).eq.0.d0).and.(basis_wave_s(outer_point,i_function).eq.0.d0) )
          outer_point = outer_point - 1
       enddo
     else
       do while ( ( outer_point.ge.1 ) .and.(basis_wave(outer_point,i_function).eq.0.d0) )
          outer_point = outer_point - 1
       enddo
     endif
     if (outer_point.lt.n_grid(fn_species(i_function))) then
        outer_point = outer_point+1
     end if

     ! create radial function spline
     basis_wave_spl(1:n_max_spline,1:n_max_grid,i_function) = 0.d0
     call cubic_spline ( basis_wave(1,i_function), outer_point, basis_wave_spl(1,1,i_function) )
     ! fudge outermost segment of spline so that it goes to zero smoothly as
     ! a parabola (but with a discontinuous forst derivative at the first non-zero
     ! logarithmic grid point)
     basis_wave_spl(2,outer_point-1,i_function) = -2.d0 * basis_wave_spl(1,outer_point-1,i_function)
     basis_wave_spl(3,outer_point-1,i_function) = basis_wave_spl(1,outer_point-1,i_function)
     basis_wave_spl(4,outer_point-1,i_function) = 0.d0
     if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
       basis_wave_s_spl(1:n_max_spline,1:n_max_grid,i_function) = 0.d0
       call cubic_spline ( basis_wave_s(1,i_function), outer_point, basis_wave_s_spl(1,1,i_function) )
       basis_wave_s_spl(2,outer_point-1,i_function) = -2.d0 * basis_wave_s_spl(1,outer_point-1,i_function)
       basis_wave_s_spl(3,outer_point-1,i_function) = basis_wave_s_spl(1,outer_point-1,i_function)
       basis_wave_s_spl(4,outer_point-1,i_function) = 0.d0
     endif

     ! create radial kinetic energy part spline
     basis_kinetic_spl(1:n_max_spline,1:n_max_grid,i_function) = 0.d0
     call cubic_spline ( basis_kinetic(1,i_function), outer_point, basis_kinetic_spl(1,1,i_function) )
     ! fudge outermost segment of spline so that it goes to zero smoothly as
     ! a parabola (but with a discontinuous forst derivative at the first non-zero
     ! logarithmic grid point)
     basis_kinetic_spl(2,outer_point-1,i_function) = -2.d0 * basis_kinetic_spl(1,outer_point-1,i_function)
     basis_kinetic_spl(3,outer_point-1,i_function) = basis_kinetic_spl(1,outer_point-1,i_function)
     basis_kinetic_spl(4,outer_point-1,i_function) = 0.d0
     if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
       basis_kinetic_s_spl(1:n_max_spline,1:n_max_grid,i_function) = 0.d0
       call cubic_spline ( basis_kinetic_s(1,i_function), outer_point, basis_kinetic_s_spl(1,1,i_function) )
       basis_kinetic_s_spl(2,outer_point-1,i_function) = -2.d0 * basis_kinetic_s_spl(1,outer_point-1,i_function)
       basis_kinetic_s_spl(3,outer_point-1,i_function) = basis_kinetic_s_spl(1,outer_point-1,i_function)
       basis_kinetic_s_spl(4,outer_point-1,i_function) = 0.d0
     endif

     if (use_basis_gradients .or. flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) then
        ! create radial derivative spline
        basis_deriv_spl(1:n_max_spline,1:n_max_grid,i_function) = 0.d0
        call cubic_spline ( basis_deriv(1,i_function), outer_point, basis_deriv_spl(1,1,i_function) )
        ! fudge outermost segment of spline so that it goes to zero smoothly as
        ! a parabola (but with a discontinuous forst derivative at the first non-zero
        ! logarithmic grid point)
        basis_deriv_spl(2,outer_point-1,i_function) = -2.d0 * basis_deriv_spl(1,outer_point-1,i_function)
        basis_deriv_spl(3,outer_point-1,i_function) = basis_deriv_spl(1,outer_point-1,i_function)
        basis_deriv_spl(4,outer_point-1,i_function) = 0.d0
     end if

     if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
        ! create radial derivative spline
        basis_deriv_s_spl(1:n_max_spline,1:n_max_grid,i_function) = 0.d0
        call cubic_spline ( basis_deriv_s(1,i_function), outer_point, basis_deriv_s_spl(1,1,i_function) )
        basis_deriv_s_spl(2,outer_point-1,i_function) = -2.d0 * basis_deriv_s_spl(1,outer_point-1,i_function)
        basis_deriv_s_spl(3,outer_point-1,i_function) = basis_deriv_s_spl(1,outer_point-1,i_function)
        basis_deriv_s_spl(4,outer_point-1,i_function) = 0.d0
     endif
     
     if (flag_rel .eq. REL_atomic_zora) then
        ! create spline of kinetic energy part needed in scaled zora
        basis_kinetic_scaled_zora_spl (1:n_max_spline,1:n_max_grid,i_function) = 0.d0
        call cubic_spline ( basis_kinetic_scaled_zora(1,i_function), outer_point, basis_kinetic_scaled_zora_spl(1,1,i_function) )
        ! fudge outermost segment of spline so that it goes to zero smoothly as
        ! a parabola (but with a discontinuous forst derivative at the first non-zero
        ! logarithmic grid point)
        basis_kinetic_scaled_zora_spl(2,outer_point-1,i_function) = -2.d0 * basis_kinetic_scaled_zora_spl (1,outer_point-1,i_function)
        basis_kinetic_scaled_zora_spl(3,outer_point-1,i_function) = basis_kinetic_scaled_zora_spl(1,outer_point-1,i_function)
        basis_kinetic_scaled_zora_spl(4,outer_point-1,i_function)=0.d0
     end if
  enddo ! end loop over basis fn splines
  
  ! If requested, do a separate verification of the basis function behavior near the cutoff.
  if (force_smooth_cutoff) then
     do i_function = 1, n_basis_fns, 1
        ! determine first "strictly zero" point on logarithmic grid
        outer_point = n_grid(fn_species(i_function))
        do while ( ( outer_point.ge.1 ) .and.(basis_wave(outer_point,i_function).eq.0.d0) )
           outer_point = outer_point - 1
        enddo
        if (outer_point.lt.n_grid(fn_species(i_function))) then
           outer_point = outer_point+1
        end if
        ! radial function first
        if ( ( basis_wave(outer_point-2,i_function) .gt.smooth_cutoff_limit) .or. ( basis_wave(outer_point-1,i_function) .gt.smooth_cutoff_limit) ) then
           write(use_unit,'(1X,A)') "* After radial function orthonormalization / packing:"
           write(use_unit,'(1X,A,A,I5,A)') &
                "* Warning: You requested strict checking of the radial ", "function cutoff, and function # ", i_function, &
                " exceeds the threshold at its two outermost points."
           stop
        end if
        ! kinetic function next
        if ( ( basis_kinetic(outer_point-2,i_function) .gt.smooth_cutoff_limit) .or. &
             ( basis_kinetic(outer_point-1,i_function) .gt.smooth_cutoff_limit) ) then
           
           write(use_unit,'(1X,A)') "* After radial function orthonormalization / packing:"
           write(use_unit,'(1X,A,A,I5,A)') &
                "* Warning: You requested strict checking of the radial ", "function cutoff, and kinetic # ", i_function, &
                " exceeds the threshold at its two outermost points."
           stop
        end if
        ! radial derivative next
        if (use_basis_gradients) then
           if ( ( basis_deriv(outer_point-2,i_function) .gt.smooth_cutoff_limit) .or. &
                ( basis_deriv(outer_point-1,i_function) .gt.smooth_cutoff_limit) ) then
              write(use_unit,'(1X,A)') "* After radial function orthonormalization / packing:"
              write(use_unit,'(1X,A,A,I5,A)') &
                   "* Warning: You requested strict checking of the radial ", "function cutoff, and deriv # ", i_function, &
                   " exceeds the threshold at its two outermost points."
              stop
           end if
        end if
     enddo
  end if
  ! end cutoff verification.
  
  ! debugging
  !write(use_unit,*) "TESTING TESTING WE RESET L_WAVE_MAX and l_shell_max! "
  !write(use_unit,*) "Original value: ", l_wave_max
  ! Hopefully l_wave_max still works. 
  l_wave_max = maxval(fn_l(1:n_basis_fns))

  ! fix l_shell_max to correspond normal basis set
  ! and save the original l_shell_max to ext_l_shell_max
  ext_l_shell_max = l_shell_max (1:n_species)

  l_shell_max = basis_l_max (1:n_species)

  ! Generate basis function indexing
  basisfn_l = fn_l(1:n_basis_fns)
  basisfn_species = fn_species(1:n_basis_fns)
  basisfn_type = fn_type(1:n_basis_fns)
  basisfn_n = fn_n(1:n_basis_fns)
  basisfn_k = fn_k(1:n_basis_fns)
  ! (Rundong) The Lsp_ arrays here may need to be further modified for ralativistic cases. 
  ! Currently I don't know what they are used for.
  if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
     call generate_full_bas_rel(n_basis_fns, basisfn_n, basisfn_l, basisfn_k, basisfn_species, &
          max_basis_L, max_n_basis_fnLsp, n_basis, n_basis_small, basis_atom, basis_small_atom,&
          basis_l, basis_small_l, basis_k, basis_small_k, basis_m, basis_small_m, basis_fn, basis_small_fn, &
          Lsp2n_basis_fnLsp, Lsp2basis_fn, Lsp2basis_sp, atom2basis_off, sp2n_basis_sp)
  else
     call generate_full_bas(n_basis_fns, basisfn_l, basisfn_species, &
     &                      max_basis_L, max_n_basis_fnLsp, n_basis, &
     &                      basis_atom, basis_l, basis_m, basis_fn, &
     &                      Lsp2n_basis_fnLsp, Lsp2basis_fn, Lsp2basis_sp, &
     &                      atom2basis_off, sp2n_basis_sp)
  endif
  max_n_basis_sp = maxval(sp2n_basis_sp)

! ATOM BSSE:
!  create the basis_mapping here:
  if (calculate_atom_bsse) then

  if (current_atom_num_bsse==1) then
      write(info_str,'(2X,A)') "Creating the basis mapping for fresh BSSE geometry"
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      if (myid.eq.0) then
        write(use_unit,'(2X,A,I8)') "| Atom number for atom-based counterpoise correction: ", current_atom_num_bsse
      endif
      do i_basis = 1, n_basis, 1
        basis_mapping(i_basis) = i_basis
      enddo
  endif !loop for current_atom_num_bsse<1
  
  if (current_atom_num_bsse>1) then
      write(info_str,'(2X,A)') "Creating the basis mapping for fresh BSSE geometry"
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      if (myid.eq.0) then
        write(use_unit,'(2X,A,I8)') "| Atom number for atom-based counterpoise correction: ", current_atom_num_bsse
      endif

      j_basis=0
      basis_counter=0
      empty_basis=0
      
  ! for the 1st atoms previous to the _atom_num_bsse : just update the basis counter
   do i_atom = 2, current_atom_num_bsse, 1
     do i_function = species_first_function(species(i_atom)), species_last_function(species(i_atom))
        i_l = fn_l(i_function)
        do i_m = -i_l, i_l
           empty_basis = empty_basis+1
        end do
     end do
   enddo      

   do i_function = species_first_function(species(1)), species_last_function(species(1))
     i_l = fn_l(i_function)
     do i_m = -i_l, i_l
        empty_basis = empty_basis+1
        basis_counter= basis_counter+1
        basis_mapping(basis_counter) = empty_basis
     end do
   end do  
      
   do i_atom = 2, current_atom_num_bsse, 1
     do i_function = species_first_function(species(i_atom)), species_last_function(species(i_atom))
        i_l = fn_l(i_function)
        do i_m = -i_l, i_l
           j_basis = j_basis+1
           basis_counter= basis_counter+1
           basis_mapping(basis_counter) = j_basis
         end do
     end do
   enddo
    
   do i_atom = current_atom_num_bsse+1,n_atoms, 1
     do i_function = species_first_function(species(i_atom)), species_last_function(species(i_atom))
        i_l = fn_l(i_function)
        do i_m = -i_l, i_l
           empty_basis = empty_basis+1
           basis_counter= basis_counter+1
           basis_mapping(basis_counter) = empty_basis
        end do
     end do
   enddo
    
  endif ! for current_atom_num_bsse>1    

  endif !  end stuff for atom_bsse

  ! All basis functions are stored and indexed. Now verify consistency.
  
  ! determine the outer radius of each basis function u(r)
  ! [i.e. the radius outside of which all integrations may be skipped
  !  because u(r) is practically zero]
  ! MUST CHECK BOTH ACTUAL FN AND SECOND DERIVATIVE
  do i_function = 1, n_basis_fns, 1
     i_grid = n_grid(fn_species(i_function))
     if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
       do while ( (abs(basis_wave(i_grid,i_function)).le.wave_threshold) .and. &
                (abs(basis_wave_s(i_grid,i_function)).le.wave_threshold) .and. &
               (abs(basis_kinetic(i_grid,i_function)).le.wave_threshold) .and. &
             (abs(basis_kinetic_s(i_grid,i_function)).le.wave_threshold) .and. (i_grid.gt.1) )
          i_grid = i_grid-1
       enddo
     else
       do while ( (abs(basis_wave(i_grid,i_function)).le.wave_threshold) .and. &
            (abs(basis_kinetic(i_grid,i_function)).le.wave_threshold) .and.(i_grid.gt.1) )
          i_grid = i_grid-1
       enddo
     endif
     if (i_grid.le.1) then
        
        if (myid.eq.0) then
           write(use_unit,'(1X,A,A)') "* Warning - a basis function is ", "lower that the requested ", "threshold value for integrations everywhere."
           write(use_unit,'(1X,A,A)') "Species : ", trim(species_name(fn_species(i_function)))
           write(use_unit,'(1X,A,A)') "Type    : ", fn_type(i_function)
           write(use_unit,'(1X,A,I3,A,I3)') "(n,l)   : ", fn_n(i_function), ",", fn_l(i_function)
        end if
        
        stop
     end if
     outer_radius(i_function) = r_grid(i_grid,fn_species(i_function))
  enddo
  

  if (use_basis_gradients .or. flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) then
     !       also check first derivative
     do i_function = 1, n_basis_fns, 1
        i_grid = n_grid(fn_species(i_function))
        do while ( (abs(basis_deriv(i_grid,i_function)).le.wave_threshold) .and. (i_grid.gt.1) )
           i_grid = i_grid-1
        enddo
        if (i_grid.eq.n_grid(fn_species(i_function))) then
           if ( r_grid(i_grid,fn_species(i_function)) .gt.outer_radius(i_function) ) then
              outer_radius(i_function) = r_grid(i_grid,fn_species(i_function))
           end if
        else
           if ( r_grid(i_grid+1,fn_species(i_function)) .gt.outer_radius(i_function) ) then
              outer_radius(i_function) = r_grid(i_grid,fn_species(i_function))
           end if
        end if
     enddo
  end if
  

  
  
  !------------ end order of basis functions----------
  
  atom_radius_sq = 0.d0  
  do i_function = 1, n_basis_fns, 1
     outer_radius_sq(i_function) = outer_radius(i_function)**2
     atom_radius_sq(fn_species(i_function)) = max (atom_radius_sq(fn_species(i_function)), outer_radius_sq(i_function) )
  enddo

  ! Make sure that multipole_radius_free >= atom_radius := max(outer_radius)
  !
  ! (Rundong) multipole_radius_free_sq was generated in get_free_atoms.f90. But I think that the present scheme is not proper 
  ! for fully-relativistic cases. We may need to perform an inward integration from the outermost grid, and compare it to 
  ! the threshold. Currently I set multipole_radius_free to atom_radius. To me, this makes sense.
  multipole_radius_free_sq = max(multipole_radius_free_sq, atom_radius_sq)
  if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) multipole_radius_free_sq = atom_radius_sq
 !multipole_radius_free_sq = atom_radius_sq  !!!!! (Rundong) For tests between rel and nonrel calculations
  multipole_radius_free = sqrt(multipole_radius_free_sq)
  atom_radius = sqrt(atom_radius_sq)

  ! Construct arrays to reconstruct fn<->species relations without
  ! basisfn_species.
  basis_fn_atom = .false.
  do i_basis_2 = 1, n_basis, 1
     basis_fn_atom(basis_fn(i_basis_2),basis_atom(i_basis_2)) = .true.
  enddo
  i_offset_spl = 0
  species_offset = 0
  do i_species = 1, n_species, 1
     basis_fn_start_spl(i_species) = i_offset_spl + 1
     i_function2 = 0
     i_species_index = 0
     do i_basis_2 = 1, n_basis_fns, 1
        if (fn_species(i_basis_2).eq.i_species) then
           i_function2 = i_function2 + 1
           i_species_index(i_function2) = i_basis_2
        end if
     enddo
     n_function = i_function2
     n_basis_fn_species(i_species) = n_function
     do i_basis_2 = 1, n_function, 1
        temp_radius(i_basis_2) = outer_radius_sq(i_species_index(i_basis_2))
     enddo
     call insertionsort( temp_radius, n_function, perm_basis_fns, perm_basis_fns_inv )
     ! Index array linking the actual order of the spline array to the original order of radial functions
     ! for Hamiltonian evaluation
     do i_spline = 1, n_function, 1
        i_offset = i_spline + species_offset
        ! store index of spline function as a function of the radial fn index
        perm_basis_fns_spl(i_offset) = species_offset + perm_basis_fns_inv(i_spline)
        ! store the inverse, i.e. the index of the radial function as a
        ! function of the index used in the array of splined radial functions
        i_radial_fn(i_offset) = species_offset + perm_basis_fns(i_spline)
     enddo
     ! store species offset and increment it for the next species
     spline_offset(i_species) = species_offset
     species_offset = species_offset + n_function
     ! and store all basis function splines in arrays in order of
     ! increasing outer_radius per species
     ! in principle, this loop also runs over _spline_ functions, not over
     ! radial functions in their original order ...
     do i_basis_2 = 1, n_function, 1
        i_offset_spl = i_offset_spl + 1
        basis_wave_ordered(i_offset_spl,1:4,1:n_grid(i_species)) = &
             basis_wave_spl(1:4,1:n_grid(i_species), i_species_index(perm_basis_fns_inv(i_basis_2)))
        basis_kinetic_ordered(i_offset_spl,1:4,1:n_grid(i_species)) = &
             basis_kinetic_spl(1:4,1:n_grid(i_species), i_species_index(perm_basis_fns_inv(i_basis_2)))
        if(flag_rel==REL_atomic_zora)then
           basis_kinetic_scaled_zora_ordered (i_offset_spl,1:4,1:n_grid(i_species)) = &
                basis_kinetic_scaled_zora_spl (1:4,1:n_grid(i_species), i_species_index(perm_basis_fns_inv(i_basis_2)))
        end if
        if (use_basis_gradients .or. flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) then
           basis_deriv_ordered(i_offset_spl,1:4,1:n_grid(i_species)) = &
                basis_deriv_spl(1:4,1:n_grid(i_species), i_species_index(perm_basis_fns_inv(i_basis_2)))
        end if
        if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
           basis_wave_s_ordered(i_offset_spl,1:4,1:n_grid(i_species)) = &
                basis_wave_s_spl(1:4,1:n_grid(i_species), i_species_index(perm_basis_fns_inv(i_basis_2)))
           basis_kinetic_s_ordered(i_offset_spl,1:4,1:n_grid(i_species)) = &
                basis_kinetic_s_spl(1:4,1:n_grid(i_species), i_species_index(perm_basis_fns_inv(i_basis_2)))
           basis_deriv_s_ordered(i_offset_spl,1:4,1:n_grid(i_species)) = &
                basis_deriv_s_spl(1:4,1:n_grid(i_species), i_species_index(perm_basis_fns_inv(i_basis_2)))
        endif
     enddo
  enddo
  
  if (myid.eq.0) then
     write(use_unit,*)
     write(use_unit,'(2X,A)') "Basis size parameters after reduction:"
     write (use_unit,'(2X,A,I8)') "| Total number of radial functions: ", n_basis_fns
     if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
        write (use_unit,'(2X,A,I8)') "| Total number of scalar large component basis functions : ", n_basis
        write (use_unit,'(2X,A,I8)') "| Total number of scalar small component basis functions : ", n_basis_small
     else
        write (use_unit,'(2X,A,I8)') "| Total number of basis functions : ", n_basis
     endif
     write(use_unit,*)
  end if
  if ( allocated(basis_wave) ) then
     deallocate ( basis_wave )
  end if
  if ( allocated(basis_wave_s) ) then
     deallocate ( basis_wave_s )
  end if
  if ( allocated(basis_deriv) ) then
     deallocate ( basis_deriv )
  end if
  if ( allocated(basis_deriv_s) ) then
     deallocate ( basis_deriv_s )
  end if
  if ( allocated(basis_kinetic) ) then
     deallocate ( basis_kinetic )
  end if
  if ( allocated(basis_kinetic_s) ) then
     deallocate ( basis_kinetic_s )
  end if
  if(allocated ( basis_kinetic_scaled_zora))then
     deallocate ( basis_kinetic_scaled_zora)
  end if
  return
end subroutine shrink_fixed_basis_phi_thresh
!******
