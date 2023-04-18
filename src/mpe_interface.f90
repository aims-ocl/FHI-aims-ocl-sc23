!****h* FHI-aims/mpe_interface
!  NAME
!    mpe_interface - this module contains all public variables for the 
!                    MPE implicit solvent model and collects necessary
!                    information from aims modules
!  SYNOPSIS

module mpe_interface

!  PURPOSE
!    This module builds the interface between the MPE continuum solvation
!    module and FHI-aims.
!
!    For more information regarding the MPE implicit solvent model
!    please refer to
!     [1] M. Sinstein, H. Oberhofer, ..., V. Blum, and K. Reuter, in preparation
!     [2] D. Rinaldi, A. Bouchy, J.L. Rivail, and V. Dillet, 
!           J. Chem. Phys. 120, 2343 (2004)
!     [3] V. Dillet, D. Rinaldi, J.G. Angyan, and J.L. Rivail, 
!           Chem. Phys. Lett. 202, 18 (1993)
!     [4] J.L. Rivail and D. Rinaldi, Chem. Phys. 18, 233-242 (1976)
!
!    Short description of public module variables:
!
!     o isc_interfaces:
!       This is an array of a derived type which encapsulated all
!       necessary information for the coupling in the MPE module
!       (number of boundary conditions, dielectric constants, as well as
!       coordinates and coord.-system for all sampling points)
!
!  USES

   use types, only: dp
   use constants, only: bohr, pi4, hartree
   use dimensions, only: &
         l_pot_max, &
         n_atoms, &
         n_occ_atoms, &
         n_empty_atoms, &
         n_my_batches, &
         n_full_points
   use species_data, only: &
         species_z, &
         species_pseudoized
   use geometry, only: coords, species, empty
   use grids, only: batches
   use mpi_tasks, only: myid, use_mpi, mpi_comm_global, &
         f_stop => aims_stop, &
         aims_stop_coll
   use aims_memory_tracking, only: &
         f_allocate => aims_allocate, &
         f_deallocate => aims_deallocate
   use localorb_io, only: localorb_info, OL_norm
   use timing, only: get_timestamps, get_times, output_times
   use applicable_citations, only: cite_reference
   use xml_write, only: &
         Xml_file_t, xml_attr, tostr, &
         xml_open_file, xml_close_file, &
         xml_elem, xml_open, xml_close
   use runtime_choices, only: &
         use_scalapack, &
         mpe_xml_logfile, &
         mpe_xml_loglevel, &
         mpe_factorization_type, &
         mpe_solvent_permittivity, &
         mpe_bulk_permittivity, &
         mpe_n_centers_ep, &
         mpe_lmax_rf, &
         mpe_lmax_ep, &
         mpe_n_boundary_conditions, &
         mpe_degree_of_determination, &
         mpe_f_sparsity_threshold, &
         mpe_skip_first_n_scf_steps, &
         mpe_nonel_model, &
         mpe_nonel_model_bulk, &
         mpe_nonel_alpha, &
         mpe_nonel_beta, &
         mpe_nonel_alpha_if, &
         mpe_nonel_alpha_bulk, &
         mpe_nonel_beta_bulk, &
         mpe_tol_adjR2, &
         mpe_tol_adjR2_wait_scf, &
         isc_cavity_type, &
         isc_isodensity_value, &
         ifp, &
         ifp_dist, &
         ifp_normal, &
         restart_read, &
         isc_cavity_restart_read, &
         isc_cavity_restart_read_file, &
         mpe_lmax_ep_from_basis, &
         mpe_lmax_ep_additional_order, &
         mpe_lmax_ep_increase_cluster
   use isc_implicit_solvent_cavity, only: &
         isc_surface_area, &
         isc_cavity_volume, &
         plane_origin
   use mpe_dielectric_interfaces, only: &
         create_dielectric_interfaces, &
         update_dielectric_interfaces, &
         cleanup_dielectric_interfaces, &
         never_updated_before
   use mpe_constants, only: MPE_CONST, ISC_CONST
   use mpe_types, only: &
         initialize_mpe_types, cleanup_mpe_types, &
         DielectricInterface, InterfacePoint, &
         DielectricContinuum, &
         Basis, SolHarmBasis, &
         BasisCenter, SolHarmBasisCenter, &
         SpVecTuple, &
            SpVecTuple_vector_from_dense_vector
   use mpe_reaction_field, only: &
         MPESettings, &
         calculate_reaction_field_coefficients, &
         calculate_potential_at_points, &
         initialize_mpe_solver, &
         adjust_basis_functions, &
         cleanup_mpe_solver
   use mpe_dielectric_interfaces_common, only: &
         calculate_distance_to_plane, &
         determine_species_radius, &
         get_atoms_and_radii, &
         R_atom

   implicit none

!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180 (2009), 2175-2196.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   private  ! make everything private by default


   ! **************              <PUBLIC>              ************************
   ! ROUTINES (subroutines and functions)
   public   mpe_main
   public   free_mpe_interface
   public   mpe_check_initialization
   public   mpe_calculate_reaction_field
   public   mpe_calculate_interaction_energy_with_nuclei
   public   mpe_calculate_nonelectrostatic_energy
   public   mpe_store_energy_contributions
   public   mpe_read_energy_contributions
   public   mpe_get_continuum_name

   ! VARIABLES
   public   mpe_dc_indices
   public   mpe_dc_indices_atoms
   public   mpe_solved
   public   mpe_lmax_ep_nucl_dep_arr
   public   mpe_charge_per_dc

   ! CONSTANTS
   public   IND_DC

   ! *************              </PUBLIC>               ***********************


   ! PARAMETERS

   real(dp), parameter :: ZERO=0.e0_dp, ONE=1.e0_dp
   ! dielectric permittivity within the cavity
   !MS: Setting anything else than 1 here will require many changes
   real(dp), parameter :: EPSILON_IN = ONE 
   ! number of boundary conditions between continua
   integer, parameter :: mpe_n_boundary_conditions_2 = 2
   type :: ContinuumIndices
      ! Actual assignment will happen in define_dc_and_if_indices
      integer, allocatable :: NA(:), R(:), Q(:), O(:)
   end type
   type(ContinuumIndices) :: IND_DC
   type :: InterfaceIndices
      ! Actual assignment will happen in define_dc_and_if_indices
      integer, allocatable :: RQ(:), RO(:), QO(:), QO_ghost(:)
   end type
   type(InterfaceIndices) :: IND_IF

   !type(BasisCenter), parameter :: CENTER_CONST = BasisCenter( &
   !               coord=(/ZERO,ZERO,ZERO/), lmin=0, lmax=0 )


   ! VARIABLES

   type, public :: MPEEnergyContributions
      ! Meaning of abbreviations:
      !  V: reaction field, i.e. additional electrostatic potential
      !  rho: electron density
      !  nuc: nuclear charge density
      real(dp) :: rho_V
      real(dp) :: nuc_V
      real(dp) :: nonel
   end type
   type(MPEEnergyContributions), public :: mpe_energy_contributions

   type(MPESettings) :: settings

   ! dielectric continua and interfaces in the problem
   type(DielectricContinuum), allocatable :: isc_continua(:)
   type(DielectricInterface), allocatable :: isc_interfaces(:)

   ! partitioning of space into dielectric continuum regions 
   integer, allocatable, target :: mpe_dc_indices(:)
   integer, allocatable, target :: mpe_dc_indices_atoms(:)

   ! Note: module variables are initialized in set_mpe_defaults()
   type(Xml_file_t) :: mpe_xml_file
   integer, public :: scf_step_of_first_call

   logical, public :: module_initialized

   logical :: mpe_solved

   ! Array of nucleus dependent expansion orders
   integer, allocatable :: mpe_lmax_ep_nucl_dep_arr(:)

   real(dp), allocatable :: mpe_charge_per_dc(:)

contains



!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/v_hartree_callback
!  NAME
!    v_hartree_callback
!  SYNOPSIS

subroutine v_hartree_callback(n, points, v_h, v_h_gradient)

!  PURPOSE
!    This routine is a wrapper for the corresponding FHI-aims procedure
!    capable of evaluating the Hartree potential and, if needed, its 
!    gradient at arbitrary points in real space.
!    The definition here has to match the one defined in the
!    mpe_reaction_field module.
!
!  USES
   use v_hartree_multipole_evaluation, only: &
         get_v_hartree_multipole_and_gradient
   implicit none

!  ARGUMENTS
   integer, intent(in) :: n
   real(dp), intent(in) :: points(3,n)
   real(dp), intent(out) :: v_h(n)
   real(dp), intent(out), optional :: v_h_gradient(3,n)

!  INPUTS
!   o n -- number of points to evaluate
!   o points -- array of point coordinates
!  OUTPUT
!   o v_h -- Hartree potential at points
!   o v_h_gradient -- gradient of Hartree potential at points
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   ! Note that in the aims routine, the gradient is also optional
   call get_v_hartree_multipole_and_gradient(n, points, v_h, v_h_gradient)

end subroutine v_hartree_callback



!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/mpe_main
!  NAME
!    mpe_main
!  SYNOPSIS

subroutine mpe_main( cur_scf_step )

!  PURPOSE
!    When solvent effect calculations are enabled, this subroutine is called
!    in every scf step after the Hartree potential has been splined but
!    before the potential is summed up.
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: cur_scf_step

!  INPUTS
!   o cur_scf_step -- current step in the scf-loop
!  OUTPUT
!    writes to standard output
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'mpe_implicit_solvent_main_routine'
   character(*), parameter :: deffmt = '2X'
   integer, parameter :: defprio = OL_norm

   logical :: geometry_changed
   real(dp) :: cpu_calc, clock_calc

   integer, parameter :: RHS_cols = 1
   real(dp) :: adjR2(RHS_cols)

   character(132) :: info_str

   call cite_reference("MPE_implicit_solvation")

   ! memorize initial scf step
   if ( scf_step_of_first_call.lt.0 ) &
      scf_step_of_first_call = cur_scf_step

   ! make sure public variables get allocated to avoid segfaults
   if ( .not. allocated(IND_DC%R) ) & ! This one should be allocated in all possible settings
      call define_dc_and_if_indices()
   call allocate_dc_indices()

   ! check if calculation can be skipped
   ! this is useful for convergence in case the calculation is started
   ! from an initial guess instead of restarting from converged vacuum density
   if ( (cur_scf_step - scf_step_of_first_call) &
         .lt. mpe_skip_first_n_scf_steps ) &
      return

   ! 1. INITIALIZE
   if ( .not.module_initialized ) then
      call localorb_info(' ')
      call localorb_info('  Initializing MPE interface')
      call get_timestamps(cpu_calc, clock_calc)
      call initialize_mpe_interface()
      call get_times(cpu_calc, clock_calc)
      call output_times(deffmt, 'MPE Initialization', &
             cpu_calc, clock_calc, defprio)
      call localorb_info(' ')
   endif

   if ( (myid.eq.0) .and. (mpe_xml_loglevel.gt.MPE_CONST % XML_NOLOG) ) then
      call xml_open('scf', &
         file=mpe_xml_file, &
         attrs=(/xml_attr('step',cur_scf_step)/) )
   endif

   ! 2. OPTIMIZE CAVITY
   call localorb_info('  Updating MPE interface')
   call get_timestamps(cpu_calc, clock_calc)
   call update_mpe_interface(geometry_changed)
   call get_times(cpu_calc, clock_calc)
   call output_times(deffmt, 'MPE Update', &
         cpu_calc, clock_calc, defprio)


   ! 3. SOLVE MPE EQUATIONS
   call localorb_info('  Calculating reaction field coefficients')
   call get_timestamps(cpu_calc, clock_calc)
   ! call WRAPPER routine
   call calculate_reaction_field_coefficients( settings, &
         geometry_changed, isc_continua, isc_interfaces, &
         v_hartree_callback, adjR2 )

   if (minval(adjR2) .lt. 1.e0_dp-mpe_tol_adjR2) then
      if (.not. mpe_tol_adjR2_wait_scf) then
         write(info_str, '(A)') ""
         call localorb_info(info_str)
         write(info_str, '(A)') "*******************************************************************"
         call localorb_info(info_str)
         write(info_str, '(1X,A,F8.4)') '* Could not solve MPE equations in given basis. adjR2 = ', minval(adjR2)
         call localorb_info(info_str)
         write(info_str, '(1X,A)') '* adjR2 tends to improve during SCF procedure; setting'
         call localorb_info(info_str)
         write(info_str, '(A)') ""
         call localorb_info(info_str)
         write(info_str, '(1X,A)') "          mpe_tol_adjR2_wait_scf     .true."
         call localorb_info(info_str)
         write(info_str, '(A)') ""
         call localorb_info(info_str)
         write(info_str, '(1X,A)') '* in control.in might help if adjR2 is not too far from tolerance.'
         call localorb_info(info_str)
         write(info_str, '(1X,A)') '* Otherwise, try increasing'
         call localorb_info(info_str)
         write(info_str, '(A)') ""
         call localorb_info(info_str)
         write(info_str, '(1X,A)') "          mpe_lmax_rf"
         call localorb_info(info_str)
         write(info_str, '(A)') ""
         call localorb_info(info_str)
         write(info_str, '(1X,A)') "* See user's guide for further info. Aims will abort now."
         call localorb_info(info_str)
         write(info_str, '(A)') "*******************************************************************"
         call localorb_info(info_str)
         write(info_str, '(A)') ""
         call localorb_info(info_str)
         call aims_stop_coll("Error (MPE): Adjusted R^2 below 1-tolerance. Could not solve MPE equations.")
      else
         mpe_solved = .false.
      endif
   else
      mpe_solved = .true.
   endif
   call get_times(cpu_calc, clock_calc)
   call output_times(deffmt, 'MPE Solver', &
         cpu_calc, clock_calc, defprio)

   if ( (myid.eq.0) .and. (mpe_xml_loglevel.ge.MPE_CONST % XML_DETAILED) ) then
      call write_xml_continua(mpe_xml_file, isc_continua)
   endif

   ! DONE

   ! set solvent effect as initialized
   module_initialized = .true.

   if ( (myid.eq.0) .and. (mpe_xml_loglevel.gt.MPE_CONST % XML_NOLOG) ) then
      call xml_close(file=mpe_xml_file)
   endif

end subroutine mpe_main

!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/define_dc_and_if_indices
!  NAME
!    define_dc_and_if_indices
!  SYNOPSIS

subroutine define_dc_and_if_indices()

!  PURPOSE
!    This subroutine define the indices of electrostatic regions and interfaces
!
!  USES
   implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2019).
!  SOURCE

    if ( ifp ) then
      ! Continuum indices
      allocate(IND_DC%NA(1))
      IND_DC%NA(1) = 0 ! undefined
      allocate(IND_DC%R(1))
      IND_DC%R(1)  = 1
      allocate(IND_DC%Q(1))
      IND_DC%Q(1)  = 2
      allocate(IND_DC%O(1))
      IND_DC%O(1)  = 3
      ! Interface indices
      allocate(IND_IF%RQ(1))
      IND_IF%RQ(1) = 1 ! solute (cavity) - liquid
      allocate(IND_IF%RO(1))
      IND_IF%RO(1) = 2 ! solute (cavity) - solid
      allocate(IND_IF%QO(1))
      IND_IF%QO(1) = 3 ! solid - liquid
      allocate(IND_IF%QO_ghost(1))
      IND_IF%QO_ghost (1)= 4 ! solid - liquid inside cavity
                      ! not included in MPE equations, only needed for interfacial energy term
      ! Only used in sum_up_whole_potential, but needs to be allocated here
      allocate(mpe_charge_per_dc(IND_DC%NA(1):IND_DC%O(1)))
   else
      ! Continuum indices
      allocate(IND_DC%NA(1))
      IND_DC%NA(1) = 0
      allocate(IND_DC%R(1))
      IND_DC%R(1)  = 1
      allocate(IND_DC%Q(1))
      IND_DC%Q(1)  = 2
      ! Interface indices
      allocate(IND_IF%RQ(1))
      IND_IF%RQ(1) = 1
      ! Only used in sum_up_whole_potential, but needs to be allocated here
      allocate(mpe_charge_per_dc(IND_DC%NA(1):IND_DC%Q(1)))
   endif ! ifp

end subroutine define_dc_and_if_indices


!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/allocate_dc_indices
!  NAME
!    allocate_dc_indices
!  SYNOPSIS

subroutine allocate_dc_indices()

!  PURPOSE
!    This subroutine sets the initial values for the global indices arrays
!    assigning integration grid points to dielectric continua.
!
!  USES
   implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   if (.not.allocated(mpe_dc_indices_atoms)) then
      call f_allocate(mpe_dc_indices_atoms, n_atoms, &
            name="Dielectric Indices at nuclear positions")
      mpe_dc_indices_atoms = IND_DC % NA(1)
   endif

   if (.not.allocated(mpe_dc_indices)) then
      call f_allocate(mpe_dc_indices, n_full_points, &
            name="+Dielectric Indices on integration grid")
      mpe_dc_indices = IND_DC % NA(1)
   endif

end subroutine allocate_dc_indices
!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/initialize_dc_indices
!  NAME
!    initialize_dc_indices
!  SYNOPSIS

subroutine initialize_dc_indices()

!  PURPOSE
!    This subroutine sets the initial values for the global indices arrays
!    assigning integration grid points to dielectric continua.
!
!  USES
   implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   real(dp), allocatable :: R_atom_sq(:)
   real(dp) :: coords(3), rho_dummy
   integer :: i_atom, i_species, i_full_points, i_batch, i_index

   mpe_dc_indices_atoms = IND_DC % R(1)
   !MS: This is not necessarily correct! However, if we allowed
   ! to have atoms in regions with a dielectric permittivity
   ! not equal to 1, sadly, we would run into other troubles as well...

   mpe_dc_indices = IND_DC % NA(1)
   if ( (isc_cavity_type.eq.ISC_CONST % CAVITY_OvlpSph) ) then
      call localorb_info("    initialize integration grid assignment")
      allocate(R_atom_sq(n_atoms))

      do i_atom = 1, n_atoms
         i_species = species(i_atom)
         if (empty(i_atom).or.species_pseudoized(i_species)) then
            R_atom_sq(i_atom) = 0
         else
            R_atom_sq(i_atom) = determine_species_radius(i_species, &
                              isc_isodensity_value)**2
         endif
      enddo

      i_full_points = 0
      do i_batch = 1, n_my_batches
         do i_index = 1, batches(i_batch)%size
            i_full_points = i_full_points + 1
            coords(:) = batches(i_batch) % points(i_index) % coords(:)
            mpe_dc_indices(i_full_points) = get_continuum_index_radii( &
               coords, radii_sq=R_atom_sq)
         enddo ! i_index
      enddo ! i_batch
   else
      call localorb_info("    pass integration grid assignment")
   endif

end subroutine initialize_dc_indices
!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/deallocate_dc_indices
!  NAME
!    deallocate_dc_indices
!  SYNOPSIS

subroutine deallocate_dc_indices()

!  PURPOSE
!    Free public module variables
!
!  USES
   implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   if (allocated(mpe_dc_indices)) &
      call f_deallocate(mpe_dc_indices, name="mpe_dc_indices")
   if (allocated(mpe_dc_indices_atoms)) &
      call f_deallocate(mpe_dc_indices_atoms, name="mpe_dc_indices_atoms")

end subroutine deallocate_dc_indices
!******
!----------------------------------------------------------------------------------------------
!****s* mpe_interface/mpe_check_initialization
!  NAME
!    mpe_check_initialization
!  SYNOPSIS

function mpe_check_initialization()

!  PURPOSE
!    This function return whether MPE has been initialized or not.
!
!  USES
   implicit none
!  ARGUMENTS
   logical :: mpe_check_initialization
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   mpe_check_initialization = module_initialized

end function mpe_check_initialization
!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/free_mpe_interface
!  NAME
!    free_mpe_interface
!  SYNOPSIS

subroutine free_mpe_interface()

!  PURPOSE
!    This subroutine sets all variables for the MPE interface.
!
!  USES
   implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   integer :: i_center, i_dc

   ! interfaces
   if (allocated(isc_interfaces)) &
      deallocate(isc_interfaces)

   ! continua with basis and coefficients
   if (allocated(isc_continua)) then
      do i_dc = lbound(isc_continua,1), ubound(isc_continua,1)
         select type (b => isc_continua(i_dc)%basis)
         class is (SolHarmBasis)
            if (allocated(b%centers)) then
               do i_center = lbound(b%centers,1),ubound(b%centers,1)
                  if (allocated(b%centers(i_center)%coeff))&
                     deallocate(b%centers(i_center)%coeff)
               enddo ! i_center
               deallocate(b%centers)
            endif
         end select
      enddo ! i_dc
      deallocate(isc_continua)
   endif

   call deallocate_dc_indices()

   ! deallocate implicit solvent cavity arrays
   call cleanup_dielectric_interfaces()

   ! clear derived MPI data types
   call cleanup_mpe_types()

   ! clear MPE solver
   call cleanup_mpe_solver(settings)

   ! close XML logging
   if ( (myid.eq.0) .and. (mpe_xml_loglevel.gt.MPE_CONST % XML_NOLOG) ) then
      call xml_close_file(file=mpe_xml_file)
   endif

end subroutine free_mpe_interface
!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/initialize_mpe_interface
!  NAME
!    initialize_mpe_interface
!  SYNOPSIS

subroutine initialize_mpe_interface()

!  PURPOSE
!    This subroutine sets all variables for the MPE interface.
!
!  USES
   implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(*), parameter :: func = 'initialize_mpe_interface'
   character(*), parameter :: deffmt = '4X'
   integer, parameter :: defprio = OL_norm

   real(dp) :: cpu_calc, clock_calc
   integer :: n_variables, n_centers_rf, i_dc, i_if, n_interfaces
   integer, allocatable :: n_requested(:)
   logical :: cavity_restart_file_exists

   character(132) :: info_str

   integer :: basis_size_isc, basis_size_ifp

   cavity_restart_file_exists = .false.
   if ( isc_cavity_restart_read ) then
      inquire(file=trim(isc_cavity_restart_read_file)//'.bin', exist=cavity_restart_file_exists)
   endif

   if (isc_cavity_type .eq. ISC_CONST%CAVITY_RhoMPStat .and. &
           .not. restart_read .and. &
           mpe_skip_first_n_scf_steps .le. 0 .and. &
           .not. cavity_restart_file_exists) then
      write(info_str, '(A)') ' * MPE: isc_cavity_type rho_multipole_static was specified but no meaningful density exists.'
      call localorb_info(info_str)
      write(info_str, '(A)') ' *      You have the following options:'
      call localorb_info(info_str)
      write(info_str, '(A)') ' *      1. Use overlapping_spheres (not recommended, only for debugging) or rho_free.'
      call localorb_info(info_str)
      write(info_str, '(A)') ' *      2. Use rho_multipole_dynamic.'
      call localorb_info(info_str)
      write(info_str, '(A)') ' *      3. Specify mpe_skip_first_n_scf_steps.'
      call localorb_info(info_str)
      write(info_str, '(A)') ' *      4. Specify and provide a suitable restart file, e.g. from a converged vacuum calculation.'
      call localorb_info(info_str)
      write(info_str, '(A)') ' *      5. Specify and provide a suitable cavity restart file in both .xyz and .bin format.'
      call localorb_info(info_str)
      write(info_str, '(A)') " *      For further information see the user's guide. Aims will abort now."
      call localorb_info(info_str)
      call f_stop('MPE: isc_cavity_type rho_multipole_static was specified but no meaningful density exists.')
   endif

   ! init XML logging
   if ( (myid.eq.0) .and. (mpe_xml_loglevel.gt.MPE_CONST % XML_NOLOG) ) then
      call xml_open_file( &
         filename=mpe_xml_logfile, &
         root="root", &
         file=mpe_xml_file, &
         tab_spaces=2 )
   endif

   call set_default_parameters()
   call set_mpe_settings(settings)
   call initialize_mpe_types(use_mpi)

   if (mpe_lmax_ep_from_basis) then
      call get_atoms_and_radii(mpe_lmax_ep_nucl_dep_arr)
      if (any(mpe_lmax_ep_nucl_dep_arr.lt.0)) &
         call f_stop('MPE: Expansion order cannot be reduced below 0.'//&
                     ' Choose higher mpe_lmax_ep_additional_order.')
   else
      call get_atoms_and_radii()
   endif

   call define_problem(isc_continua, isc_interfaces, basis_size_isc, basis_size_ifp)

   ! create cavity with desired number of sampling points
   call get_timestamps(cpu_calc, clock_calc)

   n_interfaces = size(isc_interfaces)
   allocate(n_requested(n_interfaces))
   n_requested(1) = int(mpe_degree_of_determination * basis_size_isc / &
             isc_interfaces(1)%n_bc )
     if (ifp) then
           n_requested(2) = 0
           n_requested(3) = int(mpe_degree_of_determination * &
                  basis_size_ifp / isc_interfaces(3)%n_bc )
           n_requested(4) = 0
   endif

   call create_dielectric_interfaces(n_requested, isc_interfaces)
   deallocate(n_requested)

   call get_times(cpu_calc, clock_calc)
   call output_times(deffmt, 'Cavity creation', &
          cpu_calc, clock_calc, defprio)

   call initialize_dc_indices()

   call get_timestamps(cpu_calc, clock_calc)
   !call adjust_basis_functions(isc_continua, isc_interfaces)
   ! Does not work here for interface case, moved to update_mpe_interface
   call get_times(cpu_calc, clock_calc)
   call output_times(deffmt, 'Basis update', &
          cpu_calc, clock_calc, defprio)

   call initialize_mpe_solver(settings, mpi_comm_global)

end subroutine initialize_mpe_interface
!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/update_mpe_interface
!  NAME
!    update_mpe_interface
!  SYNOPSIS

subroutine update_mpe_interface(geometry_changed)

!  PURPOSE
!    This subroutine sets all variables for the MPE interface.
!
!  USES
   implicit none
!  ARGUMENTS
   logical, intent(out) :: geometry_changed
!  INPUTS
!    none
!  OUTPUT
!   o geometry_changed -- did the geometry of the problem change during the update?
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(*), parameter :: deffmt = '4X'
   integer, parameter :: defprio = OL_norm
   real(dp) :: cpu_calc, clock_calc
   logical :: cavity_has_changed

   cavity_has_changed = .false.

   call get_timestamps(cpu_calc, clock_calc)
   if ( isc_cavity_type .eq. ISC_CONST%CAVITY_RhoMPDyn .or. &
           never_updated_before ) then
      if (use_mpi) then
         call update_dielectric_interfaces(isc_interfaces, cavity_has_changed, &
              mpi_comm=mpi_comm_global)
      else
         call update_dielectric_interfaces(isc_interfaces, cavity_has_changed)
      endif
   endif

   ! Moved here
   call adjust_basis_functions(isc_continua, isc_interfaces)

   call get_times(cpu_calc, clock_calc)
   call output_times(deffmt, 'Cavity update', &
          cpu_calc, clock_calc, defprio)

   if (cavity_has_changed) then
      ! set dielectric indices undefined in oder to force refresh
      call localorb_info("    set integration grid assignment for re-evaluation")
      mpe_dc_indices = IND_DC % NA(1)
   endif

   if (cavity_has_changed.and.module_initialized) then
      call get_timestamps(cpu_calc, clock_calc)
      ! update basis scaling
      !call adjust_basis_functions(isc_continua, isc_interfaces)
      call get_times(cpu_calc, clock_calc)
      call output_times(deffmt, 'Basis update', &
             cpu_calc, clock_calc, defprio)
   endif

   geometry_changed = cavity_has_changed

end subroutine update_mpe_interface



!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/set_default_parameters
!  NAME
!    set_default_parameters
!  SYNOPSIS

subroutine set_default_parameters()

!  PURPOSE
!    This subroutine checks whether MPE parameters have been specified in the
!    control.in and assigns default values if this is not the case.
!
!  USES
   implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   if (mpe_factorization_type.eq.MPE_CONST % FACTZN_UNDEF) then
      mpe_factorization_type = MPE_CONST % FACTZN_QRpSVD
   endif

   if (mpe_skip_first_n_scf_steps.lt.0) &
      mpe_skip_first_n_scf_steps = 0

   if (mpe_degree_of_determination.lt.ZERO) &
      mpe_degree_of_determination = 5.e0_dp

   if (mpe_lmax_rf.lt.0) mpe_lmax_rf = 8

   if (mpe_lmax_ep.ne.-1 .and. &
           (mpe_lmax_ep_from_basis .or. &
            mpe_lmax_ep_additional_order.ne.0) ) then
      call f_stop('Cannot define both fixed and basis set dependent mpe_lmax_ep')
   elseif (mpe_lmax_ep_additional_order.ne.0 .and. &
           .not. mpe_lmax_ep_from_basis) then
      call f_stop('Can only use mpe_lmax_ep_additional_order with mpe_lmax_ep_from_basis')
   elseif (mpe_lmax_ep.lt.0 .and. &
           .not. mpe_lmax_ep_from_basis) then
      mpe_lmax_ep = l_pot_max
   endif

   if (mpe_n_centers_ep.lt.0) mpe_n_centers_ep = n_occ_atoms+n_empty_atoms

   if (mpe_f_sparsity_threshold.lt.ZERO) &
      mpe_f_sparsity_threshold = ZERO

   ! Not really a parameter, but could not find a better place for this
   mpe_solved = .false.

end subroutine set_default_parameters


!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/set_mpe_settings
!  NAME
!    set_mpe_settings
!  SYNOPSIS

subroutine set_mpe_settings( settings )

!  PURPOSE
!    This subroutine applies runtime choices to MPE settings.
!
!  USES
   implicit none
!  ARGUMENTS
   type(MPESettings), intent(inout) :: settings
!  INPUTS
!   o settings -- contains MPE specific settings and (most) variables that
!                 persist between calls to the MPE solver
!  OUTPUT
!   o settings -- contains MPE specific settings and (most) variables that
!                 persist between calls to the MPE solver
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   settings % use_scalapack = use_scalapack
   settings % factorization_type = mpe_factorization_type
   settings % sparsity_threshold = mpe_f_sparsity_threshold

end subroutine set_mpe_settings



!******
!----------------------------------------------------------------------------------------------
!****s* mpe_interface/define_problem
!  NAME
!    define_problem
!  SYNOPSIS

subroutine define_problem( continua, interfaces, basis_size_isc, basis_size_ifp )

!  PURPOSE
!    This subroutine sets all variables for the MPE interface.
!
!  USES
   implicit none
!  ARGUMENTS
   type(DielectricContinuum), allocatable, intent(out) :: continua(:)
   type(DielectricInterface), allocatable, intent(out) :: interfaces(:)

   integer, intent(out) :: basis_size_isc, basis_size_ifp
!  INPUTS
!    none
!  OUTPUT
!   o continua -- dielectric continua with permittivities and basis sets
!   o interfaces -- dielectric interfaces with boundary conditions and
!                   corresponding continua
!   o basis_size_isc -- max number of basis functions to be evaluated at a
!                       cavity point
!   o basis_size_ifp -- max number of basis functions to be evaluated at an
!                       interface point
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   type(SolHarmBasis) :: basis_rf, basis_ep, basis_ep_q, basis_ep_o
   integer :: n_continua, n_interfaces, i_center, i_center_mirrored, n_centers_ep

   real(dp) :: dist_to_plane, center_to_center(3), mirrored_coords(3)

   logical, allocatable :: point_in_Q(:), point_in_O(:), neighbors_Q(:,:), neighbors_O(:,:)
   logical, allocatable :: clusters_Q(:,:), clusters_O(:,:)

   real(dp), allocatable :: R_half_sq(:)

   character(132) :: info_str, format_str

   ! 1. DIMENSIONS
   ! homogeneous embedding vs. embedding at dielectric interface?
   if ( ifp ) then
      n_continua = 3
      n_interfaces = 4
   else
      n_continua = 2
      n_interfaces = 1
   endif ! ifp

   basis_size_isc = 0
   basis_size_ifp = 0

   allocate(continua(n_continua))
   allocate(interfaces(n_interfaces))

   ! 2. BASIS SETS

   ! WORD OF CAUTION:
   ! The basis sets seem to give quite some flexibility.
   ! However, changes here MIGHT backfire because many things
   ! are not as flexible in the mpe_reaction_field module. Check!

   ! reaction field
   basis_rf%name = "reaction_field_R"
   basis_rf%solharm_type = MPE_CONST%BASIS_REG
   allocate(basis_rf%centers(1))
   basis_rf%centers(1)%coord = sum(coords,dim=2)/size(coords,dim=2)
   basis_rf%centers(1)%lmin = 0
   basis_rf%centers(1)%lmax = mpe_lmax_rf

   basis_size_isc = basis_size_isc + basis_rf%get_size()

   ! 3. DIELECTRIC CONTINUA
   ! first continuum: QM region inside of cavity
   continua(IND_DC%R(1))%eps = EPSILON_IN
   allocate(continua(IND_DC%R(1))%basis, source=basis_rf)

   ! external potential

   if (ifp) then
      allocate(R_half_sq(size(R_atom)))
      R_half_sq = R_atom*0.5e0_dp
      R_half_sq = R_half_sq*R_half_sq

      basis_ep%solharm_type = MPE_CONST%BASIS_IRR
      n_centers_ep = 2 * mpe_n_centers_ep
      allocate(basis_ep%centers(n_centers_ep))

      ! Centers on nuclei
      do i_center = 1, mpe_n_centers_ep
         basis_ep%centers(i_center)%coord = coords(:,i_center)
         basis_ep%centers(i_center)%lmin = 0
         if (mpe_lmax_ep_from_basis) then
            basis_ep%centers(i_center)%lmax = mpe_lmax_ep_nucl_dep_arr(i_center)
         else
            basis_ep%centers(i_center)%lmax = mpe_lmax_ep
         endif
      enddo

      allocate(point_in_Q(mpe_n_centers_ep))
      allocate(point_in_O(mpe_n_centers_ep))

      ! Image charges
      do i_center_mirrored = 1, mpe_n_centers_ep
         dist_to_plane = calculate_distance_to_plane(point=coords(:,i_center_mirrored), &
                 plane_normal = ifp_normal, plane_dist = ifp_dist)

         mirrored_coords = coords(:,i_center_mirrored) - 2 * ifp_normal * dist_to_plane

         ! Will be needed later to know to which basis to add this center
         point_in_Q(i_center_mirrored) = dist_to_plane .ge. 0
         point_in_O(i_center_mirrored) = .not. point_in_Q(i_center_mirrored)

         check_redundant_centers: do i_center = 1, mpe_n_centers_ep
            center_to_center = mirrored_coords - coords(:,i_center)
            if (sum(center_to_center * center_to_center) .lt. R_half_sq(i_center)) then
               point_in_Q(i_center_mirrored) = .false.
               point_in_O(i_center_mirrored) = .false.
               exit check_redundant_centers
            endif
         enddo check_redundant_centers

         basis_ep%centers(i_center_mirrored+mpe_n_centers_ep)%coord = mirrored_coords

         basis_ep%centers(i_center_mirrored+mpe_n_centers_ep)%lmin = 0
         if (mpe_lmax_ep_from_basis) then
            basis_ep%centers(i_center_mirrored+mpe_n_centers_ep)%lmax = &
                    mpe_lmax_ep_nucl_dep_arr(i_center_mirrored)
         else
            basis_ep%centers(i_center_mirrored+mpe_n_centers_ep)%lmax = mpe_lmax_ep
         endif
      enddo

      deallocate(R_half_sq)

      allocate(basis_ep_q%centers(mpe_n_centers_ep+count(point_in_Q)))
      basis_ep_q%centers(1:mpe_n_centers_ep) = basis_ep%centers(1:mpe_n_centers_ep)
      i_center_mirrored = 1
      do i_center = 1, mpe_n_centers_ep
         if (point_in_Q(i_center)) then
            basis_ep_q%centers(mpe_n_centers_ep + i_center_mirrored) = &
                     basis_ep%centers(mpe_n_centers_ep + i_center)
            i_center_mirrored = i_center_mirrored + 1
         endif
      enddo

      ! Center clustering, Q
      call find_neighboring_centers_ifp(basis_ep_q%centers, -1.e0_dp, neighbors_Q)
      call cluster_neighboring_centers(neighbors_Q, clusters_Q)
      deallocate(neighbors_Q)
      call reduce_clusters(clusters_Q)
      call reduce_centers(basis_ep_q%centers, clusters_Q)
      deallocate(clusters_Q)

      basis_ep_q%name = "external_potential_Q"
      basis_ep_q%solharm_type = basis_ep%solharm_type
      allocate(continua(IND_DC%Q(1))%basis, source=basis_ep_q)
      continua(IND_DC%Q(1))%eps = mpe_solvent_permittivity

      basis_size_ifp = basis_size_ifp + basis_ep_q%get_size()

      allocate(basis_ep_o%centers(mpe_n_centers_ep+count(point_in_O)))
      basis_ep_o%centers(1:mpe_n_centers_ep) = basis_ep%centers(1:mpe_n_centers_ep)
      i_center_mirrored = 1
      do i_center = 1, mpe_n_centers_ep
         if (point_in_O(i_center)) then
            basis_ep_o%centers(mpe_n_centers_ep + i_center_mirrored) = &
                     basis_ep%centers(mpe_n_centers_ep + i_center)
            i_center_mirrored = i_center_mirrored + 1
         endif
      enddo

      ! Center clustering, O
      call find_neighboring_centers_ifp(basis_ep_o%centers, +1.e0_dp, neighbors_O)
      call cluster_neighboring_centers(neighbors_O, clusters_O)
      deallocate(neighbors_O)
      call reduce_clusters(clusters_O)
      call reduce_centers(basis_ep_o%centers, clusters_O)
      deallocate(clusters_O)

      basis_ep_o%name = "external_potential_O"
      basis_ep_o%solharm_type = basis_ep%solharm_type
      allocate(continua(IND_DC%O(1))%basis, source=basis_ep_o)
      continua(IND_DC%O(1))%eps = mpe_bulk_permittivity

      basis_size_ifp = basis_size_ifp + basis_ep_o%get_size()
      basis_size_isc = basis_size_isc + max(basis_ep_q%get_size(), basis_ep_o%get_size())

      write(info_str, '(4X,A,A,A,I5,A)') 'Defined basis ', trim(continua(IND_DC%R(1))%basis%name), &
              ' with ', continua(IND_DC%R(1))%basis%get_size(), ' functions'
      call localorb_info(info_str)
      write(info_str, '(4X,A,A,A,I5,A)') 'Defined basis ', trim(continua(IND_DC%Q(1))%basis%name), &
              ' with ', continua(IND_DC%Q(1))%basis%get_size(), ' functions'
      call localorb_info(info_str)
      write(info_str, '(4X,A,A,A,I5,A)') 'Defined basis ', trim(continua(IND_DC%O(1))%basis%name), &
              ' with ', continua(IND_DC%O(1))%basis%get_size(), ' functions'
      call localorb_info(info_str)

   else ! homogeneous case
      basis_ep%name = "external_potential_Q"
      basis_ep%solharm_type = MPE_CONST%BASIS_IRR
      allocate(basis_ep%centers(mpe_n_centers_ep))
      do i_center = 1, mpe_n_centers_ep
         basis_ep%centers(i_center)%coord = coords(:,i_center)
         basis_ep%centers(i_center)%lmin = 0
         if (mpe_lmax_ep_from_basis) then
            basis_ep%centers(i_center)%lmax = mpe_lmax_ep_nucl_dep_arr(i_center)
         else
            basis_ep%centers(i_center)%lmax = mpe_lmax_ep
         endif
      enddo

      ! second continuum: solvent region
      continua(IND_DC%Q(1))%eps = mpe_solvent_permittivity
      allocate(continua(IND_DC%Q(1))%basis, source=basis_ep)

      basis_size_isc = basis_size_isc + basis_ep%get_size()

   endif ! ifp

   ! 4. DIELECTRIC INTERFACES
   !
   ! REMINDER: boundary conditions
   !  1. continuity of the potential
   !  2. continuity of electric flux normal to cavity surface
   !  3. continuity of the electric field tangential to the cavity surface
   !  4.   (2nd tangential cartesian direction)

   ! first interface: QM region with solvent
   interfaces(IND_IF%RQ(1))%n_bc = mpe_n_boundary_conditions
   interfaces(IND_IF%RQ(1))%dc_ind_pos = IND_DC%R(1)
   interfaces(IND_IF%RQ(1))%dc_ind_neg = IND_DC%Q(1)
   interfaces(IND_IF%RQ(1))%surface_energy = mpe_nonel_alpha
   interfaces(IND_IF%RQ(1))%volume_energy = mpe_nonel_beta

   if ( n_interfaces .eq. 4 ) then
      ! second interface: QM region with other continuum
      interfaces(IND_IF%RO(1))%n_bc = mpe_n_boundary_conditions
      interfaces(IND_IF%RO(1))%dc_ind_pos = IND_DC%R(1)
      interfaces(IND_IF%RO(1))%dc_ind_neg = IND_DC%O(1)
      interfaces(IND_IF%RO(1))%surface_energy = mpe_nonel_alpha_bulk
      interfaces(IND_IF%RO(1))%volume_energy = mpe_nonel_beta_bulk

      ! third interface: between the two continua
      interfaces(IND_IF%QO(1))%n_bc = mpe_n_boundary_conditions_2
      interfaces(IND_IF%QO(1))%dc_ind_pos = IND_DC%Q(1)
      interfaces(IND_IF%QO(1))%dc_ind_neg = IND_DC%O(1)
      interfaces(IND_IF%QO(1))%surface_energy = 0.0e0_dp
      interfaces(IND_IF%QO(1))%volume_energy = 0.0e0_dp

      ! fourth interface: between the two continua, but inside the cavity
      interfaces(IND_IF%QO_ghost(1))%n_bc = 0
      interfaces(IND_IF%QO_ghost(1))%dc_ind_pos = IND_DC%Q(1)
      interfaces(IND_IF%QO_ghost(1))%dc_ind_neg = IND_DC%O(1)
      if (mpe_nonel_alpha_if .eq. -1.0e0_dp) then
              interfaces(IND_IF%QO_ghost(1))%surface_energy = -abs(mpe_nonel_alpha-mpe_nonel_alpha_bulk)
      else
              interfaces(IND_IF%QO_ghost(1))%surface_energy = -mpe_nonel_alpha_if
      endif
      interfaces(IND_IF%QO_ghost(1))%volume_energy =       mpe_nonel_beta -mpe_nonel_beta_bulk
      ! normal vector facing towards solvent (continuum 1)
   endif

   ! deallocate
   deallocate(basis_rf%centers)
   deallocate(basis_ep%centers)
   if (ifp) then
      deallocate(basis_ep_q%centers)
      deallocate(basis_ep_o%centers)
      deallocate(point_in_Q)
      deallocate(point_in_O)
   endif

   if (allocated(mpe_lmax_ep_nucl_dep_arr)) deallocate(mpe_lmax_ep_nucl_dep_arr)


   contains

   subroutine find_neighboring_centers_ifp(centers, sgn, neighbors)

      implicit none

      type(SolHarmBasisCenter), intent(in) :: centers(:)
      real(dp), intent(in) :: sgn

      logical, allocatable, intent(out) :: neighbors(:,:)

      integer :: i_center, j_center, n_centers
      real(dp), allocatable :: dist_to_plane(:)

      n_centers = size(centers)

      allocate(neighbors(n_centers, n_centers))
      allocate(dist_to_plane(n_centers))

      do i_center = 1, n_centers
         dist_to_plane(i_center) = sgn * calculate_distance_to_plane(point=centers(i_center)%coord, &
                 plane_normal = ifp_normal, plane_dist = ifp_dist)
      enddo

      do i_center = 1, n_centers
         neighbors(i_center, i_center) = .true.
         if ( dist_to_plane(i_center) .lt. 0.e0_dp ) then
            ! Center i lies in the solvent, should not be reduced
            do j_center = i_center+1, n_centers
               neighbors(i_center, j_center) = .false.
               neighbors(j_center, i_center) = .false.
            enddo
         else
            ! Actually find neighbors
            do j_center = i_center+1, n_centers
               if ( dist_to_plane(j_center) .lt. 0.e0_dp ) then
                  ! Center j lies in the solvent, should not be reduced
                  neighbors(i_center, j_center) = .false.
                  neighbors(j_center, i_center) = .false.
               elseif ( sqrt(sum( (centers(i_center)%coord-centers(j_center)%coord)**2 )) &
                            .gt. min(dist_to_plane(i_center), dist_to_plane(j_center)) ) then
                  ! At least one of the centers lies closer to the interface than to the other center
                  neighbors(i_center, j_center) = .false.
                  neighbors(j_center, i_center) = .false.
               else
                  ! Both centers lie in the other solvent, and closer to each other
                  ! than to the interface
                  neighbors(i_center, j_center) = .true.
                  neighbors(j_center, i_center) = .true.
               endif ! j in solvent; i and j neighbors
            enddo ! j_center
         endif ! i in solvent
      enddo ! i_center

   end subroutine find_neighboring_centers_ifp

   subroutine cluster_neighboring_centers(neighbors, clusters_out)

      implicit none

      logical, intent(in) :: neighbors(:,:)

      logical, allocatable, intent(out) :: clusters_out(:,:)

      integer :: n_centers, n_clusters, i_center, i_cluster, j_center, k_center

      logical, allocatable :: clusters(:,:), cluster(:), cluster_new(:)

      logical :: have_common_cluster, has_subset

      n_centers = size(neighbors(1,:))

      allocate(clusters(n_centers,n_centers))
      clusters = .false.

      allocate(cluster(n_centers))
      cluster = .false.
      allocate(cluster_new(n_centers))
      cluster_new = .false.

      n_clusters = 0

      do i_center = 1, n_centers
         ! Check if center i has neighbors with which it is not in a cluster yet

         ! Iterate backwards here, or else it will check i=j ('Is there a cluster
         ! with i yet?') FIRST, but this should be checked LAST (only if i has no
         ! neighbors)
         check_neighbors: do j_center = n_centers, i_center, -1
            if (.not. neighbors(i_center,j_center)) &
               cycle check_neighbors
            have_common_cluster = .false.
            check_for_common_clusters: do i_cluster = 1, n_clusters
               if (clusters(i_cluster,i_center) .and. clusters(i_cluster,j_center)) then
                  have_common_cluster = .true.
                  exit check_for_common_clusters
               endif
            enddo check_for_common_clusters

            if (.not. have_common_cluster) then
               ! Create new cluster that has both i and j in it
               ! Start with the common neighbors of i and j
               cluster = neighbors(i_center,:) .and. neighbors(j_center,:)
               converge_cluster: do
                  cluster_new = cluster
                  remove_not_common: do k_center = 1, n_centers
                     if (.not. cluster(k_center)) &
                        cycle remove_not_common
                     ! Throw out all that are not mutually common with at least one other
                     cluster_new = cluster_new .and. neighbors(k_center,:)
                  enddo remove_not_common

                  ! If none got thrown out, add cluster to list
                  if (all(cluster .eqv. cluster_new)) then
                     exit converge_cluster
                  ! Else, throw out the first one that is in cluster,
                  ! but not in cluster_new and reiterate
                  else
                     remove_first_not_common: do k_center = 1, n_centers
                        if (cluster(k_center) .and. .not. cluster_new(k_center)) then
                           cluster(k_center) = .false.
                           exit remove_first_not_common
                        endif
                     enddo remove_first_not_common
                  endif
               enddo converge_cluster

               ! Add cluster to list
               ! Check if existing cluster is subset of new one
               has_subset = .false.
               check_subsets: do i_cluster = 1, n_clusters
                  if (all(clusters(i_cluster,:) .eqv. &
                          (cluster .and. clusters(i_cluster,:)))) then
                     ! Overwrite
                     clusters(i_cluster,:) = cluster
                     has_subset = .true.
                     exit check_subsets
                  endif ! clusters(i_cluster,:) is subset of cluster
               enddo check_subsets

               ! Else, append
               if (.not. has_subset) then
                  n_clusters = n_clusters + 1
                  clusters(n_clusters,:) = cluster
               endif ! .not. has_subset
            endif ! .not. have_common_cluster
         enddo check_neighbors
      enddo ! i_center

      allocate(clusters_out(n_clusters,n_centers))

      clusters_out = clusters(1:n_clusters,:)

      deallocate(clusters)
      deallocate(cluster)
      deallocate(cluster_new)

   end subroutine cluster_neighboring_centers

   subroutine reduce_clusters(clusters)

      implicit none

      logical, allocatable, intent(inout) :: clusters(:,:)

      logical, allocatable :: clusters_final(:,:)
      integer, allocatable :: additional_centers(:), priority(:)

      integer :: n_centers, n_clusters, i_center, n_clusters_final, i_cluster, j_cluster, prio

      n_centers = size(clusters(1,:))
      n_clusters = size(clusters(:,1))
      n_clusters_final = 0

      ! Remove (i.e. flush to .false.) clusters which are subsets of another
      ! cluster, so that they will not skew center priorities
      do i_cluster = 1, n_clusters
         do j_cluster = 1, n_clusters
            if (i_cluster .ne. j_cluster) then
               if (all(clusters(i_cluster,:) .eqv. &
                          clusters(i_cluster,:) .and. clusters(j_cluster,:))) then
                  ! i is subset of j
                  clusters(i_cluster,:) = .false.
               endif
            endif
         enddo
      enddo

      allocate(priority(n_centers))
      priority = 0

      ! Sort by priority, such that centers which are only contained in one cluster
      ! have highest priority etc.
      do i_center = 1, n_centers
         priority(i_center) = count(clusters(:,i_center))
      enddo

      allocate(clusters_final(n_clusters,n_centers))
      clusters_final = .false.

      allocate(additional_centers(n_clusters))

      ! For each center, (implicitly) sorted by priority, check if it is already in a
      ! final cluster. If not, add one
      reduce_clusters_main: do prio = 1, n_clusters
         iter_center: do i_center = 1, n_centers
            if (priority(i_center) .ne. prio) cycle iter_center
            if (any(clusters_final(:,i_center))) cycle iter_center

            additional_centers = 0
            ! For each cluster that contains i_center, see how many additional centers
            ! (which are in no cluster yet) it would add

            do j_cluster = 1, n_clusters
               if (clusters(j_cluster,i_center)) then
                  additional_centers(j_cluster) = count(clusters(j_cluster,:) &
                          .and. .not. any(clusters_final, dim=1))
               else
                  additional_centers(j_cluster) = 0
               endif
            enddo

            j_cluster = maxloc(additional_centers, dim=1)
            n_clusters_final = n_clusters_final + 1

            ! Do not add centers which are already contained in other cluster
            clusters_final(n_clusters_final,:) = clusters(j_cluster,:) &
                    .and. .not. any(clusters_final, dim=1)

            ! If all centers are already in a cluster, skip to the end
            if (all(any(clusters_final, dim=1))) exit reduce_clusters_main

         enddo iter_center
      enddo reduce_clusters_main

      deallocate(priority)
      deallocate(additional_centers)

      deallocate(clusters)
      allocate(clusters(n_clusters_final,n_centers))
      clusters = clusters_final(1:n_clusters_final,:)

      deallocate(clusters_final)

   end subroutine reduce_clusters

   subroutine reduce_centers(centers, clusters)
     implicit none

     type(SolHarmBasisCenter), allocatable, intent(inout) :: centers(:)
     logical, intent(in) :: clusters(:,:)

     type(SolHarmBasisCenter), allocatable :: centers_temp(:)

     integer :: n_centers, n_clusters, i_cluster

     real(dp), allocatable :: weights(:), weighted_coords(:,:)
     real(dp) :: coord_temp(3), log2

     character(132) :: info_str, format_str

     log2 = log(2.e0_dp)

     n_centers = size(centers)
     n_clusters = size(clusters(:,1))

     allocate(centers_temp(n_clusters))
     allocate(weights(n_centers))
     allocate(weighted_coords(n_centers,3))

     weights = 1.e0_dp / real(count(clusters, dim=1),dp)

    ! write(format_str, '(A,I4,A)') '(', n_centers, 'E15.8)'
    ! write(info_str, format_str) weights
    ! call localorb_info(info_str)

     do i_center = 1, n_centers
        weighted_coords(i_center,:) = weights(i_center) * centers(i_center)%coord(:)
     enddo

     do i_cluster = 1, n_clusters
        coord_temp = 0.e0_dp
        centers_temp(i_cluster)%lmin = 100
        centers_temp(i_cluster)%lmax = 0
        do i_center = 1, n_centers
           if (clusters(i_cluster,i_center)) then
              coord_temp = coord_temp + weighted_coords(i_center,:)
              centers_temp(i_cluster)%lmin = min(centers_temp(i_cluster)%lmin, centers(i_center)%lmin)
              if (mpe_lmax_ep_increase_cluster) then
                 centers_temp(i_cluster)%lmax = centers_temp(i_cluster)%lmax + &
                      2**centers(i_center)%lmax * weights(i_center) ! number of `poles', weighted
              else
                 centers_temp(i_cluster)%lmax = max(centers_temp(i_cluster)%lmax, centers(i_center)%lmax)
              endif
           endif
        enddo

        centers_temp(i_cluster)%coord = coord_temp / sum(weights, mask=clusters(i_cluster,:))
        if (mpe_lmax_ep_increase_cluster) &
           centers_temp(i_cluster)%lmax  = ceiling(log(real(centers_temp(i_cluster)%lmax,dp))/log2)
     enddo

     deallocate(centers)

     allocate(centers(n_clusters), source=centers_temp)

     write(info_str, '(4X,A)') 'Placing expansion centers for external potential at:'
     call localorb_info(info_str)
     do i_center = 1, n_clusters
        write(info_str, '(4X,I3,3E16.8)') centers(i_center)%lmax, centers(i_center)%coord * bohr
        call localorb_info(info_str)
     enddo

     deallocate(centers_temp)
     deallocate(weights)
     deallocate(weighted_coords)



  end subroutine reduce_centers

end subroutine define_problem
!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/mpe_get_continuum_name
!  NAME
!    mpe_get_continuum_name
!  SYNOPSIS

function mpe_get_continuum_name(i_dc, stop_unassigned) result(name)

!  PURPOSE
!    This function generates human-readable names for index values
!    regarding "dielectric continua".
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: i_dc
   logical, intent(in), optional :: stop_unassigned
   character(len=32) :: name

!  INPUTS
!   o i_dc -- continuum index
!  RETURNS
!   o name -- Human readable name
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   logical :: opt_stop_unassigned
   opt_stop_unassigned = .true.
   if (present(stop_unassigned)) opt_stop_unassigned = stop_unassigned

   if (i_dc .lt. 0) then
       call f_stop("***ValueError: dielectric index must be positive!")
   elseif (any(i_dc .eq. IND_DC % R)) then
      write(name,"(A)") "Inside Solute Cavity"
   elseif (any(i_dc .eq. IND_DC % Q)) then
      write(name,"(A)") "Dielectric 'Q'"
   elseif (any(i_dc .eq. IND_DC % O)) then
      write(name,"(A)") "Dielectric 'O'"
   elseif (any(i_dc .eq. IND_DC % NA)) then
      write(name,"(A)") "Not Assigned"
      if (opt_stop_unassigned) &
         call f_stop("***ValueError: encountered unassigned dielectric index!")
   else
      call f_stop("***NotImplementedError: unknown dielectric index")
   endif ! i_dc

end function mpe_get_continuum_name
!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/write_xml_SpVecTuple
!  NAME
!    write_xml_SpVecTuple
!  SYNOPSIS

subroutine write_xml_SpVecTuple( xml_file, inst, tag, ind )

!  PURPOSE
!    This subroutine writes an instance of type SpVecTuple to an XML file
!    using the xml_write module
!  USES
   implicit none
!  ARGUMENTS
   type(Xml_file_t), intent(inout) :: xml_file
   type(SpVecTuple), intent(in) :: inst
   character(*), intent(in), optional :: tag
   integer, intent(in), optional :: ind
!  INPUTS
!   o xml_file -- XML file handle
!   o inst -- instance
!   o tag -- tag for entry in XML file
!   o ind -- index in parent array, mutually exclusive with "tag"
!  OUTPUT
!   o xml_file -- XML file handle
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(*), parameter :: name_in_array = "element"

   if (present(ind)) then
      call xml_open(file=xml_file, tag=name_in_array, &
               attrs=(/ xml_attr("index", ind) /))
   elseif (present(tag)) then
      call xml_open(file=xml_file, tag=trim(tag), &
               attrs=(/ xml_attr("type", "SpVecTuple") /))
   else
      call xml_open(file=xml_file, tag=name_in_array)
   endif

   call xml_elem(file=xml_file, tag="ind", val=inst%ind)
   call xml_elem(file=xml_file, tag="val", val=inst%val)

   call xml_close(file=xml_file) !</inst>

end subroutine write_xml_SpVecTuple
!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/write_xml_SolHarmBasisCenter
!  NAME
!    write_xml_SolHarmBasisCenter
!  SYNOPSIS

subroutine write_xml_SolHarmBasisCenter( xml_file, inst, tag, ind )

!  PURPOSE
!    This subroutine writes an instance of type SolHarmBasisCenter to an XML file
!    using the xml_write module
!  USES
   implicit none
!  ARGUMENTS
   type(Xml_file_t), intent(inout) :: xml_file
   type(SolHarmBasisCenter), intent(in) :: inst
   character(*), intent(in), optional :: tag
   integer, intent(in), optional :: ind
!  INPUTS
!   o xml_file -- XML file handle
!   o inst -- instance
!   o tag -- tag for entry in XML file
!   o ind -- index in parent array, mutually exclusive with "tag"
!  OUTPUT
!   o xml_file -- XML file handle
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(*), parameter :: name_in_array = "element"
   integer :: i_sp

   if (present(ind)) then
      call xml_open(file=xml_file, tag=name_in_array, &
               attrs=(/ xml_attr("index", ind) /))
   elseif (present(tag)) then
      call xml_open(file=xml_file, tag=trim(tag), &
               attrs=(/ xml_attr("type", "SolHarmBasisCenter") /))
   else
      call xml_open(file=xml_file, tag=name_in_array)
   endif

   call xml_elem(file=xml_file, tag="coord", arr=inst%coord)
   call xml_elem(file=xml_file, tag="lmin", val=inst%lmin)
   call xml_elem(file=xml_file, tag="lmax", val=inst%lmax)
   call xml_elem(file=xml_file, tag="rscale", val=inst%rscale)

   call xml_open(file=xml_file, tag="coeff", &
               attrs=(/ xml_attr("type", "SpVecTuple"), &
                        xml_attr("shape", shape(inst%coeff) ) /))
   do i_sp = lbound(inst%coeff,1), ubound(inst%coeff,1)
      call write_xml_SpVecTuple( xml_file=xml_file, ind=i_sp, &
                                 inst=inst%coeff(i_sp) )
   enddo ! i_sp
   call xml_close(file=xml_file) !</coeff>

   call xml_close(file=xml_file) !</inst>

end subroutine write_xml_SolHarmBasisCenter
!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/write_xml_SolHarmBasis
!  NAME
!    write_xml_SolHarmBasis
!  SYNOPSIS

subroutine write_xml_SolHarmBasis( xml_file, inst, tag, ind )

!  PURPOSE
!    This subroutine writes an instance of type SolHarmBasis to an XML file
!    using the xml_write module
!  USES
   implicit none
!  ARGUMENTS
   type(Xml_file_t), intent(inout) :: xml_file
   type(SolHarmBasis), intent(in) :: inst
   character(*), intent(in), optional :: tag
   integer, intent(in), optional :: ind
!  INPUTS
!   o xml_file -- XML file handle
!   o inst -- instance
!   o tag -- tag for entry in XML file
!   o ind -- index in parent array, mutually exclusive with "tag"
!  OUTPUT
!   o xml_file -- XML file handle
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(*), parameter :: name_in_array = "element"
   integer :: i_c

   if (present(ind)) then
      call xml_open(file=xml_file, tag=name_in_array, &
               attrs=(/ xml_attr("index", ind) /))
   elseif (present(tag)) then
      call xml_open(file=xml_file, tag=trim(tag), &
               attrs=(/ xml_attr("type", "SolHarmBasis") /))
   else
      call xml_open(file=xml_file, tag=name_in_array)
   endif

   call xml_elem(file=xml_file, tag="name", val=trim(inst%name))

   call xml_open(file=xml_file, tag="centers", &
               attrs=(/ xml_attr("type", "SolHarmBasis"), &
                        xml_attr("shape", shape(inst%centers)) /))
   do i_c = lbound(inst%centers,1), ubound(inst%centers,1)
      call write_xml_SolHarmBasisCenter( xml_file=xml_file, ind=i_c, &
                                          inst=inst%centers(i_c) )
   enddo ! i_c
   call xml_close(file=xml_file) !</centers>

   call xml_elem(file=xml_file, tag="solharm_type", &
                     val=inst%solharm_type)

   call xml_close(file=xml_file) !</inst>

end subroutine write_xml_SolHarmBasis
!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/write_xml_DielectricContinuum
!  NAME
!    write_xml_DielectricContinuum
!  SYNOPSIS

subroutine write_xml_DielectricContinuum( xml_file, inst, tag, ind )

!  PURPOSE
!    This subroutine writes an instance of type DielectricContinuum 
!    to an XML file using the xml_write module
!  USES
   implicit none
!  ARGUMENTS
   type(Xml_file_t), intent(inout) :: xml_file
   type(DielectricContinuum), intent(in) :: inst
   character(*), intent(in), optional :: tag
   integer, intent(in), optional :: ind
!  INPUTS
!   o xml_file -- XML file handle
!   o inst -- instance
!   o tag -- tag for entry in XML file
!   o ind -- index in parent array, mutually exclusive with "tag"
!  OUTPUT
!   o xml_file -- XML file handle
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(*), parameter :: name_in_array = "element"

   if (present(ind)) then
      call xml_open(file=xml_file, tag=name_in_array, &
               attrs=(/ xml_attr("index", ind) /))
   elseif (present(tag)) then
      call xml_open(file=xml_file, tag=trim(tag), &
               attrs=(/ xml_attr("type", "DielectricContinuum") /))
   else
      call xml_open(file=xml_file, tag=name_in_array)
   endif

   call xml_elem(file=xml_file, tag="eps", val=inst%eps)

   select type (b => inst%basis)
      type is (SolHarmBasis)
         call write_xml_SolHarmBasis( xml_file=xml_file, tag="basis", &
                                       inst=b )
      class default
         !TODO throw error
   end select

   call xml_close(file=xml_file) !</inst>

end subroutine write_xml_DielectricContinuum
!******
!-------------------------------------------------------------------------------
!****s* mpe_interface/write_xml_continua
!  NAME
!    write_xml_continua
!  SYNOPSIS

subroutine write_xml_continua( xml_file, continua )

!  PURPOSE
!    This subroutine writes an instance of type DielectricContinuum 
!    to an XML file using the xml_write module
!  USES
   implicit none
!  ARGUMENTS
   type(Xml_file_t), intent(inout) :: xml_file
   type(DielectricContinuum), intent(in) :: continua(:)
!  INPUTS
!   o xml_file -- XML file handle
!   o continua -- dielectric continua in problem
!  OUTPUT
!   o xml_file -- XML file handle
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   integer :: i_dc

   call xml_open(file=xml_file, tag="continua", &
                  attrs=(/ xml_attr("type", "DielectricContinuum"), &
                           xml_attr("shape", shape(continua)) /))
   do i_dc = lbound(continua,1), ubound(continua,1)
      call write_xml_DielectricContinuum( xml_file=xml_file, ind=i_dc, &
                                          inst=continua(i_dc) )
   enddo ! i_dc
   call xml_close(file=xml_file) !</continua>

end subroutine write_xml_continua
!******
!****s* mpe_interface/mpe_calculate_reaction_field
!  NAME
!    mpe_calculate_reaction_field
!  SYNOPSIS

subroutine mpe_calculate_reaction_field( &
      coord, &
      rho_free, rho_mp, rho_total, &
      hartree_potential, &
      dc_index, &
      potential_correction )

!  PURPOSE
!    Calculate the reaction field at points of the integration grid on-the-fly.
!    If dc_index is undefined on entry, its value is determined and returned.
!
!  USES
   implicit none
!  ARGUMENTS
   real(dp), intent(in) :: coord(3), rho_free, rho_mp, rho_total
   real(dp), intent(in) :: hartree_potential
   integer, intent(inout) :: dc_index
   real(dp), intent(out) :: potential_correction
!  INPUTS
!   o coords -- the point where the reaction field is evaluated
!   o rho_free -- free electron density at this point
!   o rho_mp -- multipole-expanded electron density at this point
!   o rho_total -- total electron density at this point
!   o hartree_potential -- Hartree potential at this point
!   o dc_index -- index of dielectric continuum the point belongs to.
!                 If undefined on entry, the index will be determined and overwritten.
!  OUTPUT
!   o dc_index -- index of dielectric continuum the point belongs to.
!   o potential_correction -- the correction that needs to be added to the
!        Hartree potential in order to obtain the full electrostatic potential
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   character(*), parameter :: func = 'mpe_calculate_reaction_field'
   real(dp) :: vacuum_regularized_reaction_field(1)

   if (any(dc_index .eq. IND_DC % NA)) then
      select case (isc_cavity_type)
      case (ISC_CONST % CAVITY_RhoFree)
         dc_index = get_continuum_index_rho(coord, &
            rho=rho_free, rho_iso=isc_isodensity_value)
      case (ISC_CONST % CAVITY_RhoMPStat, &
            ISC_CONST % CAVITY_RhoMPDyn)
         dc_index = get_continuum_index_rho(coord, &
            rho=rho_total, rho_iso=isc_isodensity_value)
      case (ISC_CONST % CAVITY_OvlpSph)
         call f_stop("Unexpected cavity type", func)
      case default
         call f_stop("Not implemented error", func)
      end select
   endif
   call calculate_potential_at_points(1, coord, &
         isc_continua(dc_index) % basis, &
         vacuum_regularized_reaction_field(1) )
   ! just returning this potential would imply that the regularization was
   ! performed with the regular Hartree potential. However, we regularize
   ! with the Hartree potential divided by the dielectric permittivity, thus
   potential_correction = vacuum_regularized_reaction_field(1) + &
      (1.e0_dp/isc_continua(dc_index) % eps - 1.e0_dp) * hartree_potential

end subroutine mpe_calculate_reaction_field
!******
!****f* mpe_interface/get_continuum_index_radii
!  NAME
!    get_continuum_index_radii
!  SYNOPSIS
function get_continuum_index_radii(point, radii_sq) result(idx)
   implicit none
   real(dp), intent(in) :: point(3), radii_sq(n_atoms)
   integer :: idx

   real(dp), allocatable :: dist_tab_sq(:), dir_tab(:,:)

   allocate(dist_tab_sq(n_atoms))
   allocate(dir_tab(3,n_atoms))
   call tab_atom_centered_coords_v2( &
            point, dist_tab_sq, dir_tab)
   if (any(dist_tab_sq < radii_sq)) then
      idx = IND_DC % R(1)
   else
           if (.not. ifp) then
         idx = IND_DC % Q(1)
      else
         if (0 > calculate_distance_to_plane(point, &
                 ifp_normal, &
                 ifp_dist) ) then
            idx = isc_interfaces(IND_IF % QO(1)) % dc_ind_neg
         else
            idx = isc_interfaces(IND_IF % QO(1)) % dc_ind_pos
         endif
      endif
   endif
end function get_continuum_index_radii
!******
!****f* mpe_interface/get_continuum_index_rho
!  NAME
!    get_continuum_index_rho
!  SYNOPSIS
function get_continuum_index_rho(point, rho, rho_iso) result(idx)
   implicit none
   real(dp), intent(in) :: point(3), rho, rho_iso
   integer :: idx

   if (rho .ge. rho_iso) then
      idx = IND_DC % R(1)
   else
      if (.not. ifp) then
         idx = IND_DC % Q(1)
      else
         if (0 > calculate_distance_to_plane(point, &
                 ifp_normal, &
                 ifp_dist) ) then
            idx = isc_interfaces(IND_IF % QO(1)) % dc_ind_neg
         else
            idx = isc_interfaces(IND_IF % QO(1)) % dc_ind_pos
         endif
      endif
   endif
end function get_continuum_index_rho
!******
!****s* mpe_interface/mpe_calculate_interaction_energy_with_nuclei
!  NAME
!    mpe_calculate_interaction_energy_with_nuclei
!  SYNOPSIS

function mpe_calculate_interaction_energy_with_nuclei() &
      result(interaction_energy)

!  PURPOSE
!    Calculate the interaction energy of the nuclear charges with the
!    reaction field and store this information in the variable
!    "interaction_energy".
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp) :: interaction_energy
!  INPUTS
!    none
!  RETURNS
!   o interaction_energy -- interaction energy of nuclei with reaction field
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = &
                           'calculate_interaction_energy_reaction_field_nuclei'
   integer :: i_atom
   real(dp) :: potential_at_atom(1)

   interaction_energy = ZERO
   if (.not.module_initialized) return

   do i_atom = 1, n_atoms, 1

      ! skip ghost atoms as there is no nucleus
      if (empty(i_atom)) cycle

      ! we need the potential at the position of the atom
      call calculate_potential_at_points( &
            1, coords(1,i_atom), &
            isc_continua(mpe_dc_indices_atoms(i_atom)) % basis, &
            potential_at_atom )

      ! multiply solvent potential with nuclear charge and
      ! sum up to total interaction energy
      interaction_energy = interaction_energy + &
            species_z(species(i_atom)) * potential_at_atom(1)

   enddo ! i_atom

end function mpe_calculate_interaction_energy_with_nuclei
!******
!****s* mpe_interface/mpe_calculate_nonelectrostatic_energy
!  NAME
!    mpe_calculate_nonelectrostatic_energy
!  SYNOPSIS

function mpe_calculate_nonelectrostatic_energy() &
      result(energy)

!  PURPOSE
!    Evaluate non-electrostatic correction to total energy expression
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp) :: energy, surf, vol, energy_surf, energy_vol
!  INPUTS
!    none
!  RETURNS
!   o energy -- nonelectrostatic energy contribution
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = &
                           'mpe_calculate_nonelectrostatic_energy'
   
   integer :: i_if, i_p
   character(132) :: info_str

   intrinsic dot_product

   energy = ZERO
   if (.not.module_initialized) return

   select case (mpe_nonel_model)
   case(MPE_CONST % NONEL_UNDEF)
   case(MPE_CONST % NONEL_linOV)
      do i_if = 1, size(isc_interfaces), 1
         surf = 0.0e0_dp
         vol  = 0.0e0_dp
         do i_p = 1, size(isc_interfaces(i_if)%p), 1
            surf = surf+isc_interfaces(i_if)%p(i_p)%area
            vol  = vol +isc_interfaces(i_if)%p(i_p)%volume
         enddo
         energy_surf = surf*isc_interfaces(i_if)%surface_energy
         energy_vol  = vol *isc_interfaces(i_if)%volume_energy
         energy      = energy + energy_surf + energy_vol
      enddo
   case(MPE_CONST % NONEL_INVAL)
      call f_stop("Invalid non-electrostatic model specified", func)
   case default
      call f_stop("The requested non-electrostatic model is unknown", func)
   end select

end function mpe_calculate_nonelectrostatic_energy
!******
!****s* mpe_interface/mpe_store_energy_contributions
!  NAME
!    mpe_store_energy_contributions
!  SYNOPSIS

subroutine mpe_store_energy_contributions(rho_V, nuc_V, nonel)

!  PURPOSE
!    Store evaluated interaction energies in module variables.
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in), optional :: rho_V
   real(dp), intent(in), optional :: nuc_V
   real(dp), intent(in), optional :: nonel
!  INPUTS
!   o rho_V -- interaction energy of reaction field with electron debsity
!   o nuc_V -- interaction energy of reaction field with nuclei
!   o nonel -- total non-electrostatic contributions
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = &
                           'mpe_store_energy_contributions'

   if (present(rho_V)) then
      mpe_energy_contributions % rho_V = rho_V
   endif
   if (present(nuc_V)) then
      mpe_energy_contributions % nuc_V = nuc_V
   endif
   if (present(nonel)) then
      mpe_energy_contributions % nonel = nonel
   endif

end subroutine mpe_store_energy_contributions
!******
!****s* mpe_interface/mpe_read_energy_contributions
!  NAME
!    mpe_read_energy_contributions
!  SYNOPSIS

subroutine mpe_read_energy_contributions(rho_V, nuc_V, nonel, el, tot)

!  PURPOSE
!    Read evaluated interaction energies from module variables.
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(out), optional :: rho_V
   real(dp), intent(out), optional :: nuc_V
   real(dp), intent(out), optional :: nonel
   real(dp), intent(out), optional :: el
   real(dp), intent(out), optional :: tot
!  INPUTS
!    none
!  OUTPUT
!   o rho_V -- interaction energy of reaction field with electron debsity
!   o nuc_V -- interaction energy of reaction field with nuclei
!   o nonel -- total non-electrostatic contributions
!   o el -- total electrostatic contributions
!   o tot -- sum of electrostatic and non-electrostatic contributions
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = &
                           'mpe_read_energy_contributions'

   real(dp) :: tmp_el

   tmp_el = 0.5e0_dp * ( &
      mpe_energy_contributions % rho_V - &
      mpe_energy_contributions % nuc_V )

   if (present(rho_V)) then
      rho_V = mpe_energy_contributions % rho_V
   endif
   if (present(nuc_V)) then
      nuc_V = mpe_energy_contributions % nuc_V
   endif
   if (present(nonel)) then
      nonel = mpe_energy_contributions % nonel
   endif
   if (present(el)) then
      el = tmp_el
   endif
   if (present(tot)) then
      tot = tmp_el + mpe_energy_contributions % nonel
   endif

end subroutine mpe_read_energy_contributions


!!------------------------------------------------------------------------------
!!
!! Some test routines


subroutine check_consistency_spherical_harmonics(l_max)

   use cartesian_ylm, only: n_max_cartesian, &
                            initialize_cartesian_ylm, &
                            evaluate_onecenter_cartesians, &
                            tab_ylm_onecenter_cartesian
   use localorb_io, only: use_unit
   implicit none

   integer, intent(in) :: l_max

   integer, parameter :: npoints = 5
   real(dp), dimension(3,npoints), parameter :: points = &
        reshape( (/ &
            -5.0e0_dp,  3.0e0_dp, -7.0e0_dp, &
             3.0e0_dp,  0.0e0_dp,  0.5e0_dp, &
             4.0e0_dp,  0.0e0_dp,  0.0e0_dp, &
             3.0e0_dp, -3.0e0_dp,  2.0e0_dp, &
             0.0e0_dp, -1.0e0_dp, -1.0e0_dp  &
                /), shape(points))

   character(*), parameter :: func = 'check_consistency_spherical_harmonics'
   integer :: info

   ! cartesian_ylm
   real(dp), dimension(:,:), allocatable :: cartesians
   real(dp), dimension(:), allocatable :: rlylm

   ! ylm_tab
   real(dp), dimension(:), allocatable :: dist_tab
   real(dp), dimension(:,:), allocatable :: trigonom_tab
   real(dp), dimension(:,:), allocatable :: ylm_tab


   integer :: i_point, i_l, i_m, i_lm

   ! initialize cartesians
   call initialize_cartesian_ylm(l_max)


   ! ALLOCATION

   ! cartesian_ylm
   call f_allocate(cartesians, 1,n_max_cartesian, 0,l_max, name="cartesians")
   call f_allocate(rlylm, (l_max+1)**2, name="rlylm")

   ! ylm_tab
   call f_allocate(dist_tab, npoints, name="dist_tab")
   call f_allocate(trigonom_tab, 4, npoints, name="trigonom_tab")
   call f_allocate(ylm_tab, (l_max+1)**2, npoints, name="ylm_tab")


   write(use_unit,*) 'calculate spherical harmonics and plot difference'

   ! tab_ylm
   call tab_trigonom_p0(npoints, points, trigonom_tab)
   dist_tab = sqrt(sum(points*points,dim=1))

   do i_point = 1,1,1! npoints, 1

      ! cartesian_ylm
      call evaluate_onecenter_cartesians(points(:,i_point), l_max, cartesians)
      call tab_ylm_onecenter_cartesian(l_max, cartesians, rlylm)

      call increment_ylm &
          ( trigonom_tab(1,i_point), trigonom_tab(2,i_point), &
          trigonom_tab(3,i_point), trigonom_tab(4,i_point), &
          0, l_max, ylm_tab(1,i_point) )

      i_lm = 0
      do i_l = 0, l_max, 1
         do i_m = -i_l, i_l, 1
            i_lm = i_lm + 1
            write(use_unit,*) 'point', i_point, ' l=',i_l, ' m=', i_m
            write(use_unit,*) 'diff:', rlylm(i_lm) - ylm_tab(i_lm,i_point)*dist_tab(i_point)**i_l
         enddo
      enddo

   enddo ! i_point
   write(use_unit,*) 'spherical harmonics done'


   ! DEALLOCATION
   call f_deallocate(rlylm, name="rlylm")
   call f_deallocate(cartesians, name="cartesians")

   call f_deallocate(dist_tab, name="dist_tab")
   call f_deallocate(trigonom_tab, name="trigonom_tab")
   call f_deallocate(ylm_tab, name="ylm_tab")

end subroutine check_consistency_spherical_harmonics



end module mpe_interface

