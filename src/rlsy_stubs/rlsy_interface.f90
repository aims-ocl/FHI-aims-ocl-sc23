module rlsy_interface
    use mpi_tasks, only: aims_stop_coll
    !! This module defines the handle used to deal with symmetry
    private
    public :: rlsy_handle
    public :: rlsy_h
    public :: generate_rlsy_handle
    public :: rl_integration_grid
    public :: rl_spacegroup_operation
    public :: rl_spacegroup
    public :: rl_integration_batch_irr


    type rl_integration_batch_irr
        integer :: n_point
        real*8, dimension(:,:), allocatable :: folded_coordinate
    end type rl_integration_batch_irr

    type rl_integration_grid
        integer :: n_irr_batch
        type(rl_integration_batch_irr), dimension(:), allocatable :: irr_batch

    end type rl_integration_grid

    type rl_spacegroup_operation
        real*8, dimension(3,3) :: m
        integer, dimension(:), allocatable :: fmap

    end type rl_spacegroup_operation

    type rl_spacegroup
        type(rl_spacegroup_operation), dimension(:), allocatable :: op
        integer :: n_operation

    end type rl_spacegroup

    type rlsy_handle
        ! Some compilers could get angry at an empty type
        integer :: i
        type(rl_integration_grid) :: grid
        type(rl_spacegroup) :: spacegroup
    end type rlsy_handle


    ! Expose empty type
    type(rlsy_handle), save :: rlsy_h
contains

!> Same interface, does nothing.
!subroutine generate_rlsy_handle(rh,&
!    species,species_z,frac_coords,lattice_vector,&
!    basis_m,basis_l,basis_fn,basis_atom,outer_radius,&
!    r_grid_min,r_grid_inc,basis_wave_spl,basis_deriv_spl,outer_partition_radius,&
!    n_k_points_xyz,k_point_list,k_points_offset,&
!    atom_radius,free_rho_spl,free_drho_dr_spl,free_pot_es_spl,free_pot_es_at_zero,&
!    n_radial,n_angular,r_radial,r_angular,w_radial,w_angular,&
!    n_electrons,n_states,n_spin,&
!    elsi_solver,elsi_out_level,occupation_type,solver_method,&
!    elsi_elpa_n_single,occupation_width,occupation_acc,basis_threshold,&
!    verbosity,use_unit,mpi_world_communicator,nosym)
!    !> container that holds all symmetry information
!    class(rlsy_handle), intent(out) :: rh
!    !> structure information from AIMS
!    integer, dimension(:), intent(in) :: species
!    real*8, dimension(:), intent(in) :: species_z
!    real*8, dimension(:,:), intent(in) :: frac_coords
!    real*8, dimension(3,3), intent(in) :: lattice_vector
!    !> basis set information from AIMS
!    integer, dimension(:), intent(in) :: basis_m
!    integer, dimension(:), intent(in) :: basis_l
!    integer, dimension(:), intent(in) :: basis_fn
!    integer, dimension(:), intent(in) :: basis_atom
!    real*8, dimension(:), intent(in) :: outer_radius
!    real*8, dimension(:), intent(in) :: r_grid_min
!    real*8, dimension(:), intent(in) :: r_grid_inc
!    real*8, dimension(:,:,:), intent(in) :: basis_wave_spl
!    real*8, dimension(:,:,:), intent(in) :: basis_deriv_spl
!    real*8, dimension(:), intent(in) :: outer_partition_radius
!    !> k-mesh information from AIMS
!    integer, dimension(3), intent(in) :: n_k_points_xyz
!    real*8, dimension(:,:), intent(in) :: k_point_list
!    real*8, dimension(3), intent(in) :: k_points_offset
!    !> grid information from AIMS
!    integer, dimension(:), intent(in) :: n_radial
!    integer, dimension(:,:), intent(in) :: n_angular
!    real*8, dimension(:,:), intent(in) :: r_radial
!    real*8, dimension(:,:,:,:), intent(in) :: r_angular
!    real*8, dimension(:,:), intent(in) :: w_radial
!    real*8, dimension(:,:,:), intent(in) :: w_angular
!    ! integer, intent(in) :: use_batch_permutation
!    ! integer, intent(in) :: n_my_batches_work
!    ! type(batch_of_points), dimension(:), pointer, intent(in) :: batches_work
!    !> free atom
!    real*8, dimension(:), intent(in) :: atom_radius
!    real*8, dimension(:,:,:), intent(in) :: free_rho_spl
!    real*8, dimension(:,:,:), intent(in) :: free_drho_dr_spl
!    real*8, dimension(:,:,:), intent(in) :: free_pot_es_spl
!    real*8, dimension(:), intent(in) :: free_pot_es_at_zero
!    !> KS solution related things
!    real*8, intent(in) :: n_electrons
!    integer, intent(in) :: n_states
!    integer, intent(in) :: n_spin
!    !> ELSI related things
!    integer, intent(in) :: elsi_solver
!    integer, intent(in) :: elsi_out_level
!    integer, intent(in) :: occupation_type
!    integer, intent(in) :: solver_method
!    integer, intent(in) :: elsi_elpa_n_single
!    real*8, intent(in) :: occupation_width
!    real*8, intent(in) :: occupation_acc
!    real*8, intent(in) :: basis_threshold
!    !> talk a lot?
!    integer, intent(in) :: verbosity
!    !> which unit to talk to
!    integer, intent(in) :: use_unit
!    !> mpi world communicator
!    integer, intent(in) :: mpi_world_communicator
!    !> should I switch off all symmetries except for identity? Only for debugging.
!    logical, intent(in) :: nosym
!
!    ! Just stop the code and tell user to recompile.
!    call aims_stop_coll('You have to compile with RLSY enabled to use it.')
!
!end subroutine

!> Create all the information needed to play around with symmetry. This is likely far too much, and should get pruned down to the bare necessities once things converge.
subroutine generate_rlsy_handle(rh,&
    species,species_z,frac_coords,lattice_vector,&
    basis_m,basis_l,basis_fn,basis_atom,outer_radius,&
    r_grid_min,r_grid_inc,basis_wave_spl,basis_deriv_spl,outer_partition_radius,&
    n_k_points_xyz,k_point_list,k_points_offset,&
    atom_radius,free_rho_spl,free_drho_dr_spl,free_pot_es_spl,free_pot_es_at_zero, &
    n_radial,n_angular,r_radial,r_angular,w_radial,w_angular,scale_radial,l_hartree,&
    n_electrons,n_states,n_spin,&
    elsi_solver,elsi_out_level,occupation_type,solver_method,&
    elsi_elpa_n_single,occupation_width,occupation_acc,basis_threshold,&
    verbosity,use_unit,mpi_world_communicator,nosym)
    !> container that holds all symmetry information
    class(rlsy_handle), intent(out) :: rh
    !> structure information from AIMS
    integer, dimension(:), intent(in) :: species
    real*8, dimension(:), intent(in) :: species_z
    real*8, dimension(:,:), intent(in) :: frac_coords
    real*8, dimension(3,3), intent(in) :: lattice_vector
    !> basis set information from AIMS
    integer, dimension(:), intent(in) :: basis_m
    integer, dimension(:), intent(in) :: basis_l
    integer, dimension(:), intent(in) :: basis_fn
    integer, dimension(:), intent(in) :: basis_atom
    real*8, dimension(:), intent(in) :: outer_radius
    real*8, dimension(:), intent(in) :: r_grid_min
    real*8, dimension(:), intent(in) :: r_grid_inc
    real*8, dimension(:,:,:), intent(in) :: basis_wave_spl
    real*8, dimension(:,:,:), intent(in) :: basis_deriv_spl
    real*8, dimension(:), intent(in) :: outer_partition_radius
    !> k-mesh information from AIMS
    integer, dimension(3), intent(in) :: n_k_points_xyz
    real*8, dimension(:,:), intent(in) :: k_point_list
    real*8, dimension(3), intent(in) :: k_points_offset
    !> grid information from AIMS
    integer, dimension(:), intent(in) :: n_radial
    integer, dimension(:,:), intent(in) :: n_angular
    real*8, dimension(:,:), intent(in) :: r_radial
    real*8, dimension(:,:,:,:), intent(in) :: r_angular
    real*8, dimension(:,:), intent(in) :: w_radial
    real*8, dimension(:,:,:), intent(in) :: w_angular
    real*8, dimension(:), intent(in) :: scale_radial
    integer, dimension(:), intent(in) :: l_hartree
    ! integer, intent(in) :: use_batch_permutation
    ! integer, intent(in) :: n_my_batches_work
    ! type(batch_of_points), dimension(:), pointer, intent(in) :: batches_work
    !> free atom
    real*8, dimension(:), intent(in) :: atom_radius
    real*8, dimension(:,:,:), intent(in) :: free_rho_spl
    real*8, dimension(:,:,:), intent(in) :: free_drho_dr_spl
    real*8, dimension(:,:,:), intent(in) :: free_pot_es_spl
    real*8, dimension(:), intent(in) :: free_pot_es_at_zero
    !> KS solution related things
    real*8, intent(in) :: n_electrons
    integer, intent(in) :: n_states
    integer, intent(in) :: n_spin
    !> ELSI related things
    integer, intent(in) :: elsi_solver
    integer, intent(in) :: elsi_out_level
    integer, intent(in) :: occupation_type
    integer, intent(in) :: solver_method
    integer, intent(in) :: elsi_elpa_n_single
    real*8, intent(in) :: occupation_width
    real*8, intent(in) :: occupation_acc
    real*8, intent(in) :: basis_threshold
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> which unit to talk to
    integer, intent(in) :: use_unit
    !> mpi world communicator
    integer, intent(in) :: mpi_world_communicator
    !> should I switch off all symmetries except for identity? Only for debugging.
    logical, intent(in) :: nosym


    call aims_stop_coll('You have to compile with RLSY enabled to use it.')
end subroutine

end module rlsy_interface
