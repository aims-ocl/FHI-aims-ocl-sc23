!****s* FHI-aims/init_aims
!  NAME
!    init_aims
!  SYNOPSIS

subroutine init_aims

  !  PURPOSE
  !    This subroutine only initializes the most important data fields from
  !    control.in and geometry.in.
  !  USES

  use localorb_io
  use mpi_utilities
  use dimensions
  use runtime_choices
  use hartree_fock
  use basis
  use prodbas
  use physics, only: n_electrons

  !  ARGUMENTS
  !  INPUTS
  !    none
  !  OUTPUT
  !    none
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

  implicit none

  logical :: converged

  ! if an MPI-environment exists, initialize it - else, silently
  ! circumvent all MPI functionality using dummy subroutines
  call initialize_mpi()

  call write_version_stamp()

  ! Read all input files:
  ! * Set up initial dimensions
  ! * allocate initial data structures
  ! * read all input information from files control.in and geometry.in
  call read_input_data()

  ! Set up all fixed computational quantities:
  ! * Basic integration grid specifications for all species
  ! * Free-atom DFT solution for each species on a logarithmic grid
  ! * Radial basis functions for each species, tabulated on log. grids
  ! * Coefficients that prepare the computation of Ylm functions and
  !   their derivatives from cartesian products x^l*y^m*z^n
  call prepare_scf()

  ! Set up all administrative arrays connected to periodic boundary
  ! conditions.  These depend on the atomic geometry.
  call initialize_bc_dependent_lists()

  if (use_prodbas) call initialize_prodbas()

end subroutine init_aims
!******
