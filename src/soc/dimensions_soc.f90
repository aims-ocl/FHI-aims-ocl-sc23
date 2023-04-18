!****h* FHI-aims/dimensions_soc
!  NAME
!    dimensions_soc -- Various data-structure-related variables for
!                      second-variational spin-orbit coupling
!  SYNOPSIS
module dimensions_soc
!  PURPOSE
!    This module holds variables needed for manipulating the data structures in
!    second-variational spin-orbit coupling.  As the number implies, it is
!    analogous to the main dimensions module.
!
!    No subroutines should be attached to this module.
!
!    For the most part, variables in this module are initialized either during
!    the read-in of control.in and geometry.in or at the beginning of
!    the main SOC subroutine (currently named calculate_second_variational_soc)
!  USES
  implicit none
!  AUTHOR
!    William Huhn
!  HISTORY
!    August 2017 - Added.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  SOURCE

!     n_states_sr  :   Number of SR states to include in second-variational SOC.
!                      This is used to index the scalar-relativistic states, so
!                      it does not include spin channels.
!                      The default (and maximum) value is n_states.
!     sr_start_start : The starting index for SR states to include in
!                      second-variational SOC.
!                      This value should always be n_core_states_omit_from_soc/2 + 1,
!                      but I made it into its own variable to make the code cleaner.
!     n_states_soc :   Size of the SOC Hamiltonian matrix.
!                      The name exists for historical reasons, and this variable is not
!                      the analogue of n_states, as it is not used to index the
!                      eigenvalues or eigenvectors.  See "n_saved_states_soc" for more
!                      information
!                      The default (and maximum) value is 2*n_states.
!     n_core_states_omit_from_soc :
!                      Number of low-lying SR states to omit from SOC
!                      For compactness, we refer to them as core states, but they
!                      could also be semi-core or valence states
!                      Code will exit if this value would extend into unoccupied
!                      states (can't determine Fermi energy in this case)
!                      As of this writing, this value should be zero, because the
!                      deep core states contribute noticeably to the results.
!     n_high_states_omit_from_soc :
!                      Number of high-lying SR states to omit from SOC
!                      Code will exit if this value would extend into occupied
!                      states (can't determine Fermi energy in this case)
!     n_basis_soc :    Total number of basis functions for spin-orbit-coupled
!                      eigenvectors
!                      This value should be n_basis_soc_coll + n_basis_soc_ncoll
!                      The default case should be 2*n_basis:  n_basis many
!                      spacial basis functions times 2 spinors (spin up and down)
!     n_basis_soc_coll :
!                      The number of collinear basis elements for SOC,
!                      that is, the number that are either spin-up or spin-dn
!                      This number must be even, because it is assumed that
!                      the first n_basis_soc_coll/2 basis elements are spin-up
!                      and the second n_basis_soc_coll/2 basis elements are
!                      spin down
!                      Note that, because the SOC code currently relies heavily
!                      on the underlying basis set architecture (notably the
!                      indexing arrays), the code will break if this value is
!                      anything other than 2*n_basis, and there are many parts
!                      of the code that still use 2*n_basis directly.
!                      The default value is 2*n_basis
!     n_basis_soc_ncoll :
!                      The number of basis elements for SOC which are not
!                      assumed to be spin up or spin down
!                      (They may still be, they're just not required to be)
!                      The default value is 0
!                      As of this writing (3 October 2017) this functionality
!                      has not been implemented, so this value should always
!                      be zero
!     n_saved_states_soc :
!                      The number of SOC-perturbed eigenvalues and eigenvectors
!                      to save.  This is the proper indexing for the various
!                      SOC-perturbed quantities, allowing us to save memory
!                      by not saving low-lying states when they are not needed.
!                      The default value is n_states_soc, i.e. all possible states
!                      consistent with the SOC-perturbed Hamiltonian are stored
!     soc_saved_state_start :
!                      Index for the first saved SOC state in the space of all
!                      possible states, the former being relevant for indexing
!                      within the code and the latter being relevant for output
!                      of physical quantities.
!                      The default value is 1 (in retrospect, I wish I had made
!                      this an offset so that the value would be 0, oh well)
  integer :: n_states_sr
  integer :: sr_state_start
  integer :: n_states_soc
  integer :: n_core_states_omit_from_soc = 0
  integer :: n_high_states_omit_from_soc = 0
  integer :: n_basis_soc
  integer :: n_basis_soc_coll
  integer :: n_basis_soc_ncoll
  integer :: n_saved_states_soc
  integer :: soc_saved_state_start = 1
end module dimensions_soc
