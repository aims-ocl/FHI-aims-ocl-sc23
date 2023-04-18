!****h* FHI-aims/scalapack_soc
!*  NAME
!*    scalapack_soc -- Various ScaLAPACK-related variables for
!*                     second-variational spin-orbit coupling
!*  SYNOPSIS
module scalapack_soc
!*  PURPOSE
!*    This module holds the variables for handling the ScaLAPACK distribution(s)
!*    specific to second-variational spin-orbit coupling.  Right now there are
!*    two sets of variables, owing to the two different data structures commonly
!*    found in the SOC code:
!*    - the SOC-perturbed Hamiltonian (an n_states_soc*n_states_soc matrix) that 
!*      is diagonalized to generate the SOC-perturbed eigenvalues and
!*      eigenvectors, and
!*    - the SOC-perturbed eigenvectors (a n_basis_soc*n_saved_states_soc
!*      matrix).
!*    Similar to the ScaLAPACK descriptor for the scalar-relativistic
!*    eigenvectors, we define the ScaLAPACK descriptor for the SOC-perturbed
!*    eigenvectors for an n_basis_soc*n_basis_soc matrix, allowing use to reuse
!*    this descriptor in other contexts.
!*
!*    The subroutines in this module are wrappers around fundamental
!*    subroutines.  This is deliberate.  The fundamental subroutines may be
!*    found in scalapack_generic_wrapper, which (as the name implies) is
!*    intended to be used to generalize the ScaLAPACK distribution to generic
!*    matrices.
!*  USES
  use scalapack_wrapper, only: dlen_
  implicit none
!*  AUTHOR
!*    William Huhn (Duke University)
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!*    the terms and conditions of the respective license agreement."
!*  HISTORY
!*    November 2016, added.
!*  SOURCE
  PUBLIC

  ! ScaLAPACK quantities for n_states_soc*n_states_soc matrices
  integer :: block_size_soc
  integer :: sc_desc_soc(dlen_)
  integer :: mb_soc, nb_soc
  integer :: mxld_soc, mxcol_soc
  integer :: n_my_rows_soc, n_my_cols_soc
  ! Dimensions of n_states_soc
  integer, allocatable :: l_row_soc(:), l_col_soc(:)
  integer, allocatable :: my_row_soc(:), my_col_soc(:)

  ! ScaLAPACK quantities for n_basis_soc*n_basis_soc matrices
  integer :: block_size_soc_vec
  integer :: sc_desc_soc_vec(dlen_)
  integer :: mb_soc_vec, nb_soc_vec
  integer :: mxld_soc_vec, mxcol_soc_vec
  integer :: n_my_rows_soc_vec, n_my_cols_soc_vec
  ! Dimensions of n_basis_soc
  integer, allocatable :: l_row_soc_vec(:), l_col_soc_vec(:)
  integer, allocatable :: my_row_soc_vec(:), my_col_soc_vec(:)
contains
  !****f* scalapack_soc/initialize_scalapack_soc
  !*  NAME
  !*    initialize_scalapack_soc
  !*  SYNOPSIS
  subroutine initialize_scalapack_soc ( n_states_soc, n_basis_soc )
  !*  PURPOSE
  !*    Initializes the ScaLAPACK distributions for SOC.
  !*  USES
    use localorb_io, only: localorb_info
    use aims_memory_tracking, only: aims_allocate
    implicit none
  !*  ARGUMENTS
    integer, intent(in) :: n_states_soc
    integer, intent(in) :: n_basis_soc
  !*  INPUTS
  !*    o n_states_soc - Number of states in second-variational window for SOC
  !*    o n_basis_soc  - Size of basis set for SOC
  !*  OUTPUTS
  !*    None
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  HISTORY
  !*    Added September 2017
  !*  SOURCE
    character*300 :: info_str
    write(info_str, '(2X,A)') "Creating ScaLAPACK distribution for matrices&
                              & with leading dimension n_states_soc"
    call localorb_info( info_str )
    ! Allocate global-local indexing arrays for matrices with dimensions
    ! n_states_soc
    call aims_allocate(l_row_soc, n_states_soc, "l_row_soc")
    call aims_allocate(l_col_soc, n_states_soc, "l_col_soc")
    call aims_allocate(my_row_soc, n_states_soc, "my_row_soc")
    call aims_allocate(my_col_soc, n_states_soc, "my_col_soc")
    block_size_soc = 0
    call initialize_scalapack_descriptor_matrix( &
         n_states_soc, n_states_soc, block_size_soc, sc_desc_soc,  mb_soc, &
         nb_soc, mxld_soc, mxcol_soc, n_my_rows_soc, n_my_cols_soc, l_row_soc, &
         l_col_soc, my_row_soc, my_col_soc )
    write(info_str, *)
    call localorb_info( info_str )

    ! Allocate global-local indexing arrays for matrices with dimensions
    ! n_basis_soc
    write(info_str, '(2X,A)') "Creating ScaLAPACK distribution for matrices&
                              & with leading dimension n_basis_soc"
    call localorb_info( info_str )
    call aims_allocate(l_row_soc_vec, n_basis_soc, "l_row_soc_vec")
    call aims_allocate(l_col_soc_vec, n_basis_soc, "l_col_soc_vec")
    call aims_allocate(my_row_soc_vec, n_basis_soc, "my_row_soc_vec")
    call aims_allocate(my_col_soc_vec, n_basis_soc, "my_col_soc_vec")
    block_size_soc_vec = 0
    call initialize_scalapack_descriptor_matrix( &
         n_basis_soc, n_basis_soc, block_size_soc_vec, sc_desc_soc_vec, &
         mb_soc_vec, nb_soc_vec, mxld_soc_vec, mxcol_soc_vec, &
         n_my_rows_soc_vec, n_my_cols_soc_vec, l_row_soc_vec, l_col_soc_vec, &
         my_row_soc_vec, my_col_soc_vec )
    write(info_str, *)
    call localorb_info( info_str )
  end subroutine initialize_scalapack_soc
  !******

  !****f* scalapack_soc/finalize_scalapack_soc
  !*  NAME
  !*    finalize_scalapack_soc
  !*  SYNOPSIS
  subroutine finalize_scalapack_soc
  !*  PURPOSE
  !*    Deallocates arrays related to the ScaLAPACK distribution for SOC.
  !*  USES
    use aims_memory_tracking, only: aims_deallocate
    implicit none
  !*  ARGUMENTS
  !*    None, modifies modular variables
  !*  INPUTS
  !*    None
  !*  OUTPUTS
  !*    None
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  HISTORY
  !*    Added September 2017
  !*  SOURCE
    if(allocated(l_row_soc)) call aims_deallocate(l_row_soc, "l_row_soc")
    if(allocated(l_col_soc)) call aims_deallocate(l_col_soc, "l_col_soc")
    if(allocated(my_row_soc)) call aims_deallocate(my_row_soc, "my_row_soc")
    if(allocated(my_col_soc)) call aims_deallocate(my_col_soc, "my_col_soc")
    if(allocated(l_row_soc_vec)) call aims_deallocate(l_row_soc_vec, &
         "l_row_soc_vec")
    if(allocated(l_col_soc_vec)) call aims_deallocate(l_col_soc_vec, &
         "l_col_soc_vec")
    if(allocated(my_row_soc_vec)) call aims_deallocate(my_row_soc_vec, &
         "my_row_soc_vec")
    if(allocated(my_col_soc_vec)) call aims_deallocate(my_col_soc_vec, &
         "my_col_soc_vec")
  end subroutine finalize_scalapack_soc
  !******
end module scalapack_soc
!******
