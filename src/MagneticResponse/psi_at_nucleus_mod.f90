!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Provides arrays that contain the wavefunction values at the
!!  nuclei.
!!
!!  Usage: if the wavefunction value is required at atom i_atom, then
!!         inside an integration routine use
!!
!!    if (basis_atom(i_basis(i)) == i_atom) &
!!         & X = psi_at_nucleus(basis_fn(i_basis(i))),
!!
!!  where i iterates over n_compute (number of nonzero orbitals for
!!  the given batch) and X is the wavefunction value at the nucleus.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
module psi_at_nucleus_mod

  use aims_memory_tracking, only: aims_allocate
  use dimensions,           only: n_max_ind_fns, n_species
  use tools,                only: safe_deallocate
  use types,                only: dp

  implicit none

  public

  ! These arrays are constructed in the respective subroutines. For
  ! instance, psi_at_nucleus_hydro is computed in
  ! get_hydrogenic_basis_fns.
  real(dp), allocatable :: psi_at_nucleus_atomic(:,:)
  real(dp), allocatable :: psi_at_nucleus_hydro(:,:)
  real(dp), allocatable :: psi_at_nucleus_gauss(:,:)
  real(dp), allocatable :: psi_at_nucleus_ionic(:,:)
  ! This array contains the above arrays in a compact form and is also
  ! normalized. It is built up in shrink_fixed_basis_phi_thresh and it
  ! is the one actually used in integration.
  real(dp), allocatable :: psi_at_nucleus(:)

contains
  ! This is called in read_control.
  subroutine initialize_psi_at_nucleus()
    character(*), parameter :: &
         & THIS_SUB = 'psi_at_nucleus_mod::initialize_psi_at_nucleus::'
    call aims_allocate(psi_at_nucleus_atomic, n_max_ind_fns, n_species, &
         & THIS_SUB//'psi_at_nucleus_atomic')
    psi_at_nucleus_atomic = 0d0
    call aims_allocate(psi_at_nucleus_hydro, n_max_ind_fns, n_species, &
         & THIS_SUB//'psi_at_nucleus_hydro')
    psi_at_nucleus_hydro = 0d0
    call aims_allocate(psi_at_nucleus_ionic, n_max_ind_fns, n_species, &
         & THIS_SUB//'psi_at_nucleus_ionic')
    psi_at_nucleus_ionic = 0d0
    call aims_allocate(psi_at_nucleus_gauss, n_max_ind_fns, n_species, &
         & THIS_SUB//'psi_at_nucleus_gauss')
    psi_at_nucleus_gauss = 0d0
  end subroutine initialize_psi_at_nucleus

  ! This is called in final_deallocations.
  subroutine cleanup_psi_at_nucleus()
    call safe_deallocate(psi_at_nucleus_atomic)
    call safe_deallocate(psi_at_nucleus_hydro)
    call safe_deallocate(psi_at_nucleus_gauss)
    call safe_deallocate(psi_at_nucleus_ionic)
    call safe_deallocate(psi_at_nucleus)
  end subroutine cleanup_psi_at_nucleus
end module psi_at_nucleus_mod
