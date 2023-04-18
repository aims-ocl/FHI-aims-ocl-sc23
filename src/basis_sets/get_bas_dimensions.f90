!****s* FHI-aims/get_bas_dimensions
!  NAME
!    get_bas_dimensions
!  SYNOPSIS

subroutine get_bas_dimensions(n_bas_fns, basfn_l, basfn_species, &
&                             max_bas_L, max_bas_n_fnLsp, n_bas)

  !  PURPOSE
  !
  !    For a given set of (product) basis functions, induce array dimension
  !    paramters.
  !
  !  USES

  use dimensions,   only : n_species
  use runtime_choices, only : flag_rel, REL_x2c, REL_4c_dks
  use species_data, only : atoms_in_structure
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_bas_fns
  integer, intent(IN) :: basfn_l(n_bas_fns)
  integer, intent(IN) :: basfn_species(n_bas_fns)
  integer, intent(OUT) :: max_bas_L
  integer, intent(OUT) :: max_bas_n_fnLsp
  integer, intent(OUT) :: n_bas

  !  INPUTS
  !    o n_bas_fns -- Number of distinct radial parts
  !    o basfn_l -- Angular momentum quantum numbers of radial parts
  !    o basfn_species -- Atomic species to which they belong to
  !  OUTPUTS
  !    o max_bas_L -- Highest overall angular momentum
  !    o max_bas_n_fnLsp -- Maximum number of radial parts
  !                         with identical L, i_species
  !    o n_bas -- Total number of basis functions in system (or 0-cell).
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  integer :: L, i_species, i_bas_fn, n_fnLsp
  character(*), parameter :: func = 'get_bas_dimensions'

  max_bas_L = maxval(basfn_l)
  
  max_bas_n_fnLsp = 0
  do i_species = 1, n_species
     do L = 0, max_bas_L
        n_fnLsp = count(basfn_l == L .and. basfn_species == i_species)
        max_bas_n_fnLsp = max(max_bas_n_fnLsp, n_fnLsp)

     end do
  end do

  n_bas = 0
  do i_bas_fn = 1, n_bas_fns
     i_species = basfn_species(i_bas_fn)
     L = basfn_l(i_bas_fn)
     n_bas = n_bas + (2*L+1) * atoms_in_structure(i_species)
  end do

end subroutine get_bas_dimensions
!******
