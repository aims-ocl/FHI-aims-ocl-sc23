!****s* FHI-aims/generate_Lsp_indexing
!  NAME
!    generate_Lsp_indexing
!  SYNOPSIS

subroutine generate_Lsp_indexing(n_bas_fns, basfn_l, basfn_species, &
&                            max_bas_L, max_n_bas_fnLsp, n_bas, &
&                            Lsp2n_bas_fnLsp, Lsp2bas_fn, Lsp2bas_sp, &
&                            atom2bas_off, sp2n_bas_sp)

  !  PURPOSE
  !
  !    From a given list of radial parts, generate the bas_* arrays as well
  !    as the Lsp arrays
  !
  !  USES

  use dimensions,      only : n_species, n_atoms
  use geometry,        only : species
  use runtime_choices, only : flag_rel, REL_x2c, REL_4c_dks
  use mpi_tasks,       only : aims_stop
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_bas_fns
  integer, intent(IN) :: basfn_l(n_bas_fns)
  integer, intent(IN) :: basfn_species(n_bas_fns)
  integer, intent(IN) :: max_bas_L
  integer, intent(IN) :: max_n_bas_fnLsp
  integer, intent(IN) :: n_bas

  integer, intent(OUT) :: Lsp2n_bas_fnLsp(0:max_bas_L, n_species)
  integer, intent(OUT) :: Lsp2bas_fn(max_n_bas_fnLsp, 0:max_bas_L, n_species)
  integer, intent(OUT) :: Lsp2bas_sp(max_n_bas_fnLsp, 0:max_bas_L, n_species)
  integer, intent(OUT) :: atom2bas_off(n_atoms)
  integer, intent(OUT) :: sp2n_bas_sp(n_species)

  !  INPUTS
  !    o basfn_l -- i_bas_fn -> angular momentum quantum number L
  !    o basfn_species -- i_bas_fn -> i_species
  !    o max_bas_L -- Maximum angular momentum
  !    o max_n_bas_fnLsp -- Maximum number of radial parts within L-channel
  !    o n_bas -- Number of bas functions in system (or 0-cell)
  !  OUTPUTS
  !    o Lsp2n_bas_fnLsp -- L, i_species -> n_fnLsp
  !    o Lsp2bas_fn -- i_fnLsp, L, i_species -> i_bas_fn
  !    o Lsp2bas_sp -- i_fnLsp, L, i_species -> i_bas_sp
  !    o atom2bas_off -- i_atom -> atom_off
  !        i_bas = atom_off + i_bas_sp + M
  !    o sp2n_bas_sp -- i_species -> Number of bas functions per atom
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
  !    Release version, FHI-aims (2010).
  !  SOURCE

  integer :: i_fnLsp
  integer :: n_bas_uptonow, atom_off
  integer :: i_atom, i_species, i_bas, i_bas_sp, i_bas_fn, L, M
  character(*), parameter :: func = 'generate_Lsp_indexing'

  ! --- Set Lsp arrays

  Lsp2n_bas_fnLsp = 0
  Lsp2bas_fn = 0
  sp2n_bas_sp = 0
  do i_bas_fn = 1, n_bas_fns
     i_species = basfn_species(i_bas_fn)
     L = basfn_l(i_bas_fn)
     Lsp2n_bas_fnLsp(L, i_species) = Lsp2n_bas_fnLsp(L, i_species) + 1
     i_fnLsp = Lsp2n_bas_fnLsp(L, i_species)
     if (i_fnLsp > max_n_bas_fnLsp) then
        call aims_stop('max_n_bas_fnLsp too small', func)
     else if (any(Lsp2bas_fn == i_bas_fn)) then
        call aims_stop('Double usage of i_bas_fn', func)
     end if
     Lsp2bas_fn(i_fnLsp, L, i_species) = i_bas_fn
     Lsp2bas_sp(i_fnLsp, L, i_species) = sp2n_bas_sp(i_species) + L + 1
     sp2n_bas_sp(i_species) = sp2n_bas_sp(i_species) + (2*L+1)
     ! (Rundong) I don't modify Lsp2bas_sp here, because I don't know what it is
     ! used for in the later code. This variable may need to be modified in the
     ! future, for a ralativistic case.
     if((flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks))then
     endif
  end do
  if (sum(Lsp2n_bas_fnLsp) /= n_bas_fns) then
     call aims_stop('sum(Lsp2n_bas_fnLsp) mismatch', func)
  end if

  ! --- Set atomic offset

  n_bas_uptonow = 0
  do i_atom = 1, n_atoms
     i_species = species(i_atom)
     atom2bas_off(i_atom) = n_bas_uptonow
     n_bas_uptonow = n_bas_uptonow + sp2n_bas_sp(i_species)
  end do

end subroutine generate_Lsp_indexing
!******
