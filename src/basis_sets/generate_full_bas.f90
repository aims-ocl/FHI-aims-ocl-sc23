!****s* FHI-aims/generate_full_bas
!  NAME
!    generate_full_bas
!  SYNOPSIS

subroutine generate_full_bas(n_bas_fns, basfn_l, basfn_species, &
&                            max_bas_L, max_bas_n_fnLsp, n_bas, &
&                            bas_atom, bas_l, bas_m, bas_fn, &
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
  integer, intent(IN) :: max_bas_n_fnLsp
  integer, intent(IN) :: n_bas

  integer, intent(OUT) :: bas_atom(n_bas)
  integer, intent(OUT) :: bas_l(n_bas)
  integer, intent(OUT) :: bas_m(n_bas)
  integer, intent(OUT) :: bas_fn(n_bas)

  integer, intent(OUT) :: Lsp2n_bas_fnLsp(0:max_bas_L, n_species)
  integer, intent(OUT) :: Lsp2bas_fn(max_bas_n_fnLsp, 0:max_bas_L, n_species)
  integer, intent(OUT) :: Lsp2bas_sp(max_bas_n_fnLsp, 0:max_bas_L, n_species)
  integer, intent(OUT) :: atom2bas_off(n_atoms)
  integer, intent(OUT) :: sp2n_bas_sp(n_species)

  !  INPUTS
  !    o basfn_l -- i_bas_fn -> angular momentum quantum number L
  !    o basfn_species -- i_bas_fn -> i_species
  !    o max_bas_L -- Maximum angular momentum
  !    o max_bas_n_fnLsp -- Maximum number of radial parts within L-channel
  !    o n_bas -- Number of bas functions in system (or 0-cell)
  !  OUTPUTS
  !    o bas_atom -- i_bas -> i_atom
  !    o bas_l -- i_bas -> L
  !    o bas_m -- i_bas -> M
  !    o bas_fn -- i_bas -> i_bas_fn
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

  integer :: i_atom, i_species, L, i_fnLsp, M
  integer :: atom_off, i_bas_fn, i_bas_sp, i_bas
  character(*), parameter :: func = 'generate_full_bas'

  ! --- Set bas_ arrays

  call generate_bas_indexing(n_bas_fns, basfn_l, basfn_species, n_bas, &
  &                          bas_atom, bas_l, bas_m, bas_fn)

  ! --- Set Lsp arrays

  call generate_Lsp_indexing(n_bas_fns, basfn_l, basfn_species, &
  &                          max_bas_L, max_bas_n_fnLsp, n_bas, &
  &                          Lsp2n_bas_fnLsp, Lsp2bas_fn, Lsp2bas_sp, &
  &                          atom2bas_off, sp2n_bas_sp)


  ! --- Check consistency

  ! The following loops check consistency but also serve as reference
  ! documentation of how these arrays can be used.
  do i_atom = 1, n_atoms
     i_species = species(i_atom)
     atom_off = atom2bas_off(i_atom)
     do L = 0, max_bas_L
        do i_fnLsp = 1, Lsp2n_bas_fnLsp(L, i_species)
           i_bas_fn = Lsp2bas_fn(i_fnLsp, L, i_species)
           i_bas_sp = Lsp2bas_sp(i_fnLsp, L, i_species)
           if((flag_rel.eq.REL_x2c .or.flag_rel.eq.REL_4c_dks) .and. L.eq.0)then
              i_bas = atom_off + i_bas_sp - L
              if (bas_atom(i_bas) /= i_atom) then
                 call aims_stop('bas_atom mismatch', func)
              else if (bas_l(i_bas) /= L) then
                 call aims_stop('bas_l mismatch', func)
              else if (bas_fn(i_bas) /= i_bas_fn) then
                 call aims_stop('bas_fn mismatch', func)
              end if
           elseif((flag_rel.eq.REL_x2c .or.flag_rel.eq.REL_4c_dks) .and. L.ne.0)then
              do M = 1, 4*L+2
                 i_bas = atom_off + i_bas_sp -L
                 if (bas_atom(i_bas) /= i_atom) then
                    call aims_stop('bas_atom mismatch', func)
                 else if (bas_l(i_bas) /= L) then
                    call aims_stop('bas_l mismatch', func)
                 else if (bas_fn(i_bas) /= i_bas_fn) then
                    call aims_stop('bas_fn mismatch', func)
                 end if
              end do
           else
              do M = -L, L
                 i_bas = atom_off + i_bas_sp + M
                 if (bas_atom(i_bas) /= i_atom) then
                    call aims_stop('bas_atom mismatch', func)
                 else if (bas_l(i_bas) /= L) then
                    call aims_stop('bas_l mismatch', func)
                 else if (bas_m(i_bas) /= M) then
                    call aims_stop('bas_m mismatch', func)
                 else if (bas_fn(i_bas) /= i_bas_fn) then
                    call aims_stop('bas_fn mismatch', func)
                 end if
              end do
           endif
        end do
     end do
  end do

end subroutine generate_full_bas
!******
