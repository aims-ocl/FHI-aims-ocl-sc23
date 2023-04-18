!****s* FHI-aims/generate_bas_indexing
!  NAME
!    generate_bas_indexing
!  SYNOPSIS

subroutine generate_bas_indexing(n_bas_fns, basfn_l, basfn_species, n_bas, &
&                            bas_atom, bas_l, bas_m, bas_fn)

  !  PURPOSE
  !
  !    From a given list of radial parts, generate the bas_* arrays.
  !
  !  USES

  use dimensions,      only : n_atoms
  use geometry,        only : species
  use runtime_choices, only : flag_rel, REL_x2c, REL_4c_dks
  use mpi_tasks,       only : aims_stop
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_bas_fns
  integer, intent(IN) :: basfn_l(n_bas_fns)
  integer, intent(IN) :: basfn_species(n_bas_fns)
  integer, intent(IN) :: n_bas

  integer, intent(OUT) :: bas_atom(n_bas)
  integer, intent(OUT) :: bas_l(n_bas)
  integer, intent(OUT) :: bas_m(n_bas)
  integer, intent(OUT) :: bas_fn(n_bas)

  !  INPUTS
  !    o basfn_l -- i_bas_fn -> angular momentum quantum number L
  !    o basfn_species -- i_bas_fn -> i_species
  !    o n_bas -- Number of bas functions in system (or 0-cell)
  !  OUTPUTS
  !    o bas_atom -- i_bas -> i_atom
  !    o bas_l -- i_bas -> L
  !    o bas_m -- i_bas -> M
  !    o bas_fn -- i_bas -> i_bas_fn
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

  integer :: n_bas_uptonow
  integer :: i_atom, i_species, i_bas_fn, L, M
  character(*), parameter :: func = 'generate_bas_indexing'

  ! --- Set bas_ arrays

  n_bas_uptonow = 0
  do i_atom = 1, n_atoms
     i_species = species(i_atom)
     do i_bas_fn = 1, n_bas_fns
        if (basfn_species(i_bas_fn) == i_species) then
           L = basfn_l(i_bas_fn)
           do M = -L, L
              n_bas_uptonow = n_bas_uptonow + 1
              if (n_bas_uptonow > n_bas) then
                 call aims_stop('Consistenty error: Too many bas.', func)
              end if
              bas_atom(n_bas_uptonow) = i_atom
              bas_l(n_bas_uptonow) = L
              bas_m(n_bas_uptonow) = M
              bas_fn(n_bas_uptonow) = i_bas_fn
              ! (Rundong) The ralativistic part needs to be further modified for j-adapted cases. The current code is 
              ! only written like this to avoid bugs.
              if((flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) .and. L.ne.0)then
                n_bas_uptonow = n_bas_uptonow + 1
                bas_atom(n_bas_uptonow) = i_atom
                bas_l(n_bas_uptonow) = L
                bas_m(n_bas_uptonow) = M ! this variable makes no sense currently
                bas_fn(n_bas_uptonow) = i_bas_fn
              endif
           end do
        end if
     end do
  end do
  if (n_bas_uptonow /= n_bas) then
     call aims_stop('Consistenty error: Too few bas funcs.', func)
  endif

end subroutine generate_bas_indexing
!******
