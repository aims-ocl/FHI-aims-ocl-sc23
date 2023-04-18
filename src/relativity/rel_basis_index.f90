
!(Rundong) This is a relativistic counterpart of /src/basis_sets/generate_full_bas.f90
subroutine generate_full_bas_rel(n_bas_fns, basfn_n, basfn_l, basfn_k, basfn_species, &
  max_bas_L, max_bas_n_fnLsp, n_bas, n_bas_small, bas_atom, bas_small_atom, &
  bas_l, bas_small_l, bas_k, bas_small_k, bas_m, bas_small_m, bas_fn, bas_small_fn, &
  Lsp2n_bas_fnLsp, Lsp2bas_fn, Lsp2bas_sp, atom2bas_off, sp2n_bas_sp)

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
  integer, intent(IN) :: basfn_n(n_bas_fns), basfn_l(n_bas_fns), basfn_k(n_bas_fns)
  integer, intent(IN) :: basfn_species(n_bas_fns)
  integer, intent(IN) :: max_bas_L
  integer, intent(IN) :: max_bas_n_fnLsp
  integer, intent(IN) :: n_bas, n_bas_small

  integer, intent(OUT) :: bas_atom(n_bas), bas_small_atom(n_bas_small)
  integer, intent(OUT) :: bas_l(n_bas), bas_small_l(n_bas_small)
  integer, intent(OUT) :: bas_k(n_bas), bas_small_k(n_bas_small)
  integer, intent(OUT) :: bas_m(n_bas), bas_small_m(n_bas_small)
  integer, intent(OUT) :: bas_fn(n_bas), bas_small_fn(n_bas_small)

  integer, intent(OUT) :: Lsp2n_bas_fnLsp(0:max_bas_L, n_species)
  integer, intent(OUT) :: Lsp2bas_fn(max_bas_n_fnLsp, 0:max_bas_L, n_species)
  integer, intent(OUT) :: Lsp2bas_sp(max_bas_n_fnLsp, 0:max_bas_L, n_species)
  integer, intent(OUT) :: atom2bas_off(n_atoms)
  integer, intent(OUT) :: sp2n_bas_sp(n_species)

  !  INPUTS
  !    o basfn_l -- i_bas_fn -> angular momentum quantum number L
  !    o basfn_k -- i_bas_fn -> quantum number kappa for relativistic basis
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

  call generate_bas_indexing_rel(n_bas_fns, basfn_n, basfn_l, basfn_k, basfn_species, n_bas, n_bas_small, &
       bas_atom, bas_small_atom, bas_l, bas_small_l, bas_k, bas_small_k, bas_m, bas_small_m, bas_fn, bas_small_fn)

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
           if((flag_rel.eq.REL_x2c .or.flag_rel.eq.REL_4c_dks))then
         !    i_bas = atom_off + i_bas_sp - L
         !    if (bas_atom(i_bas) /= i_atom) then
         !       call aims_stop('bas_atom mismatch', func)
         !    else if (bas_l(i_bas) /= L) then
         !       call aims_stop('bas_l mismatch', func)
         !    else if (bas_fn(i_bas) /= i_bas_fn) then
         !       call aims_stop('bas_fn mismatch', func)
         !    end if
         ! elseif((flag_rel.eq.REL_x2c .or.flag_rel.eq.REL_4c_dks) .and. L.ne.0)then
         !    do M = 1, 4*L+2
         !       i_bas = atom_off + i_bas_sp -L
         !       if (bas_atom(i_bas) /= i_atom) then
         !          call aims_stop('bas_atom mismatch', func)
         !       else if (bas_l(i_bas) /= L) then
         !          call aims_stop('bas_l mismatch', func)
         !       else if (bas_fn(i_bas) /= i_bas_fn) then
         !          call aims_stop('bas_fn mismatch', func)
         !       end if
         !    end do
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

end subroutine generate_full_bas_rel

!(Rundong) This is a relativistic counterpart of /src/basis_sets/generate_bas_indexing.f90
subroutine generate_bas_indexing_rel(n_bas_fns, basfn_n, basfn_l, basfn_k, basfn_species, n_bas, n_bas_small, &
  bas_atom, bas_small_atom, bas_l, bas_small_l, bas_k, bas_small_k, bas_m, bas_small_m, bas_fn, bas_small_fn)

  !  PURPOSE
  !
  !    From a given list of radial parts, generate the bas_* arrays.
  !
  !  USES

  use dimensions,      only : n_atoms
  use geometry,        only : species
  use runtime_choices, only : flag_rel, REL_x2c, REL_4c_dks
  use mpi_tasks,       only : aims_stop
  use localorb_io,     only : use_unit
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_bas_fns   ! number of RADIAL basis functions
  integer, intent(IN) :: basfn_n(n_bas_fns), basfn_l(n_bas_fns), basfn_k(n_bas_fns)
  integer, intent(IN) :: basfn_species(n_bas_fns)
  integer, intent(IN) :: n_bas, n_bas_small

  integer, intent(OUT) :: bas_atom(n_bas), bas_small_atom(n_bas_small)
  integer, intent(OUT) :: bas_l(n_bas), bas_small_l(n_bas_small)
  integer, intent(OUT) :: bas_k(n_bas), bas_small_k(n_bas_small)
  integer, intent(OUT) :: bas_m(n_bas), bas_small_m(n_bas_small)
  integer, intent(OUT) :: bas_fn(n_bas), bas_small_fn(n_bas_small)

  !  INPUTS
  !    o basfn_l -- i_bas_fn -> angular momentum quantum number L
  !    o basfn_k -- i_bas_fn -> quantum number kappa for relativistic basis
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

  integer :: bas_n(n_bas), bas_small_n(n_bas_small) ! For test only
  integer :: n_bas_uptonow, n_bas_uptonow_s
  integer :: i_atom, i_species, i_bas_fn, N, L, M, J
  character(*), parameter :: func = 'generate_bas_indexing'

  ! --- Set bas_ arrays

  n_bas_uptonow = 0; n_bas_uptonow_s = 0
  do i_atom = 1, n_atoms
     i_species = species(i_atom)
     do i_bas_fn = 1, n_bas_fns
        if (basfn_species(i_bas_fn) == i_species) then
           N = basfn_n(i_bas_fn)
           L = basfn_l(i_bas_fn)
           if((flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks))then

               if(basfn_k(i_bas_fn).eq.L)then ! spin down
                   J=2*L-1
               elseif(basfn_k(i_bas_fn).eq.(-L-1))then ! spin up
                   J=2*L+1
               else
                   write(use_unit,*)'Error in /basis_sets/generate_bas_indexing: invalid kappa value.'
                   stop
               endif

               ! For large component scalar basis:
               if(L.eq.0)then ! s orbital
                  n_bas_uptonow = n_bas_uptonow + 1
                  bas_atom(n_bas_uptonow) = i_atom
                  bas_n(n_bas_uptonow) = N
                  bas_k(n_bas_uptonow) = basfn_k(i_bas_fn)
                  bas_l(n_bas_uptonow) = L
                  bas_m(n_bas_uptonow) = 0 
                  bas_fn(n_bas_uptonow) = i_bas_fn
               else
                  do M = -L, L
                    n_bas_uptonow = n_bas_uptonow + 1
                    bas_atom(n_bas_uptonow) = i_atom
                    bas_n(n_bas_uptonow) = N
                    bas_k(n_bas_uptonow) = basfn_k(i_bas_fn)
                    bas_l(n_bas_uptonow) = L
                    bas_m(n_bas_uptonow) = M 
                    bas_fn(n_bas_uptonow) = i_bas_fn
                  enddo
               endif

               ! For small component scalar basis:
               if(basfn_k(i_bas_fn).eq.L)then ! For a spin down case, the angular momentum is degraded.
                   do M= -(L-1), (L-1)
                      n_bas_uptonow_s = n_bas_uptonow_s + 1
                      bas_small_atom(n_bas_uptonow_s) = i_atom
                      bas_small_n(n_bas_uptonow_s) = N
                      bas_small_l(n_bas_uptonow_s) = L-1  ! Note this!
                      bas_small_k(n_bas_uptonow_s) = L-1  ! Note this! Both l and k are changed.
                      bas_small_m(n_bas_uptonow_s) = M
                      bas_small_fn(n_bas_uptonow_s) = i_bas_fn
                   enddo
               elseif(basfn_k(i_bas_fn).eq.(-L-1))then ! For a spin up case, the angular momentum is upgraded.
                   do M= -(L+1), (L+1)
                      n_bas_uptonow_s = n_bas_uptonow_s + 1
                      bas_small_atom(n_bas_uptonow_s) = i_atom
                      bas_small_n(n_bas_uptonow_s) = N
                      bas_small_l(n_bas_uptonow_s) = L+1  ! Note this!
                      bas_small_k(n_bas_uptonow_s) = -L-2 ! Note this! Both l and k are changed.
                      bas_small_m(n_bas_uptonow_s) = M 
                      bas_small_fn(n_bas_uptonow_s) = i_bas_fn
                   enddo
               endif

           else ! This "else" is left here only for reference as a non-rel code. In subroutine generate_bas_indexing_rel,
                ! the code will never go into this "else" branch.
               do M = -L, L
                  n_bas_uptonow = n_bas_uptonow + 1
                  if (n_bas_uptonow > n_bas) then
                     call aims_stop('Consistenty error: Too many bas.', func)
                  end if
                  bas_atom(n_bas_uptonow) = i_atom
                  bas_l(n_bas_uptonow) = L
                  bas_m(n_bas_uptonow) = M
                  bas_fn(n_bas_uptonow) = i_bas_fn
               end do
           endif
        end if
     end do
  end do

 !do i_bas_fn=1, n_bas_uptonow
 !   write(use_unit,"('bas_fn=',i3,3x,'bas_atom=',i3,3x,'bas_n=',i3,3x,'bas_l=',i3,3x,'bas_m=',i3,3x,'bas_k=',i3)")&
 !      bas_fn(i_bas_fn), bas_atom(i_bas_fn), bas_n(i_bas_fn), bas_l(i_bas_fn), bas_m(i_bas_fn), bas_k(i_bas_fn)
 !enddo

 !do i_bas_fn=1, n_bas_uptonow_s
 !   write(use_unit,"('small_fn=',i3,3x,'bas_atom=',i3,3x,'bas_n=',i3,3x,'bas_l=',i3,3x,'bas_m=',i3,3x,'bas_k=',i3)")&
 !      bas_small_fn(i_bas_fn), bas_small_atom(i_bas_fn), bas_small_n(i_bas_fn), bas_small_l(i_bas_fn), bas_small_m(i_bas_fn), bas_small_k(i_bas_fn)
 !enddo 

  if (n_bas_uptonow /= n_bas) then
     call aims_stop('Consistenty error: Too few large bas funcs.', func)
  endif
  if (n_bas_uptonow_s /= n_bas_small) then
     call aims_stop('Consistenty error: Too few small bas funcs.', func)
  endif

end subroutine generate_bas_indexing_rel
!



