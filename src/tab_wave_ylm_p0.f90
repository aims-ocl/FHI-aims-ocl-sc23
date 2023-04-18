!****s* FHI-aims/tab_wave_ylm_p0
!  NAME
!   tab_wave_ylm_p0
!  SYNOPSIS

subroutine tab_wave_ylm_p0 &
     ( n_compute_atoms, atom_index, & 
     trigonom_tab, basis_l_max, l_ylm_max, &
     ylm_tab &
     )

!  PURPOSE
!  Subroutine tab_wave_ylm computes the Y_lm functions needed 
!  for wave functions between a 3D cartesian point and all atoms
!
!  USES

  use dimensions
  use pbc_lists
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer :: n_compute_atoms
  integer :: atom_index(n_compute_atoms)
  real*8 trigonom_tab ( 4, n_compute_atoms )
  integer basis_l_max (n_species)
  integer l_ylm_max
  real*8 ylm_tab( (l_ylm_max+1)**2, n_compute_atoms )

!  INPUTS
!    o n_compute_atoms -- number of relevant atoms
!    o atom_index -- list of relevant atoms
!    o trigonom_tab -- trigonometric functions from subroutine tab_trigonom_p0 
!    o basis_l_max -- maximum l index of basis functions
!    o l_ylm_max -- maximum l index of wave functions
!    
!  OUTPUT
!    o ylm_tab -- Y_lm functions
!
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
!



  !  counters

  integer i_atom
  integer l_counter

  !  begin work

  do i_atom = 1, n_compute_atoms, 1

     !       bootstrap table of ylm's
     call increment_ylm &
          ( trigonom_tab(1,i_atom), trigonom_tab(2,i_atom),  &
          trigonom_tab(3,i_atom), trigonom_tab(4,i_atom),  &
          0, basis_l_max(species_center(atom_index(i_atom))), &
          ylm_tab(1,i_atom) ) 

  enddo

  if(flag_rel.eq.REL_4c_dks.or.flag_rel.eq.REL_x2c)then
    ! For fully relativistic cases, the convention of spherical harmonics should
    ! be converted to fit the scalar to spinar integration transformation
     do i_atom = 1, n_compute_atoms, 1
        if(l_ylm_max.ge.1)then ! for p orbital
           ylm_tab(4,i_atom) = -ylm_tab(4,i_atom)
        endif
        if(l_ylm_max.ge.2)then ! for d orbital
           ylm_tab(8,i_atom) = -ylm_tab(8,i_atom)
        endif
        if(l_ylm_max.ge.3)then ! for f orbital
           ylm_tab(14,i_atom) = -ylm_tab(14,i_atom)
           ylm_tab(16,i_atom) = -ylm_tab(16,i_atom)
        endif
        if(l_ylm_max.ge.4)then ! for g orbital
           ylm_tab(22,i_atom) = -ylm_tab(22,i_atom)
           ylm_tab(24,i_atom) = -ylm_tab(24,i_atom)
        endif
        if(l_ylm_max.ge.5)then ! for h orbital
           ylm_tab(32,i_atom) = -ylm_tab(32,i_atom)
           ylm_tab(34,i_atom) = -ylm_tab(34,i_atom)
           ylm_tab(36,i_atom) = -ylm_tab(36,i_atom)
        endif
        if(l_ylm_max.ge.6)then ! for the upgraded orbital of h
           ylm_tab(44,i_atom) = -ylm_tab(44,i_atom)
           ylm_tab(46,i_atom) = -ylm_tab(46,i_atom)
           ylm_tab(48,i_atom) = -ylm_tab(48,i_atom)
        endif
     enddo
  endif

  !  that's all folks

  return
end subroutine tab_wave_ylm_p0

!----------------------------------------------------------------------
!******
