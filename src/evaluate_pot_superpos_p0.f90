!****s* FHI-aims/evaluate_pot_superpos_p0
!  NAME
!    evaluate_pot_superpos_p0
!  SYNOPSIS

      subroutine evaluate_pot_superpos_p0 &
      ( i_r, potential, n_atom_list, atom_list &
      )

!  PURPOSE
!  Subroutine evaluate_pot_superpos tabulates a superposition of free-atom
!  potentials. This is used in ZORA potential calculations.
! 
!  USES

      use dimensions
      use runtime_choices
      use grids
      use pbc_lists
      use spline
      use free_atoms
      use constants
      implicit none

!  ARGUMENTS

      integer:: n_atom_list
      integer, dimension(n_atom_list) :: atom_list
      real*8,  dimension(n_atom_list) :: i_r
      real*8, dimension(n_spin) :: potential

!  INPUTS
!   o n_atom_list -- number of relevant atoms
!   o atom_list -- list of relevant atoms
!   o i_r -- tabulates the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!  OUTPUT
!   o potential -- a a superposition of free-atom  potentials
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




!  local variables

!     counters

      integer :: i_center_2, i_center_L
      integer :: i_spin

!  begin work

!     only superposition of free-atom potentials

      potential = 0.d0

!     Contributions from all atoms

        do i_center_L = 1, n_atom_list, 1
           i_center_2 = atom_list(i_center_L)
           
          do i_spin = 1,n_spin,1

            potential(i_spin) =  &
            potential(i_spin) +  &
            val_spline &
            ( i_r(i_center_L), free_potential_spl(1,1,species_center(i_center_2)), &
              n_grid(species_center(i_center_2)) )

          enddo

        enddo

      end subroutine evaluate_pot_superpos_p0
!---------------------------------------------------------------------
!******
