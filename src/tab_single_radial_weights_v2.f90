!****s* FHI-aims/tab_single_radial_weights_v2
!  NAME
!   tab_single_radial_weights_v2
!  SYNOPSIS

      subroutine tab_single_radial_weights_v2 &
      ( atom_in, dist_tab_in, i_r, &
        log_weight, radial_weight &
      )

!  PURPOSE
!  Subroutine tab_radial_weights tabulates radial weights di/dr which are
!  required to compute the derivative f'(r) of a tabulated splined function f(i)
!  on a real-space grid r(i). (where i is the index of the grid point)
!
!  USES

        use dimensions
        use geometry
        use grids
        implicit none

!  ARGUMENTS

      integer :: atom_in
      real*8 dist_tab_in
      real*8 i_r
      real*8 log_weight
      real*8 radial_weight

!  INPUTS
!   o atom_in -- atom which weights we are calculating      
!   o dist_tab_in -- distance to atom
!   o i_r -- the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!    
!  OUTPUT
!   o log_weight -- integration weight for log grid
!   o radial_weight -- integration weight for radial grid
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


      integer :: i_species
      integer :: i_atom_in
      integer :: i_coord

!  begin work


        ! This is the actual index of the current atom i_atom_in
        ! in the list of all n_atoms atoms
      i_species = species(atom_in)

      radial_weight = &
          ( dble(n_radial(i_species)+1) / i_r )**2

        radial_weight = &
          i_r / (2.d0 * scale_radial(i_species) ) &
          * (radial_weight - 1)

        log_weight = &
          log_r_grid_inc(i_species) * dist_tab_in

        log_weight = 1.d0/log_weight

!  that's it

      end subroutine tab_single_radial_weights_v2
!----------------------------------------------------------------------
!******
