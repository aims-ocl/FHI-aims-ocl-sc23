!****s* FHI-aims/free_atoms_out
!  NAME
!   free_atoms_out
!  SYNOPSIS
      subroutine free_atoms_out (potential, rho, n_wave, wave, wave_n, &
           wave_l, eigenval, descript )
! PURPOSE
!  Subroutine free_atoms_out is for debugging purposes only.
!  It's a wrapper around atom_out(), which takes data for an individual atom and
!  writes it out in a defined format.
! USES
      use constants,    only : pi
      use dimensions,   only : n_max_grid, n_species, n_max_ind_fns
      use grids,        only : n_grid, r_grid
      use species_data, only : species_name
      implicit none
!  ARGUMENTS
      real*8      :: potential(n_max_grid, n_species)
      real*8      :: rho(n_max_grid, n_species)
      integer     :: n_wave(n_species)
      real*8      :: wave(n_max_grid, n_species, n_max_ind_fns)
      integer     :: wave_n(n_species, n_max_ind_fns)
      integer     :: wave_l(n_species, n_max_ind_fns)
      real*8      :: eigenval(n_species, n_max_ind_fns)
      character*4 :: descript
! INPUTS
!  o potential -- total potential
!  o rho -- density
!  o n_wave -- number of wave functions
!  o wave -- wave functions
!  o wave_n -- quantum number n for all wave functions
!  o wave_l -- quantum number l for all wave functions
!  o eigenval -- eigenvalues
!  o descript -- file name descriptor for output
! OUTPUTS
!   none
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
!  * The suffix "_single" indicates that these quantities each refer
!    only to one individual atom, and there is no extra index (such as
!    i_species, i_atom)
!  * prefix dim_ denotes variable array dimensions

      real*8 four_pi_rsq_rho(n_max_grid)
      real*8 wave_single(n_max_grid,n_max_ind_fns)
      integer wave_n_single(n_max_ind_fns)
      integer wave_l_single(n_max_ind_fns)
      real*8 eigenval_single(n_max_ind_fns)

      integer dim_grid, dim_fns

!  counters

      integer i_species, i_grid, i_fn

!  begin work

      do i_species = 1, n_species, 1

!       first, copy everything which is species-dependent onto individual arrays

        do i_fn = 1, n_wave(i_species), 1
          wave_n_single(i_fn) = wave_n(i_species,i_fn)
          wave_l_single(i_fn) = wave_l(i_species,i_fn)
!test
!          write(use_unit,*) i_fn, wave_n_single(i_fn), wave_l_single(i_fn)
!test end
          eigenval_single(i_fn) = eigenval(i_species,i_fn)
          do i_grid = 1, n_grid(i_species),1
            wave_single(i_grid,i_fn) = wave(i_grid,i_species,i_fn)
          enddo
        enddo

        do i_grid = 1, n_grid(i_species), 1
          four_pi_rsq_rho(i_grid) = rho(i_grid, i_species) &
            * 4.d0 * pi * (r_grid(i_grid,i_species)**2.d0)
        enddo

!  next, set variable grids for atom_out

        dim_grid    = n_max_grid
        dim_fns = n_max_ind_fns

!test
!        do i_fn = 1, n_wave(i_species), 1
!          write(use_unit,*) i_fn, wave_n_single(i_fn), wave_l_single(i_fn)
!        enddo
!test end

!  now, call atom_out() to write out individual atoms

        call atom_out &
          ( species_name(i_species), n_grid(i_species), &
            r_grid(1,i_species), potential(1,i_species), &
            four_pi_rsq_rho, n_wave(i_species), wave_single, &
            wave_n_single, wave_l_single, eigenval_single, &
            descript, dim_grid, dim_fns )

      enddo

!  that's all folks

      return
    end subroutine free_atoms_out
!******
