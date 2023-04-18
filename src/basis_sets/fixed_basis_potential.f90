!****s* FHI-aims/fixed_basis_potential
!  NAME
!   fixed_basis_potential
!  SYNOPSIS

      subroutine fixed_basis_potential ( )

!  PURPOSE
!   Subroutine fixed_basis_potential pastes together the basis-defining radial
!   potential from
!    (i)  the free-atom effective potential, inside r_cutoff
!    (ii) a rapidly rising cutoff potential, added outside r_cutoff
!  USES
      use dimensions,   only : n_species
      use grids,        only : n_grid, r_grid
      use species_data, only : basis_potential, include_min_basis, r_cutoff, &
                               w_cutoff, cutoff_type, scale_cutoff
      use free_atoms,   only : free_potential
      use localorb_io,  only : use_unit, localorb_info
      implicit none
! INPUT
! none
! OUTPUT
! none
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

      character*100 :: info_str
      integer       :: i_species, i_grid
      real*8        :: cutoff_pot

      write (info_str,'(2X,A,A)') &
        "Adding cutoff potential to free-atom ", &
        "effective potential."

      call localorb_info('',use_unit)
      call localorb_info(info_str,use_unit,'(A)')

      do i_species = 1, n_species,1
      if (include_min_basis(i_species)) then

        do i_grid = 1, n_grid(i_species), 1

          basis_potential (i_grid, i_species) = &
          free_potential (i_grid, i_species)

          if ( (r_cutoff(i_species)+w_cutoff(i_species)).gt.0. ) then
            basis_potential (i_grid, i_species) = &
            basis_potential (i_grid, i_species) + &
            cutoff_pot &
            ( r_grid(i_grid,i_species), cutoff_type(i_species), &
              r_cutoff(i_species), w_cutoff(i_species), &
              scale_cutoff(i_species) )
          end if

        enddo

      end if
      enddo

      return
    end subroutine fixed_basis_potential
!******
