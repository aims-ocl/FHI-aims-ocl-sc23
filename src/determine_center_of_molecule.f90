!****s* FHI-aims/determine_center_of_molecule
!  NAME
!   determine_center_of_molecule
!  SYNOPSIS

      subroutine determine_center_of_molecule &
      ( coord_of_center, atom_coord_wrt_center)

!  PURPOSE
!  Determine the geomerical center of a given molecule.
!  If the molecule has an inversion symmetry, then find out the
!  inversion center, and evalute the relative coordination of
!  all atoms with repect to this center. Otherwise print out
!  an error message stating that there is no inversion center.
!
!  USES

      use dimensions
      use basis
      use geometry
      use constants

      implicit none

!  ARGUMENTS
      real*8 coord_of_center(3)
      real*8 atom_coord_wrt_center(3,n_atoms)

!  INPUTS 
!   o none  
!  OUTPUTS
!   o coord_of_center -- the coordinates of the inversion center
!   o atom_coord_wrt_center -- the relative position of all the atoms    
!           with respect to the inversion center 
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

!  local variables

!  counters

      integer i_atom


      coord_of_center(:) = 0
      do i_atom = 1, n_atoms, 1 
           coord_of_center(:) = coord_of_center(:) + coords(:,i_atom)
      enddo
      coord_of_center(:) = coord_of_center(:)/dble(n_atoms) 

      do i_atom = 1, n_atoms, 1 
         atom_coord_wrt_center(:,i_atom) = coords(:,i_atom) -  &
                                              coord_of_center(:)
      enddo


!      write(use_unit,'(2X, A50)') "Inversion center of the molecule :"
!      write(use_unit,'(2X,3f16.3)') coord_of_center(:)*bohr
!      write(use_unit,'(2X, A80)') &
!          "Relative poisiton of  the individual atoms :"
!      do i_atom = 1, n_atoms, 1
!       write(use_unit,'(2X,4f16.3)') atom_coord_wrt_center(:,i_atom)*bohr, &
!        sqrt(dot_product(atom_coord_wrt_center(:,i_atom),atom_coord_wrt_center(:,i_atom)))*bohr
!      enddo

      return

      end subroutine determine_center_of_molecule
!******
