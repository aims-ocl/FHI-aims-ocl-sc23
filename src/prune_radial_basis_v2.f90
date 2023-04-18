!****s* FHI-aims/prune_radial_basis_v2
!  NAME
!   prune_radial_basis_v2
!  SYNOPSIS

subroutine prune_radial_basis_v2 &
     ( dist_tab_sq, &
     n_compute_atoms, atom_index, atom_index_inv, &
     n_compute_fns, i_basis_fns, i_basis_fns_inv, &
     i_atom_fns, spline_array_start, spline_array_end )

!  PURPOSE
!     Subroutine prune_basis_v2
!     reduces full set of species, atoms, basis functions and basis to
!     relevant ones only, for a given set of n_points integration points
!
!  USES
!
  use dimensions
  use basis
  use grids
  use geometry
  implicit none
    
!  ARGUMENTS

  real*8 dist_tab_sq(n_atoms)
  integer :: n_compute_atoms
  integer :: atom_index(n_atoms)
  integer :: atom_index_inv(n_atoms)
  
  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_atoms)
  integer :: i_atom_fns(n_basis_fns*n_atoms)
  integer :: i_basis_fns_inv(n_basis_fns,n_atoms)
  
  integer :: spline_array_start(n_atoms)
  integer :: spline_array_end(n_atoms)


!  INPUTS
!    o dist_tab_sq -- (distance to atoms)**2
!  OUTPUT
!    o n_compute_atoms --  number of relevant atoms 
!    o atom_index -- global indeces of relevant atoms
!    o atom_index_inv -- inverse index list of relevant atoms
!    o n_compute_fns -- total number of radial functions relevant at these points
!                     note that each atom (not species) contributes to this number
!                     separately, so that in most cases n_basis_fns < n_compute_fns < n_basis
!    o i_basis_fns --  radial basis function indeces (up to n_basis_fns, so only per species)
!    o i_atom_fns -- gives the global atom number associated with the radial basis function
!    o i_basis_fns_inv -- inverse index list of i_basis_fns
!    o spline_array_start -- starting point of splined data information
!    o spline_array_end --  ending point of splined data information
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

	






!     local variables
  
  integer :: i_offset_spl
  logical :: spline_start_found
  
!     counters

  integer :: i_basis_1
  integer :: i_compute
  integer :: i_compute_1
  
  integer :: i_atom
  integer :: i_atom_compute
  integer :: i_index_2
  integer :: current_atom

  integer :: i_spline

  !     begin work

  !     check for atoms which contribute basis functions at the current integration point at all
  do i_atom = 1, n_atoms, 1

     ! if this atom is relevant at this point
     if ( dist_tab_sq(i_atom) .le. &
          atom_radius_sq(species(i_atom)) .and. &
          (dist_tab_sq(i_atom).gt.0.d0) ) then

        n_compute_atoms = n_compute_atoms + 1
        atom_index(n_compute_atoms) = i_atom
        atom_index_inv(i_atom) = n_compute_atoms

     end if

     !     loop over atoms
  enddo

  !     next, check for radial basis functions
  do i_atom_compute = 1, n_compute_atoms, 1

     i_atom = atom_index(i_atom_compute)

     ! This is the index where the splines for the present atom
     ! begin within the condensed array of all splined radial functions:
     i_offset_spl = basis_fn_start_spl(species(i_atom))
     spline_start_found = .false.

     ! simply run through array of splined functions to determine
     ! the beginning of the section of the spline array which is
     ! needed from the present atom at the present integration point
     do i_spline = 1, n_basis_fns, 1

        ! original radial function index associated with current spline section
        i_basis_1 = perm_basis_fns_spl(i_spline)

        !          if current basis function is associated with current atom
        if (basis_fn_atom(i_basis_1, i_atom)) then

           if ( dist_tab_sq(i_atom).le. &
                outer_radius_sq(i_basis_1)) then

              spline_array_start(i_atom_compute) = &
                   (i_spline-1) + 1

              spline_start_found = .true.
              exit

           end if

        end if
     enddo

     ! Now, determine separately the correspondence between radial functions
     ! used at the current point, and radial functions as they were originally stored
     ! in basis_fn(i_basis)
     do i_basis_1 = 1, n_basis_fns, 1
        !          if current radial function is associated with current atom
        if (basis_fn_atom(i_basis_1, i_atom)) then

           !            if this function is relevant at this point
           if ( dist_tab_sq(i_atom) .le. &
                outer_radius_sq(i_basis_1) ) then

              n_compute_fns = n_compute_fns + 1
              i_basis_fns(n_compute_fns) = i_basis_1
              i_atom_fns(n_compute_fns) = i_atom
              i_basis_fns_inv(i_basis_1,i_atom) = n_compute_fns

           end if

           i_offset_spl = i_offset_spl + 1

        end if

        !        loop over basis functions
     enddo

     ! The end of the splined array is simply the end of the block where
     ! the splined functions of the current species are stored, which is
     ! what we have after the loop over basis functions ...
     spline_array_end(i_atom_compute) = i_offset_spl - 1

     !     loop over atoms
  enddo


end subroutine prune_radial_basis_v2
!---------------------------------------------------------------------
!******
