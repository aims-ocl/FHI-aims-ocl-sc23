!****s* FHI-aims/prune_radial_basis_p0
!  NAME
!   prune_radial_basis_p0
!  SYNOPSIS

subroutine prune_radial_basis_p0 &
     ( dist_tab_sq, dist_tab, &
     dir_tab, &
     n_compute_atoms, atom_index, atom_index_inv, &
     n_compute_fns, i_basis_fns, i_basis_fns_inv, &
     i_atom_fns, spline_array_start, spline_array_end, &
     n_atom_list, atom_list)
  
!  PURPOSE
!     reduces full set radian basis functions and atom centers to
!     relevant ones only, for a given set of n_points integration points
!
!  USES

  use dimensions
  use basis
  use grids
  use geometry
  use pbc_lists
  implicit none

!  ARGUMENTS

  integer:: n_atom_list
  integer:: atom_list(n_atom_list)  
  real*8 dist_tab_sq(n_atom_list)
  real*8 dist_tab(n_atom_list)
  real*8, dimension(3, n_atom_list) :: dir_tab

  integer :: n_compute_atoms
  integer :: atom_index(n_atom_list)
  integer :: atom_index_inv(n_centers)

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_atom_list)
  integer :: i_atom_fns(n_basis_fns*n_atom_list)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: spline_array_start(n_atom_list)
  integer :: spline_array_end(n_atom_list)
      

!  INPUTS
!    o n_atom_list -- total number of atom centers, including periodic mirror images
!    o atom_list -- list of the total number of atom centers
!    o dist_tab_sq -- (distance to atom centers)**2
!    o dir_tab -- direction to atom centers
!
!  OUTPUT
!    o dir_tab -- direction to atom centers (normalized)
!    o dist_tab -- distance to atoms (normalized)
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
  integer :: i_center, i_center_L
  integer :: i_atom_compute


  integer :: i_spline

  !     begin work

  !     check for atoms which contribute basis functions at the current integration point at all

  do i_center_L = 1, n_atom_list, 1
     i_center = atom_list(i_center_L)

     ! if this atom is relevant at this point
     if ( dist_tab_sq(i_center_L) .le. &
          atom_radius_sq(species_center(i_center)) .and. &
          (dist_tab_sq(i_center_L).gt.0.d0) ) then

        n_compute_atoms = n_compute_atoms + 1
        atom_index(n_compute_atoms) = i_center
        atom_index_inv(i_center) = n_compute_atoms


        dist_tab_sq(n_compute_atoms) = dist_tab_sq(i_center_L)
        dist_tab(n_compute_atoms) = sqrt(dist_tab_sq(n_compute_atoms))


        dir_tab(1:3,n_compute_atoms) = dir_tab(1:3,i_center_L)/dist_tab(n_compute_atoms) 


     end if
     !     loop over atoms
  enddo

  !     next, check for radial basis functions
  do i_atom_compute = 1, n_compute_atoms, 1
     i_center = atom_index(i_atom_compute)

!  DB: if we hit one of the pseudoatoms, we have to CYCLE, since the pseudoatom don't have 
!      this kind of basiswaves - the pseudowaves are treated differently (see pseudodata.f90)
!      if (use_embedding_pp.and.(center_to_atom(i_center).gt.n_real_atoms)) cycle

     ! This is the index where the splines for the present atom
     ! begin within the condensed array of all splined radial functions:
     i_offset_spl = basis_fn_start_spl(species_center(i_center))
     spline_start_found = .false.

     ! simply run through array of splined functions to determine 
     ! the beginning of the section of the spline array which is
     ! needed from the present atom at the present integration point
     do i_spline = 1, n_basis_fns, 1

        ! original radial function index associated with current spline section
        i_basis_1 = perm_basis_fns_spl(i_spline)

        !          if current basis function is associated with current atom
        if (basis_fn_atom(i_basis_1, center_to_atom(i_center))) then

           if ( dist_tab_sq(i_atom_compute).le. &
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
        if (basis_fn_atom(i_basis_1, center_to_atom(i_center))) then


           !            if this function is relevant at this point
           if ( dist_tab_sq(i_atom_compute) .le. &
                outer_radius_sq(i_basis_1) ) then


              n_compute_fns = n_compute_fns + 1


              i_basis_fns(n_compute_fns) = i_basis_1
              i_atom_fns(n_compute_fns) = i_center
              i_basis_fns_inv(i_basis_1,i_center) = n_compute_fns

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



end subroutine prune_radial_basis_p0
!---------------------------------------------------------------------
!******
