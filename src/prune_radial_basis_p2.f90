!****s* FHI-aims/prune_radial_basis_p2
!  NAME
!   prune_radial_basis_p2
!  SYNOPSIS

subroutine prune_radial_basis_p2 &
     ( dim_atoms, dim_fns, &
     dist_tab_sq, dist_tab, dir_tab, &
     n_compute_atoms, atom_index, atom_index_inv, &
     n_compute_fns, i_basis_fns, i_basis_fns_inv, &
     i_atom_fns, spline_array_start, spline_array_end, &
     n_atom_list, atom_list, n_compute, i_basis, &
     n_batch_centers, batch_center, &
     one_over_dist_tab, rad_index, wave_index, l_index, l_count, &
     fn_atom, n_zero_compute, zero_index_point &
     )

!  PURPOSE
!
!     Subroutine prune_radial_basis_p2
!
!     for a given integration point, determines and indexes:
!     *  those atoms from which at least one (radial) wave function is needed at this point 
!     *  those radial functions which are non-zero at this point (# n_compute_fns)
!     *  among _all_ wave functions relevant in current batch of points (# n_compute), those 
!         which are zero at current point. (# n_zero_compute)
!
!  USES

  use dimensions
  use basis
  use grids
  use geometry
  use pbc_lists
  use localorb_io, only: use_unit
  implicit none
    
!  ARGUMENTS

  integer :: dim_atoms
  integer :: dim_fns
  integer:: n_atom_list
  integer:: atom_list(n_atom_list)

  real*8 :: dist_tab_sq(n_atom_list)
  real*8 :: dist_tab(n_atom_list)
  real*8, dimension(3, n_atom_list) :: dir_tab

  integer :: n_batch_centers
  integer :: batch_center(n_atom_list)

  integer :: n_compute_atoms
  integer :: atom_index(n_atom_list)
  integer :: atom_index_inv(n_centers)

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_atom_list)
  integer :: i_atom_fns(n_basis_fns*n_atom_list)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)

  integer :: spline_array_start(n_atom_list)
  integer :: spline_array_end(n_atom_list)

  integer :: n_compute
  integer :: i_basis(n_compute)

  real*8  :: one_over_dist_tab(n_atom_list)

  integer :: rad_index(dim_atoms)
  integer :: wave_index(dim_fns)
  integer :: l_index(dim_fns)
  integer :: l_count(dim_fns)
  integer :: fn_atom(dim_fns)
  
  integer :: n_zero_compute
  integer :: zero_index_point(n_compute)


!  INPUTS
! 
!   o dim_atoms -- dimensions that bound n_compute_atoms
!   o dim_fns -- dimensions that bound n_compute_fns
!   o n_atom_list -- total number of atom centers
!   o atom_list -- list of atom centers
!   o dist_tab_sq -- (distance to atom centers)**2
!   o dir_tab -- direction to atoms
!   o n_batch_centers -- number of batch centers
!   o batch_center -- list of batch centers
!
!  OUTPUT

!   o dist_tab -- distance to relevant atoms
!   o dir_tab -- direction to relevant atoms (normalized)
!   o n_compute_atoms -- number of relevant atoms
!   o atom_index --  list of relevant atoms
!   o atom_index_inv -- inverse citing list of atom_index
!   o n_compute_fns --  total number of radial functions relevant at these points
!                     note that each atom (not species) contributes to this number
!                     separately, so that in most cases n_basis_fns < n_compute_fns < n_basis
!   o i_basis_fns -- list of relevant radial basis functions
!   o i_basis_fns_inv --  inverse citing list of i_basis_fns
!   o spline_array_start -- starting point of splined data information
!   o spline_array_end --  ending point of splined data information
!   o n_compute -- number of relevant basis functions ( For indexing of full basis, to avoid any reindexing later)
!   o i_basis -- list of relevant basis functions ( For indexing of full basis, to avoid any reindexing later)
!   o one_over_dist_tab -- 1/(distance to atoms)
!   o  rad_index -- in the list of n_compute_fns radial functions, 
!                     the END of the radial function list associated with atom i_compute_atom
!   o wave_index -- in the list of non-zero wave functions, the START of
!                     those wave functions associated with current radial function i_compute_fn
!   o l_index -- in the list of ylm functions, the START of
!                     those ylm functions associated with current radial function i_compute_fn
!   o l_count -- in the list of ylm functions, the NUMBER of
!                     ylm functions associated with current radial function i_compute_fn, MINUS 1
!   o i_atom_fns -- in the list of non-zero atoms at the current point, the index
!                     of the atom at which radial function i_compute_fn is centered
!   o fn_atom -- ??????
!   o n_zero_compute -- number of  known zero basis functions at current point
!   o zero_index_point -- list of  known zero basis functions at current point
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
  integer :: l_aux

  !     counters

  integer :: i_batch_center
  integer :: i_basis_1
  integer :: i_compute
  integer :: i_fn
  integer :: i_lm
  integer :: i_center, i_center_L
  integer :: i_atom_compute
  integer :: i_spline




  !test
  !      logical :: t_out

  !     begin work

  ! Memo for prune_basis / prune_basis_radial:
  !
  ! * prune_basis needs to evaluate (and store!!) all relevant atoms in the current batch, so that prune_radial only
  ! needs to rescan those atoms, no others; not sure if this is already done, don't think so.


  !     check for atoms which contribute basis functions at the current integration point at all
  !     VB: This is now reduced to the list of only the atoms that are active in current batch of
  !         integration points. Thus goes the last O(N^2) part of the integration.
  do i_batch_center = 1, n_batch_centers, 1

     i_center_L = batch_center(i_batch_center)

!test
     if (.not.(dist_tab_sq(i_center_L).gt.0.d0)) then
       write(use_unit,*) "zero: ", i_center_L, dist_tab_sq(i_center_L)
     end if
!test end

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

        one_over_dist_tab(n_compute_atoms) = 1.d0/dist_tab(n_compute_atoms) 

        dir_tab(1:3,n_compute_atoms) = & 
             dir_tab(1:3,i_center_L) * one_over_dist_tab(n_compute_atoms) 
     else

        ! Some routines, like small component needs this one.
        atom_index_inv(i_center) = 0

     end if

     !     loop over atoms
  enddo


  !     next, check for radial basis functions
  do i_atom_compute = 1, n_compute_atoms, 1

     i_center = atom_index(i_atom_compute)

     ! This is the index where the splines for the present atom
     ! begin within the condensed array of all splined radial functions:
     i_offset_spl = basis_fn_start_spl(species_center(i_center))

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
              ! Found the start of the non-zero radial spline functions for current atom;
              ! store for future use.
              spline_array_start(i_atom_compute) = &
                   (i_spline-1) + 1

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

              fn_atom(n_compute_fns) = i_atom_compute

           end if

           i_offset_spl = i_offset_spl + 1

        end if

        !        loop over radial functions
     enddo
     ! store the last fn associated with current atom in functions list above
     rad_index(i_atom_compute) = n_compute_fns

     ! The end of the splined array is simply the end of the block where
     ! the splined functions of the current species are stored, which is
     ! what we have after the loop over basis functions ...
     spline_array_end(i_atom_compute) = i_offset_spl - 1

     !     loop over atoms
  enddo

  ! Now prepare reverse indexing for the benefit of pointwise wave functions later.
  ! This avoids a lot of cumbersome if statements.

  n_zero_compute = 0

  ! Loop over all basis functions u(r)/r * Ylm(Omega) that are active in current 
  ! BATCH of integration points. Sort out those basis functions that are ZERO at
  ! current point, and index those that are non-zero for later use.
  ! NOTE that this sequence relies on the fact that all Ylm functions associated
  ! with a single radial function appear DIRECTLY after one another in the array of
  ! active basis functions.

  wave_index(1:n_compute_fns) = 0

  !      if (t_out) then
  !        write(use_unit,*) "n_compute: ", n_compute
  !      end if

  i_compute = 1

  do while (i_compute .le. n_compute)
     i_basis_1 = i_basis(i_compute)
     ! radial function in list of i_compute_fns that corresponds to 
     ! current basis function
     i_fn = i_basis_fns_inv(basis_fn(Cbasis_to_basis(i_basis_1)), Cbasis_to_center(i_basis_1))
     l_aux = basis_l(Cbasis_to_basis(i_basis(i_compute)))
     if ( i_fn .eq. 0 ) then
        ! that basis function was listed as zero before and does not exist in
        ! list of radial functions

        ! the next ( 2*basis_l(i_basis(i_compute))+1 ) elements of wave
        ! will always be zero

        do i_lm = 0, 2*l_aux, 1
           n_zero_compute = n_zero_compute + 1
           zero_index_point(n_zero_compute) = i_compute + i_lm
        enddo

     else

        if (wave_index(i_fn).eq.0) then
           ! first time that this function shows up in the list - collect all
           ! related indices here.
           wave_index(i_fn) = i_compute

           ! starting index for current l in ylm_tab
           l_index(i_fn) = l_aux**2 + 1

           ! number of Ylm fns (minus one) for current l
           l_count(i_fn) = 2*l_aux
        end if

     end if
     i_compute = i_compute + 2*l_aux + 1
  enddo
!  This is a very useful debug sequence for strange crashes in the integrations.
!  If a radial function that _should_ be nonzero has no wave function pendant that
!  is non-zero, something has gone wrong, and the following check catches it.
!
!  do i_fn = 1, n_compute_fns
!     if( wave_index(i_fn) < 1)then
!        write(use_unit,*)  wave_index(1:n_compute_fns)
!        write(use_unit,*) i_fn, wave_index(i_fn)
!        stop
!     end if
!  end do
end subroutine prune_radial_basis_p2
!---------------------------------------------------------------------
!******
