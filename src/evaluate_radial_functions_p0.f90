!****s* FHI-aims/evaluate_radial_functions_p0
!  NAME
!   evaluate_radial_functions_p0
!  SYNOPSIS


subroutine evaluate_radial_functions_p0( &
     spline_array_start, spline_array_end, &
     n_compute_atoms, n_compute_fns, dist_tab, & 
     i_r, atom_index, i_basis_fns_inv,  &
     spline_data, wave_aux, derivative, n_compute, n_basis_list &
     )

!  PURPOSE
!     This subroutine evaluates wave arrays radial_wave, kinetic_wave, and
!     radial_wave_deriv (if needed) form splines using vectorisation
!     over the spline index.
!
!     The data organisation follows from shrink_fixed_basis.
!
!     The spline data is organised in the global array basis_all_spl in the
!     following order:
!     o The species constitute the main blocks of this array
!     o Within each species the spline are arranged into blocks of increasing order with
!        respect to their cutoff radiuses.
!     o Within each spline block there are two or three components to be used without or with
!        density gradients, respectively. They are ordered as: radial wave, kinetic wave, 
!        (derivative of radial wave)
!
!  USES

  use dimensions
  use grids
  use pbc_lists
  use spline
  use basis, only : perm_basis_fns_spl 

  implicit none

!  ARGUMENTS

  integer :: n_compute
  integer :: n_compute_atoms
  integer :: n_compute_fns
  integer:: n_basis_list

  integer :: spline_array_start(n_compute_atoms)
  integer :: spline_array_end(n_compute_atoms)

  real*8 :: dist_tab(n_compute_atoms)
  real*8 :: i_r(n_compute_atoms)

  integer :: atom_index(n_compute_atoms)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  real*8 :: spline_data(n_basis_fns,n_max_spline, n_max_grid)
  logical :: derivative
  real*8 :: wave_aux(n_basis_list)

! INPUTS
! o n_compute -- not use anymore
! o n_compute_atoms -- number of relevant atoms
! o n_compute_fns -- not used anymore
! o n_basis_list
! o spline_array_start -- starting point of radial functions spline arrow
! o spline_array_end -- ending point of radial functions spline arrow
! o dist_tab  -- distance to atoms
! o i_r(n_compute_atoms)
! o i_r -- the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
! o atom_index -- list of relevant atoms
! o i_basis_fns_inv -- inverse citing list of radial functions
! o spline_data -- spline data for radial functions
! o derivative -- ?????????????
!
! OUTPUT
! o wave_aux -- radial part of basis functions <- MR: it can also be the radial wave derivative if forces are on.
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





  !     locals

  real*8 :: r_point

  integer :: n_spline
  integer :: spline_start
  integer :: spline_end
  real*8 :: spline_array_aux(n_basis_fns)

  integer :: current_atom
  integer :: current_species
  integer :: current_basis_fn
  integer :: current_basis_fn_comp

  real*8 :: distance_from_atom

  !     counters

  integer :: i_atom_1
  integer :: i_spline
  integer :: i_rad
  integer :: spline_counter

  !      write(use_unit,*) 'r', n_compute
  !  write(use_unit,*) i_basis_fns_inv

  do i_atom_1 = 1, n_compute_atoms, 1

     current_atom = atom_index(i_atom_1)
     current_species = species_center(current_atom)

     spline_start = spline_array_start(i_atom_1)
     spline_end = spline_array_end(i_atom_1)

     !         print *, current_atom, spline_start, spline_end


     n_spline = spline_end - spline_start + 1

     r_point = i_r(i_atom_1)
     distance_from_atom = dist_tab(i_atom_1)

     call spline_vector_waves( r_point, spline_data, &
          n_max_grid, n_basis_fns, spline_start, &
          spline_end, n_grid(current_species), &
          n_spline, spline_array_aux )

     i_rad = (spline_start-1)

     ! now sort back splined radial functions into original
     ! array of computed radial functions at present point
     do i_spline = spline_start, spline_end, 1

        i_rad = i_rad + 1
        spline_counter = i_spline - (spline_start-1)

        current_basis_fn = perm_basis_fns_spl(i_rad)
        current_basis_fn_comp = i_basis_fns_inv(current_basis_fn,current_atom)
        if(current_basis_fn_comp.eq.0)cycle ! (Rundong) if we don't use basis permutation, this line is a must.

        if (derivative) then
           if (distance_from_atom.gt.0) then
              wave_aux(current_basis_fn_comp) = &
                   spline_array_aux(spline_counter)
           else
              wave_aux(current_basis_fn_comp) = 0.0d0
           end if
        else
           wave_aux(current_basis_fn_comp) = &
                spline_array_aux(spline_counter)
        end if
     enddo
  enddo


end subroutine evaluate_radial_functions_p0
!****** 
