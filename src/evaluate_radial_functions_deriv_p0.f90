!****s* FHI-aims/evaluate_radial_functions_deriv_p0
!  NAME
!    evaluate_radial_functions_deriv_p0
!  SYNOPSIS

subroutine evaluate_radial_functions_deriv_p0( &
     spline_array_start, spline_array_end, &
     n_compute_atoms, n_compute_fns, dist_tab, &
     i_r, atom_index, i_basis_fns_inv, &
     spline_data, wave_aux, derivative, n_size &
     )

!  PURPOSE
!     This subroutine evaluates wave arrays radial_wave_deriv form splines using vectorisation
!     over the spline index.
!
!     The data organisation follows from shrink_fixed_basis.
!
!     The spline data is organised in the global array basis_all_spl in the
!     following order:
!     o The species constitute the main blocks of this array
!     o Within each species the spline are arranged into blocks of increasing order with
!       respect to their cutoff radiuses.
!     o Within each spline block there are two or three components to be used with density gradients, respectively. 
!
!  USES

  use dimensions
  use grids
  use geometry
  use spline
  use basis, only : perm_basis_fns_spl
  use pbc_lists
  implicit none

!  ARGUMENTS


  integer :: n_compute_atoms
  integer :: n_compute_fns
  integer :: spline_array_start(n_compute_atoms)
  integer :: spline_array_end(n_compute_atoms)
  real*8 :: dist_tab(n_compute_atoms)
  real*8 :: i_r(n_compute_atoms)
  integer :: atom_index(n_compute_atoms)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  real*8 :: spline_data(n_max_spline,n_max_grid,n_basis_fns)
  logical :: derivative
  integer :: n_size
  real*8 :: wave_aux(n_size)


! INPUTS
! o spline_array_start -- starting point of spline array
! o spline_array_end -- ending point of spline array
! o n_compute_atoms -- number of relevant atoms 
! o n_compute_fns -- not used anymore
! o dist_tab -- distance to atoms
! o i_r -- the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
! o atom_index -- list of relevant atoms
! o i_basis_fns_inv -- inverse citing list of radial basis.
! o spline_data -- radial functions spline data
! o derivative -- not used anymore
! o n_size -- dimension of wave aux
! 
! OUTPUT
! o wave_aux -- derivative of radial basis functions
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
  real*8 :: scale

  !     counters

  integer :: i_atom_1
  integer :: i_spline
  integer :: i_rad
  integer :: spline_counter
  do i_atom_1 = 1, n_compute_atoms, 1

     current_atom = atom_index(i_atom_1)
     current_species = species_center(current_atom)
     distance_from_atom = dist_tab(i_atom_1)

     scale = 1.d0 / (log(r_grid_inc(current_species)) &
          * distance_from_atom)

     spline_start = spline_array_start(i_atom_1)
     spline_end = spline_array_end(i_atom_1)

     !         print *, current_atom, spline_start, spline_end

     n_spline = spline_end - spline_start + 1

     r_point = i_r(i_atom_1)

     call spline_vector_waves_deriv( r_point, spline_data, &
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
        current_basis_fn_comp = &
             i_basis_fns_inv(current_basis_fn,current_atom)

        if (distance_from_atom.gt.0) then
           ! scale derivative to get derivative w.r.t. to space
           ! instead of spline parameter (* di/dr)
           wave_aux(current_basis_fn_comp) = &
                spline_array_aux(spline_counter) * scale
        else
           wave_aux(current_basis_fn_comp) = 0.0d0
        end if
     enddo
  enddo

end subroutine evaluate_radial_functions_deriv_p0
!******	
