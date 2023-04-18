!****s* FHI-aims/evaluate_radial_functions_v2
!  NAME
!   evaluate_radial_functions_v2
!  SYNOPSIS

subroutine evaluate_radial_functions_v2( &
     spline_array_start, spline_array_end, &
     n_compute_atoms, n_compute_fns, dist_tab, &
     i_r, atom_index, i_basis_fns_inv, &
     spline_data, wave_aux, derivative &
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
!     * The species constitute the main blocks of this array
!     * Within each species the spline are arranged into blocks of increasing order with
!        respect to their cutoff radiuses.
!     * Within each spline block there are two or three components to be used without or with
!        density gradients, respectively. They are ordered as: radial wave, kinetic wave,
!        (derivative of radial wave)
!
!  USES

  use dimensions
  use grids
  use geometry
  use spline
  use basis, only : perm_basis_fns_spl
  implicit none

!  ARGUMENTS

  integer :: n_compute_atoms
  integer :: n_compute_fns

  integer :: spline_array_start(n_atoms)
  integer :: spline_array_end(n_atoms)

  real*8 :: dist_tab(n_atoms)
  real*8 :: i_r(n_atoms)

  integer :: atom_index(n_atoms)
  integer :: i_basis_fns_inv(n_basis_fns,n_atoms)
  real*8 :: spline_data(n_max_spline,n_max_grid,n_basis_fns)

  logical :: derivative
  real*8 :: wave_aux(n_basis)


! INPUTS
! o n_compute_atoms -- number of relevant atoms
! o n_compute_fns -- number of relevant radial functions
! o spline_array_start -- starting point of radial functions spline arrow
! o spline_array_end -- ending point of radial functions spline arrow
! o dist_tab -- distance to atoms
! o i_r -- the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
! o atom_index -- list of relevant atoms
! o i_basis_fns_inv -- inverse citing list of radial functions
! o spline_data -- spline data for radial functions
! o derivative -- ?????????????
!
! OUTPUT
! o wave_aux -- radial part of basis functions.
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
  integer :: i_point
  integer :: i_spline
  integer :: i_rad
  integer :: spline_counter

  do i_atom_1 = 1, n_compute_atoms, 1

     current_atom = atom_index(i_atom_1)
     current_species = species(current_atom)

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
        current_basis_fn_comp = &
             i_basis_fns_inv(current_basis_fn,current_atom)

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

end subroutine evaluate_radial_functions_v2
!******	
