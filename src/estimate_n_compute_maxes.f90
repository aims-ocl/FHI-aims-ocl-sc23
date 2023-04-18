!****s* FHI-aims/estimate_n_compute_maxes
!  NAME
!    estimate_n_compute_maxes
!  SYNOPSIS 
subroutine estimate_n_compute_maxes( )
!  PURPOSE
!    Estimates the threadwise maximal number of non-zero basis functions based
!    on a geometrical analysis of the structure and the basis functions.
!  USES
  use dimensions
  use grids
  use runtime_choices
  use localorb_io
  use pbc_lists
  use mpi_utilities
  use geometry, only: species

  implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    n_max_compute_ham and n_max_compute_fns_ham are set on exit
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
! SOURCE

  
  ! locals

  integer :: i_atom, i_coord, i_angular, i_grid, i_division
  integer :: i_radial
  integer :: division_low, division_high

  real*8 :: coord_current(3), coord_2(3)
  real*8 :: dist_tab_sq( n_centers_basis_integrals, n_max_angular )
  real*8 :: dist_tab( n_centers_basis_integrals )
  real*8 :: dir_tab( 3, n_centers_basis_integrals, n_max_angular )

  integer :: n_compute, n_compute_2
  integer :: i_basis(n_centers_basis_I)

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_centers_basis_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_basis_integrals)

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_basis_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_centers_basis_integrals)
  integer :: spline_array_end(n_centers_basis_integrals)

  character*100 :: info_str



  n_max_compute_ham = 0
  n_max_compute_fns_ham = 0

  do i_atom = 1, n_atoms,1

     do i_radial = 1, n_radial(species(i_atom)), 1

        if (myid.eq.radial_task_list(i_radial,i_atom)) then

           do i_grid = 1, n_grids, 1

              do i_division = 1, n_division_lebedev(i_grid), 1

                 division_low = division_boundaries_lebedev(i_division,i_grid) + 1
                 division_high = division_boundaries_lebedev(i_division+1,i_grid)

                 n_compute = 0
                 n_compute_2 = 0
                 i_basis = 0

                 do i_angular = division_low, division_high, 1

                    coord_current(:) = coords_center(:,i_atom) + &
                    r_angular_lebedev(:,i_angular,i_grid) * &
                    r_radial(i_radial, species(i_atom))

                    if (n_periodic > 0) then
                       call map_to_center_cell( coord_current )
                    end if

                    ! compute atom-centered coordinates of current integration point,
                    ! as viewed from all atoms
                    call tab_atom_centered_coords_p0 &
                    ( coord_current,  &
                    dist_tab_sq(1, i_angular),  &
                    dir_tab(1,1,i_angular), &
                    n_centers_basis_integrals, &
                    centers_basis_integrals )

                    ! determine which basis functions are relevant at current integration point,
                    ! and tabulate their indices

                    ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
                    call prune_basis_p0 &
                    ( dist_tab_sq(1,i_angular), &
                    n_compute, n_compute_2, i_basis,  &
                    n_centers_basis_I, n_centers_integrals, &
                    inv_centers_basis_integrals  )

                 end do

                 n_max_compute_ham = MAX(n_compute,n_compute_2,n_max_compute_ham)

                 if (n_compute.gt.0) then

                    do i_angular = division_low, division_high, 1

                       n_compute_atoms = 0
                       n_compute_fns = 0

                       call prune_radial_basis_p0 &
                       ( dist_tab_sq(1,i_angular), &
                       dist_tab, &
                       dir_tab(1,1,i_angular), &
                       n_compute_atoms, atom_index, atom_index_inv, &
                       n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                       i_atom_fns, spline_array_start, spline_array_end, &
                       n_centers_basis_integrals, centers_basis_integrals)

                       n_max_compute_fns_ham = MAX(n_compute_fns,n_max_compute_fns_ham)

                    end do

                 end if

              end do

           end do
        end if
     end do

  end do
  !     n_max_compute_fns_ham = n_max_compute_ham

  write(info_str,'(2X,A,I8,A,I3)') "| Estimated number of non-zero basis functions for the Hamiltonian : ", &
  n_max_compute_ham, " in task ", myid
  call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)
  write(info_str,'(2X,A,I8,A,I3)') "| Estimated number of non-zero radial functions for the Hamiltonian: ", &
  n_max_compute_fns_ham, " in task ", myid
  call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)

end subroutine estimate_n_compute_maxes
!******
