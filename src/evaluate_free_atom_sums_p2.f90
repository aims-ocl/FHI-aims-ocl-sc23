!****s* FHI-aims/evaluate_free_atom_sums_p2
!  NAME
!   evaluate_free_atom_sums_p2
!  SYNOPSIS

subroutine evaluate_free_atom_sums_p2 &
     ( i_r, dir_tab, &
     free_hartree_superpos, free_rho_superpos, &
     free_rho_gradient_superpos, &
     rho, rho_gradient, n_atom_list, atom_list, & 
     n_compute_atoms, i_compute2i_atom, temp_free_rho, &
     i_full_points )
     !SR: to get laplacian of free atom density add this args: (everythere where evaluate_free_atom_sums_p2 appears... )
     !, &
     !dist_tab_sq,free_rho_laplace_superpos)
     !the laplacian like implemented here works (tested for small molecules)

!  PURPOSE
!  Subroutine evaluate_free_atom_sums tabulates the superposition of free-atom
!  densities and hartree potentials on the entire integration grid, for use 
!  in the construction of the Hartree potential.
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use pbc_lists
  use species_data
  use spline
  use free_atoms
  use constants
  use heat_flux
  implicit none

!  ARGUMENTS

  integer :: n_atom_list
  integer :: atom_list(n_atom_list)

  integer:: n_compute_atoms
  integer:: i_compute2i_atom(n_compute_atoms)
  real*8 :: temp_free_rho(n_compute_atoms)

  real*8, dimension(n_atom_list) :: i_r
  real*8, dimension(3,n_atom_list) :: dir_tab
!  real*8, dimension(n_atom_list) :: dist_tab_sq
  
  real*8 :: free_hartree_superpos 
  real*8 :: free_rho_superpos 
  real*8, dimension(3) :: free_rho_gradient_superpos 
  real*8, dimension(n_spin) :: rho
  real*8, dimension(3,n_spin) :: rho_gradient
  integer,optional :: i_full_points
!SR for laplacian
!  real*8 :: free_rho_laplace_superpos
!  real*8 :: atomic_laplacian  
!  real*8 :: log_weight, radial_weight

!  INPUTS
!    o n_atom_list -- number of atoms
!    o atom_list -- list of atoms
!    o n_compute_atoms -- number of relevant atoms
!    o i_compute2i_atom -- this list contains the indices of relevant atom 
!                           centers in atom_list
!    o temp_free_rho -- this data is added to free_rho_superpos
!    o i_r  the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!    o dir_tab -- direction to atoms
!   
!  OUTPUT
!    o free_hartree_superpos -- superposition of free atoms Hartree potential
!    o free_rho_superpos --  superposition of free atoms electronic charge times 4*pi
!    o free_rho_gradient_superpos  -- superposition of free atoms gradient of electronic charge
!    o rho -- electronic charge is initialized to be free_rho_superpos/( 4*pi)
!    o rho_gradient -- gradient of electronic charge is initialized to be free_rho_gradient_superpos
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

    

  !  local variables

  real*8 atomic_gradient
  !      real*8 spin_sign

  !     counters

  integer :: i_center, i_center_L
  integer :: i_coord
  integer :: i_spin
  integer :: i_atom
  integer :: i_atom_2

  integer :: i_compute_atom
!  real :: i_r_rad
  !  begin work


  !       initialize
  free_rho_superpos = 0.d0
  free_hartree_superpos = 0.d0
  if (use_density_gradient) then
     free_rho_gradient_superpos = 0.d0
!SR laplacian
!      free_rho_laplace_superpos = 0d0
  end if
  rho = 0.d0

  !       Contributions from all atoms
  ! By way of the lists above, we only count occupied sites here, no ghost atoms
  do i_compute_atom = 1, n_compute_atoms, 1
! We need to exlude free rho of pseudoized species
    if (species_pseudoized(species_center(atom_list(i_compute_atom)))) cycle

!DB 01. october 13: VB, actually initialize_grid_storage excludes 
!   all points associated to empty sites to be considered here

     i_atom = i_compute2i_atom(i_compute_atom)
     i_center = atom_list(i_atom)

     ! electrostatic contribution
     free_hartree_superpos =  &
          free_hartree_superpos +  &
          val_spline &
          ( i_r(i_atom), free_pot_es_spl(1,1,species_center(i_center)), &
          n_grid(species_center(i_center)) )

     !            density contribution
     free_rho_superpos = free_rho_superpos +  &
          temp_free_rho(i_compute_atom)
     ! CC: Store per atom contribution to free atom rho
     ! at each grid point
     if (compute_heat_flux) then
       HF_rho_free_per_atom(i_full_points,Center_to_atom(atom_list(i_compute2i_atom(i_compute_atom)))) = &
         HF_rho_free_per_atom(i_full_points,Center_to_atom(atom_list(i_compute2i_atom(i_compute_atom)))) +  temp_free_rho(i_compute_atom)
     end if


     if (use_density_gradient) then
        ! density gradient contribution
        ! factor pi4 already divided out!! different to free_rho_superpos!!!

        atomic_gradient = &
             val_spline( i_r(i_atom), &
             free_drho_dr_spl(1,1,species_center(i_center)), &
             n_grid(species_center(i_center)) )

        atomic_gradient = atomic_gradient * pi4_inv

        do i_coord = 1,3,1
           free_rho_gradient_superpos(i_coord) = &
                free_rho_gradient_superpos(i_coord) + &
                atomic_gradient * dir_tab(i_coord, i_atom)
        enddo
!SR: for laplacian of electron density. commented everything out, since I do not need
!it anymore, but maybe it is useful for someone
! 	i_r_rad =& !i_r(i_atom)
! 	  invert_radial_grid &
!        ( sqrt(dist_tab_sq(i_atom)), &
!        n_radial(species_center(i_center)), &
!        scale_radial(species_center(i_center)))

!         call tab_single_radial_weights_v2 &
!            ( i_atom, sqrt(dist_tab_sq(i_atom)), i_r_rad, &
!             log_weight, radial_weight )        
!          atomic_laplacian = &
!               val_spline_deriv( i_r(i_atom), &
!               free_drho_dr_spl(1,1,species_center(i_center)), &
!               n_grid(species_center(i_center)) )* log_weight           
!         atomic_laplacian = atomic_laplacian * pi4_inv
!         free_rho_laplace_superpos = &
! 	    free_rho_laplace_superpos + &
! 	    atomic_laplacian + 2d0*atomic_gradient/sqrt(dist_tab_sq(i_atom))
     end if

  enddo

  if (.not.use_initial_rho) then
     ! non-polarized case
     ! (Rundong) For fully-relativistic cases, we currently use this branch.

     i_spin = 1
     rho(i_spin) = pi4_inv * free_rho_superpos

  else

     do i_spin = 1, n_spin, 1

        do i_center_L = 1, n_atom_list, 1

           i_center = atom_list(i_center_L)

           i_atom_2 = center_to_atom(i_center)
           if (.not.species_pseudoized(species_center(atom_list(i_center_L)))) then


              if (.not.(empty(i_atom_2))) then
                 rho(i_spin) = rho(i_spin) + val_spline( i_r(i_center_L), &
                      initial_rho_spl(1,1,atom_type(i_atom_2),i_spin), &
                      n_grid(species(i_atom_2)) )     
              endif

           endif 

        enddo

     enddo

  end if

  !       if required, obtain density gradient also
  if (use_density_gradient) then
     rho_gradient = 0.d0

     if (.not.use_initial_rho) then

        i_spin = 1
        rho_gradient(:, i_spin) = free_rho_gradient_superpos(:)

     else

        !           Contributions from all atoms
        do i_center_L = 1, n_atom_list, 1

           i_center = atom_list(i_center_L)

           i_atom_2 = center_to_atom(i_center)

           if (.not.empty(i_atom_2)) then
              do i_spin = 1, n_spin, 1

                 atomic_gradient = val_spline( i_r(i_center_L), &
                      initial_drho_dr_spl(1,1,atom_type(i_atom_2), &
                      i_spin), n_grid(species(i_atom_2)) )

                 do i_coord = 1,3,1

                    rho_gradient(i_coord,i_spin) = &
                         rho_gradient(i_coord,i_spin) + &
                         atomic_gradient * dir_tab(i_coord, i_center_L)

                 enddo

              enddo
           endif
        enddo

     end if

  end if

end subroutine evaluate_free_atom_sums_p2
!---------------------------------------------------------------------
!******
