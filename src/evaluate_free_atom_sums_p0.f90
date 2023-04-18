!****s* FHI-aims/evaluate_free_atom_sums_p0
!  NAME
!   evaluate_free_atom_sums_p0
!  SYNOPSIS

subroutine evaluate_free_atom_sums_p0 &
     ( dist_tab, i_r, dir_tab, &
     free_hartree_superpos, free_rho_superpos, &
     free_rho_gradient_superpos, &
     rho, rho_gradient, n_atom_list, atom_list &
     )

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
  implicit none

!  ARGUMENTS

  integer:: n_atom_list
  integer:: atom_list(n_atom_list)
  real*8, dimension(n_atom_list) :: dist_tab
  real*8, dimension(n_atom_list) :: i_r
  real*8, dimension(3,n_atom_list) :: dir_tab
  real*8 :: free_hartree_superpos 
  real*8 :: free_rho_superpos 
  real*8, dimension(3) :: free_rho_gradient_superpos 
  real*8, dimension(n_spin) :: rho
  real*8, dimension(3,n_spin) :: rho_gradient


!  INPUTS
!    o n_atom_list -- number of atoms
!    o atom_list -- list of atoms
!    o dist_tab -- distance to atoms
!    o i_r -- the distance from current integration point to all atoms
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

  integer :: i_center_2, i_center_L
  integer :: i_coord
  integer :: i_spin
  integer :: i_atom_2

  !  begin work

  !       initialize
  free_rho_superpos = 0.d0
  free_hartree_superpos = 0.d0
  if (use_density_gradient) then
     free_rho_gradient_superpos = 0.d0
  end if
  rho = 0.d0

  !       Contributions from all atoms

  do i_center_L = 1, n_atom_list, 1
! We need to exlude free rho of pseudoized species
     if (species_pseudoized(species(atom_list(i_center_L)))) cycle

     i_center_2 = atom_list(i_center_L)

     if (.not.empty(center_to_atom(i_center_2))) then
        !         electrostatic potential contribution

        if (dist_tab(i_center_L).le. &
             multipole_radius_free(species_center(i_center_2))) then

           free_hartree_superpos =  &
                free_hartree_superpos +  &
                val_spline &
                ( i_r(i_center_L), free_pot_es_spl(1,1,species_center(i_center_2)), &
                n_grid(species_center(i_center_2)) )

           !            density contribution
           free_rho_superpos = free_rho_superpos +  &
                val_spline &
                ( i_r(i_center_L), free_rho_spl(1,1,species_center(i_center_2)), &
                n_grid(species_center(i_center_2)) )


           if (use_density_gradient) then
              ! density gradient contribution
              ! factor pi4 already divided out!! different to free_rho_superpos!!!

              atomic_gradient = &
                   val_spline( i_r(i_center_L), &
                   free_drho_dr_spl(1,1,species_center(i_center_2)), &
                   n_grid(species_center(i_center_2)) )

              atomic_gradient = atomic_gradient * pi4_inv

              do i_coord = 1,3,1
                 free_rho_gradient_superpos(i_coord) = &
                      free_rho_gradient_superpos(i_coord) + &
                      atomic_gradient * dir_tab(i_coord, i_center_L)
              enddo

           end if

        end if
     end if
  enddo

  if (.not.use_initial_rho) then
     ! non-polarized case

     i_spin = 1
     rho(i_spin) = pi4_inv * free_rho_superpos

  else

     do i_spin = 1, n_spin, 1

        do i_center_L = 1, n_atom_list, 1

           i_center_2 = atom_list(i_center_L)

           i_atom_2 = center_to_atom(i_center_2)
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

           i_center_2 = atom_list(i_center_L)

           i_atom_2 = center_to_atom(i_center_2)

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

end subroutine evaluate_free_atom_sums_p0
!---------------------------------------------------------------------
!******
