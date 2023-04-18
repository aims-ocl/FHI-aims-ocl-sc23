!****s* FHI-aims/initialize_grid_storage_p1
!  NAME
!   initialize_grid_storage_p1
!  SYNOPSIS

    subroutine initialize_cube_grid_storage( &
          n_points, cube_points, &
          my_free_hartree_superpos, my_free_rho_superpos, &
          my_free_rho_gradient_superpos, &
          my_rho, my_rho_gradient  &
    )

!  PURPOSE
!   calculates relevant quantities for integration on a supplied grid, as well as initializing
!   a number of free atom references.
!  WARNING: This routine has been written with free_hartree_superpos in mind
!  WARNING: If you intend to use anything else, first make sure it is correctly calculated
!  WARNING: This is particular true for rho_gradient, where only a dummy is supplied.
!
!  USES
!
      use dimensions
      use grids
      use runtime_choices
      use geometry
      use spline
      use free_atoms
      use pbc_lists
      use species_data
      use localorb_io

implicit none

!  ARGUMENTS
       integer,intent(in) :: n_points !number of points
       real*8,intent(in) :: cube_points(3,n_points) !coordinates of grid points

      real*8, dimension(n_points) :: my_free_hartree_superpos 
      real*8, dimension(n_points) :: my_free_rho_superpos
      real*8, dimension(3, n_points) :: my_free_rho_gradient_superpos
      real*8, dimension(n_spin,n_points) :: my_rho
      real*8, dimension(3, n_spin, n_points) :: my_rho_gradient

! INPUTS
!  o n_points !number of points
!  o cube_points ! coordinates of points
! OUTPUTS
! o partition_tab -- grid integration weight
! o hartree_potential -- overall hartree potential (including any external ions)
! o free_hartree_superpos -- superposition of hartree potential of free atoms; reference
! o free_rho_superpos -- superposition of free atom densities
! o pot_ion_embed -- ion embedding potential
! o rho -- density
! o free_rho_gradient_superpos -- free atom density gradient
! o rho_gradient -- density gradient
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

!Auxilary variable
  character*100 :: info_str
  character*30, parameter :: func = 'initialize_cube_grid'

      real*8 :: dist_tab_sq(n_centers_basis_integrals)
      real*8 :: dir_tab(3,n_centers_basis_integrals)
      real*8 :: dist_tab(n_centers_basis_integrals)
      real*8 :: dir_tab_norm(3,n_centers_basis_integrals)
      real*8 :: i_r(n_centers_basis_integrals)
      integer :: n_compute_atoms, n_compute_occ_atomS
      integer, dimension(n_centers_basis_integrals) :: center_index
      integer, dimension(n_centers_basis_integrals) :: i_occ2i_compute
      real*8, dimension(n_centers_basis_integrals) :: temp_occ_rho
      real*8 :: local_rho_gradient(3,n_spin)
!       real*8 :: dummy !SR: if the laplacian of rho_free is wished uncomment this and the comment below


!local counter
   integer :: i_point, i_atom
   real*8 :: coord_current(3)


!here begins the work
  
    do i_point=1,n_points,1
          coord_current(:) = cube_points(:,i_point)
        !              tabulate current integration point as it appears from spherical
        !              coordinates centered at each atom
            call tab_atom_centered_coords_p0 &
              ( coord_current, &
              dist_tab_sq, &
              dir_tab, &
              n_centers_basis_integrals, centers_basis_integrals )


          n_compute_atoms=n_centers_basis_integrals !Uppermost estimate
          n_compute_occ_atoms= n_centers_basis_integrals !Also uppermost estimate
          do i_atom = 1,n_centers_basis_integrals,1
             center_index(i_atom)=i_atom
             i_occ2i_compute(i_atom) = i_atom
          enddo
          

 
         !Get direction and normalized distances
           call tab_global_geometry_p0 &
                ( dist_tab_sq,         &
                  dir_tab,             &
                  dist_tab,            &
                  i_r,                 &
                  dir_tab_norm,        &
                  n_compute_atoms,     &
                  center_index )



         !Initialize free atom variables
            call evaluate_free_atom_sums_p2  &
                ( i_r, dir_tab_norm,  &
                my_free_hartree_superpos(i_point),  &
                my_free_rho_superpos(i_point),  &
                my_free_rho_gradient_superpos(1, i_point), &
                my_rho(1,i_point),  &
                local_rho_gradient, n_compute_atoms, center_index, &
                n_compute_occ_atoms, i_occ2i_compute, temp_occ_rho) !,dist_tab_sq,&
                !dummy)              



    enddo

end subroutine initialize_cube_grid_storage
