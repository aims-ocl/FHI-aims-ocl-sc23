!****s* FHI-aims/output_whole_potential_on_arbitrary_points
!  NAME
!    output_whole_potential_on_arbitrary_points
!  SYNOPSIS

    subroutine output_whole_potential_on_arbitrary_points (n_points,cube_points,cube_potential)

!  PURPOSE
!  High-level wrapper to prepare all variables needed 
!  for output of potential.
!  To be abandoned and incorporated into calculate_whole_potential later on
!
!  USES

!TODO: Check which are needed and cleanup
      use localorb_io
      use dimensions
      use runtime_choices
      use physics
      implicit none

!  ARGUMENTS
!
!    none
!
!
!  INPUTS
!     integer :: n_points !number of points 
!     real*8 :: cube_points(3,n_points) !coordiantes of the points
     integer :: n_points !number of points 
     real*8 :: cube_points(3,n_points) !coordiantes of the points
!  OUTPUT
     real*8 :: cube_potential(*) !potential on given coordinates
!    none
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




! Aux variables
      character*100 :: info_str
      character(*), parameter :: func = 'scf_solver'
      integer :: i_point, i_coord
      real*8 :: my_free_hartree_superpos(n_points)
      real*8, dimension(n_points) :: my_free_rho_superpos
      real*8, dimension(n_spin,n_points) :: my_rho
      real*8, dimension(3, n_points) :: my_free_rho_gradient_superpos
      real*8, dimension(3, n_spin, n_points) :: my_rho_gradient
      real*8, dimension(n_points) :: my_partition_tab

! Initialisation of auxiliary variables
      i_point=0
      i_coord=0
      my_free_hartree_superpos=0
      my_free_rho_superpos=0
      my_rho=0
      my_free_rho_gradient_superpos=0
      my_rho_gradient=0
      my_partition_tab(1:n_points)=1.00

!Initializing grid and relevant quantities, 
!in particular my_free_hartree_superpos
      call  initialize_cube_grid_storage( &
                   n_points, cube_points, &
                   my_free_hartree_superpos, my_free_rho_superpos,  &
                   my_free_rho_gradient_superpos, &
                   my_rho, my_rho_gradient  &
                                                )
!Calculate potential on grid
      call  calculate_whole_potential_on_arbitrary_points &
                (n_points,cube_points,cube_potential, &
                 delta_v_hartree_part_at_zero, &
                delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
                my_partition_tab, my_rho, &
                my_free_hartree_superpos, my_free_rho_superpos,  &
                hartree_delta_energy, en_elec_delta, hartree_multipole_correction, &
                  en_density_embed, &
                 multipole_radius_sq, &
                l_hartree_max_far_distance, & 
                outer_potential_radius )

!Safeguard against infinite values - however, these point at a
! bug in the implementation! 
      do i_point = 1,n_points,1
        if (cube_potential(i_point)*hartree.lt.-1e20) then 
             write(use_unit,*) 'Too small value found at: ', i_point
             write(use_unit,'(A,3F12.3)') 'Coordinate', cube_points(1:3,i_point)
             cube_potential(i_point)=-1e20/hartree
        endif
        if (cube_potential(i_point)*hartree.gt.1e20) then
             cube_potential(i_point)=1e20/hartree
             write(use_unit,*) 'Too large value found'
        endif
             
      enddo

end subroutine  output_whole_potential_on_arbitrary_points
