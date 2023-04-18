!****s* FHI-aims/initialize_fixed_grids
!  NAME
!   initialize_fixed_grids
!  SYNOPSIS

subroutine initialize_fixed_grids()

!  PURPOSE
!  The subroutine initializes the user defined grids, so called fixed grids.
!  These grids are not adaptivily optimised, but are identical in all same species
!  atoms.
!
!  USES

  use dimensions
  use grids
  use species_data
  use mpi_tasks
  use localorb_io, only : use_unit
  implicit none

!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
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


	

  integer :: i_species, i_radial, i_grid_shell

  integer, dimension(:,:), allocatable :: n_lebedev


  if (.not.allocated(n_lebedev)) then
     allocate( n_lebedev(n_max_radial, n_species),stat=i_species )
     call check_allocation(i_species, 'n_lebedev                     ')
  end if

  lebedev_grid_index = 0


  do i_species = 1, n_species, 1

     do i_radial = 1, n_radial(i_species), 1
        
        i_grid_shell = 1
        do while ( (i_grid_shell.lt.n_ang_shells(i_species)) .and. &
        &          (r_ang_shell(i_grid_shell,i_species) .lt.&
        &           r_radial(i_radial, i_species)) )
           i_grid_shell = i_grid_shell + 1
        enddo
        
        n_lebedev(i_radial,i_species) = grid_ceil(n_ang_points(i_grid_shell,i_species))
        
        n_angular(i_radial,i_species) = &
             n_angular_lebedev(n_lebedev(i_radial,i_species))
        r_angular(:,:,i_radial,i_species) = &
             r_angular_lebedev(:,:,n_lebedev(i_radial,i_species))
        w_angular(:,i_radial,i_species) = &
             w_angular_lebedev(:,n_lebedev(i_radial,i_species))
        
        if (n_angular(i_radial,i_species).gt.n_max_angular) then
           ! smallest grid is already too large, abort calculation
           write(use_unit,'(1X,A,I5,A,I5,A)') "* Species ", i_species, " radial shell ",    i_radial, ":"
           write(use_unit,'(1X,A)')           "* No suitable angular grid available. "
           write(use_unit,'(1X,A,A)')         "* Adjust the prescribed number of grid points in control.in."
           stop
        end if

        lebedev_grid_index(i_radial,i_species) = n_lebedev(i_radial,i_species)

     end do
  end do

  if (allocated(n_lebedev)) then
     deallocate( n_lebedev )
  end if

end subroutine initialize_fixed_grids
!******	
