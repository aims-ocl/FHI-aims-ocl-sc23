!****h* FHI-aims/read_fixed_grid
!  NAME
!   read_fixed_grid
!  SYNOPSIS

module read_fixed_grid

!  PURPOSE
!  The module contains routines for reading so called fixed grids.
!
!  USES

  use constants  
  implicit none

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

!******	
	


contains
!----------------------------------------------------------------------------------
!****s* read_fixed_grid/read_fixed_grid_radial_data
!  NAME
!   read_fixed_grid_radial_data
!  SYNOPSIS

  subroutine read_fixed_grid_radial_data(species_name,  angular_limit, n_radial, radius)

!  PURPOSE
!  The subroutine reads radial grid part information for fixed grids.
!
!  USES
!  ARGUMENTS

    real*8:: radius
    character*20 :: species_name
    character*25 :: file_name
    integer:: n_radial, angular_limit

! INPUTS
!  o species_name -- species where grid belongs
!   
!  OUTPUT
!  o angular_limit -- maximum angular grid size
!  o n_radial -- number of radial points
!  o radius -- maximum radius of spherical grid
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



    species_name = adjustl(species_name)

    write(file_name,*) trim(species_name), '.grid'
    file_name = adjustl(file_name)
    open(88, file=trim(file_name))


    !    write(88,*) 'testi'
    !    close(88)
    !    stop
    !  open(88, file=species_name)

    read(88,*) n_radial, radius

    read(88,*) angular_limit

    !  do i_radial = 1, n_radial,1
    !     read(88,*) angular_limit(i_species)
    !  end if



    close(88)

    radius = radius/bohr

  end subroutine read_fixed_grid_radial_data
!******
!----------------------------------------------------------------------------
!****s* read_fixed_grid/fixed_grids_open_for_angular_data
!  NAME
!   fixed_grids_open_for_angular_data
!  SYNOPSIS

  subroutine fixed_grids_open_for_angular_data(species_name, file_number)

!  PURPOSE
!  The subroutine opens the file and reads radial grid part information for fixed grids from file.
!  The file is not closed.
!
!  USES
!  ARGUMENTS

    character*20 :: species_name
    character*25 :: file_name


!  INPUTS
!  o species_name -- species where grid belongs
!  o file_number -- the chanel where the file is opened and not closed.
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE






    real*8:: radius
    integer:: n_radial, angular_limit
    integer:: file_number


    species_name = adjustl(species_name)
    write(file_name,*) trim(species_name), '.grid'
    file_name = adjustl(file_name)
    open(88, file=trim(file_name))

    !  open(file_number, file=species_name)

    read(file_number,*) n_radial, radius

    read(file_number,*) angular_limit




  end subroutine fixed_grids_open_for_angular_data
!******





end module read_fixed_grid
