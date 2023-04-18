!****s* FHI-aims/write_cube_header
!  NAME
!    write_cube_header
!  SYNOPSIS

      subroutine write_cube_header(descriptor,offset,num,i_cube)

!  PURPOSE
!   Writes a header to cube format file.
!
!  USES
!
      use dimensions
      use runtime_choices
      use grids
      use geometry
      use species_data
      use mpi_utilities
      use constants
      use plot
      use free_atoms, only: average_free_es_pot   
      implicit none

!  ARGUMENTS

      integer  descriptor
      integer  num(3)
      real*8   offset(3)
!      real*8   scaling(3)
      integer :: i_cube
      


!  INPUTS
!   o descriptor -- file number where the result are written
!   o num -- number of grid points
!   o offset -- offset of the coordinates
!   o scaling -- what is a scale between grid poitns and coordinates
!
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
!
!     local
      integer i_coord
      integer i_atom
      real*8  local_offset(3)
      real*8  coords_temp(3,n_atoms)
      real*8, dimension(3,n_atoms) :: periodic_unit_cell_translations_cart
  

    coords_temp = coords
    if ((output_in_original_unit_cell).and.(n_periodic>0)) then

       call frac2cart (lattice_vector,                    &
                       periodic_unit_cell_translations,   &
                       periodic_unit_cell_translations_cart)

       coords_temp(:,:) = coords(:,:)+periodic_unit_cell_translations_cart(:,:)

       call cart2frac (lattice_vector,                       &
                       periodic_unit_cell_translations_cart, &
                       periodic_unit_cell_translations)
    end if
    
    if (flag_out_elec_real) then 

      write (descriptor,'(A)') "#Real-sapce electrostatic potential written by FHI-AIMS in a cube file format"
      write (descriptor,'(2(A,F15.8))') "#Numerical average free-atom electrostatic potential in [eV]:", average_free_es_pot*hartree,& 
            ";   Average real-space part of the electrostatic potential in [eV]:", average_potential*hartree 
    else
      write (descriptor,*) "CUBE FILE written by FHI-AIMS"
      write (descriptor,*) "*****************************"
    endif      
      local_offset = -offset

      write (descriptor,fmt='(1X,I4,3F12.6)') &
       n_atoms, (local_offset(i_coord), i_coord=1,3,1)

         write(descriptor,fmt='(1X,I4,3F12.6)') &
           cube_edge_steps(1,i_cube), cube_edge_unit(1:3,1,i_cube)
         write(descriptor,fmt='(1X,I4,3F12.6)') &
           cube_edge_steps(2,i_cube), cube_edge_unit(1:3,2,i_cube)
         write(descriptor,fmt='(1X,I4,3F12.6)') &
           cube_edge_steps(3,i_cube), cube_edge_unit(1:3,3,i_cube)


      do i_atom = 1,n_atoms,1
         write (descriptor,fmt='(1X,I4,4F12.6)') &
         nint(species_z(species(i_atom))), 0.0d0, &
         (coords_temp(i_coord,i_atom),i_coord=1,3,1)
      enddo

      end subroutine write_cube_header
!****************************************
!---------------------------------------------------------------
!****s* FHI-aims/write_cube_header_gOpenMol
!  NAME
!    write_cube_header_gOpenMol
!  SYNOPSIS

subroutine write_cube_header_gOpenMol(descriptor, cube_units,i_cube)

!  PURPOSE
!  Writes cube-files header for gOpenMol format
!
!  USES
!TODO: CHECK if all these are necessary - I highly doubt that, given the short routine.!
  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use mpi_utilities
  use constants
  use plot
  implicit none

!  ARGUMENTS

  integer  descriptor
  real*8 cube_units(3)
  integer :: i_cube

!  INPUTS
!   o descriptor -- file number where results are written
!   o cube_units -- scale of the units
!   
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
!


  !     local
  integer i_atom
  real*8 coords_temp(3,n_atoms)
  real*8, dimension(3,n_atoms) :: periodic_unit_cell_translations_cart

    coords_temp = coords
    if ((output_in_original_unit_cell).and.(n_periodic>0)) then

       call frac2cart (lattice_vector,                    &
                       periodic_unit_cell_translations,   &
                       periodic_unit_cell_translations_cart)

       coords_temp(:,:) = coords(:,:)+periodic_unit_cell_translations_cart(:,:)

       call cart2frac (lattice_vector,                       &
                       periodic_unit_cell_translations_cart, &
                       periodic_unit_cell_translations)
    end if


  write (descriptor,*) '3  2'

  write (descriptor,'(3I4)')   cube_edge_steps(1,i_cube), cube_edge_steps(2,i_cube), cube_edge_steps(3,i_cube)

  write (descriptor,fmt='(6F12.6)') &
       (-cube_units(1)* cube_edge_steps(1,i_cube)/2 + cube_origin(1,i_cube))*bohr, &
       (cube_units(1)* cube_edge_steps(1,i_cube)/2 + cube_origin(1,i_cube))*bohr, &
       (-cube_units(2)* cube_edge_steps(2,i_cube)/2 + cube_origin(2,i_cube))*bohr, &
       (cube_units(2)* cube_edge_steps(2,i_cube)/2 + cube_origin(2,i_cube))*bohr, &
       (-cube_units(3)* cube_edge_steps(3,i_cube)/2 + cube_origin(3,i_cube))*bohr, &
       (cube_units(3)* cube_edge_steps(3,i_cube)/2 + cube_origin(3,i_cube))*bohr

  open(88,file='coords.xyz')

  write(88,*) n_atoms
  write(88,*) 'cell'

     do i_atom = 1, n_atoms
        write(88,'(A,3F12.6)') trim(species_name(species(i_atom))), &
           coords_temp(3,i_atom),coords_temp(2,i_atom),coords_temp(1,i_atom)
     end do


  close(88)

end subroutine write_cube_header_gOpenMol
!*******************************************************
!---------------------------------------------------------------
!****s* FHI-aims/write_cube_header_xsf
!  NAME
!    write_cube_header_xsf
!  SYNOPSIS

subroutine write_cube_header_xsf(descriptor, i_cube)

!  PURPOSE
!  Writes cube-files header for gOpenMol format
!
!  USES
!
      use dimensions
      use runtime_choices
      use grids
      use geometry
      use species_data
      use mpi_utilities
      use constants
  use plot
  implicit none

!  ARGUMENTS

  integer ::  descriptor
  integer :: i_cube 
!  INPUTS
!   o descriptor -- file number where results are written
!   o i_cube -- number of the cube file
!   
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


!variables

!counters
integer :: i_atom
integer :: i_coord
real*8, dimension(3,n_atoms) :: coords_temp
real*8, dimension(3,n_atoms) :: periodic_unit_cell_translations_cart
    coords_temp = coords
    if ((output_in_original_unit_cell).and.(n_periodic>0)) then

       call frac2cart (lattice_vector,                    &
                       periodic_unit_cell_translations,   &
                       periodic_unit_cell_translations_cart)

       coords_temp(:,:) = coords(:,:)+periodic_unit_cell_translations_cart(:,:)

       call cart2frac (lattice_vector,                       &
                       periodic_unit_cell_translations_cart, &
                       periodic_unit_cell_translations)
    end if




! begin work
         if (n_periodic.eq.0) then
            write(descriptor,*) 'MOLECULE'
            write(descriptor,*) 'ATOMS' 
         elseif (n_periodic.eq.3) then
            write(descriptor,*) 'CRYSTAL'
            write(descriptor,*) 'PRIMVEC'
            write(descriptor,fmt='(1X,3F12.6)') lattice_vector(1:3,1)*bohr
            write(descriptor,fmt='(1X,3F12.6)') lattice_vector(1:3,2)*bohr
            write(descriptor,fmt='(1X,3F12.6)') lattice_vector(1:3,3)*bohr
            write(descriptor,*) 'CONVVEC'
            write(descriptor,fmt='(1X,3F12.6)') lattice_vector(1:3,1)*bohr
            write(descriptor,fmt='(1X,3F12.6)') lattice_vector(1:3,2)*bohr
            write(descriptor,fmt='(1X,3F12.6)') lattice_vector(1:3,3)*bohr
            write(descriptor,*) 'PRIMCOORD'
            write(descriptor,fmt='(1X,I4,A)') n_atoms, '  1'
         else
            call aims_stop('Open boundary conditions not implemented for xsf headers', 'write_cube_header_xsf') 
         endif

      do i_atom = 1,n_atoms,1
         write(descriptor,fmt='(1X,I4,4F12.6)') &
         nint(species_z(species(i_atom))), &
         (coords_temp(1:3,i_atom)*bohr)
      enddo
  
         write(descriptor,*) 'BEGIN_BLOCK_DATAGRID_3D'
         write(descriptor,*) 'g98_3D_unknown'
         write(descriptor,*) 'DATAGRID_3D_g98Cube'
         
         write(descriptor,fmt='(1X,3I4)') &
          cube_edge_steps(3,i_cube), &
          cube_edge_steps(2,i_cube), &
          cube_edge_steps(1,i_cube)

         write(descriptor,fmt='(1X,3F12.6)') &
           (cube_origin(1,i_cube)- (0.5*cube_edge_unit(1,1,i_cube)*(cube_edge_steps(1,i_cube)-1) ) &
                                 - (0.5*cube_edge_unit(1,2,i_cube)*(cube_edge_steps(2,i_cube)-1) ) &
                                 - (0.5*cube_edge_unit(1,3,i_cube)*(cube_edge_steps(3,i_cube)-1) ))*bohr , &

           (cube_origin(2,i_cube)- (0.5*cube_edge_unit(2,1,i_cube)*(cube_edge_steps(1,i_cube)-1) ) &
                                 - (0.5*cube_edge_unit(2,2,i_cube)*(cube_edge_steps(2,i_cube)-1) ) &
                                 - (0.5*cube_edge_unit(2,3,i_cube)*(cube_edge_steps(3,i_cube)-1) ))*bohr , &

           (cube_origin(3,i_cube)- (0.5*cube_edge_unit(3,1,i_cube)*(cube_edge_steps(1,i_cube)-1) ) &
                                 - (0.5*cube_edge_unit(3,2,i_cube)*(cube_edge_steps(2,i_cube)-1) ) &
                                 - (0.5*cube_edge_unit(3,3,i_cube)*(cube_edge_steps(3,i_cube)-1) ))*bohr  

         write(descriptor,fmt='(1X,3F12.6)') &
           cube_edge_unit(1:3,3,i_cube)*(cube_edge_steps(3,i_cube)-1)*bohr
         write(descriptor,fmt='(1X,3F12.6)') &
           cube_edge_unit(1:3,2,i_cube)*(cube_edge_steps(2,i_cube)-1)*bohr
         write(descriptor,fmt='(1X,3F12.6)') &
           cube_edge_unit(1:3,1,i_cube)*(cube_edge_steps(1,i_cube)-1)*bohr


end subroutine write_cube_header_xsf


