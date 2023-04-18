!****s* FHI-aims/determine_vacuum_level
!  NAME
!  determine_vacuum_level
!  SYNOPSIS
! This routine is intented to automatically determine
! the position of the vacuum level. It will try
! to find the vacuum as the biggest distance in z-direction
! between atoms and place the vacuum level in the middle of it

!Working principle: For each atom, determine the
! SMALLEST distance in z to the other atoms
!The atom pair which has the largest smallest distance
!determines the position and the size of the vacuum 
subroutine determine_vacuum_level( )
  !  PURPOSE
  !  Determine the position of the vacuum level
  !
  !  USES
  !
  use dimensions
  use runtime_choices
  use geometry
  use localorb_io
  use mpi_tasks, only: aims_stop

    !  INPUTS
  !   none
  !  OUTPUT
  !   none
  ! implicit: sets vacuum_z_level
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

 implicit none

  !local variables
  integer :: i_atom
  integer :: i_atom2
  integer :: last_atom=1 !last atom before the vacuum
  real*8 :: next_in_z(n_atoms) !distance from (i_atom) to next atom in +z-direction
  real*8 :: current_dist = 0 !working variable, current 'vacuum' distance
  integer, dimension(n_atoms) :: vacuum_atom_pair ! for debug: deterimes which atom gives the smallest distance..
  real*8 :: coord_tmp(3)
  real*8 :: coord_tmp2(3)
  character*100 :: info_str

  !begin work
  next_in_z(:)=maxval(lattice_vector(3,:))

  write(info_str,*) 'Determine position of vacuum level automatically' 
  call localorb_info(info_str)
  !OTH: Once it works, remove warning
  write(info_str,*) '**WARNING** : This feature is experimental and still requires some testing' 
  call localorb_info(info_str)
  write(info_str,*) '**PLEASE ** : If you see it fail, please report.'
  call localorb_info(info_str)

  !Determine the distance to the next atom in z-direction
  do i_atom = 1,n_atoms,1
    do i_atom2 = 1,n_atoms, 1
       if (i_atom==i_atom2) cycle !no self-distance
       
       coord_tmp(:) = coords(:,i_atom)
       coord_tmp2(:) = coords(:,i_atom2)
       call map_to_center_cell(coord_tmp)
       call map_to_center_cell(coord_tmp2)

       current_dist=coord_tmp2(3)-coord_tmp(3)
       if (current_dist.le.0) current_dist=current_dist+maxval(lattice_vector(3,:)) !we only look for distances upward
       if (current_dist.lt.next_in_z(i_atom)) then
          next_in_z(i_atom)=current_dist
          vacuum_atom_pair(i_atom)=i_atom2
       endif
    enddo
  enddo

  !Debug block
  !do i_atom=1,n_atoms
  !  write(use_unit,*) i_atom, next_in_z(i_atom)
  !enddo
  !write(use_unit,*) 'Maximum minimal distance: ', maxval(next_in_z(:))

  !Find out which atom has the largest minimal distance to the next neighout
  do i_atom=1,n_atoms
    if(next_in_z(i_atom).eq.maxval(next_in_z(:))) last_atom=i_atom
  enddo

  !Catch corner case - all atoms in one plane
  if (maxval(next_in_z(:))==0) next_in_z = maxval(lattice_vector(3,:))

  !Safeguard vs. nonsense results:
  if (maxval(next_in_z(:)).lt.10) then
    write(info_str,*) 'Vacuum is smaller than 10 Bohr.'
    call localorb_info(info_str)
    write(info_str,*) 'No sensible vacuum level can be set.'
    call localorb_info(info_str)
    write(info_str,*) 'Please set the vacuum level manually or increase the size of the vacuum level'
    call localorb_info(info_str)
    call aims_stop()
  endif

  !Now set vacuum level
  vacuum_z_level=coords(3,last_atom)+next_in_z(last_atom)/2
  if (vacuum_z_level.gt.maxval(lattice_vector(3,:))) then
      vacuum_z_level=vacuum_z_level-maxval(lattice_vector(3,:))
  endif

   write(info_str,*) 'Vacuum level determined: '
    call localorb_info(info_str)
   write(info_str,*) 'set_vacuum_level ', vacuum_z_level*bohr
   call localorb_info(info_str)

   !Debug block
   !write(info_str,*) 'Atom determining vacuum: ', last_atom
   !call localorb_info(info_str)
   !write(info_str,*) 'paired atom: ', vacuum_atom_pair(last_atom)
   !call localorb_info(info_str)
   !write(info_str,*) 'Distance: ', next_in_z(last_atom)*bohr
   !call localorb_info(info_str)
   !call aims_stop(info_str)

end subroutine determine_vacuum_level
