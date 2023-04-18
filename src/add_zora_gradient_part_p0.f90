!****s* FHI-aims/add_zora_gradient_part_p0
!  NAME
!    add_zora_gradient_part_p0
!  SYNOPSIS

 subroutine add_zora_gradient_part_p0(gradient_part, i_r,dir,dist_tab, zora, n_atom_list, atom_list)

 
!  PURPOSE
!   The subroutine adds gradients term of ZORA to gradient_part.   
!    
!  USES

   use dimensions
   use runtime_choices
   use grids
   use pbc_lists
   use spline
   use free_atoms
   use constants
   implicit none

!  ARGUMENTS

   real*8  :: gradient_part(3)
   real*8  :: i_r(n_atom_list)
   real*8  :: dir(3,n_atom_list)
   real*8  :: dist_tab(n_atom_list)
   real*8  :: zora
   integer :: n_atom_list
   integer :: atom_list(n_atom_list)

!  INPUTS
!  o i_r -- distance to atoms, in relative position of log grid.
!  o dir -- which direnction atoms are
!  o dist_tab -- distance to atoms
!  o zora -- ZORA factor, typically: (light_speed_sq / (2 * light_speed_sq - local_potential)**2)
!  o n_atom_list -- number of atoms, including periodic mirror images
!  o atom_list -- atoms, including periodic mirror images.xs
!
!  OUTPUT
!   o gradient_part -- gradient term where the ZORA part is added.
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
! SOURCE




   real*8:: V_radial_deriv
   real*8:: gradient(3)
   integer:: i_center_2, i_center_L

   gradient = 0.0

   do i_center_L = 1,n_atom_list,1
      i_center_2 = atom_list(i_center_L)
      
      V_radial_deriv =  val_spline_deriv( i_r(i_center_L), free_potential_spl(1,1,species_center(i_center_2)), &
           n_grid(species_center(i_center_2)) ) / &
           (log(r_grid_inc(species_center(i_center_2))) *  dist_tab(i_center_L))

      gradient(:) = gradient(:)  +  abs(V_radial_deriv) * dir(:,i_center_L)

   end do

!
!
!   zora:  (light_speed_sq / (2 * light_speed_sq - local_potential)**2)
!    
   gradient_part(:) = gradient_part(:) - zora *  gradient(:) 



 end subroutine add_zora_gradient_part_p0
!******
