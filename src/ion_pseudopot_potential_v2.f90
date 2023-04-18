!****s* FHI-aims/ion_pseudopot_potential_v2
!  NAME
!    ion_pseudopot_potential_v2
!  SYNOPSIS
subroutine ion_pseudopot_potential(input_coord, potential,pp_number_excluded, force)
 
! sumed up local pseudo potential of all pseudo pot except pp_number_excluded 
!  in v_2 its a simple sum over monopoles

  use dimensions
  use geometry
  use pseudodata
  use grids, only: invert_log_grid
  use spline

  
implicit none

! Imported variables

  real*8, dimension(3) :: input_coord
  real*8 :: potential
  integer :: pp_number_excluded


! local variables
  real*8, dimension(3) :: distance_vector
  real*8 :: distance, ilog
  integer :: i_pp_species, i_pp_atom, i_coord

! output
  real*8, dimension(3) :: force


  potential = 0d0
  force = 0.d0

  do i_pp_atom = 1, n_pp_atoms, 1
     if(i_pp_atom.ne.pp_number_excluded) then
     distance=0.d0
     i_pp_species = pp_species(i_pp_atom)
     distance_vector(:) = input_coord(:) - pp_coords(:,i_pp_atom) 
     distance = sqrt((distance_vector(1))**2+(distance_vector(2))**2+ &
                    (distance_vector(3))**2 )

     ! negative sign because positive charges must be attractive!
     potential = potential &
                - pp_charge(i_pp_species)/distance


     if (use_forces) then 
        do i_coord = 1,3,1
          force(i_coord) = force(i_coord) - pp_charge(i_pp_species) &
                           * distance_vector(i_coord) / distance**3.d0
        end do
     endif


     end if
  end do


end subroutine
