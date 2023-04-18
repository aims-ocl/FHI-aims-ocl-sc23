!****s* FHI-aims/ion_pseudopot_potential_v1
!  NAME
!    ion_pseudopot_potential_v1
!  SYNOPSIS
subroutine ion_pseudopot_potential_v1(input_coord, potential,pp_number_excluded)
 
! sumed up local pseudo potential of all pseudo pot except pp_number_excluded 

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
  integer :: i_pp_species, i_pp_atom


  potential = 0d0

  do i_pp_atom = 1, n_pp_atoms, 1
     if(i_pp_atom.eq.pp_number_excluded) exit
     i_pp_species = pp_species(i_pp_atom)
     distance_vector(:) = input_coord(:) - pp_coords(:,i_pp_atom) 
     distance = sqrt((distance_vector(1))**2+(distance_vector(2))**2+ &
                    (distance_vector(3))**2 )


     ilog = invert_log_grid(distance, &
                    pp_r_grid_min(i_pp_species), &
                    pp_r_grid_inc(i_pp_species))


     if (ilog.gt.n_points_pp_fn(i_pp_species)) then
        potential = -pp_charge(pp_species(i_pp_atom))/distance
     else

     potential = potential + &
           val_spline ( ilog, local_pseudopot_spl(:,:,i_pp_species), n_points_pp_fn(i_pp_species) )
     endif

  end do



end subroutine
