!****s* FHI-aims/evaluate_partition_tab_p0
!  NAME
!   evaluate_partition_tab_p0
!  SYNOPSIS

subroutine evaluate_partition_tab_p0 &
     ( i_center,  i_center_L, dist_tab, i_r, &
     radial_weight, angular_weight, &
     partition_tab,&
     n_atom_list, atom_list  &
     )

!  PURPOSE
!     Subroutine evaluate_partition
!     calculates the partition function the present integration point.
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

  integer:: n_atom_list
  integer:: atom_list(n_atom_list)
  integer i_center
  integer i_center_L
  real*8 dist_tab(n_atom_list)
  real*8 i_r(n_atom_list)
  real*8 radial_weight
  real*8 angular_weight
  real*8 partition_tab

!  INPUTS
!  o n_atom_list -- number of relevant atoms
!  o atom_list -- list of relevant atoms
!  o i_center -- center where partition function belongs, absolute index
!  o i_center_L -- center where partition function belongs, relative index
!  o dist_tab -- distance to relevant atoms
!  o i_r -- the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!  o radial_weight -- radial integration weight
!  o angular_weight -- angular integration weight
!  
!  OUTPUT
!  o partition_tab -- value of partition function
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


	


  logical :: valid

  !     variables for trial partition function
  real*8 aux_dens
  real*8 wt_dens
  real*8 partition_norm

  !     counters

  integer i_center_2, i_center_L2

  !     begin work

  ! First, create safety net: Any integration point that lies inside the innermost logarithmic
  ! grid shell of another atom must be avoided because else we would have to extrapolate to
  ! zero on the logarithmic grid (but zero is the point at -infty on a log. grid, and therefore this
  ! extrapolation would have horrible consequences ...)

  valid = .true.
  do i_center_L2 = 1, n_atom_list, 1 
     i_center_2 = atom_list(i_center_L2)

     if (dist_tab(i_center_L2).lt.r_grid_min(species_center(i_center_2))) then
        valid = .false.
        partition_tab = 0.d0
     end if
  enddo

  ! continue only if the grid point is safe
  if (valid) then

     ! initialize numerator / denominator of partition function

     if ( (flag_hartree_partition_type.ne.partition_type) ) then
       aux_dens = &
          val_spline &
          ( i_r(i_center_L), partition_rho_spl(1,1,species_center(i_center)), &
          n_grid(species_center(i_center)) )  
     else
       aux_dens = &
          val_spline &
          ( i_r(i_center_L), hartree_partition_rho_spl(1,1,species_center(i_center)), &
          n_grid(species_center(i_center)) )  
     end if

     select case(partition_type)

     case(1)
        wt_dens = &
             aux_dens / (dist_tab(i_center_L)**2.0d0)
     case(2)
        wt_dens = &
             aux_dens / abs(dist_tab(i_center_L))
     case(3)
        wt_dens = &
             aux_dens 
     case(4)
        wt_dens = &
             1.0/ ( 1+ &
             exp( (dist_tab(i_center_L)- hartree_partition_parameters(1))/ &
             hartree_partition_parameters(2) ) )
     case(6)
        wt_dens = &
             aux_dens / (dist_tab(i_center_L)**2.0d0)
     end select


     !      wt_dens = aux_dens / 
     !     +  (dist_tab(i_center))**2.0d0

     !     obtain sum over all density weights at current integration point for normalisation.
     partition_norm = wt_dens
     do i_center_L2 = 1, n_atom_list,1
        i_center_2 = atom_list(i_center_L2)

        if (i_center_2.ne.i_center) then

           !         get density of i_center_2 at current integration point
           if ( (flag_hartree_partition_type.ne.partition_type) ) then
             aux_dens = &
                val_spline &
                ( i_r(i_center_L2), partition_rho_spl(1,1,species_center(i_center_2)), &
                n_grid(species_center(i_center_2)) ) 
           else
             aux_dens = &
                val_spline &
                ( i_r(i_center_L2), hartree_partition_rho_spl(1,1,species_center(i_center_2)), &
                n_grid(species_center(i_center_2)) ) 
           end if

           !         add contribution from i_center_2 to partition_norm


           select case(partition_type)

           case(1)
              partition_norm = partition_norm + &
                   aux_dens / (dist_tab(i_center_L2)**2.0d0)
           case(2)
              partition_norm = partition_norm + &
                   aux_dens / abs(dist_tab(i_center_L2))
           case(3)
              partition_norm = partition_norm + &
                   aux_dens 
           case(4)
              partition_norm = partition_norm + &
                   1.0/ ( 1+  &
                   exp((dist_tab(i_center_L2)-hartree_partition_parameters(1))/ &
                   hartree_partition_parameters(2) ) )
           case(6)
              partition_norm = partition_norm + &
                   aux_dens / (dist_tab(i_center_L2)**2.0d0)
           end select


        end if

        !     end loop over other atoms, i_center_2
     enddo


     if (partition_norm.gt.partition_acc) then
        !       tabulate partition function
        partition_tab = &
             wt_dens / partition_norm

        !       and multiply integration weights into partition_tab
        partition_tab = &
             partition_tab * &
             radial_weight * &
             (dist_tab(i_center_L))**2.d0 * &
             angular_weight * &
             4.0d0*pi



     else
        !       this grid point is too far away - we should not be integrating here!
        partition_tab = 0.d0
     end if

     ! end if (valid)
  end if

end subroutine evaluate_partition_tab_p0

!---------------------------------------------------------------------
!******	
