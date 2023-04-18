!****s* FHI-aims/evaluate_nlcc_forces
!  NAME
!    evaluate_nlcc_forces
!  SYNOPSIS

!subroutine evaluate_nlcc_forces(nlcc_forces,local_xc_derivs, partition_tab, &
!     dist_tab_sq_global, dir_tab_global, n_points)

subroutine evaluate_nlcc_forces(local_xc_derivs, partition,&
     partial_core_rho_grad_in_batch, dist_tab_sq_global, dir_tab_global, n_points)


!  PURPOSE
!    Sums up the force term related to the nonlinear core correction
!    of KB pseudopotentials for a given batch of integration points.
!     F_i = sum_{batch points} partition_tab(i_point) * V_xc(i_point) * 
!         \frac{d rho_{pc}}{d R_i} * direction-vector
!
!  USES

  use dimensions
  use species_data
  use geometry
  use runtime_choices
  use pseudodata
  use spline
  use grids
!  ARGUMENTS

  implicit none
  real*8, dimension(n_points), intent(in) :: local_xc_derivs
  real*8, dimension(3,n_points) :: rho_grad
  real*8, dimension(3, n_atoms, n_points), intent(in) :: dir_tab_global
  real*8, dimension(n_atoms, n_points), intent(in) :: dist_tab_sq_global
  real*8, dimension(n_points), intent(in) :: partition
  real*8, dimension(3, n_points), intent(in) :: partial_core_rho_grad_in_batch
  integer :: n_points


!  INPUTS
!  o  local_xc_derivs -- exchange correlation potential
!  o  dir_tab_global -- direction to atoms
!  o  dist_tab_sq_global -- (distance to atoms)**2
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!
!  OUTPUT
!  o  nlcc_forces -- Force related to the nonlinear core correction for KB pseudopotentials.
!
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


  integer :: i_point, i_atom, n_atoms_pp
  integer :: i_pp_atom, i_coord, i_pp_species
  real*8 :: dist, i_r_log
  real*8 :: point_term
  real*8 :: d_partial_core_rho_d_R
  real*8 :: delta(3) 


  do i_point = 1, n_points, 1


     point_term = partition(i_point) * local_xc_derivs(i_point)

!  We some up over all pseudocores to get the correction core densities for
!  every single pseudocore 
        do i_pp_atom = 1, n_pp_atoms
           i_pp_species = pp_species(i_pp_atom)


           dist =  sqrt(dist_tab_sq_global(n_atoms_pp + i_pp_atom, i_point))

        d_partial_core_rho_d_R = 0

        if(dist.gt.(localpot_outer_radius(i_pp_species))) then
           d_partial_core_rho_d_R = 0.d0
        elseif (dist.lt.(pp_r_grid_min(i_pp_species))) then
           d_partial_core_rho_d_R = 0.d0
        else


           i_r_log = invert_log_grid(dist, &
               pp_r_grid_min(i_pp_species), pp_r_grid_inc(i_pp_species))
 

           d_partial_core_rho_d_R = val_spline( i_r_log, &
                        partial_core_dens_deriv_spl(1,1,i_pp_species), &
                        n_points_pp_fn(i_pp_species))

           d_partial_core_rho_d_R = d_partial_core_rho_d_R/(4*3.14159265359)
 

        endif


           
           do i_coord=1,3,1

               pp_nlcc_forces(i_coord, i_pp_atom) = pp_nlcc_forces(i_coord, i_pp_atom) +  &
                         point_term*partial_core_rho_grad_in_batch(i_coord,i_point)


           enddo

      enddo !i_pp_atoms

  enddo !i_point

  
end subroutine evaluate_nlcc_forces
!******
