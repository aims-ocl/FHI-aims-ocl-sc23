subroutine get_rho_mp(x1, rho_mp,rho_gradient_mp)

!!$ Calculates total rho and gradient from multipoles. 
!!$ Uses streamlined routines tab_singl_gradient_ylm_p2_forvdw.f90 
!!$ and increment_ylm_deriv_forvdw.f90.
!!$ 30.11.2010
!!$ Simiam Ghan (SAG), Department of Applied Physics, COMP, Aalto University, Espoo, Finland.
!!$ sag@fyslab.hut.fi

  use constants,                 only: pi4_inv
  use dimensions
  use runtime_choices 
  use grids
  use geometry
  use free_atoms
  use spline
  use localorb_io
  use pbc_lists
  use physics,                   only : multipole_radius_sq
  use hartree_potential_storage, only : get_rho_multipole_spl
  use species_data,              only : l_hartree, multipole_radius_free_sq
  implicit none


  real*8, dimension(1) :: rho_mp 
  real*8, dimension(3) :: rho_gradient_mp 
  real*8 x2(3), x1(3)  
  integer count   !,info
  real*8, dimension(1) :: rho_multipole
  real*8 dist_tab_sq
  real*8 dist_tab_in
  real*8 dist_tab_out
  real*8 dir_tab(3)
  real*8 dir_tab_in(3)
  real*8 dir_tab_out(3)
  real*8 log_weight
  real*8 radial_weight
  real*8 trigonom_tab(4)
  real*8 i_r
  real*8 i_r_log
  real*8 prefactor
  real*8 ylm_tab((l_pot_max+1)**2)
  real*8, dimension((l_pot_max+1)**2) :: delta_v_hartree_multipole_component
  real*8, dimension((l_pot_max+1)**2) :: rho_multipole_component
  integer l_h_dim
  integer :: n_coeff_hartree
  real*8, dimension( (l_pot_max+1)**2 ) :: delta_v_hartree_multipole_deriv
  real*8, dimension((l_pot_max+1)**2  ) :: rho_multipole_deriv
  real*8, dimension( (l_pot_max+1)**2, n_centers_hartree_potential   ) :: dylm_dtheta_tab
  real*8, dimension((l_pot_max+1)**2, n_centers_hartree_potential ) :: scaled_dylm_dphi_tab
  real*8 :: v_hartree_gradient_temp(3)
  real*8 :: rho_multipole_gradient_temp(3)
  real*8 :: rho_multipole_gradient(3)
  real*8 :: d_v_hartree_free_d_r
  real*8 :: d_rho_free_d_r
  real*8 rho_multipole_aux
  integer :: current_spl_atom
  integer :: current_center
  integer :: prev_atom
  integer :: l_atom_max
  integer i_center
  real*8, external :: ddot
  real*8 :: free_hartree_superpos 
  real*8 :: free_rho_superpos 
  real*8, dimension(n_spin) :: rho 
  real*8, dimension(3) :: free_rho_gradient_superpos 
  real*8, dimension(3) :: coord_current
  real*8 :: dist_tab_sq_free(n_centers_basis_integrals)
  real*8 :: dist_tab(n_centers_basis_integrals)
  real*8 :: dir_tab_free(3,n_centers_basis_integrals)
  real*8 :: dir_tab_norm(3,n_centers_basis_integrals)
  real*8 :: i_r_free(n_centers_basis_integrals)
  real*8 :: local_rho_gradient(3,n_spin)
  integer :: n_compute_atoms,i
  integer, dimension(n_centers_basis_integrals) :: center_index
  real*8, dimension(n_centers_basis_integrals) :: temp_free_rho
  integer :: n_compute_occ_atoms
  integer, dimension(n_centers_basis_integrals) :: i_occ2i_compute
  real*8, dimension(n_centers_basis_integrals) :: temp_occ_rho
  real*8 :: aux_spl
  integer :: i_atom, i_coord, i_point, i_spin !  &
  integer :: i_center_L, i_center_free, i_center_L2
  real*8, dimension(3) :: dir_current
  real*8 :: dist_current, dist_current_sq, dens_current, r_temp
  logical getdelta, getfree, freeonly
  integer :: theatom
  logical bump, moveleft
  real*8 limit
  real*8, dimension(3) :: ddd
  real*8 :: current_rho_multipole_spl((l_pot_max+1)**2, n_max_spline, n_max_radial+2)
  !SR for laplacian of free electron density
!   real*8 :: dummy
  
  !begin work
  
  x2(:) = x1(:)   
  
 
  !Apply constant density if too close to the nucleus.  Remember to zero gradients at the end!!! 
  bump = .false.      
  
  do i_atom = 1, n_atoms
     ddd(:) = abs(x2(:)-coords(:,i_atom))
     limit = r_grid_min(species(i_atom))
     
     if( ( ddd(1).le.limit) .and. (ddd(2).le.limit) .and. (ddd(3).le.limit)  ) then
        bump = .true.
        theatom = i_atom  
     endif
  enddo    !assumes cores are seperated.
  
  if(bump)then
     x2(:) = coords(:,theatom) 
     x2(1) = x2(1) + r_grid_min(species(theatom))!assumes density is spherically symmetric this close to atom.
  endif
  
  
  rho_multipole = 0.d0
  rho_multipole_gradient(:) = 0.d0
  rho_mp = 0.d0
  rho_gradient_mp(:) = 0.d0
  count = 0
  
  
!you could optionally solve for only delta or only free rho. 
  getdelta = .true.
  getfree = .true.

  if(getdelta)then
     
        do i_center = 1,  n_centers_hartree_potential, 1
           
           count = count +1

           current_center   = centers_hartree_potential(i_center)
           current_spl_atom = center_to_atom(current_center)
           rho_multipole_gradient_temp(:) = 0.d0

           call tab_single_atom_centered_coords_p0 &
                ( current_center, &
                x2,  &                           !  Coord??
                dist_tab_sq,  &
                dir_tab )
           
           l_atom_max = l_hartree(species(current_spl_atom))
           
           if (dist_tab_sq.lt.multipole_radius_sq(current_spl_atom) ) then
              
              call tab_single_atom_centered_coords_radial_log_p0 &
                   ( current_center, dist_tab_sq, dir_tab,  &
                   dist_tab_in, i_r, i_r_log, dir_tab_in )
              
                 
                 call tab_single_radial_weights_v2 &
                      ( current_spl_atom, dist_tab_in, i_r, &
                      log_weight, radial_weight )
              
              call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)
              
                 call tab_single_gradient_ylm_p2_forvdw &
                      ( trigonom_tab, l_atom_max, l_pot_max,  &
                      ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab)
            
              l_h_dim = (l_atom_max + 1)**2
              
              call get_rho_multipole_spl(current_rho_multipole_spl, current_spl_atom)
              call spline_vector_v2 &
                   ( i_r+1,  current_rho_multipole_spl, &
                   (l_pot_max+1)**2, n_max_spline, n_max_radial+2,  &
                   n_radial(species_center(current_center))+2,  &
                   l_h_dim, &
                   rho_multipole_component)
              
                 
                 rho_multipole_aux = &
                      ddot ( l_h_dim, rho_multipole_component, 1, ylm_tab, 1)
                 
                 rho_multipole = rho_multipole + rho_multipole_aux
                 
                 
                 call spline_deriv_vector_v2 &
                      ( i_r+1,  &
                      current_rho_multipole_spl,&                           
                      (l_pot_max+1)**2, n_max_spline, n_max_radial+2,  &
                      n_radial(species_center(current_center))+2,  &
                      l_h_dim, &
                      rho_multipole_deriv)
                 
                
                 rho_multipole_deriv(1:l_h_dim) = rho_multipole_deriv(1:l_h_dim) * radial_weight
                 
                 do i = 1,l_h_dim
                    
                    delta_v_hartree_multipole_deriv(i) = 0.d0  !&     !DUMMY
                 enddo
                 
                 d_v_hartree_free_d_r =  0.d0  !&   !DUMMY
                 
                 d_rho_free_d_r = pi4_inv *  &
                      val_spline_deriv(i_r_log,  &
                      renormalized_free_rho_spl(1,1,species(current_spl_atom)),  &
                      n_grid(species(current_spl_atom)))  * log_weight
                 
                 call evaluate_v_hartree_and_rho_multipole_gradient &
                      ( dist_tab_in, dir_tab_in,  &
                      trigonom_tab,  &
                      ylm_tab,  &
                      dylm_dtheta_tab,  &
                      scaled_dylm_dphi_tab,  &
                      delta_v_hartree_multipole_component,  &
                      delta_v_hartree_multipole_deriv, &
                      rho_multipole_component, &
                      rho_multipole_deriv,  &
                      d_v_hartree_free_d_r, d_rho_free_d_r, &
                      l_h_dim, &
                      v_hartree_gradient_temp, &
                      rho_multipole_gradient_temp)
               
                 do i = 1,3        !This contains delta rho AND free rho contributions! SAG
                    rho_multipole_gradient(i) = rho_multipole_gradient(i) + rho_multipole_gradient_temp(i)
                 enddo
                 
                 
              endif
           enddo   !i_centers
     
     
     rho_mp = rho_multipole        !This is delta density, mp approximation. !
     
     
     do i = 1,3
        rho_gradient_mp(i) = rho_multipole_gradient(i)    
     enddo
     
     
  endif


  !End multipole rho part. 
  
  if(getfree)then
     !Begin free atom part (as in initialize_grid_storage_p1 )
     
     coord_current(:) = x2(:)! batches(i_my_batch) % points(i_index) % coords(:)
     
     call tab_atom_centered_coords_p0 &
          ( coord_current, &
          dist_tab_sq_free, &
          dir_tab_free, &
          n_centers_basis_integrals, centers_basis_integrals )
     
     n_compute_atoms = 0
     n_compute_occ_atoms = 0
     temp_occ_rho(:) = 0.d0

     do i_center_L = 1, n_centers_basis_integrals, 1         !Worry about later. 
        i_center_free = centers_basis_integrals(i_center_L)
           
        !    if ((i_center.eq.current_atom).or. &                                   !there is no current atom now.
        if((dist_tab_sq_free(i_center_L).lt.multipole_radius_free_sq(species_center(i_center_free)))) then 
           n_compute_atoms                  = n_compute_atoms + 1
           center_index(n_compute_atoms)    = i_center_free
           dist_tab_sq_free(n_compute_atoms)     = dist_tab_sq_free(i_center_L) 
           dir_tab_free(:,n_compute_atoms)       = dir_tab_free(:,i_center_L) 
           ! This ensures consistent handling later, note n_compute_atoms <= i_center_L
           !atom_atom_index(n_compute_atoms) = i_center_L
           ! indexing for later use of atom_atom_tab, which is NOT recomputed here for speed reasons
        end if
        
     enddo

     call tab_global_geometry_p0 &
          ( dist_tab_sq_free,         &
          dir_tab_free,             &  
          dist_tab,            &
          i_r_free,                 &
          dir_tab_norm,        &
          n_compute_atoms,     & 
          center_index )
     
     do i_center_L2 = 1, n_compute_atoms

        i_center_free = center_index(i_center_L2)
        if (.not.empty(center_to_atom(i_center_free))) then
           aux_spl = val_spline &
                ( i_r_free(i_center_L2), renormalized_free_rho_spl(1,1,species_center(i_center_free)), &
                n_grid(species_center(i_center_free)) )  
           n_compute_occ_atoms = n_compute_occ_atoms + 1
           i_occ2i_compute(n_compute_occ_atoms) = i_center_L2
           temp_occ_rho(n_compute_occ_atoms) = aux_spl
        end if
        
     end do
     
     call evaluate_free_atom_sums_p2  &
          ( i_r_free, dir_tab_norm,  &
          free_hartree_superpos,  &
          free_rho_superpos,  &
          free_rho_gradient_superpos, &
          rho(1),  &
          local_rho_gradient,&
          n_compute_atoms, center_index, &
          n_compute_occ_atoms, i_occ2i_compute, temp_occ_rho) !SR: to get laplacian of free atom density add this args: ,dist_tab_sq_free,dummy )
     
     
     !free atom contributions to gradient was taken care of already above. 
     rho_mp = rho_mp + (pi4_inv * free_rho_superpos)    !Total density (multipole approximation) IN ATOMIC UNITS. 
     
     
     
  endif


if(bump)then  !cap rho near atom center
rho_gradient_mp(:) = 0d0
endif



end subroutine get_rho_mp

