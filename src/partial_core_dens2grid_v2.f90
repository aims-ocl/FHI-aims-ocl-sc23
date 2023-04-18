 subroutine partial_core_dens2grid_v2()

!  PURPOSE
!  Subroutine evaluate_free_atom_sums tabulates the superposition of free-atom
!  densities and hartree potentials on the entire integration grid, for use 
!  in the construction of the Hartree potential.
!
!  USES

   use dimensions
   use runtime_choices
   use grids
   use geometry
   use pbc_lists
   use species_data
   use spline
   use pseudodata
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
!    Release version, FHI-aims (2012).
!  SOURCE
!

    

   !  local variables
   real*8, dimension(3,n_pp_atoms) :: dir_tab

   real*8, dimension(n_pp_atoms) :: dist_tab, dist_tab_sq

   real*8 :: dir_current(3), dir_current_sq 
   real*8, dimension(3) :: coord_current, relvec
   real*8 :: dist, dist_current, dist_current_sq
   real*8 :: atomic_gradient

   real*8 :: ilog

   logical :: point_on_atom 

   integer :: current_atom

   !     counters

   integer :: i_coord, i_atom
   integer :: i_pp_atom, i_pp_species
   integer :: i_full_points, i_my_batch, i_index




   !  begin work


   !       initialize
   if(allocated(partial_core_rho)) then 
      deallocate(partial_core_rho)
   endif 

   if(allocated(partial_core_rho_grad)) then 
      deallocate(partial_core_rho_grad)
   endif 

   allocate(partial_core_rho(n_full_points))
   partial_core_rho = 0d0
!   if (use_density_gradient) then
      allocate(partial_core_rho_grad(3,n_full_points))
      partial_core_rho_grad = 0d0
!   end if

   !       Contributions from all pseudocores

   i_full_points = 0

   do i_my_batch = 1, n_my_batches, 1
     do i_index = 1, batches(i_my_batch)%size, 1 
        i_full_points = i_full_points + 1
        coord_current   = batches(i_my_batch)%points(i_index)%coords(:)
        current_atom    = batches(i_my_batch) % points(i_index) % index_atom


        dir_current(:)  = coord_current(:)- &
                          coords_center(:,centers_basis_integrals(current_atom))
        dist_current_sq = dir_current(1)*dir_current(1) &
                        + dir_current(2)*dir_current(2) &
                        + dir_current(3)*dir_current(3)
        dist_current    = sqrt(dist_current_sq)

        do i_pp_atom = 1, n_pp_atoms,1
            dir_tab(1,i_pp_atom) = coord_current(1) - &
                 pp_coords(1,i_pp_atom)
            dir_tab(2,i_pp_atom) = coord_current(2) - &
                 pp_coords(2,i_pp_atom)
            dir_tab(3,i_pp_atom) = coord_current(3) - &
                 pp_coords(3,i_pp_atom)
        
            dist_tab_sq(i_pp_atom) = &
                   dir_tab(1,i_pp_atom)*dir_tab(1,i_pp_atom) &
                 + dir_tab(2,i_pp_atom)*dir_tab(2,i_pp_atom) &
                 + dir_tab(3,i_pp_atom)*dir_tab(3,i_pp_atom) 

            dist_tab (i_pp_atom) = sqrt(dist_tab_sq(i_pp_atom))

!    Normalization of direction vector
            dir_tab(1,i_pp_atom) = dir_tab(1,i_pp_atom)/dist_tab(i_pp_atom)
            dir_tab(2,i_pp_atom) = dir_tab(2,i_pp_atom)/dist_tab(i_pp_atom)
            dir_tab(3,i_pp_atom) = dir_tab(3,i_pp_atom)/dist_tab(i_pp_atom)


        end do


        point_on_atom = .false.
           if ( dist_current.eq.0.d0) then
             point_on_atom = .true.
           end if

        if (.not. point_on_atom) then

           do i_pp_atom = 1, n_pp_atoms

               i_pp_species = pp_species(i_pp_atom)

               relvec = coord_current(:) - pp_coords(:, i_pp_atom)
               dist = sqrt(sum(relvec**2))


               ilog = invert_log_grid(dist, &
               &                    pp_r_grid_min(i_pp_species), pp_r_grid_inc(i_pp_species))


               partial_core_rho(i_full_points) = partial_core_rho(i_full_points) + &
                            val_spline ( ilog, partial_core_dens_spl(1,1,i_pp_species), &
                            n_points_pp_fn(i_pp_species) )
       
   !       if required, we need the gradient of the partial core density as well
 !              if (use_density_gradient) then
! DB 27.1.13: rausgenommen, damit die nlcc kraft damit funktioniert
                  atomic_gradient = val_spline( ilog, &
                        partial_core_dens_deriv_spl(1,1,i_pp_species), &
                        n_points_pp_fn(i_pp_species))

                  do i_coord = 1,3,1

                     partial_core_rho_grad(i_coord,i_full_points) = &
                          partial_core_rho_grad(i_coord,i_full_points) + &
                          atomic_gradient * dir_tab(i_coord, i_pp_atom)

                  enddo


  !             endif
           enddo 


        end if !point on atom

     enddo

   enddo

!  partitial_tabs for the global grid bring in a factor of 4pi
!  we have to scale here with 1 over 4pi to guarrentee normalization

   partial_core_rho(:) = partial_core_rho(:)/(4*3.14159265359)
   partial_core_rho_grad(:,:) = partial_core_rho_grad(:,:)/(4*3.14159265359)





 end subroutine partial_core_dens2grid_v2
