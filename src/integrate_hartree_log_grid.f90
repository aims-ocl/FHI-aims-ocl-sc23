!****s* FHI-aims/integrate_hartree_log_grid
!  NAME
!    integrate_hartree_log_grid
!  SYNOPSIS

subroutine integrate_hartree_log_grid &
     ( i_atom, angular_part_spl, delta_v_hartree_part_at_zero,  &
     delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
     multipole_radius, l_hartree_max_far_distance, outer_potential_radius )

!  PURPOSE
! 
! This subroutine prepares the partitioned Hartree potential (for the 
! density difference between the free-atom superposition density and the
! actual mixed density in each s.c.f. iteration) as a spline function on the 
! logarithmic integration grid for each atom.  
!
! In short, we receive as input:
!
!   Multicenter multipole density components n^\tilde_{l,m,atom}(r)
!     tabulated as a spline function on the radial shells of the 3D integration
!     grid
!
!   We produce as output the quantities which can be calculated a priori as global 
!   properties of this density:
!
!   - The electrostatic potential components at the nucleus
!   - The multipole moments of the partitioned density outside the radius 
!     where the density becomes zero. (For the l=0 multipole, this would be
!     the overall charge ... etc.)
!   - The maximum angular momentum up to which the far field components 
!     must even be considered
!   - and the actual radii of the multipole density components (i.e., the radius
!     beyond which there is not multipole density for each l,m component)
!
!   The actual partitioned Hartree components are then computed on the fly only
!   for the necessary atoms later in sum_up_whole_potential, one by one. (The required
!   storage can break a calculation, while the computational effort of subroutine
!   integrate_delta_v_hartree below is really not much.)
!
!  USES

     use dimensions
     use runtime_choices
     use grids
     use geometry
     use species_data
     use spline
     use free_atoms
     use analytic_multipole_coefficients
     use constants
     use hartree_potential_storage
     implicit none

!  ARGUMENTS

     real*8,dimension((l_pot_max+1)**2,n_max_spline,n_max_radial+2)::angular_part_spl
     real*8 delta_v_hartree_part_at_zero
     real*8::delta_v_hartree_deriv_l0_at_zero(3)
     real*8, dimension( ( l_pot_max + 1)**2) :: multipole_moments
     real*8 :: multipole_radius
     integer i_atom
     integer:: l_hartree_max_far_distance

    ! VB: Cutoff radius for Hartree potential components of present
     !     atom, as a function of l
     real*8 :: outer_potential_radius ( 0:l_pot_max)


! INPUTS
! o i_atom -- atom number to which the current potential parts belong
! o angular_part_spl -- (electron density - free atoms electron density) = delta rho in multipole expansion form.
!                        This means simply a splined version of the partitioned multipole density
!
!  OUTPUT
! o delta_v_hartree_part_at_zero -- delta Hartree potential at the centers of atoms
! o delta_v_hartree_deriv_l0_at_zero -- derivative of delta Hartree potential at centers of atoms
!                                       for l=0 component
! o multipole_moments -- multipole moments of the delta Hartree potential
! o multipole_radius --  radius of multipole components parts witch need spline
! o l_hartree_max_far_distance -- the maximum l component for periodic Hartree potential far distance part.
!                                 This means the part witch is Fourier transformed, and the part wich calculated analytically at real space
! o outer_potential_radius --  Cutoff radius for Hartree potential components of present
!                              atom, as a function of l
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

  !  angular_integral is the angular integral over
  !  Y_lm(Omega) * partition_fn(at,r,Omega) * delta_rho(r,Omega),
  !  Delley 1990, Eq. (11) 
  !  i.e. angular_integral = angular_integral(r) . For convenience

     real*8, dimension((l_pot_max+1)**2, n_max_grid  ):: angular_integral_log

     integer index_lm(-l_pot_max:l_pot_max, 0:l_pot_max)
     real*8, dimension(0:l_pot_max) :: prefactor

     real*8:: alpha
     real*8:: dr_coef
     integer:: multipole_moments_radius(0:l_pot_max)
     real*8 :: r_test
  
     integer i_grid
     integer i_l
     integer i_m
     integer i_index

!    initialize index_lm

     i_index = 0
     do i_l = 0, l_pot_max, 1
        do i_m = -i_l, i_l
           i_index = i_index + 1
           index_lm(i_m, i_l) = i_index
        enddo
     enddo

     do i_l = 0, l_pot_max, 1
        prefactor(i_l) = pi4 / ( 2.0d0*dble(i_l) + 1.0d0)
     enddo





     ! Coefficient for the integrations
     alpha = log(r_grid_inc(species(i_atom)))


     ! tabulate angular_integral_log from the spline coefficients in angular_part_spl

     call spline_angular_integral_log( angular_part_spl, angular_integral_log, i_atom )

!-----------------------------------------
!
!  Next, we occupy ourselves with analyzing the multipole density components 
!  on the logarithmic grid and derive the multipoles of a classical potential that can be used 
!
!  (i)  outside the (partitioned) charge density associated with each atom
!  (ii) in the periodic case, for a Fourier-transformed version of the long-range
!       potential later on.

     ! 1. determine the radius ("multipole_radius") outside of
     !    which the global classical potential is used

     ! VB: As it is, the multipole_radius should be the outermost radius for which the
     !     integration weights (partition_tab) on the grids of a given atom are not zero.  
     !     As the Stratmann and co. partition table is implemented right now,
     !     this radius has a well defined name: outer_partition_radius.
     !
     !     The reliance on multipole_radius_free all over the Hartree potential is
     !     historical, albeit presently correct. 

     multipole_radius =  multipole_radius_free(species(i_atom))
     
     do i_l = 0, l_hartree(species(i_atom)),1
        do i_m = - i_l, i_l, 1
           do i_grid =  1, n_grid(species(i_atom)), 1
           
              if ( abs(angular_integral_log(index_lm(i_m, i_l),i_grid)) & 
                   .gt. far_distance_hartree_multipole_radius_threshold        ) then                 
                 multipole_radius = max(multipole_radius, r_grid(i_grid,species(i_atom)))
              end if

           end do
        end do
     end do

     ! 2. In the periodic case, this radius gets reduced again. This looks very dodgy - the free-atom
     !    multipole charge density radius should NOT have been exceeded in the first place by the 
     !    partitioned charge density.

     if ( n_periodic > 0 .or. force_new_functional) then
        if(multipole_radius > extra_adding_to_hartree_potential_distance + multipole_radius_free(species(i_atom))) then
           multipole_radius = extra_adding_to_hartree_potential_distance + multipole_radius_free(species(i_atom))
        end if
     end if

     ! 3. Now we determine a safe radius where the multipole_moment associated with each lm component can be
     !    determined. Usually, the partitioned multipole charge density should still have a small finite value here.

     ! Note that this is an integer and counted in units of the logarithmic grid.
     multipole_moments_radius = 0

     do i_l = 0, l_hartree(species(i_atom)),1
        do i_m = - i_l, i_l, 1

           do i_grid =  1, n_grid(species(i_atom)), 1
              if ( abs(angular_integral_log(index_lm(i_m, i_l),i_grid)) &
                   .gt. far_distance_hartree_multipole_moment_radius_threshold ) then
                            
                 multipole_moments_radius(i_l) = max( multipole_moments_radius(i_l), i_grid)

              end if
           end do

        end do
     end do 

     ! 4. Based on multipole_moments_radius, we determine the actual multipole moments
     !    by inward integration.

     do i_l = 0, l_hartree(species(i_atom)),1
        do i_m = - i_l, i_l, 1
           multipole_moments(index_lm(i_m, i_l)) = 0.0d0
        end do
     end do

     do i_l = 0, l_hartree(species(i_atom)),1
        do i_m = - i_l, i_l, 1

           do i_grid =  multipole_moments_radius(i_l), 1, -1
              dr_coef  = r_grid(i_grid,species(i_atom))**2 * alpha &
                * r_grid(i_grid, species(i_atom))

              multipole_moments(index_lm(i_m, i_l)) = &
                multipole_moments(index_lm(i_m, i_l)) &
                + angular_integral_log(index_lm(i_m, i_l),i_grid) * dr_coef &
                *  r_grid(i_grid,species(i_atom))**(i_l)*prefactor(i_l)/4.0d0/pi
           end do

        end do
     end do

     ! VB: This is an awkward place. 
     !
     !     In principle we now know the multipole moments of
     !     the splined density components on the logarithmic grid. 
     !
     !     However, these moments
     !     will not be the same as the original multipole moments of the initial
     !     splined density components on the radial grid, which are in principle
     !     exact, and which are the only components that matter for us.
     !
     !     The splined density components on the radial grid are called angular_integral_spl, 
     !     and they are obtained by subroutine
     !     get_rho_multipole_spl in module hartree_potential_storage.f90 .
     ! 
     !     They are always recreated before use - the only thing that is ever stored is 
     !     rho_multipole itself (on the radial shells of the 3D integration grid, not yet splined). 
     !
     !     So if we want to correct the resulting potential on the splined logarithmic grid, we
     !     must actually add a small auxiliary density on the logarithmic grid only that
     !     compensates the multipole error that arose from splining onto the logarithmic grid.
     
     if (compensate_multipole_errors) then
        ! We first determine the norm of the actual compensating density 
        ! for multipole_radius as the outermost radius. We must do this here and not globally
        ! because theoretically, the multipole_radius can differ for different atoms.

        ! compensation_norm is the multipole moment of the compensating density 
        ! on the logarithmic grid.

        ! For now, only compensate the charge (zero multipole component)
        i_l = 0 

        compensation_norm(i_atom) = 0.d0
        if ( multipole_moments_radius(i_l).gt.0) then
          compensation_radius(i_atom) = r_grid(multipole_moments_radius(i_l),species(i_atom))
        else
          compensation_radius(i_atom) = 0.d0
        end if

        do i_grid =  multipole_moments_radius(i_l), 1, -1
           dr_coef  = r_grid(i_grid,species(i_atom))**2 * alpha &
             * r_grid(i_grid, species(i_atom))

           compensation_norm(i_atom) = &
             compensation_norm(i_atom) &
             + compensating_density( r_grid(i_grid,species(i_atom)), compensation_radius(i_atom), i_l ) &
             * dr_coef &
             *  r_grid(i_grid,species(i_atom))**(i_l)*prefactor(i_l)/4.0d0/pi
        end do
        ! now store the factor that must actually be multiplied to the compensating density
        ! to offset the multipole moment difference (original_multipole_moments - multipole_moments)
        if (compensation_norm(i_atom).gt.0.d0) then

          compensation_norm(i_atom) = & 
          (original_multipole_moments(1,i_atom) - multipole_moments(1)) / compensation_norm(i_atom)

          total_compensated_charge = total_compensated_charge & 
            + (original_multipole_moments(1,i_atom) - multipole_moments(1))

        end if  

        ! now add the compensation to angular_integral_log ...
        do i_grid =  multipole_moments_radius(i_l), 1, -1
           angular_integral_log(1,i_grid) = &
             angular_integral_log(1,i_grid) &
             + compensation_norm(i_atom) &
             * compensating_density( r_grid(i_grid,species(i_atom)), compensation_radius(i_atom), i_l )
        end do

        ! ... and recompute the multipole moments as above
        multipole_moments(1) = 0.d0
        do i_grid =  multipole_moments_radius(i_l), 1, -1
           dr_coef  = r_grid(i_grid,species(i_atom))**2 * alpha &
             * r_grid(i_grid, species(i_atom))

           multipole_moments(1) = &
             multipole_moments(1) &
             + angular_integral_log(1,i_grid) * dr_coef &
             *  r_grid(i_grid,species(i_atom))**(i_l)*prefactor(i_l)/4.0d0/pi
        end do

     end if


     ! 5. Determine a global maximum angular momentum for the analytical far-distance
     !    multipole potential of the present atom
     !    l_hartree_max_far_distance is additionally capped at a configurable value l_hartree_far_distance

     if ( n_periodic == 0 .and. .not. use_hartree_non_periodic_ewald ) then

       l_hartree_max_far_distance = l_hartree(species(i_atom))

     else
       ! VB: The next modification is equivalent to a discontinuous cutoff as a function of
       !     moment and can lead to discontinuities of the total energy as a function of geometry.
       !
       !     The far_distance_hartree_multipole_moment_threshold value must be low enough to make such 
       !     discontinuities numerically absolutely irrelevant.
       !
       !     Paula Havu showed very early that small components could be due to numerical noise,
       !     which is why they are cut off. 
       !
       !     However, to me (VB) it is unclear whether cutting off small multipole moments is even still
       !     necessary. At the very least, a smooth cutoff should be adopted, not a hard step.

       l_hartree_max_far_distance = 0
       do i_l = 0, l_hartree(species(i_atom)),1
          do i_m = - i_l, i_l, 1

            if ( abs(multipole_moments(index_lm(i_m, i_l))/ multipole_radius**(i_l+1)) > &
                far_distance_hartree_multipole_moment_threshold) then
              l_hartree_max_far_distance=  max(l_hartree_max_far_distance, i_l)
            else
              multipole_moments(index_lm(i_m, i_l)) = 0.d0                
            end if

          end do
       end do

     end if

     l_hartree_max_far_distance = min(l_hartree_max_far_distance, l_hartree_far_distance)

     ! 6. VB: After all the above various multipole thresholds, let's actually get to work and determine
     !    an l-dependent threshold radius for different multipole components that we can use to reduce the
     !    formal O(N^2) scaling of the Hartree potential

     if (multipole_threshold.gt.0.d0) then

        do i_l = 0, l_hartree(species(i_atom)), 1

          ! use the radius specified by multipole_moments_radius(i_l) as the inner
          ! threshold

          if (multipole_moments_radius(i_l).gt.0) then
            outer_potential_radius(i_l) = r_grid( multipole_moments_radius(i_l), species(i_atom) )
          else
            outer_potential_radius(i_l) = 0.d0
          end if

          do i_m = - i_l, i_l, 1

           if ( abs(multipole_moments( index_lm(i_m,i_l) ) ) .gt. 0.d0 ) then

              r_test = abs(multipole_moments( index_lm(i_m,i_l) )) / multipole_threshold
              r_test = log (r_test) / dble(i_l+1)
              r_test = exp(r_test)

            else
              r_test = 0.d0
            end if

            outer_potential_radius(i_l) = max( outer_potential_radius(i_l), r_test )

          enddo

        enddo

     else
       ! no threshold requested - set ridiculously high radius.

       outer_potential_radius ( 0:l_pot_max) = 100000.d0

     end if

     ! Square the outer_potential_radii for later use with dist_tab_sq in
     ! sum_up_whole_potential_p1
     do i_l = 0, l_hartree(species(i_atom)), 1
       outer_potential_radius(i_l) = outer_potential_radius(i_l)**2.d0
     enddo


     ! 7. Finally, we calculate only one piece of the Hartree potential on the nucleus
     ! by extrapolation.

     if (n_periodic .gt. 0 .or. (force_new_functional)) then
        ! "Zero limit" of the hartree potential is needed explicitly.
        ! = the potential of the monopole (l=0), at position of the nucleus.
        !   and the potential derivative of the dipole (l=1), at position of the nucleus.

        i_l = 1

        do i_m = -i_l, i_l, 1

           delta_v_hartree_deriv_l0_at_zero(i_m+2) = 0.d0
           do  i_grid = 1,  n_grid(species(i_atom)),1 
           
              dr_coef  = r_grid(i_grid,species(i_atom)) * alpha

              delta_v_hartree_deriv_l0_at_zero(i_m+2) = &
                  delta_v_hartree_deriv_l0_at_zero(i_m+2) &
                   +  angular_integral_log(index_lm(i_m, i_l),i_grid ) *  dr_coef

           end do
           delta_v_hartree_deriv_l0_at_zero(i_m+2) =  delta_v_hartree_deriv_l0_at_zero(i_m+2)*pi4/3.d0

        end do

        delta_v_hartree_part_at_zero = 0.d0
        do  i_grid = 1,  n_grid(species(i_atom)),1 
           
           dr_coef  = r_grid(i_grid,species(i_atom))**2 * alpha

           delta_v_hartree_part_at_zero = &
                delta_v_hartree_part_at_zero &
                +  angular_integral_log(index_lm(0, 0),i_grid ) *  dr_coef

        end do
        delta_v_hartree_part_at_zero =   delta_v_hartree_part_at_zero*pi4

     end if

   end subroutine integrate_hartree_log_grid


subroutine integrate_hartree_log_grid_supercell &
     ( i_atom, angular_part_spl, &  
       delta_v_hartree_part_at_zero, delta_v_hartree_deriv_l0_at_zero ) 

!  PURPOSE
! shanghui add: this is just a simplied version of integrate_hartree_log_grid, to just give
! delta_v_hartree_part_at_zero
! 
! This subroutine prepares the partitioned Hartree potential (for the 
! density difference between the free-atom superposition density and the
! actual mixed density in each s.c.f. iteration) as a spline function on the 
! logarithmic integration grid for each atom.  
!
! In short, we receive as input:
!
!   Multicenter multipole density components n^\tilde_{l,m,atom}(r)
!     tabulated as a spline function on the radial shells of the 3D integration
!     grid
!
!   We produce as output the quantities which can be calculated a priori as global 
!   properties of this density:
!
!   - The electrostatic potential components at the nucleus
!
!
!  USES

     use dimensions
     use runtime_choices
     use grids
     use geometry
     use species_data
     use spline
     use free_atoms
     use constants
     use hartree_potential_storage
     implicit none

!  ARGUMENTS

     integer i_atom
     real*8,dimension((l_pot_max+1)**2,n_max_spline,n_max_radial+2) :: angular_part_spl
     real*8 delta_v_hartree_part_at_zero
     real*8::delta_v_hartree_deriv_l0_at_zero(3)


! INPUTS
! o i_atom -- atom number to which the current potential parts belong
! o angular_part_spl -- (electron density - free atoms electron density) = delta rho in multipole expansion form.
!                        This means simply a splined version of the partitioned multipole density
!
!  OUTPUT
! o delta_v_hartree_part_at_zero -- delta Hartree potential at the centers of atoms
! o delta_v_hartree_deriv_l0_at_zero -- derivative of delta Hartree potential at centers of atoms
!                                       for l=0 component
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

  !  angular_integral is the angular integral over
  !  Y_lm(Omega) * partition_fn(at,r,Omega) * delta_rho(r,Omega),
  !  Delley 1990, Eq. (11) 
  !  i.e. angular_integral = angular_integral(r) . For convenience

     real*8, dimension((l_pot_max+1)**2, n_max_grid  ):: angular_integral_log

     integer index_lm(-l_pot_max:l_pot_max, 0:l_pot_max)

     real*8:: alpha
     real*8:: dr_coef
     real*8 :: r_test
  
     integer i_grid
     integer i_l
     integer i_m
     integer i_index

!    initialize index_lm
 
     i_index = 0
     do i_l = 0, l_pot_max, 1
        do i_m = -i_l, i_l
           i_index = i_index + 1
           index_lm(i_m, i_l) = i_index
        enddo
     enddo



     ! Coefficient for the integrations
     alpha = log(r_grid_inc(species(i_atom)))


     ! tabulate angular_integral_log from the spline coefficients in angular_part_spl
     call spline_angular_integral_log( angular_part_spl, angular_integral_log, i_atom )


     ! 7. Finally, we calculate only one piece of the Hartree potential on the nucleus
     ! by extrapolation.

     if (n_periodic .gt. 0 .or. (force_new_functional)) then
        ! "Zero limit" of the hartree potential is needed explicitly.
        ! = the potential of the monopole (l=0), at position of the nucleus.
        !   and the potential derivative of the dipole (l=1), at position of the nucleus.

        i_l = 1

        do i_m = -i_l, i_l, 1

           delta_v_hartree_deriv_l0_at_zero(i_m+2) = 0.d0
           do  i_grid = 1,  n_grid(species(i_atom)),1 
           
              dr_coef  = r_grid(i_grid,species(i_atom)) * alpha

              delta_v_hartree_deriv_l0_at_zero(i_m+2) = &
                  delta_v_hartree_deriv_l0_at_zero(i_m+2) &
                   +  angular_integral_log(index_lm(i_m, i_l),i_grid ) *  dr_coef

           end do
           delta_v_hartree_deriv_l0_at_zero(i_m+2) =  delta_v_hartree_deriv_l0_at_zero(i_m+2)*pi4/3.d0

        end do

        delta_v_hartree_part_at_zero = 0.d0
        do  i_grid = 1,  n_grid(species(i_atom)),1 
           
           dr_coef  = r_grid(i_grid,species(i_atom))**2 * alpha

           delta_v_hartree_part_at_zero = &
                delta_v_hartree_part_at_zero &
                +  angular_integral_log(index_lm(0, 0),i_grid ) *  dr_coef

        end do
        delta_v_hartree_part_at_zero =   delta_v_hartree_part_at_zero*pi4

     end if

end subroutine integrate_hartree_log_grid_supercell



!******	
!-----------------------------------------------------------------------------------------------------------
subroutine spline_angular_integral_log &
     ( angular_part_spl, angular_integral_log, i_atom )

!  PURPOSE
! 
! The subrotine calculates angular_integral_log from the spline coefficients in angular_part_spl
!
!  USES

     use dimensions
     use runtime_choices
     use grids
     use geometry
     use species_data
     use spline
     use free_atoms
     use analytic_multipole_coefficients
     use constants
     use hartree_potential_storage
     implicit none

!  ARGUMENTS

     real*8,dimension((l_pot_max+1)**2, n_max_spline, n_max_radial+2) :: angular_part_spl
     real*8,dimension((l_pot_max+1)**2, n_max_grid) :: angular_integral_log
     integer i_atom

! INPUTS
! o angular_part_spl -- (electron density - free atoms electron density) = delta rho in multipole expansion form.
!                        This means simply a splined version of the partitioned multipole density
! o i_atom -- atom index where potential parts belong
!
!  OUTPUT
! o angular_integral_log -- angular integral the logarithmic grid
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


  !  angular_integral is the angular integral over
  !  Y_lm(Omega) * partition_fn(at,r,Omega) * delta_rho(r,Omega),
  !  Delley 1990, Eq. (11) 
  !  i.e. angular_integral = angular_integral(r) . For convenience


     real*8:: i_r_radial

     integer l_h_dim
     integer i_grid

     integer :: i_l

     l_h_dim = (l_hartree(species(i_atom))+1)**2

     !  angular_integral is the angular integral over
     !  Y_lm(Omega) * partition_fn(at,r,Omega) * delta_rho(r,Omega),
     !  Delley 1990, Eq. (11) 
     !  i.e. angular_integral = angular_integral(r) . For convenience

     ! Here we spline the atom-projected and lm-decomposed charge density
     ! onto the logarithmic grid for later integration.
     ! However, the Hartree partition table is bounded by the charge density 
     ! of free atoms - outside the free-atom charge density, the partition table goes
     ! to zero, and so must the multipole charge densities here.
     ! 
     ! So:
     ! IF we are inside the radial integration grid shell that is inside the outermost 
     !    atomic radius, we spline as usual.
     ! IF we are outside the radial integration grid shell that is inside the outermost
     !    atomic radius BUT inside the outermost atomic radius, we extrapolate so that
     !    the charge density becomes zero AT the outermost radius at the latest.
     !    This is handled by the tabulated spline itself.
     ! IF we are outside the outermost atomic radius, THEN all charge densities must be 
     !    zero by definition.

     do i_grid = 1, n_grid(species(i_atom))

        if (r_grid(i_grid,species(i_atom)).lt.multipole_radius_free(species(i_atom))) then

          i_r_radial = invert_radial_grid( r_grid(i_grid,species(i_atom)), n_radial(species(i_atom)), &
             scale_radial(species(i_atom)))+1
     
          call spline_vector_v2( i_r_radial,  angular_part_spl, &
             (l_pot_max+1)**2, n_max_spline, n_max_radial+2, n_radial(species(i_atom))+2, &
             l_h_dim, angular_integral_log(1,i_grid))
    
        else

          angular_integral_log(1:l_h_dim,i_grid) = 0.d0

        end if    
     end do

     ! If we are adding an extra compensating charge to zero the charge error (monopole moment error)
     ! on the logarithmic grid as efficiently as possible, we must add that term to angular_integral_log here
     if (compensate_multipole_errors) then
       i_l = 0
       do i_grid =  1, n_grid(species(i_atom)), 1
           angular_integral_log(1,i_grid) = &
             angular_integral_log(1,i_grid) &
             + compensation_norm(i_atom) &
             * compensating_density( r_grid(i_grid,species(i_atom)), compensation_radius(i_atom), i_l )
       end do
     end if

end subroutine spline_angular_integral_log
!******	
!-----------------------------------------------------------------------------------------------------------
subroutine integrate_delta_v_hartree_internal &
     ( angular_integral_log, delta_v_hartree, n_coeff_hartree, i_atom )

!  PURPOSE
! 
! This subroutine calculates the partitioned hartree potential (for the difference density) as
! spline function on the logaritmic integration grid for each atom. 
!
! This is the internal version using angular_integral_log calculated in
! subroutine spline_angular_integral_log
!
!  USES

     use dimensions
     use runtime_choices
     use grids
     use geometry
     use species_data
     use spline
     use free_atoms
     use analytic_multipole_coefficients
     use constants
     implicit none

!  ARGUMENTS

     integer n_coeff_hartree
     real*8,dimension((l_pot_max+1)**2, n_max_grid) :: angular_integral_log
     real*8,dimension((l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid) :: delta_v_hartree
     integer i_atom

! INPUTS
! o angular_integral_log -- angular integral the logarithmic grid
!                           (the result of spline_angular_integral_log)
! o n_coeff_hartree -- number of spline coefficients for delta_v_hartree
!                      This may be 2 or 4 depending if only the value and derivative should be calculated
!                      or if all 4 spline coefficients are needed.
! o i_atom -- atom index where potential parts belong
!
!  OUTPUT
! o delta_v_hartree -- Hartree potential of delta rho = delta Hartree
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

     real*8, dimension((l_pot_max+1)**2) :: integral_zero_r
     real*8, dimension((l_pot_max+1)**2) :: integral_r_infty

     real*8, dimension(0:l_pot_max) :: prefactor

     real*8, parameter :: d_1_24 = 1.d0/24.d0

     real*8:: r_l
     real*8:: r_neg_l1
     real*8:: r_inv
     real*8:: alpha
     real*8:: dr_coef

     real*8:: r_l_AM(0:3)
     real*8:: r_neg_l1_AM(0:3)
     real*8:: r_inv_AM(0:3)
     real*8:: dr_coef_AM(0:3)

     integer l_h_dim
     integer i_grid
     integer i_index
     integer i_l
     integer i_m
     integer i_g


     do i_l = 0, l_pot_max, 1
        prefactor(i_l) = pi4 / ( 2.0d0*dble(i_l) + 1.0d0)
     enddo


     l_h_dim = (l_hartree(species(i_atom))+1)**2

     ! Coefficient for the integrations
     alpha = log(r_grid_inc(species(i_atom)))

     ! Now integrate the 1D potential components by using the Green-function solution
     ! from classical electrostatics ...
     integral_zero_r  = 0.0d0
     integral_r_infty = 0.0d0

     if(Adams_Moulton_integrator)then

        ! The Adams-Moulton linear multistep integrator, using up to 4 terms. 
        !    see http://en.wikipedia.org/wiki/Linear_multistep_method or Abramowitz/Stegun p. 896 for details. 

        ! as in integrate_hartree_log_grid, do so in two parts
        ! part (1) contains the integral from 0 up to the grid radius r
        ! part (2) is the other way around, i.e. from r to "infinity"  

        ! first integral: from zero to radius r, for all possible r

        ! first three terms: do one-by-one with increasing order of integration method. 
        i_index = 0
        do i_l = 0,  l_hartree(species(i_atom)), 1
           do i_m = - i_l, i_l, 1
              i_index = i_index+1


             ! TERM 1 : Integral_1 = h*f_1; but h = 1
              integral_zero_r(i_index) =   integral_zero_r(i_index) &
                   + angular_integral_log(i_index,1 ) *  alpha * r_grid(1,species(i_atom))**(i_l+3)

              delta_v_hartree(i_index, 1, 1) =  &
                   integral_zero_r(i_index ) / r_grid(1,species(i_atom))**(1+i_l)



             ! TERM 2 : Integral_2 = Integral_1 + h(f_2+f_1)/2 
              integral_zero_r(i_index) =   integral_zero_r(i_index) &
                   + (angular_integral_log(i_index,1 ) *  alpha * r_grid(1,species(i_atom))**(i_l+3) &
                   +  angular_integral_log(i_index,2 ) *  alpha * r_grid(2,species(i_atom))**(i_l+3)) * 0.5d0

              delta_v_hartree(i_index, 1, 2) =  &
                   integral_zero_r(i_index ) / r_grid(2,species(i_atom))**(1+i_l)


             ! TERM 3 : Integral_3 = Integral_2 + h(5f_3 + 8f_2 - f_1)/12

              integral_zero_r(i_index) =   integral_zero_r(i_index) &
                   + ( 5 * angular_integral_log(i_index, 3 ) *  alpha * r_grid(3,species(i_atom))**(i_l+3) &
                   +  8 * angular_integral_log(i_index,2 ) *  alpha * r_grid(2,species(i_atom))**(i_l+3) &
                   -  angular_integral_log(i_index,1 ) *  alpha * r_grid(1,species(i_atom))**(i_l+3))/12.d0


              delta_v_hartree(i_index, 1, 3) =  &
                   integral_zero_r(i_index ) / r_grid(3,species(i_atom))**(1+i_l)

           end do
        end do


    
      ! TERM i_grid > 4 : Integral_i = Integral_(i-1) + h[9 f_i + 19 f_(i-1) - 5 f_(i-2) + f_(i-3)]/24
        do i_grid = 4, n_grid(species(i_atom)), 1

           do i_g = 0, 3, 1           
              r_inv_AM(i_g)    = 1.d0/r_grid(i_grid-i_g,species(i_atom))
              r_l_AM (i_g)     = r_inv_AM(i_g)
              r_neg_l1_AM(i_g) = 1.d0
              dr_coef_AM(i_g)  = r_grid(i_grid-i_g,species(i_atom))**2 * alpha * r_grid(i_grid-i_g,species(i_atom))
           end do


              i_index = 0
              do i_l = 0, l_hartree(species(i_atom)), 1
      
                 do i_g = 0, 3, 1
                    r_l_AM(i_g) = r_l_AM(i_g) * r_grid(i_grid-i_g,species(i_atom))
                 end do
                 r_neg_l1_AM(0) = r_neg_l1_AM(0) * r_inv_AM(0)

                 do i_m = - i_l, i_l, 1
                    i_index = i_index+1

                    integral_zero_r(i_index) =   integral_zero_r(i_index) &
                         + ( 9 * angular_integral_log(i_index,i_grid   ) *  r_l_AM(0) * dr_coef_AM(0)  &
                         +  19 * angular_integral_log(i_index,i_grid-1 ) *  r_l_AM(1) * dr_coef_AM(1)  &
                         -  5  * angular_integral_log(i_index,i_grid-2 ) *  r_l_AM(2) * dr_coef_AM(2)  &
                         +       angular_integral_log(i_index,i_grid-3 ) *  r_l_AM(3) * dr_coef_AM(3)) * d_1_24

                    delta_v_hartree(i_index, 1, i_grid) =  integral_zero_r(i_index ) * r_neg_l1_AM(0)

                                 
                 end do
              end do
           end do        ! end radial loop for calculation of first integral
        


        
        ! start integrating from the outside in, again using the Adams-Moulton linear multistep integrator.
        ! the first three terms warrant special treatment, similar to the above integral. 
        integral_r_infty = 0d0
        i_index = 0
        do i_l = 0,  l_hartree(species(i_atom)), 1
           do i_m = - i_l, i_l, 1
              i_index = i_index+1


              ! TERM 1 : Integral_N = h*f_N; but h = 1            

              delta_v_hartree(i_index, 1, n_grid(species(i_atom))) =  delta_v_hartree(i_index, 1, n_grid(species(i_atom))) &
                   + integral_r_infty(i_index) * r_grid(n_grid(species(i_atom)),species(i_atom))**(i_l)

              integral_r_infty(i_index) = integral_r_infty(i_index) &
                   + angular_integral_log(i_index,n_grid(species(i_atom))) / & 
                     r_grid(n_grid(species(i_atom)),species(i_atom))**(i_l+1) &
                   * r_grid(n_grid(species(i_atom)),species(i_atom))**2 * alpha  * r_grid(n_grid(species(i_atom)), species(i_atom))

              delta_v_hartree(i_index, 1, n_grid(species(i_atom))) =  delta_v_hartree(i_index, 1, n_grid(species(i_atom))) &
                   * prefactor(i_l)






              ! TERM 2 : Integral_(N-1) = Integral_N + h(f_(N-1)+f_N)/2 

                 delta_v_hartree(i_index, 1, n_grid(species(i_atom))-1) =  delta_v_hartree(i_index, 1, n_grid(species(i_atom))-1)&
                      + integral_r_infty(i_index) * r_grid(n_grid(species(i_atom))-1,species(i_atom))**(i_l)


                 integral_r_infty(i_index) = integral_r_infty(i_index) &
!
                      +  (angular_integral_log(i_index,n_grid(species(i_atom))) / & 
                        r_grid(n_grid(species(i_atom)),species(i_atom))**(i_l+1) &
                      * r_grid(n_grid(species(i_atom)),species(i_atom))**2 * alpha  & 
                      * r_grid(n_grid(species(i_atom)), species(i_atom)) &
!
                      +  angular_integral_log(i_index,n_grid(species(i_atom))-1) / & 
                         r_grid(n_grid(species(i_atom))-1,species(i_atom))**(i_l+1) &
                      * r_grid(n_grid(species(i_atom))-1,species(i_atom))**2 * alpha  & 
                      * r_grid(n_grid(species(i_atom))-1, species(i_atom)))*0.5d0

                 delta_v_hartree(i_index, 1, n_grid(species(i_atom))-1) =  delta_v_hartree(i_index, 1, n_grid(species(i_atom))-1)&
                      * prefactor(i_l)


                 ! TERM 3 : Integral_(N-2) = Integral_(N-1) + h(5f_(N-2) + 8f_(N-1) - f_N)/12

                 delta_v_hartree(i_index, 1, n_grid(species(i_atom))-2) =  delta_v_hartree(i_index, 1, n_grid(species(i_atom))-2) &
                      + integral_r_infty(i_index) * r_grid(n_grid(species(i_atom))-2,species(i_atom))**(i_l)

                 integral_r_infty(i_index) = integral_r_infty(i_index) &
                      !
                      + ( -1 * angular_integral_log(i_index,n_grid(species(i_atom))) / & 
                        r_grid(n_grid(species(i_atom)),species(i_atom))**(i_l+1) &
                      * r_grid(n_grid(species(i_atom)),species(i_atom))**2 * alpha  & 
                      * r_grid(n_grid(species(i_atom)), species(i_atom)) &
                      !
                      + 8 * angular_integral_log(i_index,n_grid(species(i_atom))-1) / & 
                        r_grid(n_grid(species(i_atom))-1,species(i_atom))**(i_l+1) &
                      * r_grid(n_grid(species(i_atom))-1,species(i_atom))**2 * alpha  & 
                      * r_grid(n_grid(species(i_atom))-1, species(i_atom)) &
                      !
                      + 5 * angular_integral_log(i_index,n_grid(species(i_atom))-2) / & 
                        r_grid(n_grid(species(i_atom))-2,species(i_atom))**(i_l+1) &
                      * r_grid(n_grid(species(i_atom))-2,species(i_atom))**2 * alpha  & 
                      * r_grid(n_grid(species(i_atom))-2, species(i_atom)))/12.d0


                 delta_v_hartree(i_index, 1, n_grid(species(i_atom))-2) =  delta_v_hartree(i_index, 1, n_grid(species(i_atom))-2)&
                      * prefactor(i_l)


          end do
        end do


        ! all remaining terms
        ! Integral_i = Integral_(i+1) + h[9 f_i + 19 f_(i+1) - 5 f_(i+2) + f_(i+3)]/24
        do i_grid = n_grid(species(i_atom))-3, 1, -1

           do i_g = 0, 3, 1           
              r_inv_AM(i_g)    = 1.d0/r_grid(i_grid+i_g,species(i_atom))
              r_l_AM (i_g)     = r_inv_AM(i_g)
              r_neg_l1_AM(i_g) = 1.d0
              dr_coef_AM(i_g)  = r_grid(i_grid+i_g,species(i_atom))**2 * alpha * r_grid(i_grid+i_g,species(i_atom))
           end do

           i_index = 0 
           do i_l = 0,  l_hartree(species(i_atom)), 1
              
              do i_g = 0, 3, 1
                 r_neg_l1_AM(i_g) = r_neg_l1_AM(i_g) * r_inv_AM(i_g)
              end do
              r_l_AM(0) = r_l_AM(0) * r_grid(i_grid, species(i_atom))
              
              do i_m = - i_l, i_l, 1
                 i_index = i_index+1

                 integral_r_infty(i_index) = integral_r_infty(i_index) &
                      !
                      + ( 9 * angular_integral_log(i_index,i_grid  ) * r_neg_l1_AM(0) *  dr_coef_AM(0) &
                      +  19 * angular_integral_log(i_index,i_grid+1) * r_neg_l1_AM(1) *  dr_coef_AM(1) &
                      -   5 * angular_integral_log(i_index,i_grid+2) * r_neg_l1_AM(2) *  dr_coef_AM(2) &
                      +       angular_integral_log(i_index,i_grid+3) * r_neg_l1_AM(3) *  dr_coef_AM(3)) * d_1_24


                 delta_v_hartree(i_index, 1, i_grid) =  delta_v_hartree(i_index, 1, i_grid) &
                      + integral_r_infty(i_index) * r_l_AM(0)

                 delta_v_hartree(i_index, 1, i_grid) =  delta_v_hartree(i_index, 1, i_grid) * prefactor(i_l)

              end do
           end do
        end do ! end calculation of second integral 


       

     else !----------------------------------------------------------------------------


        ! Now to the integrations
        ! First part of the integral 0 -> r

        do i_grid = 1, n_grid(species(i_atom)), 1
           
           r_inv    = 1.d0/r_grid(i_grid,species(i_atom))
           r_l      = r_inv
           r_neg_l1 = 1.d0
           
           !       radial integration weight on the logarrithmic grid alpha*r times usual radial 
           !       integration weight from integral r^2 dr
           dr_coef  = r_grid(i_grid,species(i_atom))**2 * alpha &
                * r_grid(i_grid,species(i_atom))
           
           i_index = 0 
           do i_l = 0, l_hartree(species(i_atom)), 1
              
              
              r_l = r_l * r_grid(i_grid,species(i_atom))
              r_neg_l1 = r_neg_l1 * r_inv


              do i_m = - i_l, i_l, 1
                 i_index = i_index+1
               

                 integral_zero_r(i_index) =   integral_zero_r(i_index) &
                      + angular_integral_log(i_index,i_grid ) * r_l * dr_coef                 

                 delta_v_hartree(i_index, 1, i_grid) =  integral_zero_r(i_index ) * r_neg_l1


              end do
           end do
        end do
        

        !!  run a second time through the radial grid from outward to inward
        !!  (not the whole integration grid!)
        !!  to evaluate integral_r_infty via tabulated angular_integral


        do i_grid =  n_grid(species(i_atom)), 1, -1
           
           r_inv    = 1.d0/r_grid(i_grid,species(i_atom))
           r_l      = r_inv
           r_neg_l1 = 1.d0
           dr_coef  = r_grid(i_grid,species(i_atom))**2 * alpha &
                * r_grid(i_grid, species(i_atom))
           
           
           i_index = 0 
           do i_l = 0, l_hartree(species(i_atom)), 1
              
              r_l       = r_l * r_grid(i_grid,species(i_atom))
              r_neg_l1  = r_neg_l1 * r_inv
              
              do i_m = - i_l, i_l, 1
              i_index = i_index+1


              delta_v_hartree(i_index, 1, i_grid) =  delta_v_hartree(i_index, 1, i_grid) &
                   + integral_r_infty(i_index) * r_l

              integral_r_infty(i_index) = integral_r_infty(i_index) &
                   + angular_integral_log(i_index,i_grid) * r_neg_l1 * dr_coef

              delta_v_hartree(i_index, 1, i_grid) =  delta_v_hartree(i_index, 1, i_grid) &
                   * prefactor(i_l)

              end do
           end do
        end do

     end if !(Adams_Moulton_integrator) and else

     ! Calculate spline coefficients

     call cubic_spline_v2 ( delta_v_hartree, (l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid, &
                            n_grid(species(i_atom)), l_h_dim )

end subroutine integrate_delta_v_hartree_internal
!******	
!-----------------------------------------------------------------------------------------------------------
subroutine integrate_delta_v_hartree &
     ( angular_part_spl, delta_v_hartree, n_coeff_hartree, i_atom )

!  PURPOSE
! 
! The subrotine calculates the partitioned hartree potential (for the difference density) as
! spline function on the logaritmic integration grid for each atom. 
!
! This is just a driver routine calling spline_angular_integral_log (for tabulating angular_integral_log)
! and then integrate_delta_v_hartree_internal (for integrating delta_v_hartree)
!  USES

     use dimensions
     use runtime_choices
     use grids
     use geometry
     use species_data
     use spline
     use free_atoms
     use analytic_multipole_coefficients
     use constants
     implicit none

!  ARGUMENTS

     integer n_coeff_hartree
     real*8,dimension((l_pot_max+1)**2, n_max_spline, n_max_radial+2) :: angular_part_spl
     real*8,dimension((l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid) :: delta_v_hartree
     integer i_atom

! INPUTS
! o angular_part_spl -- (electron density - free atoms electron density) = delta rho in multipole expansion form.
!                        This means simply a splined version of the partitioned multipole density
! o n_coeff_hartree -- number of spline coefficients for delta_v_hartree
!                      This may be 2 or 4 depending if only the value and derivative should be calculated
!                      or if all spline coefficients are needed.
! o i_atom -- atom index where potential parts belong
!
!  OUTPUT
! o delta_v_hartree -- (Hartree potential of delta rho = delta Hartree
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


     real*8, dimension((l_pot_max+1)**2, n_max_grid) :: angular_integral_log

     ! tabulate angular_integral_log from the spline coefficients in angular_part_spl

     call spline_angular_integral_log( angular_part_spl, angular_integral_log, i_atom )

     ! integrate delta_v_hartree

     call integrate_delta_v_hartree_internal( angular_integral_log, delta_v_hartree, n_coeff_hartree, i_atom )

end subroutine integrate_delta_v_hartree
!******	
!-----------------------------------------------------------------------------------------------------------
subroutine integrate_average_atom_potential &
     ( delta_v_hartree, n_coeff_hartree, multipole_radius_sq, adap_outer_radius_sq, multipole_moment, &
       atom_average_es_pot, i_atom )

!  PURPOSE
! 
!    This subroutine takes the partitioned hartree potential (for the difference density) 
!    and computes the spatial average of the zero component of the real-space part
!    on the logarithmic grid.
!
!    This quantity only makes sense in periodic systems for now.
!    We also need to add the error function potential that compensates the long-range tail
!    of the Coulomb potential in Ewald's method.
!
!  USES

     use dimensions
     use runtime_choices
     use grids
     use constants
     use geometry
     use hartree_f_p_functions
     implicit none

!  ARGUMENTS

     integer n_coeff_hartree
     real*8,dimension((l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid) :: delta_v_hartree
     real*8 :: multipole_radius_sq,adap_outer_radius_sq,multipole_moment
     real*8  :: atom_average_es_pot
     integer ::i_atom

! INPUTS
! o delta_v_hartree -- (Hartree potential of delta rho = delta Hartree
! o n_coeff_hartree -- number of spline coefficients for delta_v_hartree
!                      This may be 2 or 4 depending if only the value and derivative should be calculated
!                      or if all spline coefficients are needed.
! o i_atom -- atom index where potential parts belong
!
!  OUTPUT
!
! o atom_average_es_pot -- integral of the atomic electrostatic potential over space
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals",
!    Computer Physics Communications (2009), [...].
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2013).
!  SOURCE

     integer :: i_grid
     real*8 :: alpha, dr_coef
     real*8, dimension(0:0) :: Fp

     ! Coefficient for the integrations
     alpha = log(r_grid_inc(species(i_atom)))

     ! The average that we are looking for is the actual l=0 multipole component
     ! of the atom PLUS the compensating Ewald potential that localizes the
     ! real-space component. To that end, we follow exactly the spatial division used
     ! also in sum_up_whole_potential.

     atom_average_es_pot = 0.d0
     do i_grid = 1, n_grid(species(i_atom)), 1
        dr_coef  = r_grid(i_grid,species(i_atom))**3 * alpha ! radial integration weight on logarithmic grid
        if ( (r_grid(i_grid,species(i_atom))**2) .lt. multipole_radius_sq ) then
           call F_erf(Fp,r_grid(i_grid,species(i_atom)),0)
           atom_average_es_pot = atom_average_es_pot + &
                      dr_coef * ( delta_v_hartree(1,1,i_grid) - multipole_moment * Fp(0) * pi4 )
        else if ( (r_grid(i_grid,species(i_atom))**2) .lt. adap_outer_radius_sq )then
           call F_erfc(Fp,r_grid(i_grid,species(i_atom)),0)
           atom_average_es_pot = atom_average_es_pot + &
                      dr_coef * multipole_moment * Fp(0) * pi4
        end if
     enddo
     atom_average_es_pot = atom_average_es_pot * sqrt(pi4)

   end subroutine integrate_average_atom_potential
!******	
