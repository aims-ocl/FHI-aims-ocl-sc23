!****s* FHI-aims/get_confined_basis_fns
!  NAME
!   get_confined_basis_fns
!  SYNOPSIS
 subroutine get_confined_basis_fns ( i_species )
! PURPOSE
!  Subroutine get_confined_basis_fns integrates confined Schroedinger equation
!  for one species, to provide radial functions for a given effective potential
! USES
   use constants,       only : bohr, hartree, light_speed
   use dimensions,      only : n_max_grid, use_basis_gradients
   use runtime_choices, only : atomic_solver, atomic_solver_atom_sphere, &
                               flag_KH_core_states, flag_rel, REL_atomic_zora, &
                               REL_kolning_harmon, REL_own, REL_zora, REL_x2c, &
                               REL_4c_dks, wave_threshold
   use grids,           only : n_grid, r_grid, r_grid_inc
   use species_data,    only : conf_cutoff, confined_pot, confined_kinetic, &
                               confined_wave, confined_wave_large, confined_wave_small, &
                               confined_wave_deriv, confined_outer_radius, species_z, &
                               conf_rad, conf_l, conf_n, conf_kappa, core_n_max, cutoff_type, &
                               w_cutoff, confined_eigenval, basis_dep_cutoff_thresh, &
                               scale_cutoff, r_cutoff, n_conf
   use free_atoms,      only : free_potential
   use mpi_tasks,       only : aims_stop
   use localorb_io,     only : localorb_info, use_unit
   use reigen_dftatom
   implicit none
!  ARGUMENTS
   integer i_species
!  INPUT
!  o i_species -- species index
!  OUTPUT
!  none
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
!  i_mode: This controls the mode of operation for dftseq.
!             here, i_mode = 1 : scalar relativistic full wave function
!  nuclear_z : nuclear charge of current species
!  el_mass   : electron mass in current units
!  dftseq_v1 : filled with 0.'s ; original meaning see dftseq
!  dftseq_v2 : filled with 1.'s ; original meaning see dftseq
!  log_deriv : would be logarithmic drivative for pp construction, here 0.
!  r_outer : outer radius from which to begin integration of radial eq.
!  n_grid_max : outermost grid point used
!  n_grid_turn : grid point index of classical turning point
!  wave_deriv1 : first derivative of wave function
!  wave_deriv2 : second derivative of wave function

!  dftseq input
   integer :: i_mode
   real*8 :: nuclear_z
   real*8 :: el_mass
   real*8 :: dftseq_v1(n_max_grid)
   real*8 :: dftseq_v2(n_max_grid)
   real*8 :: log_deriv
   real*8 :: r_outer

!  dftseq output
   integer :: n_grid_max
   integer :: n_grid_turn
   real*8 :: wave_deriv1 (n_max_grid)
   real*8 :: wave_deriv2 (n_max_grid)

!  other
 
   integer converged_dftatom
   real*8 :: alpha_grid
   real*8 :: int_outer, int_outer_l, int_outer_s
   real*8 :: r_prime(n_max_grid)

   character*100 :: info_str

!  counters

   integer :: i_grid, i_conf

!  FIXME - could be eliminated
   integer :: i_l, i_shell, kappa, relat

!  functions

   character :: l_to_str
   real*8 :: cutoff_pot

!  begin work

!  initialize fixed content for dftseq()
!  species-independent
   if (flag_rel.eq.0) then
     i_mode=2
   else if ((flag_rel.eq.REL_zora).or.(flag_rel.eq.REL_atomic_zora)) &
          then
     i_mode = 9
   end if
   el_mass = 1.0d0
   do i_grid = 1, n_max_grid, 1
     dftseq_v1(i_grid) = 0.0d0
     dftseq_v2(i_grid) = 1.0d0
   enddo
   log_deriv = 0.0d0
!  species-dependent
   nuclear_z = species_z(i_species)

!  grid scale factor to amend derivatives
   alpha_grid = log(r_grid_inc(i_species))
   do i_grid=1, n_grid(i_species)
     r_prime(i_grid) = alpha_grid * r_grid(i_grid,i_species)
   enddo

!  solve Schroedinger equation, shell by shell
   do i_conf = 1, n_conf(i_species), 1

     if (atomic_solver .eq. ATOMIC_SOLVER_ATOM_SPHERE) then
       write(info_str, '(2X,A)') "You have chosen the atom_sphere atomic solver and are using confined basis elements."
       call localorb_info( info_str )
       write(info_str, '(2X,A)') "Confined basis elements require the free atomic potential derived from get_free_atom"
       call localorb_info( info_str )
       write(info_str, '(2X,A)') "be valid, but this will not be the case for hybrid/HF calculations.  Since this"
       call localorb_info( info_str )
       write(info_str, '(2X,A)') "functionality has not been tested anyway, we're stopping the calculation for now."
       call localorb_info( info_str )
       ! If you are reading this, this could be very easily solved by calling atom_sphere_wrapper here, but using an LDA/GGA
       ! functional, for which we know that the resulting atomic potential will be good (similar to how we handle ionic basis
       ! functions with atom_sphere, where the same issue arises.)
       call aims_stop("Feel free to comment out this stop, but you're on your own.")
     end if
      
     conf_cutoff(i_species,i_conf) = conf_rad(i_species,i_conf) 
     i_l     = conf_l(i_species, i_conf)
     i_shell = conf_n(i_species, i_conf)

     if (flag_rel.eq.REL_own .or. &
          (flag_rel.eq.REL_KOLNING_HARMON)) then
       if (i_shell.le.core_n_max(i_l,i_species)) then
         i_mode = 10
       else
         i_mode = 9
       end if
     end if


     ! KH-separate core states:
     if (flag_KH_core_states ) then
        if (i_shell.le.core_n_max(i_l,i_species)) then
           i_mode = 10
        else
           if (flag_rel.eq.0) then
              i_mode=2
           else if ((flag_rel.eq.REL_zora).or.(flag_rel.eq.REL_atomic_zora)) &
                then
              i_mode = 9
           end if
        end if
     end if



!    amend basis-defining potential by specified cutoff
     do i_grid = 1, n_grid(i_species), 1
       confined_pot(i_grid,i_species,i_conf) = &
         free_potential(i_grid, i_species) + &
         cutoff_pot &
         ( r_grid(i_grid, i_species), cutoff_type(i_species), &
           conf_rad(i_species, i_conf), w_cutoff(i_species), &
           scale_cutoff(i_species))
     enddo

!    specify outermost integration point of Schroedinger Eqn, using
!    cutoff potential
     r_outer = conf_rad(i_species,i_conf) + w_cutoff(i_species)

     if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then

       kappa = conf_kappa(i_species, i_conf)
       if( kappa .eq. -i_l-1 )then
         relat=2  ! spin up
       elseif( kappa .eq. i_l )then
         relat=3  ! spin down
       else
         write(use_unit,*)'illegal kappa:',' get_confined_basis_fns'
         stop
       endif

       ! Solve Dirac equation, to get the large and small component radial wavefunctions
       call solve_radial_eigenproblem( n_grid(i_species), i_shell, i_l, -10.d0, 1.d-10, 100, &
         r_grid(1,i_species), r_prime, confined_pot(1,i_species,i_conf), nuclear_z, &
         light_speed, relat, .false., -1.d4, 1.d3, &
         converged_dftatom, confined_eigenval(i_species, i_conf), &
         confined_wave_large(1,i_species,i_conf), confined_wave_small(1,i_species,i_conf) )
     else ! non-rel case

       call dftseq ( i_mode, nuclear_z, n_grid(i_species), r_grid(1,i_species), &
         i_shell, i_l, el_mass, confined_pot(1,i_species,i_conf), &
         dftseq_v1, dftseq_v2, n_grid_max, n_grid_turn, log_deriv, &
         confined_eigenval(i_species, i_conf), confined_wave(1, i_species, i_conf), &
         wave_deriv1, wave_deriv2, r_outer )
     endif

     ! calculate new onset of cutoff potential if desired,
     ! rerun integrator for wave function with new potential if necessary

     if (basis_dep_cutoff_thresh(i_species).gt.0d0) then
       ! determine the shortest cutoff potential required by the threshold
       i_grid = n_grid(i_species)

       ! integrate the wave funtion from the outside to calculate if/when the set threshold for the cutoff 
       ! potential is reached
       if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
         int_outer_l = alpha_grid*r_grid(i_grid,i_species)*confined_wave_large(i_grid,i_species,i_conf)**2d0
         int_outer_s = alpha_grid*r_grid(i_grid,i_species)*confined_wave_small(i_grid,i_species,i_conf)**2d0
         do while ( int_outer_l.lt.basis_dep_cutoff_thresh(i_species) .and. &
                    int_outer_s.lt.basis_dep_cutoff_thresh(i_species) )
           i_grid = i_grid - 1
           int_outer_l = int_outer_l + &
                         alpha_grid*r_grid(i_grid,i_species)*confined_wave_large(i_grid,i_species,i_conf)**2d0
           int_outer_s = int_outer_s + &
                         alpha_grid*r_grid(i_grid,i_species)*confined_wave_small(i_grid,i_species,i_conf)**2d0
         end do
       else
         int_outer = alpha_grid*r_grid(i_grid,i_species)*confined_wave(i_grid,i_species,i_conf)**2d0
         do while (int_outer.lt.basis_dep_cutoff_thresh(i_species))
            i_grid = i_grid - 1 
            int_outer = int_outer + alpha_grid*r_grid(i_grid,i_species)*confined_wave(i_grid,i_species,i_conf)**2d0
         end do
       endif
       i_grid = i_grid + 1 ! the index determined in the last step does NOT satisfy the threshold, use one index higher
        
       conf_cutoff(i_species,i_conf) = min(r_grid(i_grid,i_species),r_cutoff(i_species))  ! set cutoff according to smaller r_cut
 
       ! recalculate basis function if necessary - same procedure as above
       if (conf_cutoff(i_species,i_conf).lt.conf_rad(i_species,i_conf)) then
          do i_grid = 1, n_grid(i_species), 1
             confined_pot(i_grid,i_species,i_conf) = free_potential(i_grid, i_species) + &
                cutoff_pot ( r_grid(i_grid, i_species), cutoff_type(i_species), &
                conf_cutoff(i_species, i_conf), w_cutoff(i_species), scale_cutoff(i_species))
          enddo
          
          ! specify outermost integration point of Schroedinger Eqn, using cutoff potential
          r_outer = conf_cutoff(i_species,i_conf) + w_cutoff(i_species)
          
         if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
           call solve_radial_eigenproblem( n_grid(i_species), i_shell, i_l, -10.d0, 1.d-10, 100, &
             r_grid(1,i_species), r_prime, confined_pot(1, i_species, i_conf), nuclear_z, &
             light_speed, relat, .false., -1.d4, 1.d3, &
             converged_dftatom, confined_eigenval(i_species, i_conf), &
             confined_wave_large(1,i_species,i_conf), confined_wave_small(1,i_species,i_conf) )
         else
           call dftseq ( i_mode, nuclear_z, n_grid(i_species), r_grid(1,i_species), &
                i_shell, i_l, el_mass, confined_pot(1,i_species,i_conf), &
                dftseq_v1, dftseq_v2, n_grid_max, n_grid_turn, log_deriv, &
                confined_eigenval(i_species, i_conf), confined_wave(1, i_species, i_conf), &
                wave_deriv1, wave_deriv2, r_outer )
         endif

       end if
     end if   ! end basis function dependent cutoff potential

!    store non-relativistic kinetic energy radial function
     if (flag_rel.eq.0) then
       do i_grid = 1, n_grid(i_species),1
         confined_kinetic(i_grid, i_species, i_conf) = &
           confined_eigenval(i_species, i_conf) - confined_pot(i_grid, i_species,i_conf)

         confined_kinetic(i_grid, i_species, i_conf) = &
           confined_kinetic(i_grid, i_species, i_conf) * confined_wave(i_grid, i_species, i_conf)
       enddo
     elseif (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) then
       do i_grid = 1, n_grid(i_species)
         confined_kinetic(i_grid,i_species,i_conf) = &
           confined_eigenval(i_species,i_conf) - confined_pot(i_grid,i_species,i_conf)
         confined_kinetic(i_grid,i_species,i_conf) = &
           confined_kinetic(i_grid,i_species,i_conf) * confined_wave_large(i_grid,i_species,i_conf)
       enddo
     else
       do i_grid = 1, n_grid(i_species),1
!        rescale wave_deriv2 [which is d^2(u)/di^2 on the logarithmic
!        grid i(r)] to become d^2(u)/dr^2

         wave_deriv2(i_grid) =  wave_deriv2(i_grid) - alpha_grid * wave_deriv1(i_grid)

         wave_deriv2(i_grid) = wave_deriv2(i_grid) / (alpha_grid * r_grid(i_grid,i_species))**2.d0

         confined_kinetic(i_grid,i_species,i_conf) = - 0.5 * wave_deriv2(i_grid) + &
           0.5 * conf_l(i_species,i_conf) * (conf_l(i_species,i_conf) + 1) &
           / (r_grid(i_grid,i_species)**2.d0) * confined_wave(i_grid,i_species,i_conf)
       enddo
     end if

     if (use_basis_gradients .and. flag_rel.ne.REL_x2c .and. flag_rel.ne.REL_4c_dks) then
     ! (Rundong) Only for non-rel case now. For x2c and 4c-dks, currently I see no necessity to 
     ! generate the basis derivative.
       do i_grid = 1, n_grid(i_species), 1
!        store derivative; rescale to give du/dr instead of du/di
!        where i=i(r) is the index on the logarithmic grid ...
         confined_wave_deriv(i_grid, i_species, i_conf) = &
           wave_deriv1(i_grid) / (alpha_grid*r_grid(i_grid,i_species))
       enddo
     end if

!test - create new outermost radius array after the fact, for each radial function
     i_grid = n_grid(i_species)
     if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
         do while (abs(confined_wave_large(i_grid,i_species,i_conf)).lt.wave_threshold .and. &
                   abs(confined_wave_small(i_grid,i_species,i_conf)).lt.wave_threshold)
            i_grid = i_grid - 1 
         end do
     else
         do while (abs(confined_wave(i_grid,i_species,i_conf)).lt.wave_threshold)
            i_grid = i_grid - 1 
         end do
     endif
     confined_outer_radius(i_species,i_conf) = r_grid(i_grid,i_species)
!end test

   enddo

!  output confined basis data

   call localorb_info('',use_unit)

   write(info_str,'(2X,A)') "List of confined basis orbitals and eigenvalues: "
   call localorb_info(info_str,use_unit,'(A)')

   if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) then
     write(info_str,'(4X,A,4X,A,4X,A,6X,A,4X,A,4X,A)') "n", "l", "k", "energy [Ha]", &
       "energy [eV]", "outer radius [A]"
     call localorb_info(info_str,use_unit,'(A)')

     do i_conf = 1, n_conf(i_species), 1
       i_l     = conf_l(i_species, i_conf)
       i_shell = conf_n(i_species, i_conf)
       kappa = conf_kappa(i_species, i_conf)
       write(info_str,'(2X,I3,2X,I3,2X,I3,2X,F15.6,F15.4,F15.6)') i_shell, i_l, kappa, &
       confined_eigenval(i_species, i_conf), confined_eigenval(i_species, i_conf)*hartree, &
       confined_outer_radius(i_species, i_conf)*bohr
       call localorb_info(info_str,use_unit,'(A)')
     enddo
     call localorb_info('',use_unit)
   else

     write(info_str,'(4X,A,4X,A,6X,A,4X,A,4X,A)') "n", "l", "energy [Ha]", &
       "energy [eV]", "outer radius [A]"
     call localorb_info(info_str,use_unit,'(A)')

     do i_conf = 1, n_conf(i_species), 1
       i_l     = conf_l(i_species, i_conf)
       i_shell = conf_n(i_species, i_conf)
       write(info_str,'(2X,I3,2X,I3,2X,F15.6,F15.4,F15.6)') i_shell, i_l, &
         confined_eigenval(i_species, i_conf), confined_eigenval(i_species, i_conf)*hartree, &
         confined_outer_radius(i_species,i_conf)*bohr
       call localorb_info(info_str,use_unit,'(A)')
     enddo
     call localorb_info('',use_unit)
   endif

!  that's all folks

   return
 end subroutine get_confined_basis_fns
!******
