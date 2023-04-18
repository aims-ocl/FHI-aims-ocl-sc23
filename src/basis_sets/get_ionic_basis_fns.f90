
!  FIXME - this is exactly the same as get_confined_basis_fns, only different
!  FIXME - free input potential format!

!****s* FHI-aims/get_ionic_basis_fns
!  NAME
!   get_ionic_basis_fns
!  SYNOPSIS

 subroutine get_ionic_basis_fns( i_species, free_ion_pot )

! PURPOSE
!  Subroutine get_ionic_basis_fns integrates Schroedinger equation
!  for one species, to provide radial functions for a given effective potential
! USES
   use constants,       only : bohr, hartree, light_speed, light_speed_sq
   use dimensions,      only : n_max_grid, use_basis_gradients, use_ext_basis
   use runtime_choices, only : flag_KH_core_states, flag_rel, REL_atomic_zora, &
                               REL_kolning_harmon, REL_own, REL_zora, &
                               REL_x2c, REL_4c_dks, wave_threshold
   use grids,           only : n_grid, r_grid, r_grid_inc
   use species_data,    only : ionic_cutoff, scale_cutoff, ionic_rad, ionic_pot, & 
                               ionic_kinetic, ionic_kinetic_small, &
                               ionic_wave, ionic_wave_large, ionic_wave_small, &
                               ionic_wave_deriv, ionic_large_deriv, ionic_small_deriv, &
                               ionic_outer_radius, n_ionic, ionic_l, ionic_n, ionic_kappa, &
                               ionic_in_large_basis, species_z, core_n_max, cutoff_type, &
                               w_cutoff, ionic_eigenval, basis_dep_cutoff_thresh
   use localorb_io,     only : use_unit, localorb_info
   use psi_at_nucleus_mod, only: psi_at_nucleus_ionic
   use reigen_dftatom
   implicit none
! ARGUMENTS
   integer :: i_species
   real*8  :: free_ion_pot(n_max_grid)
! INPUTS 
!   o i_species -- number of species to be calculated
!   o free_ion_pot -- free ion potential from prev. calculation
! OUTPUTS
!   none
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

!  local variables

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
   real*8 :: wave_deriv1 (n_max_grid), large_deriv1(n_max_grid), small_deriv1(n_max_grid)
   real*8 :: wave_deriv2 (n_max_grid)
   real*8 :: spline_coef (4*n_max_grid) ! for getting basis derivative

   ! For determining the wavefunction value at the nucleus. We do
   ! this by extrapolating the wavefunction into the nucleus.
   integer, parameter :: NUM_POINTS = 5
   ! Work variables for dgetrf and dgetrs. aA * aX = wave
   real*8 :: aA(NUM_POINTS,NUM_POINTS), aX(NUM_POINTS)
   integer :: aP(NUM_POINTS)
   integer :: i_exp

!  other

   integer converged_dftatom
   real*8 :: alpha_grid
   real*8 :: int_outer, int_outer_l, int_outer_s
   real*8 :: r_prime(n_max_grid)
   real*8 :: tmp, tmp1

   character*100 :: info_str

!  counters

   integer :: i_grid, i_ionic

!  FIXME, could eliminate
   integer :: i_l, i_n, kappa, relat

!  functions

   character :: l_to_str
   real*8 :: cutoff_pot

!  begin work

!  initialize fixed content for dftseq()
!  species-independent
   if (flag_rel.eq.0) then
     i_mode=2
   else if ((flag_rel.eq.REL_zora).or.(flag_rel.eq.REL_atomic_zora)) then
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
   ionic_cutoff(i_species,:) = ionic_rad(i_species,:)

!  grid scale factor to amend derivatives
   alpha_grid = log(r_grid_inc(i_species))
   do i_grid=1, n_grid(i_species)
     r_prime(i_grid) = alpha_grid * r_grid(i_grid,i_species)
   enddo

!  solve Schroedinger or Dirac equation, shell by shell
   do i_ionic = 1, n_ionic(i_species), 1

     i_l = ionic_l(i_species, i_ionic)
     i_n = ionic_n(i_species, i_ionic)

     if (flag_rel.eq.REL_own.or. (flag_rel.eq.REL_KOLNING_HARMON)) then
       if (i_n.le.core_n_max(i_l,i_species)) then
         i_mode = 10
       else
         i_mode = 9
       end if
     end if

     ! KH-separate core states:
     if (flag_KH_core_states ) then           
        if (i_n.le.core_n_max(i_l,i_species)) then
           i_mode = 10
        else
           if (flag_rel.eq.0) then
              i_mode=2
           else if ((flag_rel.eq.REL_zora).or.(flag_rel.eq.REL_atomic_zora))then
              i_mode = 9
           end if
        end if
     end if

!    amend basis-defining potential by specified cutoff
     do i_grid = 1, n_grid(i_species), 1
       ionic_pot(i_grid, i_species, i_ionic) = free_ion_pot(i_grid) + &
         cutoff_pot ( r_grid(i_grid, i_species), cutoff_type(i_species), &
         ionic_rad(i_species, i_ionic), w_cutoff(i_species), scale_cutoff(i_species) )
     enddo

!    specify outermost integration point of Schroedinger Eqn, using
!    cutoff potential
     r_outer = ionic_rad(i_species,i_ionic) + w_cutoff(i_species)

     if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
       kappa = ionic_kappa(i_species, i_ionic)
       if( kappa .eq. -i_l-1 )then
         relat=2  ! spin up
       elseif( kappa .eq. i_l )then
         relat=3  ! spin down
       else
         write(use_unit,*)'illegal kappa:',' get_ionic_basis_fns'
         stop
       endif
       ! Solve Dirac equation, to get the large and small component radial wavefunctions
       call solve_radial_eigenproblem( n_grid(i_species), i_n, i_l, -10.d0, 1.d-10, 100, &
         r_grid(1,i_species), r_prime, ionic_pot(1, i_species, i_ionic), nuclear_z, &
         light_speed, relat, .true., -1.d4, 10.d0, &
         converged_dftatom, ionic_eigenval(i_species, i_ionic), &
         ionic_wave_large(1,i_species,i_ionic), ionic_wave_small(1,i_species,i_ionic) )
       ! Then, get the derivative of basis:
       call fderiv(1, n_grid(i_species),r_grid(1,i_species), ionic_wave_large(1,i_species,i_ionic), large_deriv1, spline_coef)
       call fderiv(1, n_grid(i_species),r_grid(1,i_species), ionic_wave_small(1,i_species,i_ionic), small_deriv1, spline_coef)
     else ! non-rel case
       call dftseq ( i_mode, nuclear_z, n_grid(i_species), r_grid(1,i_species), &
         i_n, i_l, el_mass, ionic_pot(1, i_species, i_ionic), &
         dftseq_v1, dftseq_v2, n_grid_max, n_grid_turn, log_deriv, &
         ionic_eigenval(i_species, i_ionic), ionic_wave(1,i_species,i_ionic), &
         wave_deriv1, wave_deriv2, r_outer )
     endif

!    calculate new onset of cutoff potential if desired,
!    rerun integrator for wave function with new potential if necessary

     if (basis_dep_cutoff_thresh(i_species).gt.0d0) then
       ! determine the shortest cutoff potential required by the threshold
       i_grid = n_grid(i_species)

       ! Integrate the wave funtion from the outside to calculate if/when 
       ! the set threshold for the cutoff potential is reached.
       if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
         int_outer_l = alpha_grid*r_grid(i_grid,i_species)*ionic_wave_large(i_grid,i_species,i_ionic)**2d0
         int_outer_s = alpha_grid*r_grid(i_grid,i_species)*ionic_wave_small(i_grid,i_species,i_ionic)**2d0
         do while ( int_outer_l.lt.basis_dep_cutoff_thresh(i_species) .and. &
                    int_outer_s.lt.basis_dep_cutoff_thresh(i_species) )
           i_grid = i_grid - 1
           int_outer_l = int_outer_l + &
                         alpha_grid*r_grid(i_grid,i_species)*ionic_wave_large(i_grid,i_species,i_ionic)**2d0
           int_outer_s = int_outer_s + &
                         alpha_grid*r_grid(i_grid,i_species)*ionic_wave_small(i_grid,i_species,i_ionic)**2d0
         end do
       else
         int_outer = alpha_grid*r_grid(i_grid,i_species)*ionic_wave(i_grid,i_species,i_ionic)**2d0
         do while (int_outer.lt.basis_dep_cutoff_thresh(i_species))
           i_grid = i_grid - 1 
           int_outer = int_outer + &
                       alpha_grid*r_grid(i_grid,i_species)*ionic_wave(i_grid,i_species,i_ionic)**2d0
         end do
       endif
       i_grid = i_grid + 1 ! the index determined in the last step does NOT satisfy the threshold
       
       ionic_cutoff(i_species,i_ionic) = min(r_grid(i_grid,i_species),ionic_rad(i_species,i_ionic))
 
       ! recalculate basis function if necessary - same procedure as above
       if (ionic_cutoff(i_species,i_ionic).lt.ionic_rad(i_species,i_ionic)) then

         do i_grid = 1, n_grid(i_species), 1
           ionic_pot(i_grid, i_species, i_ionic) = free_ion_pot(i_grid) + &
              cutoff_pot ( r_grid(i_grid, i_species), cutoff_type(i_species), &
              ionic_cutoff(i_species, i_ionic), w_cutoff(i_species), &
              scale_cutoff(i_species) )
         enddo
         
         r_outer = ionic_cutoff(i_species,i_ionic) + w_cutoff(i_species)

         if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
           call solve_radial_eigenproblem( n_grid(i_species), i_n, i_l, -10.d0, 1.d-10, 100, &
             r_grid(1,i_species), r_prime, ionic_pot(1, i_species, i_ionic), nuclear_z, &
             light_speed, relat, .true., -1.d4, 10.d0, &
             converged_dftatom, ionic_eigenval(i_species, i_ionic), &
             ionic_wave_large(1,i_species,i_ionic), ionic_wave_small(1,i_species,i_ionic) )
           call fderiv(1, n_grid(i_species),r_grid(1,i_species), ionic_wave_large(1,i_species,i_ionic), large_deriv1, spline_coef)
           call fderiv(1, n_grid(i_species),r_grid(1,i_species), ionic_wave_small(1,i_species,i_ionic), small_deriv1, spline_coef)
         else 
           call dftseq ( i_mode, nuclear_z, n_grid(i_species), r_grid(1,i_species), &
                i_n, i_l, el_mass, ionic_pot(1,i_species,i_ionic), &
                dftseq_v1, dftseq_v2, n_grid_max, n_grid_turn, log_deriv, &
                ionic_eigenval(i_species,i_ionic), ionic_wave(1,i_species,i_ionic), &
                wave_deriv1, wave_deriv2, r_outer )
         endif

       end if
     end if   ! end basis function dependent cutoff potential

!    store non-relativistic kinetic energy radial function
     if (flag_rel.eq.0) then
       do i_grid = 1, n_grid(i_species)
         ionic_kinetic(i_grid, i_species, i_ionic) = &
           ionic_eigenval(i_species, i_ionic) - ionic_pot(i_grid, i_species,i_ionic)

         ionic_kinetic(i_grid, i_species, i_ionic) = &
           ionic_kinetic(i_grid, i_species, i_ionic) * ionic_wave(i_grid, i_species, i_ionic)
       enddo
     elseif (flag_rel.eq.REL_x2c) then
       do i_grid = 1, n_grid(i_species)
         tmp = ionic_eigenval(i_species,i_ionic) - ionic_pot(i_grid,i_species,i_ionic)
         tmp1 = tmp/light_speed
         ionic_kinetic(i_grid,i_species,i_ionic) = ( tmp + 0.5d0*tmp1*tmp1 ) * ionic_wave_large(i_grid,i_species,i_ionic)
         ionic_kinetic(i_grid,i_species,i_ionic) = ionic_eigenval(i_species,i_ionic) - ionic_pot(i_grid,i_species,i_ionic)
         ionic_kinetic(i_grid,i_species,i_ionic) = ionic_kinetic(i_grid,i_species,i_ionic) * ionic_wave_large(i_grid,i_species,i_ionic)
         ionic_kinetic_small(i_grid,i_species,i_ionic) = ionic_eigenval(i_species,i_ionic) - ionic_pot(i_grid,i_species,i_ionic) + 2.d0*light_speed_sq
         ionic_kinetic_small(i_grid,i_species,i_ionic) = ionic_kinetic_small(i_grid,i_species,i_ionic) * ionic_wave_small(i_grid,i_species,i_ionic)
       enddo
     elseif (flag_rel.eq.REL_4c_dks) then
       ionic_kinetic(:,i_species,i_ionic) = ionic_eigenval(i_species,i_ionic) - ionic_pot(:,i_species,i_ionic)
       ionic_kinetic(:,i_species,i_ionic) = ionic_kinetic(:,i_species,i_ionic) * ionic_wave_large(:,i_species,i_ionic)
       if(.true.)then ! use atomic balance condition
       ionic_kinetic_small(:,i_species,i_ionic) = ionic_eigenval(i_species,i_ionic) - ionic_pot(:,i_species,i_ionic) + 2.d0*light_speed_sq
       ionic_kinetic_small(:,i_species,i_ionic) = ionic_kinetic_small(:,i_species,i_ionic) * ionic_wave_small(:,i_species,i_ionic)
     ! else ! use restricted kinentic balance condtion
     ! !(Rundong) in this way, the kinetic energy term is actually in a nonrelativistic form, viz. the second derivative of the basis. 
     ! ! The precision is low. Thus, I don't recommend this strategy.
     ! call sigma_dot_p_large_wave( n_grid(i_species), r_grid(1,i_species), ionic_l(i_species,i_ionic), ionic_kappa(i_species,i_ionic),&
     !      ionic_wave_large(1,i_species,i_ionic), large_deriv1, ionic_wave_small(1,i_species,i_ionic), ionic_kinetic_small(1,i_species,i_ionic) )
     ! call fderiv(1,n_grid(i_species),r_grid(1,i_species), ionic_wave_small(1,i_species,i_ionic),small_deriv1,spline_coef)
     !!call fderiv(2,n_grid(i_species),r_grid(1,i_species), ionic_wave_large(1,i_species,i_ionic),ionic_kinetic(1,i_species,i_ionic),spline_coef)

     !! if use zora kind of atomic balance
     ! ionic_wave_small(:,i_species,i_ionic) = 2*light_speed_sq/(2*light_speed_sq-ionic_pot(:,i_species,i_ionic)) * ionic_wave_small(:,i_species,i_ionic)
     ! ionic_kinetic(:,i_species,i_ionic) = -ionic_pot(:,i_species,i_ionic) * ionic_wave_large(:,i_species,i_ionic)
     ! ionic_kinetic_small(:,i_species,i_ionic) = (2.d0*light_speed_sq-ionic_pot(:,i_species,i_ionic)) * ionic_wave_small(:,i_species,i_ionic)
       endif
     else
       do i_grid = 1, n_grid(i_species)
!        rescale wave_deriv2 [which is d^2(u)/di^2 on the logarithmic
!        grid i(r)] to become d^2(u)/dr^2

         wave_deriv2(i_grid) = wave_deriv2(i_grid) - alpha_grid * wave_deriv1(i_grid)
         wave_deriv2(i_grid) = wave_deriv2(i_grid) / (alpha_grid * r_grid(i_grid,i_species))**2.d0

         ionic_kinetic(i_grid,i_species,i_ionic) = - 0.5 * wave_deriv2(i_grid) + &
           0.5 * ionic_l(i_species,i_ionic) * (ionic_l(i_species,i_ionic) + 1) &
           / (r_grid(i_grid,i_species)**2.d0) * ionic_wave(i_grid,i_species,i_ionic)
       enddo
     end if

     if (use_basis_gradients)then
        do i_grid = 1, n_grid(i_species), 1
!         store derivative; rescale to give du/dr instead of du/di
!         where i=i(r) is the index on the logarithmic grid ...
          ionic_wave_deriv(i_grid, i_species, i_ionic) = &
            wave_deriv1(i_grid) / (alpha_grid*r_grid(i_grid,i_species))
        enddo
     end if
     if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) then
        do i_grid = 1, n_grid(i_species), 1
          ionic_large_deriv(i_grid, i_species, i_ionic) = large_deriv1(i_grid)
          ionic_small_deriv(i_grid, i_species, i_ionic) = small_deriv1(i_grid)
        enddo
     endif

!test - create new outermost radius array after the fact, for each radial function
     i_grid = n_grid(i_species)
     if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
         do while (abs(ionic_wave_large(i_grid,i_species,i_ionic)).lt.wave_threshold .and. &
                   abs(ionic_wave_small(i_grid,i_species,i_ionic)).lt.wave_threshold)
            i_grid = i_grid - 1 
         end do
     else
         do while (abs(ionic_wave(i_grid,i_species,i_ionic)).lt.wave_threshold)
            i_grid = i_grid - 1 
         end do
     endif
     ionic_outer_radius(i_species,i_ionic) = r_grid(i_grid,i_species)
!end test

     if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) goto 1000 ! (Rundong) I don't know what the following
     ! code is used for. Currently, skip it for relativistic cases. If you want to modify it for relativistic 
     ! basis, please replace the ionic_wave with ionic_wave_large (or maybe ionic_wave_small? -- I don't know).

         ! Only the s orbitals may have nonzero values at the nucleus.
         if (i_l == 0) then
            ! A is the following matrix:
            !   r1^0 r1^1 ... r1^(N-1)
            !   r2^0 r2^1 ... r2^(N-1)
            !   ...
            !   rN^0 rN^1 ... rN^(N-1)
            aA(:,1) = 1d0
            do i_exp = 2, NUM_POINTS
               aA(:,i_exp) = r_grid(:NUM_POINTS,i_species)**(i_exp-1)
            end do
            aX = ionic_wave(:NUM_POINTS,i_species,i_ionic) / r_grid(:NUM_POINTS,i_species)
            call dgetrf(size(aA,1), size(aA,2), aA, size(aA,1), aP, i_exp)
            call dgetrs('n', size(aA,1), 1, aA, size(aA,1), aP, aX, size(aX), i_exp)
            psi_at_nucleus_ionic(i_ionic,i_species) = aX(1)
         else
            psi_at_nucleus_ionic(i_ionic,i_species) = 0d0
         end if

 1000 continue
   enddo ! end of i_ionic

!  output ion basis data

   if(.not. all(ionic_in_large_basis(i_species, 1:n_ionic(i_species)))) then
     call localorb_info('',use_unit)

     write(info_str,'(2X,A)') "List of ionic basis orbitals and eigenvalues: "
     call localorb_info(info_str,use_unit,'(A)')

     if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) then
       write(info_str,'(4X,A,4X,A,4X,A,6X,A,4X,A,4X,A)') "n", "l", "k", "energy [Ha]", &
         "energy [eV]", "outer radius [A]"
       call localorb_info(info_str,use_unit,'(A)')


       do i_ionic = 1, n_ionic(i_species), 1
         if( ionic_in_large_basis(i_species,i_ionic)) cycle
         i_l = ionic_l(i_species, i_ionic)
         i_n = ionic_n(i_species, i_ionic)
         kappa = ionic_kappa(i_species, i_ionic)
         write(info_str,'(2X,I3,2X,I3,2X,I3,2X,F15.6,F15.4,F15.6)') i_n, i_l, kappa, &
         ionic_eigenval(i_species, i_ionic), ionic_eigenval(i_species, i_ionic)*hartree, &
         ionic_outer_radius(i_species, i_ionic)*bohr
         call localorb_info(info_str,use_unit,'(A)')

       enddo
       call localorb_info('',use_unit)
     else
       write(info_str,'(4X,A,4X,A,6X,A,4X,A,4X,A)') "n", "l", "energy [Ha]", &
         "energy [eV]", "outer radius [A]"
       call localorb_info(info_str,use_unit,'(A)')


       do i_ionic = 1, n_ionic(i_species), 1
         if( ionic_in_large_basis(i_species,i_ionic)) cycle
         i_l = ionic_l(i_species, i_ionic)
         i_n = ionic_n(i_species, i_ionic)

         write(info_str,'(2X,I3,2X,I3,2X,F15.6,F15.4,F15.6)') i_n, i_l, &
         ionic_eigenval(i_species, i_ionic), ionic_eigenval(i_species, i_ionic)*hartree, &
         ionic_outer_radius(i_species, i_ionic)*bohr
         call localorb_info(info_str,use_unit,'(A)')

       enddo
       call localorb_info('',use_unit)
     endif
   end if
   
   ! print out the extra basis functions if present
   if(use_ext_basis .and. any(ionic_in_large_basis(i_species,1:n_ionic(i_species)))) then
     call localorb_info('',use_unit)

     write(info_str,'(2X,A)')"List of extra ionic orbitals and eigenvalues for auxliary basis: "
     call localorb_info(info_str,use_unit,'(A)')

     write(info_str,'(4X,A,4X,A,6X,A,4X,A,4X,A)') "n", "l", "energy [Ha]", &
       "energy [eV]", "outer radius [A]"
     call localorb_info(info_str,use_unit,'(A)')


     do i_ionic = 1, n_ionic(i_species), 1
       if(.not. ionic_in_large_basis(i_species,i_ionic)) cycle
       i_l = ionic_l(i_species, i_ionic)
       i_n = ionic_n(i_species, i_ionic)

       write(info_str,'(2X,I3,2X,I3,2X,F15.6,F15.4,F15.6)') i_n, i_l, &
       ionic_eigenval(i_species, i_ionic), ionic_eigenval(i_species, i_ionic)*hartree, &
       ionic_outer_radius(i_species, i_ionic)*bohr
       call localorb_info(info_str,use_unit,'(A)')

     enddo
     call localorb_info('',use_unit)
   
   end if
!  that's all folks

   return
 end subroutine get_ionic_basis_fns

