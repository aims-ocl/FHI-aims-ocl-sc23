!****s* FHI-aims/get_species_basis_fns
!  NAME
!   get_species_basis_fns
!  SYNOPSIS
 subroutine get_species_basis_fns ( ) 
!  PURPOSE
!    Subroutine get_species_basis_fns integrates Schroedinger equation to provide
!    the atomic-like basis part for the current effective potential
!  USES
   use constants,       only : bohr, hartree, light_speed, light_speed_sq
   use dimensions,      only : n_max_grid, n_max_ind_fns, n_species, use_basis_gradients
   use runtime_choices, only : flag_KH_core_states, flag_rel, REL_kolning_harmon, &
                               REL_none, REL_own, REL_x2c, REL_4c_dks, wave_threshold
   use grids,           only : n_grid, r_grid, r_grid_inc, log_r_grid_inc
   use species_data,    only : atomic_cutoff, basis_potential, atomic_kinetic, atomic_kinetic_small, &
                               atomic_wave_deriv, atomic_large_deriv, atomic_small_deriv, &
                               atomic_outer_radius, include_min_basis, &
                               species_z, core_fn, cut_core, w_cutoff, core_r_cut, &
                               cutoff_type, r_cutoff, atomic_n, cut_atomic_basis_functions, &
                               atomic_wave, atm_wave_large, atm_wave_small, &
                               basis_dep_cutoff_thresh, scale_cutoff, &
                               atomic_l, atomic_k, atomic_eigenval, n_atomic, species_name  
   use free_atoms,      only : free_potential
   use localorb_io,     only : use_unit, localorb_info
   use reigen_dftatom
   implicit none
!  INPUT
!   none
!  OUTPUT
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
!  imported variables

!  local variables

!  we solve a 1D equation, repeatedly.
!  Out of sheer paranoia, we copy all interesting arrays onto 1D working
!  arrays.

!  n_shell : denotes maximum main quantum number for given l

   integer n_shell

!  auxiliary variables for dftseq :

!  input:
!  i_mode: This controls the mode of operation for dftseq.
!          here, i_mode = 1 : scalar relativistic full wave function
!  nuclear_z : nuclear charge of current species
!  el_mass   : electron mass in current units
!  dftseq_v1 : filled with 0.'s ; original meaning see dftseq
!  dftseq_v2 : filled with 1.'s ; original meaning see dftseq
!  log_deriv : would be logarithmic drivative for pp construction, here 0.

!  output:
!  n_grid_max : outermost grid point used
!  n_grid_turn : grid point index of classical turning point
!  wave_deriv1 : first derivative of wave function
!  wave_deriv2 : second derivative of wave function
!
!  vb variables:
!  r_outer : outer radius from which to begin integration of radial eq.

   integer :: i_mode,converged_dftatom
   real*8 :: nuclear_z
   real*8 :: el_mass
   real*8 :: dftseq_v1(n_max_grid)
   real*8 :: dftseq_v2(n_max_grid)
   real*8 :: log_deriv

   integer :: n_grid_max
   integer :: n_grid_turn
   real*8 :: wave_deriv1 (n_max_grid), large_deriv1(n_max_grid), small_deriv1(n_max_grid)
   real*8 :: wave_deriv2 (n_max_grid)
   real*8 :: spline_coef (4*n_max_grid) ! for getting basis derivative

 ! eigenvalue and wave function obtained from schrodinger equation with relativistic radial potential:
   real*8 :: nr_eigenval(n_species,n_max_ind_fns)
   real*8 :: nr_wave(n_max_grid,n_species,n_max_ind_fns)
   real*8 :: tmp, tmp1, tmp2

   real*8 :: r_outer

   real*8 :: r_prime(n_max_grid)
   real*8 :: radial_potential(n_max_grid)

   real*8 :: alpha_grid
   real*8 :: int_outer

   character*100 :: info_str

!  functions

   real*8 :: cutoff_pot

!  counters

   integer i,l, i_species, i_fn, i_grid, relat

!  begin work

   write (info_str,'(2X,A,A)') &
    "Creating atomic-like basis functions for current", &
    " effective potential."
   call localorb_info(info_str,use_unit,'(A)')

!  initialize fixed content for dftseq

   if (flag_rel.eq.0) then
     i_mode=2
   else if ((flag_rel.eq.1).or.(flag_rel.eq.2)) then
     i_mode=9
   end if

   el_mass = 1.0d0
   do i_grid = 1, n_max_grid, 1
     dftseq_v1(i_grid) = 0.0d0
     dftseq_v2(i_grid) = 1.0d0
   enddo
   log_deriv = 0.0d0

!  begin loop over all species

   do i_species = 1, n_species, 1
   if (include_min_basis(i_species)) then

!  for dftseq

     nuclear_z = species_z(i_species)

!  grid scale factor to amend derivatives

     alpha_grid = log(r_grid_inc(i_species))

     if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
       ! grid derivatives
       do i_grid=1, n_grid(i_species)
         r_prime(i_grid) = alpha_grid * r_grid(i_grid,i_species)
       enddo
     endif

!  now attempt to integrate each wavefunction n, l individually.

     do i_fn = 1, n_atomic(i_species), 1

        if (flag_rel.eq.REL_own.or.(flag_rel.eq.REL_KOLNING_HARMON)) &
             then
           if (core_fn(i_species,i_fn)) then
              i_mode = 10
           else
              i_mode = 9
           end if
        end if
        
        if (flag_KH_core_states ) then
           if (core_fn(i_species,i_fn)) then
              i_mode = 10
           else
              if (flag_rel.eq.0) then
                 i_mode=2
              else if ((flag_rel.eq.1).or.(flag_rel.eq.2)) then
                 i_mode=9
              end if
              
           end if
        end if


        ! hacked handling of core cutoff parameter for basis-defining potential.
        ! This collides with fixed_basis_potential earlier. I am only separating this
        ! here because being able to print out basis_potential later has proven
        ! tremendously useful - I do not want to destroy that possibility.


        if (core_fn(i_species,i_fn).and.cut_core(i_species)) then
          r_outer = core_r_cut(i_species) + w_cutoff(i_species)
          atomic_cutoff(i_species,i_fn) = core_r_cut(i_species)
          do i_grid = 1, n_grid(i_species), 1

            radial_potential (i_grid) = free_potential (i_grid, i_species)

            if ( (core_r_cut(i_species)+w_cutoff(i_species)).gt.0. ) then
              radial_potential (i_grid) = radial_potential (i_grid) + cutoff_pot( r_grid(i_grid,i_species), cutoff_type(i_species), &
                core_r_cut(i_species), w_cutoff(i_species), scale_cutoff(i_species) )
            end if

          enddo
        else
          ! use the standard cutoff potential
          r_outer = r_cutoff(i_species) + w_cutoff(i_species)
          atomic_cutoff(i_species,i_fn) = r_cutoff(i_species)
          radial_potential(:) = basis_potential(:,i_species)
        end if


        if(flag_rel.eq.REL_4c_dks .or. flag_rel.eq.REL_x2c)then
          ! Solve radial Dirac equation, to generate relativisitc kind basis for each shell: for fully-relativistic 4c-DKS, and X2C methods.
          if(atomic_k(i_species,i_fn).eq.-atomic_l(i_species,i_fn)-1)then
            relat=2
          elseif(atomic_k(i_species,i_fn).eq.atomic_l(i_species,i_fn))then
            relat=3
          else
            write(use_unit,*)'get_species_basis_fn: illegal kappa'
            stop
          endif
          call solve_radial_eigenproblem(n_grid(i_species), &
            atomic_n(i_species,i_fn), atomic_l(i_species,i_fn), -10.d0, 1.d-10, 100, r_grid(1,i_species), r_prime, radial_potential, nuclear_z, &
            light_speed, relat, .true., -1.d4, 10.d0, converged_dftatom, &
            atomic_eigenval(i_species,i_fn), atm_wave_large(1,i_species,i_fn), atm_wave_small(1,i_species,i_fn))
          ! For GGA, we need to get the derivative of large and small comp. basis:
          call fderiv(1, n_grid(i_species),r_grid(1,i_species), atm_wave_large(1,i_species,i_fn), large_deriv1, spline_coef)
          call fderiv(1, n_grid(i_species),r_grid(1,i_species), atm_wave_small(1,i_species,i_fn), small_deriv1, spline_coef)
          ! write(6,"('n,l,k:',i3,i3,i3,5x,e20.13)")atomic_n(i_species,i_fn), atomic_l(i_species,i_fn),atomic_k(i_species,i_fn),atomic_eigenval(i_species,i_fn)
        else
          ! Generate non-relativistic basis for each shell. Basis for ZORA are also generated here, instead in above.
          call dftseq (i_mode, nuclear_z,n_grid(i_species),r_grid(1,i_species), atomic_n(i_species,i_fn), atomic_l(i_species,i_fn), &
            el_mass, radial_potential, dftseq_v1, dftseq_v2, n_grid_max, n_grid_turn, log_deriv, &
            atomic_eigenval(i_species, i_fn), atomic_wave(1,i_species,i_fn), wave_deriv1, wave_deriv2, r_outer)
        endif


        ! calculate new onset of cutoff potential if desired,
        !   rerun integrator for wave function with new potential if necessary

        if ((basis_dep_cutoff_thresh(i_species).gt.0d0).and.cut_atomic_basis_functions(i_species)) then
           ! determine the shortest cutoff potential required by the threshold
           i_grid = n_grid(i_species)
           
           ! integrate the wave funtion from the outside to calculate if/when the set threshold for the cutoff 
           ! potential is reached
           int_outer = alpha_grid*r_grid(i_grid,i_species)*atomic_wave(i_grid,i_species,i_fn)**2d0
           do while (int_outer.lt.basis_dep_cutoff_thresh(i_species))
              i_grid = i_grid - 1 
              int_outer = int_outer + alpha_grid*r_grid(i_grid,i_species)*atomic_wave(i_grid,i_species,i_fn)**2d0
           end do
           i_grid = i_grid + 1 ! the index determined in the last step does NOT satisfy the threshold, use one index higher
           
           atomic_cutoff(i_species,i_fn) = min(r_grid(i_grid,i_species),atomic_cutoff(i_species,i_fn))  ! set cutoff according to smaller r_cut
           r_outer = atomic_cutoff(i_species,i_fn)+w_cutoff(i_species)

           ! prepare new atomic potential for entire basis function
           do i_grid = 1, n_grid(i_species), 1                  
              radial_potential (i_grid) = free_potential (i_grid, i_species) + cutoff_pot ( r_grid(i_grid,i_species), cutoff_type(i_species), &
                                          atomic_cutoff(i_species,i_fn), w_cutoff(i_species), scale_cutoff(i_species) )                  
           enddo

          if(flag_rel.eq.REL_4c_dks .or. flag_rel.eq.REL_x2c)then
            call solve_radial_eigenproblem(n_grid(i_species), &
              atomic_n(i_species,i_fn), atomic_l(i_species,i_fn), -10.d0, 1.d-10, 100, r_grid(1,i_species), r_prime, radial_potential, nuclear_z, &
              light_speed, relat, .true., -1.d4, 10.d0, converged_dftatom, &
              atomic_eigenval(i_species,i_fn), atm_wave_large(1,i_species,i_fn), atm_wave_small(1,i_species,i_fn))
            call fderiv(1, n_grid(i_species),r_grid(1,i_species), atm_wave_large(1,i_species,i_fn), large_deriv1, spline_coef)
            call fderiv(1, n_grid(i_species),r_grid(1,i_species), atm_wave_small(1,i_species,i_fn), small_deriv1, spline_coef)
            write(use_unit,"('n,l,k:',i3,i3,i3,5x,e20.13)")atomic_n(i_species,i_fn), atomic_l(i_species,i_fn),atomic_k(i_species,i_fn),atomic_eigenval(i_species,i_fn)
          else 
            call dftseq (i_mode, nuclear_z,n_grid(i_species),r_grid(1,i_species), atomic_n(i_species,i_fn), atomic_l(i_species,i_fn), &
              el_mass, radial_potential, dftseq_v1, dftseq_v2, n_grid_max, n_grid_turn, log_deriv, &
              atomic_eigenval(i_species, i_fn), atomic_wave(1,i_species,i_fn), wave_deriv1, wave_deriv2, r_outer)
          endif
           
        end if   ! end basis function dependent cutoff potential

!       store kinetic energy radial function
        if (flag_rel.eq.REL_none) then
          do i_grid = 1, n_grid(i_species),1
            atomic_kinetic(i_grid, i_species, i_fn) = atomic_eigenval(i_species, i_fn) - basis_potential(i_grid, i_species)
            atomic_kinetic(i_grid, i_species, i_fn) = atomic_kinetic(i_grid, i_species, i_fn) * atomic_wave(i_grid, i_species, i_fn)
          enddo
        else if (flag_rel.eq.REL_x2c) then
         ! (Rundong) Note: the kinetic energy operator used in DKS equation is the nonrelativistic T operator (second derivative of the wave function).
         ! Now, calculate the 2nd derivative:
         !call fderiv(2, n_grid(i_species),r_grid(1,i_species), atm_wave_large(1,i_species,i_fn), atomic_kinetic(1,i_species,i_fn), spline_coef)
         !l=atomic_l(i_species,i_fn)
         !do i_grid = 1, n_grid(i_species)
         !  !p^2/2m
         !   atomic_kinetic(i_grid,i_species,i_fn) = -0.5d0*( atomic_kinetic(i_grid,i_species,i_fn) - &
         !                                           atomic_kinetic(i_grid,i_species,i_fn)/r_grid(i_grid,i_species)**2*l*(l+1) )
         !enddo
          do i_grid = 1, n_grid(i_species),1

            tmp = atomic_eigenval(i_species, i_fn) - basis_potential(i_grid, i_species)
            tmp1 = tmp/light_speed
            atomic_kinetic(i_grid, i_species, i_fn) = ( tmp + 0.5d0*tmp1*tmp1 ) * atm_wave_large(i_grid,i_species,i_fn)
            atomic_kinetic_small(i_grid, i_species, i_fn) = atomic_eigenval(i_species,i_fn) - basis_potential(i_grid,i_species) + 2.d0*light_speed_sq
            atomic_kinetic_small(i_grid, i_species, i_fn) = atomic_kinetic_small(i_grid, i_species, i_fn) * atm_wave_small(i_grid, i_species, i_fn)
          enddo

        else if (flag_rel.eq.REL_4c_dks) then ! Q4C for NAO

          do i_grid = 1, n_grid(i_species),1
            atomic_kinetic(i_grid, i_species, i_fn) = atomic_eigenval(i_species,i_fn) - basis_potential(i_grid,i_species)
            atomic_kinetic(i_grid, i_species, i_fn) = atomic_kinetic(i_grid, i_species, i_fn) * atm_wave_large(i_grid, i_species, i_fn)
            atomic_kinetic_small(i_grid, i_species, i_fn) = atomic_eigenval(i_species,i_fn) - basis_potential(i_grid,i_species) + 2.d0*light_speed_sq
            atomic_kinetic_small(i_grid, i_species, i_fn) = atomic_kinetic_small(i_grid, i_species, i_fn) * atm_wave_small(i_grid, i_species, i_fn)
          enddo

        else
          do i_grid = 1, n_grid(i_species),1
!           rescale wave_deriv2 [which is d^2(u)/di^2 on the logarithmic
!           grid i(r)] to become d^2(u)/dr^2
            wave_deriv2(i_grid) = wave_deriv2(i_grid) - alpha_grid * wave_deriv1(i_grid)
            wave_deriv2(i_grid) = wave_deriv2(i_grid) / (alpha_grid * r_grid(i_grid,i_species))**2.d0
            atomic_kinetic(i_grid,i_species,i_fn) = - 0.5 * wave_deriv2(i_grid) + &
             0.5 * atomic_l(i_species,i_fn) * (atomic_l(i_species,i_fn) + 1) / (r_grid(i_grid,i_species)**2.d0) * atomic_wave(i_grid,i_species,i_fn)
          enddo
        end if

        if (use_basis_gradients)then
           do i_grid=1, n_grid(i_species), 1
!            store derivative; rescale to give du/dr instead of du/di
!            where i=i(r) is the index on the logarithmic grid ...
             atomic_wave_deriv(i_grid, i_species, i_fn) =  wave_deriv1(i_grid) / (alpha_grid*r_grid(i_grid,i_species))
           enddo
        end if
        if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) then
           do i_grid=1, n_grid(i_species), 1
             atomic_large_deriv(i_grid, i_species, i_fn) = large_deriv1(i_grid)
             atomic_small_deriv(i_grid, i_species, i_fn) = small_deriv1(i_grid)
           enddo
        endif

!test - create new outermost radius array after the fact, for each radial function
        i_grid = n_grid(i_species)
        if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) then
          do while (abs(atm_wave_large(i_grid,i_species,i_fn)).lt.wave_threshold .and. abs(atm_wave_small(i_grid,i_species,i_fn)).lt.wave_threshold)
             i_grid = i_grid - 1 
          end do
        else
          do while (abs(atomic_wave(i_grid,i_species,i_fn)).lt.wave_threshold)
             i_grid = i_grid - 1 
          end do
        endif
        atomic_outer_radius(i_species,i_fn) = r_grid(i_grid,i_species)
!end test

       !if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
       !   call get_small_comp_sigmadotp(i_species, n_grid(i_species), i_fn, r_grid(1,i_species), atm_wave_large(1,i_species,i_fn), large_deriv1, &
       !     atomic_eigenval(i_species,i_fn), basis_potential(1,i_species), &
       !     atm_wave_small(1,i_species,i_fn), atomic_kinetic_small(1,i_species,i_fn) ) !!!!!!!!!!!
       !endif
     enddo ! end loop of i_fn

                       !do i_fn=1, n_atomic(i_species)
                       !   write(6,*)'i_fn=',i_fn
                       !  !write(6,*)'large:'
                       !  !write(6,"(20f13.7)")atm_wave_large(:,i_species,i_fn)
                       !   write(6,*)'kinetic:'
                       !   write(6,"(15f16.7)")atomic_kinetic(:,i_species,i_fn)
                       !   write(6,*)'nr_wave:'
                       !   write(6,"(15f16.7)")nr_wave(:,i_species,i_fn)
                       !enddo

!  output atomic basis data

     call localorb_info('',use_unit)

     write(info_str,'(2X,A,A,A)') "Species ", species_name(i_species),":"
     call localorb_info(info_str,use_unit,'(A)')

     call localorb_info('',use_unit)

     write(info_str,'(2X,A)') "List of atomic basis orbitals and eigenvalues: "
     call localorb_info(info_str,use_unit,'(A)')

     if(flag_rel.eq.REL_4c_dks .or. flag_rel.eq.REL_x2c)then
       write(info_str,'(4X,A,4X,A,4X,A,6X,A,4X,A,4X,A)') "n", "l", "k", "energy [Ha]", "energy [eV]", "outer radius [A]"
       call localorb_info(info_str,use_unit,'(A)')

       do i_fn = 1, n_atomic(i_species), 1
         write(info_str,'(2X,I3,2X,I3,2X,I3,2X,F15.6,F15.4,F15.6)') &
           atomic_n(i_species, i_fn), atomic_l(i_species, i_fn), atomic_k(i_species, i_fn), &
           atomic_eigenval(i_species, i_fn), atomic_eigenval(i_species, i_fn)*hartree, &
           atomic_outer_radius(i_species,i_fn)*bohr
           call localorb_info(info_str,use_unit,'(A)')
       enddo
       call localorb_info('',use_unit)
     else
       write(info_str,'(4X,A,4X,A,6X,A,4X,A,4X,A)') "n", "l", "energy [Ha]", "energy [eV]", "outer radius [A]"
       call localorb_info(info_str,use_unit,'(A)')

       do i_fn = 1, n_atomic(i_species), 1
         write(info_str,'(2X,I3,2X,I3,2X,F15.6,F15.4,F15.6)') &
           atomic_n(i_species, i_fn), atomic_l(i_species, i_fn), &
           atomic_eigenval(i_species, i_fn), atomic_eigenval(i_species, i_fn)*hartree, &
           atomic_outer_radius(i_species,i_fn)*bohr
           call localorb_info(info_str,use_unit,'(A)')
       enddo
       call localorb_info('',use_unit)
     endif

   end if
   enddo  ! end loop of i_species

!  that's all folks

   return
 end subroutine get_species_basis_fns
!******
