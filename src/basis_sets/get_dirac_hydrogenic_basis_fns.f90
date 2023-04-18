
!----------------------------------------------------------------------
! get_dirac_hydrogenic_basis_fns() constructs hydro_wave_large() and hydro_kinetic() 
! from the analytical solutions of an artificial hydrogen-like atom.
! This subroutine has currently not been debugged. I do not use hydro basis
! without cutoff, currently.
! -- Rundong Zhao May 16, 2018
!
!----------------------------------------------------------------------
 subroutine get_dirac_hydrogenic_basis_fns ( i_species )
!  PURPOSE
!    Subroutine get_dirac_hydrogenic_basis_fns tabulates Dirac hydrogen-like
!    wave functions as polarisation functions for one species
!    For now for debugging purposes, only does the non-relativistic 
!    wave functions but uses relativistic infrastructure
!  USES
      use constants,    only : bohr, hartree
      use dimensions,   only : n_max_ind_fns, use_basis_gradients, use_ext_basis
      use grids,        only : n_grid, r_grid
      use species_data, only : n_hydro, hydro_wave_large, hydro_wave_small, &
                               hydro_wave_deriv, hydro_kinetic, hydro_in_large_basis, &
                               hydro_n, hydro_l, hydro_kappa, &
                               hydro_scale, hydro_outer_radius 
      use mpi_tasks,    only : myid, aims_stop_coll
      use localorb_io,  only : use_unit
      implicit none

      integer :: i_species

      integer :: l_shell, n_shell
      real*8, dimension(:), allocatable :: hermite_coeff
      real*8 :: radius, term, poly, deriv_term, poly_deriv
      real*8 :: z_eff(n_max_ind_fns)
      real*8 :: r_outer_max(n_max_ind_fns)
      real*8 :: r_inner_max(n_max_ind_fns)
      real*8 :: norm
      real*8 :: eigenval(n_max_ind_fns)

!  counters

      integer i_grid, i_dirac_hydro_upper

      integer i_l, i_kappa, i_shell, i_k

!  functions

      character l_to_str
      real*8 int_log_mesh

!  begin work

      if ( any(n_hydro > 0) ) then
        if (myid.eq.0) then
          write(use_unit,'(A)') "* You are attempting to use analytic forms for the Dirac hydrogenic polarization basis "
          write(use_unit,'(A)') "* functions (i.e. you have not specified a cut-off potential.)  These have not yet been"
          write(use_unit,'(A)') "* implemented.  Exiting."
        end if
        call aims_stop_coll("","get_dirac_hydrogenic_basis_fns")
      end if

      ! This is old debug code from the non-relativistic hydro, where I
      ! copy/pasted everything but used the Dirac upper component framework.

!     Tabulate Hermite polynomial times exp decay, shell by shell
      do i_dirac_hydro_upper = 1, n_hydro(i_species), 1

        l_shell = hydro_l(i_species, i_dirac_hydro_upper)
        n_shell = hydro_n(i_species, i_dirac_hydro_upper)

        if (allocated(hermite_coeff)) then
          deallocate(hermite_coeff)
        end if

!       allocate one too many, to avoid zero allocation
        allocate (hermite_coeff(n_shell-l_shell))

!       tabulate polynomial coefficients
        do i_k = 1, (n_shell-l_shell-1), 1
          hermite_coeff(i_k) = &
          - (2.d0 * (n_shell-l_shell-i_k)) / &
          (i_k * (i_k+2.d0*l_shell+1.d0) * n_shell)
        enddo

        z_eff(i_dirac_hydro_upper) = hydro_scale(i_species, i_dirac_hydro_upper)

!       tabulate wave fn for z = z_eff
        do i_grid = 1, n_grid(i_species), 1
          radius = r_grid(i_grid,i_species)*z_eff(i_dirac_hydro_upper)
!         calculate polynomial contribution
          term = 1.0d0
          poly = 1.0d0
          do i_k = 1, (n_shell-l_shell-1), 1
            term = term * radius * hermite_coeff(i_k)
            poly = poly+term
          enddo
!         calculate radial function at r_grid
          hydro_wave_large(i_grid, i_species, i_dirac_hydro_upper) = &
            exp(-radius/n_shell) * (radius**(l_shell+1)) * poly
        enddo

!       normalize wave function using int_log_grid ...
        norm = &
          int_log_mesh ( &
          hydro_wave_large(1,i_species,i_dirac_hydro_upper), &
          hydro_wave_large(1,i_species,i_dirac_hydro_upper), &
          n_grid(i_species), r_grid(1,i_species) &
                       )

        norm = sqrt(norm)

        do i_grid = 1, n_grid(i_species), 1
          hydro_wave_large(i_grid, i_species, i_dirac_hydro_upper) = &
          hydro_wave_large(i_grid, i_species, i_dirac_hydro_upper)/norm
        enddo

        if (use_basis_gradients) then
!         tabulate wave fn derivative for z = z_eff
          do i_grid = 1, n_grid(i_species), 1

            radius = r_grid(i_grid,i_species)*z_eff(i_dirac_hydro_upper)

!           calculate polynomial contribution
            term = 1.0d0
            poly = 1.0d0
            poly_deriv = 0.d0
            do i_k = 1, (n_shell-l_shell-1), 1
              deriv_term = term * hermite_coeff(i_k) * i_k
              term = term * radius * hermite_coeff(i_k)
              poly = poly+term
              poly_deriv = poly_deriv + deriv_term
            enddo

!           calculate radial derivative at r_grid
!           safeguard against l_shell = 0 ... some compilers might not like it!
            if (l_shell.eq.0) then
              hydro_wave_deriv (i_grid, i_species, i_dirac_hydro_upper) = &
                exp(-radius/n_shell) * &
                ( - radius * poly / n_shell &
                  + poly &
                  + radius * poly_deriv &
                )
            else
              hydro_wave_deriv (i_grid, i_species, i_dirac_hydro_upper) = &
                exp(-radius/n_shell) * &
                ( - (radius**(l_shell+1)) * poly / n_shell &
                  + (l_shell+1) * (radius**(l_shell)) * poly &
                  + (radius**(l_shell+1)) * poly_deriv &
                )
            end if

!           account for the transformation r -> zr
            hydro_wave_deriv (i_grid, i_species, i_dirac_hydro_upper) = &
              hydro_wave_deriv (i_grid, i_species, i_dirac_hydro_upper) &
            * z_eff(i_dirac_hydro_upper)

!           normalize derivative
            hydro_wave_deriv(i_grid,i_species,i_dirac_hydro_upper) = &
            hydro_wave_deriv(i_grid,i_species,i_dirac_hydro_upper)/norm

          enddo
        end if

!       find innermost maximum

        i_grid = 1

        do while &
          (abs(hydro_wave_large(  i_grid+1, i_species, i_dirac_hydro_upper)).ge. &
           abs(hydro_wave_large(i_grid, i_species, i_dirac_hydro_upper))  )
           i_grid = i_grid+1
           if (i_grid.eq.n_grid(i_species)) then
              if(myid.eq.0) then
                 write(use_unit,*) "* Dirac hydrogenic (upper component) wave function no. ", i_dirac_hydro_upper, ":"
                 write(use_unit,*) "* No innermost maximum - bad log. grid?"
              end if
              stop
           end if
        enddo

        r_inner_max(i_dirac_hydro_upper) = r_grid(i_grid, i_species)

!       find outermost maximum

        i_grid = n_grid(i_species)

        do while &
          (abs(hydro_wave_large(  i_grid, i_species, i_dirac_hydro_upper)).le. &
           abs(hydro_wave_large(i_grid-1, i_species, i_dirac_hydro_upper))  )
           i_grid = i_grid-1
           if (i_grid.eq.1) then
              if(myid.eq.0) then
                 write(use_unit,*) "* Dirac hydrogenic (upper component) wave function no. ", i_dirac_hydro_upper, ":"
                 write(use_unit,*) "* No outermost maximum - bad log. grid?"
              end if
              stop
           end if
        enddo

        r_outer_max(i_dirac_hydro_upper) = r_grid(i_grid, i_species)

!       now obtain kinetic energy term for same wave function ...
!       Potential v(r)=-(z_eff/r)
!       H-like eigenvalue e=-0.5*[(z_eff)^2/n^2]
!         note our convention that Ha = 1. i.e. Ry = 0.5 !!
!       We store [e-v(r)]*u(r) to calculate kinetic energy later.

        eigenval(i_dirac_hydro_upper) = - 0.5d0 * z_eff(i_dirac_hydro_upper)**2./n_shell**2.

        do i_grid = 1, n_grid(i_species), 1

          hydro_kinetic(i_grid, i_species, i_dirac_hydro_upper) = eigenval(i_dirac_hydro_upper) + z_eff(i_dirac_hydro_upper) / r_grid(i_grid,i_species)

          hydro_kinetic(i_grid, i_species, i_dirac_hydro_upper) = hydro_kinetic(i_grid, i_species, i_dirac_hydro_upper) * hydro_wave_large(i_grid, i_species, i_dirac_hydro_upper)

        enddo

      enddo

!     output Dirac hydrogenic basis data

      if(myid.eq.0 .and. (.not. all(hydro_in_large_basis(i_species, 1:n_hydro(i_species))))) then
        write(use_unit,*)
        write(use_unit,*) " List of Dirac hydrogenic (upper component) basis orbitals: "
        write(use_unit,'(4X,A,4X,A,4X,A,6X,A,6X,A,2X,A,2X,A)') "n", "l", "kappa", "effective z", "eigenvalue [eV]", "inner max. [bohr]", "outer max. [bohr]"
      end if

      do i_dirac_hydro_upper = 1, n_hydro(i_species), 1
        i_l     = hydro_l(i_species, i_dirac_hydro_upper)
        i_kappa = hydro_kappa(i_species, i_dirac_hydro_upper)
        i_shell = hydro_n(i_species, i_dirac_hydro_upper)

        if(myid.eq.0) then
          if( hydro_in_large_basis(i_species,i_dirac_hydro_upper)) cycle
          
          write(use_unit,'(2X,I3,2X,I3,2X,I3,2X,F15.6,F15.4,4X,F15.6,4X,F15.6)') i_shell, i_l, i_kappa, &
            z_eff(i_dirac_hydro_upper), eigenval(i_dirac_hydro_upper)*hartree, r_inner_max(i_dirac_hydro_upper),r_outer_max(i_dirac_hydro_upper)
        end if

      enddo

      if(use_ext_basis .and. any(hydro_in_large_basis(i_species, 1:n_hydro(i_species)))) then
        if(myid.eq.0) then
          write(use_unit,*)
          write(use_unit,*)" List of extra Dirac hydrogenic (upper component) orbitals for auxiliary basis: "
          write(use_unit,'(4X,A,4X,A,4X,A,6X,A,6X,A,2X,A,2X,A,2X,A)')"n", "l", "kappa", "effective z", "eigenvalue [eV]", &
            "inner max. [A]   ", "outer max. [A]   ", "outer radius [A]   "
        end if

        do i_dirac_hydro_upper = 1, n_hydro(i_species), 1
          if(.not. hydro_in_large_basis(i_species,i_dirac_hydro_upper)) cycle
          i_l     = hydro_l(i_species, i_dirac_hydro_upper)
          i_kappa = hydro_kappa(i_species, i_dirac_hydro_upper)
          i_shell = hydro_n(i_species, i_dirac_hydro_upper)

          if(myid.eq.0) then
            write(use_unit,'(2X,I3,2X,I3,2X,I3,2X,F15.6,F15.4,4X,F15.6,4X,F15.6,4X,F15.6)') &
              i_shell, i_l, i_kappa, z_eff(i_dirac_hydro_upper), eigenval(i_dirac_hydro_upper)*hartree, &
              r_inner_max(i_dirac_hydro_upper)*bohr,r_outer_max(i_dirac_hydro_upper)*bohr, &
              hydro_outer_radius(i_species,i_dirac_hydro_upper)*bohr
          end if

        enddo
      end if

      if(myid.eq.0) then
         write(use_unit,*)
      end if

 end subroutine get_dirac_hydrogenic_basis_fns

!-------------------------------------------------------------------------------------------------------
! int_dirac_hydrogenic_basis_fns constructs hydro_wave_large(), hydro_wave_small() and hydro_kinetic() 
! from the numerical solutions of an artificial hydrogen-like atom.
! The basic framework is similar to that of get_species_basis_fns which generates the minimal basis set.
! -- Rundong Zhao May 16, 2018
!-------------------------------------------------------------------------------------------------------
 subroutine int_dirac_hydrogenic_basis_fns ( i_species )
 !  Generate hydrogen-like basis for relativistic cases.
 !  Similar to get_ionic_basis_fns and get_conf_basis_fns.
   use constants,       only : hartree, bohr, light_speed, light_speed_sq
   use dimensions,      only : n_max_grid, n_max_ind_fns, use_ext_basis
   use runtime_choices, only : wave_threshold, flag_rel, REL_x2c, REL_4c_dks
   use grids,           only : n_grid, r_grid, r_grid_inc
   use species_data,    only : hydro_wave_large, hydro_wave_small, hydro_large_deriv, hydro_small_deriv, &
                               hydro_kinetic, hydro_kinetic_small, n_hydro, hydro_n, hydro_l, hydro_kappa, hydro_scale, &
                               hydro_in_large_basis, hydro_cutoff, hydro_outer_radius, r_cutoff, w_cutoff, &
                               cutoff_type, basis_dep_cutoff_thresh, scale_cutoff
   use mpi_tasks,       only : myid
   use localorb_io,     only : use_unit
   use reigen_dftatom
   implicit none
  
   integer,intent(in) :: i_species
  
   integer l_shell, kappa, n_shell, relat
   real*8 z_eff(n_hydro(i_species))   ! the artifical nuclear charge Z
   real*8 basis_pot(n_max_grid)  ! radial potential
   real*8 eigenval(n_hydro(i_species))
   real*8 wave_deriv1 (n_max_grid)
   real*8 wave_deriv2 (n_max_grid)
   real*8 r_outer_max(n_hydro(i_species))
   real*8 r_inner_max(n_hydro(i_species))
   real*8 r_outer
  
   integer converged_dftatom
   real*8 :: alpha_grid, tmp, tmp1
   real*8 :: int_outer_l, int_outer_s
   real*8 :: r_prime(n_max_grid)
  
   real*8, dimension(n_max_grid,n_hydro(i_species)) :: upper_comp
   real*8, dimension(n_max_grid,n_hydro(i_species)) :: upper_comp_deriv
   real*8, dimension(n_max_grid,n_hydro(i_species)) :: lower_comp
   real*8, dimension(n_max_grid,n_hydro(i_species)) :: lower_comp_deriv
   real*8, dimension(4,n_max_grid)     :: spline_coef
  
   ! counters
   integer i_grid, i_hydro
   integer i_l, i_kappa, i_shell, i,j,k
  
   ! functions
   real*8 cutoff_pot

   ! debug:
   real*8 :: upper(n_grid(i_species)),lower(n_grid(i_species))

  
   hydro_cutoff(i_species,:) = r_cutoff(i_species)
  
   ! Generate grid derivatives:
   alpha_grid = log(r_grid_inc(i_species))
   do i_grid=1, n_grid(i_species)
     r_prime(i_grid) = alpha_grid * r_grid(i_grid,i_species)
   enddo
  
   do i_hydro = 1, n_hydro(i_species), 1
  
     n_shell = hydro_n(i_species, i_hydro)
     l_shell = hydro_l(i_species, i_hydro)
     kappa   = hydro_kappa(i_species, i_hydro)
     if( kappa .eq. -l_shell-1 )then
       relat=2  ! spin up
     elseif( kappa .eq. l_shell )then
       relat=3  ! spin down
     else
       write(use_unit,*)'illegal kappa:',' get_dirac_hydrogenic_basis_fns'
       stop
     endif
  
     ! outer cutoff radius
     r_outer = r_cutoff(i_species) + w_cutoff(i_species)
  
     z_eff(i_hydro) = hydro_scale(i_species, i_hydro)
  
     ! generate radial potential with confinement potential added
     do i_grid = 1, n_grid(i_species), 1
       basis_pot(i_grid) = - z_eff(i_hydro)/r_grid(i_grid,i_species)
       basis_pot(i_grid) = basis_pot(i_grid) + cutoff_pot( r_grid(i_grid,i_species), &
         cutoff_type(i_species), r_cutoff(i_species), w_cutoff(i_species), scale_cutoff(i_species))
     enddo

     ! Solve Dirac equation, to get the large and small component radial wavefunctions
     call solve_radial_eigenproblem( n_grid(i_species), n_shell, l_shell, -10.d0, 1.d-10, 100, &
       r_grid(1,i_species), r_prime, basis_pot, z_eff(i_hydro), light_speed, relat, .true., &
       -1.d4, 10.d2, converged_dftatom, eigenval(i_hydro), upper_comp(1,i_hydro), lower_comp(1,i_hydro) )

   ! write(6,"('i_hydro:',i3,5x,'n:',i3,5x,'l:',i3,5x,'k:',i3,5x,'E:',f12.6)") i_hydro,n_shell,l_shell,kappa,eigenval(i_hydro)
   ! write(6,"(20f12.6)")upper_comp(:,i_hydro)

     ! Then calculate wavefunction derivatives
     call fderiv(1,n_grid(i_species),r_grid(1,i_species), upper_comp(1,i_hydro),upper_comp_deriv(1,i_hydro),spline_coef)
     call fderiv(1,n_grid(i_species),r_grid(1,i_species), lower_comp(1,i_hydro),lower_comp_deriv(1,i_hydro),spline_coef)

     ! Finally, update the basis function variables to the upper component values
     hydro_wave_large(:,i_species,i_hydro)  = upper_comp(:,i_hydro)
     hydro_wave_small(:,i_species,i_hydro)  = lower_comp(:,i_hydro)
     hydro_large_deriv(:,i_species,i_hydro) = upper_comp_deriv(:,i_hydro)
     hydro_small_deriv(:,i_species,i_hydro) = lower_comp_deriv(:,i_hydro)
     
     ! Calculate new onset of cutoff potential if desired. Rerun the integrator with new potential.
     if (basis_dep_cutoff_thresh(i_species).gt.0d0) then
       ! determine the shortest cutoff potential required by the threshold
       i_grid = n_grid(i_species)
  
       ! Integrate the wave funtion from the outside to calculate if/when 
       ! the set threshold for the cutoff potential is reached.
       int_outer_l = alpha_grid*r_grid(i_grid,i_species)*upper_comp(i_grid,i_hydro)**2d0
       int_outer_s = alpha_grid*r_grid(i_grid,i_species)*lower_comp(i_grid,i_hydro)**2d0
       do while (int_outer_l.lt.basis_dep_cutoff_thresh(i_species) .and. &
                int_outer_s.lt.basis_dep_cutoff_thresh(i_species))
         i_grid = i_grid - 1 
         int_outer_l = int_outer_l + alpha_grid*r_grid(i_grid,i_species)*upper_comp(i_grid,i_hydro)**2d0
         int_outer_s = int_outer_s + alpha_grid*r_grid(i_grid,i_species)*lower_comp(i_grid,i_hydro)**2d0
       end do
       i_grid = i_grid + 1 ! the index of the last step does NOT satisfy the threshold, use the former one
       
       hydro_cutoff(i_species,i_hydro) = min(r_grid(i_grid,i_species),r_cutoff(i_species))
     
       ! recalculate basis function if necessary - same procedure as above
       if (hydro_cutoff(i_species,i_hydro).lt.r_cutoff(i_species)) then
         ! retabulate wave fn for z = z_eff
         do i_grid = 1, n_grid(i_species)
            basis_pot(i_grid) = - z_eff(i_hydro)/r_grid(i_grid,i_species)
            basis_pot(i_grid) = basis_pot(i_grid) + &
              cutoff_pot ( r_grid(i_grid,i_species), cutoff_type(i_species), &
              hydro_cutoff(i_species,i_hydro), w_cutoff(i_species), scale_cutoff(i_species) )
         enddo              
  
         r_outer = hydro_cutoff(i_species,i_hydro) + w_cutoff(i_species)
       
         call solve_radial_eigenproblem( n_grid(i_species), n_shell, l_shell, -10.d0, 1.d-10, 100, &
           r_grid(1,i_species), r_prime, basis_pot, z_eff(i_hydro), light_speed, relat, .true., &
           -1.d4, 10.d2, converged_dftatom, eigenval(i_hydro), upper_comp(1,i_hydro), lower_comp(1,i_hydro) )

         ! Then calculate wavefunction derivatives
         call fderiv(1,n_grid(i_species),r_grid(1,i_species), upper_comp(1,i_hydro),upper_comp_deriv(1,i_hydro),spline_coef)
         call fderiv(1,n_grid(i_species),r_grid(1,i_species), lower_comp(1,i_hydro),lower_comp_deriv(1,i_hydro),spline_coef)

         ! Finally, update the basis function variables to the upper component values
         hydro_wave_large(:,i_species,i_hydro)  = upper_comp(:,i_hydro)
         hydro_wave_small(:,i_species,i_hydro)  = lower_comp(:,i_hydro)
         hydro_large_deriv(:,i_species,i_hydro) = upper_comp_deriv(:,i_hydro)
         hydro_small_deriv(:,i_species,i_hydro) = lower_comp_deriv(:,i_hydro)
       end if
  
     end if   ! end basis function dependent cutoff potential
  
     ! Now obtain kinetic energy term.
     if(flag_rel.eq.REL_x2c)then
       do i_grid = 1, n_grid(i_species)
         tmp = eigenval(i_hydro) - basis_pot(i_grid)
         tmp1 = tmp/light_speed
         hydro_kinetic(i_grid, i_species, i_hydro) = ( tmp + 0.5d0*tmp1*tmp1 ) * hydro_wave_large(i_grid,i_species,i_hydro)
         hydro_kinetic_small(i_grid, i_species, i_hydro) = eigenval(i_hydro) - basis_pot(i_grid) + 2.d0*light_speed_sq
         hydro_kinetic_small(i_grid, i_species, i_hydro) = hydro_kinetic_small(i_grid, i_species, i_hydro) * hydro_wave_small(i_grid, i_species, i_hydro)
       enddo
     elseif (flag_rel.eq.REL_4c_dks)then ! Q4C for NAO
       hydro_kinetic(:,i_species,i_hydro) = eigenval(i_hydro) - basis_pot(:)
       hydro_kinetic(:,i_species,i_hydro) = hydro_kinetic(:,i_species,i_hydro) * hydro_wave_large(:,i_species,i_hydro)
       if(.true.)then ! use atomic balance condition
       hydro_kinetic_small(:,i_species,i_hydro) = eigenval(i_hydro) - basis_pot(:) + 2.d0*light_speed_sq
       hydro_kinetic_small(:,i_species,i_hydro) = hydro_kinetic_small(:,i_species,i_hydro) * hydro_wave_small(:,i_species,i_hydro)
     ! else ! use restricted kinetic balance condition
     ! !(Rundong) in this way, the kinetic energy term is actually in a nonrelativistic form, viz. the second derivative of the basis.
     !!! The precision is low. Thus, I don't recommend this strategy.
     ! call sigma_dot_p_large_wave( n_grid(i_species), r_grid(1,i_species), hydro_l(i_species,i_hydro), hydro_kappa(i_species,i_hydro),&
     !      hydro_wave_large(1,i_species,i_hydro), hydro_large_deriv(1,i_species,i_hydro), hydro_wave_small(1,i_species,i_hydro), &
     !      hydro_kinetic_small(1,i_species,i_hydro) )
     ! call fderiv(1,n_grid(i_species),r_grid(1,i_species), hydro_wave_small(:,i_species,i_hydro),hydro_small_deriv(:,i_species,i_hydro),spline_coef)
     !!call fderiv(2,n_grid(i_species),r_grid(1,i_species), upper_comp(1,i_hydro),hydro_kinetic(1,i_species,i_hydro),spline_coef)

     !! if use zora kind of atomic balance
     ! hydro_wave_small(:,i_species,i_hydro) = 2*light_speed_sq/(2*light_speed_sq-basis_pot(:)) * hydro_wave_small(:,i_species,i_hydro)
     ! hydro_kinetic(:,i_species,i_hydro) = -basis_pot(:) * hydro_wave_large(:,i_species,i_hydro)
     ! hydro_kinetic_small(:,i_species,i_hydro) = (2.d0*light_speed_sq-basis_pot(:)) * hydro_wave_small(:,i_species,i_hydro)
       endif
     endif
  
     !test - create new outermost radius array after the fact, for each radial function
       i_grid = n_grid(i_species)
       do while (abs(hydro_wave_large(i_grid,i_species,i_hydro)).lt.wave_threshold .and. &
                 abs(hydro_wave_small(i_grid,i_species,i_hydro)).lt.wave_threshold)
         i_grid = i_grid - 1 
       end do
       hydro_outer_radius(i_species,i_hydro) = r_grid(i_grid,i_species)
     !end test
  
     ! find innermost maximum
     i_grid = 1
  
     do while (abs(hydro_wave_large( i_grid+1, i_species, i_hydro)) .ge. &
               abs(hydro_wave_large(i_grid, i_species, i_hydro)) )
        i_grid = i_grid+1
        if (i_grid.eq.n_grid(i_species)) then
           if(myid.eq.0) then
              write(use_unit,*) "* Dirac hydrogenic (upper component) wave function no. ", i_hydro, ":"
              write(use_unit,*) "* No innermost maximum - bad log. grid?"
           end if
          stop
        end if
     enddo
  
     r_inner_max(i_hydro) = r_grid(i_grid, i_species)
  
     ! find outermost maximum
  
     i_grid = n_grid(i_species)
  
     do while &
       (abs(hydro_wave_large(  i_grid, i_species, i_hydro)).le. &
        abs(hydro_wave_large(i_grid-1, i_species, i_hydro))  )
        i_grid = i_grid-1
        if (i_grid.eq.1) then
           if(myid.eq.0) then
              write(use_unit,*) "* Dirac hydrogenic (upper component) wave function no. ", i_hydro, ":"
              write(use_unit,*) "* No outermost maximum - bad log. grid?"
           end if
           stop
        end if
     enddo
  
     r_outer_max(i_hydro) = r_grid(i_grid, i_species)
  
   enddo ! end of i_hydro
  

   ! output Dirac hydrogenic basis data
  
   if(myid.eq.0 .and. (.not. all(hydro_in_large_basis(i_species, 1:n_hydro(i_species))))) then
     write(use_unit,*)
     write(use_unit,*)" List of Dirac hydrogenic (upper component) basis orbitals: "
     write(use_unit,'(4X,A,4X,A,4X,A,6X,A,6X,A,2X,A,2X,A,2X,A)') &
       "n", "l", "k", "effective z", "eigenvalue [eV]", &
       "inner max. [A]   ", "outer max. [A]   ","outer radius [A]"
   end if
  
   do i_hydro = 1, n_hydro(i_species), 1
     if( hydro_in_large_basis(i_species,i_hydro)) cycle
     i_l     = hydro_l(i_species, i_hydro)
     i_kappa = hydro_kappa(i_species, i_hydro)
     i_shell = hydro_n(i_species, i_hydro)
  
     if(myid.eq.0) then
       write(use_unit,'(2X,I3,2X,I3,2X,I3,2X,F15.6,F15.4,4X,F15.6,4X,F15.6,4X,F15.6)') &
         i_shell, i_l, i_kappa, z_eff(i_hydro), eigenval(i_hydro)*hartree, r_inner_max(i_hydro)*bohr, &
         r_outer_max(i_hydro)*bohr, hydro_outer_radius(i_species,i_hydro)*bohr
     end if
  
   enddo
   ! output extra basis function if present
   if(use_ext_basis .and. &
     any(hydro_in_large_basis(i_species, 1:n_hydro(i_species)))) then
     if(myid.eq.0) then
       write(use_unit,*)
       write(use_unit,*)" List of extra Dirac hydrogenic (upper component) orbitals for auxiliary basis:"
       write(use_unit,'(4X,A,4X,A,4X,A,6X,A,6X,A,2X,A,2X,A,2X,A)') "n", "l", "k", &
         "effective z", "eigenvalue [eV]", "inner max. [A]   ", "outer max. [A]   ","outer radius [A]"
     end if
  
     do i_hydro = 1, n_hydro(i_species), 1
       if(.not. hydro_in_large_basis(i_species,i_hydro)) cycle
       i_l     = hydro_l(i_species, i_hydro)
       i_kappa = hydro_kappa(i_species, i_hydro)
       i_shell = hydro_n(i_species, i_hydro)
  
       if(myid.eq.0) then
         write(use_unit,'(2X,I3,2X,I3,2X,I3,2X,F15.6,F15.4,4X,F15.6,4X,F15.6,4X,F15.6)') &
           i_shell, i_l, i_kappa, z_eff(i_hydro), eigenval(i_hydro)*hartree, r_inner_max(i_hydro)*bohr,&
           r_outer_max(i_hydro)*bohr, hydro_outer_radius(i_species,i_hydro)*bohr
       end if
  
     enddo
   end if
  
   if(myid.eq.0) then
      write(use_unit,*)
   end if

 end subroutine int_dirac_hydrogenic_basis_fns






