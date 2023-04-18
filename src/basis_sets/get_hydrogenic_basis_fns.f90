!----------------------------------------------------------------------
!
!  Provides:
!  * get_hydrogenic_basis_fns():
!    constructs hydro_wave() and hydro_kinetic() from analytic results for H atom
!  * int_hydrogenic_basis_fns()
!    integrates hydro_wave() and hydro_kinetic() numerically, for given Coulomb
!    potential
!
!----------------------------------------------------------------------

!****s* FHI-aims/get_hydrogenic_basis_fns
!  NAME
!   get_hydrogenic_basis_fns
!  SYNOPSIS
  subroutine get_hydrogenic_basis_fns ( i_species )
!  PURPOSE
!    Subroutine get_hydrogenic_basis_fns tabulates non-relativistic hydrogen-like
!    wave functions as polarisation functions for one species
!  USES
      use constants,       only : bohr, hartree
      use dimensions,      only : n_max_ind_fns, use_basis_gradients, use_ext_basis
      use grids,           only : n_grid, r_grid
      use species_data,    only : hydro_wave, hydro_wave_deriv, hydro_kinetic, &
                                  n_hydro, hydro_in_large_basis, hydro_l, hydro_n, &
                                  hydro_scale, hydro_outer_radius
      use mpi_tasks,       only : myid
      use localorb_io,     only : use_unit
      use psi_at_nucleus_mod, only: psi_at_nucleus_hydro
      implicit none
!  ARGUMENTS
      integer :: i_species
! INPUTS 
!   o i_species -- species index
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
      integer :: l_shell, n_shell
      real*8, dimension(:), allocatable :: hermite_coeff
      real*8 :: radius, term, poly, deriv_term, poly_deriv
      real*8 :: z_eff(n_max_ind_fns)
      real*8 :: r_outer_max(n_max_ind_fns)
      real*8 :: r_inner_max(n_max_ind_fns)
      real*8 :: norm
      real*8 :: eigenval(n_max_ind_fns)

!  counters

      integer i_grid, i_hydro

      integer i_l, i_shell, i_k

!  functions

      character l_to_str
      real*8 int_log_mesh

!  begin work
  
!     Tabulate Hermite polynomial times exp decay, shell by shell
      do i_hydro = 1, n_hydro(i_species), 1

        l_shell = hydro_l(i_species, i_hydro)
        n_shell = hydro_n(i_species, i_hydro)

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

        z_eff(i_hydro) = hydro_scale(i_species, i_hydro)

!       tabulate wave fn for z = z_eff
        do i_grid = 1, n_grid(i_species), 1
          radius = r_grid(i_grid,i_species)*z_eff(i_hydro)
!         calculate polynomial contribution
          term = 1.0d0
          poly = 1.0d0
          do i_k = 1, (n_shell-l_shell-1), 1
            term = term * radius * hermite_coeff(i_k)
            poly = poly+term
          enddo
!         calculate radial function at r_grid
          hydro_wave(i_grid, i_species, i_hydro) = &
            exp(-radius/n_shell) * (radius**(l_shell+1)) * poly
        enddo

!       normalize wave function using int_log_grid ...
        norm = &
          int_log_mesh ( &
          hydro_wave(1,i_species,i_hydro), &
          hydro_wave(1,i_species,i_hydro), &
          n_grid(i_species), r_grid(1,i_species) &
                       )

        norm = sqrt(norm)

        do i_grid = 1, n_grid(i_species), 1
          hydro_wave(i_grid, i_species, i_hydro) = &
          hydro_wave(i_grid, i_species, i_hydro)/norm
        enddo

        if (use_basis_gradients) then
!         tabulate wave fn derivative for z = z_eff
          do i_grid = 1, n_grid(i_species), 1

            radius = r_grid(i_grid,i_species)*z_eff(i_hydro)

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
              hydro_wave_deriv (i_grid, i_species, i_hydro) = &
                exp(-radius/n_shell) * &
                ( - radius * poly / n_shell &
                  + poly &
                  + radius * poly_deriv &
                )
            else
              hydro_wave_deriv (i_grid, i_species, i_hydro) = &
                exp(-radius/n_shell) * &
                ( - (radius**(l_shell+1)) * poly / n_shell &
                  + (l_shell+1) * (radius**(l_shell)) * poly &
                  + (radius**(l_shell+1)) * poly_deriv &
                )
            end if

!           account for the transformation r -> zr
            hydro_wave_deriv (i_grid, i_species, i_hydro) = &
              hydro_wave_deriv (i_grid, i_species, i_hydro) &
            * z_eff(i_hydro)

!           normalize derivative
            hydro_wave_deriv(i_grid,i_species,i_hydro) = &
            hydro_wave_deriv(i_grid,i_species,i_hydro)/norm

          enddo
        end if

!       find innermost maximum

        i_grid = 1

        do while &
          (abs(hydro_wave(  i_grid+1, i_species, i_hydro)).ge. &
           abs(hydro_wave(i_grid, i_species, i_hydro))  )
           i_grid = i_grid+1
           if (i_grid.eq.n_grid(i_species)) then
              if(myid.eq.0) then
                 write(use_unit,*) "* Hydrogenic wave function no. ", &
                      i_hydro, ":"
                 write(use_unit,*) "* No innermost maximum - bad log. grid?"
              end if
              stop
           end if
        enddo

        r_inner_max(i_hydro) = r_grid(i_grid, i_species)

!       find outermost maximum

        i_grid = n_grid(i_species)

        do while &
          (abs(hydro_wave(  i_grid, i_species, i_hydro)).le. &
           abs(hydro_wave(i_grid-1, i_species, i_hydro))  )
           i_grid = i_grid-1
           if (i_grid.eq.1) then
              if(myid.eq.0) then
                 write(use_unit,*) "* Hydrogenic wave function no. ", &
                      i_hydro, ":"
                 write(use_unit,*) "* No outermost maximum - bad log. grid?"
              end if
              stop
           end if
        enddo

        r_outer_max(i_hydro) = r_grid(i_grid, i_species)

!       now obtain kinetic energy term for same wave function ...
!       Potential v(r)=-(z_eff/r)
!       H-like eigenvalue e=-0.5*[(z_eff)^2/n^2]
!         note our convention that Ha = 1. i.e. Ry = 0.5 !!
!       We store [e-v(r)]*u(r) to calculate kinetic energy later.

        eigenval(i_hydro) = - 0.5d0 * z_eff(i_hydro)**2./n_shell**2.

        do i_grid = 1, n_grid(i_species), 1

          hydro_kinetic(i_grid, i_species, i_hydro) = &
          eigenval(i_hydro) &
          + z_eff(i_hydro) / r_grid(i_grid,i_species)

          hydro_kinetic(i_grid, i_species, i_hydro) = &
          hydro_kinetic(i_grid, i_species, i_hydro) * &
          hydro_wave(i_grid, i_species, i_hydro)

        enddo

        ! Only the s orbitals may have nonzero values at the nucleus.
        if (l_shell == 0) psi_at_nucleus_hydro(i_hydro, i_species) = &
             & 2*(z_eff(i_hydro)/n_shell)**(1.5d0)
      enddo

!     output hydrogenic basis data

      if(myid.eq.0 .and. (.not. all(hydro_in_large_basis(i_species, 1:n_hydro(i_species))))) then
         write(use_unit,*)
         write(use_unit,*) " List of hydrogenic basis orbitals: "
         write(use_unit,'(4X,A,4X,A,6X,A,6X,A,2X,A,2X,A)') "n", "l", &
              "effective z", "eigenvalue [eV]", "inner max. [bohr]", &
              "outer max. [bohr]"
      end if

      do i_hydro = 1, n_hydro(i_species), 1
         i_l     = hydro_l(i_species, i_hydro)
         i_shell = hydro_n(i_species, i_hydro)

         if(myid.eq.0) then
            if( hydro_in_large_basis(i_species,i_hydro)) cycle
            
            write(use_unit,'(2X,I3,2X,I3,2X,F15.6,F15.4,4X,F15.6,4X,F15.6)') &
                 i_shell, i_l, &
                 z_eff(i_hydro), eigenval(i_hydro)*hartree, &
                 r_inner_max(i_hydro),r_outer_max(i_hydro)
         end if

      enddo

      if(use_ext_basis .and. any(hydro_in_large_basis(i_species, 1:n_hydro(i_species)))) then
        if(myid.eq.0) then
            write(use_unit,*)
            write(use_unit,*) " List of extra hydrogenic orbitals for auxiliary basis: "
            write(use_unit,'(4X,A,4X,A,6X,A,6X,A,2X,A,2X,A,2X,A)') "n", "l", &
              "effective z", "eigenvalue [eV]", "inner max. [A]   ", &
              "outer max. [A]   ","outer radius [A]   "
        end if

        do i_hydro = 1, n_hydro(i_species), 1
            if(.not. hydro_in_large_basis(i_species,i_hydro)) cycle
            i_l     = hydro_l(i_species, i_hydro)
            i_shell = hydro_n(i_species, i_hydro)

            if(myid.eq.0) then
                write(use_unit,'(2X,I3,2X,I3,2X,F15.6,F15.4,4X,F15.6,4X,F15.6,4X,F15.6)') &
                    i_shell, i_l, &
                    z_eff(i_hydro), eigenval(i_hydro)*hartree, &
                    r_inner_max(i_hydro)*bohr,r_outer_max(i_hydro)*bohr, &
                    hydro_outer_radius(i_species,i_hydro)*bohr
            end if

        enddo
      end if

      if(myid.eq.0) then
         write(use_unit,*)
      end if

!  that's all folks

      return
  end subroutine get_hydrogenic_basis_fns
!******
!---------------------------------------------------------------------------------------------------

!****s* FHI-aims/int_hydrogenic_basis_fns
!  NAME
!   int_hydrogenic_basis_fns
!  SYNOPSIS
      subroutine int_hydrogenic_basis_fns ( i_species )
! PURPOSE
!  Subroutine int_hydrogenic_basis_fns tabulates non-relativistic hydrogen-like
!  wave functions as polarisation functions for one species
!  We integrate the radial equation numerically, in full analogy to get_ionic_basis_fns
!  and get_conf_basis_fns
! USES
      use constants,       only : bohr, hartree
      use dimensions,      only : n_max_grid, n_max_ind_fns, use_basis_gradients, use_ext_basis
      use runtime_choices, only : wave_threshold
      use grids,           only : n_grid, r_grid, r_grid_inc
      use species_data,    only : hydro_cutoff, hydro_kinetic, hydro_wave_deriv, &
                                  hydro_outer_radius, n_hydro, hydro_in_large_basis, r_cutoff, &
                                  hydro_l, hydro_n, w_cutoff, hydro_scale, cutoff_type, &
                                  hydro_wave, scale_cutoff, basis_dep_cutoff_thresh
      use mpi_tasks,       only : aims_stop, myid
      use localorb_io,     only : use_unit
      use psi_at_nucleus_mod, only: psi_at_nucleus_hydro
      implicit none
!  ARGUMENTS
      integer :: i_species
! INPUTS 
!   o i_species -- species number in question
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
      integer l_shell, n_shell
      real*8 z_eff(n_max_ind_fns)
      real*8 r_outer_max(n_max_ind_fns)
      real*8 r_inner_max(n_max_ind_fns)

!     dftseq input
      real*8 basis_pot(n_max_grid)
      integer i_mode
      real*8 el_mass
      real*8 dftseq_v1(n_max_grid)
      real*8 dftseq_v2(n_max_grid)
      real*8 log_deriv
      real*8 r_outer

!     dftseq output
      integer n_grid_max
      integer n_grid_turn
      real*8 eigenval(n_max_ind_fns)
      real*8 wave_deriv1 (n_max_grid)
      real*8 wave_deriv2 (n_max_grid)

!     other

      real*8 :: alpha_grid
      real*8 :: int_outer
      ! Value of the analytical wavefunction at the first grid
      ! point. The analytical function is computed without the cutoff
      ! potential.
      real*8 :: wave_an_1

!  counters

      integer i_grid, i_hydro

      integer i_l, i_shell

!  functions

      character l_to_str
      real*8 cutoff_pot
      real*8 :: laguerre_l0

!  begin work

      hydro_cutoff(i_species,:) = r_cutoff(i_species)

!     initialize general input for dftseq()

!     non-relativistic version
      i_mode=2

!     other stuff
      el_mass = 1.0d0
      do i_grid = 1, n_max_grid, 1
        dftseq_v1(i_grid) = 0.0d0
        dftseq_v2(i_grid) = 1.0d0
      enddo
      log_deriv = 0.0d0

!     grid scale factor to amend derivatives
      alpha_grid = log(r_grid_inc(i_species))

      do i_hydro = 1, n_hydro(i_species), 1

        l_shell = hydro_l(i_species, i_hydro)
        n_shell = hydro_n(i_species, i_hydro)

!       outer cutoff radius
        r_outer = r_cutoff(i_species) + w_cutoff(i_species)

        z_eff(i_hydro) = hydro_scale(i_species, i_hydro)

!       retabulate wave fn for z = z_eff
        do i_grid = 1, n_grid(i_species), 1
          basis_pot(i_grid) = - z_eff(i_hydro)/r_grid(i_grid,i_species)
          basis_pot(i_grid) = basis_pot(i_grid) + &
            cutoff_pot &
            ( r_grid(i_grid,i_species),cutoff_type(i_species), &
              r_cutoff(i_species), w_cutoff(i_species), &
              scale_cutoff(i_species))
        enddo
!       initialize eigenvalue ...
        eigenval(i_hydro) = - (z_eff(i_hydro)**2.d0) / (n_shell**2.0d0)

        call dftseq &
        ( i_mode, z_eff(i_hydro), n_grid(i_species), &
          r_grid(1,i_species), &
          n_shell, l_shell, el_mass, basis_pot, dftseq_v1, &
          dftseq_v2, n_grid_max, n_grid_turn, log_deriv, &
          eigenval(i_hydro), hydro_wave(1,i_species,i_hydro), &
          wave_deriv1, wave_deriv2, r_outer &
        )
!       calculate new onset of cutoff potential if desired,
!           rerun integrator for wave function with new potential if necessary

        if (basis_dep_cutoff_thresh(i_species).gt.0d0) then
           ! determine the shortest cutoff potential required by the threshold
           i_grid = n_grid(i_species)

           ! integrate the wave funtion from the outside to calculate if/when the set threshold for the cutoff 
           ! potential is reached
           int_outer = alpha_grid*r_grid(i_grid,i_species)*hydro_wave(i_grid,i_species,i_hydro)**2d0
           do while (int_outer.lt.basis_dep_cutoff_thresh(i_species))
              i_grid = i_grid - 1 
              int_outer = int_outer + alpha_grid*r_grid(i_grid,i_species)*hydro_wave(i_grid,i_species,i_hydro)**2d0
           end do
           i_grid = i_grid + 1 ! the index determined in the last step does NOT satisfy the threshold, use one index higher
           
           hydro_cutoff(i_species,i_hydro) = min(r_grid(i_grid,i_species),r_cutoff(i_species))  ! set cutoff according to smaller r_cut
 
           ! recalculate basis function if necessary - same procedure as above
           if (hydro_cutoff(i_species,i_hydro).lt.r_cutoff(i_species)) then
              !       retabulate wave fn for z = z_eff
              do i_grid = 1, n_grid(i_species), 1
                 basis_pot(i_grid) = - z_eff(i_hydro)/r_grid(i_grid,i_species)
                 basis_pot(i_grid) = basis_pot(i_grid) + &
                      cutoff_pot &
                      ( r_grid(i_grid,i_species),cutoff_type(i_species), &
                      hydro_cutoff(i_species,i_hydro), w_cutoff(i_species), &
                      scale_cutoff(i_species))
              enddo              

              r_outer = hydro_cutoff(i_species,i_hydro) + w_cutoff(i_species)
              
              call dftseq &
                   ( i_mode, z_eff(i_hydro), n_grid(i_species), &
                   r_grid(1,i_species), &
                   n_shell, l_shell, el_mass, basis_pot, dftseq_v1, &
                   dftseq_v2, n_grid_max, n_grid_turn, log_deriv, &
                   eigenval(i_hydro), hydro_wave(1,i_species,i_hydro), &
                   wave_deriv1, wave_deriv2, r_outer &
                   )              
           end if

        end if   ! end basis function dependent cutoff potential

!       store non-relativistic kinetic energy radial function
!       hydrogenic wave functions are always generated using the
!       non-relativistic version of dftseq, i.e. no changes necessary if
!       relativistic overall calculation

!       now obtain kinetic energy term for same wave function ...
!       Potential v(r)=-(z_eff/r)
!       H-like eigenvalue e=-[(z_eff)^2/n^2]
!       We store [e-v(r)]*u(r) to calculate kinetic energy later.
        do i_grid = 1, n_grid(i_species),1
          hydro_kinetic(i_grid, i_species, i_hydro) = &
          eigenval(i_hydro) - basis_pot(i_grid)

          hydro_kinetic(i_grid, i_species, i_hydro) = &
          hydro_kinetic(i_grid, i_species, i_hydro) * &
          hydro_wave(i_grid, i_species, i_hydro)
        enddo

        if (use_basis_gradients) then
          do i_grid = 1, n_grid(i_species), 1
            hydro_wave_deriv(i_grid, i_species, i_hydro) = &
              wave_deriv1(i_grid) &
            / (alpha_grid*r_grid(i_grid,i_species))
          enddo
        end if

!test - create new outermost radius array after the fact, for each radial function
            i_grid = n_grid(i_species)
            do while (abs(hydro_wave(i_grid,i_species,i_hydro)).lt.wave_threshold)
               i_grid = i_grid - 1 
            end do
            hydro_outer_radius(i_species,i_hydro) = r_grid(i_grid,i_species)
!end test

!       find innermost maximum
        i_grid = 1

        do while &
          (abs(hydro_wave( i_grid+1, i_species, i_hydro)).ge. &
           abs(hydro_wave(i_grid, i_species, i_hydro))  )
           i_grid = i_grid+1
           if (i_grid.eq.n_grid(i_species)) then
              if(myid.eq.0) then
                 write(use_unit,*) "* Hydrogenic wave function no. ", &
                      i_hydro, ":"
                 write(use_unit,*) "* No innermost maximum - bad log. grid?"
              end if
             stop
           end if
        enddo

        r_inner_max(i_hydro) = r_grid(i_grid, i_species)

!       find outermost maximum

        i_grid = n_grid(i_species)

        do while &
          (abs(hydro_wave(  i_grid, i_species, i_hydro)).le. &
           abs(hydro_wave(i_grid-1, i_species, i_hydro))  )
           i_grid = i_grid-1
           if (i_grid.eq.1) then
              if(myid.eq.0) then
                 write(use_unit,*) "* Hydrogenic wave function no. ", &
                      i_hydro, ":"
                 write(use_unit,*) "* No outermost maximum - bad log. grid?"
              end if
              stop
           end if
        enddo

        r_outer_max(i_hydro) = r_grid(i_grid, i_species)

        ! Only the s orbitals may have nonzero values at the nucleus.
        if (l_shell == 0) then
           ! Analytically compute the uncorrected wavefunction value
           ! at the first grid point.
           if (n_shell >= 10) call aims_stop('Hydrogenic orbitals with &
                &n_shell>=10 not supported', 'int_hydrogenic_basis_fns::')
           wave_an_1 = 2d0/n_shell*(z_eff(i_hydro)/n_shell)**1.5d0* &
                & laguerre_l0(n_shell-1,2*z_eff(i_hydro)* &
                & r_grid(1,i_species)/n_shell)* &
                & exp(-z_eff(i_hydro)*r_grid(1,i_species)/n_shell)
           ! Next, analytically compute the uncorrected wavefunction
           ! value at the nucleus.
           psi_at_nucleus_hydro(i_hydro, i_species) = &
                & 2*(z_eff(i_hydro)/n_shell)**(1.5d0)
           ! Correct the previous value by taking into account the cutoff
           ! potential. The correction factor is Psi1/Psi_an1, where Psi1
           ! is the value at the first grid point as determined by dftseq
           ! above, and Psi_an1 is the analytical uncorrected value at
           ! the first grid point.
           psi_at_nucleus_hydro(i_hydro, i_species) = &
                & psi_at_nucleus_hydro(i_hydro, i_species) * &
                & hydro_wave(1,i_species,i_hydro)/r_grid(1,i_species)/wave_an_1
        end if

      enddo

!     output hydrogenic basis data

      if(myid.eq.0 .and. (.not. all(hydro_in_large_basis(i_species, 1:n_hydro(i_species))))) then
         write(use_unit,*)
         write(use_unit,*) " List of hydrogenic basis orbitals: "
         write(use_unit,'(4X,A,4X,A,6X,A,6X,A,2X,A,2X,A,2X,A)') "n", "l", &
              "effective z", "eigenvalue [eV]", "inner max. [A]   ", &
              "outer max. [A]   ","outer radius [A]   "
      end if

      do i_hydro = 1, n_hydro(i_species), 1
        if( hydro_in_large_basis(i_species,i_hydro)) cycle
        i_l     = hydro_l(i_species, i_hydro)
        i_shell = hydro_n(i_species, i_hydro)

        if(myid.eq.0) then
           write(use_unit,'(2X,I3,2X,I3,2X,F15.6,F15.4,4X,F15.6,4X,F15.6,4X,F15.6)') &
                i_shell, i_l, &
                z_eff(i_hydro), eigenval(i_hydro)*hartree, &
                r_inner_max(i_hydro)*bohr,r_outer_max(i_hydro)*bohr, &
                hydro_outer_radius(i_species,i_hydro)*bohr
        end if

      enddo
!    output extra basis function if present
      if(use_ext_basis .and. any(hydro_in_large_basis(i_species, 1:n_hydro(i_species)))) then
        if(myid.eq.0) then
            write(use_unit,*)
            write(use_unit,*) " List of extra hydrogenic orbitals for auxiliary basis: "
            write(use_unit,'(4X,A,4X,A,6X,A,6X,A,2X,A,2X,A,2X,A)') "n", "l", &
              "effective z", "eigenvalue [eV]", "inner max. [A]   ", &
              "outer max. [A]   ","outer radius [A]   "
        end if

        do i_hydro = 1, n_hydro(i_species), 1
            if(.not. hydro_in_large_basis(i_species,i_hydro)) cycle
            i_l     = hydro_l(i_species, i_hydro)
            i_shell = hydro_n(i_species, i_hydro)

            if(myid.eq.0) then
                write(use_unit,'(2X,I3,2X,I3,2X,F15.6,F15.4,4X,F15.6,4X,F15.6,4X,F15.6)') &
                    i_shell, i_l, &
                    z_eff(i_hydro), eigenval(i_hydro)*hartree, &
                    r_inner_max(i_hydro)*bohr,r_outer_max(i_hydro)*bohr, &
                    hydro_outer_radius(i_species,i_hydro)*bohr
            end if

        enddo
      end if

      if(myid.eq.0) then
         write(use_unit,*)
      end if

!  that's all folks

      return
    end subroutine int_hydrogenic_basis_fns
!******

  ! First few associated Laguerre functions for l=0
  real*8 pure function laguerre_l0(n, x) result(y)
    integer, intent(in) :: n
    real*8, intent(in) :: x
    select case(n)
    case(0)
       y = 1
    case(1)
       y = -x+2
    case(2)
       y = 1d0/2*x**2-3*x+3
    case(3)
       y = -1d0/6*x**3+2*x**2-6*x+4
    case(4)
       y = 1d0/24*x**4 - 5d0/6*x**3 + 5*x**2 - 10*x + 5
    case(5)
       y = -1d0/120*x**5 + 1d0/4*x**4 - 5d0/2*x**3 + 10*x**2 - 15*x + 6
    case(6)
       y = 1d0/720*x**6 - 7d0/120*x**5 + 7d0/8*x**4 - 35d0/6*x**3 + &
            & 35d0/2*x**2 - 21d0*x + 7
    case(7)
       y = -1d0/5040*x**7 + 1d0/90*x**6 - 7d0/30*x**5 + 7d0/3*x**4 - &
            & 35d0/3*x**3 + 28*x**2 - 28*x + 8
    case(8)
       y = 1d0/40320*x**8 - 1d0/560*x**7 + 1d0/20*x**6 - 7d0/10*x**5 + &
            & 21d0/4*x**4 - 21d0*x**3 + 42*x**2 - 36*x + 9
    case(9)
       y = -1d0/362880*x**9 + 1d0/4032*x**8 - 1d0/112*x**7 + 1d0/6*x**6 - &
            & 7d0/4*x**5 + 21d0/2*x**4 - 35*x**3 + 60*x**2 - 45*x + 10
    end select
  end function laguerre_l0
