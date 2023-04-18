!****s* FHI-aims/evaluate_v_times_waves
!  NAME
!   evaluate_v_times_waves
!  SYNOPSIS

      subroutine evaluate_v_times_waves &
      ( l_ylm_max, ylm_tab, dist_tab, index_lm, &
        n_compute, i_basis, &
        v_times_radialwaves_spl, v_times_waves &
      )

!  PURPOSE
!     Subroutine evaluate_v_times_waves
!     gives v_times_radial_waves * Y_lm,
!     with the former calculated by the subroutine integrate_v_times_radialwaves
!
!  USES

      use dimensions
      use basis
      use grids
      use geometry
      use spline
      use prodbas
      use mpi_tasks
      use constants, only: pi
      implicit none

!  ARGUMENTS

      integer :: l_ylm_max
      real*8  :: ylm_tab ( (l_ylm_max+1)**2, n_atoms )
      real*8  :: dist_tab(n_atoms)
      integer :: index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max )

      integer :: n_compute
      integer :: i_basis(n_compute)

      real*8  ::  v_times_radialwaves_spl &
                  (n_max_spline, n_hartree_grid, n_loc_prodbas )
      real*8  ::  v_times_waves(n_loc_prodbas)

!  INPUTS
!  o  l_ylm_max -- the maximal l component of the spherical harmonics
!  o  ylm_tab -- tabulated values for the spherical hamiltonics
!  o  dist_tab -- distances to all the atoms
!  o  index_lm -- order (index) for all the different l,m channels
!  o  n_compute -- number of non-zero basis functions
!  o  i_basis -- collection of non-zero  basis functions
!  o  v_times_radialwaves_spl -- the radial part of the integral over bare Coulomb
!           interaction times basis wave functions
!  OUTPUTS
!  o  v_times_waves -- the value of the integral over the bare Coulomb interaction 
!           times basis functions. Now the angular part is also included
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

!  local variables

      real*8   i_r(n_atoms)
      real*8   v_times_radial_waves ( n_loc_prodbas )
!     counters

      integer :: i_loc_prodbas
      integer :: i_compute
      integer :: i_atom, i_species
      integer :: i_basbas
      character(*), parameter :: func = 'evaluate_v_times_waves'

!     begin work

!     tabulate the integral over the coulomb potential and basis  wave function
!        value for each basis function

      do  i_atom = 1, n_atoms

!        if(dist_tab(i_atom).gt.0.d0) then
           i_r(i_atom) = invert_log_grid ( &
                     dist_tab(i_atom), &
                     r_grid_min(species(i_atom)), &
                     r_grid_inc(species(i_atom)) &
                   )

!         endif

      enddo

      do i_compute = 1, n_compute, 1
         i_loc_prodbas = i_basis(i_compute)

         i_basbas = map_prodbas(i_loc_prodbas, myid+1)

         if(i_basbas.gt.0) then
            i_atom = basbas_atom(i_basbas)
            i_species = species(i_atom)

            ! determine the radial part of the integral over the Coulomb
            ! potential and wave function value from basis_atom at current
            ! integration point
            if (dist_tab(i_atom) <= r_grid_max(i_species)) then
               v_times_radial_waves(i_compute) = &
               val_spline &
               ( i_r(i_atom), &
               v_times_radialwaves_spl(1,1,i_loc_prodbas), &
               n_grid(i_species))

            else
               v_times_radial_waves(i_compute) &
               & = 4*pi * multipole_basbas_fn(basbas_fn(i_basbas)) / &
               &   dist_tab(i_atom)**(basbas_l(i_basbas)+1)
            end if



            ! now tabulate the full integral value by multiplying the angular part
            ! v_time_radial_waves_{i_basis}(r) * Y_{atom(i),l(i),m(i)}(theta, phi)

            v_times_waves(i_compute) = &
            ylm_tab(index_lm(basbas_m(i_basbas),basbas_l(i_basbas)), &
            i_atom ) *v_times_radial_waves(i_compute)

            !          if(dist_tab(basbas_atom(i_loc_prodbas)).lt.1.e-16) then
            !            v_times_waves(i_compute) = 0.d0
            !          endif

            !         else

            !          v_times_waves(i_compute) =0.d0
            !         endif

            ! end of if (i_basbas.gt.0)
         endif
         ! end of loop over i_compute
      enddo

      end subroutine evaluate_v_times_waves
!---------------------------------------------------------------------
!******
