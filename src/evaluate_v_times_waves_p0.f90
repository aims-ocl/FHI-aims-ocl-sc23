!****s* FHI-aims/evaluate_v_times_waves_p0
!  NAME
!   evaluate_v_times_waves_p0
!  SYNOPSIS

      subroutine evaluate_v_times_waves_p0 &
      (i_task,n_atoms_list,i_r,l_ylm_max, ylm_tab, dist_tab, &
       index_lm, n_loc_compute_current, n_loc_compute, i_loc_prodbas, &
       v_times_radialwaves_spl, v_times_waves &
      )

!  PURPOSE
!     Subroutine evaluate_v_times_waves
!     gives v_times_radial_waves * Y_lm,
!     with the former calculated by the subroutine integrate_v_times_radialwaves
!
!  USES

      use constants, only: pi
      use dimensions
      use basis
      use grids
      use geometry
      use spline
      use prodbas
      use mpi_tasks

      implicit none

!  ARGUMENTS

      integer :: i_task
      integer :: n_atoms_list
      integer :: l_ylm_max
      real*8  :: ylm_tab ( (l_ylm_max+1)**2, n_atoms_list )
      real*8  :: dist_tab(n_atoms_list)
      integer :: index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max )

      integer :: n_loc_compute_current
      integer :: n_loc_compute
      integer :: i_loc_prodbas(n_loc_compute)

      real*8  ::  i_r(n_atoms_list)
      real*8  ::  v_times_radialwaves_spl &
                  (n_max_spline, n_hartree_grid, n_loc_prodbas )
      real*8  ::  v_times_waves(n_loc_compute)

!  INPUTS
!  o  i_task -- current task number, = myid+1
!  o  n_atoms_list -- in practice equals to n_atoms+1, namely the number
!  o       of atoms in a remote unit cell plus the current atom in the
!  o       referene unit cell under consideration
!  o  i_r tabulates the distance from current integration point to all
!  o       atoms in units of the logarithmic grid, i(r)
!  o  l_ylm_max -- the maximal l component of the spherical harmonics
!  o  ylm_tab -- tabulated values for the spherical hamiltonics
!  o  dist_tab -- distances to all the atoms
!  o  index_lm -- order (index) for all the different l,m channels
!  o  n_loc_compute -- number of non-zero basis functions on n_atoms+1
!       atoms (one atom in unit cell at the origin and all atoms in one
!              remote unit cell)
!  o  n_loc_compute_current -- number of non-zero basis functions on 
!       current atom in unit cell at the origin.
!  o  i_loc_prodbas -- collection of non-zero  basis functions
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

      real*8   v_times_radial_waves ( n_loc_compute )
!     counters

      integer :: i_loc_prodbas_1
      integer :: i_compute
      integer :: i_atom
      integer :: i_index

!     begin work

!     tabulate the integral over the coulomb potential and basis  wave function
!        value for each basis function

      do i_compute = 1, n_loc_compute, 1

          i_loc_prodbas_1 = i_loc_prodbas(i_compute)
          i_index = map_prodbas(i_loc_prodbas_1,i_task)
          if(i_index.gt.0 ) then

          if(i_compute .le. n_loc_compute_current) then
            i_atom = 1
          else
            i_atom = basbas_atom(i_index) + 1
          endif

!         if( dist_tab(basbas_atom(i_loc_prodbas_1)).gt.0.d0) then

!         determine the radial part of the integral over the Coulomb potential and
!         wave function value from basis_atom at current integration point

          if (i_r(i_atom) <= n_grid(species(basbas_atom(i_index)))) then
             v_times_radial_waves(i_compute) = &
             val_spline &
             ( i_r(i_atom), &
             v_times_radialwaves_spl(1,1,i_loc_prodbas_1), &
             n_grid(species(basbas_atom(i_index))))
          else
             v_times_radial_waves(i_compute) &
             & = 4*pi * multipole_basbas_fn(basbas_fn(i_index)) / &
             &   dist_tab(i_atom)**(basbas_l(i_index)+1)
          end if

!         now tabulate the full integral value by multiplying the angular part
!         v_time_radial_waves_{i_loc_prodbas}(r) * Y_{atom(i),l(i),m(i)}(theta, phi)

          v_times_waves(i_compute) = &
          ylm_tab(index_lm(basbas_m(i_index),basbas_l(i_index)), &
          i_atom) *v_times_radial_waves(i_compute)

!          if(dist_tab(basbas_atom(i_loc_prodbas_1)).lt.1.e-16) then
!            v_times_waves(i_compute) = 0.d0
!          endif

!         else

!          v_times_waves(i_compute) =0.d0
!         endif

! end of if myid
        endif
! end of loop over i_compute
       enddo

!test
!      write(use_unit,*) "   evaluate_wave: i_loc_prodbas(1) = ", i_loc_prodbas(1)
!      write(use_unit,*) "   evaluate_wave:    wave(1) = ", wave(1)
!test end

      end subroutine evaluate_v_times_waves_p0
!---------------------------------------------------------------------
!******
