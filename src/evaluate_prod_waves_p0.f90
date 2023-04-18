!****s* FHI-aims/evaluate_prod_waves_p0
!  NAME
!   evaluate_prod_waves_p0
!  SYNOPSIS

      subroutine evaluate_prod_waves_p0 &
      ( n_atoms_list, i_r, l_ylm_max, ylm_tab, dist_tab, index_lm, &
        n_compute_current, n_compute, i_prodbas, &
        prod_wave &
      )

!  PURPOSE
!     Subroutine evaluate_waves
!     prepares only those wave function components which are later needed
!     to evaluate the Hamiltonian or density
!
!  USES

      use dimensions
      use basis
      use grids
      use geometry
      use spline
      use prodbas

      implicit none

!  imported variables

!     i_r tabulates the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!     wave is defined for only one point here, but may be stored on the outside
!         of the subroutine

! ARGUMENTS 

      integer :: n_atoms_list
      real*8  :: i_r(n_atoms_list)
      integer :: l_ylm_max
      real*8  :: ylm_tab ( (l_ylm_max+1)**2, n_atoms_list )
      real*8  :: dist_tab(n_atoms_list)
      integer :: index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max )

      integer :: n_compute
      integer :: n_compute_current
      integer :: i_prodbas(n_compute)

!     output

      real*8  :: prod_wave(n_compute)

!  INPUTS
!  o  n_atoms_list -- in practice equals to n_atoms+1, namely the number
!  o       of atoms in a remote unit cell plus the current atom in the
!  o       referene unit cell under consideration
!  o  i_r tabulates the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!  o  l_ylm_max -- maximum of l component for auxiliary basis
!  o  ylm_tab -- tabulated values for the spherical hamiltonics
!  o  dist_tab -- distance to all atoms
!  o  index_lm -- order (index) for all the different l,m channels of the 
!                auxiliary asis
!  o  n_compute -- number of nonzero auxiliary basis for the n_atoms+1
!             many atoms
!  o  n_compute_current -- number of nonzero auxiliary basis functions
!             on the current reference atom (in the original unit cell)
!  o  i_prodbas -- list all the nonzero auxiliary basis
!  OUTPUT
!  o  prod_wave -- the values for all the nonzero auxiliary basis in the current 
!         integration batch
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

!     counters

      integer :: i_basbas
      integer :: i_compute
      integer :: i_atom

      real*8 radial_prod_wave
!     begin work

!     tabulate total wave function value for each basis function

      do i_compute = 1, n_compute, 1

          i_basbas = i_prodbas(i_compute)
          if(i_compute .le. n_compute_current) then
             i_atom = 1 
          else
             i_atom = basbas_atom(i_basbas) + 1
          endif

!         determine radial wave function value from basis_atom at current integration point

          radial_prod_wave = &
          val_spline &
          ( i_r(i_atom), &
            basbas_wave_spl(1,1,basbas_fn(i_basbas)), &
            n_grid(species(basbas_atom(i_basbas))) )

!         now tabulate wave function value
!         u_{atom(i),n(i),l(i)}(r) * Y_{atom(i),l(i),m(i)}(theta, phi)

          prod_wave(i_compute) = &
          ylm_tab(index_lm(basbas_m(i_basbas), &
                  basbas_l(i_basbas)), i_atom) * &
          radial_prod_wave / &
          dist_tab(i_atom)

!          if(dist_tab(basbas_atom(i_basbas)).lt.1.e-16) then
!            prod_wave(i_compute) = 0.d0
!          endif 
!test
!          if (i_compute.eq.134) then
!            write(use_unit,*) "i_basis(1): ", i_basbas
!            write(use_unit,*) "wave(1): ", prod_wave(1)
!            write(use_unit,*) "radial_wave(1): ", radial_prod_wave(1)
!            write(use_unit,*) "dist_tab: ", dist_tab(basbas_atom(i_basbas))
!            write(use_unit,*) "ylm: ",basbas_l(i_basbas),
!     +            basbas_m(i_basbas),
!     +      ylm_tab(index_lm(basbas_m(i_basbas),basbas_l(i_basbas)),
!     +    basbas_atom(i_basbas) )
!          end if
!test end

      enddo

!test
!      write(use_unit,*) "   evaluate_wave: i_basis(1) = ", i_basis(1)
!      write(use_unit,*) "   evaluate_wave:    wave(1) = ", wave(1)
!test end

      end subroutine evaluate_prod_waves_p0
!---------------------------------------------------------------------
!******
