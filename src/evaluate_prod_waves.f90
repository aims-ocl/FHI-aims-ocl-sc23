!****s* FHI-aims/evaluate_prod_waves
!  NAME
!   evaluate_prod_waves
!  SYNOPSIS

      subroutine evaluate_prod_waves &
      ( i_r, l_ylm_max, ylm_tab, dist_tab, index_lm, &
        n_basbas_iatom, basbas_iatom, &
        wave_prod &
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

      real*8 i_r(n_atoms)
      integer :: l_ylm_max
      real*8 ylm_tab ( (l_ylm_max+1)**2, n_atoms )
      real*8 dist_tab(n_atoms)
      integer index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max )

      integer :: n_basbas_iatom
      integer :: basbas_iatom(n_basbas_iatom)

!     output

      real*8 wave_prod(n_basbas_iatom)

!  INPUTS
!  o  i_r tabulates the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!  o  l_ylm_max -- maximum of l component for auxiliary basis
!  o  ylm_tab -- tabulated values for the spherical hamiltonics
!  o  dist_tab -- distance to all atoms
!  o  index_lm -- order (index) for all the different l,m channels of the auxiliary
!        basis
!  o  n_basbas_iatom -- number of nonzero auxiliary basis
!  o  basbas_iatom -- list all the nonzero auxiliary basis
!  OUTPUT
!  o  wave_prod -- the values for all the nonzero auxiliary basis in the current 
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

      integer :: i_basis_1
      integer :: i_basbas

      real*8 radial_wave_prod
!     begin work

!     tabulate total wave function value for each basis function

      do i_basbas = 1, n_basbas_iatom, 1

        i_basis_1 = basbas_iatom(i_basbas)

!         determine radial wave function value from basis_atom at current integration point

          radial_wave_prod = &
          val_spline &
          ( i_r(basbas_atom(i_basis_1)), &
            basbas_wave_spl(1,1,basbas_fn(i_basis_1)), &
            n_grid(species(basbas_atom(i_basis_1))) )

!         now tabulate wave function value
!         u_{atom(i),n(i),l(i)}(r) * Y_{atom(i),l(i),m(i)}(theta, phi)

          wave_prod(i_basbas) = &
          ylm_tab(index_lm(basbas_m(i_basis_1), &
                  basbas_l(i_basis_1)), &
          basbas_atom(i_basis_1) ) * &
          radial_wave_prod / &
          dist_tab(basbas_atom(i_basis_1))

!          if(dist_tab(basbas_atom(i_basis_1)).lt.1.e-16) then
!            wave_prod(i_basbas) = 0.d0
!          endif 
!test
!          if (i_basbas.eq.134) then
!            write(use_unit,*) "i_basis(1): ", i_basis_1
!            write(use_unit,*) "wave(1): ", wave_prod(1)
!            write(use_unit,*) "radial_wave(1): ", radial_wave_prod(1)
!            write(use_unit,*) "dist_tab: ", dist_tab(basbas_atom(i_basis_1))
!            write(use_unit,*) "ylm: ",basbas_l(i_basis_1),
!     +            basbas_m(i_basis_1),
!     +      ylm_tab(index_lm(basbas_m(i_basis_1),basbas_l(i_basis_1)),
!     +    basbas_atom(i_basis_1) )
!          end if
!test end

      enddo

!test
!      write(use_unit,*) "   evaluate_wave: i_basis(1) = ", i_basis(1)
!      write(use_unit,*) "   evaluate_wave:    wave(1) = ", wave(1)
!test end

      end subroutine evaluate_prod_waves
!---------------------------------------------------------------------
!******
