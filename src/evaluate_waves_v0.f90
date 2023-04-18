!****s* FHI-aims/evaluate_waves_v0
!  NAME
!   evaluate_waves_v0
!  SYNOPSIS

      subroutine evaluate_waves_v0 &
           ( i_r, l_ylm_max, ylm_tab, dist_tab, index_lm, &
           n_compute, i_basis, radial_wave, wave &
           )

!  PURPOSE
!     Prepares only those wave function components which are later needed
!     to evaluate the Hamiltonian or density.
!     There is later updates for this routines: v0 -> v1 -> p0 -> p2
!
!  USES

      use dimensions
      use basis
      use grids
      use geometry, only: species
      use spline
      implicit none

!  ARGUMENTS

      real*8, intent(IN) :: i_r(n_atoms)
      integer, intent(IN) :: l_ylm_max
      real*8, intent(IN) :: ylm_tab ( (l_ylm_max+1)**2, n_atoms )
      real*8, intent(IN) :: dist_tab(n_atoms)
      integer, intent(IN) :: index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max )
      integer, intent(IN) :: n_compute
      integer, intent(IN) :: i_basis(n_compute)
      real*8, intent(OUT) :: radial_wave(n_compute)
      real*8, intent(OUT) :: wave(n_compute)

!  INPUTS
!    o i_r -- tabulates the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!    o l_ylm_max -- maximum of l
!    o ylm_tab -- spherical harmonic functions
!    o dist_tab -- distance to atoms
!    o index_lm -- order of l and m components
!    o n_compute -- number of non-zero basis functions      
!    o i_basis -- non-zero  basis functions
!
!  OUTPUT
!    o radial_wave -- radial part of the basis functions.
!    o wave -- total basis functions. wave is defined for only one point here, but may be stored on the outside
!         of the subroutine
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







!  imported variables


! 

!     counters

      integer :: i_basis_1
      integer :: i_compute

!     begin work

!     tabulate total wave function value for each basis function

      do i_compute = 1, n_compute, 1
         i_basis_1 = i_basis(i_compute)

!         determine radial wave function value from basis_atom at current integration point
!        if ( i_basis_fns_inv(basis_fn(Cbasis_to_basis(i_basis_1)),
!     +        Cbasis_to_center(i_basis_1)) .gt. 0 ) then


         radial_wave(i_compute) = &
              val_spline &
              ( i_r(basis_atom(i_basis_1)), &
              basis_wave_spl(1,1,basis_fn(i_basis_1)), &
              n_grid(species(basis_atom(i_basis_1))) )

!         now tabulate wave function value
!         u_{atom(i),n(i),l(i)}(r) * Y_{atom(i),l(i),m(i)}(theta, phi)

         wave(i_compute) = &
             ylm_tab(index_lm(basis_m(i_basis_1),basis_l(i_basis_1)), &
             basis_atom(i_basis_1) ) * radial_wave(i_compute) / &
             dist_tab(basis_atom(i_basis_1))

!test
!          if (i_compute.eq.1) then
!            write(use_unit,*) i_r(:)
!            write(use_unit,*) "radial_wave ", radial_wave(1)
!            write(use_unit,*) "dist_tab: ", dist_tab(basis_atom(i_basis_1))
!            write(use_unit,*) "ylm : ",
!     +      ylm_tab(index_lm(basis_m(i_basis_1),basis_l(i_basis_1)),
!     +      basis_atom(i_basis_1) )
!            write(use_unit,*) "wave ", wave(1)
!          end if
!test end
!       else

!         wave(i_compute) = 0.d0

!       end if

      enddo

!test
!      write(use_unit,*) "   evaluate_wave: i_basis(1) = ", i_basis(1)
!      write(use_unit,*) "   evaluate_wave:    wave(1) = ", wave(1)
!test end

      end subroutine evaluate_waves_v0
!---------------------------------------------------------------------
!
!******
