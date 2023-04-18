!****s* FHI-aims/evaluate_waves_v2
!  NAME
!   evaluate_waves_v2
!  SYNOPSIS

      subroutine evaluate_waves_v2 &
           ( l_ylm_max, ylm_tab, dist_tab, index_lm, &
           n_compute, i_basis, radial_wave, wave, &
           n_compute_atoms, atom_index_inv, &
            i_basis_fns_inv &
           )


!  PURPOSE
!     Prepares only those wave function components which are later needed
!     to evaluate the Hamiltonian or density.
!     There is later updates for this routines: v0 -> v1 -> p0 -> p2
!
!  USES


      use dimensions
      use basis
      implicit none


!  ARGUMENTS



      integer, intent(IN) :: n_compute_atoms
      integer, intent(IN) :: l_ylm_max
      real*8, intent(IN) :: ylm_tab ((l_ylm_max+1)**2, n_compute_atoms )
      real*8, intent(IN) :: dist_tab(n_compute_atoms)
      integer, intent(IN) :: index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max )
      integer, intent(IN) :: atom_index_inv(n_atoms)
      integer, intent(IN) :: n_compute
      integer, intent(IN) :: i_basis(n_basis)
      integer, intent(IN) :: i_basis_fns_inv(n_basis_fns,n_atoms)
      real*8, intent(IN) :: radial_wave(n_basis)
      real*8, intent(OUT) :: wave(n_basis)



!  INPUTS
!    o n_compute_atoms -- number of relevant atoms
!    o l_ylm_max -- maximum of l
!    o ylm_tab -- spherical harmonic functions
!    o dist_tab -- distance to atoms
!    o index_lm -- order of l and m components
!    o atom_index_inv -- inverse of atom indexes
!    o n_compute -- number of non-zero basis functions      
!    o i_basis -- non-zero  basis functions
!    o i_basis_fns_inv -- inverse of non-zero  basis functions indexs
!    o radial_wave -- radial part of the basis functions.
!
!  OUTPUT
!   o wave -- total basis functions. wave is defined for only one point here, but may be stored on the outside
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


!

!  local variables

!     counters

      integer :: i_basis_1
      integer :: i_compute




!     begin work

!     tabulate total wave function value for each basis function

      do i_compute = 1, n_compute, 1
         i_basis_1 = i_basis(i_compute)

!        now tabulate wave function value
!        u_{atom(i),n(i),l(i)}(r) * Y_{atom(i),l(i),m(i)}(theta, phi)

         if ( &
           i_basis_fns_inv(basis_fn(i_basis_1), basis_atom(i_basis_1)) &
           .gt. 0 ) then

           wave(i_compute) = &
              ylm_tab( &
              index_lm(basis_m(i_basis_1),basis_l(i_basis_1)), &
              atom_index_inv(basis_atom(i_basis_1)) ) * &
              radial_wave( &
              i_basis_fns_inv(basis_fn(i_basis_1), &
              basis_atom(i_basis_1)) ) / &
              dist_tab( &
              atom_index_inv(basis_atom(i_basis_1)) )

         else
           wave(i_compute) = 0.d0
         end if

      enddo

      end subroutine evaluate_waves_v2
!---------------------------------------------------------------------
!******
