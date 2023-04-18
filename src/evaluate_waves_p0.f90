!****s* FHI-aims/evaluate_waves_p0
!  NAME
!   evaluate_waves_p0
!  SYNOPSIS

      subroutine evaluate_waves_p0 &
           ( l_ylm_max, ylm_tab, dist_tab, index_lm,  &
           n_compute, i_basis, radial_wave, wave, &
           n_compute_atoms, atom_index_inv, &
           n_compute_fns, i_basis_fns_inv, n_basis_list  &
           )

!  PURPOSE
!     Prepares only those wave function components which are later needed
!     to evaluate the Hamiltonian or density.
!     There is later updates for this routines: v0 -> v1 -> p0 -> p2
!
!  USES

      use pbc_lists
      use dimensions, only: n_basis_fns, n_centers
      use basis, only: basis_fn, basis_m, basis_l

      implicit none


!  ARGUMENTS

      integer, intent(IN) :: n_compute_atoms
      integer, intent(IN) :: n_compute_fns
      integer, intent(IN) :: n_compute
      integer, intent(IN) :: n_basis_list
      integer, intent(IN) :: l_ylm_max
      real*8, intent(IN) :: ylm_tab ((l_ylm_max+1)**2, n_compute_atoms )
      real*8, intent(IN) :: dist_tab(n_compute_atoms)
      integer, intent(IN) :: index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max )
      integer, intent(IN) :: atom_index_inv(n_centers)
      integer, intent(IN) :: i_basis(n_compute)
      integer, intent(IN) :: i_basis_fns_inv(n_basis_fns,n_centers)
      real*8, intent(IN) :: radial_wave(n_basis_list)
      real*8, intent(OUT) :: wave(n_compute)

!  INPUTS
!    o n_compute_atoms -- number of relevant atoms
!    o n_compute_fns -- number of non-zero basis fns
!    o n_compute -- number of non-zero basis functions 
!    o n_basis_list -- total number of basis functions in whole grid
!    o l_ylm_max -- maximum of l
!    o ylm_tab -- spherical harmonic functions
!    o dist_tab -- distance to atoms
!    o index_lm -- order of l and m components
!    o atom_index_inv -- inverse of atom indexes
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







!  local variables

!     counters

      integer :: i_basis_1
      integer :: i_compute
      integer :: i_point
      integer :: i_compute_fn
      integer :: i_basis_fn_1

!     begin work

!     tabulate total wave function value for each basis function

      do i_compute = 1, n_compute, 1
         i_basis_1 = i_basis(i_compute)
        
!        now tabulate wave function value
!        (1/r)*u_{atom(i),n(i),l(i)}(r) * Y_{atom(i),l(i),m(i)}(theta, phi)

         if (  &
           i_basis_fns_inv(basis_fn(Cbasis_to_basis(i_basis_1)), Cbasis_to_center(i_basis_1)) &
           .gt. 0 ) then

           wave(i_compute) =  &
              ylm_tab(  &
              index_lm(basis_m(Cbasis_to_basis(i_basis_1)),basis_l(Cbasis_to_basis(i_basis_1))),  &
              atom_index_inv(Cbasis_to_center(i_basis_1)) ) *  &
              radial_wave(  &
              i_basis_fns_inv(basis_fn(Cbasis_to_basis(i_basis_1)), &
              Cbasis_to_center(i_basis_1)) ) / &
              dist_tab(  &
              atom_index_inv(Cbasis_to_center(i_basis_1)) )

         else
           wave(i_compute) = 0.d0
         end if

      enddo

    end subroutine evaluate_waves_p0
!---------------------------------------------------------------------
!******
