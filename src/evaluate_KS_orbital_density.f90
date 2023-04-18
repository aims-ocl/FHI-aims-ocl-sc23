!****s* FHI-aims/evaluate_KS_orbital_density
!  NAME
!    evaluate_KS_orbital_density
!  SYNOPSIS

subroutine evaluate_KS_orbital_density &
     (orbital_number, n_points, wave, n_compute, i_basis, &
     KS_eigenvector,KS_ev_compute, max_num_calc_orb,KS_orbital,rho &
     )

!  PURPOSE
!  Subroutine evaluate_KS_orbital_density takes the local geometry and basis
!  functions at a given integration point and evaluates the charge density
!  for a single KS orbital
!
!  USES

  use dimensions
  implicit none

!  ARGUMENTS

  integer :: n_points
  real*8, dimension(n_basis, n_points) :: wave

  integer :: n_compute
  integer :: i_basis(n_compute)
  integer :: orbital_number
  integer :: max_num_calc_orb

  real*8, dimension(n_basis, n_states) :: KS_eigenvector

  real*8, dimension(max_num_calc_orb, n_points) :: KS_orbital
  real*8, dimension(max_num_calc_orb, n_compute) :: KS_ev_compute
  real*8 :: rho(n_points)


! INPUTS
!  o n_points -- number of grid points
!  o wave -- basis functions
!  o n_compute -- number of relavant basis functions
!  o i_basis -- list of relevant basis functions
!  o orbital_number -- number of the Kohn-Sham orbital to be plotted
!  o max_num_calc_orb -- ?????????
!  o KS_eigenvector -- Kohn-Sham eigenvectors
!
!  OUTPUT
!  o KS_orbital -- Kohn-Sham orbitals
!  o KS_ev_compute -- ??????
!  o rho -- electron density
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


  !     local variables
  real*8, external :: ddot

  !     counters
  integer :: i_state
  integer :: i_compute
  integer :: i_point
  real*8, dimension(n_compute, n_points) :: wave_compute


  do i_compute = 1, n_compute, 1
     KS_ev_compute(1,i_compute) = KS_eigenvector( &
          i_basis(i_compute), orbital_number)
  enddo

  call dgemv('T', n_compute, n_points, 1.0d0, wave, n_basis, &
       KS_ev_compute, max_num_calc_orb, 0.0d0, KS_orbital(1,1) &
       ,max_num_calc_orb)

  do i_point = 1, n_points, 1
     rho(i_point) = KS_orbital(1,i_point)**2.d0
  enddo

  !     end work
end subroutine evaluate_KS_orbital_density
!******	
