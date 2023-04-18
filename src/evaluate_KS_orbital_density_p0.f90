!****s* FHI-aims/evaluate_KS_orbital_density
!  NAME
!    evaluate_KS_orbital_density_p0
!  SYNOPSIS

subroutine evaluate_KS_orbital_density_p0 &
     (i_k_point, n_points, wave, n_compute, n_basis_compute, &
     i_basis, KS_eigenvector_complex, KS_orbital, rho)

!  PURPOSE
!  Evaluates specified Kohn-Sham orbital and corresponding electron
!  density at cube grid points (for periodic calculations)
!
!  USES

  use dimensions

  implicit none

!  ARGUMENTS

  integer :: i_k_point
  integer :: n_points
  integer :: n_compute
  integer :: n_basis_compute
  integer, dimension(n_compute) :: i_basis
  real*8, dimension(n_basis_compute, n_points) :: wave
  complex*16, dimension(n_centers_basis_T) :: KS_eigenvector_complex
  
  real*8, dimension(n_points) :: KS_orbital
  real*8, dimension(n_points) :: rho


!  INPUTS
!   o i_k_point -- k point
!   o n_points -- number of cube grid points
!   o n_compute -- number of relevant basis fns.
!   o n_basis_compute -- maximal number of relevant basis fns.
!   o orbital_number -- label of KS orbital under consideration
!   o i_basis -- list of relevant basis fnd.
!   o wave -- basis functions
!   o KS_eigenvector_complex -- KS eigenvectors
!
!  OUTPUT
!   o KS_orbital -- specified KS orbital at Gamma point
!   o rho -- electron density arising from specified KS orbital
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

! local variables

  integer :: i_compute
  integer :: i_point
  real*8, dimension(n_compute) :: KS_ev_real_part
  
! begin work

! evaluate KS orbital/density at Gamma point only
  if (i_k_point .ne. 1) return
  
  
  do i_compute = 1, n_compute, 1
     KS_ev_real_part(i_compute) = dble(KS_eigenvector_complex(i_basis(i_compute)))
  end do

  
  call dgemv('T', n_compute, n_points, 1.0d0, wave, n_basis_compute, KS_ev_real_part, &
       1, 0.0d0, KS_orbital, 1)
  
  
  do i_point = 1, n_points, 1
     rho(i_point) = KS_orbital(i_point) ** 2.d0
  end do
  
  
end subroutine evaluate_KS_orbital_density_p0
!******
