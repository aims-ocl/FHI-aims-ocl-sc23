!****s* FHI-aims/evaluate_wave_psi
!  NAME
!   evaluate_wave_psi
!  SYNOPSIS

subroutine evaluate_wave_psi &
     (n_points, h_times_wave, n_compute, &
     KS_ev_compute, max_occ_number, KS_eigenvalue, &
     KS_orbital, i_spin, h_minus_e_times_psi)

!  PURPOSE
!  ??????????????????????????
!
!  USES

  use dimensions
  use constraint
  implicit none

!  ARGUMENTS

  integer :: n_points
  real*8, dimension(n_max_compute_dens, n_points) :: h_times_wave 

  integer :: n_compute
  integer :: max_occ_number
  real*8, dimension(max_occ_number, n_compute) :: KS_ev_compute
  real*8, dimension(n_states) :: KS_eigenvalue
  real*8, dimension(max_occ_number, n_points) :: KS_orbital
  integer :: i_spin

  real*8, dimension(max_occ_number, n_points) ::  h_minus_e_times_psi

! INPUTS
! o n_points -- number of grid points
! o h_times_wave -- Hamiltonian times basis functions
! o n_compute -- number of relevant basis functions
! o KS_ev_compute -- ???????
! o max_occ_number -- number of occupated states
! o KS_eigenvalue -- Kohn-Sham eigenvalues
! o KS_orbital -- Kohn-Sham orbitals
! o i_spin -- spin index
!
! OUTPUT
! o h_minus_e_times_psi -- ??????
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
  ! in case of spin constraint, eigenvalues have to be shifted back
  real*8 :: shift

  !     counter
  integer :: i_state
  integer :: i_region

  !  begin work



  call dgemm('N','N', max_occ_number, n_points, &
     n_compute, 1.0d0, KS_ev_compute, max_occ_number, & 
     h_times_wave, n_max_compute_dens, &
     0.0d0, h_minus_e_times_psi, max_occ_number)




end subroutine evaluate_wave_psi
!---------------------------------------------------------------------
!******	
