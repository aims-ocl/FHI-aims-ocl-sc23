!****s* FHI-aims/evaluate_h_minus_e_times_psi_v2
!  NAME
!   evaluate_h_minus_e_times_psi_v2
!  SYNOPSIS

subroutine evaluate_h_minus_e_times_psi_v2 &
     (n_points, h_times_wave, n_compute, &
     KS_ev_compute, max_occ_number, KS_eigenvalue, &
     KS_orbital, i_spin, h_minus_e_times_psi)

!  PURPOSE
!  Evaluates ( Hamiltonian - eigenvalue) * eigenfunction.
!  This is needed in Pulay forces.  
!
!  USES

  use dimensions
  use constraint
  implicit none

!  ARGUMENTS

  integer :: n_points
  real*8, dimension(n_basis, n_points) :: h_times_wave

  integer :: n_compute
  integer :: max_occ_number
  real*8, dimension(max_occ_number, n_compute) :: KS_ev_compute
  real*8, dimension(n_states) :: KS_eigenvalue
  real*8, dimension(max_occ_number, n_points) :: KS_orbital
  integer :: i_spin

  real*8, dimension(max_occ_number, n_points) :: h_minus_e_times_psi

!  INPUTS
!   o n_points -- number of grid points
!   o h_times_wave -- Hamiltonian * (radial part of basis function)
!   o n_compute -- number of relevant basis functions
!   o max_occ_number -- maximum occupartion number
!   o KS_ev_compute -- Kohn-Sham eigenvalues
!   o KS_eigenvalue -- Kohn-Sham eigenvectors
!   o KS_orbital -- Kohn-Sham orbital
!   o i_spin -- spin index
!
!  OUTPUT
!   o h_minus_e_times_psi -- ( Hamiltonian - eigenvalue) * eigenfunction
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

	


  !     local variables
  ! in case of spin constraint, eigenvalues have to be shifted back
  real*8 :: shift

  !     counter
  integer :: i_state
  integer :: i_region

  !  begin work


  call dgemm('N','N', max_occ_number, n_points, &
       n_compute, 1.0d0, KS_ev_compute, max_occ_number, &
       h_times_wave, n_basis, &
       0.0d0, h_minus_e_times_psi, max_occ_number)


  ! subtract - epsilon_i x psi
  ! unfortunately not a linear iteration
  shift = 0.d0
  do i_state = 1, max_occ_number, 1

     if (use_constraint) then

        if (n_active_regions .gt. 1) then
           shift = 0.d0
           do i_region = 1, n_active_regions, 1
              shift = shift + constraint_proj_final &
                   (i_state, i_region, i_spin) * &
                   constraint_potential(i_region, i_spin)
           end do
        else
           shift = constraint_potential(1, i_spin)
        end if
     end if

     call daxpy(n_points, - KS_eigenvalue(i_state) + shift, &
          KS_orbital(i_state, 1), max_occ_number, &
          h_minus_e_times_psi(i_state, 1), max_occ_number)
  end do



  !  end work


end subroutine evaluate_h_minus_e_times_psi_v2
!---------------------------------------------------------------------
!******	
