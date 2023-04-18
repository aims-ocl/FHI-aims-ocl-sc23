!****s* FHI-aims/evaluate_nuclear_gradients_p0
!  NAME
!    evaluate_nuclear_gradients_p0
!  SYNOPSIS

subroutine evaluate_nuclear_gradients_p0 &
     (n_points, gradient_basis_wave, n_compute, &
     KS_ev_compute, max_occ_number, i_atom_2, &
     basis_offset, n_local_compute, &
     nuclear_gradient, gga_forces_on, KS_orbital_gradient)

!  PURPOSE
!    Calculates the gradient of the eigenstates with respect to nuclear displacements. 
!    If required, the spatial gradients of the eigenstates are calculated as well which
!    are obtained by the sum of the nuclear gradients.
!  USES  

  use dimensions
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer :: n_points
  real*8, dimension(n_max_compute_dens, 3, n_points) :: gradient_basis_wave

  integer :: n_compute
  integer :: max_occ_number
  real*8, dimension(max_occ_number, n_compute) :: KS_ev_compute
  integer, dimension(n_atoms+1) :: basis_offset
  integer :: n_local_compute
  integer :: i_atom_2

  logical :: gga_forces_on

  real*8, dimension(n_states, n_max_batch_size, 3) :: nuclear_gradient
  real*8, dimension(n_states, n_max_batch_size, 3) ::  KS_orbital_gradient

!  INPUTS
!  o n_points -- number of grid points
!  o gradient_basis_wave -- ????????
!  o n_compute -- number of relevant basis functions
!  o max_occ_number -- maximum occupation number
!  o KS_ev_compute -- ??????
!  o basis_offset -- ????????
!  o n_local_compute -- ????????
!  o i_atom_2 -- ???????????
!
!  OUTPUT
!  o nuclear_gradient -- ???????
!  o KS_orbital_gradient -- gradient of KS orbitals
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
  real*8, dimension(n_local_compute,n_points) :: aux_wave_gradient
  ! number of relevant basis functions sitting on this atom

  !     counter
  integer :: i_coord
  integer :: i_local_compute
  integer :: i_point
  integer :: i_state
  integer :: i_compute

  !  begin work

  !      write(use_unit,*) "in eval_nuclear_gradients...", 
  !     +  max_occ_number, n_states

  do i_coord = 1, 3, 1
     !     condense and reorganize wave_gradient (this is expensive!)
     !     maybe outside the i_atom_2-loop for all i_compute ?
     do i_point = 1, n_points, 1
        do i_local_compute = 1, n_local_compute, 1
           aux_wave_gradient(i_local_compute, i_point) = &
                gradient_basis_wave(basis_offset(i_atom_2) + &
                i_local_compute-1,i_coord,i_point)
        end do
     end do

     call dgemm('N','N', max_occ_number, n_points, &
          n_local_compute, 1.0d0, &
          KS_ev_compute(1,basis_offset(i_atom_2)), &
          max_occ_number, aux_wave_gradient, &
          n_local_compute, 0.0d0,  &
          nuclear_gradient(1,1,i_coord), &
          n_states)

     if (use_density_gradient.and..not.gga_forces_on) then
        ! evaluate whole KS_orbital_gradient if necessary
        do i_point = 1, n_points, 1
           call daxpy(max_occ_number, 1.d0, &
                nuclear_gradient(1,i_point,i_coord), 1, &
                KS_orbital_gradient(1,i_point,i_coord), 1)
        end do
     end if
  end do

  !  end work

end subroutine evaluate_nuclear_gradients_p0
!---------------------------------------------------------------------
!******	
