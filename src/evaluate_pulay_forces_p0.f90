!****s* FHI-aims/evaluate_pulay_forces_p0
!  NAME
!    evaluate_pulay_forces_p0
!  SYNOPSIS

subroutine evaluate_pulay_forces_p0 &
     (i_atom, h_minus_e_times_psi, nuclear_gradient, &
     partition_tab, max_occ_number, n_points, n_compute_force_atoms, &
     global_atom, pulay_forces)

!  PURPOSE
!    Subroutine evaluate_pulay_forces evaluates all force contributions
!    due to basis function gradients alone, acting on the full Hamiltonian
!  USES 

  use dimensions
  use basis
  use species_data
  use geometry
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer :: max_occ_number
  integer, intent(in) :: n_points
  real*8, dimension(n_points),                               intent(in)    :: partition_tab
  real*8, dimension(n_states, n_max_batch_size, 3), intent(in)             :: nuclear_gradient
  real*8, dimension(max_occ_number, n_points),               intent(in)    :: h_minus_e_times_psi
  integer :: n_compute_force_atoms
  integer, dimension(n_atoms),                               intent(in)    :: global_atom 
  real*8, dimension(3, n_atoms),                             intent(inout) :: pulay_forces
!  INPUTS
!  o max_occ_number -- number of states with non-zero occupation
!  o n_points -- number of grid points if the batch
!  o partition_tab -- values of partition function
!  o nuclear_gradient -- gradients respect atom nucleus
!  o h_minus_e_times_psi -- ( hamiltonian - eigenvalue) * wave
!  o n_compute_force_atoms -- number of relevant atoms in this grid batch
!  o global_atom -- atom origin respect the force is calculated
!  OUTPUT
!  o  pulay_forces -- the Pulay force components are added here. 
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
  integer :: compute_force_atom

  real*8, dimension(n_points) :: grad_times_h 

  ! functions
  real*8, external :: ddot

  ! counters

  integer :: i_coord
  integer :: i_atom
  integer :: i_point
  integer :: i_state

  ! begin work

!test
!  write(use_unit,*) "n_compute_force_atoms: ", n_compute_force_atoms
!  write(use_unit,*) "global_atom : ", global_atom
!  write(use_unit,*) "n_points : ", n_points
!  write(use_unit,*) "max_occ_number : ", max_occ_number
!  write(use_unit,*) "nuclear_gradient(2) ", nuclear_gradient(1:max_occ_number,1:n_points,:,2)
!  write(use_unit,*) "nuclear_gradient(3) ", nuclear_gradient(1:max_occ_number,1:n_points,:,3)
!  write(use_unit,*) "h_minus_e_times_psi ", h_minus_e_times_psi(:,:)
!test end
 
!  do i_atom = 1, n_compute_force_atoms, 1
     compute_force_atom = global_atom(i_atom)
!     write(use_unit,*) constraint_potential(constraint_region(compute_force_atom), i_spin)
     do i_coord = 1, 3, 1
        do i_point = 1, n_points, 1
!           write(use_unit,*) "nuclear gradient"
!           write(use_unit,*) nuclear_gradient(:,i_point,i_coord,i_atom)
!           write(use_unit,*) "h_minus_e_times_psi"
!           write(use_unit,*) h_minus_e_times_psi(:,i_point)
	    grad_times_h(i_point) = &
		  ddot(max_occ_number, nuclear_gradient(1,i_point,i_coord),1,&
		  h_minus_e_times_psi(1,i_point),1)
        end do

       ! NOTE that the Pulay forces still lack a factor of 2, which is added at the very end of update_density
        pulay_forces(i_coord, compute_force_atom) = &
          pulay_forces(i_coord, compute_force_atom) + &
          ddot(n_points,partition_tab,1,grad_times_h,1)

     end do
!  end do

  ! end work

end subroutine evaluate_pulay_forces_p0
!******
!---------------------------------------------------------------------
