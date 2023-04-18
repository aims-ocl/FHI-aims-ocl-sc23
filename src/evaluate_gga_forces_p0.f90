!****s* FHI-aims/evaluate_gga_forces_p0
!  NAME
!    evaluate_gga_forces_p0
!  SYNOPSIS

subroutine evaluate_gga_forces_p0 &
     (i_atom, KS_orbital, nuclear_gradient, orb_grad_dot_rho_grad, &
     orb_hess_dot_rho_grad, &
     partition_tab, max_occ_number, n_points, &
     global_atom, gga_forces, &
     orb_hess_dot_orb_grad, meta_gga_forces_on )

!  PURPOSE
!  Subroutine evaluates the gga forces in set of integration grid points.
!
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
  real*8, dimension(max_occ_number, n_points),               intent(in)    :: KS_orbital
  real*8, dimension(n_points),                               intent(in)    :: partition_tab
  real*8, dimension(n_states, n_max_batch_size, 3),          intent(in)    :: nuclear_gradient
  real*8, dimension(n_states ,n_points),                     intent(in)    :: orb_grad_dot_rho_grad
  real*8, dimension(n_states, n_max_batch_size, 3),          intent(in)    :: orb_hess_dot_rho_grad
  integer, dimension(n_atoms),                               intent(in)    :: global_atom 
  real*8, dimension(3, n_atoms),                             intent(inout) :: gga_forces
! Contributions for meta-GGAs. AJL
  real*8, dimension(3, n_max_batch_size),                    intent(in)    :: orb_hess_dot_orb_grad
  logical,                                                   intent(in)    :: meta_gga_forces_on

!  INPUTS
!   o max_occ_number -- maximum number of states with non-zore occupation
!   o n_points -- number of grid points in the batch
!   o KS_orbital -- Kohn Sman orbitals
!   o partition_tab -- values of partition functions
!   o nuclear_gradient -- gradients respects atom nucleus position
!   o orb_grad_dot_rho_grad -- orbital gradient times density gradient
!   o orb_hess_dot_rho_grad -- orbital hessian matrix times density gradient
!   o global_atom -- atom which forces is currently calculated 
!
!  OUTPUT
!   o gga_forces -- gga force component is added here.
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




  integer :: compute_force_atom

  real*8, dimension(n_points) :: gga_k_term

  ! functions
  real*8, external :: ddot

  ! counters

  integer :: i_coord
  integer :: i_atom
  integer :: i_point

  ! begin work
  compute_force_atom = global_atom(i_atom)
  do i_coord = 1, 3, 1
     do i_point = 1, n_points, 1
        gga_k_term(i_point) = &
           ddot(max_occ_number, nuclear_gradient(1,i_point,i_coord),1,&
           orb_grad_dot_rho_grad(1,i_point),1) + &
           ddot(max_occ_number, KS_orbital(1,i_point),1,&
           orb_hess_dot_rho_grad(1,i_point,i_coord),1)

        if (meta_gga_forces_on) then
          ! Here I need to add in the meta-GGA additions to the gga_k_term,
          ! where it can then be weighted wrt to the partition_tab

          ! The alternative is to keep mGGA forces seperate, and then to
          ! weight them seperately. After this we can look at them separately.
          ! Not really necessary except for aesthetics, so not worrying about this. AJL

          ! We need to include a factor of 0.5 as the meta_gga_forces have a pre-factor of only 2,
          ! where as gga_forces are multiplied by 4 outside of this subroutine

          gga_k_term(i_point) = &
             gga_k_term(i_point) + &
             0.5d0 * orb_hess_dot_orb_grad(i_coord, i_point)
        endif
     end do

     ! NOTE that the gga forces still lack a factor of 4, which is added at the very end of update_density
     gga_forces(i_coord, compute_force_atom) = &
       gga_forces(i_coord, compute_force_atom) + &
       ddot(n_points, partition_tab,1,gga_k_term,1)
  end do

  ! end work

end subroutine evaluate_gga_forces_p0
!---------------------------------------------------------------------
!******
