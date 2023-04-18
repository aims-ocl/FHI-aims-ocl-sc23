!****s* FHI-aims/evaluate_orb_hess_dot_rho_grad_p0
!  NAME
!    evaluate_orb_hess_dot_rho_grad_p0
!  SYNOPSIS

subroutine evaluate_orb_hess_dot_rho_grad_p0(hessian_basis_wave, n_compute, &
     KS_ev_compute, xc_gradient_deriv, n_points, max_occ_number, i_atom_2, basis_offset, &
     n_local_compute, orb_hess_dot_rho_grad, KS_orbital_gradient, xc_tau_deriv, &
     orb_hess_dot_orb_grad, meta_gga_forces_on)


! PURPOSE
! The subroutine evaluates dot product between KS_orbital_gradients and
! xc_gradient_deriv (including gradient rho respectively)
! needed for gga-forces.
!
! AJL: We also compute orb_hess_dot_orb_grad in here for the Meta-GGA forces
! as all the terms are constructed i.e. this is the most efficient place
!
!  USES

  use dimensions
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer, intent(in) :: n_points
  integer, intent(in) :: n_compute
  integer, intent(in) :: n_local_compute
  integer, intent(in) :: max_occ_number
  integer, intent(in) :: i_atom_2

  real*8,  dimension(n_max_compute_dens, 6, n_points), intent(in) :: hessian_basis_wave
  real*8,  dimension(max_occ_number, n_compute),       intent(in) :: KS_ev_compute
  real*8,  dimension(3, n_points),                     intent(in) :: xc_gradient_deriv
  integer, dimension(n_atoms+1),                       intent(in) :: basis_offset
  real*8, dimension(n_states, n_max_batch_size, 3),   intent(out) :: orb_hess_dot_rho_grad

! Meta GGA requirements
  real*8, dimension(n_states, n_max_batch_size, 3),    intent(in) :: KS_orbital_gradient
  real*8,  dimension(n_points),                        intent(in) :: xc_tau_deriv
  real*8,  dimension(3, n_max_batch_size),            intent(out) :: orb_hess_dot_orb_grad 
  logical,                                             intent(in) :: meta_gga_forces_on

! INPUTS
! o hessian_basis_wave -- ???????????
! o n_compute -- number of relevant basis functions
! o KS_ev_compute -- relevant Kohn-Sham eigenvectors
! o xc_gradient_deriv -- gradient of xc energy
! o n_points -- number of grid points
! o max_occ_number -- number of occupated states
! o i_atom_2 -- atom index
! o basis_offset -- starting point of the spline
! o n_local_compute -- ??????????/
! o KS_orbital_gradient -- relevant KS orbital gradients
! o xc_tau_deriv -- gradient of the xc energy wrt tau
! o meta_gga_forces_on -- flag as to whether we compute the meta-gga forces
!
! OUTPUT
! o orb_hess_dot_rho_grad -- Hessian matrix times gradient of electron density
! o orb_hess_dot_orb_grad -- Hessian matrix times gradient of KS orbitals
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
  real*8, dimension(n_local_compute,n_points) :: aux_wave_hessian
  real*8, dimension(max_occ_number, n_points) :: nuclear_hessian
  ! mGGA local variable
  real*8, dimension(max_occ_number) :: xc_tau_deriv_times_nuclear_hessian

  ! counter
  integer :: i_state
  integer :: i_coord_1
  integer :: i_coord_2
  integer :: i_point
  integer :: i_counter
  integer :: i_local_compute

  ! Function
  real*8, external :: ddot  

  orb_hess_dot_rho_grad = 0.d0
  orb_hess_dot_orb_grad = 0.d0

  i_counter = 0
  do i_coord_1 = 1, 3, 1
     do i_coord_2 = i_coord_1, 3, 1
        i_counter = i_counter + 1 
        ! condense and reorganize wave_gradient (this is expensive!)
        ! maybe outside the i_atom_2-loop for all i_compute ?
        do i_point = 1, n_points, 1
           do i_local_compute = 1, n_local_compute, 1
              aux_wave_hessian(i_local_compute, i_point) = &
                   hessian_basis_wave(basis_offset(i_atom_2) + &
                   i_local_compute-1,i_counter,i_point)
           end do
        end do
        
        call dgemm('N','N', max_occ_number, n_points, n_local_compute, 1.0d0, &
             KS_ev_compute(1, basis_offset(i_atom_2)), max_occ_number, aux_wave_hessian, &
             n_local_compute, 0.0d0, nuclear_hessian, max_occ_number)

        ! calculate dot product from wave_hessian and rho_gradient
        do i_point = 1, n_points, 1
           call daxpy(max_occ_number, xc_gradient_deriv(i_coord_2, i_point), &
                nuclear_hessian(1,i_point), 1, &
                orb_hess_dot_rho_grad(1,i_point,i_coord_1), 1)
        end do
        
        if (i_coord_1 .ne. i_coord_2) then
           do i_point = 1, n_points, 1
              call daxpy(max_occ_number, xc_gradient_deriv(i_coord_1, i_point), &
                   nuclear_hessian(1,i_point), 1, &
                   orb_hess_dot_rho_grad(1,i_point,i_coord_2), 1)
           end do
        end if

        if (meta_gga_forces_on) then

          ! AJL: Here I need to implement changes for the meta-gga
          ! The term we are looking to compute, overall, is:
          !
          ! 2 * [ d(f_xc) / d(tau) ] * grad_at(grad(phi)) . grad(phi)
          !
          ! so by doing this here we can use the already calculated value
          ! of nuclear_hessian = grad(at)grad(phi).
          !
          ! So what do I need to do then? Pass in ther outstanding variables
          ! and add the the terms to some new data object seems to be the necessary.
          ! Can I just then dump the results in gga_forces? Or should we create
          ! meta_gga_forces? TBC. But if it goes in gga_forces remember the different
          ! pre-factors to be applied at the end of update_density_and_forces
          !
          ! For now, we will pass the necessary result into a new object that will
          ! be combined with gga_forces in evaluate_gga_forces. This new object,
          ! orb_hess_dot_orb_grad, contains then:
          !
          ! [ d(f_xc) / d(tau) ] * grad_at(grad(phi)) . grad(phi)

          do i_point = 1, n_points, 1

            xc_tau_deriv_times_nuclear_hessian = 0.d0

            ! Required if we are explicitly doing the vector dot products instead of using dgemm.
            call daxpy(max_occ_number, xc_tau_deriv(i_point), &
                       nuclear_hessian(1,i_point), 1, &
                       xc_tau_deriv_times_nuclear_hessian(1), 1)

            ! Now xc_tau_deriv_times_nuclear_hessian = [ d(f_xc) / d(tau) ] * grad_at(grad(phi))
            ! We just need to perform the dot product with KS_orbital_gradients
            ! Need to make sure the shapes for this are correct.... eurgh.

            orb_hess_dot_orb_grad(i_coord_1,i_point) = &
              orb_hess_dot_orb_grad(i_coord_1,i_point) + &
              ddot(max_occ_number, KS_orbital_gradient(1,i_point,i_coord_2),1, &
                   xc_tau_deriv_times_nuclear_hessian(1),1)

            if (i_coord_1 .ne. i_coord_2) then
               orb_hess_dot_orb_grad(i_coord_2,i_point) = &
                 orb_hess_dot_orb_grad(i_coord_2,i_point) + &
                 ddot(max_occ_number, KS_orbital_gradient(1,i_point,i_coord_1),1, &
                      xc_tau_deriv_times_nuclear_hessian(1),1)

            end if 
          end do
        end if

     end do
  end do

end subroutine evaluate_orb_hess_dot_rho_grad_p0
!******	
