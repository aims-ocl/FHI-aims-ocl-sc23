!****s* FHI-aims/evaluate_first_order_gradient_rho_dielectric
!  NAME
!   evaluate_first_order_gradient_rho_dielectric
!  SYNOPSIS

subroutine evaluate_first_order_gradient_rho_dielectric( & 
           n_points,n_compute_c, i_basis_index, & 
           wave, gradient_basis_wave, &
           first_order_density_matrix_compute, & 
           first_order_gradient_rho)

!  PURPOSE


!  USES

  use dimensions
  use runtime_choices
  use basis ! basis_atom
  use pbc_lists
  implicit none

!  ARGUMENTS

  integer :: n_points
  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
 
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8 ,  intent(in) :: first_order_density_matrix_compute(n_compute_c,n_compute_c)
 
  real*8, dimension(3, n_points) :: first_order_gradient_rho


! INPUTS
! o n_points -- number of grid points
! o first_oredr_density_matrix -- relevant numbers of density matrix

!
!  OUTPUT
! o first_order_rho : first_order_rho
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



  !     counters
  integer :: i_coord
  integer :: i_point
  integer :: i_state 
  integer :: i_basis,j_basis,i_compute,j_compute
  real*8, external :: ddot



  !  begin work
  first_order_gradient_rho(1:3,1:n_points)    = 0.0d0


       do i_point=1, n_points,1

         do i_compute=1,n_compute_c
         do j_compute=1,n_compute_c
         
          first_order_gradient_rho(1:3, i_point)= &
          first_order_gradient_rho(1:3, i_point)+ & 
          first_order_density_matrix_compute(i_compute, j_compute) * & 
          ( gradient_basis_wave(i_compute,1:3,i_point)*wave(j_compute,i_point) + & 
            gradient_basis_wave(j_compute,1:3,i_point)*wave(i_compute,i_point) )           
         enddo
         enddo

      end do  ! i_point


end subroutine evaluate_first_order_gradient_rho_dielectric
!---------------------------------------------------------------------
!******	 
