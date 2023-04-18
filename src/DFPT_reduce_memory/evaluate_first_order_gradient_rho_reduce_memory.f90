!****s* FHI-aims/evaluate_first_order_gradient_rho_reduce_memory
!  NAME
!   evaluate_first_order_gradient_rho_reduce_memory
!  SYNOPSIS

subroutine evaluate_first_order_gradient_rho_reduce_memory( & 
           n_points,n_compute_c, i_basis_index, index_hessian, & 
           wave, gradient_basis_wave, hessian_basis_wave, &
           first_order_density_matrix, density_matrix, & 
           first_order_gradient_rho, & 
           j_atom, j_coord ) 

!  PURPOSE
!  Evaluates first_order_gradient_rho = first_order_DM*gradient(wave*wave) + 
!                                       DM * first_order(gradient(wave*wave))
!  shanghui  2017.03.36

!  USES

  use dimensions
  use runtime_choices
  use basis ! basis_atom
  implicit none

!  ARGUMENTS

  integer :: n_points
  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  integer , intent(in) :: index_hessian(3,3)
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham, 6, n_points),intent(in) :: hessian_basis_wave

  real*8, dimension(n_basis, n_basis), intent(in) :: first_order_density_matrix
  real*8, dimension(n_basis, n_basis), intent(in) :: density_matrix
  real*8, dimension(3, n_points), intent(inout) :: first_order_gradient_rho
  integer, intent(in) :: j_atom, j_coord

! INPUTS
! o n_points -- number of grid points
! o first_oredr_density_matrix -- relevant numbers of density matrix
! o wave -- basis functions
! o gradient_basis_wave -- gradient of basis functions

!
!  OUTPUT
! o first_order_rho : first_order_rho= rho_gradien+ first_order_DM*wave*wave
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
  integer :: i_point
  integer :: i_state 
  integer :: i_basis,j_basis,i_compute,j_compute
  integer :: i_coord
  real*8 :: first_order_gradient_basis_basis(3)


 !     external functions
  real*8, external ::  ddot

  !  begin work
  first_order_gradient_rho(1:3,1:n_points) = 0.0d0
  

  do i_point=1, n_points,1
    

    do i_compute=1,n_compute_c
    do j_compute=1,n_compute_c
    
     i_basis=i_basis_index(i_compute)
     j_basis=i_basis_index(j_compute)          

     first_order_gradient_rho(1:3, i_point)= &
     first_order_gradient_rho(1:3, i_point)+ & 
     first_order_density_matrix(i_basis, j_basis) * & 
     ( gradient_basis_wave(i_compute,1:3,i_point)*wave(j_compute,i_point) + & 
       gradient_basis_wave(j_compute,1:3,i_point)*wave(i_compute,i_point) ) 

      first_order_gradient_basis_basis = 0.0d0 ! initial to 0.0 for safe. 

      if(j_atom.eq.basis_atom(i_basis)) then 
        do i_coord = 1,3 !loop over gradient coord
        first_order_gradient_basis_basis(i_coord) =  &
        first_order_gradient_basis_basis(i_coord)    & 
                -hessian_basis_wave(i_compute,index_hessian(i_coord,j_coord),i_point) & 
                 *wave(j_compute,i_point)  &
                -gradient_basis_wave(i_compute,j_coord,i_point) & 
                 *gradient_basis_wave(j_compute,i_coord,i_point)  
        enddo
      endif 
       
      if(j_atom.eq.basis_atom(j_basis)) then 
        do i_coord = 1,3 !loop over gradient coord
        first_order_gradient_basis_basis(i_coord) =  & 
        first_order_gradient_basis_basis(i_coord)    & 
                -hessian_basis_wave(j_compute,index_hessian(i_coord,j_coord),i_point) & 
                 *wave(i_compute,i_point)  &
                -gradient_basis_wave(i_compute,i_coord,i_point) & 
                 *gradient_basis_wave(j_compute,j_coord,i_point)  
        enddo
      endif 
  
     first_order_gradient_rho(1:3, i_point)= &
     first_order_gradient_rho(1:3, i_point)+ & 
         density_matrix(i_basis, j_basis) * & 
         first_order_gradient_basis_basis(1:3)
       
    enddo
    enddo



 end do  ! i_point


end subroutine evaluate_first_order_gradient_rho_reduce_memory
!---------------------------------------------------------------------
!******	 
