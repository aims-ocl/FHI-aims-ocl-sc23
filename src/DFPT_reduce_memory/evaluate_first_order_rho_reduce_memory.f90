!****s* FHI-aims/evaluate_first_order_rho
!  NAME
!   evaluate_first_order_rho
!  SYNOPSIS

subroutine evaluate_first_order_rho_reduce_memory( & 
           n_points,max_occ_number,  & 
           KS_ev_compute,n_compute_c, i_basis_index, & 
           wave,gradient_basis_wave, &
           density_matrix,first_order_density_matrix, first_order_DM_rho, first_order_rho, & 
           j_atom, j_coord)

!  PURPOSE
!  Evaluates first_order_rho= rho_gradien+ first_order_DM*wave*wave
!  shanghui  2012.05.02

!  Now I add output also first_order_DM_rho
!  shanghui 2012.05.28
!  debug time we can refer to evaluate_density_gradient_denmat.f90
!  USES

  use dimensions
  use runtime_choices
  use basis ! basis_atom
  implicit none

!  ARGUMENTS

  integer :: n_points
  integer :: max_occ_number
  
  real*8, dimension(n_basis,n_states) :: KS_ev_compute

  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave


  real*8, dimension(n_basis,n_basis), intent(in) :: density_matrix
  real*8, dimension(n_basis,n_basis), intent(in) :: first_order_density_matrix
  real*8, dimension(n_points), intent(out) :: first_order_DM_rho
  real*8, dimension(n_points), intent(out) :: first_order_rho
  
  integer , intent(in) :: j_atom, j_coord

! INPUTS
! o n_points -- number of grid points
! o first_oredr_density_matrix -- relevant numbers of density matrix
! o wave -- basis functions
! o gradient_basis_wave -- gradient of basis functions
! o KS_ev_compute -- KS eigenvector

!
!  OUTPUT
! o first_order_DM_rho : first_order_DM*wave*wave
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
  integer :: i_basis,j_basis,i_compute,j_compute, i_coord

  real*8, allocatable :: rho_gradient_atom(:)
  
  real*8, external ::  ddot

  !-------begin allocatation-----------------
  allocate(rho_gradient_atom(n_points))


  !-------begin work-------------------------
  rho_gradient_atom(1:n_points) = 0.0d0
  first_order_DM_rho(1:n_points) = 0.0d0
  first_order_rho(1:n_points) = 0.0d0
  

  !--------------------first: rho_gradient_atom----------------------

   do i_compute = 1,n_compute_c
      i_basis   = i_basis_index(i_compute)
   do j_compute = 1,n_compute_c
      j_basis   = i_basis_index(j_compute) 

      if(j_atom.eq.basis_atom(i_basis)) then
         do i_point=1, n_points,1
          rho_gradient_atom(i_point)= &
          rho_gradient_atom(i_point)+ & 
          density_matrix(i_basis, j_basis) * & 
          (-gradient_basis_wave(i_compute,j_coord,i_point)*wave(j_compute,i_point)) 
         enddo
      endif 

      if(j_atom.eq.basis_atom(j_basis)) then
         do i_point=1, n_points,1
          rho_gradient_atom(i_point)= &
          rho_gradient_atom(i_point)+ & 
          density_matrix(i_basis, j_basis) * & 
          (-gradient_basis_wave(j_compute,j_coord,i_point)*wave(i_compute,i_point)) 
         enddo
      endif 

   enddo
   enddo

  !-------------second: first_order_rho---------------------------------
   do i_compute=1,n_compute_c
      i_basis=i_basis_index(i_compute)
   do j_compute=1,n_compute_c
      j_basis=i_basis_index(j_compute)          

      do i_point=1, n_points,1
          first_order_DM_rho(i_point)= &
          first_order_DM_rho(i_point)+ & 
          first_order_density_matrix(i_basis, j_basis) * & 
          wave(i_compute,i_point)*wave(j_compute,i_point)           
      end do  ! i_point

   enddo
   enddo


   do i_point = 1,n_points 
      first_order_rho(i_point) =                     &
      rho_gradient_atom(i_point) +                    &
      first_order_DM_rho(i_point)
   enddo 



  !-------begin deallocatation-----------------
  deallocate(rho_gradient_atom)

end subroutine evaluate_first_order_rho_reduce_memory
!---------------------------------------------------------------------
!******	 
