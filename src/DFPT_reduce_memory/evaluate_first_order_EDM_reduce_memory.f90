!****s* FHI-aims/evaluate_first_order_EDM
!  NAME
!    evaluate_first_order_EDM
!  SYNOPSIS

subroutine evaluate_first_order_EDM_reduce_memory( & 
           first_order_S, first_order_H, & 
           KS_eigenvector, KS_eigenvalue, &
           occ_numbers, max_occ_number,  &
           first_order_U,first_order_energy_density_matrix)
               

!  PURPOSE
!    calculate the first-order EDM 

!  DM(1)=
!  shanghui,2012.05.02
!  USES

  use dimensions

!  ARGUMENTS

  implicit none
  

  real*8, dimension(n_basis, n_states, n_spin,n_k_points), intent(IN) :: KS_eigenvector
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers
  integer, dimension(n_spin), intent(IN) :: max_occ_number
  real*8, dimension(n_basis, n_basis), intent(IN) :: first_order_S
  real*8, dimension(n_basis, n_basis), intent(IN) :: first_order_H
  real*8, dimension(n_basis, n_basis), intent(IN) :: first_order_U
  real*8, dimension(n_basis,n_basis),intent(INOUT) :: first_order_energy_density_matrix

!  INPUTS
!   o KS_eigenvector -- Kohn-Sham eigenvectors (real format)
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  first_order_density_matrix


  integer :: i_basis,j_basis 
  integer :: i_state,j_state,n_occ_states
  real*8, dimension(n_basis, n_basis) ::temp_first_order
  real*8, dimension(n_basis, n_basis) ::temp_1
  real*8, dimension(n_basis, n_basis) ::temp_2
  real*8, dimension(n_basis, n_basis) ::temp_S
  real*8, dimension(n_basis, n_basis) ::temp_H
  real*8, dimension(n_basis,n_states) :: KS_eigen_times_occnumber

  !--------------initialize----------------------
  first_order_energy_density_matrix(1:n_basis,1:n_basis) = 0.0d0
  n_occ_states=max_occ_number(1) 


!----------------the first term-------------------------------
  !-----(1) CtS(1)C-----------------

  do i_basis=1,n_basis
  temp_first_order(i_basis,1:n_basis)=          &
  first_order_S( i_basis,1:n_basis)
  enddo

  temp_1 = 0.0d0
  CALL dgemm("T","N",n_basis,n_basis,n_basis,&
             1.0d0, KS_eigenvector, n_basis,temp_first_order, n_basis,&
             0.0d0,temp_1,n_basis)
  temp_S = 0.0d0
  CALL dgemm("N","N",n_basis,n_basis,n_basis,&
              1.0d0,temp_1, n_basis,KS_eigenvector, n_basis,&
              0.0d0,temp_S, n_basis)


  !------(2) Ct*H(1)*C ------------------
  temp_first_order(1:n_basis,1:n_basis)=          &
  first_order_H(1:n_basis,1:n_basis)
  temp_1=0.0d0
  CALL dgemm("T","N",n_basis,n_basis,n_basis,&
              1.0d0, KS_eigenvector, n_basis,temp_first_order, n_basis,&
              0.0d0,temp_1,n_basis)
  temp_H=0.0d0
  CALL dgemm("N","N",n_basis,n_basis,n_basis,&
              1.0d0,temp_1, n_basis,KS_eigenvector, n_basis,&
              0.0d0,temp_H,n_basis)


!--------------------DM_oo--------------------------
  do i_state=1, n_states
     KS_eigen_times_occnumber(1:n_basis,i_state)=   &
     KS_eigenvector(1:n_basis,i_state,1,1)*occ_numbers(i_state,1,1)
  enddo

  temp_first_order = 0.0d0
  do i_state=1, n_occ_states
  do j_state=1, n_occ_states
     temp_first_order(i_state,j_state) =  &
   ( temp_H(i_state,j_state) - &
     temp_S(i_state,j_state)*(KS_eigenvalue(i_state,1,1)+KS_eigenvalue(j_state,1,1)) )
  enddo 
  enddo  

  temp_1 = 0.0d0
  CALL dgemm("N","N",n_basis,n_basis,n_basis,&
             1.0d0, KS_eigen_times_occnumber, n_basis,temp_first_order, n_basis,&
             0.0d0,temp_1,n_basis)
  temp_2 = 0.0d0
  CALL dgemm("N","T",n_basis,n_basis,n_basis,&
              1.0d0,temp_1, n_basis,KS_eigenvector, n_basis,&
              0.0d0,temp_2, n_basis)

  first_order_energy_density_matrix = temp_2


!--------------------DM_ov+ DM_vo--------------------------
  do i_state=1, n_states
     KS_eigen_times_occnumber(1:n_basis,i_state)=   &
     KS_eigenvector(1:n_basis,i_state,1,1)*occ_numbers(i_state,1,1)* & 
     KS_eigenvalue(i_state,1,1)
  enddo

  temp_first_order = 0.0d0
  do i_state=1, n_occ_states
  do j_state=n_occ_states+1,n_states
     temp_first_order(j_state,i_state) = first_order_U(j_state,i_state ) 
  enddo 
  enddo  

  temp_1=0.0d0
  temp_2=0.0d0 
  CALL dgemm("N","N",n_basis,n_basis,n_basis,&
             1.0d0, KS_eigenvector, n_basis,temp_first_order, n_basis,&
             0.0d0,temp_1,n_basis)
  CALL dgemm("N","T",n_basis,n_basis,n_basis,&
              1.0d0,temp_1, n_basis,KS_eigen_times_occnumber, n_basis,&
              0.0d0,temp_2, n_basis)
  
 first_order_energy_density_matrix = first_order_energy_density_matrix +  temp_2 

  temp_1=0.0d0
  temp_2=0.0d0 
  CALL dgemm("N","T",n_basis,n_basis,n_basis,&
             1.0d0, KS_eigen_times_occnumber, n_basis,temp_first_order, n_basis,&
             0.0d0,temp_1,n_basis)
  CALL dgemm("N","T",n_basis,n_basis,n_basis,&
              1.0d0,temp_1, n_basis,KS_eigenvector, n_basis,&
              0.0d0,temp_2, n_basis)

 first_order_energy_density_matrix = first_order_energy_density_matrix +  temp_2 


 
end subroutine evaluate_first_order_EDM_reduce_memory
!******
