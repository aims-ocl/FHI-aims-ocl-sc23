!****s* FHI-aims/evaluate_first_zero_order_DM_EDM
!  NAME
!    evaluate_first_zero_order_DM_EDM
!  SYNOPSIS

subroutine evaluate_first_zero_order_DM_EDM_reduce_memory( first_order_S,first_order_H,  &
           KS_eigenvector, KS_eigenvalue, occ_numbers, max_occ_number,  &
           first_order_U,density_matrix,first_order_density_matrix, &  
           energy_density_matrix,first_order_energy_density_matrix)
               

!  PURPOSE
!    calculate the zero and first-order DM 
!    calculate the zero and first-order EDM 
! 
!  shanghui,2013.02.25 @ Berlin 
!  USES

  use dimensions

!  ARGUMENTS

  implicit none
  

  real*8, dimension(n_basis, n_states, n_spin,n_k_points), intent(IN) :: KS_eigenvector
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers
  integer, dimension(n_spin), intent(IN) :: max_occ_number
  real*8, dimension( n_basis, n_basis), intent(IN) :: first_order_S
  real*8, dimension( n_basis, n_basis), intent(IN) :: first_order_U
  real*8, dimension( n_basis, n_basis), intent(IN) :: first_order_H

  real*8, dimension(n_basis,n_basis),intent(OUT) :: density_matrix
  real*8, dimension(n_basis,n_basis),intent(OUT) :: energy_density_matrix
  real*8, dimension( n_basis,n_basis),intent(OUT) :: first_order_density_matrix
  real*8, dimension( n_basis,n_basis),intent(OUT) :: first_order_energy_density_matrix


!  INPUTS
!   o KS_eigenvector -- Kohn-Sham eigenvectors (real format)
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  density_matrix
!  energy_density_matrix
!  first_order_density_matrix
!  first_order_energy_density_matrix

  integer :: i_basis,j_basis 
  integer :: i_state,j_state,n_occ_states
  real*8, dimension(n_basis, n_basis) ::temp_first_order
  real*8, dimension(n_basis, n_basis) ::temp_1, temp_2, temp_3, temp_4
  real*8, dimension(n_basis,n_states) :: KS_eigen_times_occnumber

       
!--------------initialize----------------------       
  n_occ_states=max_occ_number(1) 

  density_matrix=0.0d0
  first_order_density_matrix=0.0d0
  energy_density_matrix=0.0d0
  first_order_energy_density_matrix=0.0d0

  do i_state=1, n_states
     KS_eigen_times_occnumber(1:n_basis,i_state)=   &
     KS_eigenvector(1:n_basis,i_state,1,1)*occ_numbers(i_state,1,1)
  enddo

  
  do i_basis=1,n_basis
    do j_basis=1,n_basis
  
    density_matrix(i_basis,j_basis)=0.0d0
    energy_density_matrix(i_basis,j_basis)=0.0d0
 
      do i_state=1, n_occ_states  
        density_matrix(i_basis,j_basis)= density_matrix(i_basis,j_basis) &
        + KS_eigen_times_occnumber(i_basis,i_state) &
        * KS_eigenvector(j_basis,i_state,1,1)

        energy_density_matrix(i_basis,j_basis)= energy_density_matrix(i_basis,j_basis) &
        + KS_eigenvalue(i_state,1,1)*KS_eigen_times_occnumber(i_basis,i_state) &
        * KS_eigenvector(j_basis,i_state,1,1)
      enddo

    enddo
  enddo





!----------------the occ-occ term-------------------------------
  !temp2=CtS(1)C----->

  do i_basis=1,n_basis
  temp_first_order(i_basis,1:n_basis)=          &
  -first_order_S(i_basis,1:n_basis)
  enddo
 
  temp_1=0.0d0
  CALL dgemm("T","N",n_basis,n_basis,n_basis,&
             1.0d0, KS_eigenvector, n_basis,temp_first_order, n_basis,&
             0.0d0,temp_1,n_basis)
  temp_2=0.0d0
  CALL dgemm("N","N",n_basis,n_basis,n_basis,&
              1.0d0,temp_1, n_basis,KS_eigenvector, n_basis,&
              0.0d0,temp_2, n_basis)


  !for DM(1) =  Cui * ( temp2)ij * Cvj----->
   do i_basis=1,n_basis
     do j_basis=1,n_basis
      temp_1(i_basis,j_basis)=0.0d0
      do i_state=1, n_occ_states
      do j_state=1, n_occ_states
 
        temp_1(i_basis,j_basis)= temp_1(i_basis,j_basis) &
        + KS_eigen_times_occnumber(i_basis,i_state) &
        * temp_2(i_state,j_state)              &
        * KS_eigenvector(j_basis,j_state,1,1) 
 
      enddo
      enddo

 first_order_density_matrix(i_basis,j_basis)=   &
       temp_1(i_basis,j_basis)

     enddo
   enddo

  !Ct*H(1)*C 
  temp_first_order(1:n_basis,1:n_basis)=          &
  first_order_H(1:n_basis,1:n_basis)
  temp_3=0.0d0
  CALL dgemm("T","N",n_basis,n_basis,n_basis,&
              1.0d0, KS_eigenvector, n_basis,temp_first_order, n_basis,&
              0.0d0,temp_3,n_basis)
  temp_4=0.0d0
  CALL dgemm("N","N",n_basis,n_basis,n_basis,&
              1.0d0,temp_3, n_basis,KS_eigenvector, n_basis,&
              0.0d0,temp_4,n_basis)


  ! for EDM(1)= E(0)*Cui * ( temp2)ij * Cvj----->
   do i_basis=1,n_basis
     do j_basis=1,n_basis
      temp_1(i_basis,j_basis)=0.0d0    
      do i_state=1, n_occ_states
      do j_state=1, n_occ_states
 
        temp_1(i_basis,j_basis)= temp_1(i_basis,j_basis) & 
        + KS_eigen_times_occnumber(i_basis,i_state) &
        * ( temp_4(i_state,j_state)     & 
           +temp_2(i_state,j_state)*(KS_eigenvalue(i_state,1,1)+KS_eigenvalue(j_state,1,1)) )    &
        * KS_eigenvector(j_basis,j_state,1,1)  
 
      enddo
      enddo


    first_order_energy_density_matrix(i_basis,j_basis)=   &
         temp_1(i_basis,j_basis)

     enddo
   enddo


 
!-----------------------end the occ-occ term----------------------------

  
!-----------------------the occ-unocc  term-------------------------------
  do i_basis=1,n_basis
     do j_basis=1,n_basis
      temp_first_order(i_basis,j_basis)=0.0d0
      temp_1(i_basis,j_basis)=0.0d0

      do i_state=1, n_occ_states
      do j_state=n_occ_states+1,n_states
        temp_first_order(i_basis,j_basis)= temp_first_order(i_basis,j_basis) &
        +  KS_eigen_times_occnumber(i_basis,i_state) &
        * first_order_U(j_state,i_state )                  &
        * KS_eigenvector(j_basis,j_state,1,1) + KS_eigenvector(i_basis,j_state,1,1) &
        * first_order_U(j_state,i_state )                  &
        *  KS_eigen_times_occnumber(j_basis,i_state) 

        temp_1(i_basis,j_basis)= temp_1(i_basis,j_basis) &
        +  KS_eigen_times_occnumber(i_basis,i_state) &
        * first_order_U(j_state,i_state )                  &
        * KS_eigenvector(j_basis,j_state,1,1)*KS_eigenvalue(i_state,1,1)    & 
        +  KS_eigenvector(i_basis,j_state,1,1) &
        * first_order_U(j_state,i_state )                  &
        *  KS_eigen_times_occnumber(j_basis,i_state)*KS_eigenvalue(i_state,1,1) 
      enddo
      enddo 
              

 first_order_density_matrix(i_basis,j_basis)=   &
 first_order_density_matrix(i_basis,j_basis)+   &
        temp_first_order(i_basis,j_basis)

 first_order_energy_density_matrix(i_basis,j_basis)=   &
 first_order_energy_density_matrix(i_basis,j_basis)+   &
        temp_1(i_basis,j_basis)
      enddo
  enddo
!-----------------------end the occ-unocc term----------------------------



end subroutine evaluate_first_zero_order_DM_EDM_reduce_memory
!******
