!****s* FHI-aims/evaluate_first_order_DM
!  NAME
!    evaluate_first_order_DM
!  SYNOPSIS

subroutine evaluate_first_order_DM( first_order_S,  &
           KS_eigenvector, KS_eigenvalue, occ_numbers, max_occ_number,  &
           first_order_U,first_order_density_matrix)
               

!  PURPOSE
!    calculate the first-order DM 

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
  real*8, dimension(3, n_atoms, n_basis, n_basis), intent(IN) :: first_order_S
  real*8, dimension(3, n_atoms, n_basis, n_basis), intent(IN) :: first_order_U
  real*8, dimension(3, n_atoms, n_basis,n_basis),intent(INOUT) :: first_order_density_matrix

!  INPUTS
!   o KS_eigenvector -- Kohn-Sham eigenvectors (real format)
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  first_order_density_matrix


  integer :: i_basis,j_basis ,i_coord, i_atom
  integer :: i_state,j_state,n_occ_states
  real*8, dimension(n_basis, n_basis) ::temp_first_order
  real*8, dimension(n_basis, n_basis) ::temp_1
  real*8, dimension(n_basis, n_basis) ::temp_2
  real*8, dimension(n_basis,n_states) :: KS_eigen_times_occnumber

  n_occ_states=max_occ_number(1) 

  do i_coord=1,3 
  do i_atom=1,n_atoms
     


!----------------the first term-------------------------------
  !temp2=CtS(1)CE----->

  do i_basis=1,n_basis
  temp_first_order(i_basis,1:n_basis)=          &
  -first_order_S(i_coord, i_atom, i_basis,1:n_basis)
  enddo

  CALL dgemm("T","N",n_basis,n_basis,n_basis,&
             1.0d0, KS_eigenvector, n_basis,temp_first_order, n_basis,&
             0.0d0,temp_1,n_basis)
  CALL dgemm("N","N",n_basis,n_basis,n_basis,&
              1.0d0,temp_1, n_basis,KS_eigenvector, n_basis,&
              0.0d0,temp_2, n_basis)


  do i_state=1, n_states
     KS_eigen_times_occnumber(1:n_basis,i_state)=   &
     KS_eigenvector(1:n_basis,i_state,1,1)*occ_numbers(i_state,1,1)
  enddo

  ! Cui * ( temp2)ij * Cvj----->
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

 first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)=   &
       temp_1(i_basis,j_basis)

      !write(use_unit,*) i_basis, j_basis, temp_1(i_basis,j_basis) 
     enddo
   enddo


!--------------------this part is not right for C(Uij+Uji)*CT------------
!--------------------so we block it and use the above ones---------------
!  temp_2(1:4,3)=0.0d0
!  temp_2(1:4,4)=0.0d0
!  KS_eigen_times_occnumber(1:4,3)=0.0d0
!  KS_eigen_times_occnumber(1:4,4)=0.0d0
!  CALL dgemm("N","N",n_basis,n_occ_states,n_occ_states,&
!             1.0d0, KS_eigen_times_occnumber, n_basis,temp_2, n_occ_states,&
!             0.0d0,temp_1,n_basis)
!
!  do i_basis=1,n_basis
!     write(use_unit,*) temp_1(i_basis,1:n_basis)   
!  enddo

!  CALL dgemm("N","T",n_basis,n_basis,n_occ_states,&
!              1.0d0,temp_1, n_basis,temp_first_order, n_basis,&
!              0.0d0,temp_2, n_basis)
!
!  stop
!  first_order_density_matrix(i_coord,i_atom,1:n_basis,1:n_basis)=   &
!       temp_2(1:n_basis,1:n_basis)
!--------------------------end dgemm part-------------------------
 
!-----------------------end the first term----------------------------

  
!-----------------------the second term-------------------------------
  do i_basis=1,n_basis
     do j_basis=1,n_basis
      temp_first_order(i_basis,j_basis)=0.0d0
      do i_state=1, n_occ_states
      do j_state=n_occ_states+1,n_states
        temp_first_order(i_basis,j_basis)= temp_first_order(i_basis,j_basis) &
        +  KS_eigen_times_occnumber(i_basis,i_state) &
        * first_order_U(i_coord, i_atom, j_state,i_state )                  &
        * KS_eigenvector(j_basis,j_state,1,1) + KS_eigenvector(i_basis,j_state,1,1) &
        * first_order_U(i_coord, i_atom, j_state,i_state )                  &
        *  KS_eigen_times_occnumber(j_basis,i_state) 
      enddo
      enddo 
              
      

 first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)=   &
 first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)+   &
        temp_first_order(i_basis,j_basis)
      enddo
  enddo
!-----------------------end the second term----------------------------



   enddo ! i_atom
   enddo ! i_coord

end subroutine evaluate_first_order_DM
!******
