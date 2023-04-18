!****s* FHI-aims/evaluate_first_order_U_polarizability
!  NAME
!    evaluate_first_order_U_polarizability
!  SYNOPSIS

subroutine evaluate_first_order_U_polarizability(first_order_H,   &
           KS_eigenvector, KS_eigenvalue, occ_numbers,max_occ_number,   &
           first_order_U,first_order_E)
               

!  PURPOSE
!    calculate the first-order U with its occupied-unoccupied part.

!  Uia(1)=(C(0)S(1)C(0)E(0)-C(0)H(1)C(0))ia/(Eii(1)-Eaa(1))
!  shanghui,2012.04.27

! note: KS_eigenvalue now is only for nspin=1 and ikpoint=1 
! shanghui, 2012.05.30
!  USES

  use dimensions

!  ARGUMENTS

  implicit none
  

  real*8, dimension(n_basis, n_states, n_spin,n_k_points), intent(IN) :: KS_eigenvector
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers
  integer, dimension(n_spin), intent(IN) :: max_occ_number
  real*8, dimension(3, n_basis, n_basis,n_spin), intent(IN) :: first_order_H
  real*8, dimension(3, n_states, n_states,n_spin), intent(OUT) :: first_order_U
  real*8, dimension(3, n_states,n_spin), intent(OUT) ::  first_order_E

!  INPUTS
!   o KS_eigenvector -- Kohn-Sham eigenvectors (real format)
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  first_order_U 


  integer :: i_state, j_state, i_coord, i_spin
  integer :: n_occ_states

  real*8, dimension(n_basis, n_basis) ::temp_first_order
  real*8, dimension(n_states, n_basis) ::temp
  real*8, dimension(n_states, n_states) ::right_2


do i_spin = 1, n_spin 

   n_occ_states=max_occ_number(i_spin)   

  do i_coord=1,3 

  !Ct*H(1)*C 
  temp_first_order(1:n_basis,1:n_basis)=          &
  first_order_H(i_coord,1:n_basis,1:n_basis,i_spin)
  CALL dgemm("T","N",n_states,n_basis,n_basis,&
              1.0d0, KS_eigenvector(:,:,i_spin,1), n_basis,temp_first_order, n_basis,&
              0.0d0,temp,n_states)
  CALL dgemm("N","N",n_states,n_states,n_basis,&
              1.0d0,temp, n_states,KS_eigenvector(:,:,i_spin,1), n_basis,&
              0.0d0,right_2,n_states)


   !Uia
   do i_state=1,n_occ_states
      do j_state=n_occ_states+1, n_states
!----------------------we do not need this term, we calculate here just for debug--------
     first_order_U(i_coord, i_state,j_state,i_spin)=         &
    (-right_2(i_state,j_state)) / &
    (KS_eigenvalue(i_state,i_spin,1)-KS_eigenvalue(j_state,i_spin,1))

!----------------in fact, we only neet U(unocc,occ)=U(a,i) like this:----------
     first_order_U(i_coord, j_state,i_state,i_spin)= &
    (-right_2(j_state,i_state)) / &
    (KS_eigenvalue(j_state,i_spin,1)-KS_eigenvalue(i_state,i_spin,1))
      enddo
   enddo
 
    do i_state=1, n_states
    first_order_E(i_coord,i_state,i_spin)=right_2(i_state,i_state)
    enddo

   enddo ! i_coord

enddo ! i_spin
end subroutine evaluate_first_order_U_polarizability
!******
