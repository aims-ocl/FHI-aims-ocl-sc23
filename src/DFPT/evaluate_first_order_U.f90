!****s* FHI-aims/evaluate_first_order_U
!  NAME
!    evaluate_first_order_U
!  SYNOPSIS

subroutine evaluate_first_order_U(first_order_H, first_order_S,  &
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
  real*8, dimension(3, n_atoms, n_basis, n_basis), intent(IN) :: first_order_H
  real*8, dimension(3, n_atoms, n_basis, n_basis), intent(IN) :: first_order_S
  real*8, dimension(3, n_atoms, n_basis, n_basis), intent(OUT) :: first_order_U
  real*8, dimension(3, n_atoms, n_basis), intent(OUT) ::  first_order_E

!  INPUTS
!   o KS_eigenvector -- Kohn-Sham eigenvectors (real format)
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  first_order_U 


  integer :: i_basis, j_basis, i_coord, i_atom
  integer :: n_occ_states

  real*8, dimension(n_basis, n_basis) ::temp_first_order
  real*8, dimension(n_basis, n_basis) ::temp
  real*8, dimension(n_basis, n_basis) ::right_1,right_2

  n_occ_states=max_occ_number(1)  ! shanghui need to change this in the future.  


  do i_coord=1,3 
  do i_atom=1,n_atoms
     

  !Ct*S(1)*C*E
  do i_basis=1,n_basis
  temp_first_order(i_basis,1:n_basis)=          &
  first_order_S(i_coord, i_atom, i_basis,1:n_basis)
  enddo

  CALL dgemm("T","N",n_basis,n_basis,n_basis,&
             1.0d0, KS_eigenvector, n_basis,temp_first_order, n_basis,&
             0.0d0,temp,n_basis)
  CALL dgemm("N","N",n_basis,n_basis,n_basis,&
              1.0d0,temp, n_basis,KS_eigenvector, n_basis,&
              0.0d0,right_1,n_basis)

  ! Cij=Aij*Bjj
  do i_basis=1, n_basis
     right_1(1:n_basis,i_basis)=    &
     right_1(1:n_basis,i_basis)*KS_eigenvalue(i_basis,1,1)
  enddo



  !Ct*H(1)*C 
  temp_first_order(1:n_basis,1:n_basis)=          &
  first_order_H(i_coord, i_atom,1:n_basis,1:n_basis)
  CALL dgemm("T","N",n_basis,n_basis,n_basis,&
              1.0d0, KS_eigenvector, n_basis,temp_first_order, n_basis,&
              0.0d0,temp,n_basis)
  CALL dgemm("N","N",n_basis,n_basis,n_basis,&
              1.0d0,temp, n_basis,KS_eigenvector, n_basis,&
              0.0d0,right_2,n_basis)


   !Uia
   do i_basis=1,n_occ_states
      do j_basis=n_occ_states+1, n_basis
!----------------------we do not need this term, we calculate here just for debug--------
     first_order_U(i_coord, i_atom, i_basis,j_basis)=         &
    (right_1(i_basis,j_basis)-right_2(i_basis,j_basis))/ &
    (KS_eigenvalue(i_basis,1,1)-KS_eigenvalue(j_basis,1,1))

!----------------in fact, we only neet U(unocc,occ)=U(a,i) like this:----------
     first_order_U(i_coord, i_atom, j_basis,i_basis)= &
    (right_1(j_basis,i_basis)-right_2(j_basis,i_basis))/ &
    (KS_eigenvalue(j_basis,1,1)-KS_eigenvalue(i_basis,1,1))
      enddo
   enddo
 
    do i_basis=1, n_basis
    first_order_E(i_coord, i_atom,i_basis)=right_2(i_basis,i_basis)-right_1(i_basis,i_basis)
    enddo

   enddo ! i_atom
   enddo ! i_coord

end subroutine evaluate_first_order_U
!******
