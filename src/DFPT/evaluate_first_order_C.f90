!****s* FHI-aims/evaluate_first_order_C
!  NAME
!    evaluate_first_order_C
!  SYNOPSIS

subroutine evaluate_first_order_C(first_order_H, first_order_S,  &
           KS_eigenvector, KS_eigenvalue, occ_numbers,max_occ_number,   &
           first_order_density_matrix)
               

!  PURPOSE
!    calculate the first_order_C----> only occ state.
!  shanghui,2012.08.30
! note: KS_eigenvalue now is only for nspin=1 and ikpoint=1 
! shanghui, 2012.08.30
!  USES

  use dimensions
  use physics, only : hamiltonian,overlap_matrix

!  ARGUMENTS

  implicit none
  

  real*8, dimension(n_basis, n_states, n_spin,n_k_points), intent(IN) :: KS_eigenvector
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers
  integer, dimension(n_spin), intent(IN) :: max_occ_number
  real*8, dimension(3, n_atoms, n_basis, n_basis), intent(IN) :: first_order_H
  real*8, dimension(3, n_atoms, n_basis, n_basis), intent(IN) :: first_order_S
  real*8, dimension(3, n_atoms, n_basis, n_basis), intent(OUT) :: first_order_density_matrix

!  INPUTS
!   o KS_eigenvector -- Kohn-Sham eigenvectors (real format)
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  first_order_C


  integer :: i_state,i_basis, j_basis, i_coord, i_atom,j_state
  integer :: n_occ_states,i_index, m,n
  integer :: info,ipiv(n_basis)

  real*8, dimension(n_basis, n_basis) ::temp_first_order
  real*8, dimension(n_basis, n_basis) ::temp
  real*8, dimension(n_basis, n_basis) ::right_1,right_2
  real*8, dimension(n_basis) ::  first_order_E
  real*8, dimension(n_basis,n_basis) :: H_ma, S_ma, P,P_occ,P_unocc
  real*8, dimension(n_basis) :: H1C, SE1C, S1EC
  real*8, dimension(n_basis,n_basis) :: first_order_C

  n_occ_states=max_occ_number(1)  ! shanghui need to change this in the future.  
 

!---------------first:prepare H0,S0,P_occ=PS,P_unocc=I-PS------------------------ 
   do i_basis=1,n_basis
     do j_basis=i_basis,n_basis
        i_index = i_basis + (j_basis-1)*j_basis/2
        H_ma(i_basis,j_basis)= hamiltonian(i_index,1)
        H_ma(j_basis,i_basis)= hamiltonian(i_index,1)
        S_ma(i_basis,j_basis)= overlap_matrix(i_index)
        S_ma(j_basis,i_basis)= overlap_matrix(i_index)
     enddo
   enddo 
  
 
    do i_basis=1, n_basis
      do j_basis=1, n_basis
       P(i_basis,j_basis)=0.0d0
       do i_state=1,n_occ_states
       P(i_basis,j_basis)=P(i_basis,j_basis)+ KS_eigenvector(i_basis,i_state,1,1)* & 
       KS_eigenvector(j_basis,i_state,1,1)
       enddo 
      enddo
    enddo  

    call  matrix_multiplication(n_basis,P,S_ma,P_occ) 

    do i_basis=1,n_basis
      do j_basis=1,n_basis
         if(i_basis.eq.j_basis) then 
         P_unocc(i_basis,j_basis)=1-P_occ(i_basis,j_basis)  
         else
         P_unocc(i_basis,j_basis)=-P_occ(i_basis,j_basis)
         endif
      enddo
    enddo


!-------------------second: calculate first_order_DM for i_coor,i_atom

  do i_coord=1,3 
  do i_atom=1,n_atoms

  
  !++++++++++++++++++++(1) calcuate first_order_E and DM1_oo ++++++++++++++   
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


   do i_basis=1,n_basis
   do j_basis=1,n_basis
     first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)= 0.0d0
      do i_state=1, n_occ_states
      do j_state=1, n_occ_states
 
     first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)=  &
     first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)-  &
         2.0d0*KS_eigenvector(i_basis,i_state,1,1) &
        * right_1(i_state,j_state)              &
        * KS_eigenvector(j_basis,j_state,1,1) 
      enddo
      enddo
   enddo
   enddo


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

 
    do i_basis=1, n_basis
    first_order_E(i_basis)=right_2(i_basis,i_basis)-right_1(i_basis,i_basis)
    !write(use_unit,*) i_basis, first_order_E(i_basis),KS_eigenvector(i_basis,1,1,1)
    enddo



 

   !+++++++++++++++++++++++++++(2) solve Ax=b ++++++++++++++   
   do i_state=1,n_occ_states

    !-------A(n,n): H-S*E-----------------------------
    do i_basis=1, n_basis
    do j_basis=1, n_basis
    temp(i_basis,j_basis)=H_ma(i_basis,j_basis)  &  
               -S_ma(i_basis,j_basis)*KS_eigenvalue(i_state,1,1)  
    enddo
    enddo

    !------B(n,n_occ): -H1*C + S*E1*C + S1*E*C------------ 
    H1C=0.0d0
    SE1C=0.0d0
    S1EC=0.0d0
    do i_basis=1, n_basis   ! our u_index
       do j_basis=1,n_basis ! our v_index
       H1C(i_basis)=H1C(i_basis)+  & 
       first_order_H(i_coord, i_atom, i_basis,j_basis)* &
       KS_eigenvector(j_basis,i_state,1,1)
       
       SE1C(i_basis)=SE1C(i_basis)+  & 
       S_ma(i_basis,j_basis)*first_order_E(i_state)* &
       KS_eigenvector(j_basis,i_state,1,1)
       
       S1EC(i_basis)=S1EC(i_basis)+  & 
       first_order_S(i_coord, i_atom, i_basis,j_basis)* &
       KS_eigenvalue(i_state,1,1)* KS_eigenvector(j_basis,i_state,1,1)
       enddo
    enddo

    right_1(1:n_basis,1)=-H1C(1:n_basis)   & 
    +SE1C(1:n_basis)+S1EC(1:n_basis)

    !--------Ax=B------------------------------------
     call dgesv( n_basis, 1, temp, n_basis,ipiv, right_1, n_basis, info )
     if (info /= 0) write (0, *) 'Error occured in dgesv!'

     
   !+++++++++++++++++++++++++++(3) new_C(1)= P_unocc * C(1)  ++++++++++++++   
     do i_basis=1,n_basis
       first_order_C(i_basis,i_state)=0.0d0
       do j_basis=1,n_basis
      first_order_C(i_basis,i_state)= first_order_C(i_basis,i_state) & 
      +P_unocc(i_basis,j_basis)* right_1(j_basis,1)
       enddo
     enddo


   enddo !i_state


   !++++++++++++++++++++++++++(4) DM1= DM1_oo+ DM1_ov +DM1_vo  ++++++++++++++   
      do i_basis=1,n_basis
        do j_basis=1,n_basis
          do i_state=1,n_occ_states
       first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)=  & 
       first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)+   &
          2*first_order_C(i_basis,i_state)*KS_eigenvector(j_basis,i_state,1,1)+ &
          2*first_order_C(j_basis,i_state)*KS_eigenvector(i_basis,i_state,1,1)
          enddo
        enddo
      enddo

   enddo ! i_atom
   enddo ! i_coord

end subroutine evaluate_first_order_C
!******

subroutine matrix_multiplication(n,A,B,C)
!---------shanghui write for C=A*B---------------
implicit none
!input: n, A,B 
!output: C
integer n,i,j,m
real*8 A(n,n),B(n,n),C(n,n)
 
   do i=1,n
      do j=1,n 
      C(i,j)=0.0d0
        do m=1,n
        C(i,j)=C(i,j)+A(i,m)*B(m,j) 
        enddo
      enddo
   enddo 

end subroutine matrix_multiplication
