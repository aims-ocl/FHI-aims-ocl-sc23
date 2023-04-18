!****s* FHI-aims/evaluate_first_order_U
!  NAME
!    evaluate_first_order_U_p0
!  SYNOPSIS

subroutine evaluate_first_order_U_p0(first_order_H_complex, first_order_S_complex,  &
           KS_eigenvector_complex, KS_eigenvalue, occ_numbers,   &
           first_order_U_complex,first_order_E)
               

!  PURPOSE
!    calculate the first-order U with its occupied-unoccupied part.

!  Uia(1)=(C(0)S(1)C(0)E(0)-C(0)H(1)C(0))ia/(Eii(1)-Eaa(1))
!  shanghui,2012.04.27

! note: KS_eigenvalue now is only for nspin=1 and ikpoint=1 
! shanghui, 2012.05.30

! shanghui now extended it to pbc case:evaluate_first_order_U_p0.f90
!  USES

  use dimensions

!  ARGUMENTS

  implicit none
  

  complex*16, dimension(n_basis, n_states, n_spin),  intent(IN) :: KS_eigenvector_complex
  real*8, dimension(n_states, n_spin), intent(IN) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin), intent(IN) :: occ_numbers

  complex*16, dimension(3,n_atoms, n_basis, n_basis),intent(IN) :: first_order_S_complex
  complex*16, dimension(3,n_atoms, n_basis, n_basis), intent(IN) :: first_order_H_complex
  complex*16, dimension(3,n_atoms, n_basis, n_basis), intent(OUT) :: first_order_U_complex
  real*8,dimension(3,n_atoms,n_basis), intent(OUT) :: first_order_E

!  INPUTS
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  first_order_U_complex 
!  first_order_E 


  integer :: i_basis, j_basis, i_coord, i_atom, i_spin
  integer :: max_occ_number(n_spin)
  integer :: n_occ_states,i_state


  complex*16 , dimension(n_basis, n_basis) ::temp_first_order
  !real*8, dimension(n_basis, n_basis) ::temp_first_order
  !real*8, dimension(n_basis, n_basis) ::temp
  !real*8, dimension(n_basis, n_basis) ::right_1,right_2

  complex*16 , dimension(n_basis, n_basis) ::temp_1
  complex*16 , dimension(n_basis, n_basis) ::temp_S
  complex*16 , dimension(n_basis, n_basis) ::temp_H
  complex*16 , dimension(n_basis, n_basis) ::temp_eigenvector

  do i_spin = 1, n_spin, 1
      max_occ_number(i_spin) = 0
      do i_state = n_states, 1, -1
       if (dabs(occ_numbers(i_state,i_spin)).gt.0.d0) then
        max_occ_number(i_spin) = i_state
        exit
       endif
      enddo
  enddo

  n_occ_states=max_occ_number(1)  
  ! shanghui need to change this in the future for n_spin


  do i_coord=1,3 
  do i_atom=1,n_atoms
     

    !---> first: temp_S: Ct*S(1)*C*E
    do i_basis=1,n_basis
     temp_first_order(i_basis,1:n_basis)=          &
     first_order_S_complex(i_coord, i_atom, i_basis,1:n_basis)
  
     !write(use_unit,*) 'S(1)'
     !write(use_unit,*) i_basis,temp_first_order(i_basis,1:n_basis)
     !write(use_unit,*) ''

     temp_eigenvector(i_basis,1:n_basis)=          &
     KS_eigenvector_complex(i_basis,1:n_basis,1)

     !write(use_unit,*) 'C'
     !write(use_unit,*) i_basis,temp_eigenvector(i_basis,1:n_basis)
     !write(use_unit,*) ''

     !write(use_unit,*) 'E'
     !write(use_unit,*) i_basis,KS_eigenvalue(i_basis,1)
    enddo

    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
 
    CALL zgemm("C","N",n_basis,n_basis,n_basis,&
               1.0d0, temp_eigenvector, n_basis,temp_first_order, n_basis,&
               0.0d0,temp_1,n_basis)

    temp_S(1:n_basis,1:n_basis)=(0.0d0,0.0d0)

    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                1.0d0,temp_1, n_basis,temp_eigenvector, n_basis,&
                0.0d0,temp_S, n_basis)


    ! Cij=Aij*Bjj  temp_S=CT*S1*C*E0
    do i_basis=1, n_basis
     temp_S(1:n_basis,i_basis)=    &
     temp_S(1:n_basis,i_basis)*KS_eigenvalue(i_basis,1)
    enddo


     !write(use_unit,*) 'E*C*S1*C', temp_S 

    ! Ct*H(1)*C 
    temp_first_order(1:n_basis,1:n_basis)=          &
    first_order_H_complex(i_coord, i_atom,1:n_basis,1:n_basis)

     !write(use_unit,*) 'H(1)',temp_first_order

    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
 
    CALL zgemm("C","N",n_basis,n_basis,n_basis,&
               1.0d0, temp_eigenvector, n_basis,temp_first_order, n_basis,&
               0.0d0,temp_1,n_basis)

    temp_H(1:n_basis,1:n_basis)=(0.0d0,0.0d0)

    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                1.0d0,temp_1, n_basis,temp_eigenvector, n_basis,&
                0.0d0,temp_H, n_basis)
 
     !write(use_unit,*) 'C*H1*C', temp_H

   !Uia
   do i_basis=1,n_occ_states
      do j_basis=n_occ_states+1, n_basis
!----------------------we do not need this term, we calculate here just for debug--------
     first_order_U_complex(i_coord, i_atom, i_basis,j_basis)=         &
    (temp_S(i_basis,j_basis)-temp_H(i_basis,j_basis))/ &
    (KS_eigenvalue(i_basis,1)-KS_eigenvalue(j_basis,1))

!----------------in fact, we only neet U(unocc,occ)=U(a,i) like this:----------
     first_order_U_complex(i_coord, i_atom, j_basis,i_basis)= &
    (temp_S(j_basis,i_basis)-temp_H(j_basis,i_basis))/ &
    (KS_eigenvalue(j_basis,1)-KS_eigenvalue(i_basis,1))
      enddo
   enddo
 
    do i_basis=1, n_basis
    first_order_E(i_coord, i_atom,i_basis)=temp_H(i_basis,i_basis)-temp_S(i_basis,i_basis)
    enddo

   enddo ! i_atom
   enddo ! i_coord

end subroutine evaluate_first_order_U_p0
!******
