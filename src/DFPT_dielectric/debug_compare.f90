!****s* FHI-aims/debug_compare
!  NAME
!    debug_compare
!  SYNOPSIS

subroutine debug_compare( & 
           momentum_matrix, &
           r_matrix,   &
           KS_eigenvector_complex, KS_eigenvalue,occ_numbers) 
 
               

!  PURPOSE
!    compare if <i|-r|j> = <i|grad|j>/(Ei-Ej)

! shanghui, 2016.02.25

!  USES

  use dimensions
  use localorb_io, only: use_unit

!  ARGUMENTS

  implicit none
  

  complex*16, dimension( n_basis, n_basis), intent(IN) :: momentum_matrix
  complex*16, dimension( n_basis, n_basis), intent(IN) :: r_matrix

  real*8, dimension(n_states, n_spin), intent(IN) :: occ_numbers
  complex*16, dimension(n_basis, n_states, n_spin),  intent(IN) :: KS_eigenvector_complex
  real*8, dimension(n_states, n_spin), intent(IN) :: KS_eigenvalue

!  INPUTS
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  first_order_U_complex 
!  first_order_E 


  integer :: i_basis, j_basis, i_spin
  integer :: max_occ_number(n_spin)
  integer :: n_occ_states,i_state

  complex*16, allocatable     ::  temp_first_order(:,:)
  complex*16, allocatable     ::  temp_1(:,:)
  complex*16, allocatable     ::  temp_H(:,:)
  complex*16, allocatable     ::  temp_M(:,:)
  complex*16, allocatable     ::  temp_eigenvector(:,:)
  complex*16     ::  zero, one

  zero = (0.d0,0.d0)
  one = (1.d0,0.d0)


  allocate(temp_first_order(n_basis,n_basis))
  allocate(temp_1(n_basis,n_basis))
  allocate(temp_H(n_basis,n_basis))
  allocate(temp_M(n_basis,n_basis))
  allocate(temp_eigenvector(n_basis,n_basis))

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



   do i_basis=1,n_basis

        temp_eigenvector(i_basis,1:n_basis)=          &
        KS_eigenvector_complex(i_basis,1:n_basis,1)

   enddo

     

 !---------(1) Ct*H(1)*C ----------------------------------------
    temp_first_order(1:n_basis,1:n_basis)=          &
    r_matrix(1:n_basis,1:n_basis)

    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
 
    CALL zgemm("C","N",n_basis,n_basis,n_basis,&
               one, temp_eigenvector, n_basis,temp_first_order, n_basis,&
               zero,temp_1,n_basis)

    temp_H(1:n_basis,1:n_basis)=(0.0d0,0.0d0)

    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                one,temp_1, n_basis,temp_eigenvector, n_basis,&
                zero,temp_H, n_basis)


 !---------(2) Ct*momentum*C ----------------------------------------
    temp_first_order(1:n_basis,1:n_basis)=          &
    momentum_matrix(1:n_basis,1:n_basis)

    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
 
    CALL zgemm("C","N",n_basis,n_basis,n_basis,&
               one, temp_eigenvector, n_basis,temp_first_order, n_basis,&
               zero,temp_1,n_basis)

    temp_M(1:n_basis,1:n_basis)=(0.0d0,0.0d0)

    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                one,temp_1, n_basis,temp_eigenvector, n_basis,&
                zero,temp_M, n_basis)

    !if(i_k_point.eq.1) then
    !write(use_unit,*) ' ' ! only when I add this output the complie -O3 at thnec could get right. 
    !endif 

   !Uia : we add H^(1) (electric filed part) here : 
   do i_basis=1,n_occ_states
      do j_basis=n_occ_states+1, n_basis
!----------------------we do not need this term, we calculate here just for debug--------
      write(use_unit,*) 'r matrix'
      write(use_unit,*) i_basis, j_basis, -temp_H(i_basis,j_basis)

      write(use_unit,*) 'moment matrix'     
      write(use_unit,*) i_basis, j_basis, & 
         temp_M(i_basis,j_basis)/(KS_eigenvalue(i_basis,1)-KS_eigenvalue(j_basis,1))

      enddo
   enddo

  deallocate(temp_first_order)
  deallocate(temp_1)
  deallocate(temp_H)
  deallocate(temp_M)
  deallocate(temp_eigenvector)

end subroutine debug_compare
!*****
