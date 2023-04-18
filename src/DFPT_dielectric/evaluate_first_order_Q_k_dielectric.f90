!****s* FHI-aims/evaluate_first_order_Q_k_dielectric
!  NAME
!    evaluate_first_order_Q_k_dielectric
!  SYNOPSIS

subroutine evaluate_first_order_Q_k_dielectric( & 
           first_order_S_k_complex, &
           first_order_H_k_complex,   &
           first_order_Q_k_complex)
               

!  PURPOSE
!    calculate the first-order Q with its occupied-unoccupied part.

!  Qij(1)=-0.5*[C(0)S(1)C(0)]ij
!  Qia(1)=[C(0)S(1)C(0)E-C(0)H(1)C(0)]ia/(Eii(1)-Eaa(1))
!  C(1)=C(0)Q(1) ===> C(1)ui =\sum_q C(0)_uq * Q(1)_qi (q = 1,n_basis)


!  USES

  use dimensions
  use runtime_choices, only : real_eigenvectors
  use mpi_tasks
  use physics, only : occ_numbers, KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue  
  use localorb_io, only: use_unit

!  ARGUMENTS

  implicit none
  

  complex*16, dimension( n_basis, n_basis, n_k_points_task), intent(IN) :: first_order_S_k_complex
  complex*16, dimension( n_basis, n_basis, n_k_points_task), intent(IN) :: first_order_H_k_complex
  complex*16, dimension( n_basis, n_basis, n_k_points_task), intent(OUT) :: first_order_Q_k_complex

!  INPUTS
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  first_order_Q_k_complex 


  integer :: i_basis, j_basis, i_spin
  integer :: i_k_task, i_k_point
  integer :: max_occ_number(n_spin)
  integer :: n_occ_states,i_state

  complex*16, allocatable     ::  temp_first_order(:,:)
  complex*16, allocatable     ::  temp_1(:,:)
  complex*16, allocatable     ::  temp_S(:,:)
  complex*16, allocatable     ::  temp_SE(:,:)
  complex*16, allocatable     ::  temp_H(:,:)
  complex*16, allocatable     ::  temp_eigenvector(:,:)
  complex*16     ::  zero, one

  zero = (0.d0,0.d0)
  one = (1.d0,0.d0)


  allocate(temp_first_order(n_basis,n_basis))
  allocate(temp_1(n_basis,n_basis))
  allocate(temp_S(n_basis,n_basis))
  allocate(temp_SE(n_basis,n_basis))
  allocate(temp_H(n_basis,n_basis))
  allocate(temp_eigenvector(n_basis,n_basis))

  i_k_task = 0
  do i_k_point = 1,n_k_points, 1
  if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
        i_k_task = i_k_task + 1

    first_order_Q_k_complex(1:n_basis,1:n_basis, i_k_task) = (0.0d0, 0.0d0)

       do i_spin = 1, n_spin, 1
          max_occ_number(i_spin) = 0
          do i_state = n_states, 1, -1
             !if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.0.d0) then
             ! The following line implies the system does not have any fractional occupation number
             if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.1.e-6) then
             max_occ_number(i_spin) = i_state
             exit
             endif
          enddo
       enddo

       n_occ_states=max_occ_number(1) 
  ! shanghui need to change this in the future for n_spin


   if(real_eigenvectors) then
          temp_eigenvector(:,:)=          &
          cmplx(KS_eigenvector(:,:,1,i_k_task))
   else
          temp_eigenvector(:,:)=          &
          !dconjg(KS_eigenvector_complex(:,:,1,i_k_task))
          KS_eigenvector_complex(:,:,1,i_k_task)
   endif

 !---------(1) Ct*S(1)*C*E ----------------------------------------
    temp_first_order(1:n_basis,1:n_basis)=          &
    first_order_S_k_complex(1:n_basis,1:n_basis, i_k_task)

    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)

    CALL zgemm("C","N",n_basis,n_basis,n_basis,&
               one, temp_eigenvector, n_basis,temp_first_order, n_basis,&
               zero,temp_1,n_basis)

    temp_S(1:n_basis,1:n_basis)=(0.0d0,0.0d0)

    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                one,temp_1, n_basis,temp_eigenvector, n_basis,&
                zero,temp_S, n_basis)

   ! Cij=Aij*Bjj  temp_S=CT*S1*C*E0
    do i_basis=1, n_basis
       temp_SE(1:n_basis,i_basis)=    &
       temp_S(1:n_basis,i_basis)*KS_eigenvalue(i_basis,1,i_k_point)
    enddo
    
 
 !---------(2) Ct*H(1)*C ----------------------------------------
    temp_first_order(1:n_basis,1:n_basis)=          &
    first_order_H_k_complex(1:n_basis,1:n_basis, i_k_task)

    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
 
    CALL zgemm("C","N",n_basis,n_basis,n_basis,&
               one, temp_eigenvector, n_basis,temp_first_order, n_basis,&
               zero,temp_1,n_basis)

    temp_H(1:n_basis,1:n_basis)=(0.0d0,0.0d0)

    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                one,temp_1, n_basis,temp_eigenvector, n_basis,&
                zero,temp_H, n_basis)


    if(i_k_point.eq.1) then
    write(use_unit,*) ' ' ! only when I add this output the complie -O3 at thnec could get right. 
    endif 

!!-----------------------(1) occ-occ terms---------------
!   do i_basis=1,n_occ_states
!      do j_basis=1, n_occ_states
!
!      first_order_Q_k_complex( i_basis,j_basis, i_k_task)=         &
!          -0.5d0*temp_S(i_basis,j_basis) 
!
!      enddo
!   enddo

!-----------------------(2) occ-unocc + unocc-occ terms---------------
   do i_basis=1,n_occ_states
      do j_basis=n_occ_states+1, n_basis
!----------------------we do not need this term, we calculate here just for debug--------
      first_order_Q_k_complex( i_basis,j_basis, i_k_task)=         &
          (temp_SE(i_basis,j_basis)-temp_H(i_basis,j_basis))/ & 
          (KS_eigenvalue(i_basis,1,i_k_point)-KS_eigenvalue(j_basis,1,i_k_point))  

!----------------in fact, we only neet U(unocc,occ)=U(a,i) like this:----------
      first_order_Q_k_complex( j_basis,i_basis,i_k_task)= & 
          (temp_SE(j_basis,i_basis)-temp_H(j_basis,i_basis))/ & 
          (KS_eigenvalue(j_basis,1,i_k_point)-KS_eigenvalue(i_basis,1,i_k_point))

      enddo
   enddo

  endif ! i_k_task   
  enddo ! i_k_point


  deallocate(temp_first_order)
  deallocate(temp_1)
  deallocate(temp_S)
  deallocate(temp_SE)
  deallocate(temp_H)
  deallocate(temp_eigenvector)

end subroutine evaluate_first_order_Q_k_dielectric
!*****
