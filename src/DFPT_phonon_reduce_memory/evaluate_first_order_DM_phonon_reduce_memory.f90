!****s* FHI-aims/evaluate_first_order_DM_phonon_reduce_memory
!  NAME
!    evaluate_first_order_DM_phonon_reduce_memory
!  SYNOPSIS

subroutine evaluate_first_order_DM_phonon_reduce_memory( &
           first_order_S_complex, & 
           KS_eigenvector, KS_eigenvector_complex, occ_numbers,  &
           first_order_U_complex,first_order_density_matrix_sparse)
               

!  PURPOSE
!    calculate the first-order DM sparse for PBC case 
!    shanghui, 2015.07.30

!  USES

  use dimensions
  use pbc_lists
  use synchronize_mpi
  use scalapack_wrapper
  use runtime_choices

!  ARGUMENTS

  implicit none
  
  real*8, dimension(n_basis, n_states, n_spin,n_k_points_task),  intent(IN) :: KS_eigenvector
  complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task),  intent(IN) :: KS_eigenvector_complex
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers

  complex*16, dimension(n_basis, n_basis,n_k_points_task),intent(IN) :: first_order_S_complex
  complex*16, dimension(n_basis, n_basis,n_k_points_task),intent(IN) :: first_order_U_complex
  complex*16, dimension(n_hamiltonian_matrix_size),intent(INOUT) :: first_order_density_matrix_sparse

!  INPUTS
!   o KS_eigenvector -- Kohn-Sham eigenvectors (real format)
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  first_order_density_matrix_sparse


  integer :: i_basis,j_basis ,i_spin , i_cell, i_index
  integer :: i_state,j_state,i_k_point,  i_k_task
  integer, dimension(n_spin, n_k_points_task)   :: max_occ_number

  complex*16 , dimension(n_basis, n_basis) ::temp_oo
  complex*16 , dimension(n_basis, n_basis) ::temp_1
  complex*16 , dimension(n_basis, n_basis) ::temp_2
  complex*16 , dimension(n_basis, n_basis) ::temp_eigenvector
  complex*16 , dimension(n_basis, n_basis) ::temp_eigenvector_occ
  complex*16 , dimension(n_basis, n_basis) ::temp_C1
  complex*16     ::  zero, one

  zero = (0.d0,0.d0)
  one = (1.d0,0.d0)

 
!-----------------initialize----------------------
  first_order_density_matrix_sparse(1:n_hamiltonian_matrix_size)= (0.0d0,0.0d0)

 if(use_scalapack)then

   do i_spin = 1, n_spin 
  !-----(1) calcualte first_order_dm at every k point------------
  ! DM1(my_k_point) =   C(C^+ S1 C)C^+ occ_number  
  !                   + (CU1) C^+ occ_number 
  !                   +  C   (CU1)^+ occ_number   
    call construct_first_order_dm_scalapack(occ_numbers, i_spin)
   
  !-----(2) add the k point phase ------------------------------
  ! DM1_sparse =  [ DM1(my_k_point) * exp(-ikR) ]
    call get_first_order_dm_complex_sparse_matrix_scalapack(first_order_density_matrix_sparse,i_spin)
  
  !-----(3) sum over {all k points and all matrix} at every core------------------------------   
    call sync_vector_complex(first_order_density_matrix_sparse, n_hamiltonian_matrix_size) 
   enddo

 else   
 !---------------------start lapack version for DM1-----------------------
 ! DM1(i_k_point) =   C(C^+ S1 C)C^+ occ_number  
 !                   + (CU1) C^+ occ_number 
 !                   +  C   (CU1)^+ occ_number   
 !find the max_occ_number
  i_k_task = 0
  do i_k_point = 1, n_k_points,1
     if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
        i_k_task = i_k_task + 1
        do i_spin = 1, n_spin, 1
           max_occ_number(i_spin,i_k_task) = 0
        do i_state = n_states, 1, -1
             if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.0.d0) then
              max_occ_number(i_spin,i_k_task) = i_state
             exit
             endif
        enddo
        enddo
     endif
  enddo


   !-----here DM = Cui^(*) Cvi= Cui_aims * Cvi_aims^(*)-------------- 
   i_k_task = 0
   do i_k_point = 1, n_k_points,1
   if (myid.eq.MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
       i_k_task = i_k_task + 1

  !=================the first term: DM^(1)_oo= C (-C^+ S C) C^+  occ_number=======================
    !          -Ct          S(1) C 
    !        = -Ct_aims^(*) S(1) C_aims^(*) ----->
   do i_basis=1,n_basis
    temp_oo(i_basis,1:n_basis)=          &
    -first_order_S_complex(i_basis,1:n_basis, i_k_task)

    ! Cui = Ciu_aims^(*)
    temp_eigenvector(i_basis,1:n_basis)=          &
    dconjg(KS_eigenvector_complex(i_basis,1:n_basis,1,i_k_task))
    enddo

    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("C","N",n_basis,n_basis,n_basis,&
               one, temp_eigenvector, n_basis,temp_oo, n_basis,&
               zero,temp_1,n_basis)

    temp_2(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                one,temp_1, n_basis,temp_eigenvector, n_basis,&
                zero,temp_2, n_basis)

   !---------[1] the loop version for DM1_oo------------------- 
   ! !temp_oo =  Cui^(*)  ( temp_2)ij  Cvj 
   ! !       =  Cui_aims ( temp_2)ij  Cvj_aims^(*) ----->
   !  temp_oo(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
   !  do i_basis=1,n_basis
   !    do j_basis=1,n_basis
 
   !     do i_state=1, max_occ_number(1,i_k_task)
   !     do j_state=1, max_occ_number(1,i_k_task)
 
   !      if (real_eigenvectors) then
   !         temp_oo(i_basis,j_basis)=temp_oo(i_basis,j_basis)       &
   !       + KS_eigenvector(i_basis,i_state,1,i_k_task) &
   !       * occ_numbers(i_state,1,i_k_point)    & ! if occ_i=occ_j, we are safe.
   !       * temp_2(j_state,i_state)                                     &
   !       * KS_eigenvector(j_basis,j_state,1,i_k_task )
 
   !      else 
   !         temp_oo(i_basis,j_basis)=temp_oo(i_basis,j_basis)       &
   !       + KS_eigenvector_complex(i_basis,i_state,1,i_k_task) &
   !       * occ_numbers(i_state,1,i_k_point)    & ! if occ_i=occ_j, we are safe.
   !       * temp_2(j_state,i_state)                                     &
   !       * dconjg(KS_eigenvector_complex(j_basis,j_state,1,i_k_task ))
   !      endif
 
   !     enddo
   !     enddo
 
   !    enddo
   !  enddo


   !-------------------[2] the Blas version for DM1_oo---------------------------------------------
    ! DM^(1)_oo = C (-C^+ S C) C^T  occ_number
    do i_basis=1,n_basis
       temp_eigenvector_occ(1:n_basis,i_basis)=          &
       KS_eigenvector_complex(1:n_basis,i_basis,1,i_k_task) & 
       *occ_numbers(i_basis,1,i_k_point) 
    enddo

    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("N","T", n_basis, n_basis, max_occ_number(1,i_k_task), &  !transa, transb, m, n, k
               one, temp_eigenvector_occ, n_basis, temp_2, n_basis,   &  ! A, LDA , B, LDB 
               zero,temp_1,n_basis)

    temp_oo(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("N","T", n_basis, n_basis, max_occ_number(1,i_k_task), &
                one,temp_1, n_basis,temp_eigenvector, n_basis, &
                zero,temp_oo, n_basis)
   

 


   !=================the second term  DM1_ov + DM1_vo=============================
   !! DM1_ov + DM1_vo = Uai^(*) Cua^(*)  Cvi    + Cui^(*)  Cva Uai

   !--------[1] the loop version for DM1_ov+DM1_vo------------------- 
   !!----------------the second term-------------------------------
   !! temp_2 = Uai^(*) Cua^(*)  Cvi          + Cui^(*)  Cva          Uai
   !!        = Uai^(*) Cua_aims Cvi_aims^(*) + Cui_aims Cva_aims^(*) Uai 
 
   ! do i_basis=1,n_basis
   !    do j_basis=1,n_basis
   !     temp_2(i_basis,j_basis)=(0.0d0,0.0d0)
 
   !     do i_state=1, max_occ_number(1,i_k_task)
   !     do j_state=max_occ_number(1,i_k_task)+1,n_states
 
   !      if (real_eigenvectors) then
   !       temp_2(i_basis,j_basis)= temp_2(i_basis,j_basis) &
   !       + KS_eigenvector(i_basis,j_state,1,i_k_task)    &
   !       * first_order_U_complex( j_state,i_state,i_k_task)       &
   !       * occ_numbers(i_state,1,i_k_point)                               &
   !       * KS_eigenvector(j_basis,i_state,1,i_k_task )         &
   !       + KS_eigenvector(i_basis,i_state,1,i_k_task)   &
   !       * occ_numbers(i_state,1,i_k_point)                               &
   !       * KS_eigenvector(j_basis,j_state,1,i_k_task )          &
   !       * first_order_U_complex( j_state,i_state,i_k_task)
   !      else 
   !       temp_2(i_basis,j_basis)= temp_2(i_basis,j_basis) &
   !       + KS_eigenvector_complex(i_basis,j_state,1,i_k_task)    &
   !       * dconjg(first_order_U_complex( j_state,i_state,i_k_task))       &
   !       * occ_numbers(i_state,1,i_k_point)                               &
   !       * dconjg(KS_eigenvector_complex(j_basis,i_state,1,i_k_task ))           &
   !       + KS_eigenvector_complex(i_basis,i_state,1,i_k_task)   &
   !       * occ_numbers(i_state,1,i_k_point)                               &
   !       * dconjg(KS_eigenvector_complex(j_basis,j_state,1,i_k_task ))           &
   !       * first_order_U_complex( j_state,i_state,i_k_task)
   !      endif 
 
   !     enddo
   !     enddo
 
   !     enddo
   ! enddo
  


 !---------[2] the BLAS version for DM1_ov+DM1_vo------------------- 
   temp_C1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
               one, temp_eigenvector, n_basis, & 
               first_order_U_complex(1:n_basis,1:n_basis,i_k_task), n_basis,&
               zero,temp_C1,n_basis)
    temp_C1=dconjg(temp_C1) !=  Cua_aims Uai^(*)

    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("N","C",n_basis,n_basis, max_occ_number(1,i_k_task) ,&
               one, temp_C1, n_basis,temp_eigenvector_occ, n_basis,&
               zero,temp_1,n_basis)

    temp_2(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("N","C",n_basis,n_basis, max_occ_number(1,i_k_task),&
                one,temp_eigenvector_occ, n_basis, temp_C1, n_basis,&
                zero,temp_2, n_basis)
    temp_2 = temp_2 + temp_1



 

     !-----------------add k_phase ----------------
     do i_basis = 1, n_basis
     do i_cell = 1,n_cells_in_hamiltonian-1

             if (index_hamiltonian(1,i_cell, i_basis) > 0) then
             do i_index = index_hamiltonian(1, i_cell, i_basis), &
                &         index_hamiltonian(2, i_cell, i_basis)
                   j_basis =  column_index_hamiltonian(i_index)

                  first_order_density_matrix_sparse(i_index)= &
                  first_order_density_matrix_sparse(i_index)  &
                          + dble( (temp_oo(i_basis,j_basis)  &
                                 + temp_2(i_basis,j_basis)) &
                           * k_phase(i_cell,i_k_point) ) &  ! DM(uR,v) = temp(u,v)* exp(-ik(-R))
                           * k_weights(i_k_point)
             end do
             end if

     end do ! i_cell 
     end do ! i_basis

    endif ! i_k_task
    enddo ! i_k_point




  !-------shanghui begin parallel------
   call  sync_vector_complex(first_order_density_matrix_sparse, n_hamiltonian_matrix_size) 
  !-------shanghui end parallel------

 endif 

end subroutine evaluate_first_order_DM_phonon_reduce_memory
!******
