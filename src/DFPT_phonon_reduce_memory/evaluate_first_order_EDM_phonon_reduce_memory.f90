!****s* FHI-aims/evaluate_first_order_EDM_phonon_reduce_memory
!  NAME
!    evaluate_first_order_EDM_phonon_reduce_memory
!  SYNOPSIS

subroutine evaluate_first_order_EDM_phonon_reduce_memory( &
           first_order_S_complex, &   
           first_order_H_complex,  &
           KS_eigenvector_complex, KS_eigenvalue, occ_numbers,  &
           first_order_U_complex,first_order_energy_density_matrix_sparse)
               

!  PURPOSE
!  (2) calculate zero and first-order EDM for PBC case for phonon
!  EDM()(0)= 
!  EDM()(1)=
!  shanghui,2013.11.28 


!  (3) paralle over k
!  n_k_points : occ_numbers, KS_eigenvalue,k_phase
!  n_k_task   :  KS_eigenvector_complex 

!    calculate the first-order EDM sparse for PBC case 
!    shanghui, 2015.08

!  USES

  use dimensions
  use runtime_choices
  use pbc_lists
  use synchronize_mpi
  use scalapack_wrapper

!  ARGUMENTS

  implicit none
  
  complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task),  intent(IN) :: KS_eigenvector_complex
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers

  complex*16, dimension(n_basis, n_basis,n_k_points_task),intent(IN) :: first_order_S_complex
  complex*16, dimension(n_basis, n_basis,n_k_points_task),intent(IN) :: first_order_H_complex
  complex*16, dimension(n_basis, n_basis,n_k_points_task), intent(IN) :: first_order_U_complex
  complex*16, dimension(n_hamiltonian_matrix_size),intent(INOUT) :: first_order_energy_density_matrix_sparse

!  INPUTS
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o first_order_S_complex -- first order S at every K  
!   o first_order_H_complex -- first order H at every K
!   o first_order_U_complex -- first order U at every K

!
!  OUTPUT
!   o density_matrix
!   o energy_density_matrix
!   o first_order_density_matrix
!   o first_order_energy_density_matrix


  integer :: i_basis,j_basis ,i_spin, i_cell, i_index
  integer :: i_state,j_state,i_k_point,i_k_task
  integer, dimension(n_spin, n_k_points_task)   :: max_occ_number

  complex*16 , dimension(n_basis, n_basis) :: temp_1, temp_H1,temp_S1
  complex*16 , dimension(n_basis, n_basis) :: temp_eigenvector
  complex*16 , dimension(n_basis, n_basis) :: Ck_S1_Ck,Ck_H1_Ck
  complex*16 , dimension(n_basis, n_basis) :: ij_term_Ei,ia_term_Ei
  complex*16     ::  zero, one

  zero = (0.d0,0.d0)
  one = (1.d0,0.d0)


!-----------------initialize----------------------
  first_order_energy_density_matrix_sparse(1:n_hamiltonian_matrix_size)= (0.0d0,0.0d0)

 if(use_scalapack)then

   do i_spin = 1, n_spin 
  !-----(1) calcualte first_order_dm at every k point------------
  ! EDM1(my_k_point) =   C(C^+ H1 C-(Ei+Ej)(C^+ S1 C) )C^+ occ_number  
  !                   + (CU1) C^+ occ_number E  
  !                   +  C  (CU1)^+ occ_number E  
    call construct_first_order_edm_scalapack(occ_numbers, KS_eigenvalue, i_spin)  
 
  !-----(2) add the k point phase ------------------------------
  ! EDM1_sparse =  [ EDM1(my_k_point) * exp(-ikR) ]
    call get_first_order_edm_complex_sparse_matrix_scalapack(first_order_energy_density_matrix_sparse,i_spin)
  
  !-----(3) sum over {all k points and all matrix} at every core------------------------------   
    call sync_vector_complex(first_order_energy_density_matrix_sparse, n_hamiltonian_matrix_size) 
   enddo

 else  
 !---------------------start lapack version for EDM1-----------------------

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

  !temp2=-CtS(1)C----->
 
  i_k_task = 0
  do i_k_point = 1, n_k_points,1
  if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
       i_k_task = i_k_task + 1

   !----------------the first term-------------------------------
    do i_basis=1,n_basis
    temp_S1(i_basis,1:n_basis)=          &
    -first_order_S_complex(i_basis,1:n_basis, i_k_task)

    temp_H1(i_basis,1:n_basis)=          &
    first_order_H_complex(i_basis,1:n_basis, i_k_task)

    ! Cui_in_my_paper = Ciu_aims^(*)
    temp_eigenvector(i_basis,1:n_basis)=          &
    dconjg(KS_eigenvector_complex(i_basis,1:n_basis,1,i_k_task))
    enddo

!----------for Ck_S1_Ck--------------------
!       = -ct          s(1) c 
!       = -ct_aims^(*) s(1) c_aims^(*) ----->

    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("C","N",n_basis,n_basis,n_basis,&
               one, temp_eigenvector, n_basis,temp_S1, n_basis,&
               zero,temp_1,n_basis)

    Ck_S1_Ck(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                one,temp_1, n_basis,temp_eigenvector, n_basis,&
                zero,Ck_S1_Ck, n_basis)

!----------for Ck_H1_Ck--------------------
!       = ct          h(1) c 
!       = ct_aims^(*) h(1) c_aims^(*) ----->
    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("C","N",n_basis,n_basis,n_basis,&
               one, temp_eigenvector, n_basis,temp_H1, n_basis,&
               zero,temp_1,n_basis)

    Ck_H1_Ck(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                one,temp_1, n_basis,temp_eigenvector, n_basis,&
                zero,Ck_H1_Ck, n_basis)


!-----------------------for ij_term_Ei------------------
!      =  Cui^(*)  (Ck_H1_Ck + (ei+ej)Ck_S1_Ck )ji  Cvj 
!      =  Cui_aims (Ck_H1_Ck + (ei+ej)Ck_S1_Ck )ji  Cvj_aims^(*) ----->

     do i_basis=1,n_basis
     do j_basis=1,n_basis

        ij_term_Ei(i_basis,j_basis)=(0.0d0,0.0d0) 

        do i_state=1, max_occ_number(1,i_k_task)
        do j_state=1, max_occ_number(1,i_k_task)
         
            ij_term_Ei(i_basis,j_basis)=ij_term_Ei(i_basis,j_basis)       &
          + KS_eigenvector_complex(i_basis,i_state,1,i_k_task) & 
          * occ_numbers(i_state,1,i_k_point)                            &
          * ( Ck_H1_Ck(j_state,i_state) + & 
              Ck_S1_Ck(j_state,i_state)  & 
              *( KS_eigenvalue(i_state,1,i_k_point)  &    
                +KS_eigenvalue(j_state,1,i_k_point)) ) &
          * dconjg(KS_eigenvector_complex(j_basis,j_state,1,i_k_task))              

        enddo
        enddo
   
     enddo
     enddo

!-----------------------end ij_term----------------------------

  
!-----------------------the ia_term-------------------------------
     do i_basis=1,n_basis
        do j_basis=1,n_basis

         ia_term_Ei(i_basis,j_basis)=(0.0d0,0.0d0)

         do i_state=1, max_occ_number(1,i_k_task)
         do j_state=max_occ_number(1,i_k_task)+1,n_states
              
           ia_term_Ei(i_basis,j_basis)= ia_term_Ei(i_basis,j_basis) &
           + KS_eigenvector_complex(i_basis,j_state,1,i_k_task)    &
           * dconjg( first_order_U_complex(j_state,i_state,i_k_task) )   &
           * occ_numbers(i_state,1,i_k_point)*KS_eigenvalue(i_state,1,i_k_point)  &
           * dconjg( KS_eigenvector_complex(j_basis,i_state,1,i_k_task ))           &
           +  KS_eigenvector_complex(i_basis,i_state,1,i_k_task)   &
           * occ_numbers(i_state,1,i_k_point)*KS_eigenvalue(i_state,1,i_k_point)  &
           * dconjg(KS_eigenvector_complex(j_basis,j_state,1,i_k_task ))           &
           * first_order_U_complex( j_state,i_state,i_k_task)              
         enddo
         enddo 
    
         enddo
     enddo
!-----------------------end the ia_term----------------------------



!-----------------------add k_phase ----------------
     do i_basis = 1, n_basis
     do i_cell = 1,n_cells_in_hamiltonian-1

        if (index_hamiltonian(1,i_cell, i_basis) > 0) then
        do i_index = index_hamiltonian(1, i_cell, i_basis), &
         &         index_hamiltonian(2, i_cell, i_basis)
            j_basis =  column_index_hamiltonian(i_index)

           first_order_energy_density_matrix_sparse(i_index)= &
           first_order_energy_density_matrix_sparse(i_index)  &
                   + dble( (ij_term_Ei(i_basis,j_basis)  &
                          + ia_term_Ei(i_basis,j_basis)) &
                    * k_phase(i_cell,i_k_point) ) &  ! DM(uR,v) = temp(u,v)* exp(-ik(-R))
                    * k_weights(i_k_point)
        end do
        end if

     end do ! i_cell 
     end do ! i_basis

 
    endif
    enddo ! i_k_point


  !-------shanghui begin parallel------
   call  sync_vector_complex(first_order_energy_density_matrix_sparse, n_hamiltonian_matrix_size) 
  !-------shanghui end parallel------

 endif 
end subroutine evaluate_first_order_EDM_phonon_reduce_memory
!******
