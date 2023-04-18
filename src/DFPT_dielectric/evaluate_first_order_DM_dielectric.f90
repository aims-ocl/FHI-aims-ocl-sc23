!****s* FHI-aims/evaluate_first_order_DM_dielectric
!  NAME
!    evaluate_first_order_DM_dielectric
!  SYNOPSIS

subroutine evaluate_first_order_DM_dielectric( &
           KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers,  &
           first_order_U_complex,first_order_density_matrix_sparse)
               

!  PURPOSE
!    calculate the first-order DM sparse for PBC case 
!    shanghui, 2015.07.30

!  USES

  use dimensions
  use pbc_lists
  use runtime_choices
  use synchronize_mpi
  use scalapack_wrapper

!  ARGUMENTS

  implicit none
  
  real*8, dimension(n_basis, n_states, n_spin,n_k_points_task),  intent(IN) :: KS_eigenvector
  complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task),  intent(IN) :: KS_eigenvector_complex
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers

  complex*16, dimension(n_states,n_states,n_k_points_task),intent(IN) :: first_order_U_complex
  real*8, dimension(n_hamiltonian_matrix_size),intent(INOUT) :: first_order_density_matrix_sparse

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

  real*8, allocatable         ::  temp_first_order(:,:)
  real*8, allocatable         ::  temp_1(:,:)
  real*8, allocatable         ::  temp_2(:,:)
  real*8, allocatable         ::  KS_eigen_times_occnumber(:,:)
  real*8, allocatable         ::  temp_real(:,:)
  
  complex*16 , allocatable :: temp_1_complex(:,:)
  complex*16 , allocatable :: temp_eigenvector(:,:)
  complex*16 , allocatable :: temp_eigenvector_occ(:,:)
  complex*16 , allocatable :: temp_C1(:,:)
  complex*16 , allocatable :: temp_complex(:,:) 
  complex*16     ::  zero, one

  zero = (0.d0,0.d0)
  one = (1.d0,0.d0)
 
!-----------------initialize----------------------
  first_order_density_matrix_sparse(1:n_hamiltonian_matrix_size)= 0.0d0

 if(use_scalapack)then

   do i_spin = 1, n_spin
  !-----(1) calcualte first_order_dm at every k point------------
  ! DM1(my_k_point) =    
  !                   + (CU1) C^+ occ_number 
  !                   +  C   (CU1)^+ occ_number   
    call construct_first_order_dm_dielectric_scalapack(occ_numbers, i_spin)

  !-----(2) add the k point phase ------------------------------
  ! DM1_sparse =  [ DM1(my_k_point) * exp(-ikR) ]
    call get_first_order_dm_complex_sparse_matrix_dielectric_scalapack(first_order_density_matrix_sparse,i_spin)

  !-----(3) sum over {all k points and all matrix} at every core------------------------------   
    call sync_vector(first_order_density_matrix_sparse, n_hamiltonian_matrix_size)
   enddo

 else

!---------------------start lapack version for DM1-----------------------
  if(real_eigenvectors) then 
     allocate(temp_first_order(n_basis,n_basis)) 
     allocate(temp_1(n_basis,n_states)) 
     allocate(temp_2(n_basis,n_basis)) 
     allocate(KS_eigen_times_occnumber(n_basis,n_states)) 
     allocate(temp_real(n_basis,n_basis)) 
  else 
     allocate(temp_1_complex(n_basis,n_basis))
     allocate(temp_eigenvector(n_basis,n_basis))
     allocate(temp_eigenvector_occ(n_basis,n_basis))
     allocate(temp_C1(n_basis,n_basis))
     allocate(temp_complex(n_basis,n_basis))
  endif

  !find the max_occ_number
  i_k_task = 0
  do i_k_point = 1, n_k_points,1

     if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
        i_k_task = i_k_task + 1
        do i_spin = 1, n_spin, 1
           max_occ_number(i_spin,i_k_task) = 0

          do i_state = n_states, 1, -1
             ! The following line implies the system does not have any fractional occupation number
             if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.1.e-6) then
              max_occ_number(i_spin,i_k_task) = i_state
             exit
             endif
          enddo

        enddo ! i_spin
     endif
  enddo

   i_k_task = 0
   do i_k_point = 1, n_k_points,1
   if (myid.eq.MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
       i_k_task = i_k_task + 1

    !===================begin real_eigenvectors=========================
    if (real_eigenvectors) then
        temp_real = 0.0d0

        do i_state=1, n_states
           KS_eigen_times_occnumber(1:n_basis,i_state)=   &
           KS_eigenvector(1:n_basis,i_state,1,i_k_task)*occ_numbers(i_state,1,i_k_point)
        enddo

        temp_first_order = 0.0d0
        do i_state=1, max_occ_number(1,i_k_task)
        do j_state=max_occ_number(1,i_k_task)+1,n_states
           temp_first_order(j_state,i_state) = dble(first_order_U_complex(j_state,i_state,i_k_task)) 
        enddo 
        enddo  
       
        temp_1=0.0d0
        temp_2=0.0d0 
        CALL dgemm("N","N",n_basis,n_states,n_states,&
                   1.0d0, KS_eigenvector, n_basis,temp_first_order, n_states,&
                   0.0d0,temp_1,n_basis)
        CALL dgemm("N","T",n_basis,n_basis,n_states,&
                    1.0d0,temp_1, n_basis,KS_eigen_times_occnumber, n_basis,&
                    0.0d0,temp_2, n_basis)
        
        temp_real = temp_real + temp_2 
       
        temp_1=0.0d0
        temp_2=0.0d0 
        CALL dgemm("N","T",n_basis,n_states,n_states,&
                   1.0d0, KS_eigen_times_occnumber, n_basis,temp_first_order, n_states,&
                   0.0d0,temp_1,n_basis)
        CALL dgemm("N","T",n_basis,n_basis,n_states,&
                    1.0d0,temp_1, n_basis,KS_eigenvector, n_basis,&
                    0.0d0,temp_2, n_basis)
       
        temp_real = temp_real + temp_2 
    !===================end real_eigenvectors=========================

    else 

    !===================begin complex_eigenvectors=========================

    ! DM   =   Ciu_aims * Cvi_aims^(*) 
    ! DM1  =   Cua_aims Uai Cvi_aims^(*) + Cui_aims Cva_aims^(*) Uai^(*) 
       do i_basis=1,n_basis
           temp_eigenvector(1:n_basis,i_basis)=          &
           KS_eigenvector_complex(1:n_basis,i_basis,1,i_k_task)

           temp_eigenvector_occ(1:n_basis,i_basis)=          &
           KS_eigenvector_complex(1:n_basis,i_basis,1,i_k_task) &
           *occ_numbers(i_basis,1,i_k_point)
       enddo

       temp_C1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
       CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                  one, temp_eigenvector, n_basis, & 
                  first_order_U_complex(1:n_basis,1:n_basis,i_k_task), n_basis,&
                  zero,temp_C1,n_basis)

       temp_1_complex(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
       CALL zgemm("N","C",n_basis,n_basis,max_occ_number(1,i_k_task),&
                  one, temp_C1, n_basis,temp_eigenvector_occ, n_basis,&
                  zero,temp_1_complex,n_basis)

       temp_complex(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
       CALL zgemm("N","C",n_basis,n_basis,max_occ_number(1,i_k_task),&
                   one,temp_eigenvector_occ, n_basis,temp_C1, n_basis,&
                   zero,temp_complex, n_basis)
       temp_complex = temp_complex + temp_1_complex

    !-------this is slow loop version, to be removed.------------
    ! do i_basis=1,n_basis
    ! do j_basis=1,n_basis
    !     temp_complex(i_basis,j_basis)=(0.0d0,0.0d0)

    !     do i_state=1, max_occ_number(1,i_k_task)
    !     do j_state=1,n_states

    !      temp_complex(i_basis,j_basis)= temp_complex(i_basis,j_basis) &
    !       + occ_numbers(i_state,1,i_k_point)          & 
    !       * KS_eigenvector_complex(i_basis,j_state,1,i_k_task)    &
    !       * first_order_U_complex( j_state,i_state,i_k_task)       &
    !       * dconjg(KS_eigenvector_complex(j_basis,i_state,1,i_k_task ))           &
    !       + occ_numbers(i_state,1,i_k_point)          & 
    !       * KS_eigenvector_complex(i_basis,i_state,1,i_k_task)   &
    !       * dconjg(KS_eigenvector_complex(j_basis,j_state,1,i_k_task ))           &
    !       * dconjg(first_order_U_complex( j_state,i_state,i_k_task))

    !     enddo
    !     enddo

    ! enddo
    ! enddo
    !===================end complex_eigenvectors=========================
    endif !real/complex


    !-----------------add k_phase ----------------
     do i_basis = 1, n_basis
     do i_cell = 1,n_cells_in_hamiltonian-1

             if (index_hamiltonian(1,i_cell, i_basis) > 0) then
             do i_index = index_hamiltonian(1, i_cell, i_basis), &
                &         index_hamiltonian(2, i_cell, i_basis)
                   j_basis =  column_index_hamiltonian(i_index)

                if (real_eigenvectors) then

                  first_order_density_matrix_sparse(i_index)= &
                  first_order_density_matrix_sparse(i_index)  &
                           + temp_real(j_basis,i_basis)  &
                           * dble( k_phase(i_cell,i_k_point) ) &  
                           * k_weights(i_k_point)
                else

                  first_order_density_matrix_sparse(i_index)= &
                  first_order_density_matrix_sparse(i_index)  &
                          + dble( (temp_complex(j_basis,i_basis) ) &
                           * dconjg(k_phase(i_cell,i_k_point) ) ) &  
                           * k_weights(i_k_point)
                endif

             end do
             end if

     end do ! i_cell 
     end do ! i_basis

    endif ! i_k_task
    enddo ! i_k_point




  !-------shanghui begin parallel------
   call sync_sparse_matrix(first_order_density_matrix_sparse) 
  !-------shanghui end parallel------
 
  if(real_eigenvectors) then 
    deallocate(temp_first_order) 
    deallocate(temp_1) 
    deallocate(temp_2) 
    deallocate(KS_eigen_times_occnumber) 
    deallocate(temp_real) 
  else 
    deallocate(temp_1_complex)
    deallocate(temp_eigenvector)
    deallocate(temp_eigenvector_occ)
    deallocate(temp_C1)
    deallocate(temp_complex)
  endif

 endif


end subroutine evaluate_first_order_DM_dielectric
!******
