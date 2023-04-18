!****s* FHI-aims/evaluate_zero_order_EDM_phonon_reduce_memory
!  NAME
!    evaluate_zero_order_EDM_phonon_reduce_memory
!  SYNOPSIS

subroutine evaluate_zero_order_EDM_phonon_reduce_memory &
           (KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers,energy_density_matrix_sparse)
               

!  PURPOSE
!    calculate the zero-order EDM for PBC case 

!  DM(uR1, vR2)(0)=sum{k} sum{i} Cui*(k) Cvi(k) exp(-ik(R2-R1))
!  shanghui,2013.12.30
!  USES

  use dimensions
  use pbc_lists
  use synchronize_mpi
  use scalapack_wrapper
  use physics, only : overlap_matrix, hamiltonian, n_electrons, &
                      chemical_potential, chemical_potential_spin
  use runtime_choices

!  ARGUMENTS

  implicit none
  
  real*8,     dimension(n_basis, n_states, n_spin,n_k_points_task),  intent(IN) :: KS_eigenvector
  complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task),  intent(IN) :: KS_eigenvector_complex
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers

  real*8, dimension(n_hamiltonian_matrix_size),intent(INOUT) :: energy_density_matrix_sparse

!  INPUTS
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates

!
!  OUTPUT
!  density_matrix


  integer :: i_basis,j_basis ,i_coord, i_atom,i_spin,i_index,i_cell
  integer :: i_state,j_state,i_k_point, i_k_task
  integer, dimension(n_spin, n_k_points_task)   :: max_occ_number

  real*8 :: occu
  integer :: n_max_occupied,info

  real*8, dimension(n_basis,n_basis) :: temp
  real*8,  dimension(n_basis, n_states) :: KS_scaled
  complex*16 , dimension(n_basis, n_basis)  :: temp_cmplx
  complex*16,  dimension(n_basis, n_states) :: KS_scaled_cmplx

 
  real*8, dimension(n_states, n_spin, n_k_points) ::  occ_numbers_KS_eigenvalue


  energy_density_matrix_sparse(1:n_hamiltonian_matrix_size) = 0.0d0

  if(use_scalapack)then

      do i_spin = 1, n_spin 
          !------note(1):occ_numbers and KS_eigenvalue now have all k_points information----------
          !------note(2):in fact we need occ_numbers*KS_eigenvalue, which is a negtive value
          !              but in orde to use it in scalapack_wrapper which have sqrt operater,
          !              so here we use - occ_numbers_KS_eigenvalue, and we use change back 
          !              after sync_density_matrix_sparse.  
          do i_state = 1, n_states
          do i_k_point = 1,n_k_points 
                occ_numbers_KS_eigenvalue(i_state,i_spin,i_k_point) =  & 
               -occ_numbers(i_state,i_spin,i_k_point) * &  
                KS_eigenvalue(i_state,i_spin,i_k_point)    
          enddo
          enddo 

         call construct_dm_scalapack(occ_numbers_KS_eigenvalue, i_spin)

         call get_sparse_matrix_scalapack(energy_density_matrix_sparse,i_spin)

         call sync_density_matrix_sparse(energy_density_matrix_sparse)
              energy_density_matrix_sparse = - energy_density_matrix_sparse ! see note(2) above.
      enddo 

  else 
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

      
       if (real_eigenvectors) then
          temp = 0.0d0
       else 
          temp_cmplx = (0.d0, 0.d0)
       endif


!       n_max_occupied = 1
!       do i_state = 1, n_states, 1
!          occu =  occ_numbers(i_state, 1, i_k_point)
!          if (occu .gt. 0.d0) then
!             n_max_occupied = i_state
!             if (real_eigenvectors) then
!                KS_scaled(:, i_state) = KS_eigenvector(:, i_state, 1, i_k_task) * sqrt(occu)
!             else 
!                KS_scaled_cmplx(:, i_state) = KS_eigenvector_complex(:, i_state, 1, i_k_task) &
!                                             * sqrt(occu)
!             endif
!          endif
!       end do

!---------------------------------(1) use lapacek------------------------------------
!      !-----here DM = Cui^(*) Cvi = Ci_aims * Ci_aims^(diagger)  the same as PM_none-----  
!       if (real_eigenvectors) then
!          call dsyrk('U', 'N', n_basis, n_max_occupied, 1.d0, KS_scaled, n_basis, &
!          &          0.d0, temp, n_basis)
!       else 
!          call zherk('U', 'N', n_basis, n_max_occupied, 1.d0, KS_scaled_cmplx, n_basis, &
!          &          0.d0, temp_cmplx, n_basis)
!       endif
! 
!      !----zherk only give upder part, so we need add the down park
!       if (real_eigenvectors) then
!          temp = temp + transpose(temp)
!          do i_basis = 1, n_basis
!             temp(i_basis, i_basis) = temp(i_basis, i_basis)*0.5d0
!          end do
!       else 
!          temp_cmplx = temp_cmplx + transpose(dconjg(temp_cmplx))
!          do i_basis = 1, n_basis
!             temp_cmplx(i_basis, i_basis) = temp_cmplx(i_basis, i_basis)*0.5d0
!          end do
!       endif
!---------------------------------(1) end use lapacek------------------------------------


!---------------------------------(2) use loop------------------------------------
!     !-----here DM = Cui^(*) Cvi= Ciu_aims * Cvi_aims^(*)-------------- 
     do i_basis=1,n_basis
       do j_basis=1,n_basis
        do i_state=1,  max_occ_number(1,i_k_task) 
         
          if(real_eigenvectors) then 
            temp(i_basis,j_basis)=temp(i_basis,j_basis)       &
          + KS_eigenvector(i_basis,i_state,1,i_k_task) & 
          * occ_numbers(i_state,1,i_k_point)*KS_eigenvalue(i_state,1,i_k_point) &  
          * KS_eigenvector(j_basis,i_state,1,i_k_task )              
          else
            temp_cmplx(i_basis,j_basis)=temp_cmplx(i_basis,j_basis)       &
          + KS_eigenvector_complex(i_basis,i_state,1,i_k_task) & 
          * occ_numbers(i_state,1,i_k_point)*KS_eigenvalue(i_state,1,i_k_point) &  
          * dconjg(KS_eigenvector_complex(j_basis,i_state,1,i_k_task ))                   
          endif 

       enddo
       enddo
     enddo  
!----------------------------(2) end loop---------------------------------------




  

      do i_basis = 1, n_basis
          do i_cell = 1,n_cells_in_hamiltonian-1
             if (index_hamiltonian(1,i_cell, i_basis) > 0) then
                do i_index = index_hamiltonian(1, i_cell, i_basis), &
                &            index_hamiltonian(2, i_cell, i_basis)
                   j_basis =  column_index_hamiltonian(i_index)

                if (real_eigenvectors) then
                  energy_density_matrix_sparse(i_index) &
                  & = energy_density_matrix_sparse(i_index) &
                  & + temp(i_basis, j_basis) *dble(k_phase(i_cell,i_k_point)) &
                    * k_weights(i_k_point)
                else
                  energy_density_matrix_sparse(i_index) &
                  & = energy_density_matrix_sparse(i_index) &
                  & + dble(temp_cmplx(i_basis, j_basis) *k_phase(i_cell,i_k_point)) &
                    * k_weights(i_k_point)
                  ! DM(uR,v) = temp(u,v)* exp(-ik(-R))
                endif

                

                end do
             end if
          end do
       end do




      endif 
    enddo ! i_k_point

    call sync_sparse_matrix(energy_density_matrix_sparse) 

    endif



end subroutine evaluate_zero_order_EDM_phonon_reduce_memory
!******
