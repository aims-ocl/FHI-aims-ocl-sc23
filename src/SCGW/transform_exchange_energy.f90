      subroutine transform_exchange_energy &
         (overlap_matrix, n_high_state,aux_exchange_self_energy,&
           exchange_self_energy)

! MODULES
      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use poles_fit
      use localorb_io, only: use_unit

      implicit none
!ARGUMENTS
       integer n_high_state
       real*8  overlap_matrix (n_states, n_basis)
       real*8  aux_exchange_self_energy(n_high_state,n_spin)
!OUTPUT
       real*8  exchange_self_energy(n_basis,n_basis)
!counters
       integer i_state
       integer i_basis
       integer j_basis
       integer i_index

       exchange_self_energy(:,:) = 0.d0

       write(use_unit,*) "   --- Transforming the exchange self energy in the", &
                          " NAO basis"
       write(use_unit,*) "   "

!       do i_spin =1, n_spin,1 
         do i_state=1, n_states, 1
          do i_basis= 1, n_basis, 1
           do j_basis= 1, n_basis, 1

              exchange_self_energy (i_basis, j_basis) = &
              exchange_self_energy (i_basis, j_basis) + &
              overlap_matrix (i_state, i_basis) * &
              overlap_matrix (i_state, j_basis) * &
              aux_exchange_self_energy(i_state,1)

!                if(.not.(i_basis .eq. j_basis) ) then
!                  exchange_self_energy (j_basis, i_basis) = &
!                  exchange_self_energy (i_basis, j_basis)
!                endif
           enddo
          enddo
         enddo 
!       enddo

        open(44,file="exchange_energy_KS.dat")
        do i_state = 1, n_states, 1
          write(44,*) i_state,  aux_exchange_self_energy(i_state, 1)
        enddo
        close(44)


        open(44,file="exchange_energy_old.dat")
        i_index= 0
        do i_basis = 1, n_basis, 1
         do j_basis = 1, n_basis, 1
          i_index= i_index+1
          write(44,*) i_index,  exchange_self_energy(i_basis, j_basis)
         enddo
        enddo
        close(44)


      end subroutine transform_exchange_energy
