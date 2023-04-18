      subroutine test_exchange_self_energy  (self_energy_freq, &
              n_homo, KS_eigenvector, self_energy_freq_KS)
! this trasform the self_energy from the NAO basis to the
! KS basis

!  USES

      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      
      implicit none
!ARGUMENTS
      integer n_homo(n_spin)
!      real*8  ovlp_NAO_KS(n_states,n_basis,n_spin)
      real*8 self_energy_freq (n_basis,n_basis,n_spin) 
      real*8 self_energy_freq_KS (n_states, n_spin)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin)
      
! OTHERS
      complex*16 aux_self_energy_freq_KS (n_states)

      
! counters
      integer i_state
      integer i_freq
      integer i_basis
      integer j_basis
      integer counter
      integer i_spin

      self_energy_freq_KS (:,:)=0.d0

      do i_spin = 1, n_spin
       do i_basis=1, n_basis, 1
         do j_basis=1, n_basis, 1
           do i_state=1, n_states, 1      

               self_energy_freq_KS (i_state,i_spin ) =&
               self_energy_freq_KS (i_state,i_spin) +&
               KS_eigenvector (i_basis,i_state,i_spin)*&
               self_energy_freq (i_basis, j_basis,i_spin)*&
               KS_eigenvector (j_basis,i_state,i_spin)  
 
           enddo 
         enddo
       enddo
      enddo

      if (myid.eq.0 .and. .false.)then
        open (123,file="exchange_self_energy_KS.dat")
        do i_state = 1, n_states, 1
          write(123,*) i_state, self_energy_freq_KS (i_state, 1)
        enddo
        close(123)
      endif

      end subroutine test_exchange_self_energy

