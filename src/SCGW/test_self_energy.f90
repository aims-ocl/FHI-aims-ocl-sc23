      subroutine test_self_energy  (self_energy_freq, &
              n_homo,nomega,omega, KS_eigenvector, &
             self_energy_freq_KS )
! this trasform the self_energy from the NAO basis to the
! KS basis

!  USES

      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi

!ARGUMENTS
      implicit none
      integer n_homo (n_spin)
      integer nomega
      real*8  omega(nomega)
!      real*8  ovlp_NAO_KS(n_states,n_basis, n_spin)
      complex*16 self_energy_freq (n_basis,n_basis, nomega, n_spin) 
      complex*16 self_energy_freq_KS (n_states, nomega, n_spin)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin)
!      complex*16  self_energy_omega (n_low_state:n_high_state, nomega, n_spin)
! OTHERS
      complex*16 aux_self_energy_freq_KS (n_states, nomega)
! counters
      integer i_state
      integer i_freq
      integer i_spin
      integer i_basis
      integer j_basis
      integer counter

      self_energy_freq_KS (:,:,:)=0.d0

      do i_spin = 1, n_spin, 1
        do i_basis=1, n_basis, 1
          do j_basis=1, n_basis, 1
            do i_state=1, n_states, 1      

               self_energy_freq_KS (i_state,:,i_spin) =&
               self_energy_freq_KS (i_state,:,i_spin) +&
               KS_eigenvector (i_basis,i_state,i_spin)*&
               self_energy_freq (i_basis, j_basis, :, i_spin)*&
               KS_eigenvector (j_basis,i_state,i_spin)  
 
            enddo 
          enddo
        enddo
      enddo

      if (myid.eq.0.and..false.)then
        open (123,file="self_energy_KS.dat")
        do i_freq=1, nomega, 1
!OCCUPIED STATES
        write(123,*) omega(i_freq), real(self_energy_freq_KS (n_homo,i_freq,1)),&
                 aimag(self_energy_freq_KS (n_homo,i_freq,1))
!        write(123,*) omega(i_freq), real(self_energy_freq_KS (n_homo-1,i_freq)),&
!                 aimag(self_energy_freq_KS (n_homo-1,i_freq))
!        write(123,*) omega(i_freq), real(self_energy_freq_KS (n_homo-2,i_freq)),&
!                 aimag(self_energy_freq_KS (n_homo-2,i_freq))
!EMPTY STATES
          write(123,*) omega(i_freq), real(self_energy_freq_KS (n_homo+1,i_freq,1)),&
                           aimag(self_energy_freq_KS (n_homo+1,i_freq,1))
          write(123,*) omega(i_freq), real(self_energy_freq_KS (n_homo+2,i_freq,1)),&
                           aimag(self_energy_freq_KS (n_homo+2,i_freq,1))
          write(123,*) omega(i_freq), real(self_energy_freq_KS (n_homo+3,i_freq,1)),&
                           aimag(self_energy_freq_KS (n_homo+3,i_freq,1))  

!        write(123,*) omega(i_freq), real(self_energy_freq_KS (1,i_freq)),&
!                              aimag(self_energy_freq_KS (1,i_freq))
!        write(123,*) omega(i_freq), real(self_energy_freq_KS (n_homo,i_freq)),& 
!                              aimag(self_energy_freq_KS (n_homo,i_freq))
!        write(123,*) omega(i_freq), real(self_energy_freq_KS (n_states,i_freq)),& 
!                              aimag(self_energy_freq_KS (n_states,i_freq))
        enddo
        close(123)
      endif

      end subroutine test_self_energy

