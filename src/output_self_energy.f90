!****s* FHI-aims/output_self_energy
!  NAME
!   output_self_energy
!  SYNOPSIS

subroutine output_self_energy &
        (n_freq,n_low_state,n_high_state, &
         omega,self_energy_freq)

!  PURPOSE
!  Print out GW or beyond-GW self energies
!  

!  USES

      use dimensions
      use runtime_choices
      use mpi_tasks

      implicit none

!  ARGUMENTS

      integer  :: n_freq
      integer  :: n_low_state
      integer  :: n_high_state
      real*8   :: omega(n_freq)
      complex*16 ::  self_energy_freq &
            (n_low_state:n_high_state, n_freq, n_spin )

! INPUTS
! o  n_freq -- integer number, the number of frequency points for the GW self-energy
! o  n_low_state  -- integer number,
!            the lowest KS/HF eigenstate for self-energy correction
! o  n_high_state -- integer number,
!            the highest KS/HF eigenstate for self-energy correction
! o  self_energy_freq  -- complex array, 
!           the calculated GW self-energy on the imaginary axis
! OUTPUTS
! o  none
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!  local variables

      character*50  filename 
      
!     counters


      integer :: i_state
      integer :: j_state

      integer :: i_freq

      integer :: i_spin

!     begin work

!    printing out  
       if (myid.eq.0) then
        do i_spin = 1, n_spin
         do j_state =  n_low_state, n_high_state, 1

          if(n_spin.eq.2) then
           write(filename,'(A,I0,A,I0,A)')"self_energy/Sigma.omega.n_", &
                        j_state,".s_",i_spin,".dat"
           else
            write(filename,'(A,I0,A)')"self_energy/Sigma.omega.n_", &
                        j_state,".dat"
          endif

           open(102, file=filename)

           do i_freq = 1, n_freq, 1

              write (102, '(2X, f16.8, 7X, 2F18.10)') &  
                 omega(i_freq), &
                 self_energy_freq(j_state,i_freq,i_spin) 
            enddo
           
            close(102)
        enddo
       enddo

       open (103, file='self_energy_omega.dat')
       write(103,*) n_freq, n_high_state - n_low_state +1
       write(103,*) (omega(i_freq),i_freq=1,n_freq)
       do i_spin = 1, n_spin 
       do i_state=n_low_state,n_high_state
          write(103,*)i_state
          write(103,*) &
             (self_energy_freq(i_state,i_freq,i_spin), &
                    i_freq=1,n_freq)
       enddo
       enddo 
       close(103)
! end of if output_self_energy
      endif

      return 
      end subroutine output_self_energy
