!****s* FHI-aims/analy_continute_self_energy
!  NAME
!   analy_continute_self_energy
!  SYNOPSIS
!  
      subroutine analy_continue_self_energy_dmft &
          (anacon_type,&
           n_low_state,n_high_state,n_freq, &
           n_max_par, &
           chemical_potential, &
           sigma_par,omega, &
           self_energy_omega,i_k_point)

!  PURPOSE
!  Subroutine analy_continue_self_energy  analytically continue 
!  the self energy from imaginary energy axis to real energy axis.
!  Two schemes can be used here: the two-pole fitting or Pade approximation 
!
! USES

      use dimensions
      use localorb_io, only: use_unit
      use mpi_tasks

      implicit none

! ARGUMENTS

      integer ::  anacon_type
      integer ::  n_low_state
      integer ::  n_high_state
      integer ::  n_freq
      integer ::  n_max_par 
      integer ::  i_k_point 

      real*8 ::   omega(n_freq)
      real*8 ::   chemical_potential
      complex*16 ::   self_energy_omega  &
               (n_low_state:n_high_state,n_low_state:n_high_state, n_freq )
               !(n_basis,n_basis, n_freq )

      complex*16 :: sigma_par(n_max_par,n_states,n_states)
      !complex*16 :: sigma_par(n_max_par,n_basis,n_basis)

! INPUTS 
!         self-energy comes from the calling subroutine
!  o  anacon_type -- integer number, if 0, the two-pole fitting for analytical
!          continuation; if 1, using Pade approximation for ana. cont.
!  o  n_low_state  -- integer number,
!          the lowest KS/HF eigenstate for self-energy correction
!  o  n_high_state -- integer number,
!          the highest KS/HF eigenstate for self-energy correction
!  o  n_freq -- integer number, the number of frequency points for the GW self-energy
!  o  n_max_par -- the number of parameters used for analytical continuation
!          For anacon_type = 0, recommended n_max_par is  4 or 5. If 4, this will
!          be the normal two-pole fitting, else if 5, it will be two-pole plus a
!          (small) constant number
!          For anacon_type = 1, recommended n_max_par is the half of n_freq
!  o  omega(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the self-energy
!  o  chemical_potential -- real number, the chemical potential of the system
!  o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
!  o  self_energy_omega  -- complex array, the calculated self-energy on imaginary
!            frequency grid for each state and each spin channel
! OUTPUT
!  o  sigma_par -- complex array, the fitting parameters from analytical continuation

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
      integer n_freq_tmp
      integer number_of_states

      real*8 en
      complex*16 selfe, dselfe
      
!  fitting omega grid range 
      
      complex*16, allocatable ::  self_energy_imgy_omega(:,:)
      complex*16, allocatable ::  self_energy_real_omega(:,:)

      complex*16, allocatable :: par(:)
      

!     counters


      integer :: i_state
      integer :: j_state
      integer :: i_KS_state
      integer :: i_freq
      integer :: i_spin

!   external function  

      integer determine_indx

!      complex*16 sig_f 

!     begin work

      

!   write out 
    
      if(myid.eq.0) then
        write(use_unit,*) 
        write(use_unit,'(2X,A)') "Analytical continuation for GW self-energy starts..."
        write(use_unit,*) "                  for k-point #", i_k_point
      endif 

      !number_of_states = n_high_state-n_low_state+1        
      number_of_states = n_high_state        

      allocate(self_energy_imgy_omega(n_freq,number_of_states))
!      call check_allocation(i_spin, 'self_energy_imgy_omega        ')


!      allocate (self_energy_real_omega(-n_freq+1:n_freq, &
!                                       number_of_states,n_spin)) 

!        do i_spin = 1, n_spin, 1
         do i_freq = 1, n_freq, 1
          do i_state = 1, number_of_states, 1
             i_KS_state= i_state + n_low_state-1
             self_energy_imgy_omega(i_freq,i_state) = &
             self_energy_omega(i_KS_state,i_KS_state,i_freq)
          enddo
         enddo
!        enddo
      
      allocate (par(n_max_par))
!      call check_allocation(i_spin, 'par                           ')

! perform analytical continuation, and get the fitting parameters 
!      do i_spin = 1, n_spin, 1
!         write(use_unit,*) " | i_spin   ", i_spin
         !do i_state = n_low_state, number_of_states, 1
         do i_state = 1,n_high_state, 1
          do j_state = 1,n_high_state, 1
!             write(use_unit,*) " |      i_state   ", i_state

            i_KS_state = i_state+n_low_state-1

            if(anacon_type .eq. 0) then
             call mr_min_lsq(n_freq,dcmplx(0.d0,omega), & 
                        self_energy_imgy_omega(:,i_state), &
                        !self_energy_omega(i_state,j_state,:), &
                        n_max_par,par) 
            elseif(anacon_type .eq. 1) then
             call get_pade_approx_para(n_freq,omega, &
                        self_energy_imgy_omega(:,i_state), &
                        n_max_par,par)
            endif

            sigma_par(1:n_max_par,i_state,j_state)=par(1:n_max_par)



         enddo
       enddo


      deallocate (par)
      deallocate (self_energy_imgy_omega)
!       deallocate (self_energy_real_omega)


      return 
      end subroutine analy_continue_self_energy_dmft
!******
!---------------------------------------------------------------------
