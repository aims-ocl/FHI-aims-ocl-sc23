!****s* FHI-aims/analy_continute_self_energy
!  NAME
!   analy_continute_self_energy
!  SYNOPSIS
!  
      subroutine analy_continue_self_energy &
          (flag_qpe, anacon_type,&
           n_low_state,n_high_state,n_freq, &
           n_max_par, &
           KS_eigenvalue,chemical_potential, &
           sigma_par,omega, &
           self_energy_omega)

!  PURPOSE
!  Subroutine analy_continue_self_energy  analytically continue 
!  the self energy from imaginary energy axis to real energy axis.
!  Two schemes can be used here: the two-pole fitting or Pade approximation 
!
! USES

      use dimensions
      use mpi_tasks
      use localorb_io,             only :use_unit
      use synchronize_mpi,         only: sync_vector_complex
      use gw_para,                 only: anacon_type_defined
      implicit none

! ARGUMENTS

      logical ::  flag_qpe
      integer ::  anacon_type
      integer ::  n_low_state
      integer ::  n_high_state
      integer ::  n_freq
      integer ::  n_max_par 

      real*8 ::   omega(n_freq)
      real*8 ::   chemical_potential
      real*8 ::   KS_eigenvalue(n_states,n_spin)
      complex*16 ::   self_energy_omega  &
               (n_low_state:n_high_state, n_freq,n_spin )

      complex*16 :: sigma_par(n_max_par,n_states,n_spin)

! INPUTS 
!  o  flag_qpe -- logical, if true, read the self-energy from files, if false, the
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
      
      complex*16, allocatable ::  self_energy_imgy_omega(:,:,:)
      complex*16, allocatable ::  self_energy_real_omega(:,:,:)

      complex*16, allocatable :: par(:)
      

!     counters


      integer :: i_state
      integer :: i_KS_state
      integer :: i_freq
      integer :: i_spin
!   test
!      integer :: i_par
!   end test

!   external function  

      integer determine_indx

!      complex*16 sig_f 

!     begin work

      

!   write out 

      ! Error trap in case the analytical continuation type was not defined in control.in.
      if (.not.anacon_type_defined) then
         if (myid.eq.0) then
            write(use_unit,'(2X,A)') "Analytical continuation of the self-energy starts..."
            write(use_unit,'(1X,A)') "* Oops - error. The analytical continuation of the self-energy"
            write(use_unit,'(1X,A)') "* from the imaginary to the real frequency axis was invoked"
            write(use_unit,'(1X,A)') "* although the keyword 'anacon_type' was not defined in control.in."
            write(use_unit,'(1X,A)') "* This should not happen (this omission should have been caught when"
            write(use_unit,'(1X,A)') "* first reading control.in in read_control.f90). Please notify us"
            write(use_unit,'(1X,A)') "* if you do encounter this error message - normally, it should"
            write(use_unit,'(1X,A)') "* never trigger."
            call aims_stop('anacon_type not defined in control.in.','analy_continue_self_energy.f90')
         end if
      end if
      
      if(myid.eq.0) then
        write(use_unit,*) 
        if(anacon_type .eq. 0) then
          write(use_unit,'(2X,A)') "Analytical continuation (two-pole fitting) starts..."
        elseif (anacon_type .eq. 1) then
          write(use_unit,'(2X,A)') "Analytical continuation (Pade approximation) starts..."
        endif
      endif 

      number_of_states = n_high_state-n_low_state+1        

      allocate(self_energy_imgy_omega(n_freq,number_of_states,n_spin),stat=i_spin)
      call check_allocation(i_spin, 'self_energy_imgy_omega        ')
      self_energy_imgy_omega = 0.0d0

      !allocate (self_energy_real_omega(-n_freq+1:n_freq, &
      !                                 number_of_states,n_spin)) 

      if(flag_qpe) then 
        if(myid == 0) then
          open (104, file='self_energy_omega.dat')
          read(104,*) n_freq_tmp, number_of_states
          
          if(n_freq_tmp.ne.n_freq) then
             write(use_unit,*) "Number of frequencies is wrong! Stop!"
          endif 

          read(104, *) (omega(i_freq), i_freq =1, n_freq)
          do i_spin = 1, n_spin, 1
             do i_state = 1, number_of_states
               read(104,*) 
               read(104, *) &
                   (self_energy_imgy_omega(i_freq,i_state,i_spin), &
                        i_freq=1,n_freq)
            enddo
          enddo
          close(104)
        endif
        call sync_vector_complex(self_energy_imgy_omega,n_freq*number_of_states*n_spin)
    
      else  
        do i_spin = 1, n_spin, 1
         do i_freq = 1, n_freq, 1
          do i_state = 1, number_of_states, 1
             i_KS_state= i_state + n_low_state-1
             self_energy_imgy_omega(i_freq,i_state,i_spin) = &
             self_energy_omega(i_KS_state,i_freq,i_spin)
          enddo
         enddo
        enddo
!  end of if flag_qpe
      endif
      
      allocate (par(n_max_par),stat=i_spin)
      call check_allocation(i_spin, 'par                           ')


! perform analytical continuation, and get the fitting parameters 
      do i_spin = 1, n_spin, 1
!         write(use_unit,*) " | i_spin   ", i_spin
         do i_state = 1, number_of_states, 1 
!             write(use_unit,*) " |      i_state   ", i_state

            i_KS_state = i_state+n_low_state-1
            en = KS_eigenvalue(i_KS_state,i_spin)-chemical_potential

            if(anacon_type .eq. 0) then
             call mr_min_lsq(n_freq,dcmplx(0.d0,omega), & 
                        self_energy_imgy_omega(:,i_state,i_spin), &
                        n_max_par,par) 
            elseif(anacon_type .eq. 1) then
             call get_pade_approx_para(n_freq,omega, &
                        self_energy_imgy_omega(:,i_state,i_spin), &
                        n_max_par,par)
            endif

            sigma_par(1:n_max_par,i_state,i_spin)=par(1:n_max_par)
! test
!            write(use_unit,*)"i_state =", i_state
!!            if(i_state .eq. 20) then
!               write(use_unit,*)myid, i_spin, i_state 
!               do i_par = 1, n_max_par, 1
!                  write(use_unit,'(I4,2f24.15)') i_par, sigma_par(i_par,i_state,i_spin)
!               enddo
!            endif
! test
           !do i_freq = n_freq, 1, -1 

           !  en=-1.d1+(n_freq-i_freq)*0.25d0
           !  call get_real_selfenergy(anacon_type,n_freq,omega, &
           !           dcmplx(en,0.d0), n_max_par, &
           !           sigma_par(1:n_max_par,i_state,i_spin),selfe)
!          !   call get_real_selfenergy(anacon_type,n_freq,omega, &
!          !            dcmplx(-omega(i_freq),0.d0), n_max_par, &
!          !            sigma_par(1:n_max_par,i_state,i_spin),selfe)
 
           !  self_energy_real_omega(-i_freq+1,i_state,i_spin) = selfe
           !enddo

           !do i_freq = 1, n_freq, 1 

           ! en=(i_freq-1)*0.25d0
           ! call get_real_selfenergy(anacon_type,n_freq,omega,&
           !           dcmplx(en,0.d0),n_max_par,&
           !           sigma_par(1:n_max_par,i_state,i_spin),selfe)
!          !  call get_real_selfenergy(anacon_type,n_freq,omega,&
!          !            dcmplx(omega(i_freq),0.d0),n_max_par,&
!          !            sigma_par(1:n_max_par,i_state,i_spin),selfe)

!          !  call get_real_selfenergy(anacon_type,n_freq,omega,&
!          !            dcmplx(0.d0,omega(i_freq)),n_max_par,&
!          !            sigma_par(1:n_max_par,i_state,i_spin),selfe)
 
           !  self_energy_real_omega(i_freq,i_state,i_spin) = selfe
           !enddo

        enddo
      enddo

      !if(myid.eq.0) then 
      !open (103, file= 'selfenergy_real_omega_new.dat') 
!     ! open (104, file= 'selfenergy_imgy_omega_new.dat') 
      !
      !     do i_freq = n_freq, 1, -1

      !       en=-1.d1+(n_freq-i_freq)*0.25d0
      !       write (103, '(2X, 200F18.10)') &
      !             en, &
      !            ((self_energy_real_omega(-i_freq+1,i_state,i_spin), &
      !             i_spin=1,n_spin),i_state =1, 30)   
!
!     !        write (104, '(2X, 200F18.10)') &
!     !             -omega(i_freq), &
!     !             (( conjg(self_energy_imgy_omega(i_freq,i_state, &
!     !               i_spin)), i_spin=1,n_spin), i_state =1, 10)   
      !    enddo

      !    do i_freq = 1, n_freq, 1

      !        en=(i_freq-1)*0.25
      !        write (103, '(2X, 200F18.10)') &
      !           en, &
!     !            omega(i_freq), &
      !           ((self_energy_real_omega(i_freq,i_state,i_spin), &
      !           i_spin=1,n_spin),i_state = 1, 30)

!     !        write (104, '(2X, 200F18.10)') &
!     !             omega(i_freq), &
!     !             ((self_energy_imgy_omega(i_freq,i_state,i_spin), &
!     !             i_spin=1,n_spin),i_state =1, 10)   
      !    enddo

      ! close(103)
!     !   close(104)
      ! endif

      deallocate (par)
      deallocate (self_energy_imgy_omega)
!       deallocate (self_energy_real_omega)


      return 
      end subroutine analy_continue_self_energy
!******
!---------------------------------------------------------------------
