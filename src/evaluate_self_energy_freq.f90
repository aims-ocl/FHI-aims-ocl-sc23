!****h* FHI-aims/evaluate_self_energies
!  NAME
!   Routines for evaluating self energies
!  SYNOPSIS

module evaluate_self_energy_freq

   use dimensions
   use runtime_choices
   use constants
   use prodbas        
   use mpi_tasks            
   use synchronize_mpi,              only: sync_timing
   use localorb_io,                  only: use_unit
   use contour_def_gw_types,         only: cd_environment_type
   use contour_def_gw_environment,   only: contour_def_integrate,&
                                           contour_def_add_residues
   use matrix_inversion,             only: invert_real_matrix,&
                                           invert_cmplx_matrix,&
                                           set_matrix_elements
   use evaluate_polarisability_freq, only: evaluate_polarisability_freq_1, &
                                           evaluate_polarisability_freq_2
   use synchronize_mpi,              only: sync_vector, &
                                           sync_vector_complex,&
                                           sync_complex_number,&
                                           sync_matrix
   use scalapack_matrix_type,        only: matrix_structure_type,&
                                           create_matrix_structure
   use timing,                       only: time_WPQ_realfreq_cd,&
                                           clock_time_WPQ_realfreq_cd,&
                                           time_Wmn_realfreq_cd,&
                                           clock_time_Wmn_realfreq_cd,&
                                           time_polar_realfreq_cd,&
                                           clock_time_polar_realfreq_cd,&
                                           time_self_energy_cd,&
                                           clock_time_self_energy_cd,&
                                           time_polar_imagfreq,&
                                           clock_time_polar_imagfreq,&
                                           time_Wmn_imagfreq, clock_time_Wmn_imagfreq,&
                                           get_timestamps,&
                                           get_times

   implicit none

   private

   public :: evaluate_self_energy_freq_1, evaluate_self_energy_freq_2,&
             evaluate_self_energy_cd, evaluate_Wmn_imag_freq_for_cd
             
 
contains

! **************************************************************************************************

  subroutine evaluate_self_energy_freq_1 &
            (n_low_state, n_high_state, &
            n_homo,  &
            occ_numbers,  &
            n_freq,n_full_freq, &
            omega, womega, &
            omega_full, womega_full, &
            chemical_potential_spin, &
            KS_eigenvalue, KS_eigenvector, &
            partition_tab,  basis_l_max, &
            ovlp_3KS, &
            self_energy_freq, &
            flag_coh, tot_time_sigma, tot_time_polar, &
            KS_eigenvalue_scf &
            )

!  PURPOSE
!  Subroutine evaluate_self_energy_freq evaluates the correlated part of 
!  the GW self-nergy on the imaginary frequency axis. 
!  
!  Sigma^_ii(iw) = i \sum_lmn O_ilm O_iln \int G_ll(iw+iw') W _mn(iw') dw'

!  USES

!  ARGUMENTS

      integer :: n_freq
      integer :: n_full_freq
      integer :: n_low_state
      integer :: n_high_state
      integer :: basis_l_max (n_species)
      integer :: n_homo(n_spin)

      real*8  :: occ_numbers(n_states,n_spin)
      real*8  :: omega(n_freq)
      real*8  :: omega_full(n_full_freq)
      real*8  :: womega(n_freq)
      real*8  :: womega_full(n_full_freq)
      real*8  :: chemical_potential_spin(n_spin)   
      real*8  :: KS_eigenvalue(n_states,n_spin)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin)
      real*8  :: ovlp_3KS(n_loc_prodbas, n_states, n_high_state,n_spin)
      real*8  :: partition_tab(n_max_angular, n_max_radial, n_atoms)

      logical flag_coh

      real*8  tot_time_polar
      real*8  tot_time_sigma

      complex*16 ::  self_energy_freq &
            ( n_low_state:n_high_state, n_freq, n_spin )
      real*8, intent(in), optional        :: KS_eigenvalue_scf(n_states,n_spin)

! INPUTS
! o  n_freq -- integer number, the number of frequency points for the GW self-energy
! o  n_full_freq -- integer number,
!            the number of frequency points for the screened Coulomb interaction W
! o  n_low_state  -- integer number,
!            the lowest KS/HF eigenstate for self-energy correction
! o  n_high_state -- integer number,
!            the highest KS/HF eigenstate for self-energy correction
! o  basis_l_max -- integer array of length n_species, 
!            the maximal angular momentum for the basis functions for each species
! o  n_homo -- integer array of length of n_spin,
!            the HOMO level for each spin channel
! o  n_lumo -- integer array of length of n_spin,
!            the LUMO level for each spin channel
! o  n_electrons -- real number
!            the total number of electrons in the system
! o  occ_numbers -- real 2-dimentianal array of length (n_states, n_spin)
!            the occupation number of the electrons for each eigenstate and each spin
! o  omega(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the self-energy
! o  omega_full(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the screened Coulomb interaction
! o  womega(n_freq) -- real array
!            the weight of the Gauss-Legendre frequency grid for the self-energy
! o  womega_full(n_freq) -- real array
!            the weigth of the Gauss-Legendre frequency grid for the screened Coulomb interaction
! o  chemical_potential -- real number, the chemical potential of the system
! o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF 
!            eigenvalue
! o  KS_eigenvector -- real array,
!            the eigenvector of the single-particle calculation
! o  ovlp_3KS -- real array
!            this is the transformed 3-cener overlap integration. Now two orbitals of
!            them are KS ones, and one is the auxiliary basis.
!            Note: for parallel calculations, the auxiliary basis are distribuated
!            among the different processors. 
! o  out_self_energy -- logic, if ture, output the self-energ on the imagniary axis 
! o  flag_coh -- logic, if ture, use the Coulomb-hole accelerator 
!             (not functioning for the time being)

! OUTPUTS
! o  self_energy_freq  -- complex array, 
!           the calculated GW self-energy on the imaginary axis
! o  tot_time_polar    -- total time for calculating the polarisability
! o  tot_time_sigma    -- total time for calculating the self-energy (except for the
!                          polarisability)
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

!     auxiliary matrices for Level 3 Blas matrix multiplications
!     n_first : the first orbital which is NOT fully occupied

      logical :: get_gw_corr
      real*8  :: gw_corr
      integer :: n_first(n_spin)
      complex*16  :: zeta
      complex*16  ::  arg

      real*8, dimension(:,:), allocatable :: polar_freq
      real*8, dimension(:,:), allocatable :: screened_coulomb

      real*8, dimension(:,:), allocatable :: aux_ovlp3KS
      real*8, dimension(:,:,:,:), allocatable :: aux_coulomb_matr
      real*8, dimension(:), allocatable :: screened_exchange
      real*8, dimension(:,:), allocatable :: coulomb_hole
      real*8, dimension(:,:), allocatable :: screened_coulomb_f
      real*8, dimension(:,:,:,:), allocatable :: ovlp_3KS_redist(:,:,:,:)
      real*8  temp_matr(n_states)
      real*8, dimension(:,:), allocatable :: my_KS_eigenvalue

!     timing
      real*8  :: time_self_energy 
      real*8  :: time_polar
      real*8  :: time_start, time_end, time_product
      real*8  :: rtime

!     parameters of the fitting tails
!       real*8  s1, s2, omega_1, omega_2
!       real*8  alpha, beta, a, b

!     counters


      integer :: i_state
      integer :: j_state

      integer :: i_basbas
      integer :: j_basbas

      integer :: i_index

      integer :: i_freq
      integer :: i_freq_1

      integer :: i_spin
      integer :: n_states_loc, i_state_loc

      integer :: info
      character(*), parameter :: func = 'evaluate_self_energy_freq_1'

!     begin work

      call cpu_time(tot_time_polar)

      if(myid.eq.0) then  
        write(use_unit,*) 
        write(use_unit,*)"-------------------------------------------------" 
        write(use_unit,'(2X,A)') "Start to calculate the self energy ... "
      endif

      get_gw_corr = .false.
      if(get_gw_corr)then
         n_low_state=1
         n_high_state=n_states
      endif
 
!     determine the highest occupied orbital level
!     such complication occurs when some of the orbitals are 
!     either not fully occupied or not fully empty
      n_first(:) = 1
      do i_spin = 1, n_spin
       do i_state = 1, n_states
        if (abs(occ_numbers(i_state,i_spin)-dble(2/n_spin)) &
                        .lt.1.d-8) then
         n_first(i_spin)= i_state + 1
        endif 
       enddo 
       if(n_first(i_spin) .gt. n_states) then
         n_first(i_spin) = n_states
       endif
      enddo

      temp_matr(:) = 0.d0

      if(myid.eq.0) then
        write(use_unit,'(2X, A,A,4I5)') &
                   "HOMO and first non-fully-occupied", &
                   " orbitals:", n_homo(:), n_first(:)
        write(use_unit,*) 
      endif
!      if (chemical_potential .lt. KS_eigenvalue (n_homo) .or.
!     +    chemical_potential .gt. KS_eigenvalue (n_lumo) ) then
!       write(use_unit,'(2X,A)') "Warning! chemical potential is not sitting"
!       write(use_unit,'(2X,A)') "in between HOMO and LUMO!"
!      endif

      n_states_loc = (n_high_state-1)/n_tasks + 1

      allocate(polar_freq(n_basbas, n_loc_prodbas), stat=info)
      call check_allocation(info, 'polar_freq', func)
      allocate(screened_coulomb( n_basbas, n_loc_prodbas ), stat=info)
      call check_allocation(info, 'screened_coulomb', func)

      allocate(aux_ovlp3KS(n_basbas, n_states ), stat=info)
      call check_allocation(info, 'aux_ovlp3KS', func)
      allocate(aux_coulomb_matr &
             (n_states_loc,n_states,n_full_freq,n_spin), stat=info)
      call check_allocation(info, 'aux_coulomb_matr', func)

      allocate(screened_exchange(n_high_state), stat=info)
      call check_allocation(info, 'screened_exchange', func)
      allocate(coulomb_hole(n_high_state,n_spin), stat=info)
      call check_allocation(info, 'coulomb_hole', func)

      allocate(ovlp_3KS_redist(n_basbas, n_states, n_states_loc, n_spin), stat=info)
      call check_allocation(info, 'ovlp_3KS_redist', func)
      allocate(screened_coulomb_f( n_basbas, n_basbas ), stat=info)
      call check_allocation(info, 'screened_coulomb_f', func)

      allocate(my_KS_eigenvalue(n_states,n_spin)) 
      if(present(KS_eigenvalue_scf)) then
        my_KS_eigenvalue(:,:) = KS_eigenvalue_scf(1:n_states,1:n_spin)
      else
        my_KS_eigenvalue(:,:) = KS_eigenvalue(1:n_states,1:n_spin)
      endif 

      ! Redistribute ovlp_3KS, the parallelization is over the 3rd parameter in ovlp_3KS_redist

      ovlp_3KS_redist = 0
      do i_spin = 1, n_spin
        do i_state = 1,n_high_state
          aux_ovlp3KS = 0.
          do j_state = 1, n_states
            do j_basbas = 1, n_loc_prodbas
              i_index = map_prodbas(j_basbas,myid+1)
              aux_ovlp3KS(i_index,j_state) = ovlp_3KS(j_basbas,j_state,i_state,i_spin)
            enddo
          enddo
          call sync_matrix(aux_ovlp3KS, n_basbas, n_states)
          if(MOD(i_state-1,n_tasks) == myid) then
            i_state_loc = (i_state-1)/n_tasks + 1
            ovlp_3KS_redist(:,:,i_state_loc,i_spin) = aux_ovlp3KS(:,:)
          endif
        enddo
      enddo

      self_energy_freq=0.d0
       
      call cpu_time (rtime) 
      tot_time_polar = rtime - tot_time_polar

      tot_time_sigma = 0.d0
      time_product = 0.d0

      aux_coulomb_matr =0.d0

      do i_freq = 1, n_full_freq, 1

          if(myid.eq.0) then
            write(use_unit,'(A,I4,f20.6)') " | i_freq ", i_freq, &
                   omega_full(i_freq)
          endif

          call cpu_time (time_polar)

!    evaluate the polarisability at frequency point i_freq 
!          call  evaluate_polarisability_freq &
!             ( 1, n_homo, n_first, n_high_state, occ_numbers, &
!               omega_full(i_freq), &
!               KS_eigenvalue, ovlp_3KS, polar_freq &
!             )
          call  evaluate_polarisability_freq_1 &
             ( 1, n_homo, n_first, n_high_state, occ_numbers, &
               omega_full(i_freq), &
               my_KS_eigenvalue, ovlp_3KS_redist, screened_coulomb_f &
             )
          do j_basbas = 1, n_loc_prodbas, 1
            i_index = map_prodbas(j_basbas,myid+1)
            do i_basbas = 1, n_basbas, 1
               polar_freq(i_basbas, j_basbas) = screened_coulomb_f(i_basbas, i_index)
             enddo
          enddo
          call cpu_time (rtime)
          time_polar = rtime - time_polar
          tot_time_polar = tot_time_polar + time_polar

          call cpu_time(time_start)
!    now calculate the screened Coulomb interaction W
          call screened_coulomb_interaction &
             ( polar_freq, screened_coulomb &
             )
          call cpu_time(time_end)
          time_product = time_product + time_end-time_start

!    subtract the bare Coulomb potential from W, Wc = W-V 
          do j_basbas = 1, n_loc_prodbas, 1
            
            i_index = map_prodbas(j_basbas,myid+1)
            do i_basbas = 1, n_basbas, 1

              if(i_index.eq.i_basbas) then 
                 screened_coulomb(i_basbas, j_basbas) = &
                   screened_coulomb(i_basbas, j_basbas) - &
                    1.d0
              endif
             enddo
          enddo

!    get the full matrix on all PEs
          screened_coulomb_f = 0
          do j_basbas = 1, n_loc_prodbas, 1
            i_index = map_prodbas(j_basbas,myid+1)
            do i_basbas = 1, n_basbas, 1
               screened_coulomb_f(i_basbas, i_index) = screened_coulomb(i_basbas, j_basbas)
             enddo
          enddo
          call sync_matrix(screened_coulomb_f,n_basbas,n_basbas)
       
!    calculate the contribution of the Coulomb hole term:
!    Sig_COH =<i|W(r,r,0)|i>
          if(i_freq.eq.1 .and. flag_coh) then 
             call integrate_coulombhole &
             ( n_high_state, &
               partition_tab, basis_l_max, &
               screened_coulomb, KS_eigenvector, &
               coulomb_hole)
          endif

!     loop over KS states which are of interest.
          do i_spin = 1, n_spin
            do i_state = n_low_state,n_high_state, 1

               if(MOD(i_state-1,n_tasks) /= myid) cycle
               i_state_loc = (i_state-1)/n_tasks + 1
 
               aux_ovlp3KS =0.d0

               call dgemm('N', 'N', n_basbas, n_states, &
                    n_basbas, 1.0d0, &
                    screened_coulomb_f, & 
                    n_basbas, &
                    ovlp_3KS_redist(:,:,i_state_loc,i_spin), &
                    n_basbas, 0.d0, &
                    aux_ovlp3KS, n_basbas &
                  )

               do j_state = 1, n_states
                  aux_coulomb_matr(i_state_loc,j_state,i_freq,i_spin) = &
                     dot_product(aux_ovlp3KS(:,j_state),ovlp_3KS_redist(:,j_state,i_state_loc,i_spin))
               enddo
!      end of i_state
            enddo
!        end of i_spin
          enddo 
          call cpu_time (time_self_energy)
          time_self_energy = time_self_energy - rtime
          tot_time_sigma = tot_time_sigma + time_self_energy

!      end of mpi task distribution
!        endif

!      end of i_freq 
       enddo


       call cpu_time (time_self_energy)
!      omega_1=omega(n_freq)
!      omega_2=omega(n_freq-1)
      do i_spin = 1, n_spin
       do i_state = n_low_state, n_high_state, 1 

         if(MOD(i_state-1,n_tasks) /= myid) cycle
         i_state_loc = (i_state-1)/n_tasks + 1

         do j_state = 1, n_states, 1
          
!         s1= aux_coulomb_matr(i_state,j_state,n_freq)
!         s2= aux_coulomb_matr(i_state,j_state,n_freq-1)
!         alpha=(omega_2*omega_2 - omega_1*omega_1)*s1*s2/(s1-s2)
!         beta=(omega_2*omega_2*s2 - omega_1*omega_1*s1)/(s1-s2)

!          if(beta.le.0) then
!           alpha=0
!           beta=1.d0
!          endif

          do i_freq_1= 1, n_freq 
              zeta = (0.d0,1.d0)*omega(i_freq_1)-  &
                    KS_eigenvalue(j_state,i_spin)+chemical_potential_spin(i_spin)

              arg = 2.d0/(1.d0-(0.d0,1.d0)*zeta/omega(n_freq-1))-1.d0

!              prefactor = alpha*zeta/(zeta*zeta-beta) 
           
              do i_freq=1, n_full_freq

               if(i_freq_1.ne.i_freq) then
                self_energy_freq(i_state,i_freq_1,i_spin) = &
                self_energy_freq (i_state,i_freq_1,i_spin) - &
               ( aux_coulomb_matr(i_state_loc,j_state,i_freq,i_spin) - &
                 aux_coulomb_matr(i_state_loc,j_state, i_freq_1,i_spin) )* &
                 zeta/(zeta*zeta + &
                      omega_full(i_freq)*omega_full(i_freq)) &
                *womega_full(i_freq)/pi 
               endif

!             if(i_state .eq. 10 .and. i_freq_1 .eq. 10) then
!                write(use_unit,'(3I4,5f12.5)') i_spin, j_state, i_freq, aux_coulomb_matr(i_state_loc,j_state,i_freq,i_spin), &
!                zeta, self_energy_freq(i_state,i_freq_1,i_spin)
!             endif

!   end of i_freq    
             enddo

            if(j_state .le. n_homo(i_spin) ) then
              self_energy_freq(i_state,i_freq_1,i_spin) = &
              self_energy_freq(i_state,i_freq_1,i_spin) - &
              aux_coulomb_matr(i_state_loc, j_state, i_freq_1,i_spin)/2.d0

!               if(use_cohsex) then
!                self_energy_freq(i_state, i_freq_1,i_spin) =
!     +            self_energy_freq (i_state, i_freq_1,i_spin) +
!     +            aux_coulomb_matr(i_state, j_state, 1, i_spin)/2.d0
!               elseif(use_screx) then
!                self_energy_freq(i_state, i_freq_1,i_spin) =
!     +            self_energy_freq (i_state, i_freq_1,i_spin) +
!     +            aux_coulomb_matr(i_state, j_state, 1, i_spin)
!                if(i_freq_1.eq.1) then
!                  temp_matr(i_state) = temp_matr(i_state) +
!     +              aux_coulomb_matr(i_state,j_state,1,i_spin) /2.d0
!                endif
!               endif

!                self_energy_freq(i_state, i_freq_1) =
!     +           self_energy_freq (i_state, i_freq_1) -
!     +           aux_coulomb_matr(i_state, j_state, i_freq_1)/2.d0*
!     +            (occ_numbers(j_state) - 1.d0)
            else
              self_energy_freq(i_state,i_freq_1,i_spin) = &
              self_energy_freq (i_state,i_freq_1,i_spin) + &
              aux_coulomb_matr(i_state_loc, j_state, i_freq_1,i_spin)/2.d0

!               if(use_cohsex) then
!                 self_energy_freq(i_state, i_freq_1,i_spin) =
!     +             self_energy_freq (i_state, i_freq_1,i_spin) -
!     +             aux_coulomb_matr(i_state, j_state, 1, i_spin)/2.d0
!!                if(i_freq_1.eq.1) then
!!                  temp_matr(i_state) = temp_matr(i_state) -
!!     +              aux_coulomb_matr(i_state,j_state,1,i_spin) /2.d0
!!                endif
!               endif
!
            endif
             
!    end of i_freq_1 loop
          enddo 

!    end of loop j_state
        enddo
         
        if(flag_coh) then
           self_energy_freq(i_state,:,i_spin) = &
            self_energy_freq (i_state,:,i_spin) - &
             coulomb_hole(i_state,i_spin)/2.d0
        endif
!        write(use_unit,*) i_state,  temp_matr(i_state)*27.21138

!     +           screened_exchange(i_state)/pi
!        write(use_unit,*) i_state, coulomb_hole(i_state) 
! end of i_state loop
       enddo
! end of i_spin loop
      enddo

      call sync_vector_complex(self_energy_freq, (n_high_state-n_low_state+1)*n_freq*n_spin)
   
 
      !evaluate Coulomb correlation energy here
      if(get_gw_corr)then
        gw_corr = 0.d0
        do i_spin = 1, n_spin
          do i_state = n_low_state , n_high_state
            do i_freq = 1, n_freq 
              gw_corr = gw_corr + self_energy_freq (i_state,i_freq, i_spin) * &
              1.d0/((0.d0,1.d0)*omega(i_freq)-(KS_eigenvalue(i_state,i_spin)&
              -chemical_potential_spin(i_spin)))*womega(i_freq) *&
              (2.d0/n_spin)/2.d0/pi
            enddo
          enddo 
        enddo 
        if(myid.eq.0) then
           write(use_unit,*) "GW correlation:  ", gw_corr*hartree
        endif
      endif
!    printing out  
!       if (myid.eq.0 .and. out_self_energy) then
!        do i_spin = 1, n_spin
!        do j_state =  n_low_state, n_high_state, 1

!         if(n_spin.eq.2) then
!          write(filename,'(A,I0,A,I0,A)')"self_energy/Sigma.omega.n_", &
!                       j_state,".s_",i_spin,".dat"
!          else
!           write(filename,'(A,I0,A)')"self_energy/Sigma.omega.n_", &
!                       j_state,".dat"
!         endif

!          open(102, file=filename)

!          do i_freq_1 = 1, n_freq, 1

!             write (102, '(2X, f16.8, 7X, 2F18.10)') &  
!                omega(i_freq_1), &
!                self_energy_freq(j_state,i_freq_1,i_spin) 
!           enddo
!          
!           close(102)
!       enddo
!      enddo

!      open (103, file='self_energy_omega.dat')
!      write(103,*) n_freq, n_high_state - n_low_state +1
!      write(103,*) (omega(i_freq_1),i_freq_1=1,n_freq)
!      do i_spin = 1, n_spin 
!      do i_state=n_low_state,n_high_state
!         write(103,*)i_state
!         write(103,*) &
!            (self_energy_freq(i_state,i_freq_1,i_spin), &
!                   i_freq_1=1,n_freq)
!      enddo
!      enddo 
!      close(103)
! end of if out_self_energy
!     endif

!      if(myid.eq.0)then
!      open (103, file='self_energy_old.dat')
!      do i_freq_1 = 1, n_freq, 1
!EMPTY STATES
!        write(103,*) omega(i_freq_1), &
!                real(self_energy_freq(n_homo(1)+1,i_freq_1,1))
!        write(103,*) omega(i_freq_1), &
!                real(self_energy_freq(n_homo(1)+2,i_freq_1,1))
!        write(103,*) omega(i_freq_1), &
!                real(self_energy_freq(n_homo(1)+3,i_freq_1,1))
!OCCUPIED STATES
!        write(103,*) omega(i_freq_1), &
!             real(self_energy_freq(n_homo,n_homo,i_freq_1,1)), &
!            aimag(self_energy_freq(n_homo,n_homo,i_freq_1,1))
!        write(103,*) omega(i_freq_1), &
!             real(self_energy_freq(n_homo-1,n_homo-1,i_freq_1,1)), &
!            aimag(self_energy_freq(n_homo-1,n_homo-1,i_freq_1,1))
!        write(103,*) omega(i_freq_1), &
!             real(self_energy_freq(n_homo-2,n_homo-2,i_freq_1,1)), &
!            aimag(self_energy_freq(n_homo-2,n_homo-2,i_freq_1,1))

!        write(103,*) omega(i_freq_1), &
!             real(self_energy_freq(1,1,i_freq_1,1)), aimag(self_energy_freq(1,1,i_freq_1,1))
!        write(103,*) omega(i_freq_1), &
!             real(self_energy_freq(n_homo,n_homo,i_freq_1,1)), aimag(self_energy_freq(n_homo,n_homo,i_freq_1,1))
!        write(103,*) omega(i_freq_1), &
!             real(self_energy_freq(n_states,n_states,i_freq_1,1)), aimag(self_energy_freq(n_states,n_states,i_freq_1,1))

!      enddo
!      close(103)
!      endif



      call cpu_time (rtime)
      time_self_energy = rtime - time_self_energy 
      tot_time_sigma = tot_time_sigma + time_self_energy


      if (allocated (aux_ovlp3KS)) then
        deallocate (aux_ovlp3KS)
      endif
      if (allocated (aux_coulomb_matr)) then
        deallocate (aux_coulomb_matr)
      endif
      if (allocated (polar_freq)) then
        deallocate (polar_freq)
      endif
      if (allocated (screened_coulomb)) then
        deallocate (screened_coulomb)
      endif
      if (allocated (screened_exchange)) then
        deallocate (screened_exchange)
      endif
      if (allocated (coulomb_hole)) then
        deallocate (coulomb_hole)
      endif
      if (allocated (ovlp_3KS_redist)) then
        deallocate (ovlp_3KS_redist)
      endif
      if (allocated (screened_coulomb_f)) then
        deallocate (screened_coulomb_f)
      endif
      if (allocated (my_KS_eigenvalue)) then
        deallocate (my_KS_eigenvalue)
      endif
      return 
  end subroutine evaluate_self_energy_freq_1

! **************************************************************************************************
  subroutine evaluate_self_energy_freq_2 &
            (n_low_state, n_high_state, &
             n_homo, occ_numbers,  &
             n_freq,n_full_freq, &
             omega, omega_full, womega_full, &
             chemical_potential_spin, &
             KS_eigenvalue, KS_eigenvector, &
             partition_tab,  basis_l_max, &
             ovlp_3KS, &
             self_energy_freq, &
             out_self_energy, flag_coh, &
             tot_time_sigma, tot_time_polar, &
             Wmn_matrix_freq, &
             KS_eigenvalue_scf &
            )

!  PURPOSE
!  Subroutine evaluate_self_energy_freq evaluates the correlated part of 
!  the GW self-nergy on the imaginary frequency axis. 
!  
!  Sigma^_ii(iw) = i \sum_lmn O_ilm O_iln \int G_ll(iw+iw') W _mn(iw') dw'


!  ARGUMENTS

      integer :: n_freq
      integer :: n_full_freq
      integer :: n_low_state
      integer :: n_high_state
      integer :: basis_l_max (n_species)
      integer :: n_homo(n_spin)

      real*8  :: occ_numbers(n_states,n_spin)
      real*8  :: omega(n_freq)
      real*8  :: omega_full(n_full_freq)
      real*8  :: womega_full(n_full_freq)
      real*8  :: chemical_potential_spin(n_spin)   
      real*8  :: KS_eigenvalue(n_states,n_spin)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin)
      real*8  :: ovlp_3KS(n_basbas, ndim1_o3KS, ndim2_o3KS, n_spin)
      real*8  :: partition_tab(n_max_angular, n_max_radial, n_atoms)

      logical out_self_energy
      logical flag_coh

      real*8  tot_time_polar
      real*8  tot_time_sigma

      complex*16 ::  self_energy_freq &
            ( n_low_state:n_high_state, n_freq, n_spin )
      real*8, intent(inout), optional     :: Wmn_matrix_freq &
                                            ( ndim2_o3KS,ndim1_o3KS,n_full_freq,n_spin)
      real*8, intent(in), optional        :: KS_eigenvalue_scf(n_states,n_spin)

! INPUTS
! o  n_freq -- integer number, the number of frequency points for the GW self-energy
! o  n_full_freq -- integer number,
!            the number of frequency points for the screened Coulomb interaction W
! o  n_low_state  -- integer number,
!            the lowest KS/HF eigenstate for self-energy correction
! o  n_high_state -- integer number,
!            the highest KS/HF eigenstate for self-energy correction
! o  basis_l_max -- integer array of length n_species, 
!            the maximal angular momentum for the basis functions for each species
! o  n_homo -- integer array of length of n_spin,
!            the HOMO level for each spin channel
! o  occ_numbers -- real 2-dimentianal array of length (n_states, n_spin)
!            the occupation number of the electrons for each eigenstate and each spin
! o  omega(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the self-energy
! o  omega_full(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the screened Coulomb interaction
! o  womega_full(n_freq) -- real array
!            the weigth of the Gauss-Legendre frequency grid for the screened Coulomb interaction
! o  chemical_potential -- real number, the chemical potential of the system
! o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF 
!            eigenvalue
! o  KS_eigenvector -- real array,
!            the eigenvector of the single-particle calculation
! o  ovlp_3KS -- real array
!            this is the transformed 3-cener overlap integration. Now two orbitals of
!            them are KS ones, and one is the auxiliary basis.
!            Note: for parallel calculations, the auxiliary basis are distribuated
!            among the different processors. 
! o  out_self_energy -- logic, if ture, output the self-energ on the imagniary axis 
! o  flag_coh -- logic, if ture, use the Coulomb-hole accelerator 
!             (not functioning for the time being)

! OUTPUTS
! o  self_energy_freq  -- complex array, 
!           the calculated GW self-energy on the imaginary axis
! o  tot_time_polar    -- total time for calculating the polarisability
! o  tot_time_sigma    -- total time for calculating the self-energy (except for the
!                          polarisability)
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

!     auxiliary matrices for Level 3 Blas matrix multiplications
!     n_first : the first orbital which is NOT fully occupied

      integer :: n_first(n_spin)
      complex*16  :: zeta
      complex*16  ::  arg

      real*8, dimension(:,:), allocatable :: polar_freq
      real*8, dimension(:,:), allocatable :: screened_coulomb

      real*8, dimension(:,:,:,:), allocatable :: aux_coulomb_matr
      real*8, dimension(:), allocatable :: screened_exchange
      real*8, dimension(:,:), allocatable :: coulomb_hole
      real*8  temp_matr(n_states)
      real*8, dimension(:,:), allocatable :: my_KS_eigenvalue

!     timing
      real(kind=8), dimension(2) :: tmp_time_pol, tmp_time_Wmn
      real*8  :: time_self_energy 
      real*8  :: time_start, time_end, time_product
      real*8  :: rtime

!     parameters of the fitting tails
!       real*8  s1, s2, omega_1, omega_2
!       real*8  alpha, beta, a, b

      character*50  filename 
      
!     counters


      integer :: i_state
      integer :: j_state

      integer :: i_freq
      integer :: i_freq_1

      integer :: i_spin
      integer :: i_state_loc, j_state_loc

      integer :: info
      character(*), parameter :: func = 'evaluate_self_energy_freq_2'

!     begin work

      if(myid.eq.0) then  
        write(use_unit,*) 
        write(use_unit,*)"-------------------------------------------------" 
        write(use_unit,'(2X,A)') "Start to calculate the self energy ... "
      endif

      ! Block size for matrix multiplies, must be a multiple of nb_aux_2d

!     determine the highest occupied orbital level
!     such complication occurs when some of the orbitals are 
!     either not fully occupied or not fully empty
      n_first(:) = 1
      do i_spin = 1, n_spin
       do i_state = 1, n_states
        if (abs(occ_numbers(i_state,i_spin)-dble(2/n_spin)) &
                        .lt.1.d-8) then
         n_first(i_spin)= i_state + 1
        endif 
       enddo 
       if(n_first(i_spin) .gt. n_states) then
         n_first(i_spin) = n_states
       endif
      enddo

      temp_matr(:) = 0.d0

      if(myid.eq.0) then
        write(use_unit,'(2X, A,A,4I5)') &
                   "HOMO and first non-fully-occupied", &
                   " orbitals:", n_homo(:), n_first(:)
        write(use_unit,*) 
      endif
!      if (chemical_potential .lt. KS_eigenvalue (n_homo) .or.
!     +    chemical_potential .gt. KS_eigenvalue (n_lumo) ) then
!       write(use_unit,'(2X,A)') "Warning! chemical potential is not sitting"
!       write(use_unit,'(2X,A)') "in between HOMO and LUMO!"
!      endif


      allocate(polar_freq(max(1,max_row_2d), max(1,max_col_2d)), stat=info)
      call check_allocation(info, 'polar_freq', func)
      allocate(screened_coulomb(max(1,max_row_2d), max(1,max_col_2d)), stat=info)
      call check_allocation(info, 'screened_coulomb', func)

      allocate(aux_coulomb_matr(ndim2_o3KS,ndim1_o3KS,n_full_freq,n_spin), stat=info)
      call check_allocation(info, 'aux_coulomb_matr', func)

      allocate(screened_exchange(n_high_state), stat=info)
      call check_allocation(info, 'screened_exchange', func)
      allocate(coulomb_hole(n_high_state,n_spin), stat=info)
      call check_allocation(info, 'coulomb_hole', func)

      !
      allocate(my_KS_eigenvalue(n_states,n_spin))  
      if(present(KS_eigenvalue_scf)) then
        my_KS_eigenvalue(:,:) = KS_eigenvalue_scf(1:n_states,1:n_spin)
      else
        my_KS_eigenvalue(:,:) = KS_eigenvalue(1:n_states,1:n_spin)
      endif 


      self_energy_freq=0.d0
       
      tot_time_sigma = 0.d0
      time_product = 0.d0

      aux_coulomb_matr = 0.d0

      tmp_time_pol(1) = time_polar_imagfreq
      tmp_time_pol(2) = clock_time_polar_imagfreq
      tmp_time_Wmn(:) = 0.d0 
      call get_timestamps(tmp_time_Wmn(1),tmp_time_Wmn(2))
      do i_freq = 1, n_full_freq, 1

          if(myid.eq.0) then
            write(use_unit,'(A,I4,f20.6)') " | i_freq ", i_freq, &
                   omega_full(i_freq)
          endif

!    evaluate the polarisability at frequency point i_freq 
          call  evaluate_polarisability_freq_2 &
             ( 1, n_homo, n_first, occ_numbers, &
               omega_full(i_freq), &
               my_KS_eigenvalue, ovlp_3KS, polar_freq, time_polar=tmp_time_pol &  
             )

          call cpu_time(time_start)
!    now calculate the screened Coulomb interaction W
          call PDLASET( 'Full', n_basbas, n_basbas, 0.d0, 1.d0, screened_coulomb, 1, 1, aux_sc_desc_2d )
          screened_coulomb(:,:) = screened_coulomb(:,:) - polar_freq(:,:)*2.d0/dble(n_spin)

          !if(use_scalapack) then
             call power_auxmat_scalapack_2d(screened_coulomb, -1.d0, '')
          !else
          !  print *,'power_auxmat_lapack_2d not yet implemented!'
          !  call aims_stop
          !  call power_auxmat_lapack_2d(screened_coulomb, -1.d0, '')
          !endif

          call cpu_time(time_end)
          time_product = time_product + time_end-time_start

!    subtract the bare Coulomb potential from W, Wc = W-V 

          ! polar_freq isn't needed any more, set it to unity matrix
          call PDLASET( 'Full', n_basbas, n_basbas, 0.d0, 1.d0, polar_freq, 1, 1, aux_sc_desc_2d )
          screened_coulomb(:,:) = screened_coulomb(:,:) - polar_freq(:,:)
       
!    calculate the contribution of the Coulomb hole term:
!    Sig_COH =<i|W(r,r,0)|i>
          if(i_freq.eq.1 .and. flag_coh) then 
             if(n_tasks>1) then
                if(myid==0) print *,'ERROR: integrate_coulombhole works only for 1 proc!'
                call aims_stop
             endif
             call integrate_coulombhole &
             ( n_high_state, &
               partition_tab, basis_l_max, &
               screened_coulomb, KS_eigenvector, &
               coulomb_hole)
          endif

          !*** calculate matrix element Wmn from (1-Pi(i\omega))^-1
          call calc_Wmn_from_WPQ_imag_freq(screened_coulomb, ovlp_3KS, aux_coulomb_matr, i_freq, n_low_state,&
                                  n_high_state)
          call cpu_time (time_self_energy)
          time_self_energy = time_self_energy - time_start
          tot_time_sigma = tot_time_sigma + time_self_energy

!      end of i_freq 
       enddo
       call get_times(tmp_time_Wmn(1), tmp_time_Wmn(2), time_Wmn_imagfreq, clock_time_Wmn_imagfreq) 

      call cpu_time (time_self_energy)
      do i_spin = 1, n_spin
        do i_state = n_low_state, n_high_state, 1 

          if(own_dim2_o3KS(i_state) /= myp2_o3KS) cycle
          i_state_loc = loc_dim2_o3KS(i_state)

          do j_state = 1, n_states, 1

            if(own_dim1_o3KS(j_state) /= myp1_o3KS) cycle
            j_state_loc = loc_dim1_o3KS(j_state)

            do i_freq_1= 1, n_freq 
              zeta = (0.d0,1.d0)*omega(i_freq_1)-  &
                    KS_eigenvalue(j_state,i_spin)+chemical_potential_spin(i_spin)

              arg = 2.d0/(1.d0-(0.d0,1.d0)*zeta/omega(n_freq-1))-1.d0

!             prefactor = alpha*zeta/(zeta*zeta-beta) 
           
              do i_freq=1, n_full_freq

                if(i_freq_1.ne.i_freq) then
                  self_energy_freq(i_state,i_freq_1,i_spin) = &
                  self_energy_freq (i_state,i_freq_1,i_spin) - &
                 ( aux_coulomb_matr(i_state_loc,j_state_loc,i_freq,i_spin) - &
                   aux_coulomb_matr(i_state_loc,j_state_loc, i_freq_1,i_spin) )* &
                   zeta/(zeta*zeta + &
                        omega_full(i_freq)*omega_full(i_freq)) &
                  *womega_full(i_freq)/pi 
                endif

!   end of i_freq    
              enddo
            
              if( j_state .le. n_homo(i_spin) ) then
                self_energy_freq(i_state,i_freq_1,i_spin) = &
                self_energy_freq(i_state,i_freq_1,i_spin) - &
                aux_coulomb_matr(i_state_loc, j_state_loc, i_freq_1,i_spin)/2.d0

              else
                self_energy_freq(i_state,i_freq_1,i_spin) = &
                self_energy_freq (i_state,i_freq_1,i_spin) + &
                aux_coulomb_matr(i_state_loc, j_state_loc, i_freq_1,i_spin)/2.d0

              endif
             
!    end of i_freq_1 loop
            enddo 

!    end of loop j_state
          enddo
         
          if(flag_coh) then
            self_energy_freq(i_state,:,i_spin) = &
             self_energy_freq (i_state,:,i_spin) - &
              coulomb_hole(i_state,i_spin)/2.d0
          endif
! end of i_state loop
        enddo
! end of i_spin loop
      enddo

      call sync_vector_complex(self_energy_freq, (n_high_state-n_low_state+1)*n_freq*n_spin)

      if(present(Wmn_matrix_freq)) then
        call check_allocation(info, 'Wmn_matrix_freq', func)
        Wmn_matrix_freq(:,:,:,:) = aux_coulomb_matr 
      endif

!    printing out  
       if (myid.eq.0 .and. out_self_energy) then
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

           do i_freq_1 = 1, n_freq, 1

              write (102, '(2X, f16.8, 7X, 2F18.10)') &  
                 omega(i_freq_1), &
                 self_energy_freq(j_state,i_freq_1,i_spin) 
            enddo
           
            close(102)
        enddo
       enddo

       open (103, file='self_energy_omega.dat')
       write(103,*) n_freq, n_high_state - n_low_state +1
       write(103,*) (omega(i_freq_1),i_freq_1=1,n_freq)
       do i_spin = 1, n_spin 
       do i_state=n_low_state,n_high_state
          write(103,*)i_state
          write(103,*) &
             (self_energy_freq(i_state,i_freq_1,i_spin), &
                    i_freq_1=1,n_freq)
       enddo
       enddo 
       close(103)
! end of if out_self_energy
      endif

      call cpu_time (rtime)
      time_self_energy = rtime - time_self_energy 
      tot_time_sigma = tot_time_sigma + time_self_energy


      if (allocated (aux_coulomb_matr)) then
        deallocate (aux_coulomb_matr)
      endif
      if (allocated (polar_freq)) then
        deallocate (polar_freq)
      endif
      if (allocated (screened_coulomb)) then
        deallocate (screened_coulomb)
      endif
      if (allocated (screened_exchange)) then
        deallocate (screened_exchange)
      endif
      if (allocated (coulomb_hole)) then
        deallocate (coulomb_hole)
      endif
      if (allocated (my_KS_eigenvalue)) then
        deallocate (my_KS_eigenvalue)
      endif

      
      call sync_timing(tmp_time_pol(1)) 
      time_polar_imagfreq = tmp_time_pol(1) 
      clock_time_polar_imagfreq = tmp_time_pol(2)
      tot_time_polar = time_polar_imagfreq

      return 
  end subroutine evaluate_self_energy_freq_2

! **************************************************************************************************
!> brief evalutes self energy for the contour deformation (CD) correction
!  o gw_cd  --  contour deformation environment
!  o aux_coulomb_matr_freq -- Wmn(i\omega) matrix elements
!  o i_state -- ith KS state
!  o i_spin -- ith spin
!  o n_homo -- the HOMO level, i.e., the number of occupied state 
!  o occ_number -- occupation number for each state and each spin
!  o n_full_freq -- the number of frequency points for numerical integration
!  o omega_full(n_freq) -- the Gauss-Legendre frequency grid for numerical integration
!  o womega_full(n_freq) -- the weigth of the Gauss-Legendre frequency grid for the numerical
!                           interaction
!  o chemical_potential -- the chemical potential of the system
!  o real_omega -- real frequency, i.e. we have Sigma_istate(real_omega)
!  o KS_eigenvalue -- KS/HF eigenvalues of the single-particle calculation
!  o ovlp_3KS -- transformed 3-center integrals
! **************************************************************************************************
  subroutine evaluate_self_energy_cd(gw_cd, aux_coulomb_matr_freq, i_state, i_spin, n_homo, n_first,&
                                     occ_numbers, n_full_freq, omega_full, womega_full, &
                                     chemical_potential_spin, real_omega, KS_eigenvalue, &
                                     KS_eigenvalue_last, ovlp_3KS, do_scalapack) 
 

      type(cd_environment_type)                                :: gw_cd
      complex(kind=8), dimension(ndim2_o3KS,ndim1_o3KS,&
        n_full_freq, n_spin), intent(in)                       :: aux_coulomb_matr_freq
      integer, intent(in)                                      :: i_state
      integer, intent(in)                                      :: i_spin
      integer, dimension(n_spin), intent(in)                   :: n_homo
      integer, dimension(n_spin), intent(in)                   :: n_first
      real(kind=8), dimension(n_states,n_spin), intent(in)     :: occ_numbers
      integer, intent(in)                                      :: n_full_freq
      real(kind=8), dimension(n_full_freq), intent(in)         :: omega_full
      real(kind=8), dimension(n_full_freq), intent(in)         :: womega_full
      real(kind=8), dimension(n_spin), intent(in)              :: chemical_potential_spin   
      real(kind=8), intent(in)                                 :: real_omega
      real(kind=8), dimension(n_states,n_spin), intent(in)     :: KS_eigenvalue
      real(kind=8), dimension(n_states,n_spin), intent(in)     :: KS_eigenvalue_last
      real(kind=8), &
        dimension(n_basbas, ndim1_o3KS, ndim2_o3KS, n_spin),&
        intent(in)                                             :: ovlp_3KS
      logical, intent(in), optional                            :: do_scalapack

      character(*), parameter :: func = 'evaluate_self_energy_cd'
 
      logical                                                  :: my_do_scalapack
      integer                                                  :: j_state, i_freq, i_freq_real,&
                                                                  i_state_loc, j_state_loc,&
                                                                  info, index_contour_def 
      real(kind=8), dimension(:,:), allocatable                :: tmp_KS
      real(kind=8), dimension(2)                               :: tmp_time
      complex(kind=8), dimension(:), allocatable               :: aux_coulomb_residue_cmplx

      tmp_time(:) = 0.d0 
      call get_timestamps(tmp_time(1), tmp_time(2))

      my_do_scalapack = .true.
      if(present(do_scalapack)) then
        my_do_scalapack = do_scalapack
      endif

      allocate(aux_coulomb_residue_cmplx(gw_cd%num_residues),stat=info)
      call check_allocation(info, 'aux_coulomb_residue_cmplx', func)
      aux_coulomb_residue_cmplx = 0.d0

      !*** calculate matrix element W_mn(e_m-e_n) for all residues e_m-e_n
      if(use_ev_scgw) then
         call calc_Wmn_for_real_frequencies(gw_cd, aux_coulomb_residue_cmplx,occ_numbers,&
                                            KS_eigenvalue_last,ovlp_3KS,i_state,i_spin, &
                                            n_homo,n_first,do_scalapack)
      else
         call calc_Wmn_for_real_frequencies(gw_cd,aux_coulomb_residue_cmplx,occ_numbers,&
                                            KS_eigenvalue,ovlp_3KS,i_state,i_spin,n_homo,& 
                                            n_first,do_scalapack)
      endif 

      
      !*** calculate now the self-energy from W_mn(i\omega) and W_mn(e_m-e_n)
      do i_freq = 1, n_full_freq + gw_cd%num_residues
         if(own_dim2_o3KS(i_state) == myp2_o3KS) then
           i_state_loc = loc_dim2_o3KS(i_state)

           do j_state = 1, n_states, 1

              if(own_dim1_o3KS(j_state) /= myp1_o3KS) cycle
              j_state_loc = loc_dim1_o3KS(j_state)
              if(i_freq <= n_full_freq) then
                call contour_def_integrate(gw_cd, aux_coulomb_matr_freq(i_state_loc,j_state_loc,i_freq,i_spin),&
                                           real_omega,KS_eigenvalue_last,i_spin,i_state,j_state,&
                                           i_freq,omega_full,womega_full,chemical_potential_spin(i_spin))
              else
                i_freq_real = i_freq - n_full_freq
                call contour_def_add_residues(gw_cd,gw_cd%self_energy%complx,&
                                              aux_coulomb_residue_cmplx(i_freq_real),&
                                              real_omega,chemical_potential_spin(i_spin),i_spin,&
                                              i_state,j_state,i_freq_real)
                
              endif

           enddo !j_state
         endif
      enddo ! i_freq

      index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
      if(my_do_scalapack) then 
        call sync_complex_number(gw_cd%self_energy%complx(index_contour_def,i_spin))
      endif
      gw_cd%self_energy%re(index_contour_def,i_spin) = &
                          real(gw_cd%self_energy%complx(index_contour_def,i_spin), kind=8)
      gw_cd%self_energy%im(index_contour_def,i_spin) = &
                          aimag(gw_cd%self_energy%complx(index_contour_def,i_spin))

     if (allocated (aux_coulomb_residue_cmplx)) then
       deallocate (aux_coulomb_residue_cmplx)
     endif

     call get_times(tmp_time(1), tmp_time(2), time_self_energy_cd, clock_time_self_energy_cd,&
                    unsynced=.TRUE.)

  end subroutine evaluate_self_energy_cd

! **************************************************************************************************
!> brief construct screened coulomb interaction Wc in the KS basis for all residues, 
!>       i.e. W_mn(e_m-e_n)
!  o gw_cd  --  contour deformation environment
!  o aux_coulomb_residue_cmplx -- W_mn matrix (also with imaginary part if needed) 
!  o occ_number -- occupation number for each state and each spin, dim: n_states,n_spin
!  o KS_eigenvalue -- KS/HF eigenvalues of the single-particle calculation, dim: n_states,n_spin
!  o ovlp_3KS -- transformed 3-center integrals, dim: (n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
!  o n_homo -- the HOMO level for each spin channel
!  o n_first -- the first orbital which is NOT fully occupied
!  o i_state -- ith KS state
!  o i_spin -- ith spin
!  o do_scalapack -- if we do scalapack or not 
! **************************************************************************************************
  subroutine calc_Wmn_for_real_frequencies(gw_cd, aux_coulomb_residue_cmplx, occ_numbers,&
                                           KS_eigenvalue, ovlp_3KS, i_state, i_spin, n_homo,&
                                           n_first, do_scalapack)

      type(cd_environment_type)                                :: gw_cd
      complex(kind=8), dimension(:), intent(inout)             :: aux_coulomb_residue_cmplx
      real(kind=8), dimension(:,:), intent(in)                 :: occ_numbers
      real(kind=8), dimension(:,:), intent(in)                 :: KS_eigenvalue
      real(kind=8), dimension(:,:,:,:), intent(in)             :: ovlp_3KS
      integer, dimension(:), intent(in)                        :: n_homo
      integer, dimension(:), intent(in)                        :: n_first
      integer, intent(in)                                      :: i_state, i_spin
      logical, intent(in), optional                            :: do_scalapack

      character(*), parameter :: func = 'calc_Wmn_for_real_frequencies'
 
      logical                                                  :: do_real_freq, my_do_scalapack
      integer                                                  :: i_freq_real, j_state, info 
      
      real(kind=8)                                             :: my_omega, my_eta
      real(kind=8), dimension(2)                               :: tmp_time
      real(kind=8), dimension(:,:), allocatable                :: screened_coulomb
      real(kind=8), dimension(:,:), allocatable                :: screened_coulomb_imag
      real(kind=8), dimension(:), allocatable                  :: aux_coulomb_residue_real
      real(kind=8), dimension(:), allocatable                  :: aux_coulomb_residue_imag
      type(matrix_structure_type)                              :: fmstruct

      tmp_time(:) = 0.d0 
      call get_timestamps(tmp_time(1), tmp_time(2))

      my_do_scalapack = .true.
      if(present(do_scalapack)) then
        my_do_scalapack = do_scalapack
      endif

      allocate(screened_coulomb(max(1,max_row_2d), max(1,max_col_2d)), stat=info)
      call check_allocation(info, 'screened_coulomb', func)
      allocate(screened_coulomb_imag(max(1,max_row_2d), max(1,max_col_2d)), stat=info)
      call check_allocation(info, 'screened_coulomb_imag', func)

      allocate(aux_coulomb_residue_real(gw_cd%num_residues),stat=info)
      call check_allocation(info, 'aux_coulomb_residue_real', func)
      allocate(aux_coulomb_residue_imag(gw_cd%num_residues),stat=info)
      call check_allocation(info, 'aux_coulomb_residue_imag', func)

      aux_coulomb_residue_real = 0.d0
      aux_coulomb_residue_imag = 0.d0
      call create_matrix_structure(fmstruct, n_basbas, n_basbas, max_row_2d, max_col_2d, &
                                   nb_aux_2d, nb_aux_2d, my_blacs_ctxt_aux_2d)

      do i_freq_real = 1, gw_cd%num_residues

         !*** contour deformation, where we have a real frequency
         my_omega = gw_cd%real_freq(i_freq_real)
         my_eta = gw_cd%eta
         do_real_freq = .TRUE.

         !*** calculate polarisability and from that the screened coulomb interaction 
         !*** in the auxiliary basis: W_PQ(e_m-e_n)
         call calc_screened_coulomb_matrix_WPQ(screened_coulomb, screened_coulomb_imag,fmstruct,&
                                              my_omega, my_eta, occ_numbers, KS_eigenvalue, ovlp_3KS,&
                                              n_homo, n_first, do_real_freq, my_do_scalapack)
 
         !*** calculate W_mn(e_m-e_n): screend coulomb interaction in MO basis for the residues 
         j_state = gw_cd%residue_from_freq(i_freq_real) 
         call calc_Wmn_from_WPQ_real_freq(screened_coulomb,ovlp_3KS,aux_coulomb_residue_real,&
                                         i_state, j_state, i_spin, i_freq_real, my_do_scalapack)
         if(flag_calc_spectral_func.or.flag_print_self_energy) then
            call calc_Wmn_from_WPQ_real_freq(screened_coulomb_imag,ovlp_3KS,aux_coulomb_residue_imag,&
                                             i_state,j_state, i_spin, i_freq_real, my_do_scalapack)
         endif 
         aux_coulomb_residue_cmplx(i_freq_real) = cmplx(aux_coulomb_residue_real(i_freq_real),&
                                                        aux_coulomb_residue_imag(i_freq_real), kind=8)

      enddo !i_freq_real

      if (allocated (aux_coulomb_residue_real)) then
        deallocate (aux_coulomb_residue_real)
      endif
      if (allocated (aux_coulomb_residue_imag)) then
        deallocate (aux_coulomb_residue_imag)
      endif
      if (allocated (screened_coulomb)) then
        deallocate (screened_coulomb)
      endif
      if (allocated (screened_coulomb_imag)) then
        deallocate (screened_coulomb_imag)
      endif

      call get_times(tmp_time(1), tmp_time(2), time_Wmn_realfreq_cd, clock_time_Wmn_realfreq_cd, &
                     unsynced=.TRUE.)

  end subroutine calc_Wmn_for_real_frequencies
! **************************************************************************************************
!> brief construct screened coulomb interaction Wc in auxiliary basisi {PQ}, i.e.(1-Pi)_PQ^-1 -V
!  o screened_coulomb -- real part of screened coulomb interaction
!  o screened_coulomb_imag -- imaginary part of screened coulomb interaction
!  o fmstruct -- structure of full distrubted matrix
!  o omega -- frequency 
!  o eta -- broadening parameter
!  o occ_number -- occupation number for each state and each spin, dim: n_states,n_spin
!  o KS_eigenvalue -- KS/HF eigenvalues of the single-particle calculation, dim: n_states,n_spin
!  o ovlp_3KS -- transformed 3-center integrals, dim: (n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
!  o n_homo -- the HOMO level for each spin channel
!  o n_first -- the first orbital which is NOT fully occupied
!  o do_real_freq -- whether we have a real or imaginary frequency
!  o do_scalapack -- if we do scalapack or not 
! **************************************************************************************************
  subroutine calc_screened_coulomb_matrix_WPQ(screened_coulomb, screened_coulomb_imag,fmstruct,&
                                              omega, eta, occ_numbers, KS_eigenvalue, ovlp_3KS,&
                                              n_homo, n_first, do_real_freq, do_scalapack)

      real(kind=8), dimension(:,:), intent(inout)              :: screened_coulomb
      real(kind=8), dimension(:,:), intent(inout)              :: screened_coulomb_imag
      type(matrix_structure_type)                              :: fmstruct
      real(kind=8), intent(in)                                 :: omega
      real(kind=8), intent(in)                                 :: eta
      real(kind=8), dimension(:,:), intent(in)                 :: occ_numbers
      real(kind=8), dimension(:,:), intent(in)                 :: KS_eigenvalue
      real(kind=8), dimension(:,:,:,:), intent(in)             :: ovlp_3KS
      integer, dimension(:), intent(in)                        :: n_homo
      integer, dimension(:), intent(in)                        :: n_first
      logical, intent(in)                                      :: do_real_freq
      logical, intent(in)                                      :: do_scalapack
      
      character(*), parameter :: func = 'calc_screened_coulomb_matrix_WPQ'

      integer                                                  :: irow, icol, info
      real(kind=8), dimension(2)                               :: tmp_time, tmp_time_pol
      real(kind=8), dimension(:,:), allocatable                :: polar_freq
      complex(kind=8), dimension(:,:), allocatable             :: screened_coulomb_cmplx

      tmp_time = 0.0d0
      call get_timestamps(tmp_time(1), tmp_time(2))

      tmp_time_pol(1) = time_polar_realfreq_cd
      tmp_time_pol(2) = clock_time_polar_realfreq_cd
      allocate(polar_freq(max(1,max_row_2d), max(1,max_col_2d)), stat=info)
      call check_allocation(info, 'polar_freq', func)
      allocate(screened_coulomb_cmplx(max(1,max_row_2d), max(1,max_col_2d)), stat=info)
      call check_allocation(info, 'screened_coulomb_cmplx', func)

      !*** evaluate the polarisability at frequency point i_freq_real 
      call  evaluate_polarisability_freq_2 (1, n_homo, n_first, occ_numbers, &
            omega, KS_eigenvalue, ovlp_3KS, polar_freq, do_real_freq,&
            eta, do_imaginary=.FALSE., do_scalapack=do_scalapack,&
            time_polar=tmp_time_pol)

      !*** now calculate the screened Coulomb interaction W
      call set_matrix_elements(screened_coulomb, fmstruct, off_val=0.0d0, diag_val=1.0d0, &
                               do_scalapack=do_scalapack)
      screened_coulomb(:,:) = screened_coulomb(:,:) - polar_freq(:,:)*2.d0/dble(n_spin)

      if(ABS(eta) > 1.0d-10) then
        do icol = 1, SIZE(screened_coulomb,2)
           do irow = 1, SIZE(screened_coulomb,1)
              screened_coulomb_cmplx(irow,icol) = cmplx(screened_coulomb(irow,icol),0.0d0,kind=8)
           enddo
        enddo
        call  evaluate_polarisability_freq_2 (1, n_homo, n_first, occ_numbers, omega, &
              KS_eigenvalue, ovlp_3KS, polar_freq, do_real_freq, eta,&
              do_imaginary=.TRUE., do_scalapack=do_scalapack, time_polar=tmp_time_pol)
        screened_coulomb(:,:) = - polar_freq(:,:)*2.d0/dble(n_spin)
        do icol = 1, SIZE(screened_coulomb,2)
           do irow = 1, SIZE(screened_coulomb,1)
             screened_coulomb_cmplx(irow,icol) = cmplx(real(screened_coulomb_cmplx(irow,icol),kind=8),  &
                                                       screened_coulomb(irow,icol),kind=8)
           enddo
        enddo
        call invert_cmplx_matrix(screened_coulomb_cmplx, fmstruct, do_scalapack)
        screened_coulomb = real(screened_coulomb_cmplx, kind=8)
        screened_coulomb_imag = aimag(screened_coulomb_cmplx)
      else 
        call invert_real_matrix(screened_coulomb, fmstruct, do_scalapack)
      endif

      !*** polar_freq isn't needed any more, set it to unity matrix
      call set_matrix_elements(polar_freq, fmstruct, off_val=0.0d0, diag_val=1.0d0, &
                               do_scalapack=do_scalapack)
      !*** subtract the bare Coulomb potential from W, Wc = W-V 
      screened_coulomb(:,:) = screened_coulomb(:,:) - polar_freq(:,:)

      !*** deallocate stuff
      if (allocated (polar_freq)) then
        deallocate (polar_freq)
      endif

      if (allocated (screened_coulomb_cmplx)) then
        deallocate (screened_coulomb_cmplx)
      endif

      time_polar_realfreq_cd = tmp_time_pol(1) 
      clock_time_polar_realfreq_cd = tmp_time_pol(2)
      call get_times(tmp_time(1), tmp_time(2), time_WPQ_realfreq_cd, clock_time_WPQ_realfreq_cd,&
                     unsynced=.TRUE.)

  end subroutine calc_screened_coulomb_matrix_WPQ

! **************************************************************************************************
!> brief calculates Wmn from (1-Pi(e_m-e_n))_{PQ}^-1, i.e. for real frequencies e_m-e_n
!  o screened_coulomb -- screened Coulomb interaction Wc = (1-Pi(e_m-e_n))_{PQ}^-1 -V
!  o ovlp_3KS -- transformed 3-center integrals, dim: (n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
!  o aux_coulomb_residue -- screened coulomb interaction in the MO basis: Wmn(e_m-e_n)
!  o i_state -- ith KS state
!  o j_state -- jth KS state
!  o i_spin -- ith spin
!  o i_freq_real -- real frequency index
!  o do_scalapack -- whether or not screened_coulomb is distributed
! **************************************************************************************************
  subroutine calc_Wmn_from_WPQ_real_freq(screened_coulomb,ovlp_3KS,aux_coulomb_residue,i_state,&
                                         j_state, i_spin, i_freq_real, do_scalapack)

      real(kind=8), dimension(:,:), intent(in)                 :: screened_coulomb
      real(kind=8), dimension(:,:,:,:), intent(in)             :: ovlp_3KS
      real(kind=8), dimension(:), intent(out)                  :: aux_coulomb_residue
      integer, intent(in)                                      :: i_state, j_state
      integer, intent(in)                                      :: i_freq_real, i_spin
      logical, intent(in)                                      :: do_scalapack

      character(*), parameter :: func = 'calc_Wmn_from_WPQ_real_freq'

      integer                                                  :: i_c, i_r, len_c, len_r, np_c,&
                                                                  np_r, np, i_c_loc, i_r_loc
      integer                                                  :: i_cs, i_rs, len_cs, len_rs, &
                                                                  n_blk, info, mpierr
      integer                                                  :: i_state_loc, j_state_loc
      real(kind=8), dimension(:,:), allocatable                :: tmp1
      real(kind=8), dimension(:,:), allocatable                :: screened_coulomb_sub
      real(kind=8), dimension(:), allocatable                  :: ovlp_3KS_ijstate
      real(kind=8), dimension(:), allocatable                  :: aux_ovlp3KS_ijstate

      !*** Block size for matrix multiplies, must be a multiple of nb_aux_2d
      n_blk = (255/nb_aux_2d+1)*nb_aux_2d

      allocate (tmp1(nb_aux_2d,nb_aux_2d), stat=info)
      call check_allocation(info, 'tmp1', func)

      allocate(screened_coulomb_sub(n_blk,n_blk), stat=info)
      call check_allocation(info, 'screened_coulomb_sub', func)

      allocate(ovlp_3KS_ijstate(n_basbas), stat=info)
      call check_allocation(info, 'ovlp3KS_ijstate', func)

      allocate(aux_ovlp3KS_ijstate(n_blk), stat=info)
      call check_allocation(info, 'aux_ovlp3KS_ijstate', func)
   
      !*** loop over the blocks of the distributed matrix screened_coulomb
      do i_c = 0, n_basbas-1, n_blk
      do i_r = 0, n_basbas-1, n_blk

        ! size of current block
        len_c = min(n_basbas-i_c, n_blk)
        len_r = min(n_basbas-i_r, n_blk)

        ! broadcast current block

        do i_rs = i_r, min(i_r+n_blk-1, n_basbas-1), nb_aux_2d
        do i_cs = i_c, min(i_c+n_blk-1, n_basbas-1), nb_aux_2d

           ! size of current block
           len_cs = min(n_basbas-i_cs, nb_aux_2d)
           len_rs = min(n_basbas-i_rs, nb_aux_2d)

           if(do_scalapack) then
             ! processor col/row of current block
             np_c = MOD(i_cs/nb_aux_2d,npcol_aux_2d)
             np_r = MOD(i_rs/nb_aux_2d,nprow_aux_2d)
             np = global_id(np_r,np_c) ! global owner of current block

             ! local offset of current block
             i_c_loc = (i_cs/(nb_aux_2d*npcol_aux_2d))*nb_aux_2d
             i_r_loc = (i_rs/(nb_aux_2d*nprow_aux_2d))*nb_aux_2d

             if(np==myid) &
                tmp1(1:len_rs,1:len_cs) = screened_coulomb(i_r_loc+1:i_r_loc+len_rs,i_c_loc+1:i_c_loc+len_cs)
             call mpi_bcast(tmp1,nb_aux_2d*len_cs,MPI_REAL8,np,mpi_comm_global,mpierr)

             screened_coulomb_sub(i_rs-i_r+1:i_rs-i_r+len_rs,i_cs-i_c+1:i_cs-i_c+len_cs) = tmp1(1:len_rs,1:len_cs)
           else
             screened_coulomb_sub(i_rs-i_r+1:i_rs-i_r+len_rs,i_cs-i_c+1:i_cs-i_c+len_cs)=&
                screened_coulomb(i_rs+1:i_rs+len_rs,i_cs+1:i_cs+len_cs)
           endif
         enddo
         enddo

         !*** loop over KS states which are of interest. For each real frequency e_m-e_n we only
         !*** multiply with the corresponding O^{nm}_P matrix element. (For imaginary frequencies
         !*** we multiply with all matrix elements)
         if(own_dim2_o3KS(i_state) == myp2_o3KS) then
            i_state_loc = loc_dim2_o3KS(i_state)
            if(own_dim1_o3KS(j_state) == myp1_o3KS) then
               j_state_loc = loc_dim1_o3KS(j_state)
               ovlp_3KS_ijstate(:) = ovlp_3KS(:,j_state_loc,i_state_loc,i_spin)
               call dgemm('N', 'N', len_r, 1, len_c, 1.0d0, &
                    screened_coulomb_sub, & 
                    ubound(screened_coulomb_sub,1), &
                    ovlp_3KS_ijstate(i_c+1), &
                    ubound(ovlp_3KS_ijstate,1), 0.d0, &
                    aux_ovlp3KS_ijstate, ubound(aux_ovlp3KS_ijstate,1) &
                    )
   
                  aux_coulomb_residue(i_freq_real) = &
                     aux_coulomb_residue(i_freq_real) + &
                     dot_product(aux_ovlp3KS_ijstate(1:len_r), &
                        ovlp_3KS_ijstate(i_r+1:i_r+len_r))
            endif
         endif

      enddo 
      enddo 

      if (allocated (aux_ovlp3KS_ijstate)) then
        deallocate (aux_ovlp3KS_ijstate)
      endif
      if (allocated (ovlp_3KS_ijstate)) then
        deallocate (ovlp_3KS_ijstate)
      endif
      if(allocated(screened_coulomb_sub)) deallocate(screened_coulomb_sub)
      if(allocated(tmp1)) deallocate(tmp1)

      return 

  end subroutine calc_Wmn_from_WPQ_real_freq

! **************************************************************************************************
!> brief calculates the Wmn matrix elements for complex frequencies i\omega fully complex, i.e.
!>       we also include a broadening parameter in the frequency integral
!  o gw_cd  --  contour deformation environment
!  o Wmn_matrix_freq_cmplx -- Wmn matrix elements, dim: (ndim2_o3KS,ndim1_o3KS,n_full_freq,n_spin)
!  o i_state -- ith KS state
!  o i_spin -- ith spin
!  o n_high_state -- highest KS/HF eigenstate for self-energy correction
!  o n_homo -- the HOMO level, i.e., the number of occupied state 
!  o occ_number -- occupation number for each state and each spin, dim: n_states,n_spin
!  o n_full_freq -- the number of frequency points for numerical integration
!  o omega_full(n_freq) -- the Gauss-Legendre frequency grid for numerical integration,
!                          dim:n_full_freq
!  o womega_full(n_freq) -- the weigth of the Gauss-Legendre frequency grid for the numerical
!                           interaction, dim: n_full_freq
!  o qp_energy --  quasi particle energies
!  o KS_eigenvalue -- KS/HF eigenvalues of the single-particle calculation, dim: n_states,n_spin
!  o ovlp_3KS -- transformed 3-center integrals, dim: (n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
! **************************************************************************************************
  subroutine calc_Wmn_for_imag_freq_fully_complex(gw_cd, Wmn_matrix_freq_cmplx, n_low_state,&
                                                  n_high_state, n_homo, occ_numbers, n_full_freq,&
                                                  omega_full, KS_eigenvalue, ovlp_3KS) 
 

      type(cd_environment_type)                                :: gw_cd
      complex(kind=8), dimension(:,:,:,:), intent(inout)       :: Wmn_matrix_freq_cmplx 
      integer, intent(in)                                      :: n_low_state, n_high_state
      integer, dimension(:), intent(in)                        :: n_homo
      real(kind=8), dimension(:,:), intent(in)                 :: occ_numbers
      integer, intent(in)                                      :: n_full_freq
      real(kind=8), dimension(:), intent(in)                   :: omega_full
      real(kind=8), dimension(:,:), intent(in)                 :: KS_eigenvalue
      real(kind=8), dimension(:, :, :, :), intent(in)          :: ovlp_3KS 

      character(*), parameter :: func = 'calc_Wmn_for_imag_freq_fully_complex'
 
      logical                                                  :: do_real_freq
      logical                                                  :: do_scalapack
      integer                                                  :: i_freq, i_state,&
                                                                  i_spin, info 
      integer, dimension(n_spin)                               :: n_first
      
      real(kind=8)                                             :: my_eta, my_omega
      real(kind=8), dimension(:,:), allocatable                :: screened_coulomb
      real(kind=8), dimension(:,:), allocatable                :: screened_coulomb_imag
      real(kind=8), dimension(:,:,:,:), allocatable            :: aux_coulomb_real
      real(kind=8), dimension(:,:,:,:), allocatable            :: aux_coulomb_imag

      type(matrix_structure_type)                              :: fmstruct


      if(myid.eq.0)then
         write(use_unit,'(2A)')" ------------------------------------------------------"
         write(use_unit,*)
         write(use_unit,"(T3,A)") "Calculate Wmn(i\omega) matrix elements fully complex:" 
         write(use_unit,*)
         write(use_unit,'(2A)')" ------------------------------------------------------"
      endif

      n_first(:) = 1
      do i_spin = 1, n_spin
       do i_state = 1, n_states
        if (abs(occ_numbers(i_state,i_spin)-dble(2/n_spin)) &
                        .lt.1.d-8) then
         n_first(i_spin)= i_state + 1
        endif 
       enddo 
       if(n_first(i_spin) .gt. n_states) then
         n_first(i_spin) = n_states
       endif
      enddo

      allocate(screened_coulomb(max(1,max_row_2d), max(1,max_col_2d)), stat=info)
      call check_allocation(info, 'screened_coulomb', func)
      screened_coulomb(:,:) = 0.0d0

      allocate(screened_coulomb_imag(max(1,max_row_2d), max(1,max_col_2d)), stat=info)
      call check_allocation(info, 'screened_coulomb_imag', func)
      screened_coulomb_imag(:,:) = 0.0d0

      allocate(aux_coulomb_real(ndim2_o3KS,ndim1_o3KS,n_full_freq,n_spin), stat=info)
      call check_allocation(info, 'aux_coulomb_real', func)
      aux_coulomb_real(:,:,:,:) = 0.0d0

      allocate(aux_coulomb_imag(ndim2_o3KS,ndim1_o3KS,n_full_freq,n_spin), stat=info)
      call check_allocation(info, 'aux_coulomb_imag', func)
      aux_coulomb_imag(:,:,:,:) = 0.0d0

      call create_matrix_structure(fmstruct, n_basbas, n_basbas, max_row_2d, max_col_2d, &
                                   nb_aux_2d, nb_aux_2d, my_blacs_ctxt_aux_2d)


      my_eta = gw_cd%eta
      do_scalapack = .TRUE.
      do_real_freq = .FALSE.

      do i_freq = 1, n_full_freq

          if(myid.eq.0) then
            write(use_unit,'(A,I4,f20.6)') " | i_freq", i_freq, &
                   omega_full(i_freq)
          endif

          my_omega = omega_full(i_freq)
          call  calc_screened_coulomb_matrix_WPQ(screened_coulomb, screened_coulomb_imag,fmstruct,&
                                                 my_omega, my_eta, occ_numbers, KS_eigenvalue, ovlp_3KS,&
                                                 n_homo, n_first, do_real_freq, do_scalapack)

          call calc_Wmn_from_WPQ_imag_freq(screened_coulomb, ovlp_3KS, aux_coulomb_real, i_freq,&
                                           n_low_state, n_high_state)
          if(ABS(my_eta) > 1.0d-10) then
            call calc_Wmn_from_WPQ_imag_freq(screened_coulomb_imag, ovlp_3KS, aux_coulomb_imag,&
                                             i_freq, n_low_state, n_high_state)
          endif
      enddo
  
      Wmn_matrix_freq_cmplx(:,:,:,:) = cmplx(aux_coulomb_real,aux_coulomb_imag,kind=8)

      if (allocated (screened_coulomb)) then
        deallocate (screened_coulomb)
      endif
      if (allocated (screened_coulomb_imag)) then
        deallocate (screened_coulomb_imag)
      endif
      if (allocated (aux_coulomb_real)) then
        deallocate (aux_coulomb_real)
      endif
      if (allocated (aux_coulomb_imag)) then
        deallocate (aux_coulomb_imag)
      endif


  end subroutine calc_Wmn_for_imag_freq_fully_complex

! **************************************************************************************************
!> brief calculates Wmn from (1-Pi(i\omega))_{PQ}^-1, i.e. for imaginary frequencies i\omega
!  o screened_coulomb -- screened Coulomb interaction Wc = (1-Pi(i\omega))_{PQ}^-1 -V
!  o ovlp_3KS -- transformed 3-center integrals, dim: (n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
!  o aux_coulomb_imag_ -- screened coulomb interaction in the MO basis: Wmn(i\omega)
!  o i_spin -- ith spin
!  o i_freq_imag -- imaginary frequency index
!  o n_low_state -- lowest KS/HF eigenstate for self-energy correction
!  o n_high_state -- highest KS/HF eigenstate for self-energy correction
! **************************************************************************************************
  subroutine calc_Wmn_from_WPQ_imag_freq(screened_coulomb,ovlp_3KS, aux_coulomb_imag_freq, i_freq_imag, &
                                n_low_state, n_high_state)

      real(kind=8), dimension(:,:), intent(in)                 :: screened_coulomb
      real(kind=8), dimension(:,:,:,:), intent(in)             :: ovlp_3KS
      real(kind=8), dimension(:,:,:,:), intent(out)            :: aux_coulomb_imag_freq
      integer, intent(in)                                      :: i_freq_imag
      integer, intent(in)                                      :: n_low_state, n_high_state

      character(*), parameter :: func = 'calc_Wmn_from_WPQ_imag_freq'

      integer                                                  :: i_c, i_r, len_c, len_r, np_c,&
                                                                  np_r, np, i_c_loc, i_r_loc
      integer                                                  :: i_cs, i_rs, len_cs, len_rs, &
                                                                  i_spin, i_state, n_blk, info, &
                                                                  mpierr
      integer                                                  :: i_state_loc, j_state_loc
      real(kind=8), dimension(:,:), allocatable                :: tmp1
      real(kind=8), dimension(:,:), allocatable                :: screened_coulomb_sub
      real(kind=8), dimension(:,:,:), allocatable              :: aux_ovlp3KS

      !*** Block size for matrix multiplies, must be a multiple of nb_aux_2d
      n_blk = (255/nb_aux_2d+1)*nb_aux_2d

      allocate (tmp1(nb_aux_2d,nb_aux_2d), stat=info)
      call check_allocation(info, 'tmp1', func)

      allocate(screened_coulomb_sub(n_blk,n_blk), stat=info)
      call check_allocation(info, 'screened_coulomb_sub', func)

      allocate(aux_ovlp3KS(n_blk,ndim1_o3KS,ndim2_o3KS), stat=info)
      call check_allocation(info, 'aux_ovlp3KS', func)

      !*** loop over the blocks of the distributed matrix screened_coulomb
      do i_c = 0, n_basbas-1, n_blk
      do i_r = 0, n_basbas-1, n_blk

        ! size of current block
        len_c = min(n_basbas-i_c, n_blk)
        len_r = min(n_basbas-i_r, n_blk)

        ! broadcast current block

        do i_rs = i_r, min(i_r+n_blk-1, n_basbas-1), nb_aux_2d
        do i_cs = i_c, min(i_c+n_blk-1, n_basbas-1), nb_aux_2d

           ! size of current block
           len_cs = min(n_basbas-i_cs, nb_aux_2d)
           len_rs = min(n_basbas-i_rs, nb_aux_2d)

           ! processor col/row of current block
           np_c = MOD(i_cs/nb_aux_2d,npcol_aux_2d)
           np_r = MOD(i_rs/nb_aux_2d,nprow_aux_2d)
           np = global_id(np_r,np_c) ! global owner of current block

           ! local offset of current block
           i_c_loc = (i_cs/(nb_aux_2d*npcol_aux_2d))*nb_aux_2d
           i_r_loc = (i_rs/(nb_aux_2d*nprow_aux_2d))*nb_aux_2d

           if(np==myid) &
              tmp1(1:len_rs,1:len_cs) = screened_coulomb(i_r_loc+1:i_r_loc+len_rs,i_c_loc+1:i_c_loc+len_cs)
           call mpi_bcast(tmp1,nb_aux_2d*len_cs,MPI_REAL8,np,mpi_comm_global,mpierr)

           screened_coulomb_sub(i_rs-i_r+1:i_rs-i_r+len_rs,i_cs-i_c+1:i_cs-i_c+len_cs) = tmp1(1:len_rs,1:len_cs)
         enddo
         enddo

         !*** loop over KS states which are of interest.
         do i_spin = 1, n_spin
           call dgemm('N', 'N', len_r, ndim1_o3KS*ndim2_o3KS, &
                 len_c, 1.0d0, &
                 screened_coulomb_sub, & 
                 ubound(screened_coulomb_sub,1), &
                 ovlp_3KS(i_c+1,1,1,i_spin), &
                 ubound(ovlp_3KS,1), 0.d0, &
                 aux_ovlp3KS, ubound(aux_ovlp3KS,1) &
               )
           do i_state = n_low_state,n_high_state, 1

              if(own_dim2_o3KS(i_state) /= myp2_o3KS) cycle
              i_state_loc = loc_dim2_o3KS(i_state)
   
              do j_state_loc = 1, ndim1_o3KS
                 aux_coulomb_imag_freq(i_state_loc,j_state_loc,i_freq_imag,i_spin) = &
                    aux_coulomb_imag_freq(i_state_loc,j_state_loc,i_freq_imag,i_spin) + &
                    dot_product(aux_ovlp3KS(1:len_r,j_state_loc,i_state_loc), &
                       ovlp_3KS(i_r+1:i_r+len_r,j_state_loc,i_state_loc,i_spin))
              enddo
           enddo ! end of i_state
         enddo ! end of i_state

      enddo 
      enddo 

      if (allocated (aux_ovlp3KS)) then
        deallocate (aux_ovlp3KS)
      endif

      if(allocated(screened_coulomb_sub)) deallocate(screened_coulomb_sub)
      if(allocated(tmp1)) deallocate(tmp1)

  end subroutine calc_Wmn_from_WPQ_imag_freq

! **************************************************************************************************
!> brief calculation of W(i\omega) needed to calculate the integral term in the CD  
!  o gw_cd  --  contour deformation environment
!  o Wmn_imag_freq_cmplx -- Wmn matrix elements, dim: (ndim2_o3KS,ndim1_o3KS,n_full_freq,n_spin)
!  o self_energy_omega -- self energy for complex frequencies i\omega
!  o n_low_state  -- the lowest KS/HF eigenstate for self-energy correction
!  o n_high_state -- the highest KS/HF eigenstate for self-energy correction
!  o occ_numbers -- the occupation number of the electrons for each eigenstate and each spin,
!                    dim: (n_states, n_spin) 
! o  omega(n_freq) -- the Gauss-Legendre frequency grid for the self-energy
! o  omega_full(n_freq) -- the Gauss-Legendre frequency grid for the screened Coulomb interaction
! o  womega_full(n_freq) -- the weigth of the Gauss-Legendre freq. grid for the scree. Coulomb inter.
! o  chemical_potential -- real number, the chemical potential of the system
! o  KS_eigenvalue -- the eigenvalues of the single-particle calculation. For DFT calculation,
!                     this is the KS eigenvalue, but for HF calculation, this is then HF eigenval.
! o  KS_eigenvector -- the eigenvector of the single-particle calculation
! o  ovlp_3KS -- this is the transformed 3-cener overlap integration. Now two orbitals of
!                them are KS ones, and one is the auxiliary basis.
! o flag_coh -- logic, if ture, use the Coulomb-hole accelerator (not functioning for the time being)
! **************************************************************************************************
   subroutine evaluate_Wmn_imag_freq_for_cd(gw_cd, Wmn_imag_freq_cmplx, self_energy_omega,&
                                            n_low_state, n_high_state, ndim2_o3KS, &
                                            ndim1_o3KS, n_spin, n_homo, occ_numbers,&
                                            n_freq, n_full_freq, omega_grid, omega_full_grid,&
                                            womega_full, chemical_potential_spin, KS_eigenvalue,&
                                            KS_eigenvector, ovlp_3KS, out_self_energy, partition_tab, &
                                            l_shell_max, flag_coh, tot_time_polar, tot_time_sigma)
                     

      type(cd_environment_type)                                :: gw_cd
      complex(kind=8), dimension(:, :, :, :), intent(inout)    :: Wmn_imag_freq_cmplx 
      complex(kind=8), dimension(:, :, :), intent(inout)       :: self_energy_omega 
      integer, intent(in)                                      :: n_low_state, n_high_state
      integer, intent(in)                                      :: ndim2_o3KS,ndim1_o3KS, n_spin
      integer, dimension(:), intent(in)                        :: n_homo
      real(kind=8), dimension(:,:), intent(in)                 :: occ_numbers
      integer, intent(in)                                      :: n_freq, n_full_freq
      real(kind=8), dimension(:), intent(in)                   :: omega_grid
      real(kind=8), dimension(:), intent(in)                   :: omega_full_grid
      real(kind=8), dimension(:), intent(in)                   :: womega_full
      real(kind=8), dimension(:), intent(in)                   :: chemical_potential_spin   
      real(kind=8), dimension(:,:), intent(in)                 :: KS_eigenvalue
      real(kind=8), dimension(:,:,:), intent(in)               :: KS_eigenvector
      real(kind=8), dimension(:,:,:,:), intent(in)             :: ovlp_3KS
      logical, intent(in)                                      :: out_self_energy
      real(kind=8), dimension(:), intent(in)                   :: partition_tab
      integer, dimension(:), intent(in)                        :: l_shell_max
      logical, intent(in)                                      :: flag_coh
      real(kind=8), intent(inout)                              :: tot_time_polar, tot_time_sigma
      

      real(kind=8), dimension(:,:,:,:), allocatable            :: Wmn_imag_freq


      if(gw_cd%full_cmplx_sigma) then
        !*** do analytic continuation for the states not calculated with CD, but
        !*** don't store the Wmn matrix elements
        call evaluate_self_energy_freq_2 &
        (n_low_state, n_high_state, &
        n_homo, occ_numbers, &
        n_freq, n_full_freq, &
        omega_grid, omega_full_grid, &
        womega_full, chemical_potential_spin, &
        KS_eigenvalue,KS_eigenvector, &
        partition_tab, l_shell_max, &
        ovlp_3KS, self_energy_omega, out_self_energy, &
        flag_coh, tot_time_sigma, tot_time_polar)
        !*** calculate the Wmn(i\omega) fully complex
        call calc_Wmn_for_imag_freq_fully_complex &
             (gw_cd, Wmn_imag_freq_cmplx, n_low_state, n_high_state,&
             n_homo, occ_numbers, n_full_freq, omega_full_grid, &
             KS_eigenvalue, ovlp_3KS)
      else
        !*** calculate self-energy for complex frequencies --> analytic continuation for
        !    states that are not corrected with CD, store also the matrix elements Wmn(i\omega) 
        allocate(Wmn_imag_freq(ndim2_o3KS,ndim1_o3KS,n_full_freq,n_spin))
        call evaluate_self_energy_freq_2 &
        (n_low_state, n_high_state, &
        n_homo, occ_numbers, &
        n_freq,n_full_freq, &
        omega_grid, omega_full_grid, &
        womega_full, chemical_potential_spin, &
        KS_eigenvalue,KS_eigenvector, &
        partition_tab, l_shell_max, &
        ovlp_3KS, self_energy_omega, out_self_energy, &
        flag_coh, tot_time_sigma, tot_time_polar, &
        Wmn_imag_freq &
        )
        Wmn_imag_freq_cmplx = cmplx(Wmn_imag_freq,0.0d0,kind=8)
        deallocate(Wmn_imag_freq)
      endif

   end subroutine evaluate_Wmn_imag_freq_for_cd
end module evaluate_self_energy_freq 
