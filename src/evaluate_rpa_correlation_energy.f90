!****s* FHI-aims/evaluate_rpa_correlation_energy
!  NAME
!   evaluate_rpa_correlation_energy
!  SYNOPSIS

      subroutine evaluate_rpa_correlation_energy &
           (n_electrons, &
            n_low_state, n_high_state, &
            n_homo, n_lumo, &
            occ_numbers, n_full_freq, &
            omega_full, womega_full, &
            chemical_potential, &
            KS_eigenvalue, KS_eigenvector, &
            ovlp_3KS, &
            rpa_c_energy &
           )

!  PURPOSE
!  Subroutine evaluate_rpa_correlation_energy evaluates the correlation
!  energy at the RPA level using the adiabatic connection fluctuation
!  dissipation theorem.
!
!  E_RPA = 1/2pi \int dw { ln(det(1-v_times_polar)) + tr(v_times_polar) }

! USES
      use runtime_choices
      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use synchronize_mpi
      use timing
      use localorb_io, only: use_unit
      use evaluate_polarisability_freq, only: evaluate_polarisability_freq_0

      implicit none

! ARGUMENTS 

      integer :: n_full_freq
      integer :: n_low_state
      integer :: n_high_state
      integer :: n_lumo(n_spin)
      integer :: n_homo(n_spin)

      real*8  :: n_electrons
      real*8  :: occ_numbers(n_states,n_spin)
      real*8  :: omega_full(n_full_freq)
      real*8  :: womega_full(n_full_freq)
      real*8  :: chemical_potential
      real*8  :: KS_eigenvalue(n_states,n_spin)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin)
      real*8  :: ovlp_3KS(n_loc_prodbas, n_states, n_high_state,n_spin)

!     output
      real*8  :: rpa_c_energy

! INPUTS
! o  n_full_freq -- integer number,
!            the number of frequency points for the screened Coulomb interaction W
! o  n_low_state  -- integer number,
!            the lowest KS/HF eigenstate taken into account in the polarisability calcua            ltions
! o  n_high_state -- integer number,
!            the highest KS/HF eigenstate. In the present case, n_high_state >= n_homo
!            should be fine. 
! o  n_electrons -- real number
!            the total number of electrons in the system
! o  occ_numbers -- real 2-dimentianal array of length (n_states, n_spin)
!            the occupation number of the electrons for each eigenstate and each spin
! o  omega_full(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the screened Coulomb interaction
! o  womega_full(n_freq) -- real array
!            the weigth of the Gauss-Legendre frequency grid for the screened Coulomb 
!            in teraction
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
!
! OUTPUT
! o  rpa_c_energy -- real number, the calculated RPA correlation energy
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

      real*8  det_v_times_polar
      real*8  trace_v_times_polar
      real*8  rpa_c_integrand
      real*8  rpa_c_integrand_along_ac_path(rpa_along_ac_path_grid)

!    local timing
      real*8  temp_time_rpa
      real*8  temp_clock_time_rpa

!     auxiliary matrices for Level 3 Blas matrix multiplications
!     n_first : the first orbital which is NOT fully occupied

      integer :: n_first(n_spin)

      real*8, dimension(:,:), allocatable :: polar_freq

!     timing

!     parameters of the fitting tails
!       real*8  s1, s2, omega_1, omega_2
!       real*8  alpha, beta, a, b

      character*50  filename

!     counters

      integer :: i_state
      integer :: i_freq
      integer :: i_spin
      integer :: i_index
      integer :: i_basis_1, i_basis_2

!     Igor, RPA potentials along AC path
      real*8  :: interval_ac_path      ! interval along AC path
      if (rpa_along_ac_path_grid .gt. 0) then
          interval_ac_path = real(rpa_along_ac_path_grid)**(-1)
          allocate(rpa_potentials_along_ac_path(rpa_along_ac_path_grid),stat=i_index)
          call check_allocation(i_index, 'rpa_potentials_along_ac_path')
          rpa_potentials_along_ac_path(:) = 0.0d0
      endif
!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"-----------------------------------", &
                  " ------------------------------------"
        write(use_unit,'(2X,A)') &
              "Start to calculate the RPA correlation energy  ... "
      endif

!     determine the highest occupied orbital level
!     such complication occurs when some of the orbitals are
!     either not fully occupied or not fully empty
      n_first(:) = 1
      do i_spin = 1, n_spin
       do i_state = 1, n_states
!        write(use_unit,'(I4,f18.10)')i_state, occ_numbers(i_state,i_spin) 
        if (abs(occ_numbers(i_state,i_spin)-dble(2/n_spin)) &
                         .lt.1.d-6) then
         n_first(i_spin)= i_state + 1
        endif
       enddo
       if(n_first(i_spin) .gt. n_states) then
         n_first(i_spin) = n_states
       endif
      enddo

      if(myid.eq.0) then
       write(use_unit,'(2X, A,A,4I5)') &
                 "HOMO and first non-fully-occupied", &
                 " orbitals:", n_homo(:), n_first(:)
       write(use_unit,*)
      endif

!      open(112, file='rpa_c_energy.dat')

      allocate(polar_freq(n_basbas, n_loc_prodbas),stat=i_index)
      call check_allocation(i_index, 'polar_freq                    ')


!   invert the bare Coulomb matrix for calculating the screened Coulomb interaction .

      time_rpa_corr = 0.d0
      clock_time_rpa_corr = 0.d0
      time_polar = 0.d0
      clock_time_polar = 0.d0

      rpa_c_energy = 0.d0 
      do i_freq = 1, n_full_freq, 1

        rpa_c_integrand = 0.d0 

        if(myid.eq.0) then
           write(use_unit,'(A,I4,f20.6)') " | i_freq ", i_freq, & 
                   omega_full(i_freq)
        endif

      call get_timestamps(temp_time_rpa, temp_clock_time_rpa )
!    evaluate the polarisability at frequency point i_freq
      call  evaluate_polarisability_freq_0 &
             ( n_low_state, n_homo, n_first, n_high_state, &
               occ_numbers, &
               omega_full(i_freq), &
               KS_eigenvalue, ovlp_3KS, polar_freq &
             )

!       if(i_freq.eq.1 .and. myid.eq.0) then
!          do i_basis_2 = 1, n_loc_prodbas, 1
!             i_basis_1 = map_prodbas(i_basis_2, myid+1) 
!             write(use_unit,'(2I4,2f20.10)') i_basis_1, i_basis_2, polar_freq(i_basis_1, i_basis_2)
!          enddo
!       endif

       call get_timestamps(rtime, clock_rtime)
       time_polar = time_polar + rtime - temp_time_rpa
       clock_time_polar = clock_time_polar + clock_rtime - temp_clock_time_rpa

       polar_freq(:,:) = polar_freq(:,:) * 2.d0/dble(n_spin)

! Hubbard correction, here simply a factor of two
! how to improve?
!        polar_freq(:,:) = polar_freq(:,:)/dble(n_spin)

!       call evaluate_rpa_integrand(polar_freq,rpa_c_integrand)

       if (rpa_along_ac_path_grid .eq. 0) then
           call evaluate_rpa_integrand(polar_freq,rpa_c_integrand)
       else
           !rpa_potentials_along_ac_path(:) = 0.0d0

           call evaluate_rpa_integrand_along_ac_path(&
               polar_freq,rpa_c_integrand,rpa_c_integrand_along_ac_path)

           !if (myid .eq. 0) then
           !   write(use_unit,'(2X,A,F16.8,F16.8)') 'debug 4',&
           !       rpa_c_integrand_along_ac_path(1), &
           !       rpa_c_integrand_along_ac_path(2)
           !endif
       endif


!       call sync_real_number(rpa_c_integrand) 

       rpa_c_energy = rpa_c_energy + &
                       rpa_c_integrand * womega_full(i_freq)

!      RPA potentials along AC path
       if (rpa_along_ac_path_grid .gt. 0) then
           do i_state = 1 , rpa_along_ac_path_grid
               rpa_potentials_along_ac_path(i_state) = &
                   rpa_potentials_along_ac_path(i_state) + &
                   rpa_c_integrand_along_ac_path(i_state) * womega_full(i_freq)
           enddo
       endif

!      write(81,'(I6,4f18.8)')i_freq, omega_full(i_freq), &
!                  womega_full(i_freq), rpa_c_integrand, rpa_c_energy

       call get_timestamps(temp_time_rpa, temp_clock_time_rpa )
       time_rpa_corr = time_rpa_corr +  temp_time_rpa - rtime
       clock_time_rpa_corr = clock_time_rpa_corr + temp_clock_time_rpa - clock_rtime
! end of loop over i_freq
      enddo
      call sync_timing(time_polar)
      call sync_timing(time_rpa_corr)

!      close(112)

      rpa_c_energy=rpa_c_energy/2.d0/pi
      if (rpa_along_ac_path_grid .gt. 0) then
          rpa_potentials_along_ac_path = rpa_potentials_along_ac_path/2.d0/pi
      endif

! Delta term correction
! This deals with the case when there are frational occupations and 
! the intra-orbtial excitations should be taken into account. The contribution
! is a delta function at zero frequency and a finite contribution when integrated
! out over frequency axis.

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"----------------------------------------------------", &
                  "-------------------------"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            "  RPA correlation energy :", rpa_c_energy, "Ha,", &
             rpa_c_energy*hartree, "eV"
!        write(use_unit,*)"----------------------------------------------------", &
!                 "-------------------------"
        write(use_unit,*)
      endif

      if (allocated (polar_freq)) then
        deallocate (polar_freq)
      endif

      return

      end subroutine evaluate_rpa_correlation_energy
!---------------------------------------------------------------------
!******
!****s* FHI-aims/evaluate_rpa_correlation_energy
!  NAME
!   evaluate_rpa_correlation_energy
!  SYNOPSIS

      subroutine evaluate_rpa_correlation_energy_2 &
           (n_electrons, &
            n_low_state, n_high_state, &
            n_homo, n_lumo, &
            occ_numbers, n_full_freq, &
            omega_full, womega_full, &
            chemical_potential, &
            KS_eigenvalue, KS_eigenvector, &
            ovlp_3KS, &
            rpa_c_energy &
           )

!  PURPOSE
!  Subroutine evaluate_rpa_correlation_energy evaluates the correlation
!  energy at the RPA level using the adiabatic connection fluctuation
!  dissipation theorem.
!
!  E_RPA = 1/2pi \int dw { ln(det(1-v_times_polar)) + tr(v_times_polar) }

! USES
      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use synchronize_mpi
      use timing
      use evaluate_polarisability_freq, only: evaluate_polarisability_freq_2
      use localorb_io, only: use_unit
      implicit none

! ARGUMENTS 

      integer :: n_full_freq
      integer :: n_low_state
      integer :: n_high_state
      integer :: n_lumo(n_spin)
      integer :: n_homo(n_spin)

      real*8  :: n_electrons
      real*8  :: occ_numbers(n_states,n_spin)
      real*8  :: omega_full(n_full_freq)
      real*8  :: womega_full(n_full_freq)
      real*8  :: chemical_potential
      real*8  :: KS_eigenvalue(n_states,n_spin)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin)
      real*8  :: ovlp_3KS(n_basbas, ndim1_o3KS, ndim2_o3KS,n_spin)

!     output
      real*8  :: rpa_c_energy

! INPUTS
! o  n_full_freq -- integer number,
!            the number of frequency points for the screened Coulomb interaction W
! o  n_low_state  -- integer number,
!            the lowest KS/HF eigenstate taken into account in the polarisability calcua            ltions
! o  n_high_state -- integer number,
!            the highest KS/HF eigenstate. In the present case, n_high_state >= n_homo
!            should be fine. 
! o  n_electrons -- real number
!            the total number of electrons in the system
! o  occ_numbers -- real 2-dimentianal array of length (n_states, n_spin)
!            the occupation number of the electrons for each eigenstate and each spin
! o  omega_full(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the screened Coulomb interaction
! o  womega_full(n_freq) -- real array
!            the weigth of the Gauss-Legendre frequency grid for the screened Coulomb 
!            in teraction
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
!
! OUTPUT
! o  rpa_c_energy -- real number, the calculated RPA correlation energy
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

      real*8  det_v_times_polar
      real*8  trace_v_times_polar
      real*8  rpa_c_integrand

!    local timing
      real*8  temp_time_rpa
      real*8  temp_clock_time_rpa

!     auxiliary matrices for Level 3 Blas matrix multiplications
!     n_first : the first orbital which is NOT fully occupied

      integer :: n_first(n_spin)

      real*8, dimension(:,:), allocatable :: polar_freq

!     timing

!     parameters of the fitting tails
!       real*8  s1, s2, omega_1, omega_2
!       real*8  alpha, beta, a, b

      character*50  filename

!     counters

      integer :: i_state
      integer :: i_freq
      integer :: i_spin
      integer :: i_index

!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"-----------------------------", &
                  "-----------------------------------------"
        write(use_unit,'(2X,A)') &
              "Start to calculate the RPA correlation energy  ... "
      endif

!     determine the highest occupied orbital level
!     such complication occurs when some of the orbitals are
!     either not fully occupied or not fully empty
      n_first(:) = 1
      do i_spin = 1, n_spin
       do i_state = 1, n_states
        if (abs(occ_numbers(i_state,i_spin)-dble(2/n_spin)) &
                         .lt.1.d-6) then
         n_first(i_spin)= i_state + 1
        endif
       enddo
       if(n_first(i_spin) .gt. n_states) then
         n_first(i_spin) = n_states
       endif
      enddo


      if(myid.eq.0) then
       write(use_unit,'(2X, A,A,4I5)') &
                    "HOMO and first non-fully-occupied", &
                    " orbitals:", n_homo(:), n_first(:)
       write(use_unit,*)
      endif

!      open(112, file='rpa_c_energy.dat')

      allocate(polar_freq(max_row_2d, max_col_2d),stat=i_index)
      call check_allocation(i_index, 'polar_freq                    ')

!   invert the bare Coulomb matrix for calculating the screened Coulomb interaction .

      time_rpa_corr = 0.d0
      clock_time_rpa_corr = 0.d0
      time_polar = 0.d0
      clock_time_polar = 0.d0

      rpa_c_energy = 0.d0 
      do i_freq = 1, n_full_freq, 1

        rpa_c_integrand = 0.d0 

        if(myid.eq.0) then
           write(use_unit,'(A,I4,f20.6)') " | i_freq ", i_freq, & 
                   omega_full(i_freq)
        endif

      call get_timestamps(temp_time_rpa, temp_clock_time_rpa )
!    evaluate the polarisability at frequency point i_freq
      call  evaluate_polarisability_freq_2 &
             ( n_low_state, n_homo, n_first, &
               occ_numbers, omega_full(i_freq), &
               KS_eigenvalue, ovlp_3KS, polar_freq &
             )


       call get_timestamps(rtime, clock_rtime)
       time_polar = time_polar + rtime - temp_time_rpa
       clock_time_polar = clock_time_polar + clock_rtime - temp_clock_time_rpa

       polar_freq(:,:) = polar_freq(:,:) * 2.d0/dble(n_spin)

! Hubaard correction, here simply a factor of two
! how to improve?
!        polar_freq(:,:) = polar_freq(:,:)/dble(n_spin)

       call evaluate_rpa_integrand_2(polar_freq,rpa_c_integrand)


!       call sync_real_number(rpa_c_integrand)     

       rpa_c_energy = rpa_c_energy + &
                       rpa_c_integrand * womega_full(i_freq)

!      write(81,'(I6,3f20.8)')i_freq, omega_full(i_freq), &
!                  womega_full(i_freq), rpa_c_integrand

       call get_timestamps(temp_time_rpa, temp_clock_time_rpa )
       time_rpa_corr = time_rpa_corr +  temp_time_rpa - rtime
       clock_time_rpa_corr = clock_time_rpa_corr + temp_clock_time_rpa - clock_rtime
! end of loop over i_freq
      enddo
      call sync_timing(time_polar)
      call sync_timing(time_rpa_corr)

!      close(112)

      rpa_c_energy=rpa_c_energy/2.d0/pi

! Delta term correction
! This deals with the case when there are frational occupations and 
! the intra-orbtial excitations should be taken into account. The contribution
! is a delta function at zero frequency and a finite contribution when integrated
! out over frequency axis.

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"--------------------------------------------", &
                  "----------------------------------"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            "  RPA correlation energy :", rpa_c_energy, "Ha,", &
             rpa_c_energy*hartree, "eV"
        write(use_unit,*)
      endif

      if (allocated (polar_freq)) then
        deallocate (polar_freq)
      endif

      return

      end subroutine evaluate_rpa_correlation_energy_2
!---------------------------------------------------------------------
!******
