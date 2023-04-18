!****** FHI-aims/evaluate_qpe_osrpa_energy
!  NAME
!   evaluate_qpe_osrpa_energy
!  SYNOPSIS

      subroutine evaluate_qpe_osrpa_energy &
           (n_electrons, &
            n_low_state, n_high_state, &
            n_homo, n_lumo, &
            occ_numbers, n_full_freq, &
            omega_full, womega_full, &
            chemical_potential, &
            KS_eigenvalue, KS_eigenvector, &
            ovlp_3KS, &
            rpa_c_energy_os, rpa_c_energy_ss, & 
            rpa_c_energy_total, rpa_c_energy_scs &
           )

!  PURPOSE
!  Subroutine evaluate_qpe_osrpa_energy evaluates the correlation
!  energy at the osPT2 level using the adiabatic connection fluctuation
!  dissipation theorem.
!
!  For the close-shell case chi_a=chi_b, the renormalized os-PT2 can be
!  calculated simply following the original RPA routine.
!  E_RPA = 1/2pi \int dw { ln(det(1-v_times_polar)) + tr(v_times_polar) }

! USES
      use runtime_choices
      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use synchronize_mpi
      use timing
      use evaluate_ospt2_polarizability_freq, only: evaluate_ospt2_polarizability_freq_0
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
      real*8  :: ovlp_3KS(n_loc_prodbas, n_states, n_high_state,n_spin)

!     output
      real*8  :: rpa_c_energy_os
      real*8  :: rpa_c_energy_ss
      real*8  :: rpa_c_energy_total
      real*8  :: rpa_c_energy_scs

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
      real*8  rpa_c_integrand_os
      real*8  rpa_c_integrand_ss
      real*8  rpa_c_integrand_total
      real*8  rpa_c_integrand_along_ac_path(rpa_along_ac_path_grid)
      real*8  tmp_renormalized_eg(10)
      !real*8  c_osrpa(30)
      real*8  c_osrpa_integrand(c_osrpa_order)

!    local timing
      real*8  temp_time_rpa
      real*8  temp_clock_time_rpa

!     auxiliary matrices for Level 3 Blas matrix multiplications
!     n_first : the first orbital which is NOT fully occupied

      integer :: n_first(n_spin)

      real*8, dimension(:,:,:), allocatable :: polar_freq

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
      integer :: i_order
      integer :: i_basis_1, i_basis_2

!     Igor, RPA potentials along AC path
      !real*8  :: interval_ac_path      ! interval along AC path
      !if (rpa_along_ac_path_grid .gt. 0) then
      !    interval_ac_path = real(rpa_along_ac_path_grid)**(-1)
      !    allocate(rpa_potentials_along_ac_path(rpa_along_ac_path_grid),stat=i_index)
      !    call check_allocation(i_index, 'rpa_potentials_along_ac_path')
      !    rpa_potentials_along_ac_path(:) = 0.0d0
      !endif
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

      allocate(polar_freq(n_basbas, n_loc_prodbas, n_spin),stat=i_index)
      call check_allocation(i_index, 'polar_freq                    ')


!     Determine the special radius of standard non-interacting response
!     matrix by Kohn-Sham SCF solutions
      tmp_renormalized_eg    = 1.0d0
      tmp_renormalized_eg(1) = 0.0d0
      call  evaluate_ospt2_polarizability_freq_0 &
             ( n_low_state, n_homo, n_first, n_high_state, &
               occ_numbers, &
               0.0d0, &
               KS_eigenvalue, ovlp_3KS, polar_freq, &
               tmp_renormalized_eg &
             )
      call evaluate_special_radius_x0(polar_freq,special_radius)
      do i_spin = 1, n_spin
        renormalized_eg(6+i_spin) = KS_eigenvalue(n_homo(i_spin),i_spin)&
                                   -KS_eigenvalue(n_lumo(i_spin),i_spin)
      end do
      do i_spin = 1, 2
        renormalized_eg(4+i_spin) = special_radius(i_spin)
        ! Renormalization based on HOMO-LUMO gaps and X0 only
        !renormalized_eg(6+i_spin) = &
        ! renormalized_eg(1)*&
        ! renormalized_eg(4+i_spin)**2.d0*&
        !   (renormalized_eg(2)+renormalized_eg(4+i_spin))**(-2.d0)*&
        ! exp(-abs(renormalized_eg(3))*&
        !   (KS_eigenvalue(n_first(i_spin),i_spin)-&
        !    KS_eigenvalue(n_homo(i_spin),i_spin)-renormalized_eg(4))**2.d0)
      end do
           !seg = eg*(1.0d0-renorm_eg(1)*&
           !                exp(-renorm_eg(3)**2.0d0*(abs(eg)-renorm_eg(4))**2.0d0)*&
           !                renorm_eg(4+i_spin)**2.d0/(renorm_eg(2)+renorm_eg(4+i_spin)**2.0d0)&
           !          )

      time_rpa_corr = 0.d0
      clock_time_rpa_corr = 0.d0
      time_polar = 0.d0
      clock_time_polar = 0.d0

      rpa_c_energy_os    = 0.d0 
      rpa_c_energy_ss    = 0.d0 
      rpa_c_energy_total = 0.d0 
      c_osrpa            = 0.d0

      ! IGOR: check if we need to check the special radius of the response matrice of spin channels
      sc_check = .true.
      if(.not.allocated(c_osrpa)) then
       allocate(c_osrpa(c_osrpa_order),stat=i_index)
       call check_allocation(i_index, 'c_osrpa        ')
       c_osrpa(:) = 0.0d0
      endif
      
      do i_freq = 1, n_full_freq, 1

        rpa_c_integrand_os = 0.d0 
        rpa_c_integrand_ss = 0.d0 
        c_osrpa_integrand  = 0.d0

        if(myid.eq.0) then
           write(use_unit,'(A,I4,f20.6)') " | i_freq ", i_freq, & 
                   omega_full(i_freq)
        endif

      call get_timestamps(temp_time_rpa, temp_clock_time_rpa )
!    evaluate the polarisability at frequency point i_freq
      call  evaluate_ospt2_polarizability_freq_0 &
             ( n_low_state, n_homo, n_first, n_high_state, &
               occ_numbers, &
               omega_full(i_freq), &
               KS_eigenvalue, ovlp_3KS, polar_freq, &
               renormalized_eg &
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

       !polar_freq(:,:,:) = polar_freq(:,:,:) * 2.d0/dble(n_spin)

! Hubbard correction, here simply a factor of two
! how to improve?

       call evaluate_osrpa_integrand(polar_freq,&
       rpa_c_integrand_os, rpa_c_integrand_ss,&
       rpa_c_integrand_total,c_osrpa_integrand)

       rpa_c_energy_os = rpa_c_energy_os + &
                       rpa_c_integrand_os * womega_full(i_freq)

       rpa_c_energy_ss = rpa_c_energy_ss + &
                       rpa_c_integrand_ss * womega_full(i_freq)

       rpa_c_energy_total = rpa_c_energy_total + &
                       rpa_c_integrand_total * womega_full(i_freq)

       do i_order = 1, c_osrpa_order
         c_osrpa(i_order) = c_osrpa(i_order) + &
                       c_osrpa_integrand(i_order) * womega_full(i_freq)
       end do

!      write(81,'(I6,4f18.8)')i_freq, omega_full(i_freq), &
!                  womega_full(i_freq), rpa_c_integrand, rpa_c_energy

       call get_timestamps(temp_time_rpa, temp_clock_time_rpa )
       time_rpa_corr = time_rpa_corr +  temp_time_rpa - rtime
       clock_time_rpa_corr = clock_time_rpa_corr + &
                             temp_clock_time_rpa - clock_rtime
! end of loop over i_freq
      enddo
      call sync_timing(time_polar)
      call sync_timing(time_rpa_corr)

!      close(112)

      rpa_c_energy_os=rpa_c_energy_os/2.d0/pi
      rpa_c_energy_ss=rpa_c_energy_ss/2.d0/pi
      rpa_c_energy_total=rpa_c_energy_total/2.d0/pi
      rpa_c_energy_scs = scsrpa_coeff(1)*rpa_c_energy_os + &
           scsrpa_coeff(2)*(rpa_c_energy_total-rpa_c_energy_os)
      do i_order = 1, c_osrpa_order
        c_osrpa(i_order) = c_osrpa(i_order)/2.0d0/pi
      end do

! Delta term correction
! This deals with the case when there are frational occupations and 
! the intra-orbtial excitations should be taken into account. The contribution
! is a delta function at zero frequency and a finite contribution when integrated
! out over frequency axis.

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*) &
              "----------------------------------------------------", &
              "-------------------------"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            "  RPA correlation energy     :", rpa_c_energy_total, "Ha,", &
             rpa_c_energy_total*hartree, "eV"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            "  OS-RPA correlation energy  :", rpa_c_energy_os, "Ha,", &
             rpa_c_energy_os*hartree, "eV"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            "  SS-RPA correlation energy  :", rpa_c_energy_ss, "Ha,", &
             rpa_c_energy_ss*hartree, "eV"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            "  SCS-RPA correlation energy :", rpa_c_energy_scs, "Ha,", &
             rpa_c_energy_scs*hartree, "eV"
        write(use_unit,'(2X,A)') &
            "Now print out os-RPA contribution in higher orders:"
        do i_order = 1, c_osrpa_order
          write(use_unit,'(2X,I6,A,f19.8)') &
               i_order," = ", c_osrpa(i_order)
        end do
        write(use_unit,*)
      endif

      if (allocated (polar_freq)) then
        deallocate (polar_freq)
      endif

      return

      end subroutine evaluate_qpe_osrpa_energy
!---------------------------------------------------------------------
!******
!****s* FHI-aims/evaluate_qpe_osrpa_energy
!  NAME
!   evaluate_qpe_osrpa_energy
!  SYNOPSIS

      subroutine evaluate_qpe_osrpa_energy_2 &
           (n_electrons, &
            n_low_state, n_high_state, &
            n_homo, n_lumo, &
            occ_numbers, n_full_freq, &
            omega_full, womega_full, &
            chemical_potential, &
            KS_eigenvalue, KS_eigenvector, &
            ovlp_3KS, &
            rpa_c_energy_os, rpa_c_energy_ss, & 
            rpa_c_energy_total, rpa_c_energy_scs &
           )

!  PURPOSE
!  Subroutine evaluate_qpe_osrpa_energy evaluates the correlation
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
      use evaluate_ospt2_polarizability_freq, only: evaluate_ospt2_polarizability_freq_2
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
      real*8  :: rpa_c_energy_os
      real*8  :: rpa_c_energy_ss
      real*8  :: rpa_c_energy_total
      real*8  :: rpa_c_energy_scs

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
      real*8  rpa_c_integrand_os
      real*8  rpa_c_integrand_ss
      real*8  rpa_c_integrand_total
      real*8  rpa_c_integrand_along_ac_path(rpa_along_ac_path_grid)
      real*8  tmp_renormalized_eg(10)
      real*8  c_osrpa_integrand(c_osrpa_order)

!    local timing
      real*8  temp_time_rpa
      real*8  temp_clock_time_rpa

!     auxiliary matrices for Level 3 Blas matrix multiplications
!     n_first : the first orbital which is NOT fully occupied

      integer :: n_first(n_spin)

      real*8, dimension(:,:,:), allocatable :: polar_freq

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
       write(use_unit,'(2X, A,A,6I5)') &
                    "HOMO and first non-fully-occupied", &
                    " orbitals:", n_homo(:), n_first(:)
       write(use_unit,*)
      endif

!      open(112, file='rpa_c_energy.dat')

      allocate(polar_freq(max_row_2d, max_col_2d, n_spin),stat=i_index)
      call check_allocation(i_index, 'polar_freq                    ')

!     Determine the special radius of standard non-interacting response
!     matrix by Kohn-Sham SCF solutions
      tmp_renormalized_eg    = 1.0d0
      tmp_renormalized_eg(1) = 0.0d0

      call  evaluate_ospt2_polarizability_freq_2 &
             ( n_low_state, n_homo, n_first, n_high_state, &
               occ_numbers, &
               0.0d0, &
               KS_eigenvalue, ovlp_3KS, polar_freq, &
               tmp_renormalized_eg &
             )
      call evaluate_special_radius_x2(polar_freq,special_radius)
      do i_spin = 1, n_spin
        renormalized_eg(6+i_spin) = KS_eigenvalue(n_homo(i_spin),i_spin)&
                                   -KS_eigenvalue(n_lumo(i_spin),i_spin)
      end do
      do i_spin = 1, 2
        renormalized_eg(4+i_spin) = special_radius(i_spin)
        ! Renormalization based on HOMO-LUMO gaps and X0 only
        !renormalized_eg(6+i_spin) = &
        ! renormalized_eg(1)*&
        ! renormalized_eg(4+i_spin)**2.d0*&
        !   (renormalized_eg(2)+renormalized_eg(4+i_spin))**(-2.d0)*&
        ! exp(-abs(renormalized_eg(3))*&
        !   (KS_eigenvalue(n_first(i_spin),i_spin)-&
        !    KS_eigenvalue(n_homo(i_spin),i_spin)-renormalized_eg(4))**2.d0)
      end do
           !seg = eg*(1.0d0-renorm_eg(1)*&
           !                exp(-renorm_eg(3)**2.0d0*(abs(eg)-renorm_eg(4))**2.0d0)*&
           !                renorm_eg(4+i_spin)**2.d0/(renorm_eg(2)+renorm_eg(4+i_spin)**2.0d0)&
           !          )

      time_rpa_corr = 0.d0
      clock_time_rpa_corr = 0.d0
      time_polar = 0.d0
      clock_time_polar = 0.d0

      rpa_c_energy_os    = 0.d0 
      rpa_c_energy_ss    = 0.d0 
      rpa_c_energy_total = 0.d0 
      c_osrpa            = 0.d0

      ! IGOR: check if we need to check the special radius of the response matrice of spin channels
      sc_check = .true.
      if(.not.allocated(c_osrpa)) then
       allocate(c_osrpa(c_osrpa_order),stat=i_index)
       call check_allocation(i_index, 'c_osrpa        ')
       c_osrpa(:) = 0.0d0
      endif

      do i_freq = 1, n_full_freq, 1

        rpa_c_integrand = 0.d0 
        rpa_c_integrand_os = 0.d0 
        rpa_c_integrand_ss = 0.d0 
        c_osrpa_integrand  = 0.d0

        if(myid.eq.0) then
           write(use_unit,'(A,I4,f20.6)') " | i_freq ", i_freq, & 
                   omega_full(i_freq)
        endif

      call get_timestamps(temp_time_rpa, temp_clock_time_rpa )
!    evaluate the polarisability at frequency point i_freq
      call  evaluate_ospt2_polarizability_freq_2 &
             ( n_low_state, n_homo, n_first, n_high_state, &
               occ_numbers, &
               omega_full(i_freq), &
               KS_eigenvalue, ovlp_3KS, polar_freq, renormalized_eg &
             )

       call get_timestamps(rtime, clock_rtime)
       time_polar = time_polar + rtime - temp_time_rpa
       clock_time_polar = clock_time_polar + clock_rtime - temp_clock_time_rpa

!       polar_freq(:,:,:) = polar_freq(:,:,:) * 2.d0/dble(n_spin)

!! Hubaard correction, here simply a factor of two
!! how to improve?
       call evaluate_osrpa_integrand_2(polar_freq,&
                          rpa_c_integrand_os, rpa_c_integrand_ss,&
                          rpa_c_integrand_total,c_osrpa_integrand &
                          )

       rpa_c_energy_os = rpa_c_energy_os + &
                       rpa_c_integrand_os * womega_full(i_freq)

       rpa_c_energy_ss = rpa_c_energy_ss + &
                       rpa_c_integrand_ss * womega_full(i_freq)

       rpa_c_energy_total = rpa_c_energy_total + &
                       rpa_c_integrand_total * womega_full(i_freq)

       call get_timestamps(temp_time_rpa, temp_clock_time_rpa )
       time_rpa_corr = time_rpa_corr +  temp_time_rpa - rtime
       clock_time_rpa_corr = clock_time_rpa_corr + temp_clock_time_rpa - clock_rtime
! end of loop over i_freq
      enddo
      call sync_timing(time_polar)
      call sync_timing(time_rpa_corr)

!      close(112)

!      rpa_c_energy=rpa_c_energy/2.d0/pi
      rpa_c_energy_os=rpa_c_energy_os/2.d0/pi
      rpa_c_energy_ss=rpa_c_energy_ss/2.d0/pi
      rpa_c_energy_total=rpa_c_energy_total/2.d0/pi
      rpa_c_energy_scs = scsrpa_coeff(1)*rpa_c_energy_os + &
           scsrpa_coeff(2)*(rpa_c_energy_total-rpa_c_energy_os)
      !do i_order = 1, c_osrpa_order
      !  c_osrpa(i_order) = c_osrpa(i_order)/2.0d0/pi
      !end do

! Delta term correction
! This deals with the case when there are frational occupations and 
! the intra-orbtial excitations should be taken into account. The contribution
! is a delta function at zero frequency and a finite contribution when integrated
! out over frequency axis.


      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*) &
              "----------------------------------------------------", &
              "-------------------------"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            "  RPA correlation energy     :", rpa_c_energy_total, "Ha,", &
             rpa_c_energy_total*hartree, "eV"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            "  OS-RPA correlation energy  :", rpa_c_energy_os, "Ha,", &
             rpa_c_energy_os*hartree, "eV"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            "  SS-RPA correlation energy  :", rpa_c_energy_ss, "Ha,", &
             rpa_c_energy_ss*hartree, "eV"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            "  SCS-RPA correlation energy :", rpa_c_energy_scs, "Ha,", &
             rpa_c_energy_scs*hartree, "eV"
        write(use_unit,*)
      endif

      if (allocated (polar_freq)) then
        deallocate (polar_freq)
      endif

      return

      end subroutine evaluate_qpe_osrpa_energy_2
!---------------------------------------------------------------------
!******
