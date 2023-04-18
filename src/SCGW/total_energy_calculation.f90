!total_energy_calculation
! Evaluate the Galitskii-Migdal total energy

      subroutine total_energy_calculation_v1 (&
                   green_fn_time,             &
                   green_fn_freq,             &
                   exchange_self_energy,      &
                   exchange_self_energy0,     &
                   hartree_pot,               &
                   hartree_pot0,              &
                   self_energy_freq,          &
                   xc_matr,                   &
                   rpa_c_energy               & 
                   )

      use dimensions
      use runtime_choices
      use species_data
      use physics
      use prodbas
      use hartree_fock
      use constants
      use mpi_tasks
      use synchronize_mpi
      use poles_fit
      use scgw_grid 
      use localorb_io, only: use_unit


      implicit none
      !ARGUMENTS
      complex*16  ::  green_fn_freq         (n_basis,n_basis,nomega,n_spin)
      complex*16  ::  self_energy_freq      (n_basis,n_basis,nomega,n_spin)
      real*8      ::  green_fn_time         (n_basis,n_basis,n_spin)
      real*8      ::  exchange_self_energy  (n_basis,n_basis,n_spin)
      real*8      ::  exchange_self_energy0 (n_basis,n_basis,n_spin)
      real*8      ::  xc_matr               (n_basis,n_basis,n_spin)
      real*8      ::  hartree_pot           (n_basis,n_basis)
      real*8      ::  hartree_pot0          (n_basis,n_basis)
      real*8      ::  Galitskii_Migdal_tot_energy
      real*8      ::  rpa_c_energy

      !INTERNAL
      real*8 dft_xc_en
      real*8 correlation_tot_energy
      real*8 exchange_energy
      real*8 exchange_energy0
      real*8 hartree_en
      real*8 hartree_en0
      real*8 scgw_kinetic_energy 
      real*8 KS_hamilt_en       

      integer i_spin
      integer i_freq
      integer i_basis, j_basis 

!start work

      call E_kin ( green_fn_time, scgw_kinetic_energy )

      call single_particle_energy_v2 ( &
         green_fn_time,           &
         hartree_pot,             &
         hartree_pot0,            &
         xc_matr ,                &
         exchange_self_energy,    &
         exchange_self_energy0,   &
         exchange_energy,         &
         exchange_energy0,        &
         hartree_en,              &
         hartree_en0,             &
         dft_xc_en,               &
         KS_hamilt_en             &
         )

      !evaluation of the correlation energy
      correlation_tot_energy = 0.d0
      do i_spin = 1, n_spin
        do i_basis = 1, n_basis, 1
          do j_basis  = 1, n_basis , 1
            do i_freq = 2, nomega, 1
              correlation_tot_energy = correlation_tot_energy +&
              self_energy_freq(i_basis, j_basis, i_freq,i_spin) *&
              green_fn_freq (i_basis, j_basis, i_freq,i_spin)* womega1 (i_freq) &
              / (2.d0 * pi)
            enddo
            correlation_tot_energy = correlation_tot_energy +&
              self_energy_freq(i_basis, j_basis,1,i_spin) *&
              green_fn_freq (i_basis, j_basis,1 ,i_spin)* womega1 (1) &
              / (2.d0 * pi)/2.d0 !the zero point should have half the 
                                 !weight for integration of the whole axis
          enddo
        enddo
      enddo
      correlation_tot_energy = correlation_tot_energy*(2.d0/n_spin)

      !sum of all the terms
      if (.not. use_hartree_fock)then
         Galitskii_Migdal_tot_energy = &
           KS_hamilt_en                & 
         - dft_xc_en                   &  
         - hartree_en0 * 2.d0          &
         + en_ion_ion                  &
         + correlation_tot_energy      &
         + hartree_en                  &
         + exchange_energy
      else
         Galitskii_Migdal_tot_energy =                &
           KS_hamilt_en                               &
         - hybrid_coeff * exchange_energy0 * 2.d0     &
         - dft_xc_en                                  &  
         - hartree_en0      * 2.d0                    & 
         + en_ion_ion                                 &
         + correlation_tot_energy                     &
         + hartree_en                                 & 
         + exchange_energy 
      endif

      !print out the total energy
      if (myid.eq.0)then
        write(use_unit,*) " "
        write(use_unit,'(A)') "   --- GW Total Energy Calculation"
        write(use_unit,'(A)') "          |"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | DFT Total Energy                  :"&
                 ,total_energy * hartree, " eV   ", total_energy , " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Galitskii-Migdal Total Energy     :",&
          Galitskii_Migdal_tot_energy* hartree, " eV   ", &
          Galitskii_Migdal_tot_energy, " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | GW Kinetic Energy                 :",&
          scgw_kinetic_energy* hartree, " eV   ", &
          scgw_kinetic_energy , " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Hartree energy from DFT density   :",&
                 (hartree_en0)* hartree, " eV   ",&
                  hartree_en0 ,   " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Hartree energy from GW density    :",&
              hartree_en * hartree, " eV   ", hartree_en ,   " Ha"
        if(.not. use_hartree_fock)then
          write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | DFT xc energy                     :", &
              dft_xc_en* hartree, " eV   ", dft_xc_en , " Ha"
        endif
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Nuclear Energy                    :", &
          en_ion_ion * hartree, " eV   " , en_ion_ion, " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | GW correlation Energy             :", &
            (correlation_tot_energy)* hartree, " eV   ", &
            correlation_tot_energy , " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | RPA correlation Energy            :", &
            (rpa_c_energy)* hartree, " eV   ", &
            rpa_c_energy , " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Initial Exchange Energy           :", &
        exchange_energy0* hartree, " eV   ", exchange_energy0, " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Exchange Energy                   :", &
        exchange_energy* hartree, " eV   ", exchange_energy, " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Sigle Particle Energy             :", &
                      KS_hamilt_en * hartree, " eV   ", &
                      KS_hamilt_en , " Ha"
        write(use_unit,*) " "
      endif

      end subroutine total_energy_calculation_v1

!------------------------------------------------
!total_energy_calculation
! Evaluate the Galitskii-Migdal total energy

      subroutine total_energy_calculation (&
                   green_fn_time, green_fn_freq,  &
               !    omega, nomega, womega,    &
                   exchange_self_energy, self_energy_freq,     &
                   xc_matr, &
                   hartree_pot, &
                   full_hartree_pot,&
                   rpa_c_energy,& 
                   n_iter, exchange_self_energy0)

      use dimensions
      use runtime_choices
      use species_data
      use physics
      use prodbas
      use hartree_fock
      use constants
      use mpi_tasks
      use synchronize_mpi
      use poles_fit
      use scgw_grid 
      use localorb_io, only: use_unit

      implicit none

      !ARGUMENTS
!      integer  nomega
!      real*8  :: omega(nomega)
!      real*8  :: womega(nomega)
      complex*16 green_fn_freq(n_basis,n_basis, nomega, n_spin)
      real*8  :: green_fn_time (n_basis,n_basis,n_spin)
      real*8  :: exchange_self_energy (n_basis,n_basis,n_spin)
      real*8  :: exchange_self_energy0 (n_basis,n_basis,n_spin)
      complex*16  self_energy_freq (n_basis,n_basis, nomega, n_spin)
      real*8  :: xc_matr (n_basis,n_basis,n_spin)
      real*8  :: hartree_pot (n_basis,n_basis,n_spin)
!      real*8  :: hartree_dft (n_basis,n_basis,n_spin)
      real*8  :: full_hartree_pot (n_basis,n_basis,n_spin)
      real*8  :: Galitskii_Migdal_tot_energy
!      real*8 ovlp_NAO_KS (n_states,n_basis,n_spin)
!      real*8 inv_overlap_matrix (n_basis,n_basis)
      real*8 rpa_c_energy
      integer n_iter

      !INTERNAL
      real*8 dft_xc_en
      real*8 correlation_tot_energy
      real*8 delta_hartree, exchange_energy, exchange_energy0
      real*8 hartree_en
!      real*8 hartree_en_ensemble
      real*8 hartree_en_GW
      real*8 sum_eigenvalues
      real*8 Galitskii_Migdal_tot_energy2
      real*8 scgw_kinetic_energy

      integer i_spin
      integer i_freq
      integer i_basis, j_basis

!start work

      call E_kin ( green_fn_time, scgw_kinetic_energy )

      call single_particle_energy ( &
         n_homo(1),&
         green_fn_time(:,:,:), &
         !inv_overlap_matrix  ,  &
         hartree_pot, &
         dft_xc_en, & 
         xc_matr,exchange_self_energy ,&
         exchange_self_energy0 ,&
         delta_hartree, exchange_energy, exchange_energy0, &
         dft_xc_en, full_hartree_pot, &
         hartree_en, hartree_en_GW )!, hartree_dft, hartree_en_ensemble)

       exchange_energy = exchange_energy/2.d0
       exchange_energy0 = exchange_energy0/2.d0

       !evaluation of the correlation energy
       correlation_tot_energy = 0.d0
       do i_spin = 1, n_spin
         do i_basis = 1, n_basis, 1
           do j_basis  = 1, n_basis , 1
             do i_freq = 2, nomega, 1
               correlation_tot_energy = correlation_tot_energy +&
               self_energy_freq(i_basis, j_basis, i_freq,i_spin) *&
               green_fn_freq (i_basis, j_basis, i_freq,i_spin)* womega (i_freq) &
               / (2.d0 * pi)
             enddo
             correlation_tot_energy = correlation_tot_energy +&
               self_energy_freq(i_basis, j_basis,1,i_spin) *&
               green_fn_freq (i_basis, j_basis,1 ,i_spin)* womega (1) &
               / (2.d0 * pi)/2.d0 !the zero point should have half the 
                                  !weight for integration of the whole axis
           enddo
         enddo
       enddo
       correlation_tot_energy = correlation_tot_energy*(2.d0/n_spin)

      !sum of all the terms
      if (.not. use_hartree_fock)then
         Galitskii_Migdal_tot_energy = &
           dft_xc_en             & 
         - dft_xc_en                   &  
         - hartree_en *2               &
         + correlation_tot_energy      &
         + en_ion_ion                  &
         + hartree_en_GW               &
         + (hartree_en_GW-hartree_en)  & 
         + exchange_energy
      else
         Galitskii_Migdal_tot_energy = &
           dft_xc_en             &
         - hartree_en                  & 
         - exchange_energy0           &
         + correlation_tot_energy      &
         + en_ion_ion                  &
         + 2.0*(hartree_en_GW-hartree_en)  &
         + (exchange_energy - exchange_energy0) 
      endif

      !print out the total energy
      if (myid.eq.0)then
        write(use_unit,*) " "
        write(use_unit,'(A)') "   --- GW Total Energy Calculation"
        write(use_unit,'(A)') "          |"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | DFT Total Energy                  :"&
                 ,total_energy * hartree, " eV   ", total_energy , " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Galitskii-Migdal Total Energy     :",&
          Galitskii_Migdal_tot_energy* hartree, " eV   ", &
          Galitskii_Migdal_tot_energy, " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | GW Kinetic Energy                 :",&
          scgw_kinetic_energy* hartree, " eV   ", &
          scgw_kinetic_energy , " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Hartree energy from DFT density   :",&
                 (hartree_en)* hartree, " eV   ",&
                  hartree_en ,   " Ha"
!        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Hartree energy for ensemble calc. :",&
!                 (hartree_en_ensemble)* hartree, " eV   ",&
!                  hartree_en_ensemble ,   " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Hartree energy from GW density    :",&
              hartree_en_GW * hartree, " eV   ", hartree_en_GW ,   " Ha"
        if(.not. use_hartree_fock)then
          write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | DFT xc energy                     :", &
              dft_xc_en* hartree, " eV   ", dft_xc_en , " Ha"
        endif
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Nuclear Energy                    :", &
          en_ion_ion * hartree, " eV   " , en_ion_ion, " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | GW correlation Energy             :", &
            (correlation_tot_energy)* hartree, " eV   ", &
            correlation_tot_energy , " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | RPA correlation Energy            :", &
            (rpa_c_energy)* hartree, " eV   ", &
            rpa_c_energy , " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Exchange Energy                   :", &
        exchange_energy* hartree, " eV   ", exchange_energy, " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Sigle Particle Energy             :", &
                      dft_xc_en* hartree, " eV   ", &
                      dft_xc_en, " Ha"
        write(use_unit,*) " "
      endif

      end subroutine total_energy_calculation
