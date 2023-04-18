!****s* FHI-aims/evaluate_C6_coeff
!  NAME
!   evaluate_C6_coeff
!  SYNOPSIS

      subroutine evaluate_C6_coeff &
           (n_electrons, &
            n_low_state, n_high_state, &
            n_homo, n_lumo, &
            occ_numbers, n_full_freq, &
            omega_full, womega_full, &
            chemical_potential, &
            KS_eigenvalue, KS_eigenvector, &
            ovlp_3KS &
           )

!  PURPOSE
!  Subroutine evaluate_C6_coeffcient at the  MP2 and RPA level
!  MP2 C6 coefficent is given by integrating over the bare polarizability,
!  while RPA C6 is given by the RPA polarizability

! USES
      use runtime_choices, only: out_polarisability
      use dimensions
      use prodbas
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use evaluate_polarisability_freq, only: evaluate_polarisability_freq_0
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

      real*8  RPA_C6_coeff
      real*8  MP2_C6_coeff
      real*8  dipole_polar
      real*8  dipole_bare_polar
!      real*8  x_orient(3)
!      real*8  xsq
!      real*8  orient_matr(3,3)
!      real*8  product_polar_matr(3,3)
      real*8  coord_of_center(3)

!     auxiliary matrices for Level 3 Blas matrix multiplications
!     n_first : the first orbital which is NOT fully occupied

      integer :: n_first(n_spin)
      integer :: n_occ
      integer :: n_unocc

      real*8, dimension(:,:), allocatable :: polar_freq
      real*8, dimension(:,:), allocatable :: atom_coord_wrt_center
!      real*8, dimension(:,:), allocatable :: dipole_mom_prodbas
      real*8, dimension(:,:,:,:), allocatable :: dipole_mom
!      real*8, dimension(:,:), allocatable :: dipole_polar_matr

!     timing

!     parameters of the fitting tails
!       real*8  s1, s2, omega_1, omega_2
!       real*8  alpha, beta, a, b

      character*50  filename

!     counters

      integer :: i_state
      integer :: i_state_1
      integer :: i_freq
      integer :: i_spin
      integer :: i_index
      integer :: i_dim, i_dim_1
      integer :: i_basbas, i_basbas_1
      integer :: i_order
real*8 ttt0
ttt0 = mpi_wtime()
!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"-------------------------------------------------"
        write(use_unit,*)
        write(use_unit,'(2X,A)') &
              "Start to calculate the C6 coefficient  ... "
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


      if(myid.eq.0) then
       write(use_unit,*)
       write(use_unit,'(2X, A,A,4I5)') "HOMO and first non-fully-occupied", &
                    " orbitals:", n_homo(:), n_first(:)
       write(use_unit,*)
      endif

      n_occ = max(n_homo(1),n_homo(n_spin))
      n_unocc =  n_states-min(n_first(1),n_first(n_spin))+1

      allocate(atom_coord_wrt_center(3,n_atoms), stat=i_index)
      call check_allocation(i_index, 'atom_coord_wrt_center         ')

!      allocate(dipole_mom_prodbas(3,n_loc_prodbas), stat=i_index)
!      call check_allocation(i_index, 'dipole_mom_prodbas            ')

!      allocate(dipole_polar_matr(3,3), stat=i_index)
!      call check_allocation(i_index, 'dipole_polar_matr             ')

!      call determine_dipole_moment_of_prodbas & 
!                (atom_coord_wrt_center,dipole_mom_prodbas)

      call determine_center_of_molecule & 
                (coord_of_center, atom_coord_wrt_center)

      allocate(dipole_mom(n_occ,n_unocc,n_spin,3), stat=i_index)
      call check_allocation(i_index, 'dipole_mom                    ')

      call integrate_dipmom_pairstates &
           (n_occ,n_unocc, n_first, coord_of_center, l_shell_max, &
            KS_eigenvector, dipole_mom)

!  determine the orientation matrix given by
!    t_ij (R_vec) = (3R_i*R_j - delta_ij |R_vec|^2)/|R_vec|^2 
!  for the time being let's assue R_i = R_j =0, R_k = |R_vec|,
!  the vector is pointing in the z direction

!      x_orient(1) = 0.d0
!      x_orient(2) = 0.d0
!      x_orient(3) = 1.d0
!      xsq = x_orient(1)*x_orient(1)+ &
!            x_orient(2)*x_orient(2)+ &
!            x_orient(3)*x_orient(3)
!
!      do i_dim = 1, 3, 1
!        do i_dim_1 = 1, 3, 1
!          if(i_dim .eq. i_dim_1) then
!            orient_matr(i_dim_1, i_dim) = 3*x_orient(i_dim_1)*x_orient(i_dim)/xsq -1
!          else
!            orient_matr(i_dim_1, i_dim) = 3*x_orient(i_dim_1)*x_orient(i_dim)/xsq 
!          endif
!        enddo
!        write(use_unit,'(2X, A, 3f16.8)') "oritention matrix: ",  &
!                (orient_matr(i_dim_1, i_dim), i_dim_1 = 1, 3)
!      enddo


      allocate(polar_freq(n_basbas, n_loc_prodbas),stat=i_index)
      call check_allocation(i_index, 'polar_freq                    ')

      if(out_polarisability) then
        open(112,file='mp2_polarizability.dat')
        open(113,file='rpa_polarizability.dat')
      endif

      RPA_C6_coeff = 0.d0
      MP2_C6_coeff = 0.d0
!if(myid==0) print *,'Start',mpi_wtime()-ttt0
      do i_freq = 1, n_full_freq, 1

!    evaluate the polarisability at frequency point i_freq
        call  evaluate_polarisability_freq_0 &
             ( n_low_state, n_homo, n_first, n_high_state, &
               occ_numbers, &
               omega_full(i_freq), &
               KS_eigenvalue, ovlp_3KS, polar_freq &
             )
!if(myid==0) print *,'evaluate_polarisability_freq',mpi_wtime()-ttt0

        polar_freq(:,:) = polar_freq(:,:) * 2.d0/dble(n_spin)


       call evaluate_dipole_polarisability &
           ( n_occ,n_unocc,n_homo,n_first, &
             n_high_state,occ_numbers, omega_full(i_freq), &
             KS_eigenvalue, ovlp_3KS, polar_freq, &
             dipole_mom, dipole_bare_polar,dipole_polar &
           )
!if(myid==0) print *,'evaluate_dipole_polarisability',mpi_wtime()-ttt0

!    get the product polarizability matrix
!    Here we assume the two fragments are the same, and this needs to be
!    extended to the more general case
!        do i_dim = 1, 3, 1
!            do i_dim_1 = 1, 3, 1
!               product_polar_matr(i_dim_1, i_dim) = &
!                 dipole_polar_matr(i_dim_1, i_dim) * dipole_polar_matr(i_dim_1, i_dim)
!            enddo
!        enddo
!
!        call dgemm( 'N', 'N', 3, 3, 3, 1.d0, orient_matr, 3, &
!                product_polar_matr, 3, 0.d0, &
!                dipole_polar_matr, 3 )   
!
!        call dgemm( 'N', 'N', 3, 3, 3, 1.d0, dipole_polar_matr, 3, & 
!                orient_matr, 3,  0.d0, &
!                product_polar_matr, 3 )   

!        dipole_bare_polar=0.d0
!        do i_dim = 1, 3, 1
!          do i_dim_1 = 1, 3, 1
!             dipole_bare_polar = dipole_bare_polar &
!                + product_polar_matr(i_dim_1,i_dim)
!          enddo
!        enddo

!        if(myid.eq.0) then
!           write(use_unit,'(2X, I4, 4f16.6)') i_freq, omega_full(i_freq), &
!                   womega_full(i_freq), dipole_bare_polar, dipole_polar
!        endif
        if(myid.eq.0 .and. out_polarisability) then
           write(112,'(2X, I4, 4f16.6)') i_freq, omega_full(i_freq), &
                   womega_full(i_freq), dipole_bare_polar
           write(113,'(2X, I4, 4f16.6)') i_freq, omega_full(i_freq), &
                   womega_full(i_freq), dipole_polar
        endif

        MP2_C6_coeff = MP2_C6_coeff +  dipole_bare_polar * dipole_bare_polar* &
                  womega_full(i_freq)

        RPA_C6_coeff = RPA_C6_coeff + dipole_polar *  &
                    dipole_polar*womega_full(i_freq)

      enddo


      RPA_C6_coeff = RPA_C6_coeff*3.d0/pi
      MP2_C6_coeff = MP2_C6_coeff*3.d0/pi

      if(out_polarisability) then
        close(112)
        close(113)
      endif

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"----------------------------------------------------", &
                  "-------------------------"
                 
        write(use_unit,'(2X,A,2X,f12.3,2X,A,f12.3,2X,A)') &
            "  RPA C6 coefficient :", RPA_C6_coeff, "Ha*Bohr^6,", &
               RPA_C6_coeff*hartree*bohr**6, "eV*Ang^6"

        write(use_unit,*)
        write(use_unit,'(2X,A,2X,f12.3,2X,A,f12.3,2X,A)') &
            "  MP2 C6 coefficient :", MP2_C6_coeff, "Ha*Bohr^6,", &
               MP2_C6_coeff*hartree*bohr**6, "eV*Ang^6"
        write(use_unit,*)"----------------------------------------------------", &
                  "-------------------------"
        write(use_unit,*)
      endif

      if (allocated (polar_freq)) then
        deallocate (polar_freq)
      endif
      if (allocated (atom_coord_wrt_center)) then
        deallocate (atom_coord_wrt_center)
      endif
!      if (allocated (dipole_polar_matr)) then
!        deallocate (dipole_polar_matr)
!      endif
      if (allocated (dipole_mom)) then
        deallocate (dipole_mom)
      endif

      return

      end subroutine evaluate_C6_coeff
!---------------------------------------------------------------------
!******

      subroutine evaluate_C6_coeff_2 &
           (n_electrons, &
            n_low_state, n_high_state, &
            n_homo, n_lumo, &
            occ_numbers, n_full_freq, &
            omega_full, womega_full, &
            chemical_potential, &
            KS_eigenvalue, KS_eigenvector, &
            ovlp_3KS &
           )

!  PURPOSE
!  Subroutine evaluate_C6_coeffcient at the  MP2 and RPA level
!  MP2 C6 coefficent is given by integrating over the bare polarizability,
!  while RPA C6 is given by the RPA polarizability

! USES
      use runtime_choices, only: out_polarisability
      use dimensions
      use prodbas
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
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
      real*8  :: ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)

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

      real*8  RPA_C6_coeff
      real*8  MP2_C6_coeff
      real*8  dipole_polar
      real*8  dipole_bare_polar
      real*8  coord_of_center(3)

!     auxiliary matrices for Level 3 Blas matrix multiplications
!     n_first : the first orbital which is NOT fully occupied

      integer :: n_first(n_spin)
      integer :: n_occ
      integer :: n_unocc

      real*8, dimension(:,:), allocatable :: polar_freq
      real*8, dimension(:,:), allocatable :: atom_coord_wrt_center
      real*8, dimension(:,:,:,:), allocatable :: dipole_mom

      character*50  filename

!     counters

      integer :: i_state
      integer :: i_state_1
      integer :: i_freq
      integer :: i_spin
      integer :: i_index
      integer :: i_dim, i_dim_1
      integer :: i_basbas, i_basbas_1
      integer :: i_order
real*8 ttt0
ttt0 = mpi_wtime()
!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"-------------------------------------------------"
        write(use_unit,*)
        write(use_unit,'(2X,A)') &
              "Start to calculate the C6 coefficient  ... "
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


      if(myid.eq.0) then
       write(use_unit,*)
       write(use_unit,'(2X, A,A,4I5)') "HOMO and first non-fully-occupied", &
                    " orbitals:", n_homo(:), n_first(:)
       write(use_unit,*)
      endif

      n_occ = max(n_homo(1),n_homo(n_spin))
      n_unocc =  n_states-min(n_first(1),n_first(n_spin))+1

      allocate(atom_coord_wrt_center(3,n_atoms), stat=i_index)
      call check_allocation(i_index, 'atom_coord_wrt_center         ')

      call determine_center_of_molecule & 
                (coord_of_center, atom_coord_wrt_center)

      allocate(dipole_mom(n_occ,n_unocc,n_spin,3), stat=i_index)
      call check_allocation(i_index, 'dipole_mom                    ')

      call integrate_dipmom_pairstates &
           (n_occ,n_unocc, n_first, coord_of_center, l_shell_max, &
            KS_eigenvector, dipole_mom)

      allocate(polar_freq(max_row_2d, max_col_2d),stat=i_index)
      call check_allocation(i_index, 'polar_freq                    ')

      if(out_polarisability) then
        open(112,file='mp2_polarizability.dat')
        open(113,file='rpa_polarizability.dat')
      endif

      RPA_C6_coeff = 0.d0
      MP2_C6_coeff = 0.d0
!if(myid==0) print *,'Start',mpi_wtime()-ttt0
      do i_freq = 1, n_full_freq, 1

!    evaluate the polarisability at frequency point i_freq
        call  evaluate_polarisability_freq_2 &
             ( n_low_state, n_homo, n_first, &
               occ_numbers, &
               omega_full(i_freq), &
               KS_eigenvalue, ovlp_3KS, polar_freq &
             )
!if(myid==0) print *,'evaluate_polarisability_freq',mpi_wtime()-ttt0

        polar_freq(:,:) = polar_freq(:,:) * 2.d0/dble(n_spin)

       call evaluate_dipole_polarisability_2 &
           ( n_occ,n_unocc,n_homo,n_first, &
             n_high_state,occ_numbers, omega_full(i_freq), &
             KS_eigenvalue, ovlp_3KS, polar_freq, &
             dipole_mom, dipole_bare_polar,dipole_polar &
           )
!if(myid==0) print *,'evaluate_dipole_polarisability',mpi_wtime()-ttt0

        if(myid.eq.0 .and. out_polarisability) then
           write(112,'(2X, I4, 4f16.6)') i_freq, omega_full(i_freq), &
                   womega_full(i_freq), dipole_bare_polar
           write(113,'(2X, I4, 4f16.6)') i_freq, omega_full(i_freq), &
                   womega_full(i_freq), dipole_polar
        endif

        MP2_C6_coeff = MP2_C6_coeff +  dipole_bare_polar * dipole_bare_polar* &
                  womega_full(i_freq)

        RPA_C6_coeff = RPA_C6_coeff + dipole_polar *  &
                    dipole_polar*womega_full(i_freq)

      enddo


      RPA_C6_coeff = RPA_C6_coeff*3.d0/pi
      MP2_C6_coeff = MP2_C6_coeff*3.d0/pi

      if(out_polarisability) then
        close(112)
        close(113)
      endif

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"----------------------------------------------------", &
                  "-------------------------"
                 
        write(use_unit,'(2X,A,2X,f12.3,2X,A,f12.3,2X,A)') &
            "  RPA C6 coefficient :", RPA_C6_coeff, "Ha*Bohr^6,", &
               RPA_C6_coeff*hartree*bohr**6, "eV*Ang^6"

        write(use_unit,*)
        write(use_unit,'(2X,A,2X,f12.3,2X,A,f12.3,2X,A)') &
            "  MP2 C6 coefficient :", MP2_C6_coeff, "Ha*Bohr^6,", &
               MP2_C6_coeff*hartree*bohr**6, "eV*Ang^6"
        write(use_unit,*)"----------------------------------------------------", &
                  "-------------------------"
        write(use_unit,*)
      endif

      if (allocated (polar_freq)) then
        deallocate (polar_freq)
      endif
      if (allocated (atom_coord_wrt_center)) then
        deallocate (atom_coord_wrt_center)
      endif
!      if (allocated (dipole_polar_matr)) then
!        deallocate (dipole_polar_matr)
!      endif
      if (allocated (dipole_mom)) then
        deallocate (dipole_mom)
      endif

      return

      end subroutine evaluate_C6_coeff_2
!---------------------------------------------------------------------
!******
