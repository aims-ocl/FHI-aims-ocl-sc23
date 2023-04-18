      subroutine get_spectrum_dmft (anacon_type, green_fn_par, n_max_par, &
             omega, nomega, filename, n_matrix ,&
             spectrum,aux_omegamax, aux_omega, aux_womega, aux_nomega,&
             DMFT_dos_LDA, DMFT_dos,new_chemical_potential)
!get the spectrum from the (analytically continued) diagonalized Green's function


      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use physics
      use localorb_io, only: use_unit

!ARGUMENTS
      implicit none
      integer n_matrix
      integer n_max_par
      integer nomega,i_e
      integer anacon_type
      complex*16 green_fn_par (n_max_par, n_matrix, n_matrix)
!      complex*16 aux_green_KS (n_max_par, n_states)
!      complex*16 inv_overlap_matrix(n_basis,n_basis)
      real*8 omega(nomega)
      !real*8,  dimension (-aux_nomega:aux_nomega) ::  spectrum
!      real*8,  dimension (:), allocatable ::  spectrum_partial
      character*15 filename
      real*8 new_chem_pot,en , de

!internal
      complex*16 diagonal_green (n_matrix)
      complex*16 aux_green (n_matrix, n_matrix)
      complex*16 denominator
      complex*16 numerator
!      real*8     w
      complex*16 w 
      complex*16 xdata (n_max_par)
      complex*16 gtmp
      real*8 n_particles
!the local grid
      integer aux_nomega
      real*8  aux_omegamax
      real*8 , dimension (aux_nomega) :: aux_omega
      real*8 , dimension (aux_nomega) :: aux_womega
      real*8 , dimension (aux_nomega) :: n_particles_freq
      real*8,  dimension (-aux_nomega:aux_nomega) ::  spectrum
      real*8, dimension(dos_n_en_points) :: DMFT_dos_LDA
      real*8, dimension(dos_n_en_points) :: DMFT_dos

!counter
      integer :: i_basis
      integer :: j_basis
      integer :: i_freq
      integer :: i_par
      integer :: i_state
      integer :: sign_freq
      integer :: n_step
      integer :: i_dat
      integer i_index

      real*8 e_fermi
      real*8 new_chemical_potential
      integer n_homo
 
      e_fermi = chemical_potential

      if (myid.eq.0)then
        open (44,file=filename)
      endif
        open (45,file='N_electrons_freq.dat')

       do sign_freq = -1, 1, 2
        do i_freq = 1, aux_nomega, 1

         w = sign_freq*aux_omega(i_freq) 

        if (myid.eq.0)then

        if ( spectrum (sign_freq*i_freq) .lt. 0) then
             !write(44, *)(real(w)+e_fermi)*hartree,&
             !          - spectrum (sign_freq*i_freq)
             write(44, *)(real(w))*hartree,&
                       - spectrum (sign_freq*i_freq)
!             write(45, *)(real(w)+e_fermi)*hartree,&
!                        spectrum_partial (sign_freq*i_freq)
         else
            !write(44, *)(real(w)+e_fermi)*hartree,&
            !            spectrum (sign_freq*i_freq)
            write(44, *)(real(w))*hartree,&
                        spectrum (sign_freq*i_freq)
         endif
         endif

!       enddo  ! i_freq
      enddo !sign_freq
      enddo !sign_freq

      !endif !anacon_type

      de= (dos_high_energy - dos_low_energy)/dble(dos_n_en_points-1)

     if(myid.eq.0)then
     ! arrange output for DOS: different for spin polarized and normal calculations
     write(use_unit,'(2X,A)') '| writing DOS (shifted by electron chemical potential)&
                          & to file DMFT_dos_total.dat'
     write(use_unit,'(2X,A)') '| writing DOS (raw data) to file DMFT_dos_total_raw.dat'
     open(88, file='DMFT_DOS_total.dat')
     open(89, file='DMFT_DOS_total_raw.dat')
        write(88,'(2A)') '# total density of states output by FHI-aims '
        write(88,'(A,F15.6,A)') '# The energy reference for this output is &
                          &the chemical potential, mu = ', &
             chemical_potential*hartree, ' eV'
        write(88,'(A)') '#    Energy (eV)      DOS '
        write(89,'(2A)') '# total density of states output by FHI-aims '
        write(89,'(A,F15.6,A)') '# The energy reference for this output is the&
                          & vacuum level'
        write(89,'(A)') '#    Energy (eV)      DOS '
        do i_e = 1, dos_n_en_points ,1
           en = dos_low_energy + dble(i_e-1)*de
           write(88,'(3F16.8)') en-chemical_potential*hartree, DMFT_dos(i_e)
           write(89,'(3F16.8)') en, DMFT_dos(i_e)
        enddo

     close(88)
     close(89)

     end if


      if(.false.)then




         n_particles = 0.d0
       do sign_freq = -1, 1, 2
         do i_freq= 1, aux_nomega, 1
          if( sign_freq*aux_omega(i_freq) .le. 0 )then
            n_particles = n_particles + spectrum(sign_freq*i_freq)*aux_womega(i_freq)
          endif

         enddo
        enddo


      endif


      return

      end subroutine get_spectrum_dmft
