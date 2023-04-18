!****s* FHI-aims/evaluate_sosex_selfenergy
!  NAME
!   evaluate_sosex_selfenergy
!  SYNOPSIS

      subroutine evaluate_sosex_selfenergy &
       ( n_low_state, n_KS_states, &
         ovlp_3KS, &
         n_freq,n_full_freq, &
         omega,&
         omega_full, womega_full, &
         chemical_potential_spin, &
         n_electrons,occ_numbers, &
         KS_eigenvalue, &
         sosex_selfenergy &
        )

!  PURPOSE 
!    Calculate the 2nd-order exchange self-energy at the imaginary frequency grid

!  USES

      use constants
      use dimensions
      use prodbas
!      use physics
      use hartree_fock
      use species_data
      use mpi_tasks
      use mpi_utilities
      use synchronize_mpi
      use runtime_choices
      use geometry
      use constants
      use timing
      use evaluate_polarisability_freq, only: evaluate_polarisability_freq_0
      use localorb_io, only: use_unit
      implicit none

!  ARGUMENTS
!     input
      integer ::   n_freq
      integer ::   n_full_freq
      integer ::   n_low_state
      integer ::   n_KS_states
      real*8  ::   omega(n_freq)
      real*8  ::   omega_full(n_full_freq)
      real*8  ::   womega_full(n_full_freq)
      real*8  ::   ovlp_3KS(n_loc_prodbas,n_states,n_states,n_spin)
      real*8  ::   chemical_potential_spin(n_spin)   
      real*8  ::   n_electrons
      real*8  ::   occ_numbers(n_states,n_spin,n_k_points)
      real*8  ::   KS_eigenvalue(n_states,n_spin,n_k_points)
!     output
      complex*16   sosex_selfenergy(n_freq,n_KS_states,n_spin)

! INPUTS
! o  n_freq -- integer number, the number of frequency points for the GW self-energy
! o  n_low_state  -- integer number, the lowest KS/HF eigenstate for self-energy correction
! o  n_KS_state -- ineger number,
!    the highest KS/HF eigenstate for self-energy correction
! o  omega(n_freq) -- the Gauss-Legendre frequency grid for the self-energy, from 0 to a maximal value
! o  omega_full(n_freq) -- real array
!           the Gauss-Legendre frequency grid for the screened coulomb matrix, from 0 to infinity
! o  womega_full(n_freq) -- real array
!           the weight for omega_full grid
! o  ovlp_3KS -- real array
!    this is the transformed 3-cener overlap integration. Now two orbitals of
!    them are KS ones, and one is the auxiliary basis.
!    Note: for parallel calculations, the auxiliary basis are distribuated
!            among the different processors.
! o  chemical_potential -- real number, the chemical potential of the system
! o  n_electrons -- number of electrons in the system
! o  occ_numbers -- occupation numbers of the eigenstates
! o  KS_eigenvalue -- real array,
!    the eigenvalues of the single-particle calculation. For DFT calculation,
!    this is the KS eigenvalue, but for HF calculation, this is then the HF 
!            eigenvalue

! OUTPUT
! o  2ox_self_energy -- complex array
!            the calculated MP2 self-energy on the imaginary frequency axis

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

      real*8, dimension(:,:), allocatable :: aux_eri
      real*8, dimension(:,:), allocatable :: aux_eri_2
      real*8, dimension(:,:), allocatable :: polar_freq
      real*8, dimension(:,:,:), allocatable :: dielec_func
      real*8, dimension(:), allocatable :: aux_ovlp3KS_vec
      real*8, dimension(:), allocatable :: ovlp_multi_w
      real*8, dimension(:,:), allocatable :: eri_over_ediff
      real*8, dimension(:), allocatable :: aux_vect_real
      real*8, dimension(:), allocatable :: aux_vect_imag
      real*8, dimension(:,:), allocatable :: ovlp_tmp
      complex*16, dimension(:), allocatable :: sosex_tmp

      complex*16 :: e_deno_1, e_deno_2
      complex*16 :: tmp
      real*8 :: ddot

      character*50 :: filename

      integer ::  n_lumo(n_spin)
      integer ::  n_lumo_min
      integer ::  n_unocc(n_spin)
      integer ::  n_unocc_max
      integer ::  n_unocc_eff

      integer ::  n_shell_aux(n_species)


!     Accuracy
!      real*8 :: nu=1.d-12
!      complex*16 :: inu=(0.d0,1.d-20)


!     Time counters

!      real*8 clock_time_2ox



!     Integer counters
      integer i_spin
      integer i_spin_1
      integer i_state
      integer i_state_2
      integer i_occ
      integer i_virt
      integer i_prodbas
      integer i_basbas
      integer i_start
      integer i_index
      integer i_index_1
      integer i_freq
      integer i_freq_1


!    determine the HOMO and LUMO level

!      if(n_KS_states .lt. n_states) then
!        if(myid.eq.0) then
!          write(use_unit,'(2X,A)') "Error: SOSEX self-energy calculation -- it", &
!           " is mandatory to set 'state_higher_limit = <n_states>' "
!        endif
!        call aims_stop('')
!      endif

      if (.not.allocated(n_homo)) then
         allocate(n_homo(n_spin))
      endif
      do i_spin = 1, n_spin
       do i_state = 1, n_states
         if(occ_numbers(i_state,i_spin,1).gt.1.d-6) then
           n_homo(i_spin)=i_state
         else
          exit
         endif
       enddo

       do i_state =  n_states, 1, -1
         if(occ_numbers(i_state,i_spin,1).lt.1.d0) then
           n_lumo(i_spin)=i_state
         endif
       enddo
      enddo

      n_homo_max = max(n_homo(1), n_homo(n_spin))
      n_lumo_min = min(n_lumo(1), n_lumo(n_spin))

      if(myid.eq.0) then
        write(use_unit,'(A)') 
        write(use_unit,'(A)') "---------------------------------------------"
        write(use_unit,'(2A)') " | Second-order screened exchange self-energy", &
                   " calculation starts ..."
        write(use_unit,'(A)') 

        write(use_unit,*) " | Homo level                             ", n_homo
        write(use_unit,*) " | Lumo level                             ", n_lumo
        write(use_unit,*) " | Total Number of states                 ", n_states

      endif
      if (n_spin.eq.2) then
       do i_spin=1,n_spin,1
        if (i_spin.eq.1) then
         if(myid.eq.0) then
          write(use_unit,'(2X,A,1X,F16.5,A,1X,F16.5,A)') &
         "| Homo-Lumo Gap of the spin up    ", &
            KS_eigenvalue(n_lumo(i_spin),i_spin,1) &
            - KS_eigenvalue(n_homo(i_spin),i_spin,1)," Ha", &
           (KS_eigenvalue(n_lumo(i_spin),i_spin,1) &
           - KS_eigenvalue(n_homo(i_spin),i_spin,1))*hartree," eV"
         endif
        else if (n_electrons.gt.1) then
          if(myid.eq.0) then
           write(use_unit,'(2X,A,1X,F16.5,A,1X,F16.5,A)') &
           "| Homo-Lumo Gap of the spin down  ", &
           KS_eigenvalue(n_lumo(i_spin),i_spin,1) &
           - KS_eigenvalue(n_homo(i_spin),i_spin,1)," Ha", &
           (KS_eigenvalue(n_lumo(i_spin),i_spin,1) &
           - KS_eigenvalue(n_homo(i_spin),i_spin,1))*hartree," eV"
          endif
        endif
       enddo
      else
        if(myid.eq.0) then
         write(use_unit,'(2X,A,1X,F16.5,A,1X,F16.5,A)') &
          "| Homo-Lumo Gap                  ", &
         KS_eigenvalue(n_lumo(n_spin),n_spin,1) &
        - KS_eigenvalue(n_homo(n_spin),n_spin,1)," Ha", &
         (KS_eigenvalue(n_lumo(n_spin),n_spin,1) &
         - KS_eigenvalue(n_homo(n_spin),n_spin,1))*hartree," eV"
        endif
      endif

      if(myid.eq.0) then
        write(use_unit,'(A)') "---------------------------------------------"
      endif


!      if (n_spin.eq.1.and.n_homo_max.eq.1) then
!       write(use_unit,*) &
!               " * Only one occupied state: SOSEX self-energy is zero! "
!       return
!      endif

      n_unocc(:) = n_states-n_lumo(:)+1
      n_unocc_max = max(n_unocc(1),n_unocc(n_spin))

      if (.not.allocated(aux_eri)) then
          allocate(aux_eri(n_unocc_max,n_homo_max))
      endif
      if (.not.allocated(aux_eri_2)) then
          allocate(aux_eri_2(n_unocc_max,n_homo_max))
      endif
      if (.not.allocated(eri_over_ediff)) then
          allocate(eri_over_ediff(n_unocc_max,2))
      endif
      if (.not.allocated(aux_vect_real)) then
          allocate(aux_vect_real(n_loc_prodbas))
      endif
      if (.not.allocated(aux_vect_imag)) then
          allocate(aux_vect_imag(n_loc_prodbas))
      endif
      if (.not.allocated(aux_eri_2)) then
          allocate(aux_eri_2(n_homo_max,n_homo_max))
      endif
      if (.not.allocated(polar_freq)) then
          allocate(polar_freq(n_basbas,n_loc_prodbas))
      endif
      if (.not.allocated(dielec_func)) then
          allocate(dielec_func(n_basbas,n_loc_prodbas,n_full_freq))
      endif
      if (.not.allocated(aux_ovlp3KS_vec)) then
          allocate(aux_ovlp3KS_vec(n_basbas))
      endif
      if (.not.allocated(ovlp_tmp)) then
          allocate(ovlp_tmp(n_loc_prodbas,n_unocc_max))
      endif
      if (.not.allocated(ovlp_multi_w)) then
          allocate(ovlp_multi_w(n_loc_prodbas))
      endif
      if (.not.allocated(sosex_tmp)) then
          allocate(sosex_tmp(n_full_freq))
      endif

      do i_freq = 1, n_full_freq, 1
        call evaluate_polarisability_freq_0 &
           ( n_low_state, n_homo, n_lumo, n_KS_states, &
             occ_numbers, omega_full(i_freq), &
             KS_eigenvalue, ovlp_3KS(:,:,1:n_KS_states,:), polar_freq &
           )

        call screened_coulomb_interaction(polar_freq,dielec_func(1,1,i_freq))
!        dielec_func(:,:,i_freq) = 0.d0
!        do i_prodbas=1, n_loc_prodbas
!          i_basbas = map_prodbas(i_prodbas, myid+1)
!          if(i_basbas.gt.0) then
!             dielec_func(i_basbas,i_prodbas,i_freq)=1.d0
!          endif
!        enddo
      enddo

      sosex_selfenergy = (0.d0,0.d0)
      do i_spin = 1, n_spin, 1

          i_index = 0
          do i_state = n_low_state, n_KS_states, 1
             if(myid.eq.0) then
               write(use_unit,'(A,I4,I6)') " | i_spin, i_state = ", i_spin, i_state
             endif
             do i_state_2 = 1, n_states, 1
!               if(i_state .le. n_homo(i_spin) .and. i_state_2.le.n_homo(i_spin) .or. &
!                 i_state .gt. n_lumo(i_spin) .and. i_state_2.gt.n_lumo(i_spin) ) cycle

              aux_ovlp3KS_vec(:)=0.d0
              do i_prodbas = 1, n_loc_prodbas, 1
                i_basbas = map_prodbas(i_prodbas,myid+1)
                if(i_basbas .gt. 0) then
                   aux_ovlp3KS_vec(i_basbas) = ovlp_3KS(i_prodbas,i_state_2,i_state,i_spin)
                endif
              enddo
              call sync_vector(aux_ovlp3KS_vec,n_basbas)

              sosex_tmp(:) = (0.d0,0.d0)
              do i_spin_1 = 1, n_spin, 1
               if(i_spin_1 .ne. i_spin) cycle
               aux_eri(:,:)=0.d0
               aux_eri_2(:,:)=0.d0
               if(n_loc_prodbas .gt. 0) then
		   call dgemm('T', 'N', n_unocc(i_spin), n_homo(i_spin), &
			 n_loc_prodbas, 1.0d0, &
			 ovlp_3KS(:,n_lumo(i_spin):n_states,i_state_2,i_spin), &
			 n_loc_prodbas, &
			 ovlp_3KS(:,1:n_homo(i_spin),i_state,i_spin), &
			 n_loc_prodbas,0.d0, &
			 aux_eri(1,1), n_unocc_max)

		   call dgemm('T', 'N', n_unocc(i_spin), n_homo(i_spin), &
			 n_loc_prodbas, 1.0d0, &
			 ovlp_3KS(:,n_lumo(i_spin):n_states,i_state,i_spin), &
			 n_loc_prodbas, &
			 ovlp_3KS(:,1:n_homo(i_spin),i_state_2,i_spin), &
			 n_loc_prodbas,0.d0, &
			 aux_eri_2(1,1), n_unocc_max)
              endif

              call sync_matrix(aux_eri,n_unocc_max,n_homo_max)
              call sync_matrix(aux_eri_2,n_unocc_max,n_homo_max)
!              aux_eri_2=aux_eri

              do i_freq = 1, n_full_freq

! Multiplying ovlp_3KS with screened Coulomb potential
                call dgemv('T', n_basbas, n_loc_prodbas, &
                      1.d0, dielec_func(1,1,i_freq), n_basbas, &
                      aux_ovlp3KS_vec, 1, 0.d0, ovlp_multi_w, 1)

                do i_occ = 1, n_homo(i_spin_1)

                  i_start = max(i_occ,n_lumo(i_spin_1))
                  n_unocc_eff = n_states - i_start +1
! this correctly accounts for fractional occupation case. When fractional occupations occur,
! n_lumo <= n_homo. The orbitals pairs involving fractional occupations should only be 
! accounted for once.
                  do i_virt = i_start, n_states, 1

                     i_index = i_virt - i_start + 1
                     i_index_1 = i_virt - n_lumo(i_spin_1) + 1

                     e_deno_1 = (0.d0,1.d0)*omega_full(i_freq)  + &
                             KS_eigenvalue(i_occ,i_spin_1,1) - &
                             KS_eigenvalue(i_virt,i_spin_1,1)

                     e_deno_2 = (0.d0,1.d0)*omega_full(i_freq)  + &
                             KS_eigenvalue(i_virt,i_spin_1,1) - &
                             KS_eigenvalue(i_occ,i_spin_1,1)
                    
                     eri_over_ediff(i_index,1) =  &
                       (occ_numbers(i_occ,i_spin_1,1) - occ_numbers(i_virt,i_spin_1,1)) * &
                       real(aux_eri(i_index_1,i_occ)/e_deno_1 - aux_eri_2(i_index_1,i_occ)/e_deno_2) &
                       * dble(n_spin)/2.d0

                     eri_over_ediff(i_index,2) =  &
                       (occ_numbers(i_occ,i_spin_1,1) - occ_numbers(i_virt,i_spin_1,1)) * &
                       imag(aux_eri(i_index_1,i_occ)/e_deno_1 - aux_eri_2(i_index_1,i_occ)/e_deno_2) &
                       * dble(n_spin)/2.d0
!                     if(i_state_2.le.n_homo(i_spin)) then
!                       eri_over_ediff(i_index,1) =  &
!                         (occ_numbers(i_occ,i_spin_1,1) - occ_numbers(i_virt,i_spin_1,1)) * &
!                          real(aux_eri(i_index_1,i_occ)/e_deno_1) 
!
!                       eri_over_ediff(i_index,2) =  &
!                         (occ_numbers(i_occ,i_spin_1,1) - occ_numbers(i_virt,i_spin_1,1)) * &
!                         imag(aux_eri(i_index_1,i_occ)/e_deno_1)
!                     else
!                       eri_over_ediff(i_index,1) =  &
!                         -(occ_numbers(i_occ,i_spin_1,1) - occ_numbers(i_virt,i_spin_1,1)) * &
!                           real(aux_eri_2(i_index_1,i_occ)/e_deno_2) 
!
!                       eri_over_ediff(i_index,2) =  &
!                         -(occ_numbers(i_occ,i_spin_1,1) - occ_numbers(i_virt,i_spin_1,1)) * &
!                          imag(aux_eri_2(i_index_1,i_occ)/e_deno_2) 
!                     endif
                  enddo
! end of loop over do i_virt

!                  ovlp_tmp(:,1:n_unocc_eff)=ovlp_3KS(:,i_start:n_states,i_occ,i_spin_1)
                  call dgemv('N', n_loc_prodbas, n_unocc_eff, &
                      1.d0, ovlp_3KS(1:n_loc_prodbas,i_start:n_states,i_occ,i_spin_1), &
!                      1.d0, ovlp_tmp, &
                      n_loc_prodbas,  eri_over_ediff(1:n_unocc_eff,1), 1, 0.d0, aux_vect_real, 1)

                  call dgemv('N', n_loc_prodbas, n_unocc_eff, &
                      1.d0, ovlp_3KS(1:n_loc_prodbas,i_start:n_states,i_occ,i_spin_1), &
!                      1.d0, ovlp_tmp, &
                      n_loc_prodbas,  eri_over_ediff(1:n_unocc_eff,2), 1, 0.d0, aux_vect_imag, 1)
!                   aux_vect(:,:) = 0.d0
!                   do i_basbas = 1, n_loc_prodbas, 1
!                     do i_virt = i_start, n_states, 1
!                         aux_vect(i_basbas,1) = aux_vect(i_basbas,1) + &
!                         ovlp_3KS(i_basbas,i_virt,i_occ,i_spin) *  &
!                         eri_over_ediff(i_virt,1)
!                         aux_vect(i_basbas,2) = aux_vect(i_basbas,2) + &
!                         ovlp_3KS(i_basbas,i_virt,i_occ,i_spin) *  &
!                         eri_over_ediff(i_virt,2)
!                    enddo
!                  enddo

                  sosex_tmp(i_freq) =  sosex_tmp(i_freq) + &
                      ddot(n_loc_prodbas,aux_vect_real,1,ovlp_multi_w,1) + &
                      (0.d0,1.d0)*ddot(n_loc_prodbas,aux_vect_imag,1,ovlp_multi_w,1)
!                  write(use_unit,*) "eri_over_ediff:", eri_over_ediff(:,:) 
!                  write(use_unit,*) "aux_vect:", aux_vect(:,:) 
!                  write(use_unit,*) "ovlp_multi_w:", ovlp_multi_w(:) 
!                  write(use_unit,*) i_freq, i_occ, i_virt, i_state_2, i_state
!                  write(use_unit,*) sosex_tmp(i_freq)


! end of loop over do i_occ
            enddo
! end of loop over i_freq
           enddo
! end of loop over i_spin_1
           enddo
!           if(i_state .eq.1) then
!             write(use_unit,*) i_state, i_state_2
!              do i_freq = 1, n_full_freq, 1
!               write(use_unit,'(2X, f16.8, 7X, 2F18.10)') omega_full(i_freq), sosex_tmp(i_freq)
!              enddo
!              write(use_unit,*) 
!           endif
!   perform the frequency convolution. i_freq_1 loops over the frequency grid for the self-energy
           do i_freq_1 = 1, n_freq, 1
             e_deno_1 = (0.d0,1.d0)*omega(i_freq_1) + chemical_potential_spin(i_spin) &
                        - KS_eigenvalue(i_state_2,i_spin,1)
             tmp = (0.d0,0.d0)
             do i_freq = 1, n_full_freq, 1

              sosex_selfenergy(i_freq_1,i_state,i_spin) =  sosex_selfenergy(i_freq_1,i_state,i_spin) + &
                      ( (sosex_tmp(i_freq)-conjg(sosex_tmp(i_freq_1)))/(e_deno_1+(0.d0,1.d0)*omega_full(i_freq)) + &
                        (conjg(sosex_tmp(i_freq))-conjg(sosex_tmp(i_freq_1)))/(e_deno_1-(0.d0,1.d0)*omega_full(i_freq)) &
                      ) * womega_full(i_freq)

!              if(i_freq_1.eq.1 .and. i_freq .eq.1 .and. i_state .ge.3 .and. i_state .le.5 ) then
!                write(use_unit,'(3I4,4f16.8)') i_spin, i_state, i_state_2, sosex_tmp(1), sosex_selfenergy(1,i_state,i_spin)
!              endif
             
!   end of loop over i_freq
            enddo

            if(KS_eigenvalue(i_state_2,i_spin,1) .le. chemical_potential_spin(i_spin)) then
              sosex_selfenergy(i_freq_1,i_state,i_spin) = sosex_selfenergy(i_freq_1,i_state,i_spin) + &
                      conjg(sosex_tmp(i_freq_1))*pi
            else
              sosex_selfenergy(i_freq_1,i_state,i_spin) = sosex_selfenergy(i_freq_1,i_state,i_spin) - &
                      conjg(sosex_tmp(i_freq_1))*pi
            endif
! end of loop over i_freq_1
          enddo
 
         

!   end of loop over i_state_2
        enddo
!   end of loop over i_state
       enddo
!    end of loop over i_spin
      enddo

      sosex_selfenergy(:,:,:)=sosex_selfenergy(:,:,:)/2.d0/pi
!      sosex_selfenergy(:,:,:)=sosex_selfenergy(:,:,:)/pi
      if(use_mpi) then
        do i_spin = 1, n_spin, 1
          call sync_matrix_complex &
             (sosex_selfenergy(:,:,i_spin),n_freq,n_KS_states)

        enddo
      endif

!      do i_spin = 1, n_spin, 1
!       do i_state = 1, 1, 1
!        write(use_unit,*) "i_state: ", i_state
!        do i_freq = 1, n_freq, 1
!          write(use_unit,'(2X, f16.8, 7X, 2F18.10)') omega(i_freq), sosex_selfenergy(i_freq,i_state,i_spin)
!        enddo
!        write(use_unit,*) 
!       enddo
!      enddo

!      Ending

      if (allocated(aux_eri)) then
        deallocate (aux_eri)
      endif
      if (allocated(aux_eri_2)) then
        deallocate (aux_eri_2)
      endif
      if (allocated(eri_over_ediff)) then
        deallocate (eri_over_ediff)
      endif
      if (allocated(aux_vect_real)) then
        deallocate (aux_vect_real)
      endif
      if (allocated(aux_vect_imag)) then
        deallocate (aux_vect_imag)
      endif
      if (allocated(polar_freq)) then
        deallocate (polar_freq)
      endif
      if (allocated(dielec_func)) then
        deallocate (dielec_func)
      endif
      if (allocated(aux_ovlp3KS_vec)) then
        deallocate (aux_ovlp3KS_vec)
      endif
      if (allocated(ovlp_multi_w)) then
        deallocate (ovlp_multi_w)
      endif
      if (allocated(ovlp_tmp)) then
        deallocate (ovlp_tmp)
      endif

      return

      end subroutine evaluate_sosex_selfenergy


