!****s* FHI-aims/evaluate_sox_selfenergy
!  NAME
!   evaluate_sox_selfenergy
!  SYNOPSIS

      subroutine evaluate_sox_selfenergy &
       ( n_freq, n_low_state, n_KS_states, &
         omega, &
         ovlp_3KS, &
         chemical_potential_spin, &
         n_electrons,occ_numbers, &
         KS_eigenvalue, &
         sox_selfenergy &
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
      use localorb_io, only: use_unit
      implicit none

!  ARGUMENTS
!     input
      integer ::   n_freq
      integer ::   n_low_state
      integer ::   n_KS_states
      real*8  ::   omega(n_freq)
      real*8  ::   ovlp_3KS(n_loc_prodbas,n_states,n_KS_states,n_spin)
      real*8  ::   chemical_potential_spin(n_spin)   
      real*8  ::   n_electrons
      real*8  ::   occ_numbers(n_states,n_spin,n_k_points)
      real*8  ::   KS_eigenvalue(n_states,n_spin,n_k_points)
!     output
      complex*16   sox_selfenergy(n_freq,n_KS_states,n_spin)

! INPUTS
! o  n_freq -- integer number, the number of frequency points for the GW self-energy
! o  n_low_state  -- integer number, the lowest KS/HF eigenstate for self-energy correction
! o  n_KS_state -- ineger number,
!    the highest KS/HF eigenstate for self-energy correction
! o  omega(n_freq) -- real array
!           the Gauss-Legendre frequency grid for the self-energy
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

      real*8 :: aux
      real*8, dimension(:,:), allocatable :: aux_eri
      real*8, dimension(:,:), allocatable :: aux_eri_2
      real*8 :: E_2ox_a, E_2ox_b, E_2ox_mpi

      complex*16 :: e_deno

      character*50 :: filename

      integer ::  n_lumo(n_spin)
      integer ::  n_lumo_min

      integer ::  n_shell_aux(n_species)

!     output
      real*8 :: E_2ox, E_2ox_1

!     Accuracy
!      real*8 :: nu=1.d-12
!      complex*16 :: inu=(0.d0,1.d-20)


!     Time counters

!      real*8 clock_time_2ox



!     Integer counters
      integer i_spin
      integer i_state
      integer i_state_2
      integer i_occ
      integer i_occ_2
      integer i_virt
      integer i_virt2
      integer i_basbas
      integer i_species
      integer i_atom
      integer term
      integer i_task
      integer i_index
      integer n_unoccupied
      integer i_start_b
      integer i_freq


!    determine the HOMO and LUMO level

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
!       if(occ_numbers(n_homo(i_spin),i_spin,1).eq.dble(2/n_spin))then
         n_lumo(i_spin) = n_homo(i_spin) +1
!       else
!         n_lumo(i_spin) = n_homo(i_spin)
!       endif
      enddo
!      endif

       if (n_spin.eq.2.and.n_lumo(n_spin).lt.n_lumo(1)) then
        n_lumo_min=n_lumo(n_spin)
        else
         n_lumo_min=n_lumo(1)
       endif


      if(myid.eq.0) then
        write(use_unit,'(A)') 
        write(use_unit,'(A)') "---------------------------------------------"
        write(use_unit,'(2A)') " | Second-order exchange self-energy", &
                   " calculation starts ..."
        write(use_unit,'(A)') "---------------------------------------------"
        if (n_spin.gt.1) then
          write(use_unit,'(A)') "  | Note : Spin-polarized system "
        else
          write(use_unit,'(A)') "  | Note : Non spin-polarized system "
        endif

      write(use_unit,'(A)') "---------------------------------------------"


!      write(use_unit,*) " | Homo level                             ", n_homo
!      write(use_unit,*) " | Lumo level                             ", n_lumo
!      write(use_unit,*) " | Total Number of states                 ", n_states
!      write(use_unit,*) " | Number of k-points                     ",n_k_points

!  endif myid
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

      n_homo_max = max(n_homo(1), n_homo(n_spin))
      if(myid.eq.0) then
        write(use_unit,'(A)') "---------------------------------------------"
      endif


!      if (n_spin.eq.1.and.n_homo_max.eq.1) then
!       write(use_unit,*) &
!              " * Only one occupied state: 2OX correction is zero! "
!      return
!     endif

      n_unoccupied = n_states-n_lumo_min+1

!      if(myid.eq.0) then
!        write(use_unit,'(A)')
!        write(use_unit,'(A)')"---------------------------------------------"
!        write(use_unit,'(A)')"| Sum of correlation terms "
!        write(use_unit,'(A)')"---------------------------------------------"
!      endif

      if (.not.allocated(aux_eri)) then
          allocate(aux_eri(n_lumo_min:n_states,n_lumo_min:n_states))
      endif
      if (.not.allocated(aux_eri_2)) then
          allocate(aux_eri_2(n_homo_max,n_homo_max))
      endif

! -------------------------------------------------------------------
!    NON-Spin polarized systems
! -------------------------------------------------------------------

      sox_selfenergy = (0.d0,0.d0)
      if (n_spin.eq.1) then

          i_index = 0
          do i_state = n_low_state, n_KS_states, 1
             if(myid.eq.0) then
               write(use_unit,*) " | i_state = ", i_state
             endif
             do i_state_2 = 1, n_homo(n_spin), 1

              call dgemm('T', 'N', n_unoccupied, n_unoccupied, &
                     n_loc_prodbas, 1.0d0, &
                     ovlp_3KS(:,n_lumo_min:n_states,i_state,n_spin), &
                     n_loc_prodbas, &
                     ovlp_3KS(:,n_lumo_min:n_states,i_state_2,n_spin), &
                     n_loc_prodbas,0.d0, &
                     aux_eri(n_lumo_min:n_states,n_lumo_min:n_states), &
                     n_unoccupied &
                    )

              call sync_matrix &
                   (aux_eri(n_lumo_min:n_states,n_lumo_min:n_states), &
                    n_unoccupied,n_unoccupied)

              i_index = i_index + 1
! MPI task distribution
              if(myid.eq.MOD(i_index,n_tasks)) then

!              if(i_state .le. n_homo(i_spin)) then

               do i_virt=n_lumo(n_spin),n_states,1
                do i_virt2=n_lumo(n_spin),n_states,1

                  E_2ox_a=aux_eri(i_virt,i_virt2)
                  E_2ox_b=aux_eri(i_virt2,i_virt)
                  E_2ox_a= -E_2ox_b* E_2ox_a

                  do i_freq = 1, n_freq
!                   if(myid.eq.0) then
!                    write(use_unit,'(A,I4,f20.6)') " | i_freq ", i_freq, & 
!                     omega(i_freq)
!                   endif

                    e_deno = (0.d0,1.d0)*omega(i_freq)  + &
                             chemical_potential_spin(n_spin) + &
                             KS_eigenvalue(i_state_2,n_spin,1) - &
                             KS_eigenvalue(i_virt,n_spin,1)- &
                             KS_eigenvalue(i_virt2,n_spin,1)

                     if (abs(e_deno).lt.1e-6)then

                       write(use_unit,'(10X,A)') &
                            "****************************************"
                       write(use_unit,'(10X,2A)') "| Warning :", &
                         " too close to degeneracy"
                       write(use_unit,'(10X,A)') &
                           "****************************************"
                     endif

                     e_deno=  1.d0/e_deno
                     sox_selfenergy (i_freq, i_state, n_spin) = &
                        sox_selfenergy (i_freq, i_state, n_spin) + &
                        E_2ox_a*e_deno

!             end of loop over i_freq
                   enddo
                    e_deno = KS_eigenvalue(i_state,n_spin,1) + &
                             KS_eigenvalue(i_state_2,n_spin,1) - &
                             KS_eigenvalue(i_virt,n_spin,1)- &
                             KS_eigenvalue(i_virt2,n_spin,1)
!                  endif
!   close b  virtual state 2
                enddo
!   close a virtual state
               enddo
              endif
             enddo

             do i_state_2 = n_lumo(n_spin), n_states, 1

              aux_eri_2(:,:) = 0.d0
              call dgemm('T', 'N', n_homo(n_spin), &
                      n_homo(n_spin), &
                      n_loc_prodbas, 1.0d0, &
                      ovlp_3KS(:,i_state,1:n_homo(n_spin),n_spin), &
                      n_loc_prodbas, &
                      ovlp_3KS(:,i_state_2,1:n_homo(n_spin),n_spin), &
                      n_loc_prodbas,0.d0, &
                      aux_eri_2,n_homo_max &
                     )

              call sync_matrix &
                     (aux_eri_2, n_homo_max,n_homo_max)

              if(myid.eq.MOD(i_state_2,n_tasks)) then

               do i_occ = 1, n_homo(n_spin), 1
                do i_occ_2 = 1, n_homo(n_spin), 1

                  E_2ox_a=aux_eri_2(i_occ, i_occ_2)
                  E_2ox_b=aux_eri_2(i_occ_2,i_occ)
                  E_2ox_a= -E_2ox_b * E_2ox_a

                  do i_freq = 1, n_freq
                    e_deno = (0.d0,1.d0)*omega(i_freq) + &
                             chemical_potential_spin(n_spin) + &
                             KS_eigenvalue(i_state_2,n_spin,1) - &
                             KS_eigenvalue(i_occ,n_spin,1)- &
                             KS_eigenvalue(i_occ_2,n_spin,1)

                     if (abs(e_deno).lt.1e-6)then

                       write(use_unit,'(10X,A)') &
                            "****************************************"
                       write(use_unit,'(10X,2A)') "| Warning :", &
                         " too close to degeneracy"
                       write(use_unit,'(10X,A)') &
                           "****************************************"
                     endif

                     e_deno=  1.d0/e_deno
                     sox_selfenergy (i_freq, i_state, n_spin) = &
                        sox_selfenergy (i_freq, i_state, n_spin) + &
                        E_2ox_a*e_deno


!             end of loop over i_freq
                   enddo

                    e_deno = KS_eigenvalue(i_state,n_spin,1) + &
                             KS_eigenvalue(i_state_2,n_spin,1) - &
                             KS_eigenvalue(i_occ,n_spin,1)- &
                             KS_eigenvalue(i_occ_2,n_spin,1)
!   end of loop over i_virt
                enddo
!   end of loop over i_occ
                enddo
!  end if i_state < n_homo
!               endif

!   end of MPI distribution
              endif
!                   if(i_state.eq.18) then
!                       write(use_unit,*) i_state, i_state_2,
!                   endif
!   end of loop over i_state_2
              enddo

!   close i state
        enddo

      else

! -------------------------------------------------------------------
!     Spin polarized systems
! -------------------------------------------------------------------

        i_index = 0
        do i_spin = 1, n_spin, 1
          do i_state = n_low_state, n_KS_states, 1
            if(myid.eq.0) then
              write(use_unit,*) " | i_spin, i_state = ", i_spin, i_state
            endif

!            do i_spin_2 = 1, n_spin, 1
                do i_state_2 = 1, n_homo(i_spin), 1

                  call dgemm('T', 'N', n_unoccupied, n_unoccupied, &
                         n_loc_prodbas, 1.0d0, &
                         ovlp_3KS(:,n_lumo_min:n_states,i_state, &
                                       i_spin), n_loc_prodbas, &
                         ovlp_3KS(:,n_lumo_min:n_states,i_state_2, &
                                  i_spin), n_loc_prodbas, 0.d0, &
                         aux_eri(n_lumo_min:n_states, &
                                 n_lumo_min:n_states), &
                         n_unoccupied &
                        )

                  call sync_matrix( aux_eri(n_lumo_min:n_states, &
                                  n_lumo_min:n_states),n_unoccupied, &
                                  n_unoccupied &
                                 )

               i_index = i_index + 1
! MPI task distribution
               if(myid.eq.MOD(i_index,n_tasks)) then

               do i_virt = n_lumo(i_spin), n_states, 1
                do i_virt2 = n_lumo(i_spin), n_states, 1
!                    if (i_state.eq.i_state_2.and.i_virt.eq.i_virt2.and.
!     +                  i_spin.eq.i_spin_2) cycle

                        E_2ox_a=aux_eri(i_virt,i_virt2)

!                        E_2ox_b=0.d0

!                        if (i_spin.eq.i_spin_2) then

                         E_2ox_b=aux_eri(i_virt2,i_virt)
                         E_2ox_a=-E_2ox_b*E_2ox_a
!                        else
!                              E_2ox_a=(E_2ox_a)*E_2ox_a

!                        endif

                        do i_freq = 1, n_freq, 1
                           e_deno = (0.d0,1.d0)*omega(i_freq)  + &
                              chemical_potential_spin(i_spin) + &
                              KS_eigenvalue(i_state_2,i_spin,1) - &
                              KS_eigenvalue(i_virt,i_spin,1)- &
                              KS_eigenvalue(i_virt2,i_spin,1)

                           e_deno=  1.d0/e_deno
                           sox_selfenergy (i_freq, i_state, i_spin) = &
                              sox_selfenergy (i_freq, i_state, i_spin) + &
                              E_2ox_a*e_deno

!      end of loop over i_freq
                      enddo
                      e_deno = (0.d0,1.d0)*omega(1)+chemical_potential_spin(i_spin) + &
                              KS_eigenvalue(i_state_2,i_spin,1) - &
                              KS_eigenvalue(i_virt,i_spin,1)- &
                              KS_eigenvalue(i_virt2,i_spin,1)

!                      write(use_unit,'(2X,3I4,4f16.6)') i_virt2, i_virt, i_spin_2,  E_2ox_a/e_deno, sox_selfenergy(1,i_state,i_spin)


!   close b  virtual state 2
                enddo
!   close a virtual state
               enddo
!   end of MPI distribution
              endif
!   close j state2
              enddo

              do i_state_2 = n_lumo(i_spin), n_states, 1

                 aux_eri_2(:,:) = 0.d0
                 call dgemm('T', 'N', n_homo(i_spin), &
                        n_homo(i_spin), &
                        n_loc_prodbas, 1.0d0, &
                        ovlp_3KS(:,1:n_homo(i_spin),i_state,i_spin), &
                        n_loc_prodbas, &
                        ovlp_3KS(:,i_state_2,1:n_homo(i_spin), &
                                 i_spin), &
                        n_loc_prodbas, 0.d0, &
                        aux_eri_2, n_homo_max &
                       )

                 call sync_matrix &
                       (aux_eri_2, n_homo_max,n_homo_max)

                 if(myid.eq.MOD(i_state_2,n_tasks)) then

                   do i_occ = 1, n_homo(i_spin), 1
                     do i_occ_2 = 1, n_homo(i_spin), 1

                       E_2ox_a=aux_eri_2(i_occ, i_occ_2)

!                       if(i_spin .eq. i_spin_2) then
                       E_2ox_b=aux_eri_2(i_occ_2,i_occ)
                       E_2ox_a= -E_2ox_b * E_2ox_a
!                       else
!                       E_2ox_a= E_2ox_a* E_2ox_a
!                       endif

                       do i_freq = 1, n_freq
                         e_deno = (0.d0,1.d0)*omega(i_freq) + &
                               chemical_potential_spin(i_spin) + &
                               KS_eigenvalue(i_state_2,i_spin,1) - &
                               KS_eigenvalue(i_occ,i_spin,1)- &
                               KS_eigenvalue(i_occ_2,i_spin,1)

                         e_deno=  1.d0/e_deno

                         sox_selfenergy (i_freq, i_state, i_spin) = &
                            sox_selfenergy (i_freq, i_state, i_spin) + &
                            E_2ox_a*e_deno

!             end of loop over i_freq
                     enddo
!                      e_deno = (0.d0,1.d0)*omega(1)+chemical_potential_spin(i_spin) + &
!                              KS_eigenvalue(i_state_2,i_spin,1) - &
!                              KS_eigenvalue(i_virt,i_spin,1)- &
!                              KS_eigenvalue(i_virt2,i_spin,1)

!                      write(use_unit,'(2X,3I4,4f16.6)') i_virt2, i_virt, i_spin_2,  E_2ox_a/e_deno, sox_selfenergy(1,i_state,i_spin)
!   end of loop over i_occ_2
              enddo
!   end of loop over i_occ
           enddo

!  MPI task distribution
          endif
!   close i state_2
          enddo

!  end of loop over i_spin2
!         enddo
!         write(use_unit,*) sox_selfenergy(1,i_state,i_spin), sox_selfenergy(2,i_state,i_spin)
!   close i state
        enddo
       enddo

      endif

      if(use_mpi) then
        do i_spin = 1, n_spin, 1
          call sync_matrix_complex &
             (sox_selfenergy(:,:,i_spin),n_freq,n_KS_states)

!          sox_selfenergy(:,:,i_spin) = 0.5d0*
!     +       sox_selfenergy(:,:,i_spin)
        enddo
      endif

!      if(myid==0) then
!        write(use_unit,*) '--------------------------------------------'
!        write(use_unit,*)
!        write(use_unit,'(10X,A,F12.3,A,F12.3,A)') &
!            "| Total time for calculating 2OX selfenergy          : ", &
!             time_2ox, " s ",clock_time_2ox," s"
!        write(use_unit,*)
!      endif
!      do i_spin = 1, n_spin, 1
!       do i_state = 1, 1, 1
!        do i_freq = 1, n_freq, 1
!          write(use_unit,'(2X, f16.8, 7X, 2F18.10)') omega(i_freq), -sox_selfenergy(i_freq,i_state,i_spin)
!        enddo
!       enddo
!      enddo


!      Ending

      if (allocated(aux_eri)) then
        deallocate (aux_eri)
      endif
      if (allocated(aux_eri_2)) then
        deallocate (aux_eri_2)
      endif

      return

      end subroutine evaluate_sox_selfenergy


!******
!---------------------------------------------------------------------------------------------------

      subroutine evaluate_sox_selfenergy_2 &
       ( n_freq, n_low_state, n_KS_states, &
         omega, &
         ovlp_3KS, &
         chemical_potential_spin, &
         n_electrons,occ_numbers, &
         KS_eigenvalue, &
         sox_selfenergy &
        )

!  PURPOSE 
!    Calculate the MP2 self-energy at the imaginary frequency grid

!  USES

      use constants
      use dimensions
      use prodbas
      use hartree_fock
      use species_data
      use mpi_tasks
      use mpi_utilities
      use synchronize_mpi
      use runtime_choices
      use geometry
      use constants
      use timing
      use localorb_io, only: use_unit
      implicit none

!  ARGUMENTS
!     input
      integer ::   n_freq
      integer ::   n_low_state
      integer ::   n_KS_states
      real*8  ::   omega(n_freq)
      real*8  ::   ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
      real*8  ::   chemical_potential_spin(n_spin)   
      real*8  ::   n_electrons
      real*8  ::   occ_numbers(n_states,n_spin,n_k_points)
      real*8  ::   KS_eigenvalue(n_states,n_spin,n_k_points)
!     output
      complex*16   sox_selfenergy(n_freq,n_KS_states,n_spin)

! INPUTS
! o  n_freq -- integer number, the number of frequency points for the GW self-energy
! o  n_low_state  -- integer number, the lowest KS/HF eigenstate for self-energy correction
! o  n_KS_state -- ineger number,
!    the highest KS/HF eigenstate for self-energy correction
! o  omega(n_freq) -- real array
!           the Gauss-Legendre frequency grid for the self-energy
! o  ovlp_3KS -- real array
!    this is the transformed 3-cener overlap integration. Now two orbitals of
!    them are KS ones, and one is the auxiliary basis.
!    Note: for parallel calculations, dimensions 2 and 3 are distribuated
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

      real*8 :: E_2ox_a, E_2ox_b

      real*8, allocatable :: tmp_o3KS(:,:)
      real*8, allocatable :: tmp_prod(:,:)
      real*8, allocatable :: tmp_eri_full(:,:,:,:)
      real*8, allocatable :: tmp_eri_tran(:,:,:)

      complex*16 :: e_deno

      integer :: n_lumo(n_spin)
      integer :: n_lumo_min

      integer :: min_lumo_loc, max_homo_loc, max_local_states1

!     Time counters


!     Integer counters
      integer i_spin
!      integer i_spin2
      integer i_state
      integer i_state2
      integer i_occ
      integer i_occ2
      integer i_virt
      integer i_virt2
      integer i_atom
      integer i_freq
      integer n_p1, n_p2
      integer i_virt_loc, i_virt2_loc, i_virt_own
      integer i_occ_loc, i_occ2_loc, i_occ_own
      integer i_state2_loc
      integer j, mpierr


!    determine the HOMO and LUMO level

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
       n_lumo(i_spin) = n_homo(i_spin) +1
      enddo

      n_lumo_min = minval(n_lumo)
      n_homo_max = maxval(n_homo)


      if(myid.eq.0) then
        write(use_unit,'(A)') 
        write(use_unit,'(A)') "---------------------------------------------"
        write(use_unit,'(A)') " | MP2 quasiparticle calculation starts ..."
        write(use_unit,'(A)') "---------------------------------------------"
        if (n_spin.gt.1) then
          write(use_unit,'(A)') "  | Note : Spin-polarized system "
        else
          write(use_unit,'(A)') "  | Note : Non spin-polarized system "
        endif
        do i_atom=1, n_atoms,1
         if (empty(i_atom)) then
           write(use_unit,'(A)') "  | Note : Counterpoise correction"
          exit
         endif
        enddo

        write(use_unit,'(A)') "---------------------------------------------"


        write(use_unit,*) " | Homo level                             ", n_homo
        write(use_unit,*) " | Lumo level                             ", n_lumo
        write(use_unit,*) " | Total Number of states                 ", n_states
        write(use_unit,*) " | Number of k-points                     ",n_k_points
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
!               " * Only one occupied state: MP2 correction is zero! "
!       return
!      endif

!      call get_timestamps(time_2ox, clock_time_2ox )


      ! ATTENTION: This code requires that ndim1_o3KS/ndim2_o3KS is globally identical
      ! and that the complete ovlp_3KS hasn't undefined values (e.g. NaN)

      sox_selfenergy = (0.d0,0.d0)

      min_lumo_loc = loc_dim1_o3ks(n_lumo_min) ! overall minimal local value corresponding to n_lumo_min
      max_homo_loc = loc_dim2_o3ks(n_homo_max) ! overall maximal local value corresponding to n_homo_max

      max_local_states1 = ndim1_o3KS - min_lumo_loc + 1

      allocate(tmp_o3KS(n_basbas,min_lumo_loc:ndim1_o3KS))
      allocate(tmp_eri_full(min_lumo_loc:ndim1_o3KS,min_lumo_loc:ndim1_o3KS,0:np1_o3KS-1,max_homo_loc))
      allocate(tmp_eri_tran(min_lumo_loc:ndim1_o3KS,min_lumo_loc:ndim1_o3KS,0:np1_o3KS-1))

      do i_spin = 1, n_spin, 1
        do i_state = n_low_state, n_KS_states, 1

          if(myid.eq.0) write(use_unit,*) " | Virt i_spin, i_state = ", i_spin, i_state

!          do i_spin2 = 1, n_spin, 1

            n_p2 = own_dim2_o3ks(i_state)

            do n_p1 = 0, np1_o3KS-1

              ! Proc (n_p1,n_p2) broadcasts its local part of
              ! ovlp_3KS_full(:,n_lumo_min:n_states,i_state,i_spin)

              if(n_p1==myp1_o3KS .and. n_p2==myp2_o3KS) &
                tmp_o3KS(:,:) = ovlp_3KS(:,min_lumo_loc:ndim1_o3KS,loc_dim2_o3ks(i_state),i_spin)

              call MPI_Bcast(tmp_o3KS,n_basbas*max_local_states1,MPI_REAL8,global_id(n_p1,n_p2),mpi_comm_global,mpierr)

              ! Multiply what we got with our local part of ovlp_3KS_full(:,n_lumo_min:n_states,1:n_homo(i_spin2),i_spin2)

              do i_state2 = 1, n_homo(i_spin), 1

                if(own_dim2_o3ks(i_state2) /= myp2_o3KS) cycle
                i_state2_loc = loc_dim2_o3ks(i_state2)

                call dgemm('T', 'N', max_local_states1, max_local_states1, n_basbas, 1.d0, &
                           tmp_o3KS, ubound(tmp_o3KS,1), &
                           ovlp_3KS(1,min_lumo_loc,i_state2_loc,i_spin), ubound(ovlp_3KS,1), &
                           0.d0, tmp_eri_full(min_lumo_loc,min_lumo_loc,n_p1,i_state2_loc),max_local_states1)
              enddo

            enddo
                

            do i_state2 = 1, n_homo(i_spin), 1

              ! The original code:
              !
              !n_unoccupied = n_states-n_lumo_min+1
              !call dgemm('T', 'N', n_unoccupied, n_unoccupied, &
              !       n_loc_prodbas, 1.0d0, &
              !       ovlp_3KS(:,n_lumo_min:n_states,i_state,n_spin), &
              !       n_loc_prodbas, &
              !       ovlp_3KS(:,n_lumo_min:n_states,i_state2,n_spin), &
              !       n_loc_prodbas,0.d0, &
              !       aux_eri(n_lumo_min:n_states,n_lumo_min:n_states), &
              !       n_unoccupied &
              !      )

              if(own_dim2_o3ks(i_state2) /= myp2_o3KS) cycle
              i_state2_loc = loc_dim2_o3ks(i_state2)

              ! Transpose the strip we have in tmp_eri_full:

!              if(i_spin .eq. i_spin2) &
              call MPI_Alltoall(tmp_eri_full(min_lumo_loc,min_lumo_loc,0,i_state2_loc), &
                                max_local_states1*max_local_states1,MPI_REAL8, &
                                tmp_eri_tran,max_local_states1*max_local_states1,MPI_REAL8, &
                                mpi_comm_rows_aux_2d,mpierr)

              do i_virt=n_lumo(i_spin),n_states,1

                do i_virt2=n_lumo(i_spin),n_states,1

                  if(own_dim1_o3ks(i_virt2) /= myp1_o3ks) cycle

                  i_virt_loc  = loc_dim1_o3ks(i_virt)
                  i_virt_own  = own_dim1_o3ks(i_virt)
                  i_virt2_loc = loc_dim1_o3ks(i_virt2)

!                  if (i_spin.eq.i_spin2) then
                    E_2ox_a=tmp_eri_full(i_virt_loc,i_virt2_loc,i_virt_own,i_state2_loc)
                    E_2ox_b=tmp_eri_tran(i_virt2_loc,i_virt_loc,i_virt_own)
                    E_2ox_a= -E_2ox_b * E_2ox_a
!                  else
!                    E_2ox_a=tmp_eri_full(i_virt_loc,i_virt2_loc,i_virt_own,i_state2_loc)
!                    E_2ox_a=(E_2ox_a)*E_2ox_a
!                     E_2ox_a=0.d0
!                  endif

                  do i_freq = 1, n_freq

                    e_deno = (0.d0,1.d0)*omega(i_freq)  + &
                             chemical_potential_spin(i_spin) + &
                             KS_eigenvalue(i_state2,i_spin,1) - &
                             KS_eigenvalue(i_virt,i_spin,1)- &
                             KS_eigenvalue(i_virt2,i_spin,1)

!                     if (abs(e_deno).lt.1e-6)then
!                       write(use_unit,'(10X,A)') "****************************************"
!                       write(use_unit,'(10X,A)') "| Warning : too close to degeneracy"
!                       write(use_unit,'(10X,A)') "****************************************"
!                     endif

                     e_deno=  1.d0/e_deno
                     sox_selfenergy (i_freq, i_state, i_spin) = &
                       sox_selfenergy (i_freq, i_state, i_spin) + &
                       E_2ox_a*e_deno

                  enddo ! i_freq
                enddo ! i_virt2
              enddo ! i_virt
            enddo ! i_state2
!          enddo ! i_spin2
        enddo ! i_state
      enddo ! i_spin

      deallocate(tmp_o3KS)
      deallocate(tmp_eri_full)
      deallocate(tmp_eri_tran)

      allocate(tmp_o3KS(n_basbas,max_homo_loc))
      allocate(tmp_prod(max_homo_loc,min_lumo_loc:ndim1_o3KS))
      allocate(tmp_eri_full(max_homo_loc,max_homo_loc,0:np2_o3KS-1,min_lumo_loc:ndim1_o3KS))
      allocate(tmp_eri_tran(max_homo_loc,max_homo_loc,0:np2_o3KS-1))

      do i_spin = 1, n_spin, 1
        do i_state = n_low_state, n_KS_states, 1

          if(myid.eq.0) write(use_unit,*) " | Occ  i_spin, i_state = ", i_spin, i_state

!          do i_spin2 = 1, n_spin, 1

            n_p1 = own_dim1_o3ks(i_state)

            do n_p2 = 0, np2_o3KS-1

              ! Proc (n_p1,n_p2) gathers and broadcasts its local part of ovlp_3KS_full(:,i_state,1:n_homo_max,n_spin)

              if(n_p1==myp1_o3KS .and. n_p2==myp2_o3KS) &
                tmp_o3KS(:,:) = ovlp_3KS(:,loc_dim1_o3ks(i_state),1:max_homo_loc,i_spin)

              call MPI_Bcast(tmp_o3KS,n_basbas*max_homo_loc,MPI_REAL8,global_id(n_p1,n_p2),mpi_comm_global,mpierr)

              ! Multiply what we got with our local part of ovlp_3KS_full(:,n_lumo_min:n_states,1:n_homo_max,n_spin)

              do j = 1, max_homo_loc
                call dgemm('T', 'N', max_homo_loc, ndim1_o3KS-min_lumo_loc+1, n_basbas, 1.d0, &
                           tmp_o3KS, ubound(tmp_o3KS,1), &
                           ovlp_3KS(1,min_lumo_loc,j,i_spin), ubound(ovlp_3KS,1), &
                           0.d0, tmp_prod, ubound(tmp_prod,1))
                
                ! Scatter the columns of the result into tmp_eri_full

                do i_state2 = n_lumo(i_spin), n_states, 1
                  if(own_dim1_o3ks(i_state2) /= myp1_o3KS) cycle
                  i_state2_loc = loc_dim1_o3ks(i_state2)
                  tmp_eri_full(:,j,n_p2,i_state2_loc) =  tmp_prod(:,i_state2_loc)
                enddo
              enddo

            enddo

            do i_state2 = n_lumo(i_spin), n_states, 1

              ! The original code:
              !
              !call dgemm('T', 'N', n_homo(n_spin), &
              !        n_homo(n_spin), &
              !        n_loc_prodbas, 1.0d0, &
              !        ovlp_3KS(:,i_state,1:n_homo(n_spin),n_spin), &
              !        n_loc_prodbas, &
              !        ovlp_3KS(:,i_state2,1:n_homo(n_spin),n_spin), &
              !        n_loc_prodbas,0.d0, &
              !        aux_eri_2,n_homo_max &
              !       )

              if(own_dim1_o3ks(i_state2) /= myp1_o3KS) cycle
              i_state2_loc = loc_dim1_o3ks(i_state2)

              ! Transpose the strip we have in tmp_eri_full:

!              if(i_spin .eq. i_spin2) &
              call MPI_Alltoall(tmp_eri_full(1,1,0,i_state2_loc),max_homo_loc*max_homo_loc,MPI_REAL8, &
                                tmp_eri_tran,max_homo_loc*max_homo_loc,MPI_REAL8, &
                                mpi_comm_cols_aux_2d,mpierr)

              do i_occ = 1, n_homo(i_spin), 1

                do i_occ2 = 1, n_homo(i_spin), 1

                  if(own_dim2_o3ks(i_occ2) /= myp2_o3ks) cycle

                  i_occ_loc  = loc_dim2_o3ks(i_occ)
                  i_occ_own  = own_dim2_o3ks(i_occ)
                  i_occ2_loc = loc_dim2_o3ks(i_occ2)

!                  if(i_spin .eq. i_spin2) then
                    E_2ox_a=tmp_eri_full(i_occ_loc,i_occ2_loc,i_occ_own,i_state2_loc)
                    E_2ox_b=tmp_eri_tran(i_occ2_loc,i_occ_loc,i_occ_own)
                    E_2ox_a= -E_2ox_b*E_2ox_a
!                  else
!                    E_2ox_a=tmp_eri_full(i_occ_loc,i_occ2_loc,i_occ_own,i_state2_loc)
!                    E_2ox_a= E_2ox_a*E_2ox_a
!                     E_2ox_a = 0.d0
!                  endif

                  do i_freq = 1, n_freq
                    e_deno = (0.d0,1.d0)*omega(i_freq) + &
                             chemical_potential_spin(i_spin) + &
                             KS_eigenvalue(i_state2,i_spin,1) - &
                             KS_eigenvalue(i_occ,i_spin,1)- &
                             KS_eigenvalue(i_occ2,i_spin,1)

!                     if (abs(e_deno).lt.1e-6)then
!                       write(use_unit,'(10X,A)') "****************************************"
!                       write(use_unit,'(10X,A)') "| Warning : too close to degeneracy"
!                       write(use_unit,'(10X,A)') "****************************************"
!                     endif

                     e_deno=  1.d0/e_deno
                     sox_selfenergy (i_freq, i_state, i_spin) = &
                        sox_selfenergy (i_freq, i_state, i_spin) + &
                        E_2ox_a*e_deno

                  enddo ! i_freq
                enddo ! i_occ2
              enddo ! i_occ
            enddo ! i_state2
!          enddo ! i_spin2
        enddo ! i_state
      enddo ! i_spin

      deallocate(tmp_o3KS)
      deallocate(tmp_prod)
      deallocate(tmp_eri_full)
      deallocate(tmp_eri_tran)

      if(use_mpi) then
        do i_spin = 1, n_spin, 1
          call sync_matrix_complex &
             (sox_selfenergy(:,:,i_spin),n_freq,n_KS_states)
        enddo
      endif


!      call get_times(time_2ox, clock_time_2ox)
!      if(myid==0) then
!        write(use_unit,*) '--------------------------------------------'
!        write(use_unit,*)
!        write(use_unit,'(10X,A,F12.3,A,F12.3,A)') &
!            "| Total time MP2 selfenergy                      : ", &
!             time_2ox, " s ",clock_time_2ox," s"
!        write(use_unit,*)
!      endif

      return

      end subroutine evaluate_sox_selfenergy_2

!******

