!****s* FHI-aims/evaluate_en_mp2_correlation_energy
!  NAME
!   evaluate_en_mp2_correlation_energy
!  SYNOPSIS 

      subroutine evaluate_en_mp2_correlation_energy &
       (n_electrons,total_energy,en_xc, &
        occ_numbers,KS_eigenvalue,KS_eigenvector,post_scf_total_energy )

!  PURPOSE
!    Calculation of the total energy perturbation terms up to the second order using Moller-Plesset theory (MP2) to the HF/DFT Energy.
!    A variant, called SCS-MP2, scales the two spin channels separately.
!    Frozen core approximation is available too. 
! USES

      use constants
      use dimensions
      use prodbas
      use species_data
      use physics, only : BSSE_full_energy, BSSE_per_atom ! for BSSE
      use mpi_tasks
      use mpi_utilities
      use synchronize_mpi
      use runtime_choices
      use geometry
      use timing
      use hartree_fock, only: coulomb_matr_lvl, coeff_3fn_ten, ovlp_3fn
      use sparse_tensor, only: dealloc_sp_ten
      use localorb_io, only: localorb_info, OL_norm, use_unit
      implicit none

!  ARGUMENTS

      real*8  :: n_electrons
      real*8  :: total_energy
      real*8  :: en_xc
      real*8  :: occ_numbers(n_states,n_spin,n_k_points)
      real*8  :: KS_eigenvalue(n_states,n_spin,n_k_points)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin,n_k_points)
      real*8  :: post_scf_total_energy
  
!  INPUTS
! o  n_homo -- HOMO level
! o  total_energy -- HF total energy
! o  en_xc -- exchange-correlation energy with a KS reference
! o  occ_numbers -- the occupation number of the electrons for each eigenstate and each spin
! o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
! o  KS_eigenvector -- real array,
!            the eigenvector of the single-particle calculation
! o  post_scf_total_energy -- the final MP2 (or scs-MP2) total energy to be printed out
!                              at the end of the output file

!  OUTPUTS

!  local variables

      real*8, dimension(:,:,:,:), allocatable :: ovlp_3KS
      real*8 :: n_bytes_o3KS
      real*8 :: aux

      real*8, dimension(:,:), allocatable :: aux_eri

      real*8 :: E_mp2_a, E_mp2_b, E_mp2_mpi
      real*8 :: E_mp2_a_scs, E_mp2_mpi_scs
      real*8 :: E_mp2_term, E_mp2_term_scs

      real*8 :: eex_energy
      real*8 :: en_total


      integer  :: n_homo(n_spin)
      integer ::  n_lumo(n_spin)
      integer ::  n_lumo_min
      integer ::  n_homo_max

      integer ::  n_shell_aux(n_species)

      character*20 :: filename

      real*8 :: E_mp2, E_mp2_1, E_mp2_scs


!     Accuracy
!      real*8 :: nu=1.d-12
!      complex*16 :: inu=(0.d0,1.d-20)


!     Time counters

      real*8 time_3basis_intg
!      real*8 time_3KS_intg
      real*8 time_mp2

      real*8 tot_time_3basis_intg
!      real*8 tot_time_3KS_intg
      real*8 tot_time_mp2
!   bsse time keeping
      real*8  temp_time_bsse
      real*8  temp_clock_time_bsse

!     Integer counters
      integer i_spin
      integer i_spin2
      integer i_state
      integer i_state2
      integer i_virt
      integer i_virt2
      integer i_basis_1
      integer i_species
      integer i_atom
      integer term
      integer i_empty(n_species)
      integer i_task
      integer i_index
      integer n_unoccupied
      integer i_start_b

!     Error variable
      integer mpierr
      integer :: info
      character*150 :: info_str
      character(*), parameter :: func = 'evaluate_mp2_correlation_energy'
!     Igor
      real*8, dimension(:,:,:), allocatable :: E_mp2_spectrum, E_mp2_spectrum_mpi
      real*8, dimension(:,:), allocatable :: E_mp2_spectrum_1D_mpi
      ! (aa|bb) (ab|bb) (ii|aa) (jj|aa) (ia|ai) (ja|aj)
      real*8, dimension(:,:), allocatable :: EN_shift_a
      real*8, dimension(:,:), allocatable :: EN_shift_b
      real*8, dimension(:,:), allocatable :: EN_shift_c
      real*8, dimension(:,:), allocatable :: EN_shift_d
      real*8, dimension(:,:), allocatable :: EN_shift_e
      real*8, dimension(:,:), allocatable :: EN_shift_f
      real*8, dimension(:,:), allocatable :: aux_eri_test
      real*8, dimension(:,:,:,:), allocatable :: ovlp_3KS_test
!  parameters
!   o n_atoms -- number of atoms
!   o n_spin -- number of spin channels
!   o n_species -- number of species
!   o n_homo -- HOMO level for each spin
!   o n_lumo -- LUMO level for each spin
!   o occ_numbers -- occupation numbers
!   o empty -- ghost atoms levels
!   o i_start_mp2 -- starting level out of the frozen core
!   o KS_eigenvalue -- Kohn-Sham eigenvalues 
!   o ovlp_3KS -- overlap three centers integrals
!   o total_energy -- energy of the ground state zero+first order
!   o eex_energy -- energy of the ground state zero+first order of exact exchange
!   o pt_mp2 -- parallel spin channels scaling in SCS approximation
!   o ps_mp2 -- opposite spin channels scaling in SCS approximation
!  OUTPUT
!   o E_mp2 -- total energy correction 
!   o tot_time_mp2 -- time to calculate MP2 correction
!  WARNING
!   SCS-MP2 feature not fully tested.
!  NOTES
!   RI-MP2 approximation is sensible to product basis cut-off potential and number of
!   basis set products, so be carefull if higher accuracy is needed.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society.
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




! Assumig the HOMO level is already correctly determined
      n_homo(:) = 0
      do i_spin = 1, n_spin
       do i_state = 1, n_states
         if(occ_numbers(i_state,i_spin,1).gt.1.d-6) then
           n_homo(i_spin)=i_state
         else
          exit
         endif
       enddo
      enddo

      if (n_spin .eq. 2 .and. use_hf_multiplicity) then
        n_homo (1) =  int((n_electrons+1)/2.d0) +  &
                     (hf_multiplicity - 1)/2
        n_homo (2) =  int((n_electrons)/2.d0) - &
                     (hf_multiplicity - 1)/2
        occ_numbers(1:n_homo(1),1,1) = 1.d0
        occ_numbers(n_homo(1)+1:n_states,1,1) = 0.d0
        occ_numbers(1:n_homo(2),2,1) = 1.d0
        occ_numbers(n_homo(2)+1:n_states,2,1) = 0.d0
      endif

      n_lumo(:) = n_homo(:) + 1

      if (n_spin.eq.2.and.n_lumo(n_spin).lt.n_lumo(1)) then
       n_lumo_min=n_lumo(n_spin)
       else
        n_lumo_min=n_lumo(1)
      endif

      if (n_lumo_min .gt. n_states) then
        e_mp2 = 0.d0
        if(myid.eq.0) then
          write(use_unit,'(2X,2A)')  &
          " * Warning: There is no unoccupied states, and hence ", &
          "the MP2 correlation energy is zero."
        endif
        return
      endif

      if(myid.eq.0) then
        write(use_unit,'(A)') "  | Note : Resolution of identity applied "
        if (use_scsmp2) then
          write(use_unit,'(A)') "  | Note : Spin scaled MP2 (SCS-MP2) "
        endif
        if (n_spin.gt.1) then
          write(use_unit,'(A)') "  | Note : Spin-polarized system "
        else
          write(use_unit,'(A)') "  | Note : Non spin-polarized system "
        endif
        if (flag_frozen_core.and.i_start_mp2.eq. &
           0.and.n_electrons.gt.2) then
          write(use_unit,'(A)') "  | Note : Frozen Core approximation  "
        else if (flag_frozen_core.and.n_electrons.gt.2) then
          write(use_unit,'(A)') "  | Note : Frozen Core approximation from level "
          write(use_unit,*) "               ",   i_start_mp2
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

!  endif myid
      endif
      if (n_spin.eq.2) then
       do i_spin=1,n_spin,1
        if (i_spin.eq.1) then
         if(myid.eq.0) then
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
         "| Homo-Lumo Gap of the spin up    ", &
            KS_eigenvalue(n_lumo(i_spin),i_spin,1) &
            - KS_eigenvalue(n_homo(i_spin),i_spin,1)," Ha", &
           (KS_eigenvalue(n_lumo(i_spin),i_spin,1) &
           - KS_eigenvalue(n_homo(i_spin),i_spin,1))*hartree," eV"
         endif
        else if (n_electrons.gt.1) then
          if(myid.eq.0) then
           write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
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
         write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
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


!     This is the implicit flag to determine whether FC treatment is used!
      if (allocated(n_fc_shell)) then

        ! if we're determining the FC orbitals automatically from last complete
        ! noble gas shell:
        if ((i_start_mp2.eq.0).and.(n_electrons.gt.2)) then

         ! Determine how many empty sites there are for this species ...

!        This here should count correctly.
         i_empty=0
         do i_atom = 1, n_atoms, 1
           if (empty(i_atom)) then
             i_empty(species(i_atom)) = i_empty(species(i_atom)) + 1
           end if
         enddo

         do i_species=1,n_species
           n_shell_aux(i_species)=n_fc_shell(i_species)-1
         enddo

         i_start_b = 0
         do i_species=1,n_species
!test
!      if (myid.eq.0) then
!        write(use_unit,*) "Species ", i_species, ":",
!     +  " atoms: ", atoms_in_structure(i_species),
!     +  " empty sites: ", i_empty(i_species)
!      end if
!test end
            i_start_mp2=n_shell_aux(i_species) &
               +3*max(0,(n_shell_aux(i_species)-1)) &
               +5*max(0,(n_shell_aux(i_species)-3)) &
               +7*max(0,(n_shell_aux(i_species)-5))
            i_start_mp2=i_start_mp2*(atoms_in_structure(i_species) &
               -i_empty(i_species))+i_start_b
            i_start_b=i_start_mp2
         enddo
         i_start_mp2=i_start_mp2+1

         if (myid.eq.0) then
          write(use_unit,'(2X,A,I5)') &
          "Number of automatically determined frozen core levels: ", &
          i_start_mp2-1
          do i_spin = 1, n_spin, 1
           if (n_spin.gt.1) then
             write(use_unit,'(2X,A,I2,A)') &
             "Spin channel", i_spin, ":"
           end if
           write(use_unit,'(2X,A,I5,1X,I5)') &
           "MP2 contribution will be calculated between levels: ", &
           i_start_mp2, n_homo(i_spin)
          enddo
          end if

         endif

      endif

      if (n_electrons.le.2) then
       i_start_mp2=1
      endif

      n_homo_max = max(n_homo(1),n_homo(n_spin))

!     check on frozen core compatibility
       if (i_start_mp2.gt.n_homo_max) then
        write(use_unit,'(1X,A)')"*--------------------------------------------"
        write(use_unit,'(1X,A)')"* Error: Requested more frozen core levels"
        write(use_unit,'(1X,A)')"* than occupied states. Aborting MP2 part."
        write(use_unit,'(1X,A)')"*--------------------------------------------"
        stop
       endif

     ! determine the post-SCF frozen core low state by igor
     if (flag_frozen_core_postSCF) then ! count the frozen core states
         call count_frozen_core_states(i_start_mp2)
         if(myid.eq.0) then
            write(use_unit,'(2X,A,I12)') &
            "Frozen core number for post-SCF calculation          :", &
            i_start_mp2 - 1
         endif
     endif



!      call cpu_time(time_3KS_intg)

      call get_timestamps(time_ovlp3fn_to_ovlp3KS, clock_time_ovlp3fn_to_ovlp3KS)
      if(.not. allocated(ovlp_3KS)) then
        if(use_hartree_fock .and. .not. sparse_o3fn) then
           allocate(ovlp_3KS(n_loc_prodbas,n_lumo_min:n_states, &
           &                 i_start_mp2:n_homo_max,n_spin), stat=info)
           call check_allocation(info, 'ovlp_3KS', func)
           n_bytes_o3KS = 8 * n_loc_prodbas * (n_states - n_lumo_min + 1)  &
                            * (n_homo_max - i_start_mp2) * n_spin
           ! igor ovlp_3KS_test
           allocate(ovlp_3KS_test(n_loc_prodbas,i_start_mp2:n_states, &
           &                 i_start_mp2:n_states,n_spin), stat=info)
           call check_allocation(info, 'ovlp_3KS_test', func)
        else
           allocate(ovlp_3KS(n_loc_prodbas, n_states, &
           &                 1:n_homo_max, n_spin), stat=info)
           n_bytes_o3KS = 8 * n_loc_prodbas * n_states * n_homo_max * n_spin
        end if

        write(info_str, "(2X,A,I7,A,I7,A,F12.3,A)") &
        & 'The ovlp_3KS matrix takes another', &
        & nint(n_bytes_o3KS / 2.d0**20), ' MiB x', n_tasks, ' procs =', &
        & n_tasks * n_bytes_o3KS / 2.d0**30, ' GiB.'
        call localorb_info(info_str, use_unit, '(A)', OL_norm)
        call localorb_info('', use_unit, '(A)', OL_norm)

        if(use_hartree_fock .and. .not. sparse_o3fn) then
            write(use_unit,'(2X,A)') "Debug : Igor 1a"
           call evaluate_ovlp_3MO(KS_eigenvector, &
                      ovlp_3KS,n_lumo_min,n_lumo)
            write(use_unit,'(2X,A)') "Debug : Igor 1b"
           call evaluate_ovlp_3MO_EN(KS_eigenvector, &
                      ovlp_3KS_test,n_lumo_min,n_lumo)
            write(use_unit,'(2X,A)') "Debug : Igor 1c"
        else
           if (sparse_o3fn) then
            write(use_unit,'(2X,A)') "Debug : Igor 2"
              call ovlp3KS_lvl_1d(n_homo_max, KS_eigenvector, &
              &                   coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
              deallocate(coulomb_matr_lvl)
              call dealloc_sp_ten(coeff_3fn_ten)
           else
            write(use_unit,'(2X,A)') "Debug : Igor 3"
              call transform_ovlp3fn(n_homo_max, KS_eigenvector, ovlp_3fn, &
              &                      ovlp_3KS)
              if (.not. calculate_atom_bsse) then
                deallocate(ovlp_3fn)
              endif
           end if
           call evaluate_eex_energy &
              ( n_homo_max,n_homo,occ_numbers, &
                ovlp_3KS,eex_energy &
              )
        endif
      endif
      call get_times(time_ovlp3fn_to_ovlp3KS, clock_time_ovlp3fn_to_ovlp3KS)

      call cpu_time(time_mp2)
      n_unoccupied = n_states-n_lumo_min+1

      if(myid.eq.0) then
        write(use_unit,'(A)')
        write(use_unit,'(A)')"---------------------------------------------"
        write(use_unit,'(A)')"| Sum of correlation terms "
        write(use_unit,'(A)')"---------------------------------------------"
      endif


      if (.not.use_scsmp2) then

      E_mp2=0.d0
      term=i_start_mp2
      ! Igor for E_mp2_spectrum
      allocate(E_mp2_spectrum(n_states,n_states,n_spin))
      allocate(E_mp2_spectrum_mpi(n_states,n_states,n_spin))
      allocate(E_mp2_spectrum_1D_mpi(n_states,n_spin))
      do i_spin=1,n_spin,1
          do i_state=1,n_states,1
              E_mp2_spectrum_1D_mpi(i_state,i_spin) = 0.0d0
              do i_state2=1,n_states,1
                  E_mp2_spectrum(i_state,i_state2,i_spin) = 0.0d0
                  E_mp2_spectrum_mpi(i_state,i_state2,i_spin) = 0.0d0
              enddo
          enddo
      enddo


      if (.not.allocated(aux_eri)) then
          allocate(aux_eri(n_lumo_min:n_states,n_lumo_min:n_states))
      endif

      if (.not.allocated(aux_eri_test)) then
          allocate(aux_eri_test(i_start_mp2:n_states,i_start_mp2:n_states))
      endif

      if (.not.allocated(EN_shift_a)) then
          allocate(EN_shift_a(n_lumo_min:n_states,n_lumo_min:n_states))
      endif
      if (.not.allocated(EN_shift_b)) then
          allocate(EN_shift_b(n_lumo_min:n_states,n_lumo_min:n_states))
      endif
      if (.not.allocated(EN_shift_c)) then
          allocate(EN_shift_c(n_lumo_min:n_states,n_lumo_min:n_states))
      endif
      if (.not.allocated(EN_shift_d)) then
          allocate(EN_shift_d(n_lumo_min:n_states,n_lumo_min:n_states))
      endif
      if (.not.allocated(EN_shift_e)) then
          allocate(EN_shift_e(n_lumo_min:n_states,n_lumo_min:n_states))
      endif
      if (.not.allocated(EN_shift_f)) then
          allocate(EN_shift_f(n_lumo_min:n_states,n_lumo_min:n_states))
      endif

! -------------------------------------------------------------------
!    NON-Spin polarized systems
! -------------------------------------------------------------------

     if (n_spin.eq.1) then
          do i_virt=n_lumo_min, n_states,1
            do i_virt2=i_virt, n_states,1
              call dgemm('T', 'N', 1, 1, &
                      n_loc_prodbas, 1.0d0, &
                      ovlp_3KS_test(:,i_virt,i_virt,n_spin), &
                      n_loc_prodbas, &
                      ovlp_3KS_test(:,i_virt2,i_virt2,n_spin),&
                      n_loc_prodbas, 0.d0, &
                      EN_shift_a(i_virt,i_virt2), &
                      n_unoccupied &
                     )
              call dgemm('T', 'N', 1, 1, &
                      n_loc_prodbas, 1.0d0, &
                      ovlp_3KS_test(:,i_virt2,i_virt,n_spin), &
                      n_loc_prodbas, &
                      ovlp_3KS_test(:,i_virt,i_virt2,n_spin), &
                      n_loc_prodbas, 0.d0, &
                      EN_shift_b(i_virt,i_virt2), &
                      n_unoccupied &
                     )
            enddo
          enddo
          call sync_matrix( EN_shift_a(n_lumo_min:n_states, &
                            n_lumo_min:n_states),n_unoccupied, &
                            n_unoccupied &
                           )
          call sync_matrix( EN_shift_b(n_lumo_min:n_states, &
                            n_lumo_min:n_states),n_unoccupied, &
                            n_unoccupied &
                           )
          i_index = 0
          do i_state=i_start_mp2, n_homo(n_spin),1
             E_mp2_term = 0.d0
             do i_state2=i_state,n_homo(n_spin),1
              !call dgemm('T', 'N', n_unoccupied, n_unoccupied, &
              !        n_loc_prodbas, 1.0d0, &
              !        ovlp_3KS(:,n_lumo_min:n_states,i_state,n_spin), &
              !        n_loc_prodbas, &
              !        ovlp_3KS(:,n_lumo_min:n_states,i_state2, &
              !        n_spin), n_loc_prodbas, 0.d0, &
              !        aux_eri(n_lumo_min:n_states, &
              !                n_lumo_min:n_states), &
              !        n_unoccupied &
              !       )

              !call sync_matrix( aux_eri(n_lumo_min:n_states, &
              !                  n_lumo_min:n_states),n_unoccupied, &
              !                  n_unoccupied &
              !                 )
              ! igor ovlp_3KS_test
              call dgemm('T', 'N', n_unoccupied, n_unoccupied, &
                      n_loc_prodbas, 1.0d0, &
                      ovlp_3KS_test(:,i_start_mp2:n_states,i_state,n_spin), &
                      n_loc_prodbas, &
                      ovlp_3KS_test(:,i_start_mp2:n_states,i_state2, &
                      n_spin), n_loc_prodbas, 0.d0, &
                      aux_eri_test(i_start_mp2:n_states, &
                              i_start_mp2:n_states), &
                      n_states &
                     )
              call sync_matrix( aux_eri_test(i_start_mp2:n_states, &
                                i_start_mp2:n_states),n_states-i_start_mp2+1, &
                                n_states-i_start_mp2+1 &
                               )

              call dgemm('T', 'N', n_unoccupied, n_unoccupied, &
                      n_loc_prodbas, 1.0d0, &
                      ovlp_3KS(:,n_lumo_min:n_states,i_state,n_spin), &
                      n_loc_prodbas, &
                      ovlp_3KS(:,n_lumo_min:n_states,i_state, &
                      n_spin), n_loc_prodbas, 0.d0, &
                      EN_shift_c(n_lumo_min:n_states, &
                              n_lumo_min:n_states), &
                      n_unoccupied &
                     )
              call sync_matrix( EN_shift_c(n_lumo_min:n_states, &
                                n_lumo_min:n_states),n_unoccupied, &
                                n_unoccupied &
                               )

              call dgemm('T', 'N', n_unoccupied, n_unoccupied, &
                      n_loc_prodbas, 1.0d0, &
                      ovlp_3KS(:,n_lumo_min:n_states,i_state2,n_spin), &
                      n_loc_prodbas, &
                      ovlp_3KS(:,n_lumo_min:n_states,i_state2, &
                      n_spin), n_loc_prodbas, 0.d0, &
                      EN_shift_d(n_lumo_min:n_states, &
                              n_lumo_min:n_states), &
                      n_unoccupied &
                     )
              call sync_matrix( EN_shift_d(n_lumo_min:n_states, &
                                n_lumo_min:n_states),n_unoccupied, &
                                n_unoccupied &
                               )
              i_index = i_index + 1
! MPI task distribution
              if(myid.eq.MOD(i_index,n_tasks)) then

               do i_virt=n_lumo(n_spin), n_states,1
                do i_virt2=n_lumo(n_spin),n_states,1
                  E_mp2_a=aux_eri_test(i_virt,i_virt2)
                  E_mp2_b=aux_eri_test(i_virt2,i_virt)
                  E_mp2_a= (2.d0*E_mp2_a- &
                             E_mp2_b)* E_mp2_a
                  !write(use_unit,'(2X,A,2I3,2F19.8)') &
                  !    "(ia|jb)",i_virt,i_virt2,&
                  !    aux_eri_test(i_virt,i_virt2),aux_eri(i_virt,i_virt2)

                  if (i_state.ne.i_state2) then
                     E_mp2_a=E_mp2_a*2.d0
                  endif

                  E_mp2_1=  (KS_eigenvalue(i_state,n_spin, 1)) &
                    +  (KS_eigenvalue(i_state2,n_spin,1) )

                  E_mp2_1= E_mp2_1-(KS_eigenvalue(i_virt,n_spin,1))

                  E_mp2_1= E_mp2_1-(KS_eigenvalue(i_virt2,n_spin,1))

                  if (en_shift_type .eq. 1) then
                      !write(use_unit,'(2X,A,2I3,2F19.8)') &
                      !    "(ii|jj) and (ij|ji)",i_state,i_state2,&
                      !    aux_eri_test(i_state,i_state2),aux_eri_test(i_state2,i_state)
                      !write(use_unit,'(2X,A,F19.8,2I3,2F19.8)') &
                      !    "E_mp2_1, (aa|bb) and (ab|ba)",E_mp2_1,i_virt,i_virt2,&
                      !    EN_shift_a(i_virt,i_virt2),EN_shift_b(i_virt,i_virt2)
                      E_mp2_1= E_mp2_1-aux_eri_test(i_state,i_state2)&
                          +aux_eri_test(i_state,i_state2)&
                          -EN_shift_a(i_virt,i_virt2) &
                          +EN_shift_b(i_virt,i_virt2) &
                          -EN_shift_c(i_virt,i_virt) &
                          -EN_shift_c(i_virt2,i_virt2) &
                          -EN_shift_d(i_virt,i_virt) &
                          -EN_shift_d(i_virt2,i_virt2)
                  elseif (en_shift_type .eq. 2) then
                      E_mp2_1= E_mp2_1+en_shift_constant
                  endif

                  if (abs(E_mp2_1).lt.1e-6)then
                     write(use_unit,'(10X,A)') &
                          "****************************************"
                     write(use_unit,'(10X,2A)') "| Warning :", &
                       " too close to degeneracy"
                     write(use_unit,'(10X,A)') &
                          "****************************************"
                  endif

                  E_mp2_1= 1.d0 / E_mp2_1
                  E_mp2_term= (E_mp2_a)*E_mp2_1 + E_mp2_term
                  ! Igor
                  E_mp2_spectrum(i_virt,i_virt2,1) = &
                      E_mp2_spectrum(i_virt,i_virt2,1) + E_mp2_a*E_mp2_1
!   close b  virtual state 2
                enddo
!   close a virtual state
               enddo
!   end of MPI distribution
              endif
!   close j state2
              enddo
              if(use_mpi) then
                  call MPI_ALLREDUCE(E_mp2_term, E_mp2_mpi, 1, &
                       MPI_DOUBLE_PRECISION, &
                      MPI_SUM, mpi_comm_global, mpierr)

                  E_mp2_term= E_mp2_mpi

              endif
              E_mp2 = E_mp2 + E_mp2_term
!      Ending

       if(myid.eq.0) then
         write(use_unit,'(2X,A,1X,I4,A,1X,F19.8,A,1X,F19.8,A)') &
                 "| Pair state ", term,"                    :" &
           , E_mp2_term * hartree, " eV", &
            E_mp2* hartree ," eV"
       endif
        term=term+1

!   close i state
        enddo

      ! Igor
      if (use_mpi) then
          call MPI_ALLREDUCE(E_mp2_spectrum, E_mp2_spectrum_mpi, &
              n_states*n_states*n_spin , &
              MPI_DOUBLE_PRECISION, &
             MPI_SUM, mpi_comm_global, mpierr)
     endif
      if(myid.eq.0) then
        write(use_unit,'(A)')
        write(use_unit,'(A)')"---------------------------------------------"
        write(use_unit,'(A)')"| Correlation energy spectrum "
        write(use_unit,'(A)')"---------------------------------------------"
        do i_virt = n_lumo(1), n_states,1
            do i_virt2 = n_lumo(1), n_states,1
                E_mp2_spectrum_1D_mpi(i_virt,1) = E_mp2_spectrum_1D_mpi(i_virt,1) &
                    + E_mp2_spectrum_mpi(i_virt,i_virt2,1)
                !write(use_unit,'(2X,A,1X,I4,I4,A,1X,F19.8,F19.8,F19.8)') &
                !    "| ",i_virt,i_virt2,'    :',KS_eigenvalue(i_virt,n_spin,1),&
                !    KS_eigenvalue(i_virt2,n_spin,1), E_mp2_spectrum_mpi(i_virt,i_virt2,1)
            enddo
            write(use_unit,'(2X,A,1X,I4,A,1X,F19.8,F19.8)') &
                "| ",i_virt,'    :',KS_eigenvalue(i_virt,n_spin,1),&
                E_mp2_spectrum_1D_mpi(i_virt,1)
        enddo
      endif

      else

! -------------------------------------------------------------------
!     Spin polarized systems
! -------------------------------------------------------------------

        i_index = 0
        do i_spin=1,n_spin,1
          do i_spin2=1,n_spin,1
             do i_state=i_start_mp2, n_homo(i_spin),1

                E_mp2_term = 0.d0
                do i_state2=i_start_mp2,n_homo(i_spin2) ,1

                  call dgemm('T', 'N', n_unoccupied, n_unoccupied, &
                         n_loc_prodbas, 1.0d0, &
                         ovlp_3KS(:,n_lumo_min:n_states,i_state, &
                                       i_spin), n_loc_prodbas, &
                         ovlp_3KS(:,n_lumo_min:n_states,i_state2, &
                                  i_spin2), n_loc_prodbas, 0.d0, &
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

               do i_virt=n_lumo(i_spin), n_states,1
                do i_virt2=n_lumo(i_spin2),n_states,1
                    if (i_state.eq.i_state2.and.i_virt.eq.i_virt2.and. &
                        i_spin.eq.i_spin2) cycle

                         E_mp2_a=aux_eri(i_virt,i_virt2)

                         E_mp2_b=0.d0


                         if (i_spin.eq.i_spin2) then

                              E_mp2_b=aux_eri(i_virt2,i_virt)

                              E_mp2_a=(E_mp2_a-E_mp2_b)*E_mp2_a
                        else
                              E_mp2_a=(E_mp2_a)*E_mp2_a

                        endif

                        E_mp2_1=  KS_eigenvalue(i_state,i_spin, 1) &
                           +  KS_eigenvalue(i_state2,i_spin2,1)

                        E_mp2_1=  E_mp2_1 - &
                                KS_eigenvalue(i_virt,i_spin,1)

                        E_mp2_1=  E_mp2_1 - &
                                KS_eigenvalue(i_virt2,i_spin2,1)




                       if (abs(E_mp2_1).lt.1e-6)then

                          write(use_unit,'(10X,A)') &
                           "****************************************"
                          write(use_unit,'(10X,2A)') "| Warning :", &
                              "too close to degeneracy"
                          write(use_unit,'(10X,A)') &
                           "****************************************"
                       endif


                       E_mp2_term= E_mp2_a/E_mp2_1+ E_mp2_term


!   close b  virtual state 2
                enddo
!   close a virtual state
               enddo
!   end of MPI distribution
              endif
!   close j state2
              enddo
              if(use_mpi) then
                   call MPI_ALLREDUCE(E_mp2_term, E_mp2_mpi, 1, &
                          MPI_DOUBLE_PRECISION, &
                          MPI_SUM, mpi_comm_global, mpierr)

                   E_mp2_term= E_mp2_mpi
              endif
              E_mp2_term= 0.5d0*E_mp2_term
              E_mp2 = E_mp2 + E_mp2_term

              if(myid.eq.0) then
               write(use_unit,'(2X,A,1X,I4,A,1X,F19.8,A,1X,F19.8,A)') &
                   "| Pair state ", term ,"                    :" &
                   , E_mp2_term * hartree, " eV" &
                  , E_mp2* hartree ," eV"
              endif

              term=term+1

!   close i state
         enddo

        enddo
       enddo

       endif

!      if(use_mpi) then
!         call MPI_ALLREDUCE(E_mp2, E_mp2_mpi, 1,
!     +        MPI_DOUBLE_PRECISION,
!     +        MPI_SUM, mpi_comm_global, mpierr)
!
!         E_mp2= E_mp2_mpi
!      endif
!      Ending

!********************************************************************
!    SCS-MP2
!********************************************************************

      else

      E_mp2=0.d0
      E_mp2_scs=0.d0
      term=i_start_mp2

      if (.not.allocated(aux_eri)) then
          allocate(aux_eri(n_lumo_min:n_states,n_lumo_min:n_states))
      endif

! -------------------------------------------------------------------
!    NON-Spin polarized systems
! -------------------------------------------------------------------

      if (n_spin.eq.1) then

          i_index = 0
          do i_state=i_start_mp2, n_homo(n_spin),1

             E_mp2_term = 0.d0
             E_mp2_term_scs = 0.d0
             do i_state2=i_state,n_homo(n_spin),1
!             do i_state2=1,n_homo(n_spin),1

              call dgemm('T', 'N', n_unoccupied, n_unoccupied, &
                      n_loc_prodbas, 1.0d0, &
                      ovlp_3KS(:,n_lumo_min:n_states,i_state,n_spin), &
                      n_loc_prodbas, &
                      ovlp_3KS(:,n_lumo_min:n_states,i_state2, &
                      n_spin), n_loc_prodbas, 0.d0, &
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

               do i_virt=n_lumo(n_spin), n_states,1
                do i_virt2=n_lumo(n_spin),n_states,1

                  E_mp2_a=aux_eri(i_virt,i_virt2)


!                  if (i_state.eq.i_state2.or.i_virt.eq.i_virt2) then
!                        E_mp2_a= E_mp2_a* E_mp2_a
!                  else
!                     E_mp2_b=0.d0

                    E_mp2_b=aux_eri(i_virt2,i_virt)

                E_mp2_a_scs= (2.d0*E_mp2_a*(pt_mp2+ps_mp2)- &
                            E_mp2_b*pt_mp2)* E_mp2_a


                E_mp2_a= (2.d0*E_mp2_a- &
                            E_mp2_b)* E_mp2_a

!                  endif

                  if (i_state.ne.i_state2) then
                     E_mp2_a=E_mp2_a*2.d0
                     E_mp2_a_scs=E_mp2_a_scs*2.d0

                  endif


                  E_mp2_1=  (KS_eigenvalue(i_state,n_spin, 1)) &
                    +  (KS_eigenvalue(i_state2,n_spin,1) )

                  E_mp2_1= E_mp2_1-(KS_eigenvalue(i_virt,n_spin,1))

                  E_mp2_1= E_mp2_1-(KS_eigenvalue(i_virt2,n_spin,1))




                   if (abs(E_mp2_1).lt.1e-6)then

                      write(use_unit,'(10X,A)') &
                           "****************************************"
                      write(use_unit,'(10X,2A)') "| Warning :", &
                        " too close to degeneracy"
                      write(use_unit,'(10X,A)') &
                           "****************************************"
                   endif

                   E_mp2_1= 1.d0 / E_mp2_1
                   E_mp2_term= (E_mp2_a)*E_mp2_1 + E_mp2_term

                  E_mp2_term_scs= (E_mp2_a_scs)*E_mp2_1 + E_mp2_term_scs

!   close b  virtual state 2
                enddo
!   close a virtual state
               enddo

!   end of MPI distribution
              endif
!   close j state2
              enddo
              if(use_mpi) then
                  call MPI_ALLREDUCE(E_mp2_term, E_mp2_mpi, 1, &
                       MPI_DOUBLE_PRECISION, &
                      MPI_SUM, mpi_comm_global, mpierr)

                  E_mp2_term= E_mp2_mpi

                  call MPI_ALLREDUCE(E_mp2_term_scs, E_mp2_mpi_scs, 1, &
                       MPI_DOUBLE_PRECISION, &
                      MPI_SUM, mpi_comm_global, mpierr)

                  E_mp2_term_scs= E_mp2_mpi_scs

              endif

              E_mp2 = E_mp2 + E_mp2_term
              E_mp2_scs = E_mp2_scs + E_mp2_term_scs
!      Ending

       if(myid.eq.0) then
         write(use_unit,'(2X,A,1X,I4,A,1X,F19.8,A,1X,F19.8,A)') &
                 "| Pair state ", term,      "                    :" &
           , E_mp2_term * hartree, " eV", &
            E_mp2* hartree ," eV"

         write(use_unit,'(2X,A,1X,I4,A,1X,F19.8,A,1X,F19.8,A)') &
                 "| SCS - Pair state ", term,"                    :" &
           , E_mp2_term_scs * hartree, " eV", &
            E_mp2_scs* hartree ," eV"

       endif
        term=term+1

!   close i state
        enddo
      else

! -------------------------------------------------------------------
!     Spin polarized systems
! -------------------------------------------------------------------

        i_index = 0
        do i_spin=1,n_spin,1
          do i_spin2=1,n_spin,1
             do i_state=i_start_mp2, n_homo(i_spin),1

                E_mp2_term = 0.d0
                E_mp2_term_scs = 0.d0
                do i_state2=i_start_mp2,n_homo(i_spin2) ,1

                  call dgemm('T', 'N', n_unoccupied, n_unoccupied, &
                         n_loc_prodbas, 1.0d0, &
                         ovlp_3KS(:,n_lumo_min:n_states,i_state, &
                                       i_spin), n_loc_prodbas, &
                         ovlp_3KS(:,n_lumo_min:n_states,i_state2, &
                                  i_spin2), n_loc_prodbas, 0.d0, &
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

               do i_virt=n_lumo(i_spin), n_states,1
                do i_virt2=n_lumo(i_spin2),n_states,1
                    if (i_state.eq.i_state2.and.i_virt.eq.i_virt2.and. &
                        i_spin.eq.i_spin2) cycle

                         E_mp2_a=aux_eri(i_virt,i_virt2)

                         E_mp2_b=0.d0


                        if (i_spin.eq.i_spin2) then

                              E_mp2_b=aux_eri(i_virt2,i_virt)

                            E_mp2_a_scs=(E_mp2_a-E_mp2_b)*E_mp2_a*pt_mp2

                              E_mp2_a=(E_mp2_a-E_mp2_b)*E_mp2_a
                        else

                              E_mp2_a_scs=(E_mp2_a)*E_mp2_a*ps_mp2

                              E_mp2_a=(E_mp2_a)*E_mp2_a


                        endif

                        E_mp2_1=  KS_eigenvalue(i_state,i_spin, 1) &
                           +  KS_eigenvalue(i_state2,i_spin2,1)

                        E_mp2_1=  E_mp2_1 - &
                                KS_eigenvalue(i_virt,i_spin,1)

                        E_mp2_1=  E_mp2_1 - &
                                KS_eigenvalue(i_virt2,i_spin2,1)




                       if (abs(E_mp2_1).lt.1e-6)then

                          write(use_unit,'(10X,A)') &
                           "****************************************"
                          write(use_unit,'(10X,2A)') "| Warning :", &
                              "too close to degeneracy"
                          write(use_unit,'(10X,A)') &
                           "****************************************"
                       endif


                       E_mp2_term= E_mp2_a/E_mp2_1+ E_mp2_term
                      E_mp2_term_scs=E_mp2_a_scs/E_mp2_1+ E_mp2_term_scs

!   close b  virtual state 2
                enddo
!   close a virtual state
               enddo
!   end of MPI distribution
              endif
!   close j state2
              enddo
              if(use_mpi) then
                   call MPI_ALLREDUCE(E_mp2_term, E_mp2_mpi, 1, &
                          MPI_DOUBLE_PRECISION, &
                          MPI_SUM, mpi_comm_global, mpierr)

                   E_mp2_term= E_mp2_mpi

                   call MPI_ALLREDUCE(E_mp2_term_scs, E_mp2_mpi_scs, 1, &
                          MPI_DOUBLE_PRECISION, &
                          MPI_SUM, mpi_comm_global, mpierr)

                   E_mp2_term_scs= E_mp2_mpi_scs

              endif
              E_mp2_term = 0.5d0*E_mp2_term
              E_mp2_term_scs = 0.5d0* E_mp2_term_scs

              E_mp2 = E_mp2 + E_mp2_term
              E_mp2_scs = E_mp2_scs + E_mp2_term_scs

              if(myid.eq.0) then
               write(use_unit,'(2X,A,1X,I4,A,1X,F19.8,A,1X,F19.8,A)') &
                   "| Pair state ", term ,"                        :" &
                   , E_mp2_term * hartree, " eV" &
                  , E_mp2* hartree ," eV"
              endif

              if(myid.eq.0) then
               write(use_unit,'(2X,A,1X,I4,A,1X,F19.8,A,1X,F19.8,A)') &
                   "| SCS Pair state ", term ,"                    :" &
                   , E_mp2_term_scs * hartree, " eV" &
                  , E_mp2_scs* hartree ," eV"
              endif

              term=term+1

!   close i state
         enddo

       enddo
      enddo

      endif
      endif

      if(.not.use_hartree_fock) then
       en_total = total_energy - en_xc + eex_energy + e_mp2
      endif




      if (allocated(ovlp_3KS)) then
        deallocate (ovlp_3KS)
      endif
      if (allocated(aux_eri)) then
        deallocate (aux_eri)
      endif

      call cpu_time(rtime)
      tot_time_mp2 = rtime- time_mp2


      if(myid.eq.0) then
       write(use_unit,'(A)')"---------------------------------------------"
       write(use_unit,'(A)') " "
       write(use_unit,'(A)') " "
       write(use_unit,'(A)') " "

       write(use_unit,'(10X,A)')"---------------------------------------------"

       write(use_unit,'(10X,A)') "| MP2 calculation comes to finish ...  "
       write(use_unit,'(A)') " "
       write(use_unit,'(10X,A)')"---------------------------------------------"

!       write(use_unit,'(10X,A,F12.3,A)') &
!          "| Total time for the KS orbital integrals                : ", &
!           tot_time_3KS_intg, " s"

       if (.not.use_scsmp2) then
        write(use_unit,'(10X,A,F12.3,A)') &
          "| Total time for calculating MP2 correction              : ", &
           tot_time_mp2, " s"
       else
        write(use_unit,'(10X,A,F12.3,A)') &
          "| Total time for calculating MP2+SCS-MP2 correction      : ", &
           tot_time_mp2, " s"
       endif

       write(use_unit,'(A)')
       write(use_unit,'(10X,A)')"---------------------------------------------"

       if (.not.use_scsmp2) then
        if(use_hartree_fock) then
          ! For DHDFT energy print out
          if(use_dftpt2) then
              write(use_unit,'(A)')
              write(use_unit,'(A)')"---------------------------------------------"
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                        "| DHDF/DFT Energy                :" &
                , total_energy, " Ha", &
                    total_energy* hartree," eV"

              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                        "| PT2 contribution               :" &
                , E_mp2* dftpt2_Ec_osPT2, " Ha" &
                , (E_mp2* dftpt2_Ec_osPT2) * hartree ," eV"
              write(use_unit,'(A)')"---------------------------------------------"

              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                        "| Total DHDF energy              :" &
                 , (total_energy + (E_mp2* dftpt2_Ec_osPT2)), " Ha" &
                 , (total_energy + (E_mp2* dftpt2_Ec_osPT2))* hartree ," eV"
              write(use_unit,'(A)')"---------------------------------------------"

              post_scf_total_energy = total_energy + E_mp2* dftpt2_Ec_osPT2
          else
              write(use_unit,'(A)')
              write(use_unit,'(A)')"---------------------------------------------"
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                        "| HF Energy                      :" &
                , total_energy, " Ha", &
                    total_energy* hartree," eV"

              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                        "| MP2 correction                 :" &
                , E_mp2, " Ha" &
                , (E_mp2) * hartree ," eV"
              write(use_unit,'(A)')"---------------------------------------------"

              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
              "| Total Energy + MP2 correction  :" &
                 , (total_energy + (E_mp2)), " Ha" &
                 , (total_energy + (E_mp2))* hartree ," eV"
              write(use_unit,'(A)')"---------------------------------------------"

              post_scf_total_energy = total_energy + E_mp2
          endif

       else
         write(use_unit,'(A)')
         write(use_unit,'(A)')"---------------------------------------------"
         write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                  "| DFT Energy                      :" &
               ,total_energy, " Ha", total_energy* hartree," eV"
         write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                  "| DFT XC Energy                   :" &
               ,en_xc, " Ha", en_xc* hartree," eV"
         write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                  "| Exact Echange Energy            :" &
               ,eex_energy, " Ha", eex_energy* hartree," eV"
         write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                  "| MP2 correction                  :" &
                 ,E_mp2, " Ha"  , (E_mp2) * hartree ," eV"

         write(use_unit,'(A)')
         write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                  "| Total (DFT+MP2) Energy          :" &
                 ,en_total, " Ha"  , en_total * hartree ," eV"

         post_scf_total_energy = en_total
       endif

      else
       write(use_unit,'(A)')
       write(use_unit,'(A)')"---------------------------------------------"
       write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                "| HF Energy                          :" &
           , total_energy, " Ha", &
               total_energy* hartree," eV"

       write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                "| MP2 correction                     :" &
           , E_mp2, " Ha" &
           , (E_mp2) * hartree ," eV"
!      write(use_unit,'(A)')"---------------------------------------------"

       write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                "| SCS-MP2 correction                 :" &
           , E_mp2_scs, " Ha" &
           , (E_mp2_scs) * hartree ," eV"
       write(use_unit,'(A)')"---------------------------------------------"

       write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
       "| Total Energy + MP2 correction      :" &
           , (total_energy + (E_mp2)), " Ha" &
           , (total_energy + (E_mp2))* hartree ," eV"
       write(use_unit,'(A)')"---------------------------------------------"
       write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
       "| Total Energy + SCS-MP2 correction  :" &
           , (total_energy + (E_mp2_scs)), " Ha" &
           , (total_energy + (E_mp2_scs))* hartree ," eV"
      write(use_unit,'(A)')"---------------------------------------------"

         post_scf_total_energy = total_energy + E_mp2_scs
      endif

!   end of if myid
      endif
      return

      end subroutine evaluate_en_mp2_correlation_energy
!****** 
