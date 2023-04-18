!****s* FHI-aims/evaluate_mp2_correlation_energy
!  NAME
!   evaluate_mp2_correlation_energy
!  SYNOPSIS 

      subroutine evaluate_mp2_correlation_energy &
       (n_electrons,total_energy,en_xc, &
        occ_numbers,KS_eigenvalue,KS_eigenvector,post_scf_total_energy )

!  PURPOSE
!    Calculation of the total energy perturbation terms up to the second order using Moller-Plesset theory (MP2) to the HF/DFT Energy.
!    A variant, called SCS-MP2, scales the two spin channels separately.
!    Frozen core approximation is available too. 
! USES
!==============================================================================================
!  Formulas behind the program, noted by Igor Ying Zhang, 2015-03-30
!
! For MP2, everything is trivial and well defined:
!
!    E_c[MP2]=1/4\sum_{i,j}^{occ}\sum_{a,b}^{vir}    |(ia||jb)|^2
!                                                 x --------------------------
!                                                      e_i+e_j-e_a-e_b
!            =1/2\sum_{i,j}^{occ}\sum_{a,b}^{vir}    (ia|jb)^2-(ia|jb)(ib|ja)
!                                                 x --------------------------
!                                                      e_i+e_j-e_a-e_b
!  In spin-polarized cases, (capital represents obtials with beta spin) 
!    (i,j;a,b) with (I,J;A,B)
!
!    E_c[MP2]=1/2\sum_{i,j}^{occ}\sum_{a,b}^{vir}    (ia|jb)^2-(ia|jb)(ib|ja)
!                                                 x --------------------------
!                                                      e_i+e_j-e_a-e_b
!            +1/2\sum_{i,J}^{occ}\sum_{a,B}^{vir}    (ia|JB)^2
!                                                 x --------------------------
!                                                      e_i+e_J-e_a-e_B
!            +1/2\sum_{I,J}^{occ}\sum_{A,B}^{vir}    (IA|JB)^2-(IA|JB)(IB|JA)
!                                                 x --------------------------
!                                                      e_I+e_J-e_A-e_B
!            +1/2\sum_{I,j}^{occ}\sum_{A,b}^{vir}    (IA|jb)^2
!                                                 x --------------------------
!                                                      e_I+e_j-e_A-e_b
!
!  Then in spin-unpolarized cases, we have
!    E_c[MP2]=\sum_{i,j}^{occ}\sum_{a,b}^{vir}    (ia|jb)^2-(ia|jb)(ib|ja)
!                                              x --------------------------
!                                                   e_i+e_j-e_a-e_b
!            +\sum_{i,j}^{occ}\sum_{a,b}^{vir}    (ia|jb)^2
!                                              x --------------------------
!                                                   e_i+e_j-e_a-e_b
!            =\sum_{i,j}^{occ}\sum_{a,b}^{vir}    2(ia|jb)^2-(ia|jb)(ib|ja)
!                                              x --------------------------
!                                                   e_i+e_j-e_a-e_b
!
!    To simplify the code, it can be further organized as:
!    E_c[MP2]=2*\sum_{i<j}^{occ}\sum_{a,b}^{vir}   2(ia|jb)^2-(ia|jb)(ib|ja) 
!                                                x -------------------------- 
!                                                     e_i+e_j-e_a-e_b         
!              +\sum_{i=j}^{occ}\sum_{a,b}^{vir}   (ia|jb)^2
!                                                x -------------------------- 
!                                                     e_i+e_j-e_a-e_b         
!
!==============================================================================================

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
      use restart
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
      character(*), parameter :: & 
      &  func = 'evaluate_mp2_correlation_energy'
!     Igor
      real*8, dimension(:,:,:), allocatable :: E_mp2_spectrum, E_mp2_spectrum_mpi
      real*8, dimension(:,:), allocatable :: E_mp2_spectrum_1D_mpi
      real*8 :: E_mp2_direct, E_mp2_exchange
      real*8 :: E_mp2_d, E_mp2_e
      real*8 :: E_mp2_term_d, E_mp2_term_e
      real*8 :: E_mp2_mpi_d, E_mp2_mpi_e

      ! for frozen virtual orbital filter
      logical :: is_fv

      E_mp2_direct = 0.d0
      E_mp2_exchange = 0.d0

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



! Initialize the array of virtual orbitals that need to be frozen
      call read_frozen_virtual_orbitals()
      if (fv_filter) then
          write(info_str,'(2X,A)') &
            '| Read in virtual orbitals that need to be frozen.'
          call localorb_info(info_str,use_unit,'(A)')
      end if

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

!       if(occ_numbers(n_homo(i_spin),i_spin,1).eq.dble(2/n_spin))then
!         n_lumo(i_spin) = n_homo(i_spin) +1
!!       else
!!         n_lumo(i_spin) = n_homo(i_spin)
!!       endif
!      enddo
!!      endif
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


!      if (n_spin.eq.1.and.n_homo_max.eq.1) then
!       E_mp2=0.d0
!       return
!      endif

!     This is the implicit flag to determine whether FC treatment is used!
      if (allocated(n_fc_shell)) then

        ! if we're determining the FC orbitals automatically from last complete
        ! noble gas shell:
        if ((i_start_mp2.eq.0).and.(n_electrons.gt.2)) then

         ! Determine how many empty sites there are for this species ...

!         This version was fundamentally broken.
!         i_empty=0
!         do i_species=1,n_species
!          do i_atom=1,atoms_in_structure(i_species)
!           if (empty(i_atom)) then
!            i_empty(i_species)= i_empty(i_species)+1
!           endif
!          enddo
!         enddo

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
      if(.not. allocated(ovlp_3KS) .or. use_dftpt2_and_hse) then
        !if(use_hartree_fock .and. .not. sparse_o3fn) then
        !   allocate(ovlp_3KS(n_loc_prodbas,n_lumo_min:n_states, &
        !   &                 i_start_mp2:n_homo_max,n_spin), stat=info)
        !   call check_allocation(info, 'ovlp_3KS', func)
        !   n_bytes_o3KS = 8 * n_loc_prodbas * (n_states - n_lumo_min + 1)  &
        !                    * (n_homo_max - i_start_mp2) * n_spin
        !else
           allocate(ovlp_3KS(n_loc_prodbas, n_states, &
           &                 1:n_homo_max, n_spin), stat=info)
           n_bytes_o3KS = 8 * n_loc_prodbas * n_states * n_homo_max * n_spin
        !end if

        write(info_str, "(2X,A,I7,A,I7,A,F12.3,A)") &
        & 'The ovlp_3KS matrix takes another', &
        & nint(n_bytes_o3KS / 2.d0**20), ' MiB x', n_tasks, ' procs =', &
        & n_tasks * n_bytes_o3KS / 2.d0**30, ' GiB.'
        call localorb_info(info_str, use_unit, '(A)', OL_norm)
        call localorb_info('', use_unit, '(A)', OL_norm)

        !if(use_hartree_fock .and. .not. sparse_o3fn) then
        !   call evaluate_ovlp_3MO(KS_eigenvector, &
        !              ovlp_3KS,n_lumo_min,n_lumo)
        !else
           if (sparse_o3fn) then
              call ovlp3KS_lvl_1d(n_homo_max, KS_eigenvector, &
              &                   coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
              deallocate(coulomb_matr_lvl)
              call dealloc_sp_ten(coeff_3fn_ten)
           else
              call transform_ovlp3fn(n_homo_max, KS_eigenvector, ovlp_3fn, &
              &                      ovlp_3KS)
              if (.not. calculate_atom_bsse .and. .not. use_dftpt2_and_lrc) then
                deallocate(ovlp_3fn)
              endif
           end if
           !call evaluate_eex_energy &
           !   ( n_homo_max,n_homo,occ_numbers, &
           !     ovlp_3KS,eex_energy &
           !   )
        !endif
      endif
      call get_times(time_ovlp3fn_to_ovlp3KS, clock_time_ovlp3fn_to_ovlp3KS)
      call evaluate_eex_energy &
         ( n_homo_max,n_homo,occ_numbers, &
           ovlp_3KS,eex_energy &
         )

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
      E_mp2_direct=0.d0
      E_mp2_exchange=0.d0
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

! -------------------------------------------------------------------
!    NON-Spin polarized systems
! -------------------------------------------------------------------

     if (n_spin.eq.1) then

          i_index = 0
          do i_state=i_start_mp2, n_homo(n_spin),1

             E_mp2_term = 0.d0
             do i_state2=i_state,n_homo(n_spin),1
!             do i_state2=1,n_homo(n_spin),1

              if (n_loc_prodbas .gt. 0) then
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
              else
                  aux_eri = 0.d0
              endif

              call sync_matrix( aux_eri(n_lumo_min:n_states, &
                                n_lumo_min:n_states),n_unoccupied, &
                                n_unoccupied &
                               )

              i_index = i_index + 1
! MPI task distribution
              if(myid.eq.MOD(i_index,n_tasks)) then

               do i_virt=n_lumo(n_spin), n_states,1
                call determine_frozen_virtual(i_virt, fv_orbs_n, fv_orbs, fv_filter, is_fv)
                if (is_fv) cycle
                do i_virt2=n_lumo(n_spin),n_states,1
                  call determine_frozen_virtual(i_virt2, fv_orbs_n, fv_orbs, fv_filter, is_fv)
                  if (is_fv) cycle

                  E_mp2_a=aux_eri(i_virt,i_virt2)


!                  if (i_state.eq.i_state2.or.i_virt.eq.i_virt2) then
!                        E_mp2_a= E_mp2_a* E_mp2_a
!                  else
!                     E_mp2_b=0.d0

                    E_mp2_b=aux_eri(i_virt2,i_virt)

                  E_mp2_a= (2.d0*E_mp2_a- &
                             E_mp2_b)* E_mp2_a

!                  endif

                  if (i_state.ne.i_state2) then
                     E_mp2_a=E_mp2_a*2.d0
                  endif


                  E_mp2_1=  (KS_eigenvalue(i_state,n_spin, 1)) &
                    +  (KS_eigenvalue(i_state2,n_spin,1) )

                  E_mp2_1= E_mp2_1-(KS_eigenvalue(i_virt,n_spin,1))

                  E_mp2_1= E_mp2_1-(KS_eigenvalue(i_virt2,n_spin,1))

                  E_mp2_1= E_mp2_1 + en_shift_constant  ! for level shift




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
                E_mp2_term_d= 0.d0
                E_mp2_term_e= 0.d0
                do i_state2=i_start_mp2,n_homo(i_spin2) ,1
                if (n_loc_prodbas .gt. 0) then
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
                else
                    aux_eri = 0.d0
                endif

                  call sync_matrix( aux_eri(n_lumo_min:n_states, &
                                  n_lumo_min:n_states),n_unoccupied, &
                                  n_unoccupied &
                                 )

               i_index = i_index + 1
! MPI task distribution
               if(myid.eq.MOD(i_index,n_tasks)) then

               do i_virt=n_lumo(i_spin), n_states,1
                call determine_frozen_virtual(i_virt, fv_orbs_n, fv_orbs, fv_filter, is_fv)
                if (is_fv) cycle
                do i_virt2=n_lumo(i_spin2),n_states,1
                    call determine_frozen_virtual(i_virt2, fv_orbs_n, fv_orbs, fv_filter, is_fv)
                    if (is_fv) cycle
                    !if (i_state.eq.i_state2.and.i_virt.eq.i_virt2.and. &
                    !    i_spin.eq.i_spin2) cycle

                         E_mp2_a=aux_eri(i_virt,i_virt2)
                         E_mp2_b=0.d0

                         E_mp2_d=aux_eri(i_virt,i_virt2)
                         E_mp2_e=0.d0

                         if (i_spin.eq.i_spin2) then

                              E_mp2_b=aux_eri(i_virt2,i_virt)
                              E_mp2_a=(E_mp2_a-E_mp2_b)*E_mp2_a

                              E_mp2_e=-E_mp2_d*aux_eri(i_virt2,i_virt)
                              E_mp2_d=E_mp2_d*E_mp2_d

                        else
                              E_mp2_a=(E_mp2_a)*E_mp2_a

                              E_mp2_d=E_mp2_d*E_mp2_d
                        endif

                        E_mp2_1=  KS_eigenvalue(i_state,i_spin, 1) &
                           +  KS_eigenvalue(i_state2,i_spin2,1)

                        E_mp2_1=  E_mp2_1 - &
                                KS_eigenvalue(i_virt,i_spin,1)

                        E_mp2_1=  E_mp2_1 - &
                                KS_eigenvalue(i_virt2,i_spin2,1)
                        
                        E_mp2_1= E_mp2_1 + en_shift_constant  ! for level shift


                       if (abs(E_mp2_1).lt.1e-6)then

                          write(use_unit,'(10X,A)') &
                           "****************************************"
                          write(use_unit,'(10X,2A)') "| Warning :", &
                              "too close to degeneracy"
                          write(use_unit,'(10X,A)') &
                           "****************************************"
                       endif


                       E_mp2_term= E_mp2_a/E_mp2_1+ E_mp2_term
                       E_mp2_term_d= E_mp2_d/E_mp2_1+ E_mp2_term_d
                       E_mp2_term_e= E_mp2_e/E_mp2_1+ E_mp2_term_e

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

                   call MPI_ALLREDUCE(E_mp2_term_d, E_mp2_mpi_d, 1, &
                          MPI_DOUBLE_PRECISION, &
                          MPI_SUM, mpi_comm_global, mpierr)

                   E_mp2_term_d= E_mp2_mpi_d

                   call MPI_ALLREDUCE(E_mp2_term_e, E_mp2_mpi_e, 1, &
                          MPI_DOUBLE_PRECISION, &
                          MPI_SUM, mpi_comm_global, mpierr)

                   E_mp2_term_e= E_mp2_mpi_e
              endif
              E_mp2_term   = 0.5d0*E_mp2_term
              E_mp2_term_d = 0.5d0*E_mp2_term_d
              E_mp2_term_e = 0.5d0*E_mp2_term_e

              E_mp2          = E_mp2 + E_mp2_term
              E_mp2_direct   = E_mp2_direct + E_mp2_term_d
              E_mp2_exchange = E_mp2_exchange + E_mp2_term_e

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

       if(myid.eq.0) then
        write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
            "| Direct term                   :" &
            , E_mp2_direct, " Ha" &
           , E_mp2_direct* hartree ," eV"
        write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
            "| Exchange term                 :" &
            , E_mp2_exchange, " Ha" &
           , E_mp2_exchange* hartree ," eV"
       endif

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
                call determine_frozen_virtual(i_virt, fv_orbs_n, fv_orbs, fv_filter, is_fv)
                if (is_fv) cycle
                do i_virt2=n_lumo(n_spin),n_states,1
                  call determine_frozen_virtual(i_virt2, fv_orbs_n, fv_orbs, fv_filter, is_fv)
                  if (is_fv) cycle

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

                  E_mp2_1= E_mp2_1 + en_shift_constant  ! for level shift




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
                call determine_frozen_virtual(i_virt, fv_orbs_n, fv_orbs, fv_filter, is_fv)
                if (is_fv) cycle
                do i_virt2=n_lumo(i_spin2),n_states,1
                    call determine_frozen_virtual(i_virt2, fv_orbs_n, fv_orbs, fv_filter, is_fv)
                    if (is_fv) cycle
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

                        E_mp2_1= E_mp2_1 + en_shift_constant  ! for level shift



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
              write(use_unit,'(A)')"---------------------------------------------"
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                        "| Exact exchange                 :" &
                , eex_energy * dftpt2_Ex_hf, " Ha", &
                    eex_energy * dftpt2_Ex_hf * hartree," eV"

              total_energy = total_energy + eex_energy * dftpt2_Ex_hf

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

      end subroutine evaluate_mp2_correlation_energy
!****** 

!****s* FHI-aims/evaluate_mp2_correlation_energy
!  NAME
!   evaluate_mp2_correlation_energy
!  SYNOPSIS 

      subroutine evaluate_mp2_correlation_energy_2 &
       (n_electrons,total_energy,en_xc, &
        occ_numbers,KS_eigenvalue,KS_eigenvector,post_scf_total_energy)

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
      use restart
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

      real*8 :: E_mp2_a, E_mp2_b
      real*8 :: E_mp2_a_scs
      real*8 :: E_mp2_term, E_mp2_term_scs

      real*8 :: eex_energy
      real*8 :: en_total


      integer :: n_homo(n_spin)
      integer :: n_lumo(n_spin)
      integer :: n_lumo_min
      integer :: n_homo_max
      integer :: min_lumo_loc, max_homo_loc, max_local_states1

      integer ::  n_shell_aux(n_species)

      real*8 :: E_mp2, E_mp2_1, E_mp2_scs

      real*8, allocatable :: tmp_o3KS(:,:)
      real*8, allocatable :: tmp_eri_full(:,:,:,:)
      real*8, allocatable :: tmp_eri_tran(:,:,:)

!     Time counters

!      real*8 time_3KS_intg
      real*8 time_mp2

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
      integer i_species
      integer i_atom
      integer term
      integer i_empty(n_species)
      integer i_start_b
      integer n_p1, n_p2
      integer i_virt_loc, i_virt2_loc, i_virt_own
      integer i_state2_start
      integer i_state2_loc

      ! for frozen virtual orbital filter
      logical :: is_fv

!     Error variable
      integer mpierr
      character*150 :: info_str
      character(*), parameter :: func = 'evaluate_mp2_correlation_energy_2'

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

! Initialize the array of virtual orbitals that need to be frozen
      call read_frozen_virtual_orbitals()
      if (fv_filter) then
          write(info_str,'(2X,A)') &
            '| Read in virtual orbitals that need to be frozen.'
          call localorb_info(info_str,use_unit,'(A)')
      end if

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

      n_lumo_min = minval(n_lumo(:))

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

!      write(use_unit,'(10X,A)')

!     This is the implicit flag to determine whether FC treatment is used!
      if (allocated(n_fc_shell)) then

        ! if we're determining the FC orbitals automatically from last complete
        ! noble gas shell:
        if ((i_start_mp2.eq.0).and.(n_electrons.gt.2)) then

         ! Determine how many empty sites there are for this species ...

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

      n_homo_max = maxval(n_homo(:))

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

        ! Attention: ndim1_o3KS/ndim2_o3KS must be the same everywhere
        ! and ovlp_3KS must not contains undefined (NaN) entries!
        ndim1_o3KS = (n_states-1)/np1_o3KS + 1
        ndim2_o3KS = (n_homo_max-1)/np2_o3KS + 1
        allocate( ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin) )
        ovlp_3KS = 0.
        n_bytes_o3KS = 8 * n_basbas * ndim1_o3KS * ndim2_o3KS * n_spin

        write(info_str, "(2X,A,I7,A,I7,A,F12.3,A)") &
        & 'The ovlp_3KS matrix takes another', nint(n_bytes_o3KS / 2**20), &
        & ' MiB x', n_tasks, ' procs =', &
        & n_tasks * n_bytes_o3KS / 2.d0**30, ' GiB.'
        call localorb_info(info_str, use_unit, '(A)', OL_norm)
        call localorb_info('', use_unit, '(A)', OL_norm)

        if (sparse_o3fn) then
           call ovlp3KS_lvl_2d(n_homo_max, KS_eigenvector, &
           &                   coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
           deallocate(coulomb_matr_lvl)
           call dealloc_sp_ten(coeff_3fn_ten)
        else
           call transform_ovlp3fn_2(n_homo_max,KS_eigenvector, &
           &                        ovlp_3fn, ovlp_3KS)
           if (.not. calculate_atom_bsse) then
             deallocate(ovlp_3fn)
           endif
        end if

        if(.not.use_hartree_fock) then
          call evaluate_eex_energy_2 &
             ( n_homo_max,n_homo,occ_numbers, &
               ovlp_3KS,eex_energy &
             )
        endif
      endif
      call get_timestamps(rtime, clock_rtime)
      time_ovlp3fn_to_ovlp3KS=rtime-time_ovlp3fn_to_ovlp3KS
      clock_time_ovlp3fn_to_ovlp3KS=clock_rtime-clock_time_ovlp3fn_to_ovlp3KS

      call cpu_time(rtime)
!      tot_time_3KS_intg = rtime- time_3KS_intg


      call cpu_time(time_mp2)

      !if(myid.eq.0) then
      !  write(use_unit,'(A)')
      !  write(use_unit,'(A)')"---------------------------------------------"
      !  write(use_unit,'(A)')"| Sum of correlation terms "
      !  write(use_unit,'(A)')"---------------------------------------------"
      !endif


      E_mp2=0.d0
      E_mp2_scs=0.d0
      term=i_start_mp2

!---------------------------------------------------------------------------------------------------

      min_lumo_loc = loc_dim1_o3ks(n_lumo_min) ! overall minimal local value corresponding to n_lumo_min
      max_homo_loc = loc_dim2_o3ks(n_homo_max) ! overall maximal local value corresponding to n_homo_max

      max_local_states1 = ndim1_o3KS - min_lumo_loc + 1

      allocate(tmp_o3KS(n_basbas,min_lumo_loc:ndim1_o3KS))
      allocate(tmp_eri_full(min_lumo_loc:ndim1_o3KS,min_lumo_loc:ndim1_o3KS,0:np1_o3KS-1,max_homo_loc))
      allocate(tmp_eri_tran(min_lumo_loc:ndim1_o3KS,min_lumo_loc:ndim1_o3KS,0:np1_o3KS-1))

      do i_spin=1,n_spin,1
        do i_spin2=i_spin,n_spin,1 ! Only terms which are not covered by symmetry
          do i_state=i_start_mp2, n_homo(i_spin),1

            !if(myid.eq.0) write(use_unit,*) " | i_spin, i_spin2, i_state = ", i_spin, i_spin2, i_state

            E_mp2_term = 0.d0
            E_mp2_term_scs = 0.d0

            if(i_spin.eq.i_spin2) then
              ! Symmetry in tmp_eri matrices may be exploited
              i_state2_start = i_state
            else
              ! No symmetry in tmp_eri (but symmetry of i_spin/i_spin2 taken into account)
              i_state2_start = i_start_mp2
            endif

            n_p2 = own_dim2_o3ks(i_state)

            do n_p1 = 0, np1_o3KS-1

              ! Proc (n_p1,n_p2) broadcasts its local part of
              ! ovlp_3KS_full(:,n_lumo_min:n_states,i_state,i_spin)

              if(n_p1==myp1_o3KS .and. n_p2==myp2_o3KS) &
                tmp_o3KS(:,:) = ovlp_3KS(:,min_lumo_loc:ndim1_o3KS,loc_dim2_o3ks(i_state),i_spin)

              call MPI_Bcast(tmp_o3KS,n_basbas*max_local_states1,MPI_REAL8,global_id(n_p1,n_p2),mpi_comm_global,mpierr)

              ! Multiply what we got with our local part of ovlp_3KS_full(:,n_lumo_min:n_states,x:n_homo(i_spin2),i_spin2)

              do i_state2=i_state2_start,n_homo(i_spin2),1

                if(own_dim2_o3ks(i_state2) /= myp2_o3KS) cycle
                i_state2_loc = loc_dim2_o3ks(i_state2)

                call dgemm('T', 'N', max_local_states1, max_local_states1, n_basbas, 1.d0, &
                           tmp_o3KS, ubound(tmp_o3KS,1), &
                           ovlp_3KS(1,min_lumo_loc,i_state2_loc,i_spin2), ubound(ovlp_3KS,1), &
                           0.d0, tmp_eri_full(min_lumo_loc,min_lumo_loc,n_p1,i_state2_loc),max_local_states1)
              enddo

            enddo

            do i_state2=i_state2_start,n_homo(i_spin2),1

              ! Original code:
              !
              ! n_unoccupied = n_states-n_lumo_min+1
              ! call dgemm('T', 'N', n_unoccupied, n_unoccupied, n_loc_prodbas, 1.0d0, &
              !        ovlp_3KS(:,n_lumo_min:n_states,i_state,i_spin), n_loc_prodbas, &
              !        ovlp_3KS(:,n_lumo_min:n_states,i_state2,i_spin2), n_loc_prodbas, 0.d0, &
              !        aux_eri(n_lumo_min:n_states,n_lumo_min:n_states), n_unoccupied)

              if(own_dim2_o3ks(i_state2) /= myp2_o3KS) cycle
              i_state2_loc = loc_dim2_o3ks(i_state2)

              ! Transpose the strip we have in tmp_eri_full:

              if(i_spin .eq. i_spin2) &
              call MPI_Alltoall(tmp_eri_full(min_lumo_loc,min_lumo_loc,0,i_state2_loc), &
                                max_local_states1*max_local_states1,MPI_REAL8, &
                                tmp_eri_tran,max_local_states1*max_local_states1,MPI_REAL8, &
                                mpi_comm_rows_aux_2d,mpierr)

              do i_virt=n_lumo(i_spin), n_states,1
                call determine_frozen_virtual(i_virt, fv_orbs_n, fv_orbs, fv_filter, is_fv)
                if (is_fv) cycle
                do i_virt2=n_lumo(i_spin2),n_states,1
                  call determine_frozen_virtual(i_virt2, fv_orbs_n, fv_orbs, fv_filter, is_fv)
                  if (is_fv) cycle

                  if(own_dim1_o3ks(i_virt2) /= myp1_o3ks) cycle

                  i_virt_loc  = loc_dim1_o3ks(i_virt)
                  i_virt_own  = own_dim1_o3ks(i_virt)
                  i_virt2_loc = loc_dim1_o3ks(i_virt2)

                  E_mp2_a=tmp_eri_full(i_virt_loc,i_virt2_loc,i_virt_own,i_state2_loc)
                  if(i_spin.eq.i_spin2) then
                    E_mp2_b=tmp_eri_tran(i_virt2_loc,i_virt_loc,i_virt_own)
                  else
                    E_mp2_b=0.d0
                  endif

                  if(n_spin==1) then
                    E_mp2_a_scs=(2.d0*E_mp2_a*(pt_mp2+ps_mp2)-E_mp2_b*pt_mp2)*E_mp2_a
                    E_mp2_a=(2.d0*E_mp2_a-E_mp2_b)*E_mp2_a
                  else
                    if(i_spin.eq.i_spin2) then
                      E_mp2_a_scs=(E_mp2_a-E_mp2_b)*E_mp2_a*pt_mp2
                      E_mp2_a=(E_mp2_a-E_mp2_b)*E_mp2_a
                    else
                      E_mp2_a_scs=(E_mp2_a)*E_mp2_a*ps_mp2
                      E_mp2_a=(E_mp2_a)*E_mp2_a
                    endif
                  endif

                  ! If symmetry is exploited, take symmetric term into account
                  if (i_spin.ne.i_spin2 .or. i_state.ne.i_state2) then
                    E_mp2_a=E_mp2_a*2.d0
                    E_mp2_a_scs=E_mp2_a_scs*2.d0
                  endif

                  E_mp2_1 = KS_eigenvalue(i_state,i_spin,1)   &
                          + KS_eigenvalue(i_state2,i_spin2,1) &
                          - KS_eigenvalue(i_virt,i_spin,1)    &
                          - KS_eigenvalue(i_virt2,i_spin2,1)

                  if (abs(E_mp2_1).lt.1e-6)then
                    write(use_unit,'(10X,A)') "****************************************"
                    write(use_unit,'(10X,A)') "| Warning : too close to degeneracy"
                    write(use_unit,'(10X,A)') "****************************************"
                  endif

                  E_mp2_term = E_mp2_term + E_mp2_a/E_mp2_1
                  E_mp2_term_scs = E_mp2_term_scs + E_mp2_a_scs/E_mp2_1

                enddo ! i_virt2
              enddo ! i_virt
            enddo ! i_state2

            if(use_mpi) then
              call sync_real_number(E_mp2_term)
            endif
            E_mp2_term = E_mp2_term / dble(n_spin) 

            if(use_scsmp2) then
              if(use_mpi) then
                call sync_real_number(E_mp2_term_scs)
              endif
              E_mp2_term_scs = E_mp2_term_scs / dble(n_spin) 
            endif

            E_mp2 = E_mp2 + E_mp2_term
            if(use_scsmp2) E_mp2_scs = E_mp2_scs + E_mp2_term_scs

            !if(myid.eq.0) then
            !  write(use_unit,'(2X,A,1X,I4,A,1X,F19.8,A,1X,F19.8,A)')  &
            !    "| Pair state ", term ,"                    :", &
            !    E_mp2_term * hartree, " eV", E_mp2* hartree, " eV"

            !  if(use_scsmp2) &
            !  write(use_unit,'(2X,A,1X,I4,A,1X,F19.8,A,1X,F19.8,A)')  &
            !    "| SCS Pair state ", term ,"                :", &
            !    E_mp2_term_scs * hartree, " eV", E_mp2_scs* hartree, " eV"
            !endif

            term = term + 1

          enddo ! i_state
        enddo ! i_spin2
      enddo ! i_spin

!---------------------------------------------------------------------------------------------------

      if(.not.use_hartree_fock) then
       en_total = total_energy - en_xc + eex_energy + e_mp2
      endif

      if (allocated(ovlp_3KS)) then
        deallocate (ovlp_3KS)
      endif

      call cpu_time(rtime)
      tot_time_mp2 = rtime- time_mp2


      if(myid.eq.0) then
       write(use_unit,'(A)')"---------------------------------------------"
       write(use_unit,'(A)') " "
       write(use_unit,'(A)') " "
       write(use_unit,'(A)') " "

       write(use_unit,'(10X,A)')"---------------------------------------------"

       write(use_unit,'(10X,A)') "| MP2 calculation comes to stop ...  "
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

              write(use_unit,'(A)')"---------------------------------------------"
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                                        "| Exact exchange                 :" &
                , eex_energy * dftpt2_Ex_hf, " Ha", &
                    eex_energy * dftpt2_Ex_hf * hartree," eV"

              total_energy = total_energy + eex_energy * dftpt2_Ex_hf

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
            post_scf_total_energy = total_energy + E_mp2*dftpt2_Ec_osPT2

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
          endif
          post_scf_total_energy = total_energy + E_mp2

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
!   record the total energy for atom_bsse 
!   and also timings
     if (calculate_atom_bsse .and. use_hartree_fock .and. use_mp2) then
!  if first structure
         if (current_atom_num_bsse==0) then
             BSSE_full_energy= (total_energy + (E_mp2))* hartree
             call get_timestamps(temp_time_bsse, temp_clock_time_bsse )
             time_bsse =temp_time_bsse-time_bsse
             clock_time_bsse =  temp_clock_time_bsse-clock_time_bsse 
             if(myid.eq.0) then
                write(use_unit,'(2X,A,F12.3,A)') &
                "Total time taken for full structure              : ", &
                time_bsse, " s"
             endif
         else       
!  if second structure and so on
             BSSE_per_atom(current_atom_num_bsse)=(total_energy + (E_mp2))* hartree
             call get_timestamps(temp_time_bsse, temp_clock_time_bsse )
             time_bsse =temp_time_bsse-time_bsse
             clock_time_bsse =  temp_clock_time_bsse-clock_time_bsse 
             if(myid.eq.0) then
                write(use_unit,'(2X,A,F12.3,A)') &
                "Total time taken for the last BSSE step              : ", &
                time_bsse, " s"
             endif
         endif
     endif
      return

      end subroutine evaluate_mp2_correlation_energy_2
!****** 

    subroutine determine_frozen_virtual(i_virt, fv_orbs_n, fv_orbs, fv_filter, is_fv)
        integer :: i_virt, fv_orbs_n
        integer :: fv_orbs(fv_orbs_n)
        logical :: fv_filter, is_fv

        integer :: i_fv

        is_fv = .false.
        if (.not. fv_filter) return
        do i_fv = 1, fv_orbs_n
            if (i_virt.eq.fv_orbs(i_fv)) then
                is_fv = .true.
                return
            end if
        end do
    end subroutine determine_frozen_virtual
