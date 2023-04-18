!****s* FHI-aims/evaluate_en_mp2_correlation_energy
!  NAME
!   evaluate_en_mp2_correlation_energy
!  SYNOPSIS 

      subroutine evaluate_iepa_mp2_correlation_energy &
       (n_electrons,total_energy,en_xc, &
        occ_numbers,KS_eigenvalue,KS_eigenvector,post_scf_total_energy )

!  PURPOSE
!    Calculation of the second-order Bethe-Goldstone equation correlation (BGE2),
!    which can also be named as the second-order IEPA in quantum chemistry (IEPA2).
!    Frozen core approximation is available too. 
!
!  AUTHOR
!    Igor Ying Zhang
!
!==============================================================================================
!  Formulas behind the program, noted by Igor Ying Zhang, 2015-03-30
!
! For IEPA2
!
!    E_c[IEPA2]=\sum_{i<j}^{occ}\sum_{a<b}^{vir}         |(ia||jb)|^2
!                                                x --------------------------
!                                                    e_i+e_j-e_a-e_b+e_{ij}
!  For each electron-pair correlation e_{ij}:
!    e_{ij} = \sum_{a<b}^{vir}          |(ia||jb)|^2
!                              x --------------------------
!                                 e_i+e_j-e_a-e_b+e_{ij}
!           = \sum_{a<b}^{vir}      |(ia|jb)-(ib|ja)|^2
!                              x --------------------------
!                                 e_i+e_j-e_a-e_b+e_{ij}
!  Then:
!   E_c[IEPA2]=\sum{i<j}^{occ}e_{ij}
!
!  In spin-polarized cases, (capital represents obtials with beta spin)         
!   (i,j;a,b) with (I,J;A,B) and                                               j
!   consider that the real index of (i,a) is smaller than (I,A):               ->                   (i=1,j=2)    
!    e_{ij} = \sum_{a<b}^{vir}    |(ia|jb)-(ib|ja)|^2                     i | |------------------|------------------|
!                               x --------------------------                V |                  |   ---------------|
!                                  e_i+e_j-e_a-e_b+e_{ij}                     |                  |  b   ------------|
!           = \sum_{a<b}^{vir}    (ia|jb)^2-2(ia|jb)(ib|ja)+(ib|ja)^2         |                  |  ->     ---------|
!                               x -----------------------------------         |                  |a|          ------|
!                                      e_i+e_j-e_a-e_b+e_{ij}                 |                  | V             ---|
!           = \sum_{a,b}^{vir}    (ia|jb)^2-(ia|jb)(ib|ja)                    |------------------|------------------|
!                               x -----------------------------------         |                  |                  |
!                                      e_i+e_j-e_a-e_b+e_{ij}                 |                  |                  |
!    e_{iJ} = \sum_{a,B}^{vir}    |(ia|JB)-(iB|Ja)|^2                         |                  |                  |
!                               x --------------------------                  |                  |                  |
!                                  e_i+e_J-e_a-e_B+e_{iJ}                     |                  |                  |
!           = \sum_{a,B}^{vir}    (ia|JB)^2                                   |------------------|------------------|
!                               x --------------------------
!                                  e_i+e_J-e_a-e_B+e_{iJ}
!    e_{IJ} = \sum_{A<B}^{vir}    |(IA|JB)-(IB|JA)|^2
!                               x --------------------------
!                                  e_I+e_J-e_A-e_B+e_{IJ}
!           = \sum_{A,B}^{vir}    (IA|JB)^2-(IA|JB)(IB|JA)
!                               x --------------------------
!                                  e_I+e_J-e_A-e_B+e_{IJ}
!    Then
!    E_c[IEPA2]=\sum_{i,J}^{occ}e_{iJ} + \sum_{i<j}^{occ}e_{ij} + \sum_{I<J}^{occ}e_{IJ} 
!              = E_c[IEPA,os] + E_c[IEPA,ss]
!    where:
!    E_c[IEPA2,os]=\sum_{i,J}^{occ}e_{iJ}
!    E_c[IEPA2,ss]=\sum_{i<j}^{occ}e_{ij} + \sum_{I<J}^{occ}e_{IJ} 
!
!  Then in spin-unpolarized cases, we have
!    e_{ij,ss} = \sum_{a,b}^{vir}    (ia|jb)^2-(ia|jb)(ib|ja)
!                                 x ---------------------------
!                                    e_i+e_j-e_a-e_b+e_{ij,ss}
!    e_{ij,os} = \sum_{a,b}^{vir}    (ia|jb)^2
!                                  x --------------------------
!                                     e_i+e_j-e_a-e_b+e_{ij,os}
!    Then
!    E_c[IEPA2]= \sum_{i,j}^{occ}e_{ij,os} + 2*\sum_{i<j}^{occ}e_{ij,ss}
!              = E_c[IEPA2,os] + E_c[IEPA2,ss]
!
!    To simpify the code, it can be further organized as:
!    E_c[IEPA2]= E_c[IEPA2,os] + E_c[IEPA2,ss]
!              = \sum_{i,j}^{occ}e_{ij,os} + 2*\sum_{i<j}^{occ}e_{ij,ss}
!              = 2*\sum_{i<j}^{occ} e_{ij,os} + 2*\sum_{i<j}^{occ}e_{ij,ss}
!                + \sum_{i}^{occ}e_{ii,os}
!              = 2*\sum_{i<j}^{occ} (e_{ij,os} + e_{ij,ss})
!                + \sum_{i}^{occ}e_{ii,os}
!
!==============================================================================================
! USES

      use constants
      use dimensions
      use prodbas
      use species_data
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
      real*8 :: E_mp2_a_ss, E_mp2_mpi_ss
      real*8 :: E_mp2_a_os, E_mp2_mpi_os
      real*8 :: E_mp2_term, E_mp2_term_ss, E_mp2_term_os

      real*8 :: eex_energy
      real*8 :: en_total


      integer  :: n_homo(n_spin)
      integer ::  n_lumo(n_spin)
      integer ::  n_lumo_min
      integer ::  n_homo_max

      integer ::  n_shell_aux(n_species)

      character*20 :: filename

      real*8 :: E_mp2, E_mp2_1, E_mp2_os, E_mp2_ss


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
      character(*), parameter :: func = &
                     'evaluate_iepa_mp2_correlation_energy'
!     Igor
      !real*8, dimension(:,:), allocatable :: numerator_a
      !real*8, dimension(:,:), allocatable :: numerator_b
      real*8, dimension(:,:), allocatable :: denominator
      real*8 :: E_iepa
      real*8 :: E_iepa_ss
      real*8 :: E_iepa_os
      real*8 :: E_iepa_term
      real*8 :: E_iepa_term_ss
      real*8 :: E_iepa_term_os
      real*8 :: E_dftpt2_hf, E_dftpt2_pt2, E_dftpt2_iepa
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
        e_mp2    = 0.d0
        e_mp2_ss = 0.d0
        e_mp2_os = 0.d0
        if(myid.eq.0) then
          write(use_unit,'(2X,2A)')  &
          " * Warning: There is no unoccupied states, and hence ", &
          "the MP2 correlation energy is zero."
        endif
        return
      endif

      if(myid.eq.0) then
        write(use_unit,'(A)') &
                "  | Note : Resolution of identity applied "
        if (n_spin.gt.1) then
          write(use_unit,'(A)') "  | Note : Spin-polarized system "
        else
          write(use_unit,'(A)') "  | Note : Non spin-polarized system "
        endif
        if (flag_frozen_core.and.i_start_mp2.eq. &
           0.and.n_electrons.gt.2) then
          write(use_unit,'(A)') "  | Note : Frozen Core approximation  "
        else if (flag_frozen_core.and.n_electrons.gt.2) then
          write(use_unit,'(A)') &
                 "  | Note : Frozen Core approximation from level "
          write(use_unit,*) "               ",   i_start_mp2
        endif
        do i_atom=1, n_atoms,1
          if (empty(i_atom)) then
            write(use_unit,'(A)') "  | Note : Counterpoise correction"
            exit
          endif
        enddo

      write(use_unit,'(A)') &
                   "---------------------------------------------"


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
       write(use_unit,'(A)') &
              "---------------------------------------------"
     endif

     n_homo_max = max(n_homo(1),n_homo(n_spin))

     ! determine the post-SCF frozen core low state by igor
     if (flag_frozen_core_postSCF) then ! count the frozen core states
         call count_frozen_core_states(i_start_mp2)
         if(myid.eq.0) then
            write(use_unit,'(2X,A,I12)') &
            "Frozen core number for post-SCF calculation          :", &
            i_start_mp2 - 1
         endif
     endif

!     check on frozen core compatibility
     if (i_start_mp2.gt.n_homo_max) then
        write(use_unit,'(1X,A)')&
               "*--------------------------------------------"
        write(use_unit,'(1X,A)')&
               "* Error: Requested more frozen core levels"
        write(use_unit,'(1X,A)')&
               "* than occupied states. Aborting MP2 part."
        write(use_unit,'(1X,A)')&
               "*--------------------------------------------"
        stop
     endif

!      call cpu_time(time_3KS_intg)

      call get_timestamps(time_ovlp3fn_to_ovlp3KS, clock_time_ovlp3fn_to_ovlp3KS)
      if(.not. allocated(ovlp_3KS)) then
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

        if(use_hartree_fock .and. .not. sparse_o3fn) then
          call evaluate_ovlp_3MO(KS_eigenvector, &
                     ovlp_3KS,n_lumo_min,n_lumo)
        else
           if (sparse_o3fn) then
              call ovlp3KS_lvl_1d(n_homo_max, KS_eigenvector, &
              &                   coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
              deallocate(coulomb_matr_lvl)
              call dealloc_sp_ten(coeff_3fn_ten)
           else
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
        write(use_unit,'(A)')&
            "---------------------------------------------"
        write(use_unit,'(A)')"| Sum of correlation terms "
        write(use_unit,'(A)')&
            "---------------------------------------------"
      endif

      E_mp2     = 0.d0
      E_mp2_ss  = 0.d0
      E_mp2_os  = 0.d0
      E_iepa    = 0.d0
      E_iepa_ss = 0.d0
      E_iepa_os = 0.d0
      term=i_start_mp2

      if (.not.allocated(aux_eri)) then
          allocate(aux_eri(n_lumo_min:n_states,n_lumo_min:n_states))
      endif
      if (.not.allocated(denominator)) then
          allocate(denominator(n_lumo_min:n_states,n_lumo_min:n_states))
      endif

! -------------------------------------------------------------------
!    NON-Spin polarized systems
! -------------------------------------------------------------------

     if (n_spin.eq.1) then
          !i_index = 0

          do i_state=i_start_mp2, n_homo(n_spin),1
            do i_state2=i_state,n_homo(n_spin),1

              denominator     = 0.d0
              E_mp2_term      = 0.d0
              E_mp2_term_ss   = 0.d0
              E_mp2_term_os   = 0.d0
              E_iepa_term     = 0.d0
              E_iepa_term_ss  = 0.d0
              E_iepa_term_os  = 0.d0

              !i_index         = i_index + 1

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

              ! Now calculate e_{ij,os} and e_{ij,ss} according to the formulas in  
              ! ==============================================================
              !  In spin-unpolarized cases:
              !    e_{ij,ss} = \sum_{a,b}^{vir}    (ia|jb)^2-(ia|jb)(ib|ja)
              !                                 x ---------------------------
              !                                    e_i+e_j-e_a-e_b+e_{ij,ss}
              !    e_{ij,os} = \sum_{a,b}^{vir}    (ia|jb)^2
              !                                  x --------------------------
              !                                     e_i+e_j-e_a-e_b+e_{ij,os}
              ! ==============================================================
              ! the begining of this subroutine
              i_index = 0
              do i_virt=n_lumo(n_spin), n_states,1
                do i_virt2=n_lumo(n_spin),n_states,1
                   i_index = i_index + 1
                   ! MPI task distribution
                   if(myid.eq.MOD(i_index,n_tasks)) then
                      E_mp2_a=aux_eri(i_virt,i_virt2)
                      E_mp2_b=aux_eri(i_virt2,i_virt)

                      ! separate parallel-spin and opposite-spin contributions

                      E_mp2_a_ss = (E_mp2_a-E_mp2_b)* E_mp2_a
                      E_mp2_a_os = E_mp2_a * E_mp2_a 

                      E_mp2_1= (KS_eigenvalue(i_virt,n_spin, 1)) &
                        +  (KS_eigenvalue(i_virt2,n_spin,1) )

                      E_mp2_1= E_mp2_1 - KS_eigenvalue(i_state,n_spin,1) &
                        - KS_eigenvalue(i_state2,n_spin,1)

                      if (abs(E_mp2_1).lt.1e-6)then
                         write(use_unit,'(10X,A)') &
                              "****************************************"
                         write(use_unit,'(10X,2A)') "| Warning :", &
                           " too close to degeneracy"
                         write(use_unit,'(10X,A)') &
                              "****************************************"
                      endif

                      denominator(i_virt,i_virt2) = E_mp2_1

                      E_mp2_1= 1.d0 / E_mp2_1

                      E_mp2_term    = (E_mp2_a_ss+E_mp2_a_os)*E_mp2_1 + E_mp2_term
                      E_mp2_term_ss = E_mp2_a_ss*E_mp2_1 + E_mp2_term_ss
                      E_mp2_term_os = E_mp2_a_os*E_mp2_1 + E_mp2_term_os
!   end of MPI distribution
                  endif
!   close b  virtual state 2
                enddo
!   close a virtual state
              enddo
              if(use_mpi) then
                  call MPI_ALLREDUCE(E_mp2_term, E_mp2_mpi, 1, &
                         MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpi_comm_global, mpierr)
                  E_mp2_term= E_mp2_mpi
                  call MPI_ALLREDUCE(E_mp2_term_ss, E_mp2_mpi_ss, 1, &
                         MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpi_comm_global, mpierr)
                  E_mp2_term_ss = E_mp2_mpi_ss
                  call MPI_ALLREDUCE(E_mp2_term_os, E_mp2_mpi_os, 1, &
                         MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpi_comm_global, mpierr)
                  E_mp2_term_os = E_mp2_mpi_os
                  call sync_matrix( denominator(n_lumo_min:n_states, &
                                    n_lumo_min:n_states),n_unoccupied, &
                                    n_unoccupied &
                                   )
              endif

              ! Now do e_{ij,os} and e_{ij,ss} iterations
              call nsolver_iepa_mp2(aux_eri,denominator, &
                  1, 1, i_state, i_state2, n_lumo,&
                  E_mp2_term, E_mp2_term_ss, E_mp2_term_os, &
                  n_lumo_min, coupling_pt2_factor, &
                  coupling_pt2_screen, coupling_pt2_shift, &
                  1.0d-7, E_iepa_term, E_iepa_term_ss, E_iepa_term_os)

              if (i_state.ne.i_state2) then
                  E_mp2_term     = E_mp2_term*2.0d0
                  E_mp2_term_os  = E_mp2_term_os*2.0d0
                  E_mp2_term_ss  = E_mp2_term_ss*2.0d0
                  E_iepa_term    = E_iepa_term*2.0d0
                  E_iepa_term_os = E_iepa_term_os*2.0d0
                  E_iepa_term_ss = E_iepa_term_ss*2.0d0
              endif
              E_mp2    = E_mp2    - E_mp2_term
              E_mp2_ss = E_mp2_ss - E_mp2_term_ss
              E_mp2_os = E_mp2_os - E_mp2_term_os

              E_iepa    = E_iepa - E_iepa_term
              E_iepa_ss = E_iepa_ss - E_iepa_term_ss
              E_iepa_os = E_iepa_os - E_iepa_term_os

            enddo !   close j state2
          enddo !   close i state

     else

! -------------------------------------------------------------------
!     Spin polarized systems
! -------------------------------------------------------------------

        !i_index = 0
        do i_spin=1,n_spin,1
          do i_spin2=i_spin,n_spin,1
            if (i_spin.eq.i_spin2) then
            !=======================================================================================
            !    e_{ij} = \sum_{a<b}^{vir}    |(ia|jb)-(ib|ja)|^2                  
            !                               x --------------------------           
            !                                  e_i+e_j-e_a-e_b+e_{ij}              
            !    e_{IJ} = \sum_{A<B}^{vir}    |(IA|JB)-(IB|JA)|^2
            !                               x --------------------------
            !                                  e_I+e_J-e_A-e_B+e_{IJ}
            !    E_c[IEPA2,ss]=\sum_{i<j}^{occ}e_{ij} + \sum_{I<J}^{occ}e_{IJ} 
            !=======================================================================================
              do i_state=i_start_mp2, n_homo(i_spin),1
                do i_state2=i_state+1,n_homo(i_spin2) ,1
                  denominator     = 0.d0
                  E_mp2_term      = 0.d0
                  E_mp2_term_ss   = 0.d0
                  E_iepa_term     = 0.d0
                  E_iepa_term_ss  = 0.d0
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
                  i_index = 0
                  do i_virt=n_lumo(i_spin), n_states,1
                    do i_virt2=i_virt+1,n_states,1
                      i_index  = i_index + 1
                      if(myid.eq.MOD(i_index,n_tasks)) then
                        E_mp2_a=aux_eri(i_virt,i_virt2)
                        E_mp2_b=aux_eri(i_virt2,i_virt)
                        E_mp2_1=  KS_eigenvalue(i_virt,i_spin, 1) &
                           +  KS_eigenvalue(i_virt2,i_spin2,1)
                        E_mp2_1=  E_mp2_1 - KS_eigenvalue(i_state,i_spin,1) &
                           -  KS_eigenvalue(i_state2,i_spin2,1)
                        if (abs(E_mp2_1).lt.1e-6)then
                           write(use_unit,'(10X,A)') &
                            "****************************************"
                           write(use_unit,'(10X,2A)') "| Warning :", &
                               "too close to degeneracy"
                           write(use_unit,'(10X,A)') &
                            "****************************************"
                        endif
                        denominator(i_virt,i_virt2) = E_mp2_1
                        E_mp2_a       = (E_mp2_a-E_mp2_b)**2
                        E_mp2_term    = E_mp2_a/E_mp2_1+ E_mp2_term
                        E_mp2_term_ss = E_mp2_a/E_mp2_1+ E_mp2_term_ss
                      endif !   end of MPI distribution
                    enddo !   close b  virtual state 2
                  enddo !   close a virtual state
                  if(use_mpi) then
                      call MPI_ALLREDUCE(E_mp2_term, E_mp2_mpi, 1, &
                             MPI_DOUBLE_PRECISION, &
                             MPI_SUM, mpi_comm_global, mpierr)
                      E_mp2_term= E_mp2_mpi
                      call MPI_ALLREDUCE(E_mp2_term_ss, E_mp2_mpi_ss, 1, &
                             MPI_DOUBLE_PRECISION, &
                             MPI_SUM, mpi_comm_global, mpierr)
                      E_mp2_term_ss = E_mp2_mpi_ss
                      call sync_matrix( denominator(n_lumo_min:n_states, &
                                        n_lumo_min:n_states),n_unoccupied, &
                                        n_unoccupied &
                                       )
                  endif

                  E_mp2    = E_mp2 - E_mp2_term
                  E_mp2_ss = E_mp2_ss - E_mp2_term_ss

                  call nsolver_iepa_mp2(aux_eri,denominator, &
                      i_spin, i_spin2, i_state, i_state2, n_lumo,&
                      E_mp2_term, E_mp2_term_ss, E_mp2_term_os, &
                      n_lumo_min, coupling_pt2_factor, &
                      coupling_pt2_screen, coupling_pt2_shift, &
                      1.0d-7, E_iepa_term, E_iepa_term_ss, E_iepa_term_os)

                  E_iepa    = E_iepa - E_iepa_term
                  E_iepa_ss = E_iepa_ss - E_iepa_term_ss

                enddo ! i_state
              enddo ! i_state2

            else ! if (i_spin.ne.i_spin2)
            !===========================================================
            !  In spin-polarized cases:
            !    e_{iJ} = \sum_{a,B}^{vir}    |(ia|JB)|^2
            !                               x --------------------------
            !                                  e_i+e_J-e_a-e_B+e_{iJ}   
            !    E_c[IEPA2,os]=\sum_{i,J}^{occ}e_{iJ}
            !===========================================================
              do i_state=i_start_mp2, n_homo(i_spin),1
                do i_state2=i_start_mp2,n_homo(i_spin2) ,1

                  denominator     = 0.d0
                  E_mp2_term      = 0.d0
                  E_mp2_term_os   = 0.d0
                  E_iepa_term     = 0.d0
                  E_iepa_term_os  = 0.d0

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

                  i_index = 0
                  do i_virt=n_lumo(i_spin), n_states,1
                    do i_virt2=n_lumo(i_spin2),n_states,1
                      i_index = i_index + 1
                      ! MPI task distribution
                      if(myid.eq.MOD(i_index,n_tasks)) then
                        E_mp2_a=aux_eri(i_virt,i_virt2)
                        E_mp2_1=  KS_eigenvalue(i_virt,i_spin, 1) &
                           +  KS_eigenvalue(i_virt2,i_spin2,1)
                        E_mp2_1=  E_mp2_1 - KS_eigenvalue(i_state,i_spin,1) &
                           -  KS_eigenvalue(i_state2,i_spin2,1)

                        if (abs(E_mp2_1).lt.1e-6)then
                           write(use_unit,'(10X,A)') &
                            "****************************************"
                           write(use_unit,'(10X,2A)') "| Warning :", &
                               "too close to degeneracy"
                           write(use_unit,'(10X,A)') &
                            "****************************************"
                        endif

                        denominator(i_virt,i_virt2) = E_mp2_1

                        E_mp2_a       = E_mp2_a*E_mp2_a
                        E_mp2_term    = E_mp2_a/E_mp2_1+ E_mp2_term
                        E_mp2_term_os = E_mp2_a/E_mp2_1+ E_mp2_term_os

                      endif !   end of MPI distribution
                    enddo !   close b  virtual state 2
                  enddo !   close a virtual state
                  if(use_mpi) then
                      call MPI_ALLREDUCE(E_mp2_term, E_mp2_mpi, 1, &
                             MPI_DOUBLE_PRECISION, &
                             MPI_SUM, mpi_comm_global, mpierr)
                      E_mp2_term= E_mp2_mpi
                      call MPI_ALLREDUCE(E_mp2_term_os, E_mp2_mpi_os, 1, &
                             MPI_DOUBLE_PRECISION, &
                             MPI_SUM, mpi_comm_global, mpierr)
                      E_mp2_term_os = E_mp2_mpi_os
                      call sync_matrix( denominator(n_lumo_min:n_states, &
                                        n_lumo_min:n_states),n_unoccupied, &
                                        n_unoccupied &
                                       )
                  endif

                  E_mp2    = E_mp2 - E_mp2_term
                  E_mp2_os = E_mp2_os - E_mp2_term_os

                  call nsolver_iepa_mp2(aux_eri,denominator, &
                      i_spin, i_spin2, i_state, i_state2, n_lumo,&
                      E_mp2_term, E_mp2_term_ss, E_mp2_term_os, &
                      n_lumo_min, coupling_pt2_factor, &
                      coupling_pt2_screen, coupling_pt2_shift, &
                      1.0d-7, E_iepa_term, E_iepa_term_ss, E_iepa_term_os)

                  E_iepa    = E_iepa - E_iepa_term
                  E_iepa_os = E_iepa_os - E_iepa_term_os
                enddo ! i_state
              enddo ! i_state2
            endif ! if (i_spin .eq. i_spin2) and (i_spin .ne. i_spin2)
          enddo !i_spin
        enddo !i_spin2

      endif ! if (n_spin .eq. 1) or (n_spin .ne. 1) i.e. close shell or open shell

      !if(.not.use_hartree_fock) then
      en_total = total_energy - en_xc + eex_energy + e_mp2
      !endif

      if (allocated(ovlp_3KS)) then
        deallocate (ovlp_3KS)
      endif
      if (allocated(aux_eri)) then
        deallocate (aux_eri)
      endif
      if (allocated(denominator)) then
        deallocate(denominator)
      endif

      call cpu_time(rtime)
      tot_time_mp2 = rtime- time_mp2

      if(myid.eq.0) then
          write(use_unit,'(A)')&
                 "---------------------------------------------"
          write(use_unit,'(A)') " "
          write(use_unit,'(A)') " "
          write(use_unit,'(A)') " "

          write(use_unit,'(10X,A)')&
                 "---------------------------------------------------"

          write(use_unit,'(10X,A)') &
                 "| sBGE2 calculation comes to finish ...  "
          write(use_unit,'(A)') " "
          write(use_unit,'(10X,A)')&
                 "---------------------------------------------------"

           write(use_unit,'(10X,A,F12.3,A)') &
        "| Total time for calculating sBGE2 correction              : ", &
              tot_time_mp2, " s"

          if(use_dftpt2) then
             E_dftpt2_hf = eex_energy*dftpt2_Ex_hf
             E_dftpt2_pt2 = E_mp2_ss*dftpt2_Ec_ssPT2+E_mp2_os*dftpt2_Ec_osPT2
             E_dftpt2_iepa = E_iepa_ss*dftpt2_Ec_ssPT2+E_iepa_os*dftpt2_Ec_osPT2
             write(use_unit,'(A)') " "
             write(use_unit,'(A)')&
         "---------------------------------------------"
             write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
        "| DHDF/DFT Energy                :" &
                   ,total_energy," Ha"  &
                   ,total_energy* hartree," eV"
             write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
        "| Exx contribution               :" &
                   ,E_dftpt2_hf, " Ha" &
                   ,E_dftpt2_hf*hartree," eV"
             write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
        "| PT2 contribution               :" &
               , E_dftpt2_pt2, " Ha" &
               , E_dftpt2_pt2 * hartree ," eV"
             write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
        "| sBGE2 contribution             :" &
               , E_dftpt2_iepa, " Ha" &
               , E_dftpt2_iepa * hartree ," eV"
            write(use_unit,'(A)')&
        "---------------------------------------------"

            if (flag_dftpt2_dft_part .eq. 33) then
                write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
        "| Total ZRPS energy              :" &
                   , (total_energy + E_dftpt2_hf + E_dftpt2_iepa), " Ha" &
                   , (total_energy + E_dftpt2_hf + E_dftpt2_iepa)* hartree ," eV"
                write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
        "| Total ZRPS(DH) energy          :" &
                   , (total_energy + E_dftpt2_hf + E_dftpt2_pt2), " Ha" &
                   , (total_energy + E_dftpt2_hf + E_dftpt2_pt2)* hartree ," eV"
                write(use_unit,'(A)')&
        "---------------------------------------------"
            endif
            post_scf_total_energy = total_energy + E_mp2*dftpt2_Ec_osPT2
        else
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
            write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
       "| MP2 correction (ss)             :" &
                    ,E_mp2_ss, " Ha"  , (E_mp2_ss) * hartree ," eV"
            write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
       "| MP2 correction (os)             :" &
                    ,E_mp2_os, " Ha"  , (E_mp2_os) * hartree ," eV"
            write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
       "| Total (DFT+MP2) Energy          :" &
                    ,en_total, " Ha"  , en_total * hartree ," eV"

            write(use_unit,'(A)')
            post_scf_total_energy = en_total
            if (en_shift_type.eq.3) then
                write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
       "| SCPT2 contribution              :" &
                  , E_iepa, " Ha" , (E_iepa) * hartree ," eV"
                write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
       "| SCPT2 contribution (ss)         :" &
                  , E_iepa_ss, " Ha" , (E_iepa_ss) * hartree ," eV"
                write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
       "| SCPT2 contribution (os)         :" &
                  , E_iepa_os, " Ha" , (E_iepa_os) * hartree ," eV"
                write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
       "| Total SCPT2 energy              :" &
                   , (en_total - E_mp2 + (E_iepa)), " Ha" &
                   , (en_total - E_mp2 + (E_iepa))* hartree ," eV"
                write(use_unit,'(A)')&
       "---------------------------------------------"
                post_scf_total_energy = en_total - E_mp2 + E_iepa
            endif
        endif
!   end of if myid
      endif

      return

      end subroutine evaluate_iepa_mp2_correlation_energy
!****** 
