!****s* FHI-aims/evaluate_dftpt2_correlation
!  NAME
!   evaluate_dftpt2_correlation
!  SYNOPSIS 

      subroutine evaluate_dftpt2_correlation &
       (n_electrons,total_energy,en_xc, &
        occ_numbers,KS_eigenvalue,KS_eigenvector,post_scf_total_energy )

!  PURPOSE
!    Calculation of the second-order PT2 or lrc-PT2 for double-hybrid DFAs.
!    Frozen core approximation is available too. 
!
!  AUTHOR
!    Igor Ying Zhang
!
!==============================================================================================
!  Formulas behind the program, noted by Igor Ying Zhang, 2015-03-30
!
! For PT2
!
!    E_c[PT2]=\sum_{i<j}^{occ}\sum_{a<b}^{vir}         |(ia||jb)|^2
!                                                x --------------------------
!                                                    e_i+e_j-e_a-e_b
!  For each electron-pair correlation e_{ij}:
!    e_{ij} = \sum_{a<b}^{vir}          |(ia||jb)|^2
!                              x --------------------------
!                                 e_i+e_j-e_a-e_b
!           = \sum_{a<b}^{vir}      |(ia|jb)-(ib|ja)|^2
!                              x --------------------------
!                                 e_i+e_j-e_a-e_b
!  Then:
!   E_c[PT2]=\sum{i<j}^{occ}e_{ij}
!
!  In spin-polarized cases, (capital represents obtials with beta spin)         
!   (i,j;a,b) with (I,J;A,B) and                                               j
!   consider that the real index of (i,a) is smaller than (I,A):               ->                   (i=1,j=2)    
!    e_{ij} = \sum_{a<b}^{vir}    |(ia|jb)-(ib|ja)|^2                     i | |------------------|------------------|
!                               x --------------------------                V |                  |   ---------------|
!                                  e_i+e_j-e_a-e_b                            |                  |  b   ------------|
!           = \sum_{a<b}^{vir}    (ia|jb)^2-2(ia|jb)(ib|ja)+(ib|ja)^2         |                  |  ->     ---------|
!                               x -----------------------------------         |                  |a|          ------|
!                                      e_i+e_j-e_a-e_b                        |                  | V             ---|
!           = \sum_{a,b}^{vir}    (ia|jb)^2-(ia|jb)(ib|ja)                    |------------------|------------------|
!                               x -----------------------------------         |                  |                  |
!                                      e_i+e_j-e_a-e_b                        |                  |                  |
!    e_{iJ} = \sum_{a,B}^{vir}    |(ia|JB)-(iB|Ja)|^2                         |                  |                  |
!                               x --------------------------                  |                  |                  |
!                                  e_i+e_J-e_a-e_B                            |                  |                  |
!           = \sum_{a,B}^{vir}    (ia|JB)^2                                   |------------------|------------------|
!                               x --------------------------
!                                  e_i+e_J-e_a-e_B       
!    e_{IJ} = \sum_{A<B}^{vir}    |(IA|JB)-(IB|JA)|^2
!                               x --------------------------
!                                  e_I+e_J-e_A-e_B
!           = \sum_{A,B}^{vir}    (IA|JB)^2-(IA|JB)(IB|JA)
!                               x --------------------------
!                                  e_I+e_J-e_A-e_B
!    Then
!    E_c[PT2]=\sum_{i,J}^{occ}e_{iJ} + \sum_{i<j}^{occ}e_{ij} + \sum_{I<J}^{occ}e_{IJ} 
!              = E_c[PT,os] + E_c[IEPA,ss]
!    where:
!    E_c[PT2,os]=\sum_{i,J}^{occ}e_{iJ}
!    E_c[PT2,ss]=\sum_{i<j}^{occ}e_{ij} + \sum_{I<J}^{occ}e_{IJ} 
!
!  Then in spin-unpolarized cases, we have
!    e_{ij,ss} = \sum_{a,b}^{vir}    (ia|jb)^2-(ia|jb)(ib|ja)
!                                 x ---------------------------
!                                    e_i+e_j-e_a-e_b+e_{ij,ss}
!    e_{ij,os} = \sum_{a,b}^{vir}    (ia|jb)^2
!                                  x --------------------------
!                                     e_i+e_j-e_a-e_b+e_{ij,os}
!    Then
!    E_c[PT2]= \sum_{i,j}^{occ}e_{ij,os} + 2*\sum_{i<j}^{occ}e_{ij,ss}
!              = E_c[PT2,os] + E_c[PT2,ss]
!
!    To simpify the code, it can be further organized as:
!    E_c[PT2]= E_c[PT2,os] + E_c[PT2,ss]
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
! o  total_energy -- The DFT-part energy 
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
      character(*), parameter :: func = 'evaluate_dftpt2_correlation'
!     Igor
      !real*8, dimension(:,:), allocatable :: numerator_a
      !real*8, dimension(:,:), allocatable :: numerator_b
      !real*8, dimension(:,:), allocatable :: denominator
      !real*8 :: E_iepa_term
      !real*8 :: E_iepa_term_ss
      !real*8 :: E_iepa_term_os
      !real*8 :: E_dftpt2_iepa
      real*8 :: E_dftpt2_hf, E_dftpt2_pt2
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
!   o total_energy -- energy of the dft part of the double-hybrid total energy
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
!    Igor Ying Zhang, FHI
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

     if (.not. lrc_pt2_started) then
       if(myid.eq.0) then
         write(use_unit,'(A)') "  | Note : Resolution of identity applied "
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
       endif !  endif myid
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
     endif !(.not. lrc_pt2_started)

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
      if(allocated(ovlp_3KS)) deallocate(ovlp_3KS)
      allocate(ovlp_3KS(n_loc_prodbas, n_states, &
      &                 1:n_homo_max, n_spin), stat=info)
      n_bytes_o3KS = 8 * n_loc_prodbas * n_states * n_homo_max * n_spin

      write(info_str, "(2X,A,I7,A,I7,A,F12.3,A)") &
      & 'The ovlp_3KS matrix takes another', &
      & nint(n_bytes_o3KS / 2.d0**20), ' MiB x', n_tasks, ' procs =', &
      & n_tasks * n_bytes_o3KS / 2.d0**30, ' GiB.'
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      call localorb_info('', use_unit, '(A)', OL_norm)

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
      endif

      if (.not. lrc_pt2_started) then
        call evaluate_eex_energy(n_homo_max,n_homo,occ_numbers, &
             ovlp_3KS,eex_energy)
      endif

      call get_times(time_ovlp3fn_to_ovlp3KS, clock_time_ovlp3fn_to_ovlp3KS)

      call cpu_time(time_mp2)
      n_unoccupied = n_states-n_lumo_min+1

      !if(myid.eq.0) then
      !  write(use_unit,'(A)')
      !  write(use_unit,'(A)')"---------------------------------------------"
      !  write(use_unit,'(A)')"| Sum of correlation terms "
      !  write(use_unit,'(A)')"---------------------------------------------"
      !endif

      E_mp2     = 0.d0
      E_mp2_ss  = 0.d0
      E_mp2_os  = 0.d0
      term=i_start_mp2

      if (allocated(aux_eri)) deallocate(aux_eri)
      allocate(aux_eri(n_lumo_min:n_states,n_lumo_min:n_states))
      !if (.not.allocated(denominator)) then
      !    allocate(denominator(n_lumo_min:n_states,n_lumo_min:n_states))
      !endif

! -------------------------------------------------------------------
!    NON-Spin polarized systems
! -------------------------------------------------------------------

     if (n_spin.eq.1) then

          do i_state=i_start_mp2, n_homo(n_spin),1
            do i_state2=i_state,n_homo(n_spin),1

              E_mp2_term      = 0.d0
              E_mp2_term_ss   = 0.d0
              E_mp2_term_os   = 0.d0

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

                      E_mp2_1= 1.d0 / E_mp2_1

                      E_mp2_term    = (E_mp2_a_ss+E_mp2_a_os)*E_mp2_1 + E_mp2_term
                      E_mp2_term_ss = E_mp2_a_ss*E_mp2_1 + E_mp2_term_ss
                      E_mp2_term_os = E_mp2_a_os*E_mp2_1 + E_mp2_term_os

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
                  call MPI_ALLREDUCE(E_mp2_term_os, E_mp2_mpi_os, 1, &
                         MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpi_comm_global, mpierr)
                  E_mp2_term_os = E_mp2_mpi_os
              endif

              if (i_state.ne.i_state2) then
                  E_mp2_term     = E_mp2_term*2.0d0
                  E_mp2_term_os  = E_mp2_term_os*2.0d0
                  E_mp2_term_ss  = E_mp2_term_ss*2.0d0
              endif
              E_mp2    = E_mp2    - E_mp2_term
              E_mp2_ss = E_mp2_ss - E_mp2_term_ss
              E_mp2_os = E_mp2_os - E_mp2_term_os

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
            !                                  e_i+e_j-e_a-e_b
            !    e_{IJ} = \sum_{A<B}^{vir}    |(IA|JB)-(IB|JA)|^2
            !                               x --------------------------
            !                                  e_I+e_J-e_A-e_B
            !    E_c[PT2,ss]=\sum_{i<j}^{occ}e_{ij} + \sum_{I<J}^{occ}e_{IJ} 
            !=======================================================================================
              do i_state=i_start_mp2, n_homo(i_spin),1
                do i_state2=i_state+1,n_homo(i_spin2) ,1
                  E_mp2_term      = 0.d0
                  E_mp2_term_ss   = 0.d0
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
                  endif

                  E_mp2    = E_mp2 - E_mp2_term
                  E_mp2_ss = E_mp2_ss - E_mp2_term_ss

                enddo ! i_state
              enddo ! i_state2

            else ! if (i_spin.ne.i_spin2)
            !===========================================================
            !  In spin-polarized cases:
            !    e_{iJ} = \sum_{a,B}^{vir}    |(ia|JB)|^2
            !                               x --------------------------
            !                                  e_i+e_J-e_a-e_B
            !    E_c[PT2,os]=\sum_{i,J}^{occ}e_{iJ}
            !===========================================================
              do i_state=i_start_mp2, n_homo(i_spin),1
                do i_state2=i_start_mp2,n_homo(i_spin2) ,1

                  E_mp2_term      = 0.d0
                  E_mp2_term_os   = 0.d0

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
                  endif

                  E_mp2    = E_mp2 - E_mp2_term
                  E_mp2_os = E_mp2_os - E_mp2_term_os

                enddo ! i_state
              enddo ! i_state2
            endif ! if (i_spin .eq. i_spin2) and (i_spin .ne. i_spin2)
          enddo !i_spin
        enddo !i_spin2

      endif ! if (n_spin .eq. 1) or (n_spin .ne. 1) i.e. close shell or open shell

      if (allocated(ovlp_3KS)) deallocate (ovlp_3KS)
      if (allocated(aux_eri)) deallocate (aux_eri)

      call cpu_time(rtime)
      tot_time_mp2 = rtime- time_mp2

      if (.not. lrc_pt2_started) then
        E_dftpt2_hf = eex_energy*dftpt2_Ex_hf
        E_dftpt2_pt2 = E_mp2_ss*dftpt2_Ec_ssPT2+E_mp2_os*dftpt2_Ec_osPT2
        post_scf_total_energy = total_energy + E_dftpt2_hf+E_dftpt2_pt2
        if(myid.eq.0) then
          write(use_unit,'(A)')&
              "---------------------------------------------"
          write(use_unit,'(A)') " "
          write(use_unit,'(A)') " "
          write(use_unit,'(2X,A)')&
              "---------------------------------------------------"
          write(use_unit,'(2X,A)') &
              "| PT2 calculation comes to finish ...  "
          write(use_unit,'(2X,A)')&
              "---------------------------------------------------"

          write(use_unit,'(2X,A,F12.3,A)') &
         "| Total time for calculating PT2 correction              : ", &
             tot_time_mp2, " s"
          write(use_unit,'(A)')&
         "---------------------------------------------"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| OS-PT2 contribution            :", &
                E_mp2_os, " Ha", E_mp2_os * hartree ," eV"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| SS-PT2 contribution            :", &
                E_mp2_ss, " Ha", E_mp2_ss * hartree ," eV"
          write(use_unit,'(A)')&
         "---------------------------------------------"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| DHDF/DFT Energy                :", &
                total_energy," Ha", total_energy* hartree," eV"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Exx contribution               :", &
                E_dftpt2_hf, " Ha", E_dftpt2_hf*hartree," eV"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| PT2 contribution               :", &
                E_dftpt2_pt2, " Ha", E_dftpt2_pt2 * hartree ," eV"
          write(use_unit,'(A)')&
         "---------------------------------------------"

          if (flag_dftpt2_dft_part .eq. 30) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total XYG3 energy              :", &
                post_scf_total_energy," Ha", &
                post_scf_total_energy*hartree ," eV"
              write(use_unit,'(A)')&
              "---------------------------------------------"
          else if (flag_dftpt2_dft_part .eq. 31) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total xDH-PBE0 energy          :", &
                post_scf_total_energy," Ha", &
                post_scf_total_energy*hartree ," eV"
              write(use_unit,'(A)')&
              "---------------------------------------------"
          else if (flag_dftpt2_dft_part .eq. 32) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total XYGJOS energy            :", &
                post_scf_total_energy," Ha", &
                post_scf_total_energy*hartree ," eV"
              write(use_unit,'(A)')&
              "---------------------------------------------"
          else if (flag_dftpt2_dft_part .eq. 34) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total XYG5 energy              :", &
                post_scf_total_energy," Ha", &
                post_scf_total_energy*hartree ," eV"
              write(use_unit,'(A)')&
              "---------------------------------------------"
          endif
        endif !   end of if myid
      else !(lrc_pt2_started)
        E_dftpt2_pt2 = E_mp2_ss*dftpt2_Ec_sslrcPT2+E_mp2_os*dftpt2_Ec_oslrcPT2
        post_scf_total_energy = total_energy + E_dftpt2_pt2
        if(myid.eq.0) then
          write(use_unit,'(A)')&
          "---------------------------------------------"
          write(use_unit,'(A)') " "
          write(use_unit,'(A)') " "
          write(use_unit,'(2X,A)')&
          "---------------------------------------------------"
          write(use_unit,'(2X,A)') &
          "| lrc-PT2 calculation comes to finish ...  "
          write(use_unit,'(2X,A)')&
          "---------------------------------------------------"

          write(use_unit,'(2X,A,F12.3,A)') &
          "| Total time for calculating lrc-PT2 correction          : ", &
             tot_time_mp2, " s"
          write(use_unit,'(A)')&
          "---------------------------------------------"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| lrc-OS-PT2 contribution         :", &
                E_mp2_os, " Ha", E_mp2_os * hartree ," eV"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| lrc-SS-PT2 contribution         :", &
                E_mp2_ss, " Ha", E_mp2_ss * hartree ," eV"
          write(use_unit,'(A)')&
          "---------------------------------------------"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| lrc-PT2 contribution            :", &
                E_dftpt2_pt2, " Ha", E_dftpt2_pt2 * hartree ," eV"
          write(use_unit,'(A)')&
          "---------------------------------------------"
          if (flag_dftpt2_dft_part .eq. 30) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total lrc-XYG3 energy          :", &
                post_scf_total_energy," Ha", &
                post_scf_total_energy*hartree ," eV"
              write(use_unit,'(A)')&
          "---------------------------------------------"
          else if (flag_dftpt2_dft_part .eq. 34) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total lrc-XYG6 energy          :", &
                post_scf_total_energy," Ha", &
                post_scf_total_energy*hartree ," eV"
              write(use_unit,'(A)')&
              "---------------------------------------------"
          endif
        endif ! end of if myid
      endif

      return

      end subroutine evaluate_dftpt2_correlation
!****** 
!****s* FHI-aims/evaluate_dftpt2_correlation_2
!  NAME
!   evaluate_dftpt2_correlation_2
!  SYNOPSIS 

      subroutine evaluate_dftpt2_correlation_2 &
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
      real*8 :: E_mp2_a_scs, E_mp2_a_os, E_mp2_a_ss
      real*8 :: E_mp2_term, E_mp2_term_scs, E_mp2_term_os, E_mp2_term_ss

      real*8 :: E_mp2_mpi
      real*8 :: E_mp2_mpi_ss
      real*8 :: E_mp2_mpi_os

      real*8 :: eex_energy
      real*8 :: en_total


      integer :: n_homo(n_spin)
      integer :: n_lumo(n_spin)
      integer :: n_lumo_min
      integer :: n_homo_max
      integer :: min_lumo_loc, max_homo_loc, max_local_states1

      integer ::  n_shell_aux(n_species)

      real*8 :: E_mp2, E_mp2_1, E_mp2_scs, E_mp2_os, E_mp2_ss

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

      real*8 :: E_dftpt2_hf, E_dftpt2_pt2

!     Error variable
      integer mpierr
      character*150 :: info_str
      character(*), parameter :: func = &
          'evaluate_mp2_correlation_energy_2'

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

      if (.not. lrc_pt2_started) then
        if(myid.eq.0) then
          write(use_unit,'(A)')&
              "  | Note : Resolution of identity applied "
          write(use_unit,'(A)') &
              "  |        2D memory distribution applied "
          write(use_unit,'(A)') &
              "  |        Second-order perturbation theory (PT2) "
          if (n_spin.gt.1) then
            write(use_unit,'(A)') "  |        Spin-polarized system "
          else
            write(use_unit,'(A)') &
                "  |        Non spin-polarized system "
          endif
          !if (flag_frozen_core.and.i_start_mp2.eq. &
          !   0.and.n_electrons.gt.2) then
          !  write(use_unit,'(A)') "  | Note : Frozen Core approximation  "
          !else if (flag_frozen_core.and.n_electrons.gt.2) then
          !  write(use_unit,'(A)') "  | Note : Frozen Core approximation from level "
          !  write(use_unit,*) "               ",   i_start_mp2
          !endif
          ! do i_atom=1, n_atoms,1
          !  if (empty(i_atom)) then
          !    write(use_unit,'(A)') "  | Note : Counterpoise correction"
          !    exit
          !  endif
          ! enddo

          write(use_unit,'(A)') &
              "---------------------------------------------"


          write(use_unit,*) &
              " | Homo level                             ", n_homo
          write(use_unit,*) &
              " | Lumo level                             ", n_lumo
          write(use_unit,*) &
              " | Total Number of states                 ", n_states
          write(use_unit,*) &
              " | Number of k-points                     ",n_k_points

          if (n_spin.eq.2) then
           do i_spin=1,n_spin,1
            if (i_spin.eq.1) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Homo-Lumo Gap of the spin up    ", &
                KS_eigenvalue(n_lumo(i_spin),i_spin,1) &
                - KS_eigenvalue(n_homo(i_spin),i_spin,1)," Ha", &
                (KS_eigenvalue(n_lumo(i_spin),i_spin,1) &
                - KS_eigenvalue(n_homo(i_spin),i_spin,1))*hartree," eV"
            else if (n_electrons.gt.1) then
               write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
               "| Homo-Lumo Gap of the spin down  ", &
               KS_eigenvalue(n_lumo(i_spin),i_spin,1) &
               - KS_eigenvalue(n_homo(i_spin),i_spin,1)," Ha", &
               (KS_eigenvalue(n_lumo(i_spin),i_spin,1) &
               - KS_eigenvalue(n_homo(i_spin),i_spin,1))*hartree," eV"
            endif
           enddo
          else
             write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
              "| Homo-Lumo Gap                  ", &
             KS_eigenvalue(n_lumo(n_spin),n_spin,1) &
             - KS_eigenvalue(n_homo(n_spin),n_spin,1)," Ha", &
             (KS_eigenvalue(n_lumo(n_spin),n_spin,1) &
             - KS_eigenvalue(n_homo(n_spin),n_spin,1))*hartree," eV"
          endif
          write(use_unit,'(A)') &
              "---------------------------------------------"
        endif !(myid.eq.0)
      endif !(.not. lrc_pt2_started)

      ! determine the post-SCF frozen core low state by igor
      if (flag_frozen_core_postSCF) then ! count the frozen core states
          call count_frozen_core_states(i_start_mp2)
          if(myid.eq.0) then
             write(use_unit,'(2X,A,I12)') &
             "Frozen core number for post-SCF calculation          :", &
             i_start_mp2 - 1
          endif
      endif

      n_homo_max = maxval(n_homo(:))
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
      if (allocated(ovlp_3KS)) deallocate(ovlp_3KS)
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

      if(.not.lrc_pt2_started) then
        call evaluate_eex_energy_2 &
           ( n_homo_max,n_homo,occ_numbers, &
             ovlp_3KS,eex_energy &
           )
      endif
      call get_timestamps(rtime, clock_rtime)
      time_ovlp3fn_to_ovlp3KS=rtime-time_ovlp3fn_to_ovlp3KS
      clock_time_ovlp3fn_to_ovlp3KS=clock_rtime-clock_time_ovlp3fn_to_ovlp3KS

      call cpu_time(rtime)
!      tot_time_3KS_intg = rtime- time_3KS_intg


      call cpu_time(time_mp2)

      E_mp2=0.d0
      E_mp2_os=0.d0
      E_mp2_ss=0.d0
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
            E_mp2_term_ss = 0.d0
            E_mp2_term_os = 0.d0

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
                tmp_o3KS(:,:) = ovlp_3KS(:,min_lumo_loc:ndim1_o3KS,&
                                         loc_dim2_o3ks(i_state),i_spin)

              call MPI_Bcast(tmp_o3KS,n_basbas*max_local_states1,&
                             MPI_REAL8,global_id(n_p1,n_p2),&
                             mpi_comm_global,mpierr)

              ! Multiply what we got with our local part of ovlp_3KS_full(:,n_lumo_min:n_states,x:n_homo(i_spin2),i_spin2)

              do i_state2=i_state2_start,n_homo(i_spin2),1

                if(own_dim2_o3ks(i_state2) /= myp2_o3KS) cycle
                i_state2_loc = loc_dim2_o3ks(i_state2)

                call dgemm('T', 'N', max_local_states1, &
                         max_local_states1, n_basbas, 1.d0, &
                         tmp_o3KS, ubound(tmp_o3KS,1), &
                         ovlp_3KS(1,min_lumo_loc,i_state2_loc,i_spin2),&
                         ubound(ovlp_3KS,1), 0.d0, &
                         tmp_eri_full(min_lumo_loc,min_lumo_loc,n_p1,i_state2_loc),&
                         max_local_states1)
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

              E_mp2_a   = 0.0d0
              E_mp2_a_os= 0.0d0
              E_mp2_a_ss= 0.0d0

              do i_virt=n_lumo(i_spin), n_states,1
                do i_virt2=n_lumo(i_spin2),n_states,1

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
                    E_mp2_a_os=E_mp2_a*E_mp2_a
                    E_mp2_a_ss=(E_mp2_a-E_mp2_b)*E_mp2_a
                    E_mp2_a=(2.d0*E_mp2_a-E_mp2_b)*E_mp2_a
                  else
                    if(i_spin.eq.i_spin2) then
                      E_mp2_a=(E_mp2_a-E_mp2_b)*E_mp2_a
                      E_mp2_a_ss=E_mp2_a
                    else
                      E_mp2_a=E_mp2_a*E_mp2_a
                      E_mp2_a_os=E_mp2_a
                    endif
                  endif

                  ! If symmetry is exploited, take symmetric term into account
                  if (i_spin.ne.i_spin2 .or. i_state.ne.i_state2) then
                    E_mp2_a=E_mp2_a*2.d0
                    E_mp2_a_ss=E_mp2_a_ss*2.d0
                    E_mp2_a_os=E_mp2_a_os*2.d0
                  endif

                  E_mp2_1 = KS_eigenvalue(i_state,i_spin,1)   &
                          + KS_eigenvalue(i_state2,i_spin2,1) &
                          - KS_eigenvalue(i_virt,i_spin,1)    &
                          - KS_eigenvalue(i_virt2,i_spin2,1)

                  if (abs(E_mp2_1).lt.1e-6)then
                    write(use_unit,'(10X,A)')&
                        "****************************************"
                    write(use_unit,'(10X,A)') &
                        "| Warning : too close to degeneracy"
                    write(use_unit,'(10X,A)')&
                        "****************************************"
                  endif

                  E_mp2_term = E_mp2_term + E_mp2_a/E_mp2_1
                  E_mp2_term_ss = E_mp2_term_ss + E_mp2_a_ss/E_mp2_1
                  E_mp2_term_os = E_mp2_term_os + E_mp2_a_os/E_mp2_1

                enddo ! i_virt2
              enddo ! i_virt
            enddo ! i_state2

            if(use_mpi) then
              call sync_real_number(E_mp2_term)
              call sync_real_number(E_mp2_term_ss)
              call sync_real_number(E_mp2_term_os)
            endif
            E_mp2_term = E_mp2_term / dble(n_spin) 
            E_mp2_term_ss = E_mp2_term_ss / dble(n_spin) 
            E_mp2_term_os = E_mp2_term_os / dble(n_spin) 

            E_mp2 = E_mp2 + E_mp2_term
            E_mp2_ss = E_mp2_ss + E_mp2_term_ss
            E_mp2_os = E_mp2_os + E_mp2_term_os

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

      if (allocated(ovlp_3KS)) then
        deallocate (ovlp_3KS)
      endif

      call cpu_time(rtime)
      tot_time_mp2 = rtime- time_mp2

      if (.not. lrc_pt2_started) then
        E_dftpt2_hf = eex_energy*dftpt2_Ex_hf
        E_dftpt2_pt2 = E_mp2_ss*dftpt2_Ec_ssPT2+E_mp2_os*dftpt2_Ec_osPT2
        post_scf_total_energy = total_energy + E_dftpt2_hf+E_dftpt2_pt2
        if(myid.eq.0) then
          write(use_unit,'(A)')&
              "---------------------------------------------"
          write(use_unit,'(A)') " "
          write(use_unit,'(A)') " "
          write(use_unit,'(2X,A)')&
              "---------------------------------------------------"
          write(use_unit,'(2X,A)')&
              "| PT2 calculation comes to finish ...  "
          write(use_unit,'(2X,A)')&
              "---------------------------------------------------"

          write(use_unit,'(2X,A,F12.3,A)') &
          "| Total time for calculating PT2 correction              : ", &
             tot_time_mp2, " s"
          write(use_unit,'(A)')&
              "---------------------------------------------"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| OS-PT2 contribution            :", &
                E_mp2_os, " Ha", E_mp2_os * hartree ," eV"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| SS-PT2 contribution            :", &
                E_mp2_ss, " Ha", E_mp2_ss * hartree ," eV"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Full-PT2 contribution          :", &
                E_mp2, " Ha", E_mp2 * hartree ," eV"
          write(use_unit,'(A)')&
              "---------------------------------------------"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| DHDF/DFT Energy                :", &
                total_energy," Ha", total_energy* hartree," eV"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Exx contribution               :", &
                E_dftpt2_hf, " Ha", E_dftpt2_hf*hartree," eV"
          write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| PT2 contribution               :", &
                E_dftpt2_pt2, " Ha", E_dftpt2_pt2 * hartree ," eV"
          write(use_unit,'(A)')&
                "---------------------------------------------"

          if (flag_dftpt2_dft_part .eq. 30) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total XYG3 energy              :", &
                post_scf_total_energy," Ha", &
                post_scf_total_energy*hartree ," eV"
              write(use_unit,'(A)')&
                "---------------------------------------------"
          else if (flag_dftpt2_dft_part .eq. 31) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total xDH-PBE0 energy          :", &
                post_scf_total_energy," Ha", &
                post_scf_total_energy*hartree ," eV"
              write(use_unit,'(A)')&
                "---------------------------------------------"
          else if (flag_dftpt2_dft_part .eq. 32) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total XYGJOS energy            :", &
                post_scf_total_energy," Ha", &
                post_scf_total_energy*hartree ," eV"
              write(use_unit,'(A)')&
                "---------------------------------------------"
          else if (flag_dftpt2_dft_part .eq. 34) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total XYG5 energy              :", &
                post_scf_total_energy," Ha", &
                post_scf_total_energy*hartree ," eV"
              write(use_unit,'(A)')&
                "---------------------------------------------"
          endif
        endif !   end of if myid
      else !(lrc_pt2_started)
        E_dftpt2_pt2 = E_mp2_ss*dftpt2_Ec_sslrcPT2+E_mp2_os*dftpt2_Ec_oslrcPT2
        post_scf_total_energy = total_energy + E_dftpt2_pt2
        if(myid.eq.0) then
          write(use_unit,'(A)')&
          "---------------------------------------------"
          write(use_unit,'(A)') " "
          write(use_unit,'(A)') " "
          write(use_unit,'(2X,A)')&
          "---------------------------------------------------"
          write(use_unit,'(2X,A)') &
          "| lrc-PT2 calculation comes to finish ...  "
          write(use_unit,'(2X,A)')&
          "---------------------------------------------------"

          write(use_unit,'(2X,A,F12.3,A)') &
          "| Total time for calculating lrc-PT2 correction          : ", &
             tot_time_mp2, " s"
          write(use_unit,'(A)')&
          "---------------------------------------------"
          write (6,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| lrc-OS-PT2 contribution         :", &
                E_mp2_os, " Ha", E_mp2_os * hartree ," eV"
          write (6,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| lrc-SS-PT2 contribution         :", &
                E_mp2_ss, " Ha", E_mp2_ss * hartree ," eV"
          write(use_unit,'(A)')&
          "---------------------------------------------"
          write (6,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| lrc-PT2 contribution            :", &
                E_dftpt2_pt2, " Ha", E_dftpt2_pt2 * hartree ," eV"
          write(use_unit,'(A)')&
          "---------------------------------------------"
          if (flag_dftpt2_dft_part .eq. 30) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total lrc-XYG3 energy          :", &
                post_scf_total_energy," Ha", &
                post_scf_total_energy*hartree ," eV"
              write(use_unit,'(A)')&
          "---------------------------------------------"
          else if (flag_dftpt2_dft_part .eq. 34) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total lrc-XYG6 energy          :", &
                post_scf_total_energy," Ha", &
                post_scf_total_energy*hartree ," eV"
              write(use_unit,'(A)')&
              "---------------------------------------------"
          endif
        endif
      endif ! end of lrc_pt2_started

      return

      end subroutine evaluate_dftpt2_correlation_2
!****** 

