!****s* FHI-aims/evaluate_en_mp2_correlation_energy_v02
!  NAME
!   evaluate_en_mp2_correlation_energy
!  SYNOPSIS 

      subroutine evaluate_osmp2_correlation_energy &
       (n_electrons,total_energy,en_xc, &
        occ_numbers,KS_eigenvalue,KS_eigenvector,post_scf_total_energy )

!  PURPOSE
!    Calculation of the total energy perturbation terms up to the second order using 
!    Moller-Plesset theory (MP2) to the HF/DFT Energy.
!    A variant, called SCS-MP2, scales the two spin channels separately.
!    Frozen core approximation is available too. 
!
! USES

      use constants
      use dimensions
      use prodbas
      use species_data
      use physics, only : partition_tab, rho, rho_gradient, kinetic_density
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

      real*8 :: ddot
      real*8, dimension(:,:,:,:), allocatable :: ovlp_3KS
      real*8 :: n_bytes_o3KS
      real*8 :: aux

      real*8, dimension(:,:), allocatable   :: aux_eri
      real*8, dimension(:,:,:), allocatable :: Ec_1st 
      real*8, dimension(:,:), allocatable   :: xc_matr
      real*8, dimension(:,:), allocatable :: HF_ex_matr
      real*8, dimension(:,:,:), allocatable :: KS_xc_matr

      real*8 :: E_mp2_a, E_mp2_b, E_mp2_mpi
      real*8 :: E_mp2_a_ss, E_mp2_mpi_ss
      real*8 :: E_mp2_a_os, E_mp2_mpi_os
      real*8 :: E_mp2_term, E_mp2_term_ss, E_mp2_term_os
      
      real*8 :: ec_1st_ij_os, ec_1st_ij_ss
      real*8 :: vee_matr_os, vee_matr_ss

      real*8 :: eex_energy
      real*8 :: en_total


      integer  :: n_homo(n_spin)
      integer ::  n_lumo(n_spin)
      integer ::  n_lumo_min
      integer ::  n_homo_max

      integer ::  n_shell_aux(n_species)
      integer ::  i_count

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
      integer i_state3
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
              'evaluate_osmp2_correlation_energy'
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

      write(use_unit,*) &
          " | Homo level                             ", n_homo
      write(use_unit,*) &
          " | Lumo level                             ", n_lumo
      write(use_unit,*) &
          " | Total Number of states                 ", n_states
      write(use_unit,*) &
          " | Number of k-points                     ",n_k_points

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
        write(use_unit,'(1X,A)') &
            "*--------------------------------------------"
        write(use_unit,'(1X,A)') &
            "* Error: Requested more frozen core levels"
        write(use_unit,'(1X,A)') &
            "* than occupied states. Aborting MP2 part."
        write(use_unit,'(1X,A)') &
            "*--------------------------------------------"
        stop
     endif

!      call cpu_time(time_3KS_intg)

      call get_timestamps(time_ovlp3fn_to_ovlp3KS, clock_time_ovlp3fn_to_ovlp3KS)
      if(.not. allocated(ovlp_3KS)) then
        allocate(ovlp_3KS(n_loc_prodbas, n_states, &
        &                 1:n_homo_max, n_spin), stat=info)
        n_bytes_o3KS = 8 * n_loc_prodbas * n_states * n_homo_max * n_spin
        call check_allocation(info, 'ovlp_3KS', func)

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
        end if
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
        write(use_unit,'(A)') &
            "---------------------------------------------"
        write(use_unit,'(A)') &
            "| Sum of correlation terms "
        write(use_unit,'(A)') &
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
      !if (.not.allocated(denominator)) then
      !    allocate(denominator(n_lumo_min:n_states,n_lumo_min:n_states))
      !endif

! -------------------------------------------------------------------
!    NON-Spin polarized systems
! -------------------------------------------------------------------

      if (n_spin.eq.1) then
        !i_index = 0
        do i_state=i_start_mp2, n_homo(n_spin),1
          do i_state2=i_state,n_homo(n_spin),1

            !denominator     = 0.d0
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
            i_count = 0
            do i_virt=n_lumo(n_spin), n_states,1
              ! ================================================
              ! Now calculate the double-excitation contribution
              ! ================================================
              do i_virt2=n_lumo(n_spin),n_states,1
                 i_count = i_count + 1
                ! MPI task distribution
                if(myid.eq.MOD(i_count,n_tasks)) then
                  E_mp2_a=aux_eri(i_virt,i_virt2)
                  !E_mp2_b=aux_eri(i_virt2,i_virt)

                  ! separate parallel-spin and opposite-spin contributions
                  !E_mp2_a_ss = (E_mp2_a-E_mp2_b)* E_mp2_a
                  E_mp2_a_os = E_mp2_a * E_mp2_a 

                  if (i_state.ne.i_state2) then
                     !E_mp2_a_ss=E_mp2_a_ss*2.d0
                     E_mp2_a_os=E_mp2_a_os*2.d0
                  endif

                  E_mp2_1= (KS_eigenvalue(i_virt,n_spin, 1)) &
                    +  (KS_eigenvalue(i_virt2,n_spin,1) )

                  E_mp2_1= E_mp2_1-(KS_eigenvalue(i_state,n_spin,1))

                  E_mp2_1= E_mp2_1-(KS_eigenvalue(i_state2,n_spin,1))

                  !E_mp2_1= E_mp2_1-Ec_1st(i_state,i_state2,n_spin)

                  if (abs(E_mp2_1).lt.1e-6)then
                     write(use_unit,'(10X,A)') &
                          "****************************************"
                     write(use_unit,'(10X,2A)') "| Warning :", &
                       " too close to degeneracy"
                     write(use_unit,'(10X,A)') &
                          "****************************************"
                  endif

                  !denominator(i_virt,i_virt2) = E_mp2_1
                  !denominator(i_virt,i_virt2) = E_mp2_1 &
                  !    - Ec_1st(i_state,i_state,1) &
                  !    - Ec_1st(i_state2,i_state2,1)

                  E_mp2_1= 1.d0 / E_mp2_1

                  !E_mp2_term    = (E_mp2_a_ss+E_mp2_a_os)*E_mp2_1 + E_mp2_term
                  !E_mp2_term_ss = E_mp2_a_ss*E_mp2_1 + E_mp2_term_ss
                  E_mp2_term_os = E_mp2_a_os*E_mp2_1 + E_mp2_term_os
!   end of MPI distribution
                endif
!   close b  virtual state 2
              enddo
!   close a virtual state
            enddo
            if(use_mpi) then
                ! sync second-order terms
                !call MPI_ALLREDUCE(E_mp2_term, E_mp2_mpi, 1, &
                !       MPI_DOUBLE_PRECISION, &
                !       MPI_SUM, mpi_comm_global, mpierr)
                !E_mp2_term = -E_mp2_mpi
                !call MPI_ALLREDUCE(E_mp2_term_ss, E_mp2_mpi_ss, 1, &
                !       MPI_DOUBLE_PRECISION, &
                !       MPI_SUM, mpi_comm_global, mpierr)
                !E_mp2_term_ss = -E_mp2_mpi_ss
                call MPI_ALLREDUCE(E_mp2_term_os, E_mp2_mpi_os, 1, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_SUM, mpi_comm_global, mpierr)
                E_mp2_term_os = -E_mp2_mpi_os
                !call sync_matrix( denominator(n_lumo_min:n_states, &
                !                  n_lumo_min:n_states),n_unoccupied, &
                !                  n_unoccupied &
                !                 )
            endif

            !if (myid.eq.0) then
            !    write(use_unit,*) "igor debug", denominator
            !endif

            !E_mp2    = E_mp2    + E_mp2_term
            !E_mp2_ss = E_mp2_ss + E_mp2_term_ss
            E_mp2_os = E_mp2_os + E_mp2_term_os

            if (myid .eq. 0) then
               write(use_unit,'(A,I5,I5)') &
                   'The electron pair: ',i_state, i_state2
            endif

          enddo !   close i state2
        enddo !   close i state

! -------------------------------------------------------------------
!     Spin polarized systems
! -------------------------------------------------------------------
      else ! (n_spin.ne.1)
        !i_index = 0
        do i_spin=1,n_spin,1
          do i_spin2=1,n_spin,1
            if (i_spin .eq. i_spin2) cycle
            do i_state=i_start_mp2, n_homo(i_spin),1
              do i_state2=i_start_mp2,n_homo(i_spin2) ,1

                  !denominator     = 0.d0
                  !E_mp2_term      = 0.d0
                  !E_mp2_term_ss   = 0.d0
                  E_mp2_term_os   = 0.d0
                  !E_iepa_term     = 0.d0
                  !E_iepa_term_ss  = 0.d0
                  !E_iepa_term_os  = 0.d0
                  !i_index         = i_index + 1

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

                  i_count = 0
                  do i_virt=n_lumo(i_spin), n_states,1
                    do i_virt2=n_lumo(i_spin2),n_states,1
                      i_count = i_count + 1
                      ! MPI task distribution
                      if(myid.eq.MOD(i_count,n_tasks)) then

                        E_mp2_a=aux_eri(i_virt,i_virt2)

                        E_mp2_1=  KS_eigenvalue(i_virt,i_spin, 1) &
                           +  KS_eigenvalue(i_virt2,i_spin2,1)
                        E_mp2_1=  E_mp2_1 - &
                                KS_eigenvalue(i_state,i_spin,1)
                        E_mp2_1=  E_mp2_1 - &
                                KS_eigenvalue(i_state2,i_spin2,1)

                        if (abs(E_mp2_1).lt.1e-6)then
                           write(use_unit,'(10X,A)') &
                            "****************************************"
                           write(use_unit,'(10X,2A)') "| Warning :", &
                               "too close to degeneracy"
                           write(use_unit,'(10X,A)') &
                            "****************************************"
                        endif

                        E_mp2_a       = E_mp2_a * E_mp2_a
                        E_mp2_term_os = E_mp2_a/E_mp2_1+ E_mp2_term_os

                      endif !   end of MPI distribution
                    enddo !   close b  virtual state 2
                  enddo !   close a virtual state
                  if(use_mpi) then
                      call MPI_ALLREDUCE(E_mp2_term_os, E_mp2_mpi_os, 1, &
                             MPI_DOUBLE_PRECISION, &
                             MPI_SUM, mpi_comm_global, mpierr)
                      E_mp2_term_os = E_mp2_mpi_os
                  endif
                  E_mp2_os = E_mp2_os - 0.5d0*E_mp2_term_os

              enddo ! close j state
            enddo ! close i state
          enddo ! close beta spin
        enddo ! close alpha spin
      endif


      if (allocated(ovlp_3KS)) then
        deallocate (ovlp_3KS)
      endif
      if (allocated(aux_eri)) then
        deallocate (aux_eri)
      endif

      call cpu_time(rtime)
      tot_time_mp2 = rtime- time_mp2

      if (use_dftpt2) then
        en_total = total_energy + eex_energy * dftpt2_Ex_hf + &
                   E_mp2_os * dftpt2_Ec_osPT2
        if(myid.eq.0) then
            write(use_unit,'(A)') &
                "---------------------------------------------"
            write(use_unit,'(A)') " "
            write(use_unit,'(A)') " "
            write(use_unit,'(A)') " "

            write(use_unit,'(10X,A)') &
                "---------------------------------------------------"

            write(use_unit,'(10X,A)') &
                "| DH-DFT calculation comes to finish ...  "
            write(use_unit,'(A)') " "
            write(use_unit,'(10X,A)') &
                "---------------------------------------------------"

             write(use_unit,'(10X,A,F12.3,A)') &
      "| Total time for calculating osPT2 correction              : ", &
                 tot_time_mp2, " s"

            write(use_unit,'(A)')
            write(use_unit,'(10X,A)') &
                "---------------------------------------------"

            write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| DHDFT/DFT Energy                :", &
                total_energy, " Ha", total_energy* hartree," eV"
            write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| DFT XC Energy                   :", &
                en_xc, " Ha", en_xc* hartree," eV"
            write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Exact Echange Energy            :", &
                eex_energy, " Ha", eex_energy* hartree," eV"
            write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| PT2 correction (os)             :", &
                E_mp2_os, " Ha"  , (E_mp2_os) * hartree ," eV"
            if (flag_dftpt2_dft_part .eq. 31) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total xDH-PBE0 Energy           :", &
                en_total, " Ha"  , en_total * hartree ," eV"
            else if (flag_dftpt2_dft_part .eq. 32) then
              write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total XYGJOS Energy             :", &
                en_total, " Ha"  , en_total * hartree ," eV"
            endif

            write(use_unit,'(A)')
            post_scf_total_energy = en_total
!   end of if myid
        endif
      else
        !if(.not.use_hartree_fock) then
        en_total = total_energy - en_xc + eex_energy + E_mp2_os
        !endif
        if(myid.eq.0) then
            write(use_unit,'(A)') &
                "---------------------------------------------"
            write(use_unit,'(A)') " "
            write(use_unit,'(A)') " "
            write(use_unit,'(A)') " "

            write(use_unit,'(10X,A)') &
                "---------------------------------------------------"

            write(use_unit,'(10X,A)') &
                "| osMP2 calculation comes to finish ...  "
            write(use_unit,'(A)') " "
            write(use_unit,'(10X,A)') &
                "---------------------------------------------------"

             write(use_unit,'(10X,A,F12.3,A)') &
      "| Total time for calculating SCPT2 correction              : ", &
                tot_time_mp2, " s"

            write(use_unit,'(A)')
            write(use_unit,'(10X,A)') &
                "---------------------------------------------"

            write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| DFT or HF Energy                :", &
                total_energy, " Ha", total_energy* hartree," eV"
            write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| DFT XC Energy                   :", &
                en_xc, " Ha", en_xc* hartree," eV"
            write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Exact Echange Energy            :", &
                eex_energy, " Ha", eex_energy* hartree," eV"
            write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| MP2 correction (os)             :", &
                E_mp2_os, " Ha"  , (E_mp2_os) * hartree ," eV"
            write(use_unit,'(2X,A,1X,F19.8,A,1X,F19.8,A)') &
                "| Total (DFT+MP2) Energy          :", &
                en_total, " Ha"  , en_total * hartree ," eV"

            write(use_unit,'(A)')
            post_scf_total_energy = en_total
!   end of if myid
        endif
      endif

      return

      end subroutine evaluate_osmp2_correlation_energy
!****** 
