!****s* FHI-aims/evaluate_screx_matrix
!  NAME
!   evaluate_screx_matrix
!  SYNOPSIS

      subroutine evaluate_screx_matrix &
           ( KS_eigenvalue, &
             occ_numbers, KS_eigenvector, &
             transformed_ovlp3fn)

!  PURPOSE
!  Subroutine evaluate_screx_matrix evaluates the screened exchange matrix
!  to be added to the noninterecting (o.k, the Hartree contribution is there)
!  Hamiltonian H0
!
!  The screx_matr contains the screened exchange matrix if use_screx = .true.,
!  and the COHSEX matrix if use_cohsex = .true.

!  USES

      use dimensions
      use prodbas
      use hartree_fock
      use gw_para
      use mpi_tasks
      use synchronize_mpi
      use runtime_choices
      use constants
      use evaluate_polarisability_freq, only: evaluate_polarisability_freq_3
      use localorb_io, only: use_unit

      implicit none

!  ARGUMENTS

!      integer :: n_homo_max
      real*8, dimension(n_states,n_spin,n_k_points) :: KS_eigenvalue
      real*8, dimension(n_states,n_spin,n_k_points) :: occ_numbers
      real*8, dimension(n_basis,n_states,n_spin,n_k_points) :: &
                       KS_eigenvector
      real*8, dimension(n_basbas,n_basis,n_loc_states,n_spin) :: &
                       transformed_ovlp3fn

!  INPUTS
!  o  number_of_loops -- integer number, the current self-consistent loop
!  o  KS_eigenvalue -- real array,
!        the eigenvalues of the single-particle calculation. For DFT calculation,
!        this is the KS eigenvalue, but for HF calculation, this is then the HF
!        eigenvalue
!  o  KS_eigenvector -- real array,
!        the eigenvector of the single-particle calculation
!  o  transformed_ovlp3fn  --   real array, the 3-point integrals over auxiliary 
!        (product) basis, a regular NAO basis, and a single-partilce (KS/HF) orbital. 
!
!  OUTPUT
!  o none
!  "screx_matr" is defined in MODULE hatree_fock 
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

       integer ::  n_first(n_spin)
!     auxiliary matrices for Level 3 Blas matrix multiplications

       real*8, dimension(:,:), allocatable ::  screened_coulomb
       real*8, dimension(:,:), allocatable ::  aux_ovlp_matr
       real*8, dimension(:,:), allocatable ::  aux_screx_matr

      real*8, parameter ::  omega = 1.d-5
!     counters

      integer :: i_state
      integer :: i_state_loc
      integer i_basis_1
      integer i_basis_2
      integer i_basis_3
      integer i_prodbas_1
      integer i_spin
      integer i_state_count
      integer i_task
      integer i_index

!     begin work

!      do i_state = 1, n_states
!       do i_spin = 1, n_spin
!         if(KS_eigenvalue(i_state, i_spin,1).le.chemical_potential) then
!           n_homo(i_spin)=i_state
!         endif
!       enddo
!      enddo


      if(myid.eq.0) then
        write(use_unit,'(2X,2A)') &
              "Evaluating polarisability at zero energy ", &
                         "and screened exchange matrix ... "
      endif

!     determine the highest occupied orbital level
!     such complication occurs when some of the orbitals are
!     either not fully occupied or not fully empty
      n_first(:) = 1
      do i_spin = 1, n_spin
       do i_state = 1, n_states
        if (abs(occ_numbers(i_state,i_spin,1)-dble(2/n_spin)) &
                         .lt.1.d-6) then
         n_first(i_spin)= i_state + 1
        endif
       enddo
       if(n_first(i_spin) .gt. n_states) then
         n_first(i_spin) = n_states
       endif
      enddo


      if(.not.allocated(screened_coulomb)) then
        allocate(screened_coulomb(n_basbas,n_basbas),stat=i_index)
        call check_allocation(i_index, 'screened_coulomb_f            ')
 
        screened_coulomb(:,:) = 0.d0
      endif

!    screened_coulomb now contains the polarizability
      call evaluate_polarisability_freq_3 &
           ( n_low_state, n_homo, n_first, n_loc_states, &
             occ_numbers, 0.d0, &
             KS_eigenvalue, KS_eigenvector,  transformed_ovlp3fn, &
             screened_coulomb &
           )

!    this is really the screened coulomb interaction W at zero frequency
      call screened_coulomb_interaction_2 &
            ( screened_coulomb )

      do i_basis_1 = 1, n_basbas, 1
         screened_coulomb(i_basis_1, i_basis_1) = &
            screened_coulomb(i_basis_1, i_basis_1) -1.d0
      enddo

      allocate(aux_ovlp_matr(n_basbas,n_basis),stat=i_index)
      call check_allocation(i_index, 'aux_ovlp_matr                 ')

      allocate(aux_screx_matr(n_basis,n_basis),stat=i_index)
      call check_allocation(i_index, ' aux_screx_matr               ')


      screx_matr=0.d0

      do i_spin = 1, n_spin, 1
        do i_state = 1, n_states, 1

          if(i_state .gt. n_homo(i_spin) .and. &
              .not.  use_cohsex)  exit

          if(myid.ne.mod(i_state-1,n_tasks)) cycle
          i_state_loc = (i_state-1)/n_tasks + 1

          call dgemm('N', 'N', n_basbas, n_basis, &
                        n_basbas, 1.0d0, &
                       screened_coulomb, n_basbas, &
                        transformed_ovlp3fn(:,:,i_state_loc,i_spin), &
                        n_basbas, 0.d0, &
                        aux_ovlp_matr, n_basbas &
                     )

          call dgemm('T', 'N', n_basis, n_basis, &
                        n_basbas, 1.0d0, &
                        transformed_ovlp3fn(:,:,i_state_loc,i_spin), &
                        n_basbas, aux_ovlp_matr, &
                        n_basbas, 0.d0, &
                        aux_screx_matr, n_basis &
                      )

          if (use_cohsex) then
!  COHSEX approximation
             if(i_state.le.n_homo(i_spin)) then
                screx_matr(:,:,i_spin) = screx_matr(:,:,i_spin) + &
                 aux_screx_matr(:,:) * &
                 occ_numbers(i_state,i_spin,1)*dble(n_spin/4.d0)
             endif
             if(i_state.ge.n_homo(i_spin)) then
               screx_matr(:,:,i_spin) = screx_matr(:,:,i_spin) - &
                aux_screx_matr(:,:) * &
                (dble(2.d0/n_spin)-occ_numbers(i_state,i_spin,1)) &
                *dble(n_spin/4.d0)
             endif
!             write(use_unit,'(2I4,2f18.6)') i_state, i_state_loc, aux_screx_matr(1,1), aux_screx_matr(1,1)
          else
! SEX approximation
           screx_matr(:,:,i_spin) = screx_matr(:,:,i_spin) + &
              aux_screx_matr(:,:) &
              *occ_numbers(i_state,i_spin,1)*dble(n_spin/2.d0)
          endif

!          write(use_unit,'(2f18.10)')screx_matr(1,1,i_spin),aux_screx_matr(1,1)

! end of i_state loop
       enddo

      enddo
!  synchronise the exchange matrix
      do i_spin = 1, n_spin

        call sync_matrix(screx_matr(:,:,i_spin), &
                 n_basis, n_basis)
      enddo

      if(allocated(aux_ovlp_matr)) then
         deallocate(aux_ovlp_matr)
      endif
      if(allocated(aux_screx_matr)) then
         deallocate(aux_screx_matr)
      endif
      if(allocated(screened_coulomb)) then
         deallocate(screened_coulomb)
      endif
      end subroutine evaluate_screx_matrix
!---------------------------------------------------------------------
!******
