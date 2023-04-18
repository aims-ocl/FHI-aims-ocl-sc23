!****s* FHI-aims/evaluate_O_2bs1HF
!  NAME
!   evaluate_O_2bs1HF
!  SYNOPSIS

      subroutine evaluate_O_2bs1HF &
           ( KS_eigenvector)

!  PURPOSE
!  This routines performs the 3-center overlap integral transformation
!  by multiplying the ovlp_3fn with the KS eigenvectors.
!  The resultant O_2bs1HF (defined in MODULE hartree-fock) matrix elements
!  are the integrals over auxiliary (product) basis, a regular NAO basis, 
!  and a single-partilce (KS/HF) orbital.
!
!  USES

      use dimensions
      use prodbas
      use hartree_fock
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: use_unit

      implicit none

!  ARGUMENTS

      real*8, dimension(n_basis,n_states,n_spin,n_k_points) :: &
              KS_eigenvector

!  INPUT
!  o  KS_eigenvector -- real array,
!       the eigenvector of the single-particle (KS/HF) self-consistent calculation
!  OUTPUT
!   none
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

!     auxiliary matrices for Level 3 Blas matrix multiplications
!      real*8   aux_ovlp_matr(n_basis,n_basis)
!      real*8   aux_ovlp_matr_1(n_basis,n_basis)
       real*8, dimension(:), allocatable :: temp_ovlp
       real*8, dimension(:,:), allocatable :: aux_ovlp_matr
       real*8, dimension(:,:,:), allocatable :: KS_eigenvector_loc

       integer ::  mpierr

!     counters

      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_basis_3
      integer :: i_prodbas_1
      integer :: i_state
      integer :: i_state_loc
      integer :: i_index
      integer :: i_spin
      integer :: i_task

!     begin work

!     first multiplication between eigenvector and basis overlap.

      if(myid.eq.0) then
       write(use_unit, *)
       write(use_unit, '(2X,A)') &
                "Evaluating transformed overlap integrals ... "
       write(use_unit, *)
      endif

      if(use_mpi) then
         call MPI_Bcast(KS_eigenvector, &
              n_basis*n_states*n_spin*n_k_points, &
              MPI_DOUBLE_PRECISION, &
              0, mpi_comm_global, mpierr)
      endif
!  n_loc_states=(n_states-1)/n_tasks + 1
      if(.not. allocated(KS_eigenvector_loc)) then
        allocate(KS_eigenvector_loc(n_basis,n_loc_states,n_spin),stat=i_index)
        call check_allocation(i_index, 'KS_eigenvector_loc           ')
      endif

      KS_eigenvector_loc(:,:,:)=0.d0
      do i_state=1, n_states, 1
         if(myid .ne. mod(i_state-1, n_tasks)) cycle
         i_state_loc = (i_state-1)/n_tasks + 1
         KS_eigenvector_loc(:,i_state_loc,:) = KS_eigenvector(:,i_state,:,1)
      enddo

      if(.not. allocated(aux_ovlp_matr)) then
        allocate(aux_ovlp_matr(n_basis,n_basis),stat=i_index)
        call check_allocation(i_index, 'aux_ovlp_matr                 ')
      endif
      if(.not. allocated(temp_ovlp)) then
        allocate(temp_ovlp(n_basis_pairs),stat=i_index)
        call check_allocation(i_index, 'temp_ovlp                     ')
      endif

      O_2bs1HF(:,:,:,:) = 0.d0
      do i_basis_1 = 1, n_max_loc_prodbas

         i_index = 0
         aux_ovlp_matr(:,:) = 0.d0

         do i_task = 1, n_tasks, 1
            if (myid == i_task -1 ) then
               if(i_basis_1 <=n_loc_prodbas) then
                 temp_ovlp(:) = ovlp_3fn(:, i_basis_1)
               else
                 temp_ovlp(:) = 0.d0
               endif
            endif
            if(use_mpi) then
              call mpi_bcast(temp_ovlp,n_basis_pairs, &
                         MPI_REAL8, i_task-1, mpi_comm_global, mpierr)
            endif

            do i_basis_2 = 1, n_basis, 1
              do i_basis_3 = 1, i_basis_2, 1
                i_index= basis_nghr(i_basis_3, i_basis_2)

                if(i_index .gt.0) then
                   aux_ovlp_matr(i_basis_2,i_basis_3) = temp_ovlp(i_index)
                   aux_ovlp_matr(i_basis_3,i_basis_2) = temp_ovlp(i_index)
                endif

              enddo
            enddo

            i_prodbas_1 = map_prodbas(i_basis_1, i_task)
            if(i_prodbas_1 == 0) cycle
            do i_spin = 1, n_spin

              call dgemm('N', 'N', n_basis, n_loc_states, &
                  n_basis, 1.0d0, aux_ovlp_matr, n_basis, &
                  KS_eigenvector_loc(1,1,i_spin), n_basis, &
                  0.d0, O_2bs1HF(i_prodbas_1,:,:,i_spin), &
                  n_basis)
            enddo

!  end of loop over i_task
         enddo


!  end of loop over i_basis_1
      enddo

!      do i_spin = 1, n_spin
!      do i_state_loc = 1, n_loc_states, 1
!         call sync_matrix(O_2bs1HF(:,:,i_state_loc,i_spin), &
!          n_basbas,n_basis)
!       enddo
!      enddo

!    printing out
!      if(myid.eq.0) then
!       do i_basis_3 = 5, 5
!        do i_basis_1 = 1, n_basbas
!         do i_basis_2 = 1, 1
!              write(use_unit, '(2X,I6, I7,2X,I7, 2X,I7, 7X, 2F18.10)')
!     +            myid,  i_basis_1, i_basis_2, i_basis_3,
!     +              O_2bs1HF(i_basis_1,i_basis_2,i_basis_3,1),
!     +              O_2bs1HF(i_basis_1,i_basis_2,i_basis_3,2)
!            enddo
!          enddo
!       enddo
!      endif
      if(allocated(temp_ovlp)) then
        deallocate(temp_ovlp)
      endif
      if(allocated(aux_ovlp_matr)) then
        deallocate(aux_ovlp_matr)
      endif
      if(allocated(KS_eigenvector_loc)) then
        deallocate(KS_eigenvector_loc)
      endif

      end subroutine evaluate_O_2bs1HF
!---------------------------------------------------------------------
!******
