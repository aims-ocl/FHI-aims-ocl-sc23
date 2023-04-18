!****s* FHI-aims/evalute_bare_ci
!  NAME
!   evaluate_bare_ci
!  SYNOPSIS

      subroutine evaluate_bare_ci &
        ()

!  PURPOSE
!  Subroutine evaluate_bare_ci evaluates and writes all bare coulomb
!  integrals for molecules
!
!  USES
        use dimensions
        use runtime_choices
        use species_data
        use physics
        use prodbas
        use hartree_fock
        use constants
        use mpi_tasks
        use synchronize_mpi
        use localorb_io, only: use_unit
        implicit none

!  ARGUMENTS
!    none
!  OUTPUT
!    none

!  local variables :
        real*8 ddot
!
!   indices
        integer :: i_state
        integer :: i_spin
!    ovlp_pb2KS : (ij|mu)_sigma where i,j are the indices of molecular orbitals and
!       mu denotes the product basis index
        real*8, dimension(:,:,:,:), allocatable :: ovlp_pb2KS
!    bare_ci : stores the bare coulomb integrals
        real*8  bare_ci(n_states,n_spin)
!    mpi error variable
        integer mpierr


        if(use_mpi) then
            call MPI_Bcast(KS_eigenvector, &
            n_basis*n_states*n_spin*n_k_points, &
            MPI_DOUBLE_PRECISION, &
            0, mpi_comm_global, mpierr)
        endif

        if (myid .eq. 0) then
          write(use_unit,'(2X,A)') &
          "------------------------------------------------------------"
          write(use_unit,'(2X,A)') &
          " Evaluation and output of molecular bare Coulomb integrals:"
        end if

        if (.not.allocated(ovlp_pb2KS)) then
            allocate(ovlp_pb2KS(n_loc_prodbas,n_states, &
            n_states,n_spin) )
        end if

        call transform_ovlp3fn (n_states,KS_eigenvector, ovlp_pb2KS)

        do i_spin = 1, n_spin
          do i_state = 1, n_states
            bare_ci (i_state,i_spin) = &
                ddot( n_loc_prodbas, &
                    ovlp_pb2KS(:, i_state, i_state, i_spin),1, &
                    ovlp_pb2KS(:, i_state, i_state, i_spin),1)
          enddo
        enddo

        call sync_matrix(bare_ci,n_states,n_spin)

        ! write results
        if (myid .eq. 0) then
            do i_spin = 1, n_spin

                write(use_unit,*)
                if (spin_treatment.eq.1) then
                    write(use_unit,'(2X,A,I2)') "Spin ", i_spin 
                end if

                write(use_unit,'(2X,A,4X,A,4X,A)') &
                "State", "molecular CI [Ha]", "molecular CI [eV]"

                do i_state = 1, n_states
                  write(use_unit,'(2X,I5,6X,F8.5,8X,F14.6)') &
                    i_state, bare_ci(i_state, i_spin), &
                    bare_ci(i_state, i_spin)*hartree
                enddo

            enddo
        end if


!       cleaning up
        if (allocated(ovlp_pb2KS)) then
          deallocate(ovlp_pb2KS)
        end if
        return
        end
