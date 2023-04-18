!****s* FHI-aims/fciqmc/ci_acc_Es.f90
!  NAME
!    initialize_ci_acc_Es
!  SYNOPSIS 

      subroutine initialize_ci_acc_Es()

!  PURPOSE
!
!  The '''initialize_ci_acc_Es''' calculate and save the eigenvalues of Fock opperator
!
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!    This file was written by Tonghao SHen
!  SOURCE

  use dimensions
  use runtime_choices
  use physics
  use mpi_tasks
  use mpi_utilities
  use fciqmc_module

    implicit none
    ! temp variables
    integer, dimension(:,:), allocatable :: occ_num_1
    ! temp indices
    integer :: i_state, j_state, a_state, b_state, i_spin
    integer :: i_ci
    integer :: errnum
    !character*128 :: info_str
    !logical :: exit_flag
    real*8  :: configuration_eigenvalue
    integer :: index_start(4)
    
    if (.not.allocated(occ_num_1)) then
        allocate(occ_num_1(n_states,2), stat=errnum)
        call check_allocation(errnum, 'occ_num_1 (1D) in CI')
    endif

    !i_ci = 0
    if (flag_single_excitation) then
        ! ==============================================================
        ! now construct matrix elements between one-electron substutions
        ! ==============================================================
        ! First for the substitution in the alpha space
        ! --------------------------------------------------------------

        i_ci = memory_table(1,2,1,myid+1)-1
        index_start = index_table(:,1,myid+1)
        do i_state = index_start(1), n_elec(1), 1
          do a_state = index_start(2), n_states, 1

            i_ci = i_ci + 1
            if (i_ci .gt. memory_table(2,2,1,myid+1)) goto 666

            occ_num_1 = occ_num_0
            occ_num_1(i_state,1) = 0
            occ_num_1(a_state,1) = 1
            call evaluate_eigenvalue_configuration(occ_num_1,configuration_eigenvalue)
            ci_acc_Es(i_ci) = configuration_eigenvalue

          enddo
          index_start(2) = n_elec(1)+1
        enddo ! one-electron excitations in the alpha spin space
666     continue
        ! --------------------------------------------------------------
        ! Second for the substitution in the beta space
        ! --------------------------------------------------------------
        i_ci = memory_table(1,2,2,myid+1)-1
        index_start = index_table(:,2,myid+1)
        do j_state = index_start(1), n_elec(2), 1
          do b_state = index_start(2), n_states, 1

            i_ci = i_ci + 1
            if (i_ci .gt. memory_table(2,2,2,myid+1)) goto 667

            occ_num_1 = occ_num_0
            occ_num_1(j_state,2) = 0
            occ_num_1(b_state,2) = 1
            call evaluate_eigenvalue_configuration(occ_num_1,configuration_eigenvalue)
            ci_acc_Es(i_ci) = configuration_eigenvalue

          enddo
          index_start(2) = n_elec(2)+1
        enddo ! one-electron excitations in the beta spin space
667     continue


    endif ! if (flag_single_excitation) then
    ! ==============================================================
    ! now construct matrix elements between two-electron substutions
    ! ==============================================================
    ! --------------------------------------------------------------
    ! First for the substitution in both spaces
    ! --------------------------------------------------------------
    i_ci = memory_table(1,2,3,myid+1)-1
    index_start = index_table(:,3,myid+1)
    do i_state = index_start(1), n_elec(1), 1
      do a_state = index_start(2), n_states, 1
        do j_state = index_start(3), n_elec(2), 1
          do b_state = index_start(4), n_states, 1

            i_ci = i_ci + 1
            if (i_ci .gt. memory_table(2,2,3,myid+1)) goto 668

            occ_num_1 = occ_num_0
            occ_num_1(i_state,1) = 0  !occ_alpha
            occ_num_1(j_state,2) = 0  !occ_beta
            occ_num_1(a_state,1) = 1  !vir_alpha
            occ_num_1(b_state,2) = 1  !vir_beta
            call evaluate_eigenvalue_configuration(occ_num_1,configuration_eigenvalue)
            ci_acc_Es(i_ci) = configuration_eigenvalue

          enddo
          index_start(4) = n_elec(2)+1
        enddo
        index_start(3) = i_start_ci
      enddo
      index_start(2) = n_elec(1)+1
    enddo ! two-electron excitations in both spin spaces
668 continue
    ! --------------------------------------------------------------
    ! Second for the substitution in the alpha spin space
    ! --------------------------------------------------------------
    i_ci = memory_table(1,2,4,myid+1)-1
    index_start = index_table(:,4,myid+1)
    do i_state = index_start(1), n_elec(1), 1
      do j_state = index_start(3), n_elec(1), 1 
        do a_state = index_start(2), n_states, 1
          do b_state = index_start(4), n_states, 1

            i_ci = i_ci + 1
            if (i_ci .gt. memory_table(2,2,4,myid+1)) goto 669

            occ_num_1 = occ_num_0
            occ_num_1(i_state,1) = 0  !occ_alpha
            occ_num_1(j_state,1) = 0  !occ_alpha
            occ_num_1(a_state,1) = 1  !vir_alpha
            occ_num_1(b_state,1) = 1  !vir_alpha
            call evaluate_eigenvalue_configuration(occ_num_1,configuration_eigenvalue)
            ci_acc_Es(i_ci) = configuration_eigenvalue

          enddo
          if(a_state .eq. n_states) then
              index_start(4) = n_elec(1)+2
          else
              index_start(4) = a_state+2
          endif
        enddo
        index_start(2) = n_elec(1)+1 
      enddo
      index_start(3) = i_state+2
    enddo ! two-electron excitations in the alpha spin space
669 continue
    ! --------------------------------------------------------------
    ! Second for the substitution in the beta spin space
    ! --------------------------------------------------------------
    i_ci = memory_table(1,2,5,myid+1)-1
    index_start = index_table(:,5,myid+1)
    do i_state = index_start(1), n_elec(2), 1
      do j_state = index_start(3), n_elec(2), 1 
        do a_state = index_start(2), n_states, 1
          do b_state = index_start(4), n_states, 1

            i_ci = i_ci + 1
            if (i_ci .gt. memory_table(2,2,5,myid+1)) goto 670

            occ_num_1 = occ_num_0
            occ_num_1(i_state,2) = 0  !occ_beta
            occ_num_1(j_state,2) = 0  !occ_beta
            occ_num_1(a_state,2) = 1  !vir_beta
            occ_num_1(b_state,2) = 1  !vir_beta
            call evaluate_eigenvalue_configuration(occ_num_1,configuration_eigenvalue)
            ci_acc_Es(i_ci) = configuration_eigenvalue

          enddo
          if(a_state .eq. n_states) then
              index_start(4) = n_elec(2)+2
          else
              index_start(4) = a_state+2
          endif
        enddo
        index_start(2) = n_elec(2)+1
      enddo
      index_start(3) = i_state+2
    enddo ! two-electron excitations in the beta spin space
670 continue


    deallocate(occ_num_1)

end subroutine initialize_ci_acc_Es

