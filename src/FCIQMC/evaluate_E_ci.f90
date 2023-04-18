!****s* FHI-aims/fciqmc/evaluate_E_ci
!  NAME
!    evaluate_E_ci
!  SYNOPSIS 

      subroutine evaluate_E_ci()

!  PURPOSE
!
!  The '''evaluate_E_ci''' evaluate truncated CI energy based on given c_vect
!   and w_vect
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
!    This file was written by Jan Kloppenburg
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
    real*8  :: E_ci_1a,E_ci_1b,E_ci_2a,E_ci_2b,E_ci_1a1b
    integer :: index_start(4)
    
    if (.not.allocated(occ_num_1)) then
        allocate(occ_num_1(n_states,2), stat=errnum)
        call check_allocation(errnum, 'occ_num_1 (1D) in CI')
    endif

    if (myid .eq. 0) then
        E_ci = c_0*(E_HF*c_0 + w_0)
    else
        E_ci = 0.0d0
    endif


! dddddd

        E_ci_1a=0.0d0
        E_ci_1b=0.0d0
        E_ci_2a=0.0d0
        E_ci_2b=0.0d0
        E_ci_1a1b=0.0d0
! eeeeee




    !i_ci = 0
    if (flag_single_excitation) then
        ! ==============================================================
        ! now construct matrix elements between one-electron substutions
        ! ==============================================================
        ! First for the substitution in the alpha space
        ! --------------------------------------------------------------
        !write(use_unit,'(2X,A,F16.8)') &
        !   'Igor debug E_CI 1st excitation before :',E_ci

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
            E_ci = E_ci + c_vect(i_ci) * &
                (&
                (configuration_eigenvalue + V_00)*c_vect(i_ci) + &
                w_vect(i_ci)&
                )

            E_ci_1a = E_ci_1a + c_vect(i_ci) * &
                (&
                (configuration_eigenvalue + V_00)*c_vect(i_ci) + &
                w_vect(i_ci)&
                )

          enddo
          index_start(2) = n_elec(1)+1
        enddo ! one-electron excitations in the alpha spin space
666     continue
        !write(use_unit,'(2X,A,F16.8)') &
        !   'Igor debug E_CI 1st excitation after :',E_ci
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
            E_ci = E_ci + c_vect(i_ci) * &
                (&
                (configuration_eigenvalue + V_00)*c_vect(i_ci) + &
                w_vect(i_ci)&
                )

            E_ci_1b = E_ci_1b + c_vect(i_ci) * &
                (&
                (configuration_eigenvalue + V_00)*c_vect(i_ci) + &
                w_vect(i_ci)&
                )

          enddo
          index_start(2) = n_elec(2)+1
        enddo ! one-electron excitations in the beta spin space
667     continue

!        print*,'E_ci_1a',E_ci_1a,'E_ci_1b',E_ci_1b

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
            E_ci = E_ci + c_vect(i_ci) * &
                (&
                (configuration_eigenvalue + V_00)*c_vect(i_ci) + &
                w_vect(i_ci)&
                )

            E_ci_1a1b = E_ci_1a1b + c_vect(i_ci) * &
                (&
                (configuration_eigenvalue + V_00)*c_vect(i_ci) + &
                w_vect(i_ci)&
                )


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
            E_ci = E_ci + c_vect(i_ci) * &
                (&
                (configuration_eigenvalue + V_00)*c_vect(i_ci) + &
                w_vect(i_ci)&
                )

            E_ci_2a = E_ci_2a + c_vect(i_ci) * &
                (&
                (configuration_eigenvalue + V_00)*c_vect(i_ci) + &
                w_vect(i_ci)&
                )

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
            E_ci = E_ci + c_vect(i_ci) * &
                (&
                (configuration_eigenvalue + V_00)*c_vect(i_ci) + &
                w_vect(i_ci)&
                )

            E_ci_2b = E_ci_2b + c_vect(i_ci) * &
                (&
                (configuration_eigenvalue + V_00)*c_vect(i_ci) + &
                w_vect(i_ci)&
                )


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

!    print*,'E_ci_2a',E_ci_2a,'E_ci_2b',E_ci_2b
!    print*,'E_ci_1a1b',E_ci_1a1b,c_0*(E_HF*c_0 + w_0)

    call sync_real_number(E_ci)

    deallocate(occ_num_1)

end subroutine evaluate_E_ci

