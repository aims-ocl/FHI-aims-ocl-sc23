!****s* FHI-aims/fciqmc/ci_calculation
!  NAME
!   ci_calculation
!  SYNOPSIS 

 subroutine ci_calculation ( )

!  PURPOSE
!    Calculation of the correlation energy based on the configuration interaction algorithm
!    Frozen core approximation is available too. 
! USES

      use physics
      use mpi_tasks
      use synchronize_mpi
      use runtime_choices
      use timing
      use fciqmc_module

      implicit none

!  ARGUMENTS

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
    ! variables
    real*8 :: dE_ci, dE_ci_prev, E_ci_prev
    real*8 :: E_corr_MP2, E_corr_MP3
    integer :: i_scf
    ! temp indices
    integer :: i_state! j_state, i_spin
    integer :: i_start
    integer :: errnum
    integer :: n_valence_a, n_valence_b
    integer :: info
    !character*128 :: info_str
    ! temp indices for timing
    real*8  :: time_CI,     clock_time_CI
    real*8  :: time_w_1a,   clock_time_w_1a
    real*8  :: time_w_1a1b, clock_time_w_1a1b
    real*8  :: time_w_2a,   clock_time_w_2a
    real*8  :: time_E_ci,   clock_time_E_ci
    real*8  :: time_c_vect, clock_time_c_vect
    real*8  :: tot_time_CI,     tot_clock_time_CI
    real*8  :: tot_time_w_1a,   tot_clock_time_w_1a
    real*8  :: tot_time_w_1a1b, tot_clock_time_w_1a1b
    real*8  :: tot_time_w_2a,   tot_clock_time_w_2a
    real*8  :: tot_time_E_ci,   tot_clock_time_E_ci
    real*8  :: tot_time_c_vect, tot_clock_time_c_vect

    ! debug
    integer :: unit_w_vect, unit_c_vect

    Logical :: flag_ci_acc_S_singular

    
    ! for timing
    call get_timestamps(time_CI, clock_time_CI)

    if (myid.eq.0) then
        write(use_unit,'(2X,A)') "--------------------------------------------------------"
        write(use_unit,'(2X,A)') "|"
        write(use_unit,'(2X,A)') "| Configuration interaction (CI) calculation starts ... "
        write(use_unit,'(2X,A)') "|"
        write(use_unit,'(2X,A)') "--------------------------------------------------------"
    endif
    
    select case(ci_truncation)
    case(-1)
      if (myid.eq.0) then
          write(use_unit,'(2X,A)') &
              "FCIQMC : Full CI calculation by QMC"
      endif
      ! prepare two- and four-index integrals
      !call ci_initialization()
      ! generate FCIDUMP for NECI
      call output_FCIDUMP_cluster_1d()
      goto 888

    case(2)
      ! prepare two- and four-index integrals
      call ci_initialization()


      ! determin the number of configurations
      n_configuration      = 0 ! total number of configurations
      n_configuration_1a   = 0 ! one-electron excitations in the alpha space
      n_configuration_1b   = 0 ! one-electron excitations in the beta space
      n_configuration_2a   = 0 ! two-electron excitations in the alpha space
      n_configuration_2b   = 0 ! two-electron excitations in the beta space
      n_configuration_1a1b = 0 ! one in alpha, one in beta
      if (myid.eq.0) then
          if (flag_single_excitation) then
              write(use_unit,'(2X,A,A)') &
                  "CISD : CI algorithm with ", &
                  "single and double excitations"
          else
              write(use_unit,'(2X,A,A)') &
                  "CID : CI algorithm with ", &
                  "double excitations"
          endif
      endif
      ! the number of valence electrons in the alpha spin space
      n_valence_a = n_elec(1) - i_start_ci + 1
      ! the number of valence electrons in the beta spin space
      n_valence_b = n_elec(2) - i_start_ci + 1
      !=================================================================
      ! the one-electron excitations has the same form for all the cases
      !=================================================================
      if (flag_single_excitation) then
          ! one substitusion in alpha space
          n_configuration_1a = n_configuration_1a + &
              n_valence_a*(n_states-n_elec(1))
          ! one substitusion in beta space
          n_configuration_1b = n_configuration_1b + &
              n_valence_b*(n_states-n_elec(2))
      endif
      !=================================================================
      ! the two-electron excitations have to de considered respectively
      !=================================================================
      ! one substitusion in alpha space, and one in beta
      n_configuration_1a1b = n_configuration_1a1b + &
          n_valence_a*(n_states-n_elec(1)) * &   ! possibility in the alpha spin space
          n_valence_b*(n_states-n_elec(2))       ! possibility in the beta spin space
      ! two substitusions in alpha space
      n_configuration_2a = n_configuration_2a + &
          (n_valence_a*(n_valence_a-1))/2 * &  ! possibility of excitation
          ((n_states-n_elec(1))*(n_states-n_elec(1)-1))/2 ! possibility to occupy
      ! two substitusions in beta space
      n_configuration_2b = n_configuration_2b + &
          (n_valence_b*(n_valence_b-1))/2 * &  ! possibility of excitation
          ((n_states-n_elec(2))*(n_states-n_elec(2)-1))/2 ! possibility of occupation
      !==================================================================
      ! NOTICE :: The convention here is that the alpha-electron number
      !           is alway larger than the beta-electron number.
      !           Therefore we don't need to consider the system where
      !           the beta spin number is larger than the alpha spin number
      !==================================================================

      n_configuration = n_configuration_1a + n_configuration_1b + &
          n_configuration_1a1b + n_configuration_2a + n_configuration_2b

      n_configuration_vector(1) = 0
      n_configuration_vector(2) = n_configuration_1a
      n_configuration_vector(3) = n_configuration_1a + n_configuration_1b
      n_configuration_vector(4) = n_configuration_1a + n_configuration_1b + n_configuration_1a1b
      n_configuration_vector(5) = n_configuration - n_configuration_2b

      if (myid .eq. 0) then
          write(use_unit,'(2X,A,I8,A)') &
              "| There are ", n_configuration, " configurations in this CI calculation, including"
          write(use_unit,'(2X,A,I8,A)') &
              "|           ", n_configuration_1a, " one-electron excitations in the alpha spin space"
          write(use_unit,'(2X,A,I8,A)') &
              "|           ", n_configuration_1b, " one-electron excitations in the beta spin space"
          write(use_unit,'(2X,A,I8,A)') &
              "|           ", n_configuration_1a1b, &
              " two-electron excitations in both spin spaces respectively"
          write(use_unit,'(2X,A,I8,A)') &
              "|           ", n_configuration_2a, " two-electron excitations in the alpha spin space"
          write(use_unit,'(2X,A,I8,A)') &
              "|           ", n_configuration_2b, " two-electron excitations in the beta spin space"
      endif

      !if (n_spin .eq. 1 .and. n_elec(1) .eq. n_elec(2)) then
      !    n_configuration = n_configuration - n_configuration_2b
      !endif

      call parallel_initialization()

      !if (.not.allocated(c_vect)) then
      !    allocate(c_vect(n_configuration), stat=errnum)
      !    call check_allocation(errnum, 'c_vect (1D) in CI')
      !endif
      !if (.not.allocated(w_vect)) then
      !    !allocate(w_vect(n_configuration+1), stat=errnum)
      !    allocate(w_vect(n_configuration), stat=errnum)
      !    call check_allocation(errnum, 'w_vect (1D) in CI')
      !endif

      c_vect = 0.0d0
      w_vect = 0.0d0
      c_0    = 1.0d0
      w_0    = 0.0d0
      norm_A = 1.0d0

    case(4)
      if (myid.eq.0) then
          write(use_unit,'(2X,A,A)') &
              "CISDTQ : CI algorithm with ", &
              "single, double, triple and quadruple excitations"
      endif
    case default
       call aims_stop('Assertion: Unknown ci_trunction', 'ci_calculation')
    end select
    
    !===============================================================
    ! Now do CI iteration
    !===============================================================
    ! initialize the correlation energy based on the HF ground state
    ! with E_corr = 0.0d0
    E_ci         = E_HF
    E_corr       = 0.0d0
    dE_ci        = 1.0d0
    dE_ci_prev   = 1.0d0
    i_scf        = 0

    ! initialization for the MP2 energy
    w_0    = 0.0d0
    w_vect = 0.0d0
    c_vect = 0.0d0
    c_0    = 1.0d0
    call evaluate_w_vector_1()
    call evaluate_w_vector_1a1b()
    call evaluate_w_vector_2a()
    if (n_elec(1) .ne. n_elec(2) .or. n_spin .eq. 2) then
        call evaluate_w_vector_2b()
    else
        !i_start = n_configuration-n_configuration_2b+1
        do i_state = memory_table(1,2,5,myid+1), memory_table(2,2,5,myid+1),1
          w_vect(i_state) = w_vect(i_state - memory_table(3,2,4,myid+1))
        enddo
    endif
    if (flag_single_excitation) then
        call evaluate_w_vector_1a()
        if (n_elec(1) .ne. n_elec(2) .or. n_spin .eq. 2) then
            call evaluate_w_vector_1b()
        else
            !i_start = n_configuration_1a + 1
            !do i_state = i_start, n_configuration_1a+n_configuration_1b ,1
            do i_state = memory_table(1,2,2,myid+1), memory_table(2,2,2,myid+1) ,1
              w_vect(i_state) = w_vect(i_state - memory_table(3,2,1,myid+1))
            enddo
        endif
    endif

    ! prepare the first-order wavefunction "c_vect"
    ! which is the same as that of the MP2 calculation

    flag_not_initialization = .false.
    call evaluate_c_vector()

    if ((ci_acc_method.eq.2).or.(ci_acc_method.eq.3)) then
      ! CI acceleration (Pople's method
      ! save wave-function (Cvect)
      ci_acc_save_c0(1)=1.0d0
      ci_acc_save_cvect(1,:)=0.0d0

      n_dim_ci_acc_matrix = 1
      ci_acc_end = 1
      Call ci_acc_save_itgl(1)

    end if



    if (myid .eq. 0) write(use_unit,'(2X,A,F16.8)') '| Norm(A) =', norm_A
    flag_not_initialization = .true.
    !if (myid .eq. 0) then
    !    write(use_unit,*) "igor debug c_vect_0   ", &
    !        c_0
    !    write(use_unit,*) "igor debug c_vect_1a  ", &
    !        c_vect(memory_table(1,2,1,1)), &
    !        memory_table(1,1,1,1)
    !    write(use_unit,*) "igor debug c_vect_1b  ", &
    !        c_vect(memory_table(1,2,2,1)), &
    !        memory_table(1,1,2,1)
    !    write(use_unit,*) "igor debug c_vect_1a1b", &
    !        c_vect(memory_table(1,2,3,1)), &
    !        memory_table(1,1,3,1)
    !    write(use_unit,*) "igor debug c_vect_2a  ", &
    !        c_vect(memory_table(1,2,4,1)), &
    !        memory_table(1,1,4,1)
    !    write(use_unit,*) "igor debug c_vect_2b  ", &
    !        c_vect(memory_table(1,2,5,1)), &
    !        memory_table(1,1,5,1)
    !endif


    ! =============================================
    ! calculate the MP2 energy based on given c_vect
    ! =============================================

    E_corr_MP2 = 0.0d0 
    !do i_state = n_configuration_1a+n_configuration_1b+1, n_configuration, 1
    !  E_corr_MP2 = E_corr_MP2 + &
    !      w_vect(i_state+1)*c_vect(i_state)
    !enddo
    do i_state = memory_table(1,2,3,myid+1), memory_table(2,2,5,myid+1), 1
      E_corr_MP2 = E_corr_MP2 + &
          w_vect(i_state)*c_vect(i_state)
    enddo
    call sync_real_number(E_corr_MP2)
    ! Just time the normalization factor, because w_vect used here is the
    ! initialized one
    E_corr_MP2 = E_corr_MP2 * norm_A


    ! ===========================================================
    ! calculate the MP3 and 2th CI energies based on given c_vect
    ! ===========================================================

    w_vect = 0.0d0
    call evaluate_w_vector_1()

    ! timing
    call get_timestamps(time_w_1a1b, clock_time_w_1a1b)

    call evaluate_w_vector_1a1b()
    
    ! timing
    call get_timestamps(time_w_2a, clock_time_w_2a)
    tot_time_w_1a1b       = time_w_2a       - time_w_1a1b
    tot_clock_time_w_1a1b = clock_time_w_2a - clock_time_w_1a1b

    call evaluate_w_vector_2a()

    ! timing
    call get_timestamps(tot_time_w_2a, tot_clock_time_w_2a)
    tot_time_w_2a         = tot_time_w_2a       - time_w_2a
    tot_clock_time_w_2a   = tot_clock_time_w_2a - clock_time_w_2a

    if (n_elec(1) .ne. n_elec(2) .or. n_spin .eq. 2) then
        call evaluate_w_vector_2b()
    else
        do i_state = memory_table(1,2,5,myid+1), memory_table(2,2,5,myid+1),1
          w_vect(i_state) = w_vect(i_state - memory_table(3,2,4,myid+1))
        enddo
    endif
    if (flag_single_excitation) then
        ! timing
        call get_timestamps(time_w_1a, clock_time_w_1a)

        call evaluate_w_vector_1a()

        ! timing
        call get_timestamps(tot_time_w_1a, tot_clock_time_w_1a)
        tot_time_w_1a         = tot_time_w_1a       - time_w_1a
        tot_clock_time_w_1a   = tot_clock_time_w_1a - clock_time_w_1a

        if (n_elec(1) .ne. n_elec(2) .or. n_spin .eq. 2) then
            call evaluate_w_vector_1b()
        else
            !i_start = n_configuration_1a + 1
            !do i_state = i_start, n_configuration_1a+n_configuration_1b ,1
            !  w_vect(i_state) = w_vect(i_state - n_configuration_1a)
            do i_state = memory_table(1,2,2,myid+1), memory_table(2,2,2,myid+1) ,1
              w_vect(i_state) = w_vect(i_state - memory_table(3,2,1,myid+1))
            enddo
        endif
    endif


    !if (myid .eq. 0) then
    !    write(use_unit,*) "igor debug w_vect_0   ", &
    !        w_0
    !    write(use_unit,*) "igor debug w_vect_1a  ", &
    !        w_vect(memory_table(1,2,1,1)), &
    !        memory_table(1,1,1,1)
    !    write(use_unit,*) "igor debug w_vect_1b  ", &
    !        w_vect(memory_table(1,2,2,1)), &
    !        memory_table(1,1,2,1)
    !    write(use_unit,*) "igor debug w_vect_1a1b", &
    !        w_vect(memory_table(1,2,3,1)), &
    !        memory_table(1,1,3,1)
    !    write(use_unit,*) "igor debug w_vect_2a  ", &
    !        w_vect(memory_table(1,2,4,1)), &
    !        memory_table(1,1,4,1)
    !    write(use_unit,*) "igor debug w_vect_2b  ", &
    !        w_vect(memory_table(1,2,5,1)), &
    !        memory_table(1,1,5,1)
    !endif

    call get_timestamps(time_E_ci, clock_time_E_ci)

    call evaluate_E_ci()

    call get_timestamps(tot_time_E_ci, tot_clock_time_E_ci)
    tot_time_E_ci       = tot_time_E_ci - time_E_ci
    tot_clock_time_E_ci = tot_clock_time_E_ci - clock_time_E_ci

    if (myid .eq. 0) then
        E_corr_MP3 = (E_ci - E_HF) * norm_A**2 - E_corr_MP2
        write(use_unit,'(2X,A,E16.8,A,E16.8,A)') '| Ec[MP2] =', &
            E_corr_MP2, ' Ha,   E[MP2] =', &
            E_HF + E_corr_MP2+en_ion_ion,' Ha'
        write(use_unit,'(2X,A,E16.8,A,E16.8,A)') '| Ec[MP3] =', &
            E_corr_MP3, &
            ' Ha,   E[MP3] =',E_HF+E_corr_MP2+E_corr_MP3+en_ion_ion,' Ha'
        write(use_unit,'(2X,A,E16.8,A,E16.8,A)') '| Ec[CI]  =', &
            E_ci - E_HF, &
            ' Ha,   E[CI]  =',E_ci+en_ion_ion,' Ha'
    endif

    !if (myid .eq. 0) then
    !    write(use_unit,'(2X,A,I5)') '| The first-order c_vect for MP2', &
    !        n_configuration
    !    write(use_unit,'(4X,A,4(A5),A21)') &
    !        '| Spin Case ', 'I','J','A','B','Value     '
    !endif
    !call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!
    !call output_w_vect(1.0d0,ci_coeff_threshold)
    !call output_CI_amplitudes(1.0d0/c_0,ci_coeff_threshold)
    !call output_CI_amplitudes_v02(1.0d0/c_0,ci_coeff_threshold)

    call get_timestamps(time_c_vect, clock_time_c_vect)

    if ((ci_acc_method.eq.2).or.(ci_acc_method.eq.3)) then

!      Call ci_acc_save_itgl(2)

      ci_acc_start = 0
      Call ci_acc_orthogonal(flag_ci_acc_S_singular)
      if (.not.(flag_ci_acc_S_singular)) then
        if (n_dim_ci_acc_matrix.lt.n_ci_acc_save) then
          n_dim_ci_acc_matrix = n_dim_ci_acc_matrix + 1
        end if
        ci_acc_start = 1
      end if

!      Call ci_acc_optwvfn()
!      n_dim_ci_acc_matrix=2
!      ci_acc_start=1
!      ci_acc_end=2

!      Call ci_acc_optwvfn()
      call evaluate_E_ci()
    end if

    call evaluate_c_vector()

    call get_timestamps(tot_time_c_vect, tot_clock_time_c_vect)
    tot_time_c_vect       = tot_time_c_vect - time_c_vect
    tot_clock_time_c_vect = tot_clock_time_c_vect - clock_time_c_vect
    if (myid .eq. 0) then
       write(use_unit,'(2X,A,f16.3,A,f16.3,A)') &
           "| Evaluate w_vect(1a1b) time  : ", tot_time_w_1a1b, " s ", tot_clock_time_w_1a1b, " s. "
       write(use_unit,'(2X,A,f16.3,A,f16.3,A)') &
           "| Evaluate w_vect(2a) time    : ", tot_time_w_2a, " s ", tot_clock_time_w_2a, " s. "
       write(use_unit,'(2X,A,f16.3,A,f16.3,A)') &
           "| Evaluate E_ci time          : ", tot_time_E_ci, " s ", tot_clock_time_E_ci, " s. "
       write(use_unit,'(2X,A,f16.3,A,f16.3,A)') &
           "| Evaluate c_vect time        : ", tot_time_c_vect, " s ", tot_clock_time_c_vect, " s. "
    endif

    if (myid .eq. 0) write(use_unit,'(2X,A,F16.8)') '| Norm(A) =', norm_A

    !if (myid .eq. 0) then
    !    write(use_unit,'(2X,A,I5)') '| The first-order c_vect', &
    !        n_configuration
    !    call output_array(n_configuration,c_vect)
    !endif

    dE_ci = E_ci - E_HF
    dE_ci_prev = dE_ci
    E_ci_prev = E_ci

    if (myid .eq. 0) then
        write(use_unit,'(2X,A,I3,A,E16.8,A)') &
            "| The CISD correlation energy in ",i_scf,"th iteration is", &
            E_ci - E_HF, ' Ha'
    endif

    SCF_LOOP: do while (abs(dE_ci) .gt. ci_accuracy_etot .and. &
       i_scf .lt. n_scf_ci)

        i_scf = i_scf + 1

        ! generate w_vect
        w_vect = 0.0d0
        !w_vect(1) = ddot(n_configuration, c_vect(:), 1, b_vect(:), 1)
        call evaluate_w_vector_1()
        call evaluate_w_vector_1a1b()
        call evaluate_w_vector_2a()
        if (n_elec(1) .ne. n_elec(2) .or. n_spin .eq. 2) then
            call evaluate_w_vector_2b()
        else
            !i_start = n_configuration-n_configuration_2b+1
            !do i_state = i_start, n_configuration,1
            !  w_vect(i_state) = w_vect(i_state - n_configuration_2a)
            do i_state = memory_table(1,2,5,myid+1), memory_table(2,2,5,myid+1),1
              w_vect(i_state) = w_vect(i_state - memory_table(3,2,4,myid+1))
            enddo
        endif
        if (flag_single_excitation) then
            call evaluate_w_vector_1a()
            if (n_elec(1) .ne. n_elec(2) .or. n_spin .eq. 2) then
                call evaluate_w_vector_1b()
            else
                !i_start = n_configuration_1a + 1
                !do i_state = i_start, n_configuration_1a+n_configuration_1b ,1
                !  w_vect(i_state) = w_vect(i_state - n_configuration_1a)
                do i_state = memory_table(1,2,2,myid+1), memory_table(2,2,2,myid+1) ,1
                  w_vect(i_state) = w_vect(i_state - memory_table(3,2,1,myid+1))
                enddo
            endif
        endif

        call evaluate_E_ci()
        if (myid.eq.0) then
          write(use_unit,*) 'E_ci = ',E_ci+en_ion_ion
        end if

        if ((ci_acc_method.eq.2).or.(ci_acc_method.eq.3)) then

          Call ci_acc_orthogonal(flag_ci_acc_S_singular)

          if (.not.(flag_ci_acc_S_singular)) then
            if (n_dim_ci_acc_matrix.lt.n_ci_acc_save) then
              n_dim_ci_acc_matrix = n_dim_ci_acc_matrix + 1
            end if
          end if

          Call ci_acc_optwvfn()
        end if

        ! calculate CI energy
        call evaluate_E_ci()
        if (myid.eq.0) then
          write(use_unit,*) 'E_ci ',E_ci+en_ion_ion
        end if
        ! generate c_vect for the next iteration
        if (ci_acc_method.eq.3) then
          call ci_acc_evaluate_c()
        else
          call evaluate_c_vector()
        end if


        if (myid .eq. 0) then
            write(use_unit,'(2X,A,F16.8)') '| Norm(A) =', norm_A
            !write(use_unit,'(2X,A,I5)') '| The first-order c_vect', &
            !    n_configuration
            !call output_array(n_configuration,c_vect)
        endif

        dE_ci_prev = dE_ci
        dE_ci      = E_ci - E_ci_prev

        !if (dE_ci*dE_ci_prev .lt. 0.0d0) then ! to damp the oscillation in truncated CI calculations
        !    E_ci = 0.5d0*(E_ci + E_ci_prev)
        !    dE_ci = E_ci - E_ci_prev
        !endif

        E_ci_prev = E_ci

        if (myid .eq. 0) then
            write(use_unit,'(2X,A,I3,A,E16.8,A,f16.8,A,E16.8,A)') &
                "| Ec[CI] in the ",i_scf,"th iteration =", E_ci - E_HF, &
                ' Ha. E[CI] =',E_ci+en_ion_ion,' Ha, with dEc[CI] =',dE_ci,' Ha'
        endif
    end do SCF_LOOP

    !if (myid .eq. 0) then
    !    write(use_unit,'(2X,A,I5)') '| Final c_vect with the dimension of', &
    !        n_configuration
    !    write(use_unit,'(4X,A,4(A5),A21)') &
    !        '| Spin Case ', 'I','J','A','B','Value     '
    !endif
    !call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!
    !call output_CI_amplitudes(1.0d0/c_0,ci_coeff_threshold)
    !call output_CI_amplitudes_v02(1.0d0/c_0,ci_coeff_threshold)

888 continue

    ! for timing
    call get_timestamps(tot_time_CI, tot_clock_time_CI)
    tot_time_CI       = tot_time_CI       - time_CI
    tot_clock_time_CI = tot_clock_time_CI - clock_time_CI
    if(myid.eq.0) then
       if (abs(dE_ci) .le. ci_accuracy_etot) then
         write(use_unit,"(2X,A53,2x,I4,2x,A11)") &
            'Configuration Interaction SCF procedure converged in',i_scf,'iterations.'
         write(use_unit,"(2X,A,F16.8,A)") 'E(CI) = ',E_ci+en_ion_ion,' Ha'
       else
         write(use_unit,"(2X,A50,2x,I4,A11)") &
            'Configuration Interaction SCF procedure failed in',i_scf,'iterations'
       end if

       write(use_unit,'(2X,A,A19,A20)') &
           "Timing analysis               : ", "max(cpu_time) ", &
           "wall_clock(cpu1) "
       write(use_unit,'(2X,A,f16.3,A,f16.3,A)') &
           "| Total CI time               : ", tot_time_CI, " s ", tot_clock_time_CI, " s. "
       if (flag_single_excitation) then
           write(use_unit,'(2X,A,f16.3,A,f16.3,A)') &
               "| Evaluate w_vect(1a) time    : ", tot_time_w_1a, " s ", tot_clock_time_w_1a, " s. "
       endif
       write(use_unit,'(2X,A,f16.3,A,f16.3,A)') &
           "| Evaluate w_vect(1a1b) time  : ", tot_time_w_1a1b, " s ", tot_clock_time_w_1a1b, " s. "
       write(use_unit,'(2X,A,f16.3,A,f16.3,A)') &
           "| Evaluate w_vect(2a) time    : ", tot_time_w_2a, " s ", tot_clock_time_w_2a, " s. "
    endif
    
    call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!
    ! clean up the memory for the CI calculation
    call cleanup_ci()

    end subroutine ci_calculation
!****** 

