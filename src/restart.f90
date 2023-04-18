!****h* FHI-aims/restart
!  NAME
!    restart
!  SYNOPSIS
module restart
!  PURPOSE
!     Module to house restart options for the scf solver
!     uses the firsty iteration of an SCF cycle to compute the last
!     known density.
!     If we use scalapack and the density matrix we restart more efficently from the density.
!     That implementation is linked in update_density_densmat.f90.
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
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
  implicit none

  logical :: restart_zero_iteration
  logical :: restart_first_warnings = .true.

!******
CONTAINS


!------------------------------------------------------------------------------
!****s* restart/determine_occupied_states
!  NAME
!    determine_occupied_states
!  SYNOPSIS
SUBROUTINE determine_occupied_states(n_states_occupied)
!  PURPOSE
!    Determines the highest occupied state for all k points of this task.
!  USES
    use dimensions, only: n_states, n_spin, n_k_points_task
    use physics, only: occ_numbers

    implicit none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    The highest occupied state for all k points of this task.
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

    ! Arguments
    INTEGER, INTENT(OUT) :: n_states_occupied

    ! Local variable
    INTEGER :: i_states, i_spin

    do i_states = n_states, 1, -1
        n_states_occupied = i_states

        ! First not non occupied state.
        if(ANY(occ_numbers(i_states, 1:n_spin, 1:n_k_points_task) > 1.d-5)) then
            exit
        end if
    enddo

END SUBROUTINE determine_occupied_states
!******

!------------------------------------------------------------------------------
!****s* restart/write_restart_info_cluster
!  NAME
!    write_restart_info_cluster
!  SYNOPSIS
  subroutine write_restart_info_cluster
!  PURPOSE
!    Write the restart information for a cluster to a file.
!  USES
    use mpi_tasks
    use localorb_io
    use runtime_choices, only: restart_write_file, fo_dft, fo_info_file, restart_for_postscf
    use dimensions, only: n_basis, n_states, n_spin
    !#use physics, only: KS_eigenvector, KS_eigenvalue, occ_numbers
    use physics

    implicit none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

    ! Local variables
    INTEGER :: i_basis, i_states, i_spin
    INTEGER :: n_states_occupied, n_states_to_save
    character*150 :: info_str

    call reshape_matrices_n_states(n_states)

    ! Only CPU 0 writes information.
    if (myid.ne.0) then
        return
    endif

    ! Write message to console.
    write(info_str,'(2X,2A)') 'Writing cluster restart information to file ', restart_write_file
    call localorb_info(info_str,use_unit,'(A)')

    ! Open the restart file for writing.
    open(file = restart_write_file, unit = 7, status = 'unknown', form = 'unformatted')

    ! Determine used states
    CALL determine_occupied_states(n_states_occupied)

    ! Save all occupied states and if possible a few additional.
    if (fo_dft .or. restart_for_postscf) then
        n_states_to_save = n_states
    else
        !n_states_to_save = min(n_states, n_states_occupied + 3)
        n_states_to_save = n_states
    endif

    write(info_str,'(2X,A,I9,A)') 'Restart file: Writing ', n_states_to_save, ' states.'
    call localorb_info(info_str,use_unit,'(A)')

    if (fo_dft) then
        open(file = fo_info_file, unit = 25, status = 'unknown', form = 'formatted')
        write(25,*) n_states_occupied
        write(25,*) n_states - n_states_occupied
        close(unit=25)
    endif

    ! Write header information
    write(7) n_basis
    write(7) n_states_to_save
    write(7) n_spin

    ! Write KS_eigenvector(n_basis,n_states,n_spin,n_k_points_task)
    do i_basis = 1, n_basis
        do i_states = 1, n_states_to_save
            do i_spin = 1, n_spin
                write(7) KS_eigenvector(i_basis,i_states,i_spin,1)
            end do
        end do
    end do

    ! Write eigenvalues and occupation numbers.
    do i_states = 1, n_states_to_save
        do i_spin = 1, n_spin
            write(7) KS_eigenvalue(i_states,i_spin,1), &
                          occ_numbers(i_states,i_spin,1)
        end do
    end do

    ! Close the file.
    close(unit = 7)

  end subroutine write_restart_info_cluster
!******

!------------------------------------------------------------------------------
!****s* restart/write_restart_info_periodic
!  NAME
!    write_restart_info_periodic
!  SYNOPSIS
  subroutine write_restart_info_periodic
!  PURPOSE
!    Write restart information for periodic systems to file.
!  USES
    use localorb_io
    use runtime_choices, only: restart_write_file, real_eigenvectors, force_single_restartfile, &
        restart_eigenvectors_periodic, restart_for_postscf
    use dimensions, only: n_basis, n_states, n_spin, n_k_points, n_k_points_task
    use physics, only: KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers
    use scalapack_wrapper, only: my_k_point, my_scalapack_id, n_scalapack_tasks
    use pbc_lists, only: k_point_list
    use mpi_tasks, only: myid, n_tasks

    implicit none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

    ! Local variables
    INTEGER :: i_basis, i_states, i_spin, i_k
    INTEGER :: n_states_occupied, n_states_to_save
    character*150 :: info_str
    character*40 :: restart_write_file_evp

    if (force_single_restartfile) then
        ! avoid error for serial runs
        if (n_tasks > 1) then
            ! CPU 1 has Gamma-K point, all other just write basic info which can be skipped
            ! Also important to avoid race conditions while writing since for force_single_restartfile
            ! we only write one file.
            if (myid.ne.1) then
                return
            endif
        endif
    endif

    ! Write message to console
    write(info_str,'(2X,2A)') 'Writing periodic restart information to file ', restart_write_file
    call localorb_info(info_str,use_unit,'(A)')


    if (restart_eigenvectors_periodic) then

        if (n_k_points.gt.n_tasks) then
            write(info_str,'(A)') '***  WARNING: restart_eigenvectors_periodic can currently not be used'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(A)') 'if n_k_points > n_tasks! Therefore, no restart files will be created!'
            call localorb_info(info_str,use_unit,'(A)')
            return
        elseif (n_scalapack_tasks.eq.0.and.n_tasks.gt.n_k_points.and. &
            myid.gt.0.and.myid.le.n_k_points) then
            my_k_point = myid
        elseif (n_scalapack_tasks.eq.0.and.n_tasks.gt.n_k_points.and. &
            (myid.eq.0.or.myid.gt.n_k_points)) then
            return
        elseif (n_scalapack_tasks.eq.0.and.n_tasks.eq.n_k_points.and.myid.gt.0) then
            my_k_point = myid
        elseif (n_scalapack_tasks.eq.0.and.n_tasks.eq.n_k_points.and.myid.eq.0) then
            my_k_point = n_k_points
        endif


        if (my_scalapack_id.eq.0) then

            ! Create dummy restart file
            open(file = restart_write_file, unit = 7, status = 'unknown', form = 'unformatted')
            close (unit = 7)

            ! Define restart file name for each k-point
            if (my_k_point.ge.1.and.my_k_point.lt. 10) then
                write(restart_write_file_evp,'(A,I1)') trim(restart_write_file), my_k_point
            else if (my_k_point.ge.10.and.my_k_point.lt.100) then
                write(restart_write_file_evp,'(A,I2)') trim(restart_write_file), my_k_point
            else if (my_k_point.ge.100.and.my_k_point.lt.1000) then
                write(restart_write_file_evp,'(A,I3)') trim(restart_write_file), my_k_point
            else if (my_k_point.ge.1000.and.my_k_point.lt.10000) then
                write(restart_write_file_evp,'(A,I4)') trim(restart_write_file), my_k_point
            else if (my_k_point.ge.10000.and.my_k_point.lt.100000) then
                write(restart_write_file_evp,'(A,I5)') trim(restart_write_file), my_k_point
            else if (my_k_point.ge.100000.and.my_k_point.lt.1000000) then
                write(restart_write_file_evp,'(A,I6)') trim(restart_write_file), my_k_point
            endif

        else
            return
        endif

    endif

    ! Open the restart file for writing
    if (restart_eigenvectors_periodic) then
        open(file = restart_write_file_evp, unit = 7, status = 'unknown', form = 'unformatted')
    else
        open(file = restart_write_file, unit = 7, status = 'unknown', form = 'unformatted')
    endif

    ! Determine used states
    CALL determine_occupied_states(n_states_occupied)

    ! Save all occupied states and if possible a few additional
    if (restart_for_postscf) then
        n_states_to_save = n_states
    else
        n_states_to_save = min(n_states, n_states_occupied + 3)
    end if

    ! Write header information
    if (restart_eigenvectors_periodic) then
        write(7) k_point_list(my_k_point,1), k_point_list(my_k_point,2), k_point_list(my_k_point,3)
    endif
    write(7) n_basis
    write(7) n_states_to_save
    write(7) n_spin
    write(7) n_k_points_task

    ! Write KS_eigenvectors. (real or complex)
    do i_k = 1, n_k_points_task
        do i_spin = 1, n_spin
            do i_states = 1, n_states_to_save
                do i_basis = 1, n_basis
                    if(real_eigenvectors)then
                        write(7) KS_eigenvector(i_basis,i_states,i_spin,i_k)
                    else ! complex eigenvector
                        write(7) KS_eigenvector_complex(i_basis,i_states,i_spin,i_k)
                    endif
                end do ! i_basis
            end do ! i_states
        end do ! i_spin
    end do ! i_k

    ! Write KS_eigenvalues occupation numbers
    do i_k = 1,n_k_points
        do i_spin = 1, n_spin
            do i_states = 1, n_states_to_save
                write(7) KS_eigenvalue(i_states,i_spin,i_k), &
                         occ_numbers(i_states,i_spin,i_k)
            end do ! i_states
        end do ! i_spin
    end do ! i_k

    ! Close the restart file
    close(unit = 7)

  end subroutine write_restart_info_periodic
!******

!------------------------------------------------------------------------------
!****s* restart/write_restart_info
!  NAME
!    write_restart_info
!  SYNOPSIS
  subroutine write_restart_info
!  PURPOSE
!    write restart information to file
!  USES
    use runtime_choices, only: use_scalapack, use_density_matrix, restart_write, &
            fo_dft, fo_info_file, force_single_restartfile, restart_eigenvectors_periodic
    use dimensions, only: n_periodic

    implicit none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
    ! If we use scalapack and the density matrix we restart more efficently from the density.
    ! That implementation is linked in update_density_densmat.f90.
    if((use_scalapack .and. use_density_matrix).and..not.(force_single_restartfile.or.fo_dft &
        .or.restart_eigenvectors_periodic)) then
        return
    end if

    if (restart_write) then
        if (n_periodic.eq.0) then
            CALL write_restart_info_cluster()
        else
            CALL write_restart_info_periodic()
        end if

        restart_first_warnings = .false.

    end if

  end subroutine write_restart_info
!******

!------------------------------------------------------------------------------
!****s* restart/invalid_restart_files
!  NAME
!    invalid_restart_files
!  SYNOPSIS
  SUBROUTINE invalid_restart_files
!  PURPOSE
!    Writes an error message and stops aims because of invalid restart files.
!  USES
    use localorb_io
    use mpi_utilities
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
    implicit none
    ! Local variables
    CHARACTER*150 :: info_str

    ! Write the error message.
    write(info_str,'(2X,A)') '* WARNING: restart file cannot belong to current system.'
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') '*          Please check'
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') '*          Aborting execution'
    call localorb_info(info_str,use_unit,'(A)')

    ! Stop aims.
    call aims_stop

END SUBROUTINE invalid_restart_files
!******

!------------------------------------------------------------------------------
!****s* restart/check_restart_info_cluster_validity
!  NAME
!    check_restart_info_cluster_validity
!  SYNOPSIS
  SUBROUTINE check_restart_info_cluster_validity(n_basis_file,  &
                                                 n_states_file, &
                                                 n_spin_file)
!  PURPOSE
!    Checks the written cluster restart file header for validity.
!  USES
    use dimensions, only: n_basis, n_states, n_spin
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    The basis set size the number of saved states and the number of spin channels saved in the
!    file header.
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
    implicit none
    ! Arguments
    INTEGER,INTENT(IN) :: n_basis_file
    INTEGER,INTENT(IN) :: n_states_file
    INTEGER,INTENT(IN) :: n_spin_file

    ! Local variable
    LOGICAL :: is_restart_file_valid

    ! A restart file is considered valid if: the basis is of the same size, the same number of
    ! spin channels are used.
    is_restart_file_valid = (n_basis  .eq. n_basis_file) &
                            .AND. &
                            (n_states .ge. n_states_file) &
                            .AND. &
                            (n_spin   .eq. n_spin_file)

    ! If restart file is invalid.
    if (.not.is_restart_file_valid) then
        CALL invalid_restart_files()
    end if

END SUBROUTINE check_restart_info_cluster_validity
!******

!------------------------------------------------------------------------------
!****s* restart/read_restart_info_cluster
!  NAME
!    read_restart_info_cluster
!  SYNOPSIS
  SUBROUTINE read_restart_info_cluster
!  PURPOSE
!    Read SCF information from last known scf cycle for clusters from files.
!  USES
    use runtime_choices, only: restart_read_file
    use physics, only: KS_eigenvector, KS_eigenvalue, occ_numbers
    use dimensions,only: n_basis, n_states, n_spin, force_occupation_projector, &
                         start_force_occ, n_force_occ
    use force_occupation,only: previous_eigenvector, &
                               allocate_previous_eigenvector, &
                               force_occ_pr_state, forced_occ_number, &
                               force_occ_spin
    use synchronize_mpi
    use localorb_io, only: localorb_info, use_unit
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
    implicit none

    ! Local variables
    INTEGER :: n_basis_file, n_states_file, n_spin_file
    INTEGER :: i_basis, i_spin, i_states, i_force_occ
    CHARACTER*150 :: info_str

    if (force_occupation_projector) then
      call allocate_previous_eigenvector
    end if

    if (myid.eq.0) then

        ! Write message to tell reading cluster restart file.
        write(info_str,'(2X,A)') ''
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,2A)') 'Reading cluster restart information from file ', restart_read_file
        call localorb_info(info_str,use_unit,'(A)')

        ! Open restart file
        open(file = restart_read_file, unit = 7, status = 'old', form = 'unformatted')

        ! Read header information
        read(7) n_basis_file
        read(7) n_states_file
        read(7) n_spin_file

        ! Check for restart file validity
        CALL check_restart_info_cluster_validity(n_basis_file,  &
                                                 n_states_file, &
                                                 n_spin_file)

        ! Read in eigenvectors
        ! If fewer eigenvectors are saved than used in this calculation just fill them up with zeros
        do i_basis = 1, n_basis
            do i_states = 1, n_states
                do i_spin = 1, n_spin
                    if (i_states <= n_states_file) then
                        read(7) KS_eigenvector(i_basis,i_states,i_spin,1)
                    else
                        KS_eigenvector(i_basis,i_states,i_spin,1) = 0.d0
                    endif
                end do
            end do
        end do

        ! Set the previous eigenvector for the force occupation projector
        if (force_occupation_projector) then
            previous_eigenvector(:,:,:,:) = KS_eigenvector(:,:,:,:)
            start_force_occ = .true.
        end if

        ! Read in the eigenvalues and the occupation numbers.
        ! If fewer eigenvales and occupation numbers are saved than used in this calculation just fill them up with zero.
        do i_states = 1, n_states
            do i_spin = 1, n_spin
                if (i_states <= n_states_file) then
                    read(7) KS_eigenvalue(i_states,i_spin,1), &
                            occ_numbers(i_states,i_spin,1)
                else
                    KS_eigenvalue(i_states,i_spin,1) = 0.d0
                    occ_numbers(i_states,i_spin,1)   = 0.d0
                end if
            end do
        end do

        ! Set the previous occupations for the force occupation projector
        if (force_occupation_projector) then
            do i_force_occ = 1, n_force_occ
                occ_numbers(force_occ_pr_state(i_force_occ), &
                            force_occ_spin(i_force_occ), 1) = &
                            forced_occ_number(i_force_occ)
            end do
        end if

        ! Close the restart file.
        close(unit = 7)
    else
        ! All other CPUs not equal 0
        KS_eigenvector = 0.d0
        KS_eigenvalue  = 0.d0
        occ_numbers    = 0.d0
    end if

    ! MPI broadcast to all other threads
    if (force_occupation_projector) then
            call sync_restart(KS_eigenvector, KS_eigenvalue, occ_numbers, &
                              n_basis, n_states, n_spin, start_force_occ, &
                              previous_eigenvector)
    else
            call sync_restart(KS_eigenvector, KS_eigenvalue, occ_numbers, &
                              n_basis, n_states, n_spin)
    end if

  END SUBROUTINE read_restart_info_cluster
!******

!------------------------------------------------------------------------------
!****s* restart/check_restart_info_periodic_validity
!  NAME
!    check_restart_info_periodic_validity
!  SYNOPSIS
  SUBROUTINE check_restart_info_periodic_validity(n_basis_file,  &
                                                  n_states_file, &
                                                  n_spin_file, &
                                                  n_k_points_task_file)
!  PURPOSE
!    Checks the written restart file header for periodic systems for validity.
!  USES
    use dimensions,only : n_basis, n_states, n_spin, n_k_points_task
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    The basis set size the number of saved states, the number of spin channels, and the number
!    of k points saved in the file header.
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
    implicit none
    ! Arguments
    INTEGER,INTENT(IN) :: n_basis_file
    INTEGER,INTENT(IN) :: n_states_file
    INTEGER,INTENT(IN) :: n_spin_file
    INTEGER,INTENT(IN) :: n_k_points_task_file

    ! Local variable
    LOGICAL :: is_restart_file_valid

    ! A restart file is considered valid if: the basis is of the same size, the same number of spin channels are used.
    is_restart_file_valid = (n_basis         .eq. n_basis_file)  &
                            .AND. &
                            (n_states        .ge. n_states_file) &
                            .AND. &
                            (n_spin          .eq. n_spin_file)   &
                            .AND. &
                            (n_k_points_task .eq. n_k_points_task_file)


    ! If restart file is invalid.
    if (.not.is_restart_file_valid) then
        CALL invalid_restart_files()
    end if

END SUBROUTINE check_restart_info_periodic_validity
!------------------------------------------------------------------------------
!****s* restart/read_restart_info_periodic
!  NAME
!    read_restart_info_periodic
!  SYNOPSIS
  SUBROUTINE read_restart_info_periodic
!  PURPOSE
!    Read SCF information from last known scf cycle for peridodic systems from files.
!  USES
    use runtime_choices, only: real_eigenvectors, restart_read_file, use_scalapack
    use synchronize_mpi
    use physics, only: KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers
    use dimensions,only: n_basis, n_states, n_spin, n_k_points, n_k_points_task, force_occupation_projector, &
                         start_force_occ, n_force_occ
    use force_occupation,only: previous_eigenvector, previous_eigenvector_complex, allocate_previous_eigenvector, &
                               force_occ_pr_state, forced_occ_number, &
                               force_occ_spin
    use runtime_choices, only: force_single_restartfile, restart_eigenvectors_periodic
    use scalapack_wrapper, only: my_k_point,  my_scalapack_id, n_scalapack_tasks
    use pbc_lists, only: k_point_list
    use localorb_io, only: localorb_info, use_unit

!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
    implicit none

    ! Local variables
    INTEGER :: n_basis_file, n_states_file, n_spin_file, n_k_points_task_file
    INTEGER :: i_basis, i_spin, i_states, i_k, i_force_occ
    CHARACTER*150 :: info_str
    character*40 :: form1
    character*40 :: restart_read_file_evp
    REAL*8,dimension(3) :: my_k_point_file
    LOGICAL :: scalapack_restart_file_exists

    if (force_occupation_projector) then
      call allocate_previous_eigenvector
    end if

    if(real_eigenvectors)then
        write(form1,'(A,I3,A)') '(',n_spin,'E30.18E4)'
    else ! complex eigenvector
        write(form1,'(A,I3,A)') '(',2*n_spin,'E30.18E4)'
    end if

    ! Write message to tell reading restart file for periodic system.
    write(info_str,'(2X,A)') ''
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,2A)') 'Reading periodic restart information from file ', restart_read_file
    call localorb_info(info_str,use_unit,'(A)')

    ! in case of rotated periodic restarts only gamma point calculations
    ! are possible. Therefore, the restart files for each process
    ! are identical. To improve transferability of restart files, only
    ! one file is necessary and considered in this case
    if (force_single_restartfile.and..not.use_scalapack) then
        !! avoid error for serial runs
        if (n_tasks > 1) then
            !! CPU 1 has Gamma-K point, all other just write basic info which can be skipped
            if (myid.ne.1) then
                !! initialize to zero
                KS_eigenvector = 0.d0
                KS_eigenvalue  = 0.d0
                occ_numbers    = 0.d0
                return
            endif
        endif
    endif


    if (restart_eigenvectors_periodic) then

        if (n_k_points.gt.n_tasks) then
            call aims_stop('restart_eigenvectors_periodic can currently not be used if n_k_points > n_tasks!')
        elseif (n_scalapack_tasks.eq.0.and.n_tasks.gt.n_k_points.and. &
            myid.gt.0.and.myid.le.n_k_points) then
            my_k_point = myid
        elseif (n_scalapack_tasks.eq.0.and.n_tasks.gt.n_k_points.and. &
            (myid.eq.0.or.myid.gt.n_k_points)) then
            return
        elseif (n_scalapack_tasks.eq.0.and.n_tasks.eq.n_k_points.and.myid.gt.0) then
            my_k_point = myid
        elseif (n_scalapack_tasks.eq.0.and.n_tasks.eq.n_k_points.and.myid.eq.0) then
            my_k_point = n_k_points
        endif


        ! Define restart file name for each k-point
        if (my_k_point.ge.1.and.my_k_point.lt. 10) then
            write(restart_read_file_evp,'(A,I1)') trim(restart_read_file), my_k_point
        else if (my_k_point.ge.10.and.my_k_point.lt.100) then
            write(restart_read_file_evp,'(A,I2)') trim(restart_read_file), my_k_point
        else if (my_k_point.ge.100.and.my_k_point.lt.1000) then
            write(restart_read_file_evp,'(A,I3)') trim(restart_read_file), my_k_point
        else if (my_k_point.ge.1000.and.my_k_point.lt.10000) then
            write(restart_read_file_evp,'(A,I4)') trim(restart_read_file), my_k_point
        else if (my_k_point.ge.10000.and.my_k_point.lt.100000) then
            write(restart_read_file_evp,'(A,I5)') trim(restart_read_file), my_k_point
        else if (my_k_point.ge.100000.and.my_k_point.lt.1000000) then
            write(restart_read_file_evp,'(A,I6)') trim(restart_read_file), my_k_point
        endif

        inquire(FILE=restart_read_file_evp,EXIST=scalapack_restart_file_exists)
        if (.not.scalapack_restart_file_exists) then
            call aims_stop('The required restart files do not exist!')
        endif

    endif


    ! Open restart files
    if (restart_eigenvectors_periodic) then
        open(file = restart_read_file_evp, unit = 7, status = 'old', form = 'unformatted')
    else
        open(file = restart_read_file, unit = 7, status = 'old', form = 'unformatted')
    endif

    ! Read header information.
    if (restart_eigenvectors_periodic) then
        read(7) my_k_point_file
        if((my_k_point_file(1).ne.k_point_list(my_k_point,1)) &
            .or. (my_k_point_file(2).ne.k_point_list(my_k_point,2)) &
            .or. (my_k_point_file(3).ne.k_point_list(my_k_point,3))) then
            call aims_stop('The k-grid of the restart files does not match the current calculation!')
        endif
    endif
    read(7) n_basis_file
    read(7) n_states_file
    read(7) n_spin_file
    read(7) n_k_points_task_file

     !Check for restart file validity.
     CALL check_restart_info_periodic_validity(n_basis_file,  &
                                               n_states_file, &
                                               n_spin_file,   &
                                               n_k_points_task_file)


    ! Read in the KS_eigenvectors (real or complex)
    ! If fewer eigenvectors are saved than used in this calculation just fill them up with zeros
    do i_k = 1,n_k_points_task
        do i_spin = 1, n_spin
            do i_states = 1, n_states
                do i_basis = 1, n_basis

                    ! If information for current state are safed in the file.
                    if (i_states <= n_states_file) then
                        if(real_eigenvectors)then
                            read(7) KS_eigenvector(i_basis,i_states,i_spin,i_k)
                        else ! complex eigenvectors
                            read(7) KS_eigenvector_complex(i_basis,i_states,i_spin,i_k)
                        end if
                    else
                        if(real_eigenvectors)then
                            KS_eigenvector(i_basis,i_states,i_spin,i_k) = 0.d0
                        else ! complex eigenvectors
                            KS_eigenvector_complex(i_basis,i_states,i_spin,i_k) = (0.d0, 0.d0)
                        end if
                    end if

                end do ! i_basis
            end do ! i_state
         end do ! i_spin
    end do ! i_k

    ! Read in the KS_eigenvalues and the occupation numbers for all k points
    ! If fewer eigenvales and occupation numbers are saved than used in this calculation just fill them up with zero.
    do i_k = 1,n_k_points
        do i_spin = 1, n_spin
            do i_states = 1, n_states

                ! If information for current state are safed in the file.
                if (i_states <= n_states_file) then
                    read(7) KS_eigenvalue(i_states, i_spin, i_k), &
                            occ_numbers(i_states, i_spin, i_k)
                else
                    KS_eigenvalue(i_states ,i_spin,i_k) = 0.d0
                    occ_numbers(i_states, i_spin, i_k)   = 0.d0
                endif

            end do ! i_states
        end do ! i_spin
    end do ! i_k

    ! Close the restart file.
    close(unit = 7)

    ! Set the previous eigenvector for the force occupation projector
    if (force_occupation_projector) then
        if(real_eigenvectors)then
            previous_eigenvector(:,:,:,:) = KS_eigenvector(:,:,:,:)
        else
            previous_eigenvector_complex(:,:,:,:) = KS_eigenvector_complex(:,:,:,:)
        end if
        start_force_occ = .true.
        do i_force_occ = 1, n_force_occ
            occ_numbers(force_occ_pr_state(i_force_occ), &
                        force_occ_spin(i_force_occ), :) = &
                        forced_occ_number(i_force_occ)
        end do
    end if

    !if (force_single_restartfile) then
       ! MPI broadcast to all other threads
       !if (force_occupation_projector) then
               !call sync_periodic_restart(KS_eigenvector, KS_eigenvalue, occ_numbers, &
                                         !n_basis, n_states, n_spin, n_k_points, start_force_occ, &
                                         !previous_eigenvector)
       !else
               !call sync_periodic_restart(KS_eigenvector, KS_eigenvalue, occ_numbers, &
                                         !n_basis, n_states, n_spin, n_k_points)
       !end if

       !if (force_occupation_projector) then
           !start_force_occ = .true.
       !end if
    !end if
  END SUBROUTINE read_restart_info_periodic
!******


!------------------------------------------------------------------------------
!****s* restart/read_restart_info
!  NAME
!    read_restart_info
!  SYNOPSIS
  subroutine read_restart_info
!  PURPOSE
!    Read scf information from last known scf cycle from file.
!  USES
    use runtime_choices, only: use_scalapack, use_density_matrix, restart_read, force_single_restartfile, &
        restart_eigenvectors_periodic
    use dimensions, only: n_periodic

!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
    implicit none

    ! If we use scalapack and the density matrix we restart more efficently from the density.
    ! That implementation is linked in update_density_densmat.f90.
    if((use_scalapack .and. use_density_matrix).and..not.(force_single_restartfile &
        .or.restart_eigenvectors_periodic)) then
        return
    end if

    if (restart_read) then

        if (n_periodic.eq.0) then
            CALL read_restart_info_cluster()
        else
            CALL read_restart_info_periodic()
        end if

    end if

  end subroutine read_restart_info
!******
  subroutine read_frozen_virtual_orbitals()
    use runtime_choices, only: fv_orbs, fv_orbs_n, fv_orbs_file, fv_filter
    !integer(kind=4) :: i_index
    !integer(kind=4), allocatable, dimension(:) :: fv_array
    integer(kind=4) :: i_fv, i_index
    character*800 :: info_str

    if (fv_filter) then
        open(file = fv_orbs_file, unit = 7, status = 'old', form = 'formatted')
        read(7,'(I8)') fv_orbs_n
        allocate(fv_orbs(fv_orbs_n),stat=i_index)
        do i_fv = 1, fv_orbs_n
          read(7,'(I8)') fv_orbs(i_fv)
        end do
        close(7)
        !write(info_str,'(2X,A)') &
        !  '| Read in virtual orbitals that need to be frozen.'
        !call localorb_info(info_str,use_unit,'(A)')
    else
        fv_orbs_n = 1
        allocate(fv_orbs(1),stat=i_index)
        fv_orbs=0
    end if
  end subroutine read_frozen_virtual_orbitals

end module restart
