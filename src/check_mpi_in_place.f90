!---------------------------------------------------------------
!****s* check_mpi_in_place
!  NAME
!    check_mpi_in_place
!  SYNOPSIS
      subroutine check_mpi_in_place( check_success )
!  PURPOSE
!    Checks the correctness of the implementation of the
!    'MPI_IN_PLACE' flag in the MPI implementation used by the
!    code. This flag is a feature that can save some memory
!    in synchronizations of arrays across MPI tasks. However, this
!    flag is found to be broken with reasonable frequency in different
!    MPI implementations. Thus, we now test if this flag is intact, and
!    if it isn't, then it will not be used. See also the
!    'use_mpi_in_place' keyword in the FHI-aims manual.        
!  USES
        ! We only use an absolute minimum of modules here - that is, mpi_tasks .
        ! All needed variables (incl. STDOUT and STDERR) may be taken from there.
        use mpi_tasks        
        implicit none
!  ARGUMENTS
      logical :: check_success
!  INPUTS
!    o none
!  OUTPUT
!    o check_success - indicates whether the 'MPI_IN_PLACE' flag works
!  AUTHOR
!    FHI-aims team
!  HISTORY
!    Volker Blum (July 2017)
!  SOURCE

    ! size of test data array
    integer :: n_data
      
    ! array that contains test data for MPI_IN_PLACE
    real*8, allocatable :: test_data(:)

    real*8, allocatable :: aux_test_data(:)
    
    logical :: aux_check_success, check_conventional_mpi
    
    integer :: i_data
    integer :: mpierr
    
    if (.not.use_mpi) return

    ! initialize
    check_success = .true.
    
    n_data = 10
    
    allocate(test_data(n_data),stat=mpierr)
    call check_allocation(mpierr,'test_data','check_mpi_in_place')

    ! seed test data array for allreduce call below
    if (myid.eq.0) then
       test_data(:) = 1.d0
    else
       test_data(:) = 0.d0
    end if

    ! Sum the test_data array over all MPI tasks
    call MPI_ALLREDUCE(MPI_IN_PLACE, &
         test_data(:), &
         n_data, &
         MPI_DOUBLE_PRECISION, &
         MPI_SUM, &
         mpi_comm_global, &
         mpierr )

    ! The value of all entries of test_data should now be 1.d0 on all MPI tasks.
    ! If that is not the case, then the MPI_IN_PLACE flag may be broken.
    ! (NB. Understood that this test can be done more elegantly. I spelled it out.)
    do i_data = 1, n_data, 1
       if (test_data(i_data).ne.1.d0) then
          check_success = .false.
       end if
    end do

    ! Synchronize check_success variable across all MPI tasks - if one is false,
    ! aux_check_success will be false. Obviously, do not use MPI_IN_PLACE .
    call MPI_ALLREDUCE(check_success, aux_check_success, 1, MPI_LOGICAL, &
         MPI_LAND, mpi_comm_global, mpierr )

    check_success = aux_check_success

    if (.not.check_success) then
       ! Not good news - MPI_IN_PLACE did not work.
       ! For good measure, check whether the intended synchronization works at all
       ! (even without MPI_IN_PLACE)

       if (myid.eq.0) then
          write(STDOUT,'(1X,A)') "* MPI_IN_PLACE does not appear to work as intended."
          write(STDOUT,'(1X,A)') "* Checking whether MPI_ALLREDUCE works at all."
       end if
       
       ! seed test data array for allreduce call below
       if (myid.eq.0) then
          test_data(:) = 1.d0
       else
          test_data(:) = 0.d0
       end if

       allocate(aux_test_data(n_data), stat=mpierr)
       call check_allocation(mpierr,'aux_test_data','check_mpi_in_place')
       aux_test_data(:) = 0.d0
       call MPI_ALLREDUCE(test_data(:), &
         aux_test_data(:), &
         n_data, &
         MPI_DOUBLE_PRECISION, &
         MPI_SUM, &
         mpi_comm_global, &
         mpierr )
       test_data(:) = aux_test_data(:)

       ! The value of all entries of test_data should now be 1.d0 on all MPI tasks.
       ! If that is not the case, then the entire MPI implementation may be broken.
       ! (NB. Understood that this test can be done more elegantly. I spelled it out.)
       check_conventional_mpi = .true.
       do i_data = 1, n_data, 1
          if (test_data(i_data).ne.1.d0) then
             check_conventional_mpi = .false.
          end if
       end do

       ! Synchronize check_conventional_mpi variable across all MPI tasks - if one is false,
       ! aux_check_success will be false. 
       call MPI_ALLREDUCE(check_conventional_mpi, aux_check_success, 1, MPI_LOGICAL, &
            MPI_LAND, mpi_comm_global, mpierr )

       check_conventional_mpi = aux_check_success

       if (.not.check_conventional_mpi) then
       ! This is unpleasant, since it would indicate a failure of a standard
       ! MPI_ALLREDUCE call. In this (unlikely) case, the code should be stopped.
          if (myid.eq.0) then
             write(STDOUT,'(1X,A)') "* It appears that a conventional MPI_ALLREDUCE call also does not work."
             write(STDOUT,'(1X,A)') "* This is unpleasant, since it indicates a broken MPI implementation."
             write(STDOUT,'(1X,A)') "* This test was performed by subroutine check_mpi_in_place. Please"
             write(STDOUT,'(1X,A)') "* verify this subroutine if you arrive here and doubt the result."
             write(STDOUT,'(1X,A)') "* Stopping the execution of the code - sorry for any inconvenience."
          end if
          call aims_stop('MPI implementation appears to be broken.','check_mpi_in_place')
       else
          if (myid.eq.0) then
             write(STDOUT,'(1X,A)') "* Without MPI_IN_PLACE, MPI_ALLREDUCE appears to work."         
          end if
       end if
          
    end if

    if (allocated(test_data)) then
      deallocate(test_data)
    end if

    if (allocated(aux_test_data)) then
      deallocate(aux_test_data)
    end if

    end subroutine check_mpi_in_place
