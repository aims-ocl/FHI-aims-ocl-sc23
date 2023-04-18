!****h* FHI-aims/check_cpu_consistency
!  NAME
!    check_cpu_consistency_matrix
!  SYNOPSIS
!
 subroutine check_cpu_consistency_matrix (matrix, matrix_name, dim1, dim2)
!  PURPOSE
!    Subroutine checks that all CPUs work on the same lattice and
!    stress tensors. It appears that for certain versions of the
!    Intel Fortran compiler a single bit 
!    may be flipped between stored values of the system geometry on
!    different MPI tasks.
! USES
   use dimensions
   use localorb_io
   use mpi_utilities
   use runtime_choices, only : output_level, cpu_consistency_threshold
   use timing, only : warn_cpu_consistency
   use synchronize_mpi_basic, only: sync_vector,sync_logical,get_max_double
   use mpi_tasks, only: SYNC_OR

   implicit none

!  ARGUMENTS

   integer, intent(in)                        :: dim1, dim2
   real*8, intent (inout), dimension(dim1,dim2) :: matrix
   character (len=*) , intent (in)              :: matrix_name

!  INPUT
!    matrix: Array of two dimensions, to be checked for consistency between
!            different MPI tasks. Example (coming from outside) would be 
!            coords(3, n_atoms).
!            Note that the requirement of 2 dimensions is not strict - just use
!            a "1" for dim2 of a one-dimensional array and the array will be
!            reshaped internally as is the standard in Fortran.
!    matrix_name: Array name of incoming array, for writeout below
!    dim1, dim2: array dimensions of the incoming array. 
!            Trying this subroutine on any integer*8 array dimension will lead to 
!            problems, but should be addressed easily by creating a modified version.
!  OUTPUT

!  AUTHOR
!    Bjoern Lange 2013 / Volker Blum 2016
!  HISTORY
!    Release version, FHI-aims (2016).
!  COPYRIGHT
!   Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  SOURCE

! local variables
    
   real*8, dimension (:,:), allocatable         :: processor_zero_matrix
   real*8                                       :: delta
   integer                                      :: i_dim1, i_dim2
   logical                                      :: found_inconsistency
   logical                                      :: stop_run
   real*8                                       :: max_delta_task, max_delta_all
    
   character*200 :: info_str
   
! begin work

   ! Check that dimensions are greater than zero
   if (dim1 .eq. 0) then
      write(info_str,'(2X,2A)') "*** Error: First dimension in ", &
         "check_cpu_consistency.f90 is ZERO!" 
      call localorb_info(info_str)
      stop
   endif

   if (dim2 .eq. 0) then
      write(info_str,'(2X,2A)') "*** Error: Second dimension in ", &
         "check_cpu_consistency.f90 is ZERO!" 
      call localorb_info(info_str)
      stop
   endif

   if(.not.(USE_MPI)) then
      return
   endif
   
   allocate (processor_zero_matrix(1:dim1,1:dim2))

   ! Set collector array to zero 
   found_inconsistency = .false.
   stop_run = .false.
   max_delta_task = 0.d0
   max_delta_all = 0.d0
   processor_zero_matrix = 0.d0


   ! Feed values from task_id 0 in collector array
   if (myid == 0) then
     processor_zero_matrix(1:dim1,1:dim2) = matrix(1:dim1,1:dim2) 
   end if

   ! Synchronize the collector arrays on all tasks
   call sync_vector(processor_zero_matrix, dim1 * dim2, &
                    mpi_comm_global)

   ! Check for inconsistencies
   do i_dim1 = 1, dim1, 1 
     do i_dim2 = 1, dim2, 1
       delta = matrix(i_dim1,i_dim2) - processor_zero_matrix(i_dim1,i_dim2)
       if ( delta .ne. 0.d0 ) then
         found_inconsistency = .true.
         if (output_level.eq.'full') then ! per-CPU output can be a lot of output - write only for 'full'
           write(use_unit,'(1X,A,I4,2X,I1,2X,I1)') &
           "* WARNING: Array mismatch between different CPUs on myID ", myid
           write(use_unit,'(1X,A,E33.16)') &
           "* Value on task 0 is: ", processor_zero_matrix(i_dim1,i_dim2)
           write(use_unit,'(1X,A,I10,A,E33.16)') &
           "* Value on MPI task number ", myid, " is: ", matrix(i_dim1,i_dim2)
           write(use_unit,'(1X,A,E16.7)') & 
           "* The detected difference DELTA is: ", delta
         end if
         if (dabs(delta).gt.max_delta_task) then
             max_delta_task = dabs(delta)
         end if
         if (dabs(delta).gt.cpu_consistency_threshold) then
            stop_run = .true.
         end if 
       endif

     enddo
   enddo

   call sync_logical(found_inconsistency,SYNC_OR)
   call sync_logical(stop_run,SYNC_OR)
   call get_max_double(max_delta_all, max_delta_task)

   ! Did we find an inconsistency ? Tell the user!
   if (found_inconsistency) then

      write(info_str,'(1X,3A)') &
         "* WARNING: The array ", &
         trim(matrix_name), &
         " had different values on different MPI tasks. This is unusual."
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(1X,A)') &
         "* The array will be set to the values stored on MPI task 0 on all other processors."
      call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '* We have observed that flips of individual bits can happen for'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '* arrays that should formally be the same, but this happens only for'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '* certain compilers and / or hardware platforms.'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A,E13.6,A)') &
         '* The numerical value of the maximal discrepancy in your case is: ',max_delta_all,' .'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A,E13.6,A)') &
         '* For comparison, the predefined cpu_consistency_threshold is   : ',cpu_consistency_threshold,' .'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         ' '
       call localorb_info(info_str,use_unit,'(A)')

       if (stop_run) then

         write (info_str,'(1X,A)') &
           '* The discrepancy was larger than the predefined threshold'
         call localorb_info(info_str,use_unit,'(A)')
         write (info_str,'(1X,A)') &
           '* cpu_consistency_threshold . In this case, we are not sure'
         call localorb_info(info_str,use_unit,'(A)')
         write (info_str,'(1X,A)') &
           '* that there is not some unwanted reason for the discrepancy'
         call localorb_info(info_str,use_unit,'(A)')
         write (info_str,'(1X,A)') &
           '* and stop the run. Apologies for any inconvenience.'
         call localorb_info(info_str,use_unit,'(A)')
         write (info_str,'(1X,A)') &
           "* You may obtain more numerical information by rerunning with 'output_level full'."
         call localorb_info(info_str,use_unit,'(A)')
         write (info_str,'(1X,A)') &
           '* Note that one can control this behavior using the'
         call localorb_info(info_str,use_unit,'(A)')
         write (info_str,'(1X,A)') &
           '* check_cpu_consistency and cpu_consistency_threshold keywords.'
         call localorb_info(info_str,use_unit,'(A)')
         write (info_str,'(1X,A)') &
           '* This warning is here to help make sure that you are aware of the problem'
         call localorb_info(info_str,use_unit,'(A)')
         write (info_str,'(1X,A)') &
           '* and can check if overriding would be safe.'
         call localorb_info(info_str,use_unit,'(A)')

         ! aims_stop_coll will work here because we synchronized 
         ! the stop_run array further up.
         call aims_stop_coll('Consistency of arrays between different MPI tasks not ensured.', & 
              'check_cpu_consistency_matrix')

       else
         ! The discepancy was smaller than the predefined threshold
         ! cpu_consistency_threshold . In this case, we will continue
         ! the run with the values stored on processor myid=0 but
         ! issue a warning at the end of the run.

         write (info_str,'(1X,A)') &
           '* Usually, this is not problematic but it would be good to check '
         call localorb_info(info_str,use_unit,'(A)')
         write (info_str,'(1X,A)') &
           '* that this is indeed all that is going on.'
         call localorb_info(info_str,use_unit,'(A)')

         ! Also output a warning at the end of the run
         warn_cpu_consistency = .true.

         ! Set data identical to process zero
         matrix(1:dim1,1:dim2) = processor_zero_matrix(1:dim1,1:dim2)

       end if

   endif
      
   ! Clean up
   deallocate (processor_zero_matrix)

 end subroutine check_cpu_consistency_matrix

