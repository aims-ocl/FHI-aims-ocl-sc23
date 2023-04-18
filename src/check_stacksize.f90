    subroutine check_stacksize()

      use localorb_io, only: localorb_allinfo, localorb_info
      use mpi_tasks, only: myid, aims_stop, mpi_comm_global
      implicit none

      integer:: err, maxstack, curstack, info
      ! We use a minimum stacksize of 60 MB as a proxy for the hard limit, as
      ! Apple in its infinite wisdom does not allow a true unlimited stack size
      ! on Mac OSX.
      integer, parameter:: min_stacksize = 60*1024*1024
      character*200 :: message

      interface
         subroutine get_stacksize(mstack, cstack, error) bind(C,name="get_stacksize_c")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(kind=C_INT) :: mstack, cstack, error
         end subroutine get_stacksize
      end interface

      ! call get_stacksize(maxstack, curstack, err)

      ! if (err .ne. 0) then
      !    write (message,'(A,I0,A,I0)') &
      !       " *** Could not determine stack size for task ", myid, &
      !       "! Error code: ", err
      !    call localorb_allinfo(message)
      ! end if

      ! if (maxstack .eq. -1 .or. maxstack .eq. 0) then
      !    write (message,'(A,I0,A)') "  | Maximum stacksize for task ", myid, &
      !       ": unlimited"
      !    call localorb_allinfo(message)
      ! else
      !    write (message,'(A,I0,A,I0,A)') "  | Maximum stacksize for task ", &
      !       myid, ": ", maxstack/1024/1024, " [Mb]"
      !    call localorb_allinfo(message)
      !    write (message,'(2A)') &
      !       "   *** Stacksize is not set to unlimited. ", &
      !       "This might cause problems..."
      !    call localorb_info(message)
      ! endif

      ! if (curstack .eq. -1 .or. curstack .eq. 0) then
      !    write (message,'(A,I0,A)') "  | Current stacksize for task ", myid, &
      !       ": unlimited"
      !    call localorb_allinfo(message)
      ! else
      !    write (message,'(A,I0,A,I0,A)') "  | Current stacksize for task ", &
      !       myid, ": ", curstack/1024/1024, " [Mb]"
      !    call localorb_allinfo(message)
      !    if (curstack < min_stacksize) then
      !       write (message,'(A)') ""
      !       call localorb_info(message)
      !       write (message,'(2A)') &
      !          " *** We have determined that your current stacksize is not ", &
      !          "unlimited or set to "
      !       call localorb_info(message)
      !       write (message,'(2A)') &
      !          " *** the hard limit, which will cause random crashes in most ",&
      !          "scientific"
      !       call localorb_info(message)
      !       write (message,'(1A)') &
      !          " *** computing packages (including FHI-aims.)"
      !       call localorb_info(message)
      !       write (message,'(1A)') &
      !          " ***"
      !       call localorb_info(message)
      !       write (message,'(2A)') &
      !          " *** For your own safety, we are stopping the calculation ", &
      !          "early, and you will "
      !       call localorb_info(message)
      !       write (message,'(1A)') &
      !          " *** need to increase the current stacksize."
      !       call localorb_info(message)
      !       write (message,'(1A)') &
      !          " ***"
      !       call localorb_info(message)
      !       write (message,'(2A)') &
      !          " *** Please run the 'ulimit -s unlimited' (Linux) or ", &
      !          "'ulimit -s hard' (Mac)"
      !       call localorb_info(message)
      !       write (message,'(1A)') &
      !          " *** command before running FHI-aims."
      !       call localorb_info(message)
      !       write (message,'(1A)') &
      !          " ***"
      !       call localorb_info(message)
      !       write (message,'(2A)') &
      !          " *** We recommend placing this command in your local ", &
      !          ".bashrc and "
      !       call localorb_info(message)
      !       write (message,'(1A)') &
      !          " *** .bash_profile files and in all queuing scripts."
      !       call localorb_info(message)
      !       write (message,'(2A)') &
      !          " ***"
      !       call localorb_info(message)
      !       call mpi_barrier(mpi_comm_global, info)
      !       call aims_stop("Current stacksize too low.  Stopping.")
      !    end if
      ! endif

    end subroutine check_stacksize
