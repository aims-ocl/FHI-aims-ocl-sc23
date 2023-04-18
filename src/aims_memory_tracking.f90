!****h* FHI-aims/aims_memory_tracking
!  NAME
!    aims_memory_tracking
!  SYNOPSIS
module aims_memory_tracking
!  PURPOSE
!    Provides a set of subroutines to aid in memory tracking
!  USES
  implicit none
!  AUTHOR
!    The original code for the allocate interface was written by Rainer
!    Johanni for Fock matrix calculations, where it was a collection of module
!    subroutines with the ability to output memory allocated to stdio.  It was
!    later expanded into a general purpose module by William Huhn to include a
!    running tally of total memory allocated through the interface.
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180, 2175 (2009).
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  NOTES
!    The usage of this module is:
!
!    1)  During initialization, call aims_mem_initialize() to set up module
!        variables.
!
!    2)  Allocate important/large arrays using aims_allocate(), an interface
!        around various allocate_* module procedures.
!        - aims_allocate has the following prototype:
!
!          call aims_allocate( variable_name, <array dims>, "+variable_name" )
!
!          where <array dims> is a list of array dimensions and the "+" is
!          omitted if the allocation should be done silently.
!        - Every allocate_* module procedure outputs the size of the allocated
!          array to stdout (if desired), allocates the array, checks the
!          allocation, and finally updates the memory tracking statistics
!        - update_when_allocating() would be called at the end of each
!          allocate_* subroutine to update the memory tracking statistics
!
!    3)  Make sure to deallocate those arrays explicitly using
!        aims_deallocate(), an interface around various deallocate_* module
!        procedures.
!        - aims_deallocate has the following prototype:
!
!          call aims_deallocate( variable_name, "+variable_name" )
!
!          where the "+" is omitted if the deallocation should be done silently.
!        - Every deallocate_* module procedure outputs the size of the
!          deallocated array to stdout (if desired), deallocates the array,
!          checks the deallocation, and finally updates the memory tracking
!          statistics
!        - update_when_deallocating() would be called at the end of each
!          deallocate_* subroutine to update the memory tracking statistics
!
!    4)  At the end of a calculation, call aims_mem_final_output() to output
!        statistics
!
!    *)  In order to measure the memory consumption of a specific part
!        of the code separately, put push_current_memory() before the
!        relevant block of code. After that block,
!        pop_memory_estimate() returns maximum memory usage (in bytes)
!        for that block relative to where push_current_memory() was
!        placed. One could get memory usage of a block within a block
!        by using these procedure in a nested way (unlimited number of
!        levels is supported).
!
!  TODO
!    *   Output more of the largest arrays, not just the largest (maybe five?)
!  HISTORY
!    June 2017 - Fock matrix and SOC versions of allocate combined into one
!                module
!  SOURCE
  private
  public :: aims_mem_initialize
  public :: aims_mem_current_output
  public :: aims_mem_final_output
  public :: aims_allocate
  public :: aims_deallocate
  public :: push_current_memory
  public :: pop_memory_estimate
  ! WPH (14 Feb 2018): update_when_allocating should not be called outside this
  !     module under normal circumstances, since the memory tracking
  !     infrastructure relies on the memory stack being strictly regulated.
  !     However, we've run into a bug where the PGI compiler (17.4 and 17.10)
  !     will rarely corrupt array bounds for select matrices allocated through
  !     the subroutine call, forcing us to inline those aims_allocate
  !     statements.
  public :: update_when_allocating
  public :: aims_mem_debug
  public :: aims_mem_sync

  integer, parameter :: max_name_length = 80
  real*8, parameter :: megabyte = 1048576.d0

  ! When .true., force all allocation/deallocation statements to be written to
  ! screen, as well as the value of current_memory_estimate every time it is
  ! modified
  logical :: aims_mem_debug

  ! When .true., synchronize memory tracking stats after EVERY SCF step.
  ! When .false., only synchronize after the entire FHI-aims run.
  logical :: aims_mem_sync

  ! The current total for tracked allocated memory
  integer*8 :: current_memory_estimate

  ! The maximum value for total allocated memory encountered so far
  integer*8 :: max_memory_estimate

  ! The largest block of tracked allocated memory encountered so far
  integer*8 :: largest_array_size

  ! The name of the largest tracked allocated array encountered so far
  character(max_name_length) :: largest_array_name

  ! The name of the tracked array allocated when the current value for
  ! max_memory_estimate was encountered (not necessarily the largest array!)
  character(max_name_length) :: when_max_memory

  interface aims_allocate
    module procedure allocate_r1
    module procedure allocate_r2
    module procedure allocate_r2_i64
    module procedure allocate_r3
    module procedure allocate_r4
    module procedure allocate_r5
    module procedure allocate_r6
    module procedure allocate_r1_lb
    module procedure allocate_r2_lb
    module procedure allocate_r3_lb
    module procedure allocate_r4_lb
    module procedure allocate_r5_lb

    module procedure allocate_c1
    module procedure allocate_c2
    module procedure allocate_c3
    module procedure allocate_c4
    module procedure allocate_c5
    module procedure allocate_c3_lb

    module procedure allocate_i1
    module procedure allocate_i2
    module procedure allocate_i3
    module procedure allocate_i4
    module procedure allocate_i5
    module procedure allocate_i1_lb
    module procedure allocate_i2_lb
    module procedure allocate_i3_lb
    module procedure allocate_i4_lb
    module procedure allocate_i5_lb

    module procedure allocate_l1
    module procedure allocate_l2
  end interface

  interface aims_deallocate
    module procedure deallocate_r1
    module procedure deallocate_r2
    module procedure deallocate_r3
    module procedure deallocate_r4
    module procedure deallocate_r5
    module procedure deallocate_r6

    module procedure deallocate_c1
    module procedure deallocate_c2
    module procedure deallocate_c3
    module procedure deallocate_c4
    module procedure deallocate_c5

    module procedure deallocate_i1
    module procedure deallocate_i2
    module procedure deallocate_i3
    module procedure deallocate_i4
    module procedure deallocate_i5

    module procedure deallocate_l1
    module procedure deallocate_l2
  end interface

  type :: linked_list
     type(linked_list), pointer :: next
     integer :: depth = 0
     integer*8 :: saved_values(2)
  end type linked_list

  type(linked_list) :: memory_stack

contains

  !===========================================================================!
  !                            Utility subroutines                            !
  !===========================================================================!

  !----------------------------------------------------------------------------
  subroutine aims_mem_initialize()
    ! Initialize the memory tracking infrastructure
    use localorb_io, only: use_unit,OL_norm,localorb_info
    implicit none

    character(256) :: info_str

    aims_mem_debug = .false.
    aims_mem_sync = .true.

    current_memory_estimate = 0
    max_memory_estimate = 0
    largest_array_size = 0
    largest_array_name = "Placeholder Text"
    when_max_memory = "Placeholder Text"

    if (aims_mem_debug) then
      write(info_str,"(A,I15)") &
        "current_memory_estimate ",current_memory_estimate
      call localorb_info(info_str,use_unit,"(4X,A)",OL_norm)
    end if

  end subroutine aims_mem_initialize
  !----------------------------------------------------------------------------
  subroutine aims_mem_current_output()
    ! Output memory tracking statistics in the middle of a calculation
    use localorb_io, only: use_unit,OL_norm,localorb_info
    use mpi_tasks, only: myid,n_tasks,mpi_min_max_all
    use synchronize_mpi_basic, only: sync_real_number,send_string,receive_string
    implicit none

    character(256) :: info_str
    character(1) :: n_digits ! For prettier output
    ! For determining minimum, maximum, and average memory usage
    real*8 :: mem_min
    real*8 :: mem_max
    real*8 :: mem_ave
    integer :: mem_min_loc
    integer :: mem_max_loc

    write(n_digits,"(I1)") int(log(real(n_tasks,8))/log(10d0),4)+1

    write(info_str,"(A)") "Partial memory accounting:"
    call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)

    if (aims_mem_sync) then
      mem_ave = current_memory_estimate/megabyte
      call mpi_min_max_all(mem_ave,mem_min,mem_min_loc,mem_max,mem_max_loc)
      call sync_real_number(mem_ave)
      mem_ave = mem_ave/n_tasks

      call localorb_info("| Current value for overall tracked memory usage:", &
        use_unit,"(2X,A)",OL_norm)

      write(info_str,"(A,F12.3,A,I"//n_digits//",A)") &
        "|   Minimum: ",mem_min," MB (on task ",mem_min_loc,")"
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)

      write(info_str,"(A,F12.3,A,I"//n_digits//",A)") &
        "|   Maximum: ",mem_max," MB (on task ",mem_max_loc,")"
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)

      write(info_str,"(A,F12.3,A)") "|   Average: ",mem_ave," MB"
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)

      mem_ave = max_memory_estimate/megabyte
      call mpi_min_max_all(mem_ave,mem_min,mem_min_loc,mem_max,mem_max_loc)
      call sync_real_number(mem_ave)
      mem_ave = mem_ave/n_tasks

      call localorb_info("| Peak value for overall tracked memory usage:", &
        use_unit,"(2X,A)",OL_norm)

      if (mem_min_loc > 0) then
        if (myid == mem_min_loc) then
          call send_string(when_max_memory,max_name_length,0)
        else if (myid == 0) then
          call receive_string(when_max_memory,max_name_length,mem_min_loc)
        end if
      end if

      write(info_str,"(A,F12.3,A,I"//n_digits//",A)") &
        "|   Minimum: ",mem_min," MB (on task ",mem_min_loc, &
        " after allocating "//trim(when_max_memory)//")"
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)

      if (mem_max_loc > 0) then
        if (myid == mem_max_loc) then
          call send_string(when_max_memory,max_name_length,0)
        else if (myid == 0) then
          call receive_string(when_max_memory,max_name_length,mem_max_loc)
        end if
      end if

      write(info_str,"(A,F12.3,A,I"//n_digits//",A)") &
        "|   Maximum: ",mem_max," MB (on task ",mem_max_loc, &
        " after allocating "//trim(when_max_memory)//")"
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)

      write(info_str,"(A,F12.3,A)") "|   Average: ",mem_ave," MB"
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)

      mem_ave = largest_array_size/megabyte
      call mpi_min_max_all(mem_ave,mem_min,mem_min_loc,mem_max,mem_max_loc)
      call sync_real_number(mem_ave)
      mem_ave = mem_ave/n_tasks

      call localorb_info("| Largest tracked array allocation so far:", &
        use_unit,"(2X,A)",OL_norm)

      if (mem_min_loc > 0) then
        if (myid == mem_min_loc) then
          call send_string(largest_array_name,max_name_length,0)
        else if (myid == 0) then
          call receive_string(largest_array_name,max_name_length,mem_min_loc)
        end if
      end if

      write(info_str,"(A,F12.3,A,I"//n_digits//",A)") &
        "|   Minimum: ",mem_min," MB ("//trim(largest_array_name)// &
        " on task ",mem_min_loc,")"
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)

      if (mem_max_loc > 0) then
        if (myid == mem_max_loc) then
          call send_string(largest_array_name,max_name_length,0)
        else if (myid == 0) then
          call receive_string(largest_array_name,max_name_length,mem_max_loc)
        end if
      end if

      write(info_str,"(A,F12.3,A,I"//n_digits//",A)") &
        "|   Maximum: ",mem_max," MB ("//trim(largest_array_name)// &
        " on task ",mem_max_loc,")"
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)

      write(info_str,"(A,F12.3,A)") "|   Average: ",mem_ave," MB"
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)
    else
      write(info_str,"(A,F12.3,A)") &
        "| Current value for overall tracked memory usage on task 0 :", &
        dble(current_memory_estimate)/megabyte," MB"
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)

      write(info_str,"(A,F12.3,A,A)") &
        "| Peak value for overall tracked memory usage on task 0    :", &
        dble(max_memory_estimate)/megabyte," MB after allocating ", &
        trim(when_max_memory)
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)

      write(info_str,"(A,F12.3,A,A)") &
        "| Largest tracked array allocation on task 0 so far        :", &
        dble(largest_array_size)/megabyte," MB  when allocating ", &
        trim(largest_array_name)
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)
    end if

    if (aims_mem_debug) then
      write(info_str,"(A,I15)") &
        "| current_memory_estimate ",current_memory_estimate
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)
    end if

    write(info_str,"(A)") "Note:  These values currently only include a"// &
      " subset of arrays which are explicitly tracked."
    call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)
    write(info_str,"(A)") 'The "true" memory usage will be greater.'
    call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)

  end subroutine aims_mem_current_output
  !----------------------------------------------------------------------------
  subroutine aims_mem_final_output()
    ! Tear down the memory tracking infrastructure and output important
    ! statistics at the end of a calculation
    use localorb_io, only: use_unit,OL_norm,localorb_info
    use mpi_tasks, only: myid,n_tasks,mpi_min_max_all
    use synchronize_mpi_basic, only: sync_long_integer,sync_real_number, &
        send_string,receive_string
    implicit none

    character(256) :: info_str
    character(1) :: n_digits ! For prettier output
    ! For determining minimum, maximum, and average memory usage
    real*8 :: mem_min
    real*8 :: mem_max
    real*8 :: mem_ave
    integer :: mem_min_loc
    integer :: mem_max_loc

    write(n_digits,"(I1)") int(log(real(n_tasks,8))/log(10d0),4)+1

    ! This needs to be zero across all cpus
    call sync_long_integer(current_memory_estimate)

    write(info_str,"(A)") "Partial memory accounting:"
    call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)

    write(info_str,"(A,F12.6,A)") &
      "| Residual value for overall tracked memory usage across tasks: ", &
      dble(current_memory_estimate)/megabyte," MB (should be 0.000000 MB)"
    call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)

    mem_ave = max_memory_estimate/megabyte
    call mpi_min_max_all(mem_ave,mem_min,mem_min_loc,mem_max,mem_max_loc)
    call sync_real_number(mem_ave)
    mem_ave = mem_ave/n_tasks

    if (abs(current_memory_estimate) < 1.d-15) then
      call localorb_info("| Peak values for overall tracked memory usage:", &
        use_unit,"(10X,A)",OL_norm)

      if (mem_min_loc > 0) then
        if (myid == mem_min_loc) then
          call send_string(when_max_memory,max_name_length,0)
        else if (myid == 0) then
          call receive_string(when_max_memory,max_name_length,mem_min_loc)
        end if
      end if

      write(info_str,"(A,F12.3,A,I"//n_digits//",A)") &
        "|   Minimum: ",mem_min," MB (on task ",mem_min_loc, &
        " after allocating "//trim(when_max_memory)//")"
      call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)

      if (mem_max_loc > 0) then
        if (myid == mem_max_loc) then
          call send_string(when_max_memory,max_name_length,0)
        else if (myid == 0) then
          call receive_string(when_max_memory,max_name_length,mem_max_loc)
        end if
      end if

      write(info_str,"(A,F12.3,A,I"//n_digits//",A)") &
        "|   Maximum: ",mem_max," MB (on task ",mem_max_loc, &
        " after allocating "//trim(when_max_memory)//")"
      call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)

      write(info_str,"(A,F12.3,A)") "|   Average: ",mem_ave," MB"
      call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)
    else
      write(info_str,"(A)") "|"
      call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)
    end if

    mem_ave = largest_array_size/megabyte
    call mpi_min_max_all(mem_ave,mem_min,mem_min_loc,mem_max,mem_max_loc)
    call sync_real_number(mem_ave)
    mem_ave = mem_ave/n_tasks

    call localorb_info("| Largest tracked array allocation:",use_unit, &
      "(10X,A)",OL_norm)

    if (mem_min_loc > 0) then
      if (myid == mem_min_loc) then
        call send_string(largest_array_name,max_name_length,0)
      else if (myid == 0) then
        call receive_string(largest_array_name,max_name_length,mem_min_loc)
      end if
    end if

    write(info_str,"(A,F12.3,A,I"//n_digits//",A)") &
      "|   Minimum: ",mem_min," MB ("//trim(largest_array_name)//" on task ", &
      mem_min_loc,")"
    call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)

    if (mem_max_loc > 0) then
      if (myid == mem_max_loc) then
        call send_string(largest_array_name,max_name_length,0)
      else if (myid == 0) then
        call receive_string(largest_array_name,max_name_length,mem_max_loc)
      end if
    end if

    write(info_str,"(A,F12.3,A,I"//n_digits//",A)") &
      "|   Maximum: ",mem_max," MB ("//trim(largest_array_name)//" on task ", &
      mem_max_loc,")"
    call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)

    write(info_str,"(A,F12.3,A)") "|   Average: ",mem_ave," MB"
    call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)

    if (abs(current_memory_estimate) > 1.d-15) then
      write(info_str,"(A)") "|"
      call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)
      write(info_str,"(A)") "| Because the residual tracked memory usage"// &
        " was not zero, the peak value for overall"
      call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)
      write(info_str,"(A)") "| tracked memory usage was not output, as it"// &
        " is probably incorrect.  This is likely due "
      call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)
      write(info_str,"(A)") "| to a missing aims_deallocate statement in"// &
        " the source code, which should be otherwise"
      call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)
      write(info_str,"(A)") "| harmless."
      call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)
      write(info_str,"(A)") "|"
      call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)
    end if

    if (aims_mem_debug) then
      write(info_str,"(A,I15)") &
        "| current_memory_estimate ",current_memory_estimate
      call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)
    end if

    write(info_str,"(A)") "Note:  These values currently only include a"// &
      " subset of arrays which are explicitly tracked."
    call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)
    write(info_str,"(A)") 'The "true" memory usage will be greater.'
    call localorb_info(info_str,use_unit,"(10X,A)",OL_norm)

    current_memory_estimate = 0
    max_memory_estimate = 0
    largest_array_size = 0
    largest_array_name = "Placeholder Text"
    when_max_memory = "Placeholder Text"

  end subroutine aims_mem_final_output
  !----------------------------------------------------------------------------
  subroutine update_when_allocating(info,name,mem)
    ! Error check the allocation and update memory statistics
    ! Should be called at the end of every allocate_* subroutine
    use localorb_io, only: use_unit,OL_norm,localorb_info
    use mpi_tasks, only: check_allocation
    implicit none

    integer, intent(in) :: info
    character(*), intent(in) :: name
    integer*8, intent(in) :: mem
    character(256) :: info_str

    call check_allocation(info,name)
    current_memory_estimate = current_memory_estimate+mem

    if (current_memory_estimate > max_memory_estimate) then
      max_memory_estimate = current_memory_estimate

      if (name(1:1) == "+") then
        when_max_memory = name(2:)
      else
        when_max_memory = name
      end if
    end if

    if (mem > largest_array_size) then
      largest_array_size = mem

      if (name(1:1) == "+") then
        largest_array_name = name(2:)
      else
        largest_array_name = name
      end if
    end if

    if (aims_mem_debug) then
      write(info_str,"(A,I15,A,I15)") &
        "mem ",mem," current_memory_estimate ",current_memory_estimate
      call localorb_info(info_str,use_unit,"(4X,A)",OL_norm)
    end if

  end subroutine update_when_allocating
  !----------------------------------------------------------------------------
  subroutine update_when_deallocating(mem)
    ! Error check the deallocation and update memory statistics
    ! Should be called at the end of every deallocate_* subroutine
    use localorb_io, only: use_unit,OL_norm,localorb_info
    implicit none

    integer*8, intent(in) :: mem
    character(256) :: info_str

    current_memory_estimate = current_memory_estimate-mem

    if (aims_mem_debug) then
      write(info_str,"(A,I15,A,I15)") &
        "mem ",mem," current_memory_estimate ",current_memory_estimate
      call localorb_info(info_str,use_unit,"(4X,A)",OL_norm)
    end if

  end subroutine update_when_deallocating
  !----------------------------------------------------------------------------
  subroutine print_before_allocating(mem,name)
    ! Print information about an allocation
    use localorb_io, only: use_unit,OL_norm,localorb_info
    implicit none

    integer*8, intent(in) :: mem
    character(*), intent(in) :: name
    integer :: name1
    character(256) :: info_str

    if (name(1:1) == "+" .or. aims_mem_debug) then
      if (name(1:1) == "+") then
        name1 = 2
      else
        name1 = 1
      end if
      write(info_str,"(A,F12.3,A,A)") &
        "Allocating ",dble(mem)/megabyte," MB for ",name(name1:)
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)
    end if

  end subroutine
  !----------------------------------------------------------------------------
  subroutine print_before_deallocating(mem,name)
    ! Print information about a deallocation
    use localorb_io, only: use_unit,OL_norm,localorb_info
    implicit none

    integer*8, intent(in) :: mem
    character(*), intent(in) :: name
    integer :: name1
    character(256) :: info_str

    if (name(1:1) == "+" .or. aims_mem_debug) then
      if (name(1:1) == "+") then
        name1 = 2
      else
        name1 = 1
      end if
      write(info_str,"(A,F12.3,A,A)") &
        "Deallocating ",dble(mem)/megabyte," MB for ",name(name1:)
      call localorb_info(info_str,use_unit,"(2X,A)",OL_norm)
    end if

  end subroutine
  !----------------------------------------------------------------------------
  recursive subroutine linked_list_push(node,depth)
    type(linked_list), intent(inout) :: node
    integer, intent(inout) :: depth
    integer, save :: level = 0

    if (level == depth) then
      allocate(node%next)
      depth = depth+1
      node%saved_values = [current_memory_estimate,max_memory_estimate]
      max_memory_estimate = 0
    else
      level = level+1
      call linked_list_push(node%next,depth)
      level = level-1
    end if
  end subroutine linked_list_push
  !----------------------------------------------------------------------------
  recursive integer*8 function linked_list_pop(node,depth) result(y)
    type(linked_list), intent(inout) :: node
    integer, intent(inout) :: depth
    integer, save :: level = 1

    if (level == depth) then
      y = max_memory_estimate-node%saved_values(1)
      max_memory_estimate = max(max_memory_estimate,node%saved_values(2))
      depth = depth-1
      deallocate(node%next)
    else
      level = level+1
      y = linked_list_pop(node%next,depth)
      level = level-1
    end if
  end function linked_list_pop
  !----------------------------------------------------------------------------
  subroutine push_current_memory()
    call linked_list_push(memory_stack,memory_stack%depth)
  end subroutine push_current_memory
  !----------------------------------------------------------------------------
  integer*8 function pop_memory_estimate() result(y)
    y = linked_list_pop(memory_stack,memory_stack%depth)
  end function pop_memory_estimate
  !----------------------------------------------------------------------------

  !===========================================================================!
  !                         Allocation subroutines                            !
  !===========================================================================!

  !----------------------------------------------------------------------------
  subroutine allocate_r1(arr,ub1,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: ub1
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(8,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_r2(arr,ub1,ub2,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:)
    integer, intent(in) :: ub1,ub2
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(8,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_r2_i64(arr,ub1,ub2,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:)
    integer*8, intent(in) :: ub1,ub2
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = ub1*ub2*int(8,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_r3(arr,ub1,ub2,ub3,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:,:)
    integer, intent(in) :: ub1,ub2,ub3
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(ub3,kind=8)*int(8,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2,ub3),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_r4(arr,ub1,ub2,ub3,ub4,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:,:,:)
    integer, intent(in) :: ub1,ub2,ub3,ub4
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(ub3,kind=8)*int(ub4,kind=8) &
      *int(8,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2,ub3,ub4),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_r5(arr,ub1,ub2,ub3,ub4,ub5,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:,:,:,:)
    integer, intent(in) :: ub1,ub2,ub3,ub4,ub5
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(ub3,kind=8)*int(ub4,kind=8) &
      *int(ub5,kind=8)*int(8,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2,ub3,ub4,ub5),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_r6(arr,ub1,ub2,ub3,ub4,ub5,ub6,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:,:,:,:,:)
    integer, intent(in) :: ub1,ub2,ub3,ub4,ub5,ub6
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(ub3,kind=8)*int(ub4,kind=8) &
      *int(ub5,kind=8)*int(ub6,kind=8)*int(8,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2,ub3,ub4,ub5,ub6),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_r1_lb(arr,lb1,ub1,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: lb1,ub1
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1-lb1+1,kind=8)*int(8,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(lb1:ub1),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_r2_lb(arr,lb1,ub1,lb2,ub2,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:)
    integer, intent(in) :: lb1,ub1,lb2,ub2
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1-lb1+1,kind=8)*int(ub2-lb2+1,kind=8)*int(8,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(lb1:ub1,lb2:ub2),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_r3_lb(arr,lb1,ub1,lb2,ub2,lb3,ub3,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:,:)
    integer, intent(in) :: lb1,ub1,lb2,ub2,lb3,ub3
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1-lb1+1,kind=8)*int(ub2-lb2+1,kind=8)*int(ub3-lb3+1,kind=8) &
      *int(8,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(lb1:ub1,lb2:ub2,lb3:ub3),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_r4_lb(arr,lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:,:,:)
    integer, intent(in) :: lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1-lb1+1,kind=8)*int(ub2-lb2+1,kind=8)*int(ub3-lb3+1,kind=8) &
      *int(ub4-lb4+1,kind=8)*int(8,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_r5_lb(arr,lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4,lb5,ub5,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:,:,:,:)
    integer, intent(in) :: lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4,lb5,ub5
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1-lb1+1,kind=8)*int(ub2-lb2+1,kind=8)*int(ub3-lb3+1,kind=8) &
      *int(ub4-lb4+1,kind=8)*int(ub5-lb5+1,kind=8)*int(8,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb5:ub5),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_c1(arr,ub1,name)
    implicit none

    complex*16, allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: ub1
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(16,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_c2(arr,ub1,ub2,name)
    implicit none

    complex*16, allocatable, intent(inout) :: arr(:,:)
    integer, intent(in) :: ub1,ub2
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(16,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_c3(arr,ub1,ub2,ub3,name)
    implicit none

    complex*16, allocatable, intent(inout) :: arr(:,:,:)
    integer, intent(in) :: ub1,ub2,ub3
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(ub3,kind=8)*int(16,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2,ub3),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_c4(arr,ub1,ub2,ub3,ub4,name)
    implicit none

    complex*16, allocatable, intent(in out) :: arr(:,:,:,:)
    integer, intent(in) :: ub1,ub2,ub3,ub4
    character(*), intent(in) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(ub3,kind=8)*int(ub4,kind=8) &
      *int(16,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2,ub3,ub4),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine allocate_c4
  !----------------------------------------------------------------------------
  subroutine allocate_c5(arr,ub1,ub2,ub3,ub4,ub5,name)
    implicit none

    complex*16, allocatable, intent(in out) :: arr(:,:,:,:,:)
    integer, intent(in) :: ub1,ub2,ub3,ub4,ub5
    character(*), intent(in) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(ub3,kind=8)*int(ub4,kind=8) &
      *int(ub5,kind=8)*int(16,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2,ub3,ub4,ub5),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine allocate_c5
  !----------------------------------------------------------------------------
  subroutine allocate_c3_lb(arr,lb1,ub1,lb2,ub2,lb3,ub3,name)
    implicit none

    complex*16, allocatable, intent(inout) :: arr(:,:,:)
    integer, intent(in) :: lb1,ub1,lb2,ub2,lb3,ub3
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1-lb1+1,kind=8)*int(ub2-lb2+1,kind=8)*int(ub3-lb3+1,kind=8) &
      *int(16,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(lb1:ub1,lb2:ub2,lb3:ub3),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_i1(arr,ub1,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: ub1
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(4,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_i2(arr,ub1,ub2,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:,:)
    integer, intent(in) :: ub1,ub2
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(4,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_i3(arr,ub1,ub2,ub3,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:,:,:)
    integer, intent(in) :: ub1,ub2,ub3
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(ub3,kind=8)*int(4,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2,ub3),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_i4(arr,ub1,ub2,ub3,ub4,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:,:,:,:)
    integer, intent(in) :: ub1,ub2,ub3,ub4
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(ub3,kind=8)*int(ub4,kind=8) &
      *int(4,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2,ub3,ub4),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_i5(arr,ub1,ub2,ub3,ub4,ub5,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:,:,:,:,:)
    integer, intent(in) :: ub1,ub2,ub3,ub4,ub5
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(ub3,kind=8)*int(ub4,kind=8) &
      *int(ub5,kind=8)*int(4,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2,ub3,ub4,ub5),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_i1_lb(arr,lb1,ub1,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: lb1,ub1
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1-lb1+1,kind=8)*int(4,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(lb1:ub1),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_i2_lb(arr,lb1,ub1,lb2,ub2,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:,:)
    integer, intent(in) :: lb1,ub1,lb2,ub2
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1-lb1+1,kind=8)*int(ub2-lb2+1,kind=8)*int(4,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(lb1:ub1,lb2:ub2),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_i3_lb(arr,lb1,ub1,lb2,ub2,lb3,ub3,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:,:,:)
    integer, intent(in) :: lb1,ub1,lb2,ub2,lb3,ub3
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1-lb1+1,kind=8)*int(ub2-lb2+1,kind=8)*int(ub3-lb3+1,kind=8) &
      *int(4,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(lb1:ub1,lb2:ub2,lb3:ub3),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_i4_lb(arr,lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:,:,:,:)
    integer, intent(in) :: lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1-lb1+1,kind=8)*int(ub2-lb2+1,kind=8)*int(ub3-lb3+1,kind=8) &
      *int(ub4-lb4+1,kind=8)*int(4,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_i5_lb(arr,lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4,lb5,ub5,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:,:,:,:,:)
    integer, intent(in) :: lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4,lb5,ub5
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1-lb1+1,kind=8)*int(ub2-lb2+1,kind=8)*int(ub3-lb3+1,kind=8) &
      *int(ub4-lb4+1,kind=8)*int(ub5-lb5+1,kind=8)*int(4,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb5:ub5),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_l1(arr,ub1,name)
    implicit none

    logical, allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: ub1
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(4,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine allocate_l2(arr,ub1,ub2,name)
    implicit none

    logical, allocatable, intent(inout) :: arr(:,:)
    integer, intent(in) :: ub1,ub2
    character(*) :: name
    integer*8 :: mem
    integer :: info

    mem = int(ub1,kind=8)*int(ub2,kind=8)*int(4,kind=8)

    call print_before_allocating(mem,name)
    allocate(arr(ub1,ub2),stat=info)
    call update_when_allocating(info,name,mem)

  end subroutine
  !----------------------------------------------------------------------------

  !===========================================================================!
  !                        Deallocation subroutines                           !
  !===========================================================================!

  !----------------------------------------------------------------------------
  subroutine deallocate_r1(arr,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:)
    character(*) :: name
    integer*8 :: mem

    mem = size(arr,kind=8)*int(8,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_r2(arr,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:)
    character(*) :: name
    integer*8 :: mem

    mem = size(arr,kind=8)*int(8,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_r3(arr,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:,:)
    character(*) :: name
    integer*8 :: mem

    mem = size(arr,kind=8)*int(8,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_r4(arr,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:,:,:)
    character(*) :: name
    integer*8 :: mem

    mem = size(arr,kind=8)*int(8,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_r5(arr,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:,:,:,:)
    character(*) :: name
    integer*8 :: mem

    mem = size(arr,kind=8)*int(8,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_r6(arr,name)
    implicit none

    real*8, allocatable, intent(inout) :: arr(:,:,:,:,:,:)
    character(*) :: name
    integer*8 :: mem

    mem = size(arr,kind=8)*int(8,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_c1(arr,name)
    implicit none

    complex*16, allocatable, intent(inout) :: arr(:)
    character(*) :: name
    integer*8 :: mem

    mem = int(16,kind=8)*size(arr,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_c2(arr,name)
    implicit none

    complex*16, allocatable, intent(inout) :: arr(:,:)
    character(*) :: name
    integer*8 :: mem

    mem = int(16,kind=8)*size(arr,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_c3(arr,name)
    implicit none

    complex*16, allocatable, intent(inout) :: arr(:,:,:)
    character(*) :: name
    integer*8 :: mem

    mem = int(16,kind=8)*size(arr,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_c4(arr,name)
    implicit none

    complex*16, allocatable, intent(inout) :: arr(:,:,:,:)
    character(*) :: name
    integer*8 :: mem

    mem = int(16,kind=8)*size(arr,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_c5(arr,name)
    implicit none

    complex*16, allocatable, intent(inout) :: arr(:,:,:,:,:)
    character(*) :: name
    integer*8 :: mem

    mem = int(16,kind=8)*size(arr,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_i1(arr,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:)
    character(*) :: name
    integer*8 :: mem

    mem = size(arr,kind=8)*int(4,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_i2(arr,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:,:)
    character(*) :: name
    integer*8 :: mem

    mem = size(arr,kind=8)*int(4,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_i3(arr,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:,:,:)
    character(*) :: name
    integer*8 :: mem

    mem = size(arr,kind=8)*int(4,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_i4(arr,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:,:,:,:)
    character(*) :: name
    integer*8 :: mem

    mem = size(arr,kind=8)*int(4,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_i5(arr,name)
    implicit none

    integer, allocatable, intent(inout) :: arr(:,:,:,:,:)
    character(*) :: name
    integer*8 :: mem

    mem = size(arr,kind=8)*int(4,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_l1(arr,name)
    implicit none

    logical, allocatable, intent(inout) :: arr(:)
    character(*) :: name
    integer*8 :: mem

    mem = size(arr,kind=8)*int(4,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
  subroutine deallocate_l2(arr,name)
    implicit none

    logical, allocatable, intent(inout) :: arr(:,:)
    character(*) :: name
    integer*8 :: mem

    mem = size(arr,kind=8)*int(4,kind=8)

    call print_before_deallocating(mem,name)
    deallocate(arr)
    call update_when_deallocating(mem)

  end subroutine
  !----------------------------------------------------------------------------
end module aims_memory_tracking
