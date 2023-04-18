!****h* FHI-aims/timing_core
!  NAME
!    timing_core
!  SYNOPSIS
module timing_core
!  PURPOSE
!    From the notes in the timing module:
!
!    WPH (2018 Jan 18):  I've split the timing module into two modules:
!    - a timing_core module which handles the time collection and output
!      formatting in an aims-independent fashion
!    - the (original) timing module which functions as a wrapper around
!      timing_core and keeps track of aims-specific timings.
!    This change was made to allow other aims-independent functionality to
!    access aims' core timing subroutines without introducing dependencies on
!    core aims-specific modules such as runtime_choices and dimensions.
!  USES
  implicit none
!  AUTHOR
!    FHI-aims team
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    January 2018:  Forked off of the timing module
!******

  character*120 , private :: info_str

  public :: output_datetime
  public :: time_convert
  public :: get_timestamps
  public :: get_times
  public :: output_timeheader
  public :: output_times
  public :: output_redundant_times
  public :: start_timer
  public :: stop_timer
  public :: output_timer

contains

!-------------------------------------------------------------------------------
!****s* timing_core/output_datetime
!  NAME
!    output_datetime
!  SYNOPSIS
  subroutine output_datetime(time_zero, clock_time_zero)
!  PURPOSE
!    Prints the current date and time to screen, as well as returning them
!    through arguments.
!  USES
    use localorb_io, only: localorb_info, use_unit
!  INPUT
!    none
!  OUTPUT
    real*8, intent(out) :: time_zero
    real*8, intent(out) :: clock_time_zero
!  AUTHOR
!    FHI-aims team.
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  HISTORY
!    January 2018 - Forked off initial_timings
!  SOURCE
    character*8  :: cdate
    character*10 :: ctime

    write(info_str,'(A)') ''
    call localorb_info( info_str, use_unit )

    call date_and_time(cdate, ctime)
    write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, &
         ", Time     :  ", ctime
    call localorb_info( info_str, use_unit )

    call cpu_time(time_zero)
    write(info_str,'(2X,A,E22.15,A)') "Time zero on CPU 1             :  ", &
         time_zero, "  s."
    call localorb_info( info_str, use_unit )

    call time_convert (cdate, ctime, clock_time_zero)
    write(info_str,'(2X,A,F22.3,A)')  "Internal wall clock time zero  :  ", &
         clock_time_zero, "  s."
    call localorb_info( info_str, use_unit )

    write(info_str,*) ''
    call localorb_info ( info_str, use_unit )
  end subroutine output_datetime
!******

!-------------------------------------------------------------------------------
!****s* timing_core/time_convert
!  NAME
!    time_convert
!  SYNOPSIS
  subroutine time_convert (cdate, ctime, seconds)
!  PURPOSE
!    Produces a wall clock time in seconds (real*8), based on the character
!    strings for date and time as returned by the Fortran intrinsic subroutine
!    date_and_time . The time zero is the beginning of 2009. Leap years are
!    accounted for every four years, i.e., NOT correctly from year 2100 on!
!  INPUT
    character*8,  intent(in)  :: cdate
    character*10, intent(in)  :: ctime
!  OUTPUT
    real*8,       intent(out) :: seconds
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
    integer :: value
    integer :: int_year

    real*8 :: year
    real*8 :: month
    real*8 :: day
    real*8 :: hour
    real*8 :: minute
    real*8 :: second
    real*8 :: millisecond

    read ( cdate(1:4), '(I4)' ) value
    int_year = value
    year = dble( value ) - 2009d0  ! year 2009 becomes an arbitrary zero

    day = year * 365 + floor(year/4d0)

    read ( cdate(5:6), '(I2)' ) value

    value = value-1
    do while (value.gt.0)
      if (value.eq.1) then
        day = day+31
      else if (value.eq.2) then
        if (mod(int_year,4).eq.0) then
          day = day+29
        else
          day = day+28
        end if
      else if (value.eq.3) then
        day = day+31
      else if (value.eq.4) then
        day = day+30
      else if (value.eq.5) then
        day = day+31
      else if (value.eq.6) then
        day = day+30
      else if (value.eq.7) then
        day = day+31
      else if (value.eq.8) then
        day = day+31
      else if (value.eq.9) then
        day = day+30
      else if (value.eq.10) then
        day = day+31
      else if (value.eq.11) then
        day = day+30
      end if

      value = value - 1
    end do

    read ( cdate(7:8), '(I2)' ) value
    day = day + dble( value ) - 1

    read ( ctime(1:2), '(I2)' ) value
    hour = dble( value )

    read ( ctime(3:4), '(I2)' ) value
    minute = dble( value )

    read ( ctime(5:6), '(I2)' ) value
    second = dble( value )

    read ( ctime(8:10), '(I3)' ) value
    millisecond = dble( value )

    seconds = day * 24d0 * 3600d0 + hour * 3600d0 + minute*60d0 + second &
         + millisecond * 0.001d0
  end subroutine time_convert
!******

!-------------------------------------------------------------------------------
!****s* timing_core/get_timestamps
!  NAME
!    get_timestamps
!  SYNOPSIS
  subroutine get_timestamps (cpu_seconds, clock_seconds)
!  PURPOSE
!    Returns CPU time stamps and wall clock time stamps
!  INPUT
!    none
!  OUTPUT
      real*8, intent(out) ::cpu_seconds
      real*8, intent(out) ::clock_seconds
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
      ! internal variables
      character*8 cdate
      character*10 ctime

      call cpu_time( cpu_seconds )

      call date_and_time(cdate, ctime)
      call time_convert (cdate, ctime, clock_seconds)
  end subroutine get_timestamps
!******

!-------------------------------------------------------------------------------
!****s* timing_core/get_times
!  NAME
!    get_times
!  SYNOPSIS
  subroutine get_times(cpu_seconds, clock_seconds, tot_cpu_seconds, &
                       tot_clock_seconds, unsynced)
!  PURPOSE
!    This subroutine is meant to be used as a timing tool.  Best to
!    explain it with examples:
!
!    ! - Time a single calculation (named calc):
!    real*8 :: cpu_calc, clock_calc
!    call get_timestamps(cpu_calc, clock_calc)
!    ! Do expensive calculation
!    call get_times(cpu_calc, clock_calc)
!    ! Now cpu_calc contains the maximum CPU time among all nodes and
!    ! clock_calc contains the clock time on myid==0 in seconds.
!
!    ! - Time an (expensive) calculation within a loop (named calc):
!    real*8 :: cpu_calc, clock_calc, tot_cpu_calc, tot_clock_calc
!    tot_cpu_calc = 0.d0; tot_clock_calc = 0.d0
!    do i = 1, N
!      ! Do other stuff
!      call get_timestamps(cpu_calc, clock_calc)
!      ! Do expensive calculation
!      call get_times(cpu_calc, clock_calc, tot_cpu_calc, tot_clock_calc)
!      ! Do other stuff
!    end do
!    ! Now cpu_calc_tot contains the accumulated maximum CPU times
!    ! among all nodes and clock_calc_tot contains the accumulated clock
!    ! time on myid==0 in seconds.
!  USES
    use synchronize_mpi_basic
    implicit none
!  ARGUMENTS
    real*8, intent(INOUT) :: cpu_seconds, clock_seconds
    real*8, intent(INOUT), optional :: tot_cpu_seconds, tot_clock_seconds
    logical, intent(IN), optional :: unsynced
!  INPUTS
!    o cpu_seconds, clock_seconds -- timestamps as given by
!                                    get_timestamps() before the
!                                    calculation
!    o tot_cpu_seconds, tot_clock_seconds -- acculuative total times
!                                            (optional)
!    o unsynced -- (optional) disable syncing of cpu_times
!                  This is useful for fine-grained time accumulations.
!                  But do not forget to sync the accumulated results.
!  OUTPUTS
!    o cpu_seconds, clock_seconds -- CPU and clock times since
!                                    get_timestamps().
!    o tot_cpu_seconds, tot_clock_seconds -- acculuative total times
!                                            (optional)
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2010).
!  SOURCE
    logical :: sync
    real*8 :: rtime, clock_rtime
    character*8 :: cdate
    character*10 :: ctime

    sync = .true.; if (present(unsynced)) sync = .not. unsynced

    ! (max) CPU time
    call cpu_time(rtime)
    cpu_seconds = rtime - cpu_seconds
    if (sync) call sync_timing(cpu_seconds)

    ! wall clock time
    call date_and_time(cdate, ctime)
    call time_convert (cdate, ctime, clock_rtime)
    clock_seconds = clock_rtime - clock_seconds

    ! optionally accumulate
    if (present(tot_cpu_seconds)) then
      tot_cpu_seconds = tot_cpu_seconds + cpu_seconds
    end if
    if (present(tot_clock_seconds)) then
      tot_clock_seconds = tot_clock_seconds + clock_seconds
    end if
  end subroutine get_times
!******

!-------------------------------------------------------------------------------
!****s* timing_core/output_timeheader
!  NAME
!    output_timeheader
!  SYNOPSIS
  subroutine output_timeheader(fmt, description, priority)
!  PURPOSE
!    Prepare header for timing table
!  USES
    use localorb_io, only: localorb_info, use_unit
    implicit none
!  ARGUMENTS
    character(*), intent(IN) :: fmt
    character(*), intent(IN) :: description
    integer, intent(IN), optional :: priority
!  INPUTS
!    o fmt -- Format descriptor for beginning of line
!    o description -- Description of timing table
!    o priority -- How important is it to output?
!  OUTPUTS
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2010).
!  SOURCE
    character*100 :: whole_fmt
    character*200 :: info_str, desc_buffer

    write(whole_fmt, "('(',A,',',A,')')") trim(fmt), "A45,': ',A14,'    ',A16"
    desc_buffer = trim(description)
    write(info_str, whole_fmt) desc_buffer(1:45), "max(cpu_time)", &
         "wall_clock(cpu1)"
    call localorb_info(info_str, use_unit, "(A)", priority)
  end subroutine output_timeheader
!******

!-------------------------------------------------------------------------------
!****s* timing_core/output_times
!  NAME
!    output_times
!  SYNOPSIS
  subroutine output_times(fmt, description, cpu_seconds, clock_seconds, &
                          priority)
!  PURPOSE
!    Output a single line of timing information.
!  USES
    use localorb_io, only: localorb_info, use_unit
    implicit none
!  ARGUMENTS
    character(*), intent(IN) :: fmt
    character(*), intent(IN) :: description
    real*8, intent(IN) :: cpu_seconds, clock_seconds
    integer, intent(IN), optional :: priority
!  INPUTS
!    o fmt -- Format descriptor for beginning of line
!    o description -- Description of what took so long.
!    o cpu_seconds, clock_seconds -- CPU and clock time it took.
!    o priority -- How important is it to output?
!  OUTPUTS
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2010).
!  SOURCE
    character*100 :: whole_fmt
    character*200 :: info_str, desc_buffer

    write(whole_fmt, "('(',A,',',A,')')") trim(fmt), &
         "'| ', A43,': ',F12.3,' s    ',F12.3,' s'"
    desc_buffer = trim(description)
    write(info_str, whole_fmt) desc_buffer(1:43), cpu_seconds, clock_seconds
    call localorb_info(info_str, use_unit, "(A)", priority)
  end subroutine output_times
!******

!-------------------------------------------------------------------------------
!****s* timing_core/output_redundant_times
!  NAME
!    output_times
!  SYNOPSIS
  subroutine output_redundant_times(fmt, description, cpu_seconds, &
                                    clock_seconds, priority)
!  PURPOSE
!    Output a single line of timing information.
!    This is an exact copy of output_times except it adds brackets around
!    the output times to indicate that they are already reflected by another
!    timer written in the same batch of times, e.g., at the end of the code.
!  USES
    use localorb_io, only: localorb_info, use_unit
    implicit none
!  ARGUMENTS
    character(*), intent(IN) :: fmt
    character(*), intent(IN) :: description
    real*8, intent(IN) :: cpu_seconds, clock_seconds
    integer, intent(IN), optional :: priority
!  INPUTS
!    o fmt -- Format descriptor for beginning of line
!    o description -- Description of what took so long.
!    o cpu_seconds, clock_seconds -- CPU and clock time it took.
!    o priority -- How important is it to output?
!  OUTPUTS
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2010).
!  SOURCE
    character*100 :: whole_fmt
    character*200 :: info_str, desc_buffer

    write(whole_fmt, "('(',A,',',A,')')") trim(fmt), &
         "'| ', A43,':(',F12.3,' s)  (',F12.3,' s)'"
    desc_buffer = trim(description)
    write(info_str, whole_fmt) desc_buffer(1:43), cpu_seconds, clock_seconds
    call localorb_info(info_str, use_unit, "(A)", priority)
  end subroutine output_redundant_times
!*****

!-------------------------------------------------------------------------------
!****s* timing_core/start_timer
!  NAME
!    start_timer
!  SYNOPSIS
  subroutine start_timer(timestamps)
!  PURPOSE
!    An alternative interface to get_timestamps().
!  USES
    implicit none
!  ARGUMENTS
    real*8, intent(OUT) :: timestamps(2)   ! In general it will be 4.
!  INPUTS
!    none
!  OUTPUTS
!    timestamps(1) -- Current CPU timestamp
!    timestamps(2) -- Current wallclock timestamp
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2010).
!  SOURCE
    character(*), parameter :: func = 'start_timer'

    call get_timestamps(timestamps(1), timestamps(2))
  end subroutine start_timer
!******

!-------------------------------------------------------------------------------
!****s* timing_core/stop_timer
!  NAME
!    stop_timer
!  SYNOPSIS
  subroutine stop_timer(timestamps, unsynced)
!  PURPOSE
!    An alternative interface to get_timestamps().
!    In order for timestamps(3:4) to be useful, they have to be reset
!    before their first use (e.g. by 'timestamps = 0.d0')
!  USES
    implicit none
!  ARGUMENTS
    real*8, intent(INOUT) :: timestamps(4)
    logical, intent(IN), optional :: unsynced
!  INPUTS
!    timestamps(1) -- CPU timestamp at start (from start_timer())
!    timestamps(2) -- Wallclock timestamp at start (from start_timer())
!    timestamps(3) -- Accumulative CPU time
!    timestamps(4) -- Accumulative Wallclock time
!    unsynced -- (optional) see get_times().
!  OUTPUTS
!    timestamps(1) -- Maximum CPU time of this period
!    timestamps(2) -- Wallclock time of this period
!    timestamps(3) -- Accumulative CPU time
!    timestamps(4) -- Accumulative Wallclock time
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2010).
!  SOURCE
    character(*), parameter :: func = 'stop_timer'
    call get_times(timestamps(1), timestamps(2), timestamps(3), timestamps(4), &
                   unsynced)
  end subroutine stop_timer
!******

!-------------------------------------------------------------------------------
!****s* timing_core/output_timer
!  NAME
!    output_timer
!  SYNOPSIS
  subroutine output_timer(description, timestamps, fmt, priority)
!  PURPOSE
!    An alternative interface to get_timestamps().
!    Can be called either
!       call output_timer('This calculation', timestamps(1:2))
!    for the result of the last timing or
!       call output_timer('All calculations', timestamps(3:4))
!    for the accumulative result of all timings.
!  USES
    implicit none
!  ARGUMENTS
    character(*), intent(IN) :: description
    real*8, intent(IN) :: timestamps(2)
    character(*), intent(IN), optional :: fmt
    integer, intent(IN), optional :: priority
!  INPUTS
!    description -- Description of what took so long
!    timestamps(1) -- CPU time
!    timestamps(2) -- Wallclock time
!    fmt -- Format to be used (default to '2X')
!    priority -- Output priority (default to infinitely high)
!  OUTPUTS
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2010).
!  SOURCE
    character*100 :: myfmt
    character(*), parameter :: func = 'output_timer'

    if (present(fmt)) then
       write(myfmt,"(A)") fmt
    else
       write(myfmt,"(A)") '2X'
    end if
    call output_times(trim(myfmt), description, timestamps(1), timestamps(2), &
          priority)
  end subroutine output_timer
!******
end module timing_core
!******
