!****h* FHI-aims/debug_output
!  NAME
!    debug_output
!  SYNOPSIS

module debug_output

  !  PURPOSE
  !    Subroutines for debug output
  !  USES

  implicit none

  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

contains

  !----------------------------------------------------------------------------
  !****s* debug_output/debug_array
  !  NAME
  !    debug_array
  !  SYNOPSIS

  subroutine debug_array(unit, array, pretext, fmt, donode, dopar, tofile)

    !  PURPOSE
    !
    !  USES

    use mpi_tasks, only: myid
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: unit
    real*8, intent(IN) :: array(:,:)
    character(*), intent(IN), optional :: pretext
    character(*), intent(IN), optional :: fmt
    logical, intent(IN), optional :: donode
    logical, intent(IN), optional :: dopar
    character(*), intent(IN), optional :: tofile

    !  INPUTS
    !    o unit -- unit number to output to
    !    o array -- array to output
    !    o pretext -- text to start all lines (def: None)
    !    o fmt -- format for a single number (def: '(ES24.16)')
    !    o donode -- output node number (myid) in front of each line (def: F)
    !    o dopar -- add parentheses (def: F)
    !    o tofile -- if present, open a temporary file with given unit
    !                to write to.
    !  OUTPUTS
    !    none (writes to disk)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character*60 :: usefmt
    logical :: use_dopar
    integer :: i, j

    if (present(fmt)) then
       if (fmt(1:1) == '(') then
          usefmt = trim(fmt)
       else
          write(usefmt, "('(',A,')')") trim(fmt)
       end if
    else
       usefmt = '(ES24.16)'
    end if
    use_dopar = .false.; if (present(dopar)) use_dopar = dopar
    if (present(tofile)) open(unit, FILE=tofile, status='replace')
    do i = 1, size(array, 1)
       if (present(donode)) then
          if (donode) write(unit, "(I2,' ')", advance='NO') myid
       end if
       if (present(pretext)) write(unit, "(A,' ')", advance='NO') pretext
       if (use_dopar) write(unit, "('(')", advance='NO')
       do j = 1, size(array, 2)
          write(unit, usefmt, advance='NO') array(i, j)
       end do
       if (use_dopar) write(unit,"(')')", advance='NO')
       write(unit, "()")
    end do
    write(unit,"()")
    if (present(tofile)) close(unit)

  end subroutine debug_array
  !******
  !----------------------------------------------------------------------------
  !****s* debug_output/debug_plot_data
  !  NAME
  !    debug_plot_data
  !  SYNOPSIS

  subroutine debug_plot_data(filename, x, y1, y2, y3, y4)

    !  PURPOSE
    !    Output into xmgrace/gnuplot format
    !  USES

    implicit none

    !  ARGUMENTS

    character(*), intent(IN)            :: filename
    real*8, intent(IN) :: x(:), y1(:)
    real*8, intent(IN), optional :: y2(:), y3(:), y4(:)

    !  INPUTS
    !    o filename -- file to write to (gets replaced)
    !  OUTPUTS
    !    o x -- x coordinate (first column)
    !    o y1 -- y coordintate (second column)
    !    o y2, y3, y4 -- optional additional columns
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i

    open(20, FILE=filename, status='replace')
    do i = 1, size(x)
       write(20,"(2E25.16E3)", advance='NO') x(i), y1(i)
       if (present(y2)) write(20,"(E25.16E3)", advance='NO') y2(i)
       if (present(y3)) write(20,"(E25.16E3)", advance='NO') y3(i)
       if (present(y4)) write(20,"(E25.16E3)", advance='NO') y4(i)
       write(20,"()")
    end do
    close(20)

  end subroutine debug_plot_data
  !******
  !--------------------------------------------------------------------------
  !****s* debug_output/debug_plot_log_data
  !  NAME
  !    debug_plot_log_data
  !  SYNOPSIS

  subroutine debug_plot_log_data(filename, ln0, lnrange, y1, y2, y3, y4)

    !  PURPOSE
    !    plot data on log grid
    !  USES

    implicit none

    !  ARGUMENTS

    character(*), intent(IN) :: filename
    real*8, intent(IN) :: ln0, lnrange
    real*8, intent(IN) :: y1(:)
    real*8, intent(IN), optional :: y2(:), y3(:), y4(:)

    !  INPUTS
    !    o filename -- file to write to (gets replaced)
    !  OUTPUTS
    !    o x -- x coordinate (first column)
    !    o y1 -- y coordintate (second column)
    !    o y2, y3, y4 -- optional additional columns
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i
    real*8, allocatable :: x(:)

    allocate(x(size(y1)))
    do i = 1, size(y1)
       x(i) = exp(ln0 + (i-1) * lnrange / size(y1))
    end do

    call debug_plot_data(filename, x, y1, y2, y3, y4)

  end subroutine debug_plot_log_data
  !******
end module debug_output
!******
