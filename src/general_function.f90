!****h* FHI-aims/general_function
!  NAME
!    general_function
!  SYNOPSIS

module general_function

  !  PURPOSE
  !
  !     Run-time definable general function in one variable.
  !    
  !  USES
  use constants
  use mpi_tasks
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

  integer, parameter :: GF_ZERO = 1
  integer, parameter :: GF_POWER = 2
  integer, parameter :: GF_SIN = 3
  integer, parameter :: GF_COS = 4

  real*8, parameter, private :: c_light_cm_ps = 2.9979246d-2


  type gen_func
     private
     integer :: func_type
     real*8 :: param
  end type gen_func
  
contains

  !****s* general_function/parse_general_function
  !  NAME
  !    parse_general_function
  !  SYNOPSIS
  subroutine parse_general_function(func_str, func_struct)

    !  PURPOSE
    !    Parse a general function described in func (see below)
    !    at given times.
    !  USES

    implicit none

    !  ARGUMENTS

    character*(*), intent(IN) :: func_str
    type(gen_func), intent(OUT) :: func_struct

    !  INPUTS
    !    o func_str -- string containing the actual function, which can be
    !        one of
    !        the keywords 'none' (a dummy), 'constant', 'linear', 'quadratic',
    !        'cubic', or 'polynomial N' where N is the order, or
    !        'sin/cos freq unit', with unit being either 'cm^-1' (wave number)
    !        or 'fs' (period length).
    !  OUTPUT
    !    o func_struct -- Internal representation of that function
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character*200 :: desc_str, unit, info_str
    real*8 :: freq
    character(*), parameter :: func = 'parse_general_function'

    read(func_str, *) desc_str
    select case (desc_str)
    case('none')
       func_struct%func_type = GF_ZERO
       func_struct%param = 0.d0
    case('constant')
       func_struct%func_type = GF_POWER
       func_struct%param = 0.d0
    case('linear')
       func_struct%func_type = GF_POWER
       func_struct%param = 1.d0
    case('quadratic')
       func_struct%func_type = GF_POWER
       func_struct%param = 2.d0
    case('cubic')
       func_struct%func_type = GF_POWER
       func_struct%param = 3.d0
    case('polynomial', 'power')
       func_struct%func_type = GF_POWER
       read(func_str, *) desc_str, func_struct%param
    case('sin', 'cos')
       if (desc_str == 'sin') then
          func_struct%func_type = GF_SIN
       else
          func_struct%func_type = GF_COS
       end if
       read(func_str, *) desc_str, freq, unit
       select case (unit)
       case('cm^-1')   ! wave number -> inverse wave length
          func_struct%param = freq * 2*pi * c_light_cm_ps
       case('fs')      ! period
          func_struct%param = 1000 * 2*pi / freq
       case default
          write(info_str, "(A,': Unknown unit ',A)") trim(func), trim(unit)
          call aims_stop(info_str)
       end select
    case default
       write(info_str, "(A,': Unknown function type ',A)") &
       & trim(func), trim(desc_str)
       call aims_stop(info_str)
    end select

  end subroutine parse_general_function
  !******
  !----------------------------------------------------------------------------
  !****s* general_function/dealloc_general_function
  !  NAME
  !    dealloc_general_function
  !  SYNOPSIS
  subroutine dealloc_general_function(func_struct)

    !  PURPOSE
    !    Dealloc any dynamic data structures within func_struct
    !  USES

    implicit none

    !  ARGUMENTS

    type(gen_func), intent(INOUT) :: func_struct

    !  INPUTS
    !    none
    !  OUTPUT
    !    o func_struct
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'dealloc_general_function'

    func_struct%func_type = GF_ZERO
    func_struct%param = 0.d0

  end subroutine dealloc_general_function
  !******
  !----------------------------------------------------------------------------
  !****s* general_function/format_general_function
  !  NAME
  !    format_general_function
  !  SYNOPSIS
  subroutine format_general_function(func_struct, func_str)

    !  PURPOSE
    !    Format a general function described in func (see below)
    !    at given times.
    !  USES

    implicit none

    !  ARGUMENTS

    type(gen_func), intent(IN) :: func_struct
    character*(*), intent(OUT) :: func_str

    !  INPUTS
    !    o func_struct -- Internal representation of a function
    !  OUTPUT
    !    o func -- String containing the function.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character*150 :: info_str
    character*20 :: func_type
    real*8 :: period, wave
    character(*), parameter :: func = 'format_general_function'

    select case (func_struct%func_type)
    case(GF_ZERO)
       write(func_str, "('none')")
    case(GF_POWER)
       write(func_str, "('power ',F10.4)") func_struct%param
    case(GF_SIN, GF_COS)
       if (func_struct%func_type == GF_SIN) then
          func_type = 'sin'
       else
          func_type = 'cos'
       end if
       period = (1000 * 2*pi) / func_struct%param
       wave = func_struct%param / (2*pi * c_light_cm_ps)
       write(func_str, "(A,' ',F16.8,' fs    # ',F16.8,' cm^-1')") &
       & trim(func_type), period, wave
    case default
       write(info_str, "(A, 'Unknown function type.')") trim(func)
       call aims_stop(info_str)
    end select

  end subroutine format_general_function
  !******
  !----------------------------------------------------------------------------
  !****s* general_function/eval_general_function
  !  NAME
  !    eval_general_function
  !  SYNOPSIS
  subroutine eval_general_function(func_struct, n_points, values, timestep)

    !  PURPOSE
    !    Evaluate a general function described in func_struct
    !    at given times.
    !  USES

    implicit none

    !  ARGUMENTS

    type(gen_func), intent(IN) :: func_struct
    integer, intent(IN) :: n_points
    real*8, intent(OUT) :: values(0:n_points)
    real*8, intent(IN) :: timestep

    !  INPUTS
    !    o func_struct -- actual function.
    !    o n_points -- number of points to evaluate the functions on.
    !    o times -- times for which to evaluate the function, in ps.
    !  OUTPUT
    !    o values -- values of the function
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character*150 :: info_str
    integer :: i_point, order
    real*8 :: times(0:n_points)
    character(*), parameter :: func = 'eval_general_function'

    do i_point = 0, n_points
       times(i_point) = -i_point
    end do

    select case (func_struct%func_type)
    case(GF_ZERO)
       values = 0.d0   ! Drop this -> least squares stabilize from noise.
    case(GF_POWER)
       values = times**func_struct%param
    case(GF_SIN)
       values = sin(func_struct%param * times * timestep)
    case(GF_COS)
       values = cos(func_struct%param * times * timestep)
    case default
       write(info_str, "(A, 'Unknown function type')") trim(func)
       call aims_stop(info_str)
    end select

  end subroutine eval_general_function
  !******
  !----------------------------------------------------------------------------
  !****s* general_function/number_of_nones
  !  NAME
  !    number_of_nones
  !  SYNOPSIS
  integer function number_of_nones(func_structs)

    !  PURPOSE
    !    Evaluate a general function described in func_struct
    !    at given times.
    !  USES

    implicit none

    !  ARGUMENTS

    type(gen_func), intent(IN) :: func_structs(:)

    !  INPUTS
    !    o func_structs -- general functions.
    !  OUTPUT
    !    o values -- values of the function
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_func

    number_of_nones = 0
    do i_func = 1, size(func_structs)
       if (func_structs(i_func)%func_type == GF_ZERO) then
          number_of_nones = number_of_nones + 1
       end if
    end do

  end function number_of_nones
  !******
end module general_function
!******
