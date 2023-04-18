!****h* FHI-aims/check_environment
!  NAME
!    check_environment
!  SYNOPSIS

module check_environment

  !  PURPOSE
  !
  !    Provide calls to check the settings of environmental variables.
  !
  !    PLEASE - no dependencies on other modules if at all possible. Pass
  !             arguments explicitly if needed.
  !
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
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  private

  public check_environment_variable

  interface check_environment_variable
     module procedure check_env_var_int, check_env_var_string
  end interface check_environment_variable

contains

  !----------------------------------------------------------------------------
  !****s* check_environmet/check_env_var_int
  !  NAME
  !    check_env_var_int
  !  SYNOPSIS

  function check_env_var_int(varname,value,myid) result(lset)

    !  PURPOSE
    !    Check ${varname} to have integer "value". 
    !  USES

    USE localorb_io, only: use_unit
    implicit none

    !  ARGUMENTS

    character(*), intent(IN) :: varname
    integer, intent(in) :: value
    integer, intent(in) :: myid

    !  INPUTS
    !    o varname -- Name of environment variable
    !    o value -- Value it should have
    !    o myid -- Process id
    !  OUTPUTS
    !    o lset -- .true. if correctly set.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    character(len=100)  :: result
    integer             :: stat, converted_string
    logical             :: lset
    character(*), parameter :: func = 'check_env_var_int'

    lset=.false.

    call get_environment_variable(varname, value=result, STATUS=stat)

    if  (stat /= 0) then
       ! Environment variable varname is not set. Notify user
       if (myid .eq. 0) then
          write(use_unit,'(2X,a,a,a)') "*** Environment variable ",trim(varname), " is not set"
          write(use_unit,'(2X,a,i0)') "*** For performance reasons you might want to set it to ",value
       endif
    else
       ! environment variable is set to result

       ! we want to compare the result (stored in a string) with
       ! the given value
       ! => convert sting to integer
       read ( result, '(i10)') converted_string

       if (converted_string .ne. value) then
          if (myid .eq. 0) then
             write(use_unit,'(2X,a,a,a)') "*** Environment variable ",trim(varname), " is set to ",trim(result)
             write(use_unit,'(2X,a,i0)') "*** For performance reasons you might want to set it to ",value
          endif
       else
          ! variable is set to correct value; report this back
          lset = .true.
          if (myid .eq. 0) then
             write(use_unit,'(2X,5a)') "| Environment variable ", trim(varname), " correctly set to ", trim(result), "."
          end if
       endif
    endif  ! stat /= 0
  end function check_env_var_int
  !******
  !----------------------------------------------------------------------------
  !****s* check_environmet/check_env_var_string
  !  NAME
  !    check_env_var_string
  !  SYNOPSIS

  function check_env_var_string(varname,value,myid) result(lset)

    !  PURPOSE
    !    Check ${varname} to have string "value". 
    !  USES

    USE localorb_io, only: use_unit
    implicit none

    !  ARGUMENTS

    character(*), intent(IN) :: varname
    character(*), intent(IN) :: value
    integer, intent(in) :: myid

    !  INPUTS
    !    o varname -- Name of environment variable
    !    o value -- Value it should have
    !    o myid -- Process id
    !  OUTPUTS
    !    o lset -- .true. if correctly set.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    character(len=100)  :: result
    integer             :: stat
    logical :: lset
    character(*), parameter :: func = 'check_env_var_string'

    lset = .false.

    call get_environment_variable(varname, value=result, STATUS=stat)

    if  (stat /=0 ) then
       ! environment variable varname is not set; notify user
       if (myid .eq. 0) then
          write(use_unit,*) " *** Environment variable ",trim(varname), " is not set"
          write(use_unit,'(a,a)') " *** For performance reasons you might want to set it to ",value
       endif
    else
       ! environment variable is set to result

       if (result .ne. value) then
          if (myid .eq. 0) then
             write(use_unit, '(a,a,a)') " *** Environment variable ",trim(varname), " is set to ",trim(result)
             write(use_unit, '(a,i0)') " *** For performance reasons you might want to set it to ",value
          endif
       else
          ! variable is set to correct value; report this back
          lset = .true.
          if (myid .eq. 0) then
             write(use_unit, '(5a)') "  | Environment variable ", trim(varname), " correctly set to ", trim(result), "."
          end if
       endif
    end if

  end function check_env_var_string
  !******

end module check_environment
!******

