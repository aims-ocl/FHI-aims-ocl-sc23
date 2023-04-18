!****h* FHI-aims/debugmanager
!  NAME
!    timing
!  SYNOPSIS
      module debugmanager
! USES
      use localorb_io
      use runtime_choices
!  PURPOSE
!    This module controls modular debugging. It stores for which modules
!    debugging has been requested in the inputfile and provides conditional
!    output functions for debug statements.
!  AUTHOR
!    Arvid Conrad Ihrig, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio molecular simulations with numeric atom-centered orbitals"
!    Computer Physics Communications 180, 2175-2196 (2009).
!    http://dx.doi.org/10.1016/j.cpc.2009.06.022
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Dev version, FHI-aims (2014).
!  SOURCE
      implicit none

      public :: init_debuglist, register_debugmodule, activate_debugging, &
                debugprint, module_is_debugged, cleanup_debugmanager

      private
      ! all module-wide symbols following this statement will not be made public

      type debugswitch
         character*150 :: tag
         logical :: active
      end type debugswitch

      type(debugswitch), allocatable :: debuglist(:)
      integer :: debugcount
      character(len=2000) :: info_str

      contains
!******

!****s* debugmanager/init_debuglist
!  NAME
!    init_debuglist
!  SYNOPSIS
      subroutine init_debuglist()
!  PURPOSE
!
!  INPUT
      use mpi_tasks, only: check_allocation
      implicit none
!     none
!  OUTPUT
!     none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  SEE ALSO
!    FHI-aims aimsclub Wiki.
!  SOURCE
      integer :: info
      character(*), parameter :: funcname = "init_debuglist"

      if (allocated(debuglist)) then
         deallocate(debuglist)
      end if
      allocate(debuglist(100), stat=info)
      call check_allocation(info, funcname)
      debugcount = 0

      end subroutine

!****s* debugmanager/register_debugmodule
!  NAME
!    register_debugmodule
!  SYNOPSIS
      subroutine register_debugmodule(tag, active)
!  PURPOSE
!     Enable debug output for the requested module
!  USES
      use mpi_tasks, only: aims_stop
      implicit none
!  INPUT
      character(len=*), intent(in) :: tag
      logical, optional :: active
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  SEE ALSO
!    FHI-aims aimsclub Wiki.
!  SOURCE
      character(*), parameter :: funcname = "register_citation"
      integer :: i
      logical :: initmode
      type(debugswitch) :: newdebug

      if (.not.present(active)) then
         initmode = .false.
      else
         initmode = active
      endif

      newdebug%tag = tag
      newdebug%active = initmode
      do i = 1, debugcount
         if (debuglist(i)%tag == tag) then
            call aims_stop("Debug Module registered twice: "//tag)
         endif
      enddo
      debuglist(debugcount+1) = newdebug
      debugcount = debugcount + 1

      end subroutine register_debugmodule

!****s* debugmanager/activate_debugging
!  NAME
!    activate_debugging
!  SYNOPSIS
      subroutine activate_debugging(tag)
!  PURPOSE
!  USE
      use mpi_tasks, only: aims_stop
      implicit none
!  INPUT
      character(len=*), intent(in) :: tag
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  SEE ALSO
!    FHI-aims aimsclub Wiki.
!  SOURCE
      integer i

      do i = 1, debugcount
         if (debuglist(i)%tag == tag) then
            debuglist(i)%active = .true.
            exit
         elseif (i == debugcount) then
            write(info_str,'(1X,A)') &
            "* ERROR - you requested to debug an unknown module."
            call localorb_info(info_str)
            write(info_str,'(1X,A,A)') "* Unknown module tag: ", tag
            call localorb_info(info_str)
            write(info_str,'(1X,A)') &
               "* If this is no typo, please add it to init_debugging.f90."
            call localorb_info(info_str)
            write(info_str,'(A)') ""
            call localorb_info(info_str)
            call aims_stop("Tried to debug unknown module!")
         endif
      enddo

      end subroutine activate_debugging

!****s* debugmanager/debugprint
!  NAME
!    debugprint
!  SYNOPSIS
      subroutine debugprint(message, tag)
!  PURPOSE
!     Wrapper around localorb_io, only called if debugging for module tag
!     has been requested
!  INPUT
      implicit none
      character(len=*), intent(in) :: message
      character(len=*), intent(in) :: tag
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  SEE ALSO
!    FHI-aims aimsclub Wiki.
!  SOURCE
      if (module_is_debugged(tag)) then
         call localorb_info(message)
      end if

      end subroutine debugprint


!****s* debugmanager/module_is_debugged
!  NAME
!    module_is_debugged
!  SYNOPSIS
     logical function module_is_debugged(tag)
!  PURPOSE
!
!  INPUT
      implicit none
      character(len=*), intent(in) :: tag
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  SEE ALSO
!    FHI-aims aimsclub Wiki.
!  SOURCE
      integer i

      module_is_debugged = .False.
      do i = 1, debugcount
         if (trim(debuglist(i)%tag) == trim(tag)) then
            module_is_debugged = debuglist(i)%active
            exit
         endif
      enddo

      end function module_is_debugged

!******

!****s* debugmanager/cleanup_debugmanager
!  NAME
!    cleanup_debugmanager
!  SYNOPSIS
      subroutine cleanup_debugmanager()
!  PURPOSE
!     Memory cleanup for debugmanager
!  USES
        implicit none
!  INPUT
!     none
!  OUTPUT
!     none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  SEE ALSO
!    FHI-aims aimsclub Wiki.
!  SOURCE
        if (allocated(debuglist)) deallocate(debuglist)
      end subroutine
!******


   end module debugmanager
