module python_interface
!****h* FHI-aims/xml_write
! PURPOSE
!   Stub files when AIMS is not linked to CFFI
! AUTHOR
!   Jan Hermann
! CREATION DATE
!   2016-04-24
! ADDITIONAL INFO
!   See python_interface.f90.
! COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e.V. Please note
!   that any use of the "FHI-aims-Software" is subject to the terms and
!   conditions of the respective license agreement.
!******
use localorb_io, only: localorb_info

implicit none

private

public :: register_python_hook, run_python_hook

type :: Hook_t
    logical, public :: registered = .false.
end type

type :: HookRegister_t
    type(Hook_t) :: post_scf
    type(Hook_t) :: post_hirshfeld
end type

type(HookRegister_t), public :: python_hooks

contains

integer function register_python_hook(label, filename, attrib) result(retcode)
    character(len=*), intent(in) :: label, filename, attrib

    call localorb_info('*** Error: python_hook requires cffi')
    retcode = -1
end function

subroutine run_python_hook(hook)
    type(Hook_t), intent(in) :: hook
end subroutine

end module
