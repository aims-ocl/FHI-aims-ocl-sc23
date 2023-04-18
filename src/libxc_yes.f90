!!****f* FHI-aims/libxc_yes_no
!!
!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Provides the function libxc_yes_no, which returns either true or
!!  false depending on whether aims was compiled with libxc or not.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
!!  SOURCE
!!
pure logical function libxc_yes_no() result(y)
  y = .true.
end function libxc_yes_no
!!***
