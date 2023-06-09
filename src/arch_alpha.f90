
!****h* FHI-aims/arch_alpha
!  NAME
!    arch_specific -- arch_alpha
!  SYNOPSIS

module arch_specific

  !  PURPOSE
  !
  ! This module contains any architecture specific 
  ! calls for the "alpha" case.
  !
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !
  !******

  implicit none

contains
  !
  !-------------------------------------------------------------------------
  !****f* arch_alpha/arch_erf
  !  NAME
  !    arch_erf
  !  SYNOPSIS

  real*8 function arch_erf(argument)

    !  PURPOSE
    !  The function calls Erf(x)  function 
    !
    implicit none
    !  ARGUMENTS

    real*8 :: argument

    !  INPUTS
    !    o argument -- argument of erf-function
    !  OUTPUT
    !    o arch_erf -- value of erf-function
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    ! SOURCE


    real*8, external :: derf

    arch_erf = derf(argument)

    return
  end function arch_erf
    !******
  !-------------------------------------------------------------------------
  !****f* arch_alpha/arch_erfc
  !  NAME
  !    arch_erfc
  !  SYNOPSIS

  real*8 function arch_erfc(argument)

    !  PURPOSE
    !  The function calls Erfc(x)  function 
    !
    implicit none
    !  ARGUMENTS

    real*8 :: argument

    !  INPUTS
    !    o argument -- argiment of erfc-function
    !  OUTPUT
    !    o arch_erf -- value of erfc-function
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    ! SOURCE


    real*8, external :: derfc

    arch_erfc = derfc(argument)

    return
  end function arch_erfc
    !******
  !-------------------------------------------------------------------------
end module arch_specific
