!****s* thermodynamic_integration/TDI_update_lambda
!  NAME
!    TDI_update_lambda
!  SYNOPSIS

subroutine TDI_update_lambda()

!  PURPOSE
!   Updates runtime value of lambda
!
!  USES

  use dimensions
  use timing    
  use runtime_choices
  use thermodynamic_integration
  implicit none

!  ARGUMENTS
!  INPUTS
!  OUTPUT
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
!  SOURCE
  
  TDI_lambda  = TDI_lambda + TDI_dlambda_dt*MD_tstep 

end subroutine TDI_update_lambda
!******
