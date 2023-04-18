!****s* FHI-aims/tab_local_ylm
!  NAME
!   tab_local_ylm
!  SYNOPSIS

subroutine tab_local_ylm &
     ( trigonom_tab, l_hartree, l_ylm_max, &
     ylm_tab &
     )

!  PURPOSE
!  Subroutine tab_local_ylm computes the Y_lm functions at the present
!  integration point ONLY for the present integration center, no other atoms
!
!  USES

  implicit none

!  ARGUMENTS

  real*8 trigonom_tab ( 4 )
  integer l_hartree
  integer l_ylm_max
  real*8 ylm_tab( (l_ylm_max+1)**2 )

!  INPUTS
!    o trigonom_tab -- trigonometric functions from subroutine tab_trigonom_p0 
!    o basis_l_max -- maximum l index of basis functions
!    o l_ylm_max -- maximum l index of wave functions
!
!  OUTPUT
!    o ylm_tab -- Y_lm functions
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
!  SOURCE



  !  counters

  !  begin work

  !       bootstrap table of ylm's
  call increment_ylm &
       ( trigonom_tab(1), trigonom_tab(2), &
       trigonom_tab(3), trigonom_tab(4), &
       0, l_hartree, ylm_tab )

  !     that's all folks

end subroutine tab_local_ylm

!----------------------------------------------------------------------
!******
