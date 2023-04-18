!****s* FHI-aims/set_aims_defaults
!  NAME
!    set_aims_defaults
!  SYNOPSIS
subroutine set_aims_defaults ( )
!  PURPOSE
!    Set an FHI-aims run to a default (uninitialized!) state by calling various
!    module subroutines which set default values for their modules
!  USES
   use dimensions, only: set_dimensions_defaults
   use geometry, only: set_geometry_defaults
   use grids, only: set_grids_defaults
   use localorb_io, only: set_localorb_io_defaults
   use pbc_lists, only: deallocate_pbc_lists
   use physics, only: set_physics_defaults
   use prodbas, only: cleanup_basbas
   use runtime_choices, only: set_runtime_choices_defaults
   use spline, only: cleanup_spline
   implicit none
!  ARGUMENTS
!    o None
!  INPUTS
!    o None
!  OUTPUT
!    o None
!  AUTHOR
!    William Huhn (Duke University)
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted and accepted and
!    published for about 11 years now.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    May 2019 - Created.
!  NOTES
!    In principle, every module which includes module variables should be have
!    a subroutine which is called here. 
!  SOURCE
   call set_dimensions_defaults ( )
   call set_geometry_defaults ( )
   call set_grids_defaults ( )
   call set_localorb_io_defaults ( )
   call deallocate_pbc_lists ( )
   call set_physics_defaults ( )
   call cleanup_basbas ( )
   call set_runtime_choices_defaults ( )
   call cleanup_spline ( )
   call set_mpe_defaults ( )
end subroutine
!******
