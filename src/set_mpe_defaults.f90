!****s* FHI-aims/set_mpe_defaults
!  NAME
!    set_mpe_defaults
!  SYNOPSIS
subroutine set_mpe_defaults ( )
!  PURPOSE
!    Sets default values for modules used in MPE solver
!  USES
   use mpe_dielectric_interfaces, only: never_updated_before
   use mpe_interface, only: scf_step_of_first_call, MPEEnergyContributions, &
       mpe_energy_contributions, mpe_interface_initialized => module_initialized
   use mpe_types, only: dummy_mpi, &
       mpe_types_initialized => module_initialized
   use types, only: dp
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
!    This ugly hack of a subroutine resets the module variables which control
!    MPE module initialization, needed when running aims as a library.
!
!    Ordinarily I would do this via explicit calls to destructor methods, but
!    the MPE module relationship is complicated enough that I've decided to
!    mimic the original logic via a subroutine instead.
!  SOURCE
   never_updated_before = .true.
   scf_step_of_first_call = -1
   mpe_interface_initialized = .false.
   dummy_mpi = .false.
   mpe_types_initialized = .false.

   mpe_energy_contributions % rho_V = 0.e0_dp
   mpe_energy_contributions % nuc_V = 0.e0_dp
   mpe_energy_contributions % nonel = 0.e0_dp
end subroutine set_mpe_defaults
!******
