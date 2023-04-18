!****s* FHI-aims/init_metadyn_
!  NAME
!   init_metadyn_
!  SYNOPSIS
subroutine  init_metadyn_(n_atoms_p, MD_tstep, masses_p, charges_p,  cell_type , boltzmann_kB_p ,  plumed_file)
! PURPOSE
!   This is a stub, i.e. a place holder. It ensures that AIMS can be 
!   compiled without a C compiler if metadynamics PLUMED is not requested.
!   As such it does nothing (except maybe to end the program, showing to a 
!   programmer that something went wrong). 
!   If you ended up here, you likely mixed up some input options.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals",
!    Computer Physics Communications (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2009).
!  SOURCE
  use dimensions
  use mpi_tasks, only: STDERR
!  use constants
  implicit none
  
  
  integer 	:: n_atoms_p, cell_type, boltzmann_kB_p
  real*8  	:: MD_tstep
  character*40 	:: plumed_file
  real*8, dimension(n_atoms)   :: masses_p , charges_p
  
  write(STDERR,*) "* Attention. Your run has called a subroutine intended for "
  write(STDERR,*) "* the metadynamics plugin called PLUMED"
  write(STDERR,*) "* This functionality must be compiled into the code using a "
  write(STDERR,*) "* subroutine written in C. The code version that you are using "
  write(STDERR,*) "* was not compiled with support PLUMED "
  write(STDERR,*) "* "
  write(STDERR,*) "* Either modify your control.in file, or build/use the code version "
  write(STDERR,*) "* including PLUMED"
  write(STDERR,*) "* Stopping the code for now. "
  stop
      
end subroutine init_metadyn_
!******

!****s* FHI-aims/meta_force_calculation
!  NAME
!   meta_force_calculation
!  SYNOPSIS 
subroutine meta_force_calculation_( lattice_vector, MD_stepcount , px, py, pz , fx, fy, fz,  plumed_energy) 
! PURPOSE
!   This is a stub, i.e. a place holder. It ensures that AIMS can be 
!   compiled without a C compiler if metadynamics PLUMED is not requested.
!   As such it does nothing (except maybe to end the program, showing to a 
!   programmer that something went wrong). 
!   If you ended up here, you likely mixed up some input options.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals",
!    Computer Physics Communications (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2009).
!  SOURCE
  use mpi_tasks, only: STDERR
  use dimensions	
  implicit none
  
  
  integer 			:: MD_stepcount 
  real*8, dimension(3,3)	:: lattice_vector
  real*8, dimension(1,n_atoms)  :: px, py, pz , fx, fy, fz
  real*8 :: plumed_energy 
  
  write(STDERR,*) "* Attention. Your run has called a subroutine intended for "
  write(STDERR,*) "* the metadynamics plugin called PLUMED"
  write(STDERR,*) "* This functionality must be compiled into the code using a "
  write(STDERR,*) "* subroutine written in C. The code version that you are using "
  write(STDERR,*) "* was not compiled with support PLUMED "
  write(STDERR,*) "* "
  write(STDERR,*) "* Either modify your control.in file, or build/use the code version "
  write(STDERR,*) "* including PLUMED"
  write(STDERR,*) "* Stopping the code for now. "
  stop
      
end subroutine meta_force_calculation_
!******
