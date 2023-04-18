!****s* FHI-aims/atom_bsse_results
!  NAME
!   atom_bsse_results
!  SYNOPSIS

  subroutine atom_bsse_results

!  PURPOSE
!  Collection of all bsse results
!
!  USES

  use localorb_io, only: use_unit
  use mpi_tasks, only: myid
  use dimensions, only: n_atoms, use_hartree_fock, use_mp2, use_rpa_ene
  use physics
  use timing

!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
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

  implicit none
  
  integer :: i_atom
  real*8 :: energy_rpa, energy_rpa_se, energy
  
  
  ! this is for PBE+RPA calculations
   if(use_rpa_ene) then
      !RPA
      energy_rpa=BSSE_full_energy_RPA
      if(myid.eq.0) then  
          write(use_unit,'(2X,A,f21.8,A)') &
        "  RPA total energy full structure: ", &
         BSSE_full_energy_RPA, " eV"
      endif
      do i_atom = 1, n_atoms
       energy_rpa =  energy_rpa - BSSE_atom_RPA(i_atom)
       if(myid.eq.0) then  
          write(use_unit,'(2X,A,I8,f21.8,A)') &
            "  RPA total energy for atom= ", i_atom, &
            BSSE_atom_RPA(i_atom), " eV"
       endif       
      enddo    
      if(myid.eq.0) then  
         write(use_unit,'(2X,A,f21.8,A)') &
           "  BSSE-corrected RPA atomization energy: ", &
               energy_rpa, " eV"
      endif
      
      ! RPA+SE
      energy_rpa_se=BSSE_full_energy_RPA_SE
      if(myid.eq.0) then  
          write(use_unit,'(2X,A,f21.8,A)') &
        "  RPA+SE total energy full structure: ", &
         BSSE_full_energy_RPA_SE, " eV"
      endif
      do i_atom = 1, n_atoms
       energy_rpa_se =  energy_rpa_se - BSSE_atom_RPA_SE(i_atom)
       if(myid.eq.0) then  
          write(use_unit,'(2X,A,I8,f21.8,A)') &
            "  RPA+SE total energy for atom= ", i_atom, &
            BSSE_atom_RPA_SE(i_atom), " eV"
       endif       
      enddo    
      if(myid.eq.0) then  
         write(use_unit,'(2X,A,f21.8,A)') &
           "  BSSE-corrected RPA+SE atomization energy: ", &
               energy_rpa_se, " eV"
      endif
    ! end of use_rpa_ene
    ! hf+mp2  
   elseif (use_hartree_fock .and. use_mp2) then 
      energy=BSSE_full_energy
      if(myid.eq.0) then  
          write(use_unit,'(2X,A,f21.8,A)') &
        "  HF + MP2  for full structure: ", &
         BSSE_full_energy, " eV"
      endif
      do i_atom = 1, n_atoms
       energy =  energy - BSSE_per_atom(i_atom)
       if(myid.eq.0) then  
          write(use_unit,'(2X,A,I8,f21.8,A)') &
            "  HF + MP2 for atom= ", i_atom, &
            BSSE_per_atom(i_atom), " eV"
       endif       
      enddo    
      if(myid.eq.0) then  
         write(use_unit,'(2X,A,f21.8,A)') &
           "  BSSE-corrected HF+ MP2 atomization energy: ", &
               energy, " eV"
      endif
   endif 
   
   ! this part is common - work still needs to be done to extract the species energy

    !  do i_species = 1, n_species
     !  energy =  energy - species_energy(species_number(i_species))
    !  enddo
      
    !  if(myid.eq.0) then  
    !      write(use_unit,'(2X,A,f21.8,A)') &
    !         "  BSSE-corrected total energy: ", &
    !           energy, " eV"
    !  endif
    
! deallocations :
      ! array that would have been deallocated if atomization bsse had not been used
    
    !  if(allocated(ovlp_3fn))then
    !    deallocate( ovlp_3fn )
    !  end if
      
      ! arrays related directly to bsse calculations:  BSSE_atom_RPA, BSSE_atom_RPA_SE,BSSE_per_atom
      
      if(allocated(BSSE_atom_RPA))then
        deallocate( BSSE_atom_RPA )
      end if
      if(allocated(BSSE_atom_RPA_SE))then
        deallocate( BSSE_atom_RPA_SE )
      end if
      if(allocated(BSSE_per_atom))then
        deallocate( BSSE_per_atom )
      end if
      
  end subroutine atom_bsse_results
!******	
