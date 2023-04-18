!****s* FHI-aims/classical_field
!  NAME
!    classical_field
!  SYNOPSIS

    subroutine classical_field ( )

!  PURPOSE
!  High-level wrapper around the classical force field calculation
!  Here, we call for the calculation of (additional) classical fields.
!
!  USES
      use physics
      use thermodynamic_integration
      use dimensions
      use localorb_io
      use runtime_choices
      use plumed_new_interface
      implicit none

!  ARGUMENTS

!  INPUTS
!    o none
!  OUTPUTS
!    o none
! 
!  AUTHOR
!    Christian Carbogno,  Fritz-Haber Institute of the Max-Planck-Society
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
!    Development version, FHI-aims (2010).
!  SOURCE
!

    !! Classical harmonic field for the thermodynamic integration
    ! Please take care that the correct forces and energy combination
    ! is specified in predict_new_geometry as well.
    ! IF NOT, the entropy corrections might enter twice
    if (use_thermodynamic_integration) then


      call localorb_info("", use_unit,'(2X,A)')
      call localorb_info("-----------------------------------------------------------", use_unit,'(2X,A)')
      if (.not. use_harmonic_pot_only) then
        call localorb_info("Thermodynamic integration:                                 ", use_unit,'(2X,A)')
      end if
      call localorb_info("| Calculating the quasi-harmonic potential for this geometry", use_unit,'(2X,A)')
      call update_QH_potential()
      call localorb_info("-----------------------------------------------------------", use_unit,'(2X,A)')

      ! Use free energy:
      TDI_ANH_energy      = total_energy + 2.d0*entropy_correction
      
      ! Setup the hybrid system for the thermodynamic integration:
      total_energy      = ( TDI_lambda*TDI_ANH_energy )    + ( (1d0 - TDI_lambda)*( TDI_atoms_energy+TDI_Segment_energy_offset ) )
      total_forces(:,:) = ( TDI_lambda*total_forces(:,:) ) + ( (1d0 - TDI_lambda)*  TDI_atoms_forces(:,:) )

    else if (use_reffree_AS) then
      ! Use free energy:
      TDI_ANH_energy      = total_energy + 2.d0*entropy_correction
      
      ! Setup the hybrid system for the thermodynamic integration:
      total_energy      = ( TDI_lambda*TDI_ANH_energy )    
      total_forces(:,:) = ( TDI_lambda*total_forces(:,:) ) 

    else if (plumed_new_plugin) then
       call plumed_forces ()  
    end if

    end subroutine classical_field
!******
