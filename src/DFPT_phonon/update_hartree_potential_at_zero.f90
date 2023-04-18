!****s* FHI-aims/update_hartree_potential_shanghui_p0
!  NAME
!    update_hartree_potential_shanghui_p0
!  SYNOPSIS

subroutine update_hartree_potential_at_zero &
     ( ) 

!  PURPOSE
!  
!  use integrate_hartree_log_grid to  give  some point at ZERO.  
!  
!  shanghui 2015
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use spline
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use constants
  use hartree_potential_storage , only :  get_rho_multipole_supercell_spl, &  
                                 delta_v_hartree_part_at_zero_supercell, & 
                                 delta_v_hartree_deriv_l0_at_zero_supercell 
  use pbc_lists, only : centers_in_hamiltonian,center_to_atom, &
                        center_in_sc_DFPT_to_atom
  implicit none

!  ARGUMENTS
!
!  OUTPUT
!   o  delta_v_hartree_part_at_zero_supercell --  multipole expansion of the Hartree potential at origin of atoms
!   o  delta_v_hartree_deriv_l0_at_zero_supercell -- derivates of  multipole expansion of the Hartree potential at origin of atoms
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


  real*8, allocatable :: current_rho_multipole_spl(:,:,:)

!  counters

  integer i_atom_multipole,i_center_multipole
  integer i_index 


  character*100 :: info_str


!  begin work


  allocate(current_rho_multipole_spl((l_pot_max+1)**2, n_max_spline, n_max_radial+2),stat=i_index)
  call check_allocation(i_index, 'current_rho_multipole_spl')


  do i_center_multipole = 1, n_centers_in_sc_DFPT
 
     i_atom_multipole = center_in_sc_DFPT_to_atom(i_center_multipole)
 
     call get_rho_multipole_supercell_spl(current_rho_multipole_spl, i_center_multipole, i_atom_multipole)
 
     call integrate_hartree_log_grid_supercell(i_atom_multipole,current_rho_multipole_spl, &
                delta_v_hartree_part_at_zero_supercell(i_center_multipole), &
                delta_v_hartree_deriv_l0_at_zero_supercell(1:3,i_center_multipole))  
 
  enddo  ! i_center_multipole


  deallocate(current_rho_multipole_spl)

end subroutine  update_hartree_potential_at_zero
!******
