!****s* FHI-aims/trans_hessian_to_dynamical_matrix
!  NAME
!    trans_hessian_to_dynamical_matrix
!  SYNOPSIS

subroutine trans_hessian_to_dynamical_matrix(hessian, & 
                                            dynamical_matrix,q_phase_band)

!  PURPOSE
!    hessian(3,n_centers_in_hamiltonian,3,n_atoms) ===> 
!    dynamical_matrix(3*n_atoms,3*n_atoms) @ i_q_point
!
!  USES

  use dimensions
  use pbc_lists
  implicit none

!  ARGUMENTS

  real*8 , dimension(3,n_centers_in_sc_DFPT,3,n_atoms), intent(in) :: hessian
  complex*16 ,dimension(3*n_atoms,3*n_atoms), intent(out) ::dynamical_matrix
  complex*16 ,dimension(n_cells_in_sc_DFPT), intent(in) :: q_phase_band

!  INPUTS
!    o hessian
!
!  OUTPUT
!   o dynamical_matrix 
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
! SOURCE


  integer :: i_center, i_cell, i_atom, j_atom
  integer :: i_coord, j_coord

  !-------------initial-------------------- 
  dynamical_matrix = (0.d0,0.0d0)


   do i_center = 1,n_centers_in_sc_DFPT
      i_atom   = center_in_sc_DFPT_to_atom(i_center)
      i_cell   = center_in_sc_DFPT_to_cell_in_sc_DFPT(i_center)

   do j_atom = 1,n_atoms
 
      do i_coord = 1,3 
      do j_coord = 1,3
 
         dynamical_matrix(3*i_atom+i_coord-3, 3*j_atom+j_coord-3) =  & 
         dynamical_matrix(3*i_atom+i_coord-3, 3*j_atom+j_coord-3) +  & 
         hessian(i_coord,i_center,j_coord,j_atom)*q_phase_band(i_cell) 

      enddo 
      enddo

   enddo ! j_atom
   enddo ! i_center


end subroutine trans_hessian_to_dynamical_matrix
!******
