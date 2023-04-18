!****s* FHI-aims/DFPT_phonon/tab_two_center_coords_PBC
!  NAME
!   tab_two_center_coords_PBC
!  SYNOPSIS

subroutine tab_two_center_coords_PBC( i_center_1, i_center_2, &
     dist_tab_sq, dir_tab )

!  PURPOSE
!  Subroutine tab_atom_centered_coords tabulates the: 
!  R(i_center_2)-R(i_center_1) under PBC. 
!
!  USES

  use dimensions
  use pbc_lists
  implicit none

!  ARGUMENTS

      integer :: i_center_1
      integer :: i_center_2
      real*8 dist_tab_sq
      real*8 dir_tab ( 3 )

!  INPUTS
!    o i_center_1 -- 
!    o i_center_2 -- 
!  OUTPUT
!    o dist_tab_sq -- (distance from atom center)**2
!    o dir_tab --R(i_center_2)-R(i_center_1) 
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
!

	

!  local variables

!     counters

      integer :: i_coord,i_cell_1,i_cell_2,i_cell_delta,i_center_2_sc_DFPT, i_center_2_new

!  begin work

!       tabulate dist_tab, dir_tab first

!     dir_tab = (u2+R2)-(u1+R1) = (u2+R2-R1) - u1
 
      i_cell_1 = center_in_sc_DFPT_to_cell_in_sc_DFPT(i_center_1)
      i_cell_2 = center_in_sc_DFPT_to_cell_in_sc_DFPT(i_center_2)
 
      i_cell_delta = cell_diff_sc_DFPT(i_cell_2,i_cell_1)

      i_center_2_sc_DFPT = cell_and_atom_to_center_sc_DFPT(i_cell_delta,  &  
                                 center_in_sc_DFPT_to_atom(i_center_2) ) 

      i_center_2_new = center_in_sc_DFPT_to_center(i_center_2_sc_DFPT) 

     do i_coord = 1,3,1
        dir_tab(i_coord) = coords_center(i_coord,i_center_2_new) - &
                           coords_center(i_coord,center_in_sc_DFPT_to_atom(i_center_1))
     end do


     dist_tab_sq = 0.d0
     do i_coord = 1,3,1
        dist_tab_sq = dist_tab_sq + dir_tab(i_coord)**2.0d0         
     enddo

!  that's it

    end subroutine tab_two_center_coords_PBC
!----------------------------------------------------------------------
!******
