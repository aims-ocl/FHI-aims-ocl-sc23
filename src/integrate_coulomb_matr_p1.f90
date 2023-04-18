!****s* FHI-aims/integrate_coulomb_matr_p1
!  NAME
!    integrate_coulomb_matr_p1
!  SYNOPSIS

subroutine integrate_coulomb_matr_p1 &
      (basis_l_max, coulomb_matr, current_cell_index)

  !  PURPOSE
  !
  !    This subroutine is intended to calculate the coulomb interaction matrix
  !    elements between every two auxiliary basis functions in a periodic
  !    system.  One basis function is sitting in the unit cell at origin,
  !    while the other sits in unit cell given by "current_cell_index".
  !
  !  USES

  use dimensions
  use sbt_overlap_aims
  use prodbas
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: basis_l_max (n_species)
  integer, intent(IN) :: current_cell_index(3)
  real*8, intent(OUT) :: coulomb_matr(n_basbas,n_loc_prodbas)

  !  INPUTS
  !    o basis_l_max -- integer array, the maximal angular momentum number
  !                     of the basis functions for each species (ignored)
  !    o current_cell_index -- unit cell of the column product basis functions
  !  OUTPUTS
  !    o coulomb_matr -- real array, the Coulomb interaction matrix within
  !                      the auxiliary basis   
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  character(*), parameter :: func = 'integrate_coulomb_matr_p1'

  call integrate_auxmat_by_atomic_sbt(coulomb_matr, ovlp_type_bare_or_hse_coul, .false.,&
  &                                   current_cell_index)

end subroutine integrate_coulomb_matr_p1
!******
