!****h* FHI-aims/DFPT_phonon/DFPT_phonon_supercell
!  NAME
!    DFPT_phonon_supercell -- module for DFPT_phonon_supercell
!  SYNOPSIS
module DFPT_phonon_supercell
!  PURPOSE
!    This module contains the supercell arrays. 
!  USES
!  use runtime_choices, only : use_scalapack

  implicit none

   real*8, allocatable :: overlap_supercell(:,:)
   real*8, allocatable :: hamiltonian_supercell(:,:)
   real*8, allocatable :: KS_eigenvalue_supercell(:)
   real*8, allocatable :: KS_eigenvector_supercell(:,:)

end module DFPT_phonon_supercell
