!****s* FHI-aims/spglib_symm_stub.f90
!  NAME
!    spglib_symm_stub.f90
!  SYNOPSIS


module spglib_symmetry

!  PURPOSE
  !
  !    Stub module in case spglib has not been compiled
  !
  !  Stub Subroutines:
  !  o get_symmetry
  !  o write_symm_info
  !  o symmetrize_forces
  !  o ...
  !
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
  !    Release version, FHI-aims (2011).

  !  USES

  use mpi_tasks, only: aims_stop
  use localorb_io

  ! Variables
  implicit none
  integer, dimension(:,:,:), allocatable :: spg_rotations
  real*8, dimension(:,:),  allocatable :: spg_shift
  integer :: num_symmetry
  integer, dimension(:),  allocatable :: map
  integer, dimension(:),  allocatable :: map_inv
  integer, dimension(:),  allocatable :: map_sym
  integer, dimension(:), allocatable :: spg_partition_tab
  integer :: n_k_points_nosym
  integer :: n_k_points_xyz_nosym(3)


contains

subroutine check_spglib(flag_spg_lib)

  !  PURPOSE
  !
  !    Check if spglib has been compiled
  !

  implicit none

  !  ARGUMENTS
  !    flag_spg_lib : boolean here set to false since we are in the stub module
  !  INPUTS
     logical, intent(INOUT) :: flag_spg_lib
  !  OUTPUTS
  !    none
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
  !    Release version, FHI-aims (2011).
  !  SOURCE

  ! Local Variables

  flag_spg_lib = .false.

end subroutine check_spglib


subroutine get_symmetry()

  implicit none

  call aims_stop("SPGlib needed","get_symmetry")

end subroutine get_symmetry
!******


subroutine write_symm_info()

  implicit none

  call aims_stop("SPGlib needed","write_symm_info")

end subroutine write_symm_info
!******

subroutine get_symmetry_for_lattice &
           (lattice_vector, n_atoms, species, coords, symprec, &
            spacegroup_found, n_operations_sym, spacegroup_number, &
            international_symbol, schoenflies)
  use, intrinsic :: iso_c_binding, only: c_char
  implicit none

  real*8 , intent(in)  :: lattice_vector(3,3)
  integer, intent(in)  :: n_atoms
  integer, intent(in)  :: species(n_atoms)
  real*8,  intent(in)  :: coords(3, n_atoms)
  real*8,  intent(out) :: symprec
  logical, intent(out) :: spacegroup_found
  integer, intent(out) :: n_operations_sym
  integer, intent(out) :: spacegroup_number
  character(len=14, kind=c_char), intent(out) :: international_symbol
  character(len=10, kind=c_char), intent(out) :: schoenflies

  call aims_stop("SPGlib needed","get_symmetry_for_lattice")

end subroutine get_symmetry_for_lattice
!******

subroutine symmetrize_forces(arg)

  implicit none

  real*8 :: arg(*)

  call aims_stop("SPGlib needed","symmetrize_forces")

end subroutine symmetrize_forces


subroutine destroy_symmetry()

  implicit none

  call aims_stop("SPGlib needed","destroy_symmetry")

end subroutine destroy_symmetry

subroutine destroy_symmats()
  implicit none

  call aims_stop("SPGlib needed", "destroy_symmats")
end subroutine destroy_symmats

!****s* FHI-aims/k_point_symmetry_check
!  NAME
!    k_point_symmetry_check
!  SYNOPSIS
subroutine k_point_symmetry_check_spg(n_k_xyz, k_off, k_number)

  !  PURPOSE
  !
  !    Idea: Figure out system symmetries and reduce the k-points accordingly.
  !
  !    Reduction is done solely done by putting positive values into k_number
  !    entries of irreducible k-points and zero out entries of redundant
  !    k-points.
  !
  !    JW: With no symmetry operations on basis functions (and real-valued
  !    spherical harmonics) at hand, making use of system symmetry is not
  !    trivial.  But time inversion is a different issue, so use only that.
  !
  !  USES
  use localorb_io
  implicit none
  !  ARGUMENTS

  integer, intent(IN) :: n_k_xyz(3)
  real*8, intent(IN) :: k_off(3)
  integer, intent(OUT) :: k_number(n_k_xyz(1), n_k_xyz(2), n_k_xyz(3))
  !  INPUTS
  !    o n_k_xyz -- number of original k-points in the three directions
  !    o k_off -- k-point offsets (parameter k_offset in control.in)
  !  OUTPUT
  !    o  k_number(i_k_xyz(1), i_k_xyz(2), i_k_xyz(3))
  !       -- contains the number of k-points that have been mapped on to the
  !       original entries, this is required to determine the proper k-weights
  !       in the end
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

  character*140 :: info_str

  write(info_str,"(10X,A)") "SPGlib not found! Using regular symmetry for reducing the k-points"
  call localorb_info(info_str)

  call k_point_symmetry_check(n_k_xyz, k_off, k_number)


end subroutine k_point_symmetry_check_spg

subroutine out_symm_mats()

  implicit none

  call aims_stop("SPGlib needed","out_symm_mats")

end subroutine out_symm_mats

end module spglib_symmetry
