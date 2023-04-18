!****s* FHI-aims/k_point_symmetry_check
!  NAME
!    k_point_symmetry_check
!  SYNOPSIS
subroutine k_point_symmetry_check(n_k_xyz, k_off, k_number)

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

  use localorb_io, only: localorb_info
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

  integer, parameter :: n_max_symmetries = 48
  real*8 :: symmats(3, 3, n_max_symmetries), origins(3, n_max_symmetries)
  integer :: n_symmetries
  character*140 :: info_str
  integer :: i
  character(*), parameter :: func = 'k_point_symmetry_check'

  write(info_str,"(2X,A)") "Using symmetry for reducing the k-points"
  call localorb_info(info_str)

  ! JW: The idea would be the following:
  ! !  call get_system_symmetries(n_max_symmetries, n_symmetries, &
  ! !  &                          symmats, origins)
  ! But for now, we can only use time reversal.
  n_symmetries = 2
  symmats(:,:,1:2) = 0.d0
  do i = 1, 3
     symmats(i, i, 1) =  1.d0   ! identity
     symmats(i, i, 2) = -1.d0   ! inversion
  end do
  ! JW: Please note that one should in general generate a group from both the
  ! system symmetries and the inversion symmetries.

  call symmetry_reduce_k_points(n_symmetries, symmats, &
  &                             n_k_xyz, k_off, k_number)

end subroutine k_point_symmetry_check
!******
