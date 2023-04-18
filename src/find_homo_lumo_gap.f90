!****s* FHI-aims/find_homo_lumo_gap
!  NAME
!   find_homo_lumo_gap
!  SYNOPSIS
subroutine find_homo_lumo_gap &
      ( KS_eigenvalue, occ_numbers, homo_level, lumo_level, homo_occ, lumo_occ,&
        i_kpt_homo, i_kpt_lumo, i_spin_homo, i_spin_lumo, found_min_direct_gap,&
        min_direct_gap, i_kpt_min_direct_gap, i_spin_min_direct_homo, &
        i_spin_min_direct_lumo )
!  PURPOSE
!    Find the HOMO level, LUMO level, and minimum direct gap and returns them.
!    Does not output to screen.  To avoid unintential modification of state, no
!    module variables should be modified.
!  USES
  use dimensions, only : n_states, n_spin, spin_degeneracy, n_k_points
  use constants, only : hartree
  implicit none

  real*8,  intent(in)  :: KS_eigenvalue(n_states, n_spin, n_k_points)
  real*8,  intent(in)  :: occ_numbers(n_states, n_spin, n_k_points)
  real*8,  intent(out) :: homo_level
  real*8,  intent(out) :: lumo_level
  real*8,  intent(out) :: homo_occ
  real*8,  intent(out) :: lumo_occ
  integer, intent(out) :: i_kpt_homo
  integer, intent(out) :: i_kpt_lumo
  integer, intent(out) :: i_spin_homo
  integer, intent(out) :: i_spin_lumo
  logical, intent(out) :: found_min_direct_gap
  real*8,  intent(out) :: min_direct_gap
  integer, intent(out) :: i_kpt_min_direct_gap
  integer, intent(out) :: i_spin_min_direct_homo
  integer, intent(out) :: i_spin_min_direct_lumo

  !  counters
  real*8  :: direct_gap
  real*8 :: current_homo_level, current_lumo_level
  real*8 :: midpoint
  integer :: current_homo_spin, current_lumo_spin
  integer :: current_homo_state, current_lumo_state
  integer :: i_state_homo, i_state_lumo
  integer :: i_state, i_spin, i_k_point

  ! Determine HOMO and LUMO values (VBM and CBM in case of periodic systems)
  ! safe initial values. If we break these, we have a problem somwehere else.
  homo_level = -10000000.0d0
  lumo_level = 10000000.0d0

  ! We also look for the minimum direct gap value at a given k-point.
  min_direct_gap = 10000000.0d0
  found_min_direct_gap = .false.

  ! Define the correct "half occupation" (with or without spin)
  midpoint = spin_degeneracy/2.0d0

  homo_occ = 2.0d0
  lumo_occ = 0.0d0
  i_kpt_homo = 0
  i_kpt_lumo = 0
  do i_k_point = 1, n_k_points, 1
    ! The "current" variables refer to HOMO, LUMO, etc. at the current k-point
    ! only.
    ! Their purpose is solely the determination of the direct gap.
    ! They are thus only meaningful in periodic systems. We do not
    ! expect any significant overhead in non-periodic systems, as searching
    ! for the direct gap is then a waste of time.
    current_homo_level = -10000000.0d0
    current_lumo_level =  10000000.0d0
    current_homo_spin = 0
    current_lumo_spin = 0

    do i_spin = 1, n_spin, 1
      do i_state = 1, n_states, 1
        ! We first search for the global HOMO and LUMO (any k-point)
        if (occ_numbers(i_state, i_spin, i_k_point) .ge. midpoint) then
          ! check if homo
          ! "HOMO" also includes Fermi level ("ge" above)
          if (KS_eigenvalue(i_state, i_spin, i_k_point) .gt. homo_level) then
            homo_level = KS_eigenvalue(i_state, i_spin, i_k_point)
            homo_occ = occ_numbers(i_state, i_spin, i_k_point)
            i_kpt_homo = i_k_point
            i_spin_homo = i_spin
            i_state_homo = i_state
          end if
        end if

        if (occ_numbers(i_state, i_spin, i_k_point) .le. midpoint) then
          ! check if lumo
          ! LUMO must also include Fermi level ("le" above), else we may get
          ! nonsensical gaps (i.e., a gap in a molecule with half-occupied
          ! orbitals)
          if (KS_eigenvalue(i_state, i_spin, i_k_point) .lt. lumo_level) then
            lumo_level = KS_eigenvalue(i_state, i_spin, i_k_point)
            lumo_occ = occ_numbers(i_state, i_spin, i_k_point)
            i_kpt_lumo = i_k_point
            i_spin_lumo = i_spin
            i_state_lumo = i_state
          end if
        end if

        if (n_k_points.gt.1) then
          ! We next do the same thing again, but this time we search the direct
          ! gap at the present k-point.
          ! Tricky enough, the direct gap could be between HOMO and LUMO on
          ! different spin channels.
          if (occ_numbers(i_state, i_spin, i_k_point) .ge. midpoint) then
            ! check if homo
            ! "HOMO" also includes Fermi level ("ge" above)
            if (KS_eigenvalue(i_state, i_spin, i_k_point) .gt. current_homo_level) then
              current_homo_level   = KS_eigenvalue(i_state, i_spin, i_k_point)
              current_homo_spin = i_spin
              current_homo_state = i_state
            end if
          end if
          if (occ_numbers(i_state, i_spin, i_k_point) .le. midpoint) then
            ! check if lumo
            ! LUMO must also include Fermi level ("le" above), else we may get
            ! nonsensical gaps (i.e., a gap in a molecule with half-occupied
            ! orbitals)
            if (KS_eigenvalue(i_state, i_spin, i_k_point) .lt. current_lumo_level) then
              current_lumo_level   = KS_eigenvalue(i_state,i_spin,i_k_point)
              current_lumo_spin = i_spin
              current_lumo_state = i_state
            end if
          end if
        end if ! ( n_k_points .gt. 1)
      end do
    end do

    ! if we have more than one k-point, check for the minimum direct gap here:
    if (n_k_points.gt.1) then
      if ( (current_lumo_spin .ne. 0) .and. (current_homo_spin .ne. 0) ) then
        direct_gap = current_lumo_level - current_homo_level
        if (direct_gap .lt. min_direct_gap) then
          min_direct_gap = direct_gap
          found_min_direct_gap = .true.
          i_spin_min_direct_homo = current_homo_spin
          i_spin_min_direct_lumo = current_lumo_spin
          i_kpt_min_direct_gap = i_k_point
        end if
      end if
    end if
  end do

  ! To avoid unintentional modification of state, estimate_low_gap is modified
  ! in find_and_output_homo_lumo_gap, not here.
end subroutine find_homo_lumo_gap
