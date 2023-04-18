!  subroutine to check whether coordinates are
!  in cavity radius
!
!  R. Gehrke (2005)
!

subroutine check_cavity(coords, cavity_radius, outside)

  use vector_array_class
  use cluster

  implicit none

!  imported variables

!  input
  type (vector_array), intent(in) :: coords
  real*8, intent(in) :: cavity_radius

!  output
  logical, intent(out) :: outside

!  counter
  integer :: i_atom

  outside = .false.
  do i_atom = 1, n_atoms, 1
     if (norm(coords%coords(i_atom)) .ge. cavity_radius) then
        outside = .true.
     end if
  end do
  
end subroutine check_cavity
