! subroutine to shift center of coordinates of an array
! to the origin
!
! R. Gehrke (2005)
!

  subroutine center_structure(this)

    use vector_array_class
    use cluster
    
    implicit none

!  imported variables

!  input
    type (vector_array), intent(inout) :: this

!  local variables
    type (vector) :: center
    real*8 :: inv_n_atoms

!  counter
    integer :: i_atom
    inv_n_atoms = 1. / n_atoms
    center%x = 0; center%y = 0; center%z = 0;
    
    do i_atom = 1, n_atoms, 1
       center%x = center%x + this%coords(i_atom)%x
       center%y = center%y + this%coords(i_atom)%y
       center%z = center%z + this%coords(i_atom)%z
    end do

    do i_atom = 1, n_atoms, 1
       this%coords(i_atom)%x = this%coords(i_atom)%x - inv_n_atoms * center%x
       this%coords(i_atom)%y = this%coords(i_atom)%y - inv_n_atoms * center%y
       this%coords(i_atom)%z = this%coords(i_atom)%z - inv_n_atoms * center%z
    end do

  end subroutine center_structure

