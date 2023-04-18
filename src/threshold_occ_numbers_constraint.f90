!-------------------------------------------------------------------------------
!  takes a given set of occupation numbers and sets all occupation numbers
!  below a certain threshold to zero; returns the resulting total charge
!  inaccuracy.
!

      subroutine threshold_occ_numbers_constraint &
        ( constraint_electrons, occ_numbers, diff_electrons_thr &
        )

      use dimensions
      use runtime_choices

      implicit none

!  imported variables

!  input
      real*8, intent(in) :: constraint_electrons

!  output
      real*8, dimension(n_states), intent(inout) :: occ_numbers
      real*8, intent(out) :: diff_electrons_thr

!  local variables
      real*8 :: temp_n_electrons

!  counter
      integer :: i_state


      temp_n_electrons   = 0.d0

      do i_state = 1, n_states, 1
          if (dabs(occ_numbers(i_state)).lt.occupation_thr) then
             occ_numbers(i_state) = 0.d0
          end if
          temp_n_electrons = temp_n_electrons + &
            occ_numbers(i_state)
      enddo



        diff_electrons_thr = temp_n_electrons - &
         constraint_electrons


      end subroutine threshold_occ_numbers_constraint

!-------------------------------------------------------------------------------
