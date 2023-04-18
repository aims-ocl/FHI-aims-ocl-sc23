!------------------------------------------------------------------------------
!****s* FHI-aims/check_occs
!  NAME
!    check_occs
!  SYNOPSIS

subroutine check_occs(caller, occ_numbers, expect_kweights)

  !  PURPOSE
  !    Check if kweights are incorporated into occ_numbers.
  !    Output error message if expectation is wrong.
  !  USES

  use physics, only: n_electrons
  use dimensions
  use localorb_io
  use pbc_lists
  implicit none

  !  ARGUMENTS

  character(*), intent(IN) :: caller
  real*8, intent(IN) :: occ_numbers(n_states, n_spin, n_k_points)
  logical, intent(IN) :: expect_kweights

  !  INPUTS
  !    o caller -- Calling subroutine for output
  !    o occ_numbers -- Current occupation numbers
  !    o expect_kweights -- Expect k_weights to be incorporated? 
  !  OUTPUTS
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  real*8 :: BZ_integral_over_one
  integer :: i_k_point
  character*100 :: info_str

  if (expect_kweights) then
     BZ_integral_over_one = sum(occ_numbers)
  else
     BZ_integral_over_one = 0.d0
     do i_k_point = 1, n_k_points
        BZ_integral_over_one &
        & = BZ_integral_over_one &
        & + sum(occ_numbers(:,:, i_k_point)) * k_weights(i_k_point)
     end do
  end if

  if (abs(BZ_integral_over_one - n_electrons) > 1d-5*n_electrons) then
     if (expect_kweights) then
        write(info_str, &
        &     "('*** ',A,' expected kweighted occ_numbers')") &
        & trim(caller)
     else
        write(info_str, &
        &     "('*** ',A,' expected pure occ_numbers')") &
        & trim(caller)
     end if
     call localorb_info(info_str)
     write(info_str, &
     &     "('*** sum(occ_numbers) =',ES24.16,' instead of',ES24.16)") &
     & sum(occ_numbers), n_electrons
     call localorb_info(info_str)
  end if

end subroutine check_occs
!******
