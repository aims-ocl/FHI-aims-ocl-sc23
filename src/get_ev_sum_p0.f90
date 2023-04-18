!****s* FHI-aims/get_ev_sum_p0
!  NAME
!    get_ev_sum_p0
!  SYNOPSIS

subroutine get_ev_sum_p0 &
     ( occ_numbers, KS_eigenvalue, av_core_shift, ev_sum, &
     ev_sum_shifted &
     )


  !  PURPOSE
  !
  !  Subroutine get_ev_sum calculates the 
  !  sum of eigenvalues for the non-spinpolarised, non-charged case.!
  !
  !  USES

  use dimensions
  use geometry
  use constants
  use pbc_lists
  implicit none

  !  ARGUMENTS

  real*8, dimension(n_states,n_spin, n_k_points) :: occ_numbers 
  real*8, dimension(n_states,n_spin, n_k_points) :: KS_eigenvalue
  real*8 :: av_core_shift
  real*8 :: ev_sum
  real*8 :: ev_sum_shifted

  !  INPUTS
  !    o occ_numbers -- occupations of the states
  !    o KS_eigenvalue -- Kohn-Sham eigenvalues
  !    o av_core_shift -- core shift
  !
  !  OUTPUT
  !    o ev_sum -- sum of eigenvalues
  !    o ev_sum_shifted -- shifted sum of eigenvalues
  !
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



  !  local variables

  !     counters

  integer :: i_state
  integer :: i_spin
  integer :: i_k_point

  call check_occs('get_ev_sum_p0', occ_numbers, .false.)

  !  begin work


  ev_sum = 0.d0
  ev_sum_shifted = 0.d0


  do i_k_point = 1, n_k_points, 1
     do i_spin = 1, n_spin, 1
        do i_state = 1, n_states, 1


           ev_sum = ev_sum + k_weights(i_k_point) * &
                occ_numbers(i_state,i_spin,i_k_point)*KS_eigenvalue(i_state,i_spin, i_k_point)

           !               ev_sum_shifted = &
           !                    ev_sum_shifted + occ_numbers(i_state,i_spin,i_k_point) * &
           !                    (KS_eigenvalue(i_state,i_spin,i_k_point) - av_core_shift) * k_weights(i_k_point)

        enddo
     enddo
  end do



end subroutine get_ev_sum_p0
!---------------------------------------------------------------------
!******
