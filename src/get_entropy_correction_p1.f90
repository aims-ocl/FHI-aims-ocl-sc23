!----------------------------------------------------------------------
! Subroutine get_entropy_correction_p1
!
! Evaluates entropy correction to total energy due to occupation smearing
!
! References:
!
! * Kresse et al., Comp. Mat. Sci. 6, 15 - 50 (1996)
!
! * Verstraete and Gonze, Phys. Rev. B 65, 035111 (2001)
!
! * Marzari, doctoral dissertation (1996)
!
! R. Gehrke (2005)
!
!----------------------------------------------------------------------

subroutine get_entropy_correction_p1(KS_eigenvalue, occ_numbers, &
              chemical_potential, chemical_potential_spin, entropy_correction)

  use constants, only: sqrt_pi, pisqrt_inv
  use dimensions, only: spin_degeneracy, n_states, n_spin, n_k_points
  use localorb_io, only: localorb_info, use_unit, OL_norm
  use mpi_tasks, only: aims_stop
  use pbc_lists, only: k_weights
  use runtime_choices, only: occupation_type, occupation_width, entropy_thr, &
      n_methfessel_paxton

  implicit none

  ! input
  real*8, dimension(n_states, n_spin, n_k_points), intent(in) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points), intent(in) :: occ_numbers
  real*8, intent(in) :: chemical_potential
  real*8, dimension(n_spin), intent(in) :: chemical_potential_spin

  ! output
  real*8, intent(out) :: entropy_correction

  ! local variables
  real*8 :: one_over_ow_sqr
  real*8 :: diff_energy
  real*8 :: prefactor
  real*8 :: hermite_arg
  real*8 :: arg
  real*8 :: gauss_weight
  real*8 :: A
  real*8 :: H_even
  real*8 :: H_odd
  real*8 :: one_over_ow
  real*8 :: one_over_spin
  real*8 :: occ_number_with_no_k_weight
  real*8 :: delta
  real*8 :: const

  ! counter
  integer :: i_state
  integer :: i_spin
  integer :: i_k_point
  integer :: i_mp

  character*200 :: info_str

  call check_occs("get_entropy_correction", occ_numbers, .false.)

  one_over_ow_sqr = 1.d0 / (occupation_width * occupation_width)
  one_over_ow = 1.d0 / occupation_width
  one_over_spin = 1.d0 / spin_degeneracy

  entropy_correction = 0.d0

  select case (occupation_type)
  case (0)
     ! gaussian smearing
     ! entropy correction corresponds to (- 0.5 * sum_n \sigma * S_n)
     prefactor = spin_degeneracy * 0.5d0 * occupation_width / (2.d0 * sqrt_pi)

     do i_k_point = 1, n_k_points, 1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1
              diff_energy = KS_eigenvalue(i_state,i_spin,i_k_point) &
                 - chemical_potential_spin(i_spin)
              entropy_correction = entropy_correction + exp(-diff_energy &
                 * diff_energy * one_over_ow_sqr) * k_weights(i_k_point)
           end do
        end do
     end do

  case (1)
     ! fermi smearing
     ! factor two in the prefactor due to spin (so 2 * 0.5 = 1 ...)
     prefactor = spin_degeneracy * 0.5d0 * occupation_width

     do i_k_point = 1, n_k_points, 1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1
              if (((1 - occ_numbers(i_state,i_spin,i_k_point) * one_over_spin) &
                 > entropy_thr) .and. ((occ_numbers(i_state,i_spin,i_k_point) &
                 * one_over_spin) > entropy_thr)) then
                 entropy_correction = entropy_correction &
                    - k_weights(i_k_point) &
                    * (occ_numbers(i_state,i_spin,i_k_point) * one_over_spin &
                    * log(occ_numbers(i_state,i_spin,i_k_point) &
                    * one_over_spin) + (1 &
                    - occ_numbers(i_state,i_spin,i_k_point) * one_over_spin) &
                    * log(1 - occ_numbers(i_state,i_spin,i_k_point) &
                    * one_over_spin))
              end if
           end do
        end do
     end do

  case (2)
     ! methfessel-paxton
     prefactor = spin_degeneracy * 0.25d0 * occupation_width

     do i_k_point = 1, n_k_points, 1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1
              hermite_arg = (KS_eigenvalue(i_state,i_spin,i_k_point) &
                 - chemical_potential_spin(i_spin)) * one_over_ow
              gauss_weight = exp(-hermite_arg * hermite_arg)
              A = pisqrt_inv
              H_even = 1.d0
              H_odd = 2.d0 * hermite_arg
              ! zero order contribution
              entropy_correction = entropy_correction + pisqrt_inv &
                 * gauss_weight * k_weights(i_k_point)

              do i_mp = 1, n_methfessel_paxton, 1
                 ! higher order contribution
                 A = -1.d0 / dble(4 * i_mp) * A
                 H_even = 2.d0 * hermite_arg * H_odd - 2.d0 * i_mp * H_even
                 H_odd = 2.d0 * hermite_arg * H_even - 2.d0 * (i_mp + 1) * H_odd
                 entropy_correction = entropy_correction + A * H_even &
                    * gauss_weight * k_weights(i_k_point)
              end do
           end do
        end do
     end do

  ! for integer occupation by igor
  case (3)
     ! for integer occupation, none entropy_correction should be considered
     entropy_correction = 0.d0
     prefactor = 0.d0

  ! cubic polynomical
  case (4)
     ! S_i  = \int_{-\infty}^{x_i} -t g(t) dt
     ! g(t) = \frac{3}{4} (t^2 - 1)
     delta = 0.75d0 * sqrt_pi * occupation_width
     prefactor = spin_degeneracy * 0.5d0 * occupation_width * 0.1875d0 &
        / (delta**4)
     const = delta**2

     do i_k_point = 1, n_k_points, 1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1
              if (KS_eigenvalue(i_state,i_spin,i_k_point) &
                 > chemical_potential_spin(i_spin) - delta &
                 .and. KS_eigenvalue(i_state,i_spin,i_k_point) &
                 < chemical_potential_spin(i_spin) + delta) then
                 diff_energy = KS_eigenvalue(i_state,i_spin,i_k_point) &
                    - chemical_potential_spin(i_spin)
                 entropy_correction = entropy_correction + (((diff_energy**2) &
                    - const)**2) * k_weights(i_k_point)
              end if
           end do
        end do
     end do

  ! Marzari-Vanderbilt (cold)
  case (5)
     prefactor = 0.5d0 * sqrt(0.5d0) * spin_degeneracy * occupation_width &
        * pisqrt_inv

     do i_k_point = 1, n_k_points, 1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1
              diff_energy = KS_eigenvalue(i_state,i_spin,i_k_point) &
                 - chemical_potential_spin(i_spin)
              arg = -diff_energy * one_over_ow - sqrt(0.5d0)
              entropy_correction = entropy_correction - arg * exp(-arg**2) &
                 * k_weights(i_k_point)
           end do
        end do
     end do

  case default
     write(info_str,"(2X,2A)") "Unknown occupation_type in", &
        " get_entropy_correction_p1."
     call localorb_info(info_str,use_unit,"(A)",OL_norm)
     write(info_str,"(2X,A)") "* Aborting."
     call localorb_info(info_str,use_unit,"(A)",OL_norm)
     call aims_stop
  end select

  entropy_correction = - (prefactor * entropy_correction)

end subroutine get_entropy_correction_p1
