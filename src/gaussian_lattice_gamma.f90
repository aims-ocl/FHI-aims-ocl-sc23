!****s* FHI-aims/gaussian_lattice_gamma
!  NAME
!    gaussian_lattice_gamma
!  SYNOPSIS

subroutine gaussian_lattice_gamma(grid, pos, exponent, max_L, max_n_LM, &
&                               exclude_zero, n_kvec, kvec, lattice_sum)

  !  PURPOSE
  !
  !   Calculate the lattice sum over a Gaussian:
  !
  !      LS(k) := \sum_R |pos+R|^{N+L} \exp(-exponent*|pos+R|^2) Y_{LM}(pos+R)
  !                      * \exp(ik(pos+R))
  !
  !            \approx V_ec^-1 * \int d^3R ...
  !
  !   Modified from gaussian_lattice_sum, but here we only need the G=0 component
  !   and no summation over the reciprocal G vectors is performed. 
  !  USES

  use bravais
  use mpi_tasks, only: check_allocation
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: grid(3, 3)
  real*8, intent(IN) :: pos(3)
  real*8, intent(IN) :: exponent
  integer, intent(IN) :: max_L
  integer, intent(IN) :: max_n_LM
  logical, intent(IN) :: exclude_zero
  !f2py integer, intent(IN) :: exclude_zero
  integer, intent(IN) :: n_kvec
  real*8, intent(IN) :: kvec(3, n_kvec)
  complex*16, intent(OUT) :: lattice_sum(max_n_LM, n_kvec)

  !  INPUTS
  !    o grid -- Lattice vectors for the lattice sum
  !    o pos -- Position in the 0-cell
  !    o exponent -- Exponent in exp(- exponent * |pos+R|**2)
  !    o max_L -- Maximum angular momentum the sum is needed for.
  !    o max_n_LM -- Array dimension of lattice_sum(:,:), max_n_LM >= (max_L+1).
  !    o exclude_zero -- The R==0. can be excluded by this flag.
  !    o n_kvec -- Number of wave vectors.
  !    o kvec -- Wave vectors.
  !  OUTPUTS
  !    o lattice_sum -- Output array.
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

  real*8, parameter :: eta = 45.d0   ! e^-45 == 2.9e-20
  integer :: num_R(3), R1, R2, R3
  integer :: n_LM, i_LM, i_kvec, L, M
  real*8, external :: double_factorial
  real*8, allocatable :: ylms(:)
!  complex*16, allocatable :: subgls(:,:,:)
  real*8 :: Rvec(3), vec(3), Rmax, vec_sq, vec_abs
  real*8 :: val_s, val_l, val_lm
  real*8 :: phi
  complex*16 :: tmp(3)
  complex*16 :: phase(n_kvec)
  integer :: info
  character*150 :: info_str
  character(*), parameter :: func = 'gaussian_lattice_gamma'

  n_LM = (max_L + 1)**2
  allocate(ylms(n_LM), stat=info)
  call check_allocation(info, 'ylms', func)

  Rmax = sqrt(eta / exponent) + sqrt(sum(pos**2))
  call get_n_supercells(3, grid, Rmax, num_R)

!  allocate(subgls(max_n_LM, n_kvec, 2), stat=info)
!  call check_allocation(info, 'subgls', func)

  ! JW: This can be done more efficient:
  ! * Calculate the Gaussian and the phases by simple complex multiplication.
  ! * Use Cartesian radial parts instead of spherical harmonics
  !   and mimic the integration procedure in
  !   [Wieferink, Krueger, Pollmann, PRB 74, 205311 (2006)].  The expansion
  !   coefficients can be found in analytic_multipole_coefficients.f90 [see
  !   reference therein].  The Cartesian terms can also be used to obtain the
  !   higher r^2n terms as used in get_gaussian_Vq(), so the conversion should
  !   be done in an external wrapper.

  tmp = (0.d0, 0.d0)
  lattice_sum = (0.d0, 0.d0)
  num_R(:)=0
  do R1 = -num_R(1), num_R(1)
!     subgls(:,:, 1) = 0.d0
     do R2 = -num_R(2), num_R(2)
!        subgls(:,:, 2) = 0.d0
        do R3 = -num_R(3), num_R(3)
!           if (exclude_zero .and. R1 == 0 .and. R2 == 0 .and. R3 == 0) cycle
           Rvec = grid(:, 1)*R1 + grid(:, 2)*R2 + grid(:, 3)*R3
           vec = pos + Rvec
           vec = pos 
           vec_sq = dot_product(vec, vec)
           vec_abs = sqrt(vec_sq)

!           val_s = exp(- exponent*vec_sq) * vec_abs**power
           val_s = exp(- exponent*vec_sq) 

           call ylm_real(vec, max_L, ylms)

           do i_kvec = 1, n_kvec
              phi = dot_product(vec, kvec(:, i_kvec))
              phase(i_kvec) = cmplx(cos(phi), sin(phi), kind(0.d0))
           end do

           do L = 0, max_L
!              val_l = val_s * vec_abs**L
              val_l = val_s 
              do M = -L, L
                 i_LM = L**2 + L + M + 1
                 val_lm = val_l * ylms(i_LM)
!                   write(use_unit,*) L, M, val_l, ylms(i_LM), val_lm
                 do i_kvec = 1, n_kvec
                    lattice_sum(i_LM, i_kvec) = lattice_sum(i_LM, i_kvec) &
                    &                         + phase(i_kvec) * val_lm
!                  if(i_kvec==1 .and. L == 1 .and. abs(power-2.d0) .lt. 1.d-10 .and. &
!                     all(abs(pos).lt.1.d-10)) then
!                   tmp(M+L+1) = tmp(M+L+1) + exp(-exponent*vec_sq)*ylms(i_LM)*ylms(i_LM)
!!                   write(use_unit,'(5I4,6f16.8)') R1, R2, R3, L, M, exponent, vec_sq, exp(-exponent*vec_sq), ylms(i_LM), &
!!                   tmp(M+L+1)
!                  endif
                 end do
              end do
           end do
        end do
!        subgls(:,:, 1) = subgls(:,:, 1) + subgls(:,:, 2)
     end do
!     lattice_sum = lattice_sum + subgls(:,:, 2)
  end do
  deallocate(ylms)

end subroutine gaussian_lattice_gamma
