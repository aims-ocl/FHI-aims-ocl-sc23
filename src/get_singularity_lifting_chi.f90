!****s* FHI-aims/get_singularity_lifting_chi
!  NAME
!    get_singularity_lifting_chi
!  SYNOPSIS

subroutine get_singularity_lifting_chi( &
&               recip_lattice_vector, n_k_points_xyz, start_Gamma, thres, chi)

  !  PURPOSE
  !
  !    Wrapper around get_singularity_lifting_chi_gamma() to perform linear
  !    extrapolation.  Stsrts with start_Gamma and scales it with 0.5 until
  !    the extrapolated results differ by less than thres.
  !
  !  USES

  use mpi_tasks
  use localorb_io
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: recip_lattice_vector(3, 3)
  integer, intent(IN) :: n_k_points_xyz(3)
  real*8, intent(IN) :: start_Gamma, thres
  real*8, intent(OUT) :: chi

  !  INPUTS
  !    o recip_lattice_vector -- Reciprocal lattice
  !    o n_k_points_xyz -- Number of Monkhorst-Pack k-points in given
  !                        reciprocal direction.
  !    o Gamma -- Starting Gaussian decay coefficients
  !    o thres -- Convergence criterion
  !  OUTPUTS
  !    o chi -- Gamma point "value" of 1/q^2
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

  real*8 :: Gamma_old, chi_old, Gamma_new, chi_new
  real*8 :: chi_extra, chi_extra_old
  real*8 :: dchi
  integer :: i_iter
  integer, parameter :: n_max_iter = 10
  character*150 :: info_str
  character(*), parameter :: func = 'get_singularity_lifting_chi'

  call localorb_info("  Calculating singularity lifting chi")


  chi_extra_old = 0.5d0 * huge(chi_extra_old)
  Gamma_old = start_Gamma
  call get_singularity_lifting_chi_gamma( &
  &        recip_lattice_vector, n_k_points_xyz, Gamma_old, chi_old)
  write(info_str,"(2X,'| Gamma:',ES10.2,'; chi:',F24.16)") &
  & Gamma_old, chi_old
  call localorb_info(info_str)
  do i_iter = 1, n_max_iter
     Gamma_new = Gamma_old / 2.d0
     call get_singularity_lifting_chi_gamma( &
     &        recip_lattice_vector, n_k_points_xyz, Gamma_new, chi_new)

     ! Get finite difference derivative
     dchi =(chi_new - chi_old) / (Gamma_new - Gamma_old)
     ! Extrapolate to zero.
     chi_extra = chi_new + dchi * (0.d0 - Gamma_new)
     write(info_str,"(2X,'| Gamma:',ES10.2,'; chi:',F24.16,'; extra:',F24.16)")&
     & Gamma_new, chi_new, chi_extra
     call localorb_info(info_str)

     ! Stop if converged
     if (abs(chi_extra - chi_extra_old) < thres) exit

     Gamma_old = Gamma_new
     chi_old = chi_new
     chi_extra_old = chi_extra
  end do
  write(info_str,"(2X,'| Using chi:',F24.16,'; last diff:',ES10.2)") &
  & chi_extra, chi_extra - chi_extra_old
  call localorb_info(info_str)
  if (i_iter > n_max_iter) call aims_stop('Not converged to thres', func)

  chi = chi_extra

end subroutine get_singularity_lifting_chi
!******
!------------------------------------------------------------------------------
!****s* FHI-aims/get_singularity_lifting_chi_gamma
!  NAME
!    get_singularity_lifting_chi_gamma
!  SYNOPSIS

subroutine get_singularity_lifting_chi_gamma( &
&                             recip_lattice_vector, n_k_points_xyz, Gamma, chi)

  !  PURPOSE
  !
  !    Calculate the chi-value for lifting the Coulomb singularity at q->0 in
  !    [1].  The prefactor is chosen such that chi can perform as a drop-in
  !    replacement for 1/q^2.  To get chi in the nomenclature of [1], set
  !    n_k_points_xyz to (/1, 1, 1/) and multiply chi by 4*pi/cell_volume.
  !
  !    From [1] it is best to take the limit Gamma->0.  This is best achieved
  !    to take two different small but finite values and do linear
  !    extrapolation.
  !
  !  USES

  use constants
  use bravais
  use mpi_tasks, only: aims_stop
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: Gamma
  real*8, intent(IN) :: recip_lattice_vector(3, 3)
  integer, intent(IN) :: n_k_points_xyz(3)
  real*8, intent(OUT) :: chi

  !  INPUTS
  !    o Gamma -- Gaussian decay coefficient (the smaller, the more "accurate")
  !    o recip_lattice_vector -- Reciprocal lattice
  !    o n_k_points_xyz -- Number of Monkhorst-Pack k-points in given
  !                        reciprocal direction.
  !  OUTPUTS
  !    o chi -- Gamma point "value" of 1/q^2
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    [1] Peter Broqvist, Audrius Alkauskas, and Alfredo Pasquarello,
  !    "Hybrid-functional calculations with plane-wave basis sets: Effect of
  !    singularity correction on total energies, energy eigenvalues, and
  !    defect energy levels", Physical Review B 80, 085114 (2009).
  !    [2] Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
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

  real*8 :: qgrid(3, 3), recip_cell_volume
  integer, parameter :: n_kvec = 1
  real*8, parameter :: kvec(3) = (/0.d0, 0.d0, 0.d0/)
  real*8, parameter :: pos(3) = (/0.d0, 0.d0, 0.d0/)
  real*8, parameter :: power = -2.d0
  logical, parameter :: exclude_zero = .true.
  integer, parameter :: max_L = 0
  integer, parameter :: max_n_LM = 1
  complex*16 :: lattice_sum(max_n_LM, n_kvec)
  real*8 :: discrete_sum_without, continuous_integral
  integer :: i, n_k_points
  character(*), parameter :: func = 'get_singularity_lifting_chi_gamma'

  ! --- Finite sum

  do i = 1, 3
     qgrid(:, i) = recip_lattice_vector(:, i) / n_k_points_xyz(i)
  end do
  ! lattice_sum := \sum_{q \in qgrid^{q/=0} q^{-2} \exp(-Gamma*q^2) Y_{00}(q)
  call gaussian_lattice_sum(qgrid, pos, power, Gamma, max_L, max_n_LM, &
  &                         exclude_zero, n_kvec, kvec, lattice_sum)
  if (aimag(lattice_sum(1, 1)) > 1d-10) then
     call aims_stop('Nonvanishing imaginary part', func)
  end if
  discrete_sum_without = dble(lattice_sum(1,1)) * sqrt(4*pi)
  ! = \sum_{q \in qgrid^{q/=0} q^{-2} \exp(-Gamma*q^2)

  ! --- Continuous integral

  continuous_integral = 4*pi / sqrt(Gamma) * sqrt(pi)/2
  ! = \int_K^3 d^3q exp(-Gamma*q**2)/q^2
  ! = 4pi/sqrt(Gamma) * \int_0^\infty dx e^{-x^2}.

  ! --- Result

  call get_cell_volume(recip_lattice_vector, recip_cell_volume)
  n_k_points = product(n_k_points_xyz)

  chi = continuous_integral * n_k_points / recip_cell_volume &
  &   - discrete_sum_without

end subroutine get_singularity_lifting_chi_gamma
!******
