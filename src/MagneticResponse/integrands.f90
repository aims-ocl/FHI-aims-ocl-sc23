!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Provides the function bodies that are used for integration in
!!  'integrate'. Each subroutine calculates H|phi> over a batch of
!!  grid points. 'integrate' is a generic subroutine in
!!  integration.f90 that can be used to integrate most functions of
!!  interest. Each integral is described by the derived type
!!  integrand_t. It provides each integrand an ID and additional
!!  information such as whether the gradient of the basis function is
!!  required for evaluation or whether to store the full matrix in a
!!  packed array after integration (default is to store only the lower
!!  triangle). By default, 'integrate' only evaluates the basic
!!  quantities such the basis function values at grid points for all
!!  integrand types. If one wants to compute, e.g., the orbital zeeman
!!  integral, the gradient of the basis functions must be requested by
!!  setting the basis_gradient (see below) flag to .true., else the
!!  computation will result in a segfault. All such flags are
!!  processed by 'integrate' before the actual integration.
!!
!!  COMMENTS
!!
!!  Various prefactors (e.g., alpha^2) have been omitted here and are
!!  found elsewhere. Also, additional operations may be performed
!!  afterwards (e.g, symmetrization or antisymmetrization).
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
module integrands

  use dimensions, only: n_spin, use_gga
  use MR_global,  only: calc_full_tensor
  use tools,      only: norm2, mul
  use types,      only: dp

  implicit none

  public
  private :: norm2

  type :: integrand_t
     integer :: id = 0
     ! Number of spatial (or generalized) directions that are
     ! simultaneously processed.
     integer :: n_dimensions = 1
     ! 0 - No explicit symmetrization required
     ! 1 - Symmetrize
     ! 2 - Antisymmetrize
     ! Do not use this with store_full_matrix=.true.
     integer :: symmetrization = 0
     ! Whether to require the evaluation of Nabla|phi>
     logical :: basis_gradient = .false.
     ! Whether to transfer the full contents or only the lower
     ! triangle to a packed matrix.
     logical :: store_full_matrix = .false.
     ! Whether to evaluate the unperturbed Hamiltonian
     logical :: hamiltonian0 = .false.
     ! Whether the electron density and density gradient (if GGA) are
     ! used (this are implied by hamiltonian0).
     logical :: rho = .false.
     ! Often the size of work arrays is proportional to
     ! n_dimensions. In case of GIAOs, the number of dimensions (e.g.,
     ! spatial directions) and the size of work arrays can be
     ! different. If so, use GIAO_work_size to explicitly request
     ! additional memory for certain work arrays. See integration.f90.
     integer :: GIAO_work_size = 0
     ! Whether to find positions relative to basis function n
     logical :: r_n_basis = .false.
     ! Whether to find differences of position vectors of basis
     ! function m and n, R_m-R_n
     logical :: r_mn_basis = .false.
     ! For a spin-polarized case, whether the up/down matrix elements
     ! need to be calculated separately.
     logical :: spin_sensitive = .false.
     ! Whether to calculate only the diagonal elements of the full tensor
     logical :: calc_full_tensor = .false.
  end type integrand_t

  ! MAGNETIC RESPONSE

  integer, parameter :: FERMI_CONTACT_SR = 0
  type(integrand_t) :: INT_FERMI_CONTACT_SR

  integer, parameter :: PARAMAGNETIC_SO = 1
  type(integrand_t) :: INT_PARAMAGNETIC_SO

  integer, parameter :: PARAMAGNETIC_SO_SR = 2
  type(integrand_t) :: INT_PARAMAGNETIC_SO_SR

  integer, parameter :: SPIN_DIPOLE_XX_YY_ZZ = 3
  type(integrand_t) :: INT_SPIN_DIPOLE_XX_YY_ZZ

  integer, parameter :: SPIN_DIPOLE_XY_XZ_YZ = 4
  type(integrand_t) :: INT_SPIN_DIPOLE_XY_XZ_YZ

  integer, parameter :: DIAMAGNETIC_SO = 5
  type(integrand_t) :: INT_DIAMAGNETIC_SO

  integer, parameter :: DIAMAGNETIC_SO_SR = 6
  type(integrand_t) :: INT_DIAMAGNETIC_SO_SR

  integer, parameter :: ORBITAL_ZEEMAN = 7
  type(integrand_t) :: INT_ORBITAL_ZEEMAN

  integer, parameter :: GIAO_PARAMAGNETIC = 8
  type(integrand_t) :: INT_GIAO_PARAMAGNETIC

  integer, parameter :: DIAMAGNETIC_SHIELDING_STD = 9
  type(integrand_t) :: INT_DIAMAGNETIC_SHIELDING_STD

  integer, parameter :: GIAO_OVERLAP1 = 10
  type(integrand_t) :: INT_GIAO_OVERLAP1

  integer, parameter :: GIAO_OVERLAP2 = 11
  type(integrand_t) :: INT_GIAO_OVERLAP2

  integer, parameter :: DIAMAGNETIC_SHIELDING = 12
  type(integrand_t) :: INT_DIAMAGNETIC_SHIELDING

  integer, parameter :: ELECTRIC_FIELD_GRADIENT = 13
  type(integrand_t) :: INT_ELECTRIC_FIELD_GRADIENT

  integer, parameter :: DIAMAGNETIC_MAGNETIZABILITY = 14
  type(integrand_t) :: INT_DIAMAGNETIC_MAGNETIZABILITY

  integer, parameter :: DIAMAGNETIC_MAGNETIZABILITY_STD = 15
  type(integrand_t) :: INT_DIAMAGNETIC_MAGNETIZABILITY_STD

  ! DFPT

  integer, parameter :: XC_KERNEL = 16
  type(integrand_t) :: INT_XC_KERNEL

  integer, parameter :: FIRST_ORDER_HARTREE = 17
  type(integrand_t) :: INT_FIRST_ORDER_HARTREE

  ! OTHER

  integer, parameter :: ELECTRON_DENSITY = 18
  type(integrand_t) :: INT_ELECTRON_DENSITY

  integer, parameter :: POSITION = 19
  type(integrand_t) :: INT_POSITION

contains

  subroutine initialize_integrands()
    INT_FERMI_CONTACT_SR%id = FERMI_CONTACT_SR
    INT_FERMI_CONTACT_SR%basis_gradient = .true.

    INT_PARAMAGNETIC_SO%id = PARAMAGNETIC_SO
    INT_PARAMAGNETIC_SO%n_dimensions   = 3
    INT_PARAMAGNETIC_SO%symmetrization = 2
    INT_PARAMAGNETIC_SO%basis_gradient = .true.

    INT_PARAMAGNETIC_SO_SR%id = PARAMAGNETIC_SO_SR
    INT_PARAMAGNETIC_SO_SR%n_dimensions   = 3
    INT_PARAMAGNETIC_SO_SR%symmetrization = 2
    INT_PARAMAGNETIC_SO_SR%basis_gradient = .true.

    INT_SPIN_DIPOLE_XX_YY_ZZ%id = SPIN_DIPOLE_XX_YY_ZZ
    INT_SPIN_DIPOLE_XX_YY_ZZ%n_dimensions = 3

    INT_SPIN_DIPOLE_XY_XZ_YZ = INT_SPIN_DIPOLE_XX_YY_ZZ
    INT_SPIN_DIPOLE_XY_XZ_YZ%id = SPIN_DIPOLE_XY_XZ_YZ

    INT_DIAMAGNETIC_SO%id = DIAMAGNETIC_SO
    INT_DIAMAGNETIC_SO%calc_full_tensor = calc_full_tensor
    if (calc_full_tensor) then
       INT_DIAMAGNETIC_SO%n_dimensions = 9
    else
       INT_DIAMAGNETIC_SO%n_dimensions = 3
    end if

    INT_DIAMAGNETIC_SO_SR%id = DIAMAGNETIC_SO_SR
    INT_DIAMAGNETIC_SO_SR%calc_full_tensor = calc_full_tensor
    if (calc_full_tensor) then
       INT_DIAMAGNETIC_SO_SR%n_dimensions = 9
    else
       INT_DIAMAGNETIC_SO_SR%n_dimensions = 3
    end if

    INT_ORBITAL_ZEEMAN%id = ORBITAL_ZEEMAN
    INT_ORBITAL_ZEEMAN%n_dimensions   = 3
    INT_ORBITAL_ZEEMAN%symmetrization = 2
    INT_ORBITAL_ZEEMAN%basis_gradient = .true.

    INT_GIAO_PARAMAGNETIC%id = GIAO_PARAMAGNETIC
    INT_GIAO_PARAMAGNETIC%n_dimensions      = 3
    INT_GIAO_PARAMAGNETIC%basis_gradient    = .true.
    INT_GIAO_PARAMAGNETIC%store_full_matrix = .true.
    INT_GIAO_PARAMAGNETIC%hamiltonian0      = .true.
    INT_GIAO_PARAMAGNETIC%r_n_basis         = .true.
    INT_GIAO_PARAMAGNETIC%r_mn_basis        = .true.
    INT_GIAO_PARAMAGNETIC%spin_sensitive    = .true.

    INT_DIAMAGNETIC_SHIELDING_STD%id = DIAMAGNETIC_SHIELDING_STD
    INT_DIAMAGNETIC_SHIELDING_STD%calc_full_tensor = calc_full_tensor
    if (calc_full_tensor) then
       INT_DIAMAGNETIC_SHIELDING_STD%n_dimensions = 9
    else
       INT_DIAMAGNETIC_SHIELDING_STD%n_dimensions = 3
    end if

    INT_GIAO_OVERLAP1%id = GIAO_OVERLAP1
    INT_GIAO_OVERLAP1%n_dimensions      = 3
    INT_GIAO_OVERLAP1%store_full_matrix = .true.
    INT_GIAO_OVERLAP1%r_mn_basis        = .true.

    INT_GIAO_OVERLAP2%id = GIAO_OVERLAP2
    INT_GIAO_OVERLAP2%GIAO_work_size   = 6
    INT_GIAO_OVERLAP2%r_mn_basis       = .true.
    INT_GIAO_OVERLAP2%calc_full_tensor = calc_full_tensor
    if (calc_full_tensor) then
       INT_GIAO_OVERLAP2%n_dimensions = 9
    else
       INT_GIAO_OVERLAP2%n_dimensions = 3
    end if

    INT_DIAMAGNETIC_SHIELDING%id = DIAMAGNETIC_SHIELDING
    INT_DIAMAGNETIC_SHIELDING%symmetrization   = 1
    INT_DIAMAGNETIC_SHIELDING%basis_gradient   = .true.
    INT_DIAMAGNETIC_SHIELDING%r_n_basis        = .true.
    INT_DIAMAGNETIC_SHIELDING%r_mn_basis       = .true.
    INT_DIAMAGNETIC_SHIELDING%calc_full_tensor = calc_full_tensor
    if (calc_full_tensor) then
       INT_DIAMAGNETIC_SHIELDING%n_dimensions   = 9
       INT_DIAMAGNETIC_SHIELDING%GIAO_work_size = 9
    else
       INT_DIAMAGNETIC_SHIELDING%n_dimensions   = 3
       INT_DIAMAGNETIC_SHIELDING%GIAO_work_size = 6
    end if

    INT_ELECTRIC_FIELD_GRADIENT%id = ELECTRIC_FIELD_GRADIENT
    INT_ELECTRIC_FIELD_GRADIENT%n_dimensions = 3

    INT_DIAMAGNETIC_MAGNETIZABILITY%id = DIAMAGNETIC_MAGNETIZABILITY
    INT_DIAMAGNETIC_MAGNETIZABILITY%symmetrization   = 1
    INT_DIAMAGNETIC_MAGNETIZABILITY%GIAO_work_size   = 9
    INT_DIAMAGNETIC_MAGNETIZABILITY%basis_gradient   = .true.
    INT_DIAMAGNETIC_MAGNETIZABILITY%hamiltonian0     = .true.
    INT_DIAMAGNETIC_MAGNETIZABILITY%r_n_basis        = .true.
    INT_DIAMAGNETIC_MAGNETIZABILITY%r_mn_basis       = .true.
    INT_DIAMAGNETIC_MAGNETIZABILITY%spin_sensitive   = .true.
    INT_DIAMAGNETIC_MAGNETIZABILITY%calc_full_tensor = calc_full_tensor
    if (calc_full_tensor) then
       INT_DIAMAGNETIC_MAGNETIZABILITY%n_dimensions = 9
    else
       INT_DIAMAGNETIC_MAGNETIZABILITY%n_dimensions = 3
    end if

    INT_DIAMAGNETIC_MAGNETIZABILITY_STD%id = DIAMAGNETIC_MAGNETIZABILITY_STD
    INT_DIAMAGNETIC_MAGNETIZABILITY_STD%calc_full_tensor = calc_full_tensor
    if (calc_full_tensor) then
       INT_DIAMAGNETIC_MAGNETIZABILITY_STD%n_dimensions = 9
    else
       INT_DIAMAGNETIC_MAGNETIZABILITY_STD%n_dimensions = 3
    end if

    INT_XC_KERNEL%id = XC_KERNEL
    ! n_dimensions must be explicitly set elsewhere!
    INT_XC_KERNEL%rho            = .true.
    INT_XC_KERNEL%spin_sensitive = .true.
    if (use_gga) then
       INT_XC_KERNEL%symmetrization = 1
       INT_XC_KERNEL%basis_gradient = .true.
    end if

    INT_FIRST_ORDER_HARTREE%id = FIRST_ORDER_HARTREE
    ! n_dimensions must be explicitly set elsewhere!

    INT_ELECTRON_DENSITY%id = ELECTRON_DENSITY
    ! n_dimensions must be explicitly set elsewhere!
    INT_ELECTRON_DENSITY%spin_sensitive = .true.
    if (use_gga) INT_ELECTRON_DENSITY%basis_gradient = .true.

    INT_POSITION%id = POSITION
    INT_POSITION%n_dimensions   = 3
    INT_POSITION%symmetrization = 1
  end subroutine initialize_integrands

  !! Fermi contact, scalar relativistic (for nonrelativistic, see
  !! integrate_FC.f90)
  !!
  !!       alpha^2   K
  !! H = - ------- ----- r_A Nabla
  !!          3    r_A^3
  !!
  pure subroutine f_fermi_contact_sr(n_compute, n_points, partition, &
       & nabla_wave, r_A, zora_factor, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: nabla_wave(n_compute,3,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(in) :: zora_factor(n_compute,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points)
    integer :: ip
    do ip = 1, n_points
       matrix_wave(:,ip) = partition(ip)*zora_factor(:,ip)* &
            & matmul(nabla_wave(:,:,ip),r_A(:,ip))/norm2(r_A(:,ip))**3
    end do
  end subroutine f_fermi_contact_sr
  !!
  !! Paramagnetic spin-orbit, nonrelativistic
  !!
  !!       alpha^2 r_A x Nabla
  !! H = -i------- -----------
  !!          2       r_A^3
  !!
  pure subroutine f_paramagnetic_SO(n_compute, n_points, partition, &
       & nabla_wave, r_A, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: nabla_wave(n_compute,3,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: ip
    do ip = 1, n_points
       matrix_wave(:,ip,1) = &
            & r_A(2,ip)*nabla_wave(:,3,ip) - r_A(3,ip)*nabla_wave(:,2,ip)
       matrix_wave(:,ip,2) = &
            & r_A(3,ip)*nabla_wave(:,1,ip) - r_A(1,ip)*nabla_wave(:,3,ip)
       matrix_wave(:,ip,3) = &
            & r_A(1,ip)*nabla_wave(:,2,ip) - r_A(2,ip)*nabla_wave(:,1,ip)
       matrix_wave(:,ip,:) = matrix_wave(:,ip,:)*partition(ip)/ &
            & norm2(r_A(:,ip))**3
    end do
  end subroutine f_paramagnetic_SO
  !!
  !! Paramagnetic spin-orbit, scalar relativistic
  !!
  !!       alpha^2   r_A x Nabla
  !! H = -i------- K -----------
  !!          2         r_A^3
  !!
  pure subroutine f_paramagnetic_SO_sr(n_compute, n_points, partition, &
       & nabla_wave, r_A, zora_factor, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: nabla_wave(n_compute,3,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(in) :: zora_factor(n_compute,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: ip
    do ip = 1, n_points
       matrix_wave(:,ip,1) = zora_factor(:,ip)* &
            & r_A(2,ip)*nabla_wave(:,3,ip) - r_A(3,ip)*nabla_wave(:,2,ip)
       matrix_wave(:,ip,2) = zora_factor(:,ip)* &
            & r_A(3,ip)*nabla_wave(:,1,ip) - r_A(1,ip)*nabla_wave(:,3,ip)
       matrix_wave(:,ip,3) = zora_factor(:,ip)* &
            & r_A(1,ip)*nabla_wave(:,2,ip) - r_A(2,ip)*nabla_wave(:,1,ip)
       matrix_wave(:,ip,:) = &
            & matrix_wave(:,ip,:)*partition(ip)/norm2(r_A(:,ip))**3
    end do
  end subroutine f_paramagnetic_SO_sr
  !!
  !! Spin-dipole
  !!
  !! H_i = Sum_j=1^3 Sigma_z H_ij
  !!
  !!                3r_A[i]r_A[j]r_A^-2 - delta_ij
  !! H_ij = alpha^2 ------------------------------
  !!                             r_A^3
  !!
  !! It is the H_ij terms that are evaluated here. In the first
  !! subroutine, delta_ij = 1 and in the second one delta_ij = 0.
  pure subroutine f_spin_dipole_xx_yy_zz(n_compute, n_points, partition, &
       & wave, r_A, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: ip
    do ip = 1, n_points
       matrix_wave(:,ip,1) = (3*r_A(1,ip)**2/norm2(r_A(:,ip))**2 - 1)*wave(:,ip)
       matrix_wave(:,ip,2) = (3*r_A(2,ip)**2/norm2(r_A(:,ip))**2 - 1)*wave(:,ip)
       matrix_wave(:,ip,3) = (3*r_A(3,ip)**2/norm2(r_A(:,ip))**2 - 1)*wave(:,ip)
       matrix_wave(:,ip,:) = &
            & matrix_wave(:,ip,:)*partition(ip)/norm2(r_A(:,ip))**3
    end do
  end subroutine f_spin_dipole_xx_yy_zz

  pure subroutine f_spin_dipole_xy_xz_yz(n_compute, n_points, partition, &
       & wave, r_A, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: ip
    do ip = 1, n_points
       matrix_wave(:,ip,1) = r_A(1,ip)*r_A(2,ip)*wave(:,ip)
       matrix_wave(:,ip,2) = r_A(1,ip)*r_A(3,ip)*wave(:,ip)
       matrix_wave(:,ip,3) = r_A(2,ip)*r_A(3,ip)*wave(:,ip)
       matrix_wave(:,ip,:) = &
            & matrix_wave(:,ip,:)*partition(ip)/norm2(r_A(:,ip))**5
    end do
    matrix_wave = 3*matrix_wave
  end subroutine f_spin_dipole_xy_xz_yz
  !!
  !! Diamagnetic spin-orbit, nonrelativistic
  !!
  !!     alpha^4 r_A r_B - r_B r_A^T
  !! H = ------- -------------------
  !!        2        r_A^3 r_B^3
  !!
  pure subroutine f_diamagnetic_SO_diag(n_compute, n_points, partition, wave, &
       & r_A, r_B, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(in) :: r_B(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: ip, i_dir
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir) = (sum(r_A(:,ip)*r_B(:,ip)) - &
               & r_B(i_dir,ip)*r_A(i_dir,ip))*wave(:,ip)
       end do
       matrix_wave(:,ip,:) = matrix_wave(:,ip,:)*partition(ip)/ &
            & norm2(r_A(:,ip))**3/norm2(r_B(:,ip))**3
    end do
  end subroutine f_diamagnetic_SO_diag

  pure subroutine f_diamagnetic_SO_full(n_compute, n_points, partition, wave, &
       & r_A, r_B, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(in) :: r_B(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,3)
    integer :: ip, i_dir
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir,i_dir) = (sum(r_A(:,ip)*r_B(:,ip)) - &
               & r_B(i_dir,ip)*r_A(i_dir,ip))*wave(:,ip)
       end do
       matrix_wave(:,ip,1,2) = -r_B(1,ip)*r_A(2,ip)*wave(:,ip)
       matrix_wave(:,ip,1,3) = -r_B(1,ip)*r_A(3,ip)*wave(:,ip)
       matrix_wave(:,ip,2,1) = -r_B(2,ip)*r_A(1,ip)*wave(:,ip)
       matrix_wave(:,ip,2,3) = -r_B(2,ip)*r_A(3,ip)*wave(:,ip)
       matrix_wave(:,ip,3,1) = -r_B(3,ip)*r_A(1,ip)*wave(:,ip)
       matrix_wave(:,ip,3,2) = -r_B(3,ip)*r_A(2,ip)*wave(:,ip)
       matrix_wave(:,ip,:,:) = matrix_wave(:,ip,:,:)*partition(ip)/ &
            & norm2(r_A(:,ip))**3/norm2(r_B(:,ip))**3
    end do
  end subroutine f_diamagnetic_SO_full
  !!
  !! Diamagnetic spin-orbit, scalar relativistic
  !!
  !!     alpha^4 r_A r_B - r_B r_A^T
  !! H = ------- ------------------- zora_factor
  !!        2        r_A^3 r_B^3
  !!
  pure subroutine f_diamagnetic_SO_sr_diag(n_compute, n_points, partition, &
       & wave, r_A, r_B, zora_factor, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(in) :: r_B(3,n_points)
    real(dp), intent(in) :: zora_factor(n_compute,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: ip, i_dir
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir) = (sum(r_A(:,ip)*r_B(:,ip)) - &
               & r_B(i_dir,ip)*r_A(i_dir,ip))*wave(:,ip)*zora_factor(:,ip)
       end do
       matrix_wave(:,ip,:) = matrix_wave(:,ip,:)*partition(ip)/ &
            & norm2(r_A(:,ip))**3/norm2(r_B(:,ip))**3
    end do
  end subroutine f_diamagnetic_SO_sr_diag

  pure subroutine f_diamagnetic_SO_sr_full(n_compute, n_points, partition, &
       & wave, r_A, r_B, zora_factor, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(in) :: r_B(3,n_points)
    real(dp), intent(in) :: zora_factor(n_compute,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,3)
    integer :: ip, i_dir
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir,i_dir) = (sum(r_A(:,ip)*r_B(:,ip)) - &
               & r_B(i_dir,ip)*r_A(i_dir,ip))*wave(:,ip)*zora_factor(:,ip)
       end do
       matrix_wave(:,ip,1,2) = -r_B(1,ip)*r_A(2,ip)*wave(:,ip)*zora_factor(:,ip)
       matrix_wave(:,ip,1,3) = -r_B(1,ip)*r_A(3,ip)*wave(:,ip)*zora_factor(:,ip)
       matrix_wave(:,ip,2,1) = -r_B(2,ip)*r_A(1,ip)*wave(:,ip)*zora_factor(:,ip)
       matrix_wave(:,ip,2,3) = -r_B(2,ip)*r_A(3,ip)*wave(:,ip)*zora_factor(:,ip)
       matrix_wave(:,ip,3,1) = -r_B(3,ip)*r_A(1,ip)*wave(:,ip)*zora_factor(:,ip)
       matrix_wave(:,ip,3,2) = -r_B(3,ip)*r_A(2,ip)*wave(:,ip)*zora_factor(:,ip)
       matrix_wave(:,ip,:,:) = matrix_wave(:,ip,:,:)*partition(ip)/ &
            & norm2(r_A(:,ip))**3/norm2(r_B(:,ip))**3
    end do
  end subroutine f_diamagnetic_SO_sr_full
  !!
  !! Electric Field Gradient, nonrelativistic
  !!
  !!      3 r_A r_A^T - r_A^T r_A   1
  !! H = ------------------------- ---
  !!                r_A^5           3
  !!
  pure subroutine f_electric_field_gradient_diag(n_compute, n_points, &
       & partition, wave, r_A, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: i_dir, ip
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir) = (r_A(i_dir,ip)*r_A(i_dir,ip) - &
               & sum(r_A(:,ip)*r_A(:,ip))/3)*wave(:,ip)
       end do
       matrix_wave(:,ip,:) = &
            & matrix_wave(:,ip,:)*partition(ip)/norm2(r_A(:,ip))**5
    end do
  end subroutine f_electric_field_gradient_diag

  pure subroutine f_electric_field_gradient_full(n_compute, n_points, &
       & partition, wave, r_A, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,3)
    integer :: i_dir, ip
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir,i_dir) = (r_A(i_dir,ip)*r_A(i_dir,ip) - &
               & sum(r_A(:,ip)*r_A(:,ip))/3)*wave(:,ip)
       end do
       matrix_wave(:,ip,1,2) = r_A(1,ip)*r_A(2,ip)*wave(:,ip)
       matrix_wave(:,ip,1,3) = r_A(1,ip)*r_A(3,ip)*wave(:,ip)
       matrix_wave(:,ip,2,1) = r_A(2,ip)*r_A(1,ip)*wave(:,ip)
       matrix_wave(:,ip,2,3) = r_A(2,ip)*r_A(3,ip)*wave(:,ip)
       matrix_wave(:,ip,3,1) = r_A(3,ip)*r_A(1,ip)*wave(:,ip)
       matrix_wave(:,ip,3,2) = r_A(3,ip)*r_A(2,ip)*wave(:,ip)
       matrix_wave(:,ip,:,:) = &
            & matrix_wave(:,ip,:,:)*partition(ip)/norm2(r_A(:,ip))**5
    end do
  end subroutine f_electric_field_gradient_full
  !!
  !! Orbital zeeman
  !!
  !!       i
  !! H = - - r x Nabla
  !!       2
  !!
  pure subroutine f_orbital_zeeman(n_compute, n_points, partition, nabla_wave, &
       & r, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: nabla_wave(n_compute,3,n_points)
    real(dp), intent(in) :: r(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: ip
    do ip = 1, n_points
       matrix_wave(:,ip,1) = &
            & r(2,ip)*nabla_wave(:,3,ip) - r(3,ip)*nabla_wave(:,2,ip)
       matrix_wave(:,ip,2) = &
            & r(3,ip)*nabla_wave(:,1,ip) - r(1,ip)*nabla_wave(:,3,ip)
       matrix_wave(:,ip,3) = &
            & r(1,ip)*nabla_wave(:,2,ip) - r(2,ip)*nabla_wave(:,1,ip)
       matrix_wave(:,ip,:) = matrix_wave(:,ip,:)*partition(ip)
    end do
  end subroutine f_orbital_zeeman
  !!
  !! Diamagnetic shielding in the standard formalism (no GIAOs)
  !!
  !!     alpha^2 r r_A - r_A r^T
  !! H = ------- ---------------
  !!        2         r_A^3
  !!
  pure subroutine f_diamagnetic_shielding_std_diag(n_compute, n_points, &
       & partition, wave, r_0, r_A, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_0(3,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: i_dir, ip
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir) = (sum(r_0(:,ip)*r_A(:,ip)) - &
               & r_A(i_dir,ip)*r_0(i_dir,ip))*wave(:,ip)
       end do
       matrix_wave(:,ip,:) = &
            & matrix_wave(:,ip,:)*partition(ip)/norm2(r_A(:,ip))**3
    end do
  end subroutine f_diamagnetic_shielding_std_diag

  pure subroutine f_diamagnetic_shielding_std_full(n_compute, n_points, &
       & partition, wave, r_0, r_A, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_0(3,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,3)
    integer :: ip, i_dir
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir,i_dir) = (sum(r_0(:,ip)*r_A(:,ip)) - &
               & r_A(i_dir,ip)*r_0(i_dir,ip))*wave(:,ip)
       end do
       matrix_wave(:,ip,1,2) = -r_A(1,ip)*r_0(2,ip)*wave(:,ip)
       matrix_wave(:,ip,1,3) = -r_A(1,ip)*r_0(3,ip)*wave(:,ip)
       matrix_wave(:,ip,2,1) = -r_A(2,ip)*r_0(1,ip)*wave(:,ip)
       matrix_wave(:,ip,2,3) = -r_A(2,ip)*r_0(3,ip)*wave(:,ip)
       matrix_wave(:,ip,3,1) = -r_A(3,ip)*r_0(1,ip)*wave(:,ip)
       matrix_wave(:,ip,3,2) = -r_A(3,ip)*r_0(2,ip)*wave(:,ip)
       matrix_wave(:,ip,:,:) = &
            & matrix_wave(:,ip,:,:)*partition(ip)/norm2(r_A(:,ip))**3
    end do
  end subroutine f_diamagnetic_shielding_std_full
  !!
  !! First-order GIAO overlap matrix
  !!
  !!          1
  !! S1_mn = --- R_mn x <m|r|n>,
  !!          2
  !!
  !! where R_mn = R_m - R_n, where R_m is the origin of basis function
  !! m. The cross-product is preformed in GIAO_mult_psi_and_pack.f90
  !!
  pure subroutine f_giao_overlap1(n_compute, n_points, partition, wave, &
       & grid_points, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: grid_points(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: i_dir, ip
    do ip = 1, n_points
       matrix_wave(:,ip,3) = partition(ip)*wave(:,ip)
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir) = grid_points(i_dir,ip)*matrix_wave(:,ip,3)
       end do
    end do
  end subroutine f_giao_overlap1
  !!
  !! Second-order GIAO overlap matrix
  !!
  !!          1
  !! S2_mn = --- R_mn x <m|rr^T|n> x R_mn.
  !!          4
  !!
  !! Same comments as with f_giao_overlap1.
  !!
  pure subroutine f_giao_overlap2(n_compute, n_points, partition, wave, r, &
       & matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,2)
    integer :: ip, i_dir
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir,1) = r(i_dir,ip)*r(i_dir,ip)*wave(:,ip)
       end do
       matrix_wave(:,ip,1,2) = r(1,ip)*r(2,ip)*wave(:,ip)
       matrix_wave(:,ip,2,2) = r(1,ip)*r(3,ip)*wave(:,ip)
       matrix_wave(:,ip,3,2) = r(2,ip)*r(3,ip)*wave(:,ip)
       matrix_wave(:,ip,:,:) = partition(ip)*matrix_wave(:,ip,:,:)
    end do
  end subroutine f_giao_overlap2
  !!
  !! Paramagnetic GIAO integrals. These are required for paramagnetic
  !! shielding and the magnetizability in the GIAO formalism.
  !!
  !!        i
  !! H = - --- (H1+H2)
  !!        2
  !!
  !! This subroutine calculates the H1 part, which does not include
  !! any R_mn terms:
  !!
  !! <m|H1|n> = <m|r_n x Nabla |n>
  !!
  pure subroutine f_giao_paramagnetic_0(n_compute, n_points, partition, &
       & nabla_wave, r_n, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: nabla_wave(n_compute,3,n_points)
    real(dp), intent(in) :: r_n(n_compute,3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: ip
    do ip = 1, n_points
       matrix_wave(:,ip,1) = &
            & r_n(:,2,ip)*nabla_wave(:,3,ip) - r_n(:,3,ip)*nabla_wave(:,2,ip)
       matrix_wave(:,ip,2) = &
            & r_n(:,3,ip)*nabla_wave(:,1,ip) - r_n(:,1,ip)*nabla_wave(:,3,ip)
       matrix_wave(:,ip,3) = &
            & r_n(:,1,ip)*nabla_wave(:,2,ip) - r_n(:,2,ip)*nabla_wave(:,1,ip)
       matrix_wave(:,ip,:) = matrix_wave(:,ip,:)*partition(ip)
    end do
  end subroutine f_giao_paramagnetic_0
  !!
  !! <m|H2|n> = - R_mn x <m|r H0|n> (LDA)
  !!
  !! Same comments as with f_giao_overlap1.
  !! Note: with pure Hartree-Fock, the xc-potential, vxc, is empty.
  !!
  pure subroutine f_giao_paramagnetic_lda(n_compute, n_points, partition, &
       & wave, kinetic_wave, hartree_pot, vxc, grid_points, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: kinetic_wave(n_compute,n_points)
    real(dp), intent(in) :: hartree_pot(n_points)
    real(dp), intent(in) :: vxc(2,n_points)
    real(dp), intent(in) :: grid_points(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,n_spin)
    integer :: i_dir, ip
    if (n_spin == 1) then
       do ip = 1, n_points
          matrix_wave(:,ip,3,1) = &
               & -kinetic_wave(:,ip) - (hartree_pot(ip) + vxc(1,ip))*wave(:,ip)
          do i_dir = 1, 3
             matrix_wave(:,ip,i_dir,1) = &
                  & grid_points(i_dir,ip)*matrix_wave(:,ip,3,1)
          end do
          matrix_wave(:,ip,:,1) = matrix_wave(:,ip,:,1)*partition(ip)
       end do
    else
       do ip = 1, n_points
          matrix_wave(:,ip,3,1) = &
               & -kinetic_wave(:,ip) - (hartree_pot(ip) + vxc(1,ip))*wave(:,ip)
          matrix_wave(:,ip,3,2) = &
               & -kinetic_wave(:,ip) - (hartree_pot(ip) + vxc(2,ip))*wave(:,ip)
          do i_dir = 1, 3
             matrix_wave(:,ip,i_dir,:) = &
                  & grid_points(i_dir,ip)*matrix_wave(:,ip,3,:)
          end do
          matrix_wave(:,ip,:,:) = matrix_wave(:,ip,:,:)*partition(ip)
       end do
    end if
  end subroutine f_giao_paramagnetic_lda
  !!
  !! Calculates the LDA term, H0, and the first GGA term H2, where
  !!
  !! <m|H2|n> = - R_mn x <m|r H0|n> - R_mn x <m|H3 + H4|n> (GGA)
  !!
  !! <m|H3|n> = 2 Int Nabla rho(r) vsigma(r) Omega_ij(r) D r
  !!
  !! H4 is calculated in a different subroutine.
  !! vsigma(r) = D Exc/D gamma(r)
  !! gamma = [ |Nabla rho_1|^2, Nabla rho_1 Nabla rho_2, |Nabla rho_2|^2 ]
  !! Omega_ij(r) = wave_i(r)*wave_j(r)
  !!
  pure subroutine f_giao_paramagnetic_gga(n_compute, n_points, partition, &
       & wave, kinetic_wave, hartree_pot, vxc, vsigma, &
       & grid_points, local_rho_gradient, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: kinetic_wave(n_compute,n_points)
    real(dp), intent(in) :: hartree_pot(n_points)
    real(dp), intent(in) :: vxc(2,n_points)
    real(dp), intent(in) :: vsigma(3,n_points)
    real(dp), intent(in) :: grid_points(3,n_points)
    real(dp), intent(in) :: local_rho_gradient(3,2,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,n_spin)
    integer :: i_dir, ip
    if (n_spin == 1) then
       do ip = 1, n_points
          ! r H0
          matrix_wave(:,ip,3,1) = &
               & -kinetic_wave(:,ip) - (hartree_pot(ip) + vxc(1,ip))*wave(:,ip)
          do i_dir = 1, 3
             matrix_wave(:,ip,i_dir,1) = &
                  & grid_points(i_dir,ip)*matrix_wave(:,ip,3,1)
          end do
          ! First GGA term [H3]
          do i_dir = 1, 3
             matrix_wave(:,ip,i_dir,1) = matrix_wave(:,ip,i_dir,1) - &
                  & (2d0*vsigma(1,ip)*local_rho_gradient(i_dir,1,ip) + &
                  & vsigma(2,ip)*local_rho_gradient(i_dir,2,ip))*wave(:,ip)
          end do
          matrix_wave(:,ip,:,1) = matrix_wave(:,ip,:,1)*partition(ip)
       end do
    else
       do ip = 1, n_points
          ! r H0
          matrix_wave(:,ip,3,1) = &
               & -kinetic_wave(:,ip) - (hartree_pot(ip) + vxc(1,ip))*wave(:,ip)
          matrix_wave(:,ip,3,2) = &
               & -kinetic_wave(:,ip) - (hartree_pot(ip) + vxc(2,ip))*wave(:,ip)
          do i_dir = 1, 3
             matrix_wave(:,ip,i_dir,:) = &
                  & grid_points(i_dir,ip)*matrix_wave(:,ip,3,:)
          end do
          ! First GGA term [H3]
          do i_dir = 1, 3
             matrix_wave(:,ip,i_dir,1) = matrix_wave(:,ip,i_dir,1) - &
                  & (2d0*vsigma(1,ip)*local_rho_gradient(i_dir,1,ip) + &
                  & vsigma(2,ip)*local_rho_gradient(i_dir,2,ip))*wave(:,ip)
             matrix_wave(:,ip,i_dir,2) = matrix_wave(:,ip,i_dir,2) - &
                  & (2d0*vsigma(3,ip)*local_rho_gradient(i_dir,2,ip) + &
                  & vsigma(2,ip)*local_rho_gradient(i_dir,1,ip))*wave(:,ip)
          end do
          matrix_wave(:,ip,:,:) = matrix_wave(:,ip,:,:)*partition(ip)
       end do
    end if
  end subroutine f_giao_paramagnetic_gga
  !!
  !! The second GGA term, H4, which includes the gradient:
  !!
  !! <m|H4|n> = 2 Int vsigma(r) Nabla rho(r) Nabla [Omega_mn(r)] D r
  !!
  !! Note: H4+H4^T must be taken after this.
  !!
  pure subroutine f_giao_paramagnetic_gga_grad(n_compute, n_points, &
       & partition, nabla_wave, vsigma, matrix_wave, grid_points, &
       & local_rho_gradient)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: nabla_wave(n_compute,3,n_points)
    real(dp), intent(in) :: grid_points(3,n_points)
    real(dp), intent(in) :: vsigma(3,n_points)
    real(dp), intent(in) :: local_rho_gradient(3,2,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,n_spin)
    integer :: i_dir, ip
    if (n_spin == 1) then
       do ip = 1, n_points
          ! Second GGA term
          matrix_wave(:,ip,3,1) = -matmul(nabla_wave(:,:,ip), &
               & 2d0*vsigma(1,ip)*local_rho_gradient(:,1,ip) + &
               & vsigma(2,ip)*local_rho_gradient(:,2,ip))
          do i_dir = 1, 3
             matrix_wave(:,ip,i_dir,1) = &
                  & grid_points(i_dir,ip)*matrix_wave(:,ip,3,1)
          end do
          matrix_wave(:,ip,:,1) = matrix_wave(:,ip,:,1)*partition(ip)
       end do
    else
       do ip = 1, n_points
          ! Second GGA term
          matrix_wave(:,ip,3,1) = -matmul(nabla_wave(:,:,ip), &
               & 2d0*vsigma(1,ip)*local_rho_gradient(:,1,ip) + &
               & vsigma(2,ip)*local_rho_gradient(:,2,ip))
          matrix_wave(:,ip,3,2) = -matmul(nabla_wave(:,:,ip), &
               & 2d0*vsigma(3,ip)*local_rho_gradient(:,2,ip) + &
               & vsigma(2,ip)*local_rho_gradient(:,1,ip))
          do i_dir = 1, 3
             matrix_wave(:,ip,i_dir,:) = &
                  & grid_points(i_dir,ip)*matrix_wave(:,ip,3,:)
          end do
          matrix_wave(:,ip,:,:) = matrix_wave(:,ip,:,:)*partition(ip)
       end do
    end if
  end subroutine f_giao_paramagnetic_gga_grad
  !!
  !! Diamagnetic shielding, with GIAOs
  !!
  !!     alpha^2
  !! H = ------- (H1+H2)
  !!        2
  !!                         r (r_A x Nabla)^T
  !! <m|H1|n> = R_mn x < m | ----------------- | n >,
  !!                               r_A^3
  !!
  !! Same comments as with f_giao_overlap1.
  !!
  pure subroutine f_diamagnetic_shielding_1_diag(n_compute, n_points, &
       & partition, nabla_wave, r, r_A, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: nabla_wave(n_compute,3,n_points)
    real(dp), intent(in) :: r(3,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,2,3)
    integer :: ip
    do ip = 1, n_points
       matrix_wave(:,ip,1,1) = &
            & (r_A(2,ip)*nabla_wave(:,3,ip) - r_A(3,ip)*nabla_wave(:,2,ip))
       matrix_wave(:,ip,1,2) = &
            & (r_A(3,ip)*nabla_wave(:,1,ip) - r_A(1,ip)*nabla_wave(:,3,ip))
       matrix_wave(:,ip,1,3) = &
            & (r_A(1,ip)*nabla_wave(:,2,ip) - r_A(2,ip)*nabla_wave(:,1,ip))
       ! Compact indexing is used here
       matrix_wave(:,ip,2,3) = r(2,ip)*matrix_wave(:,ip,1,3) ! 2 3
       matrix_wave(:,ip,1,3) = r(1,ip)*matrix_wave(:,ip,1,3) ! 1 3
       matrix_wave(:,ip,2,2) = r(3,ip)*matrix_wave(:,ip,1,2) ! 3 2
       matrix_wave(:,ip,1,2) = r(1,ip)*matrix_wave(:,ip,1,2) ! 1 2
       matrix_wave(:,ip,2,1) = r(3,ip)*matrix_wave(:,ip,1,1) ! 3 1
       matrix_wave(:,ip,1,1) = r(2,ip)*matrix_wave(:,ip,1,1) ! 2 1
       matrix_wave(:,ip,:,:) = &
            & matrix_wave(:,ip,:,:)*partition(ip)/norm2(r_A(:,ip))**3
    end do
  end subroutine f_diamagnetic_shielding_1_diag

  pure subroutine f_diamagnetic_shielding_1_full(n_compute, n_points, &
       & partition, nabla_wave, r, r_A, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: nabla_wave(n_compute,3,n_points)
    real(dp), intent(in) :: r(3,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,3)
    integer :: i_dir, j_dir, ip
    do ip = 1, n_points
       matrix_wave(:,ip,1,1) = &
            & (r_A(2,ip)*nabla_wave(:,3,ip) - r_A(3,ip)*nabla_wave(:,2,ip))
       matrix_wave(:,ip,1,2) = &
            & (r_A(3,ip)*nabla_wave(:,1,ip) - r_A(1,ip)*nabla_wave(:,3,ip))
       matrix_wave(:,ip,1,3) = &
            & (r_A(1,ip)*nabla_wave(:,2,ip) - r_A(2,ip)*nabla_wave(:,1,ip))
       do j_dir = 1, 3
          ! The i_dir=1 case should be the last one to be overwritten.
          do i_dir = 3, 1, -1
             matrix_wave(:,ip,i_dir,j_dir) = &
                  & r(i_dir,ip)*matrix_wave(:,ip,1,j_dir)
          end do
       end do
       matrix_wave(:,ip,:,:) = &
            & matrix_wave(:,ip,:,:)*partition(ip)/norm2(r_A(:,ip))**3
    end do
  end subroutine f_diamagnetic_shielding_1_full
  !!
  !!                  r_n r_A - r_A r_n^T
  !! <m|H2|n> = < m | ------------------- | n >
  !!                         r_A^3
  !!
  pure subroutine f_diamagnetic_shielding_2_diag(n_compute, n_points, &
       & partition, wave, r_A, r_n, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(in) :: r_n(n_compute,3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: i_dir, ip
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir) = (matmul(r_n(:,:,ip),r_A(:,ip)) - &
               & r_A(i_dir,ip)*r_n(:,i_dir,ip))*wave(:,ip)
       end do
       matrix_wave(:,ip,:) = &
            & matrix_wave(:,ip,:)*partition(ip)/norm2(r_A(:,ip))**3
    end do
  end subroutine f_diamagnetic_shielding_2_diag

  pure subroutine f_diamagnetic_shielding_2_full(n_compute, n_points, &
       & partition, wave, r_A, r_n, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_A(3,n_points)
    real(dp), intent(in) :: r_n(n_compute,3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,3)
    integer :: i_dir, ip
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir,i_dir) = (matmul(r_n(:,:,ip),r_A(:,ip)) - &
               & r_A(i_dir,ip)*r_n(:,i_dir,ip))*wave(:,ip)
       end do
       matrix_wave(:,ip,1,2) = -r_A(1,ip)*r_n(:,2,ip)*wave(:,ip)
       matrix_wave(:,ip,1,3) = -r_A(1,ip)*r_n(:,3,ip)*wave(:,ip)
       matrix_wave(:,ip,2,1) = -r_A(2,ip)*r_n(:,1,ip)*wave(:,ip)
       matrix_wave(:,ip,2,3) = -r_A(2,ip)*r_n(:,3,ip)*wave(:,ip)
       matrix_wave(:,ip,3,1) = -r_A(3,ip)*r_n(:,1,ip)*wave(:,ip)
       matrix_wave(:,ip,3,2) = -r_A(3,ip)*r_n(:,2,ip)*wave(:,ip)
       matrix_wave(:,ip,:,:) = &
            & matrix_wave(:,ip,:,:)*partition(ip)/norm2(r_A(:,ip))**3
    end do
  end subroutine f_diamagnetic_shielding_2_full
  !!
  !! Diamagnetic magnetizability in the standard formalism (no GIAOs)
  !!
  !!        1
  !! H = - --- (r^2 - rr^T)
  !!        4
  !!
  pure subroutine f_diamagnetic_magnetizability_std_diag(n_compute, n_points, &
       & partition, wave, r_0, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_0(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: i_dir, ip
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir) = &
               & (norm2(r_0(:,ip))**2 - r_0(i_dir,ip)**2)*wave(:,ip)
       end do
       matrix_wave(:,ip,:) = matrix_wave(:,ip,:)*partition(ip)
    end do
  end subroutine f_diamagnetic_magnetizability_std_diag

  pure subroutine f_diamagnetic_magnetizability_std_full(n_compute, n_points, &
       & partition, wave, r_0, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_0(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,3)
    integer :: i_dir, ip
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir,i_dir) = &
               & (norm2(r_0(:,ip))**2 - r_0(i_dir,ip)**2)*wave(:,ip)
       end do
       matrix_wave(:,ip,1,2) = -r_0(1,ip)*r_0(2,ip)*wave(:,ip)
       matrix_wave(:,ip,1,3) = -r_0(1,ip)*r_0(3,ip)*wave(:,ip)
       matrix_wave(:,ip,2,3) = -r_0(2,ip)*r_0(3,ip)*wave(:,ip)
       matrix_wave(:,ip,:,:) = matrix_wave(:,ip,:,:)*partition(ip)
    end do
    ! This could be further optimized but the non-GIAO version is for
    ! testing purposes anyway.
    matrix_wave(:,:,3,1) = matrix_wave(:,:,1,3)
    matrix_wave(:,:,3,2) = matrix_wave(:,:,2,3)
    matrix_wave(:,:,2,1) = matrix_wave(:,:,1,2)
  end subroutine f_diamagnetic_magnetizability_std_full
  !!
  !! Diamagnetic magnetizability, with GIAOs
  !!
  !! <m|H|n> = <m|H1|n> + <m|H2|n> + <m|H3|n> + <m|H4|n> + <m|H5|n>
  !!
  !! Five different terms (three in LDA) make up the diamagnetic
  !! magnetizability integral. This subroutine calculates the H1 part.
  !!
  !!             1
  !! <m|H1|n> = --- < m | r_n^2 - r_n r_n^T | n >
  !!             4
  !!
  !! Note: this term does not contain any cross products with R_mn.
  !!
  pure subroutine f_diamagnetic_magnetizability_0_diag(n_compute, n_points, &
       & partition, wave, r_n, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_n(n_compute,3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: i_dir, ip
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir) = &
               & (sum(r_n(:,:,ip)**2,2) - r_n(:,i_dir,ip)**2)*wave(:,ip)
       end do
       matrix_wave(:,ip,:) = partition(ip)*matrix_wave(:,ip,:)
    end do
  end subroutine f_diamagnetic_magnetizability_0_diag

  pure subroutine f_diamagnetic_magnetizability_0_full(n_compute, n_points, &
       & partition, wave, r_n, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: r_n(n_compute,3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,3)
    integer :: i_dir, ip
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir,i_dir) = &
               & (sum(r_n(:,:,ip)**2,2) - r_n(:,i_dir,ip)**2)*wave(:,ip)
       end do
       matrix_wave(:,ip,2,1) = -r_n(:,2,ip)*r_n(:,1,ip)*wave(:,ip)
       matrix_wave(:,ip,3,1) = -r_n(:,3,ip)*r_n(:,1,ip)*wave(:,ip)
       matrix_wave(:,ip,1,2) = -r_n(:,1,ip)*r_n(:,2,ip)*wave(:,ip)
       matrix_wave(:,ip,3,2) = -r_n(:,3,ip)*r_n(:,2,ip)*wave(:,ip)
       matrix_wave(:,ip,1,3) = -r_n(:,1,ip)*r_n(:,3,ip)*wave(:,ip)
       matrix_wave(:,ip,2,3) = -r_n(:,2,ip)*r_n(:,3,ip)*wave(:,ip)
       matrix_wave(:,ip,:,:) = partition(ip)*matrix_wave(:,ip,:,:)
    end do
  end subroutine f_diamagnetic_magnetizability_0_full
  !!
  !!             1
  !! <m|H2|n> = --- R_mn x < m | 2r (r_n x Nabla)^T | n >
  !!             4
  !!
  !! Note: this term includes one R_mn cross product.
  !!
  pure subroutine f_diamagnetic_magnetizability_R_mn_diag(n_compute, n_points, &
       & partition, nabla_wave, r, r_n, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: nabla_wave(n_compute,3,n_points)
    real(dp), intent(in) :: r(3,n_points)
    real(dp), intent(in) :: r_n(n_compute,3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,2,3)
    integer :: ip
    do ip = 1, n_points
       matrix_wave(:,ip,1,1) = &
            & (r_n(:,2,ip)*nabla_wave(:,3,ip) - r_n(:,3,ip)*nabla_wave(:,2,ip))
       matrix_wave(:,ip,1,2) = &
            & (r_n(:,3,ip)*nabla_wave(:,1,ip) - r_n(:,1,ip)*nabla_wave(:,3,ip))
       matrix_wave(:,ip,1,3) = &
            & (r_n(:,1,ip)*nabla_wave(:,2,ip) - r_n(:,2,ip)*nabla_wave(:,1,ip))
       ! Compact indexing is used here
       matrix_wave(:,ip,2,3) = r(2,ip)*matrix_wave(:,ip,1,3) ! 2 3
       matrix_wave(:,ip,1,3) = r(1,ip)*matrix_wave(:,ip,1,3) ! 1 3
       matrix_wave(:,ip,2,2) = r(3,ip)*matrix_wave(:,ip,1,2) ! 3 2
       matrix_wave(:,ip,1,2) = r(1,ip)*matrix_wave(:,ip,1,2) ! 1 2
       matrix_wave(:,ip,2,1) = r(3,ip)*matrix_wave(:,ip,1,1) ! 3 1
       matrix_wave(:,ip,1,1) = r(2,ip)*matrix_wave(:,ip,1,1) ! 2 1
       matrix_wave(:,ip,:,:) = 2*matrix_wave(:,ip,:,:)*partition(ip)
    end do
  end subroutine f_diamagnetic_magnetizability_R_mn_diag

  pure subroutine f_diamagnetic_magnetizability_R_mn_full(n_compute, n_points, &
       & partition, nabla_wave, r, r_n, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: nabla_wave(n_compute,3,n_points)
    real(dp), intent(in) :: r(3,n_points)
    real(dp), intent(in) :: r_n(n_compute,3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,3)
    integer :: i_dir, j_dir, ip
    do ip = 1, n_points
       matrix_wave(:,ip,1,1) = &
            & (r_n(:,2,ip)*nabla_wave(:,3,ip) - r_n(:,3,ip)*nabla_wave(:,2,ip))
       matrix_wave(:,ip,1,2) = &
            & (r_n(:,3,ip)*nabla_wave(:,1,ip) - r_n(:,1,ip)*nabla_wave(:,3,ip))
       matrix_wave(:,ip,1,3) = &
            & (r_n(:,1,ip)*nabla_wave(:,2,ip) - r_n(:,2,ip)*nabla_wave(:,1,ip))
       do j_dir = 1, 3
          ! The i_dir=1 case should be the last one to be overwritten.
          do i_dir = 3, 1, -1
             matrix_wave(:,ip,i_dir,j_dir) = &
                  & r(i_dir,ip)*matrix_wave(:,ip,1,j_dir)
          end do
       end do
       matrix_wave(:,ip,:,:) = 2*matrix_wave(:,ip,:,:)*partition(ip)
    end do
  end subroutine f_diamagnetic_magnetizability_R_mn_full
  !!
  !!             1
  !! <m|H3|n> = --- R_mn x < m | r H0 r^T | n > x R_mn
  !!             4
  !!
  !! Note: this term includes two cross products. H0 is the LDA part
  !! of the unperturbed Hamiltonian. Two GGA terms, if necessary, are
  !! calculated afterwards.
  !!
  pure subroutine f_diamagnetic_magnetizability_lda(n_compute, n_points, &
       & partition, wave, kinetic_wave, hartree_pot, vxc, r, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: kinetic_wave(n_compute,n_points)
    real(dp), intent(in) :: hartree_pot(n_points)
    real(dp), intent(in) :: vxc(2,n_points)
    real(dp), intent(in) :: r(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,2,n_spin)
    integer :: ip, i_dir
    if (n_spin == 1) then
       do ip = 1, n_points
          matrix_wave(:,ip,3,2,1) = &
               & (kinetic_wave(:,ip) + (hartree_pot(ip) + vxc(1,ip))*wave(:,ip))
          do i_dir = 1, 3
             matrix_wave(:,ip,i_dir,1,1) = &
                  & r(i_dir,ip)*r(i_dir,ip)*matrix_wave(:,ip,3,2,1)
          end do
          matrix_wave(:,ip,1,2,1) = r(1,ip)*r(2,ip)*matrix_wave(:,ip,3,2,1)
          matrix_wave(:,ip,2,2,1) = r(1,ip)*r(3,ip)*matrix_wave(:,ip,3,2,1)
          matrix_wave(:,ip,3,2,1) = r(2,ip)*r(3,ip)*matrix_wave(:,ip,3,2,1)
          matrix_wave(:,ip,:,:,1) = partition(ip)*matrix_wave(:,ip,:,:,1)
       end do
    else
       do ip = 1, n_points
          matrix_wave(:,ip,3,2,1) = &
               & (kinetic_wave(:,ip) + (hartree_pot(ip) + vxc(1,ip))*wave(:,ip))
          matrix_wave(:,ip,3,2,2) = &
               & (kinetic_wave(:,ip) + (hartree_pot(ip) + vxc(2,ip))*wave(:,ip))
          do i_dir = 1, 3
             matrix_wave(:,ip,i_dir,1,:) = &
                  & r(i_dir,ip)*r(i_dir,ip)*matrix_wave(:,ip,3,2,:)
          end do
          matrix_wave(:,ip,1,2,:) = r(1,ip)*r(2,ip)*matrix_wave(:,ip,3,2,:)
          matrix_wave(:,ip,2,2,:) = r(1,ip)*r(3,ip)*matrix_wave(:,ip,3,2,:)
          matrix_wave(:,ip,3,2,:) = r(2,ip)*r(3,ip)*matrix_wave(:,ip,3,2,:)
          matrix_wave(:,ip,:,:,:) = partition(ip)*matrix_wave(:,ip,:,:,:)
       end do
    end if
  end subroutine f_diamagnetic_magnetizability_lda
  !!
  !! For spin channel a:
  !!
  !!             1
  !! <m|H4|n> = --- R_mn x < m | xi(r) r^T | n > x R_mn,
  !!             2
  !!
  !! where xi = 2(partial Exc/partial gamma_aa) Nabla rho_a(r)
  !!           + (partial Exc/partial gamma_ab) Nabla rho_b(r)
  !!
  !! (xi is a vector.) Likewise for spin channel b.
  !! Note: the reason a factor of two is seen below is to cancel a
  !! factor of 1/4 later.
  !!
  pure subroutine f_diamagnetic_magnetizability_gga(n_compute, n_points, &
       & partition, wave, vsigma, r, local_rho_gradient, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: vsigma(3,n_points)
    real(dp), intent(in) :: r(3,n_points)
    real(dp), intent(in) :: local_rho_gradient(3,2,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,3,n_spin)
    integer :: ip, i_dir, j_dir
    if (n_spin == 1) then
       do j_dir = 1, 3
          do i_dir = 1, 3
             do ip = 1, n_points
                matrix_wave(:,ip,i_dir,j_dir,1) = 2*partition(ip)* &
                     & (2d0*vsigma(1,ip)*local_rho_gradient(i_dir,1,ip) + &
                     & vsigma(2,ip)*local_rho_gradient(i_dir,2,ip))* &
                     & r(j_dir,ip)*wave(:,ip)
             end do
          end do
       end do
    else
       do j_dir = 1, 3
          do i_dir = 1, 3
             do ip = 1, n_points
                matrix_wave(:,ip,i_dir,j_dir,1) = 2*partition(ip)* &
                     & (2d0*vsigma(1,ip)*local_rho_gradient(i_dir,1,ip) + &
                     & vsigma(2,ip)*local_rho_gradient(i_dir,2,ip))* &
                     & r(j_dir,ip)*wave(:,ip)
                matrix_wave(:,ip,i_dir,j_dir,2) = 2*partition(ip)* &
                     & (2d0*vsigma(3,ip)*local_rho_gradient(i_dir,2,ip) + &
                     & vsigma(2,ip)*local_rho_gradient(i_dir,1,ip))* &
                     & r(j_dir,ip)*wave(:,ip)
             end do
          end do
       end do
    end if
  end subroutine f_diamagnetic_magnetizability_gga
  !!
  !!             1
  !! <m|H5|n> = --- R_mn x I_mn x R_mn,
  !!             4
  !!
  !! where I_mn = Int r xi(r) r^T Nabla[phi_m(r)phi_n(r)] D r,
  !! where xi was defined in f_diamagnetic_magnetizability_gga.
  !! Note: H5+H5^T must be taken after this.
  !!
  pure subroutine f_diamagnetic_magnetizability_gga_grad(n_compute, &
       & n_points, partition, nabla_wave, vsigma, r, local_rho_gradient, &
       & matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: nabla_wave(n_compute,3,n_points)
    real(dp), intent(in) :: vsigma(3,n_points)
    real(dp), intent(in) :: r(3,n_points)
    real(dp), intent(in) :: local_rho_gradient(3,2,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,2,n_spin)
    integer :: ip, i_dir
    if (n_spin == 1) then
       do ip = 1, n_points
          matrix_wave(:,ip,3,2,1) = matmul(nabla_wave(:,:,ip), &
               & 2*vsigma(1,ip)*local_rho_gradient(:,1,ip) + &
               & vsigma(2,ip)*local_rho_gradient(:,2,ip))
          do i_dir = 1, 3
             matrix_wave(:,ip,i_dir,1,1) = &
                  & r(i_dir,ip)*r(i_dir,ip)*matrix_wave(:,ip,3,2,1)
          end do
          matrix_wave(:,ip,1,2,1) = r(1,ip)*r(2,ip)*matrix_wave(:,ip,3,2,1)
          matrix_wave(:,ip,2,2,1) = r(1,ip)*r(3,ip)*matrix_wave(:,ip,3,2,1)
          matrix_wave(:,ip,3,2,1) = r(2,ip)*r(3,ip)*matrix_wave(:,ip,3,2,1)
          matrix_wave(:,ip,:,:,1) = partition(ip)*matrix_wave(:,ip,:,:,1)
       end do
    else
       do ip = 1, n_points
          matrix_wave(:,ip,3,2,1) = matmul(nabla_wave(:,:,ip), &
               & 2*vsigma(1,ip)*local_rho_gradient(:,1,ip) + &
               & vsigma(2,ip)*local_rho_gradient(:,2,ip))
          matrix_wave(:,ip,3,2,2) = matmul(nabla_wave(:,:,ip), &
               & 2*vsigma(3,ip)*local_rho_gradient(:,2,ip) + &
               & vsigma(2,ip)*local_rho_gradient(:,1,ip))
          do i_dir = 1, 3
             matrix_wave(:,ip,i_dir,1,:) = &
                  & r(i_dir,ip)*r(i_dir,ip)*matrix_wave(:,ip,3,2,:)
          end do
          matrix_wave(:,ip,1,2,:) = r(1,ip)*r(2,ip)*matrix_wave(:,ip,3,2,:)
          matrix_wave(:,ip,2,2,:) = r(1,ip)*r(3,ip)*matrix_wave(:,ip,3,2,:)
          matrix_wave(:,ip,3,2,:) = r(2,ip)*r(3,ip)*matrix_wave(:,ip,3,2,:)
          matrix_wave(:,ip,:,:,:) = partition(ip)*matrix_wave(:,ip,:,:,:)
       end do
    end if
  end subroutine f_diamagnetic_magnetizability_gga_grad
  !!
  !! First-order contribution from the LDA potential
  !!
  !! <i|V_xc1|j> = Int xc_kernel_s(r) Omega_ij(r) D r,
  !!
  !! where Omega_ij(r) = wave_i(r)*wave_j(r) and
  !!
  pure subroutine f_xc_kernel_LDA(n_compute, n_points, partition, wave, &
       & xc_kernel_s, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: xc_kernel_s(n_spin,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,n_spin)
    integer ::  i_spin, ip
    do i_spin = 1, n_spin
       do ip = 1, n_points
          matrix_wave(:,ip,i_spin) = &
               & partition(ip)*xc_kernel_s(i_spin,ip)*wave(:,ip)
       end do
    end do
  end subroutine f_xc_kernel_LDA
  !!
  !! First-order contribution from the GGA potential
  !!
  !! <i|V_xc1|j> = Int xc_kernel_s(r) Omega_ij(r) D r
  !!               + Int xc_kernel_v(r) Nabla Omega_ij(r) D r,
  !!
  !! where Omega_ij(r) = wave_i(r)*wave_j(r) and
  !! Nabla Omega_ij(r) = 2 wave_i(r)*Nabla wave_j(r)
  !!
  pure subroutine f_xc_kernel_GGA(n_compute, n_points, partition, wave, &
         & nabla_wave, xc_kernel_s, xc_kernel_v, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: nabla_wave(n_compute,3,n_points)
    real(dp), intent(in) :: xc_kernel_s(n_spin,n_points)
    real(dp), intent(in) :: xc_kernel_v(3,n_spin,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,n_spin)
    integer :: i_spin, ip
    do i_spin = 1, n_spin
       do ip = 1, n_points
          matrix_wave(:,ip,i_spin) = &
               & partition(ip)*(xc_kernel_s(i_spin,ip)*wave(:,ip) + &
               & 2d0*matmul(nabla_wave(:,:,ip),xc_kernel_v(:,i_spin,ip)))
       end do
    end do
  end subroutine f_xc_kernel_GGA
  !!
  !! First-order Hartree potential
  !!
  pure subroutine f_first_order_hartree(n_dirs, n_compute, n_points, &
       & partition, wave, first_order_hartree, matrix_wave)
    integer, intent(in) :: n_dirs, n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: first_order_hartree(n_points,n_dirs)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,n_dirs)
    integer :: i_dir, ip
    do i_dir = 1, n_dirs
       do ip = 1, n_points
          matrix_wave(:,ip,i_dir) = &
               & partition(ip)*first_order_hartree(ip,i_dir)*wave(:,ip)
       end do
    end do
  end subroutine f_first_order_hartree
  !!
  !! Constructs a density and optionally a density gradient using the
  !! density matrix:
  !!
  !! rho(r) = Sum_ij phi_i(r) n_ij phi_j(r)
  !! rho_gradient(r) = 2*Sum_ij phi_i(r) n_ij Nabla phi_j(r)
  !!
  !! If the input is a first-order density matrix, the outputs are
  !! also first-order quantities.
  !!
  pure subroutine f_electron_density(n_compute, n_points, i_spin, &
       & wave, density_matrix, n_phi, local_rho, nabla_wave, local_rho_gradient)
    integer, intent(in) :: n_compute, n_points, i_spin
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: density_matrix(n_compute,n_compute)
    real(dp), intent(out) :: n_phi(n_compute,n_points)
    real(dp), intent(out) :: local_rho(2,n_points)
    real(dp), intent(in), optional :: nabla_wave(n_compute,3,n_points)
    real(dp), intent(out), optional :: local_rho_gradient(3,2,n_points)
    integer :: i_dir
    call mul(density_matrix, wave, n_phi, &
         & M=n_compute, N=n_points, K=n_compute, do_serial=.true.)
    local_rho(i_spin,:) = sum(wave*n_phi,1)
    if (present(nabla_wave)) then
       do i_dir = 1, 3
          local_rho_gradient(i_dir,i_spin,:) = &
               & sum(nabla_wave(:,i_dir,:)*n_phi,1)
       end do
       local_rho_gradient(:,i_spin,:) = 2d0*local_rho_gradient(:,i_spin,:)
    end if
  end subroutine f_electron_density
  !!
  !! Matrix elements of the position vector: <phi|r|phi>
  !!
  pure subroutine f_position(n_compute, n_points, partition, wave, &
       & grid_points, matrix_wave)
    integer, intent(in) :: n_compute, n_points
    real(dp), intent(in) :: partition(n_points)
    real(dp), intent(in) :: wave(n_compute,n_points)
    real(dp), intent(in) :: grid_points(3,n_points)
    real(dp), intent(out) :: matrix_wave(n_compute,n_points,3)
    integer :: ip, i_dir
    do ip = 1, n_points
       do i_dir = 1, 3
          matrix_wave(:,ip,i_dir) = grid_points(i_dir,ip)*wave(:,ip)
       end do
       matrix_wave(:,ip,:) = partition(ip)*matrix_wave(:,ip,:)
    end do
  end subroutine f_position

  ! TESTING

  ! pure subroutine f_test(n_compute, n_points, partition, &
  !      & wave, kinetic_wave, hartree_pot, vxc, grid_points, matrix_wave)
  !   integer, intent(in) :: n_compute, n_points
  !   real(dp), intent(in) :: partition(n_points)
  !   real(dp), intent(in) :: wave(n_compute,n_points)
  !   real(dp), intent(in) :: kinetic_wave(n_compute,n_points)
  !   real(dp), intent(in) :: hartree_pot(n_points)
  !   real(dp), intent(in) :: vxc(2,n_points)
  !   real(dp), intent(in) :: grid_points(3,n_points)
  !   real(dp), intent(out) :: matrix_wave(n_compute,n_points,3,n_spin)
  !   integer :: i_dir, ip
  !   do ip = 1, n_points
  !      matrix_wave(:,ip,3,1) = &
  !           & -kinetic_wave(:,ip)*0 - (hartree_pot(ip) + 0*vxc(1,ip))*wave(:,ip)
  !      do i_dir = 1, 3
  !         matrix_wave(:,ip,i_dir,1) = &
  !              & grid_points(i_dir,ip)*matrix_wave(:,ip,3,1)
  !      end do
  !      matrix_wave(:,ip,:,1) = matrix_wave(:,ip,:,1)*partition(ip)
  !   end do
  ! end subroutine f_test

  ! END TESTING

end module integrands
