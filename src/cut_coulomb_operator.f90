!****h* FHI-aims/cut_coulomb_operator
!  NAME
!    cut_coulomb_operator
!  SYNOPSIS

module cut_coulomb_operator

  !  PURPOSE
  !
  !    Provide a cut Coulomb operator similar to the one used in [2].  As we
  !    are using NAOs [1] and not Gaussians, simply cutting outside a given
  !    radius would be numerically unstable.  Instead we multiply the true 1/r
  !    operator with a shifted erfc() function in logarithmic space.  The
  !    width of the erfc() function should be adjusted to the density of
  !    sbtgrid points (which is in general extremely dense).
  !
  !  USAGE
  !
  !    Before the first call, this module needs to be initialized to the SBT
  !    grid to be used (see logsbt.f90 for grid convenstions).  Before the
  !    next initialization, it needs to be deallocated.
  !
  !    One way to use is to get the radial Coulomb potentials for different
  !    radial parts from cut_coulomb_logsbt_integrator().  The other is to
  !    directly make use of the Fourier transform of the cut Coulomb potential
  !    stored in the module variable cutCb_vvk.
  !
  !  USES

  implicit none

  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    [1] Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !        Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !        "Ab initio simulations with Numeric Atom-Centered Orbitals:
  !        FHI-aims", Computer Physics Communications 180, 2175 (2009).
  !    [2] Manuel Guidon, Juerg Hutter, and Joost VandeVondele,
  !        "Robust Periodic Hartree-Fock Exchange for Large-Scale
  !        Simulations Using Gaussian Basis Sets", Journal of Chemical Theory
  !        and Computation 5, 3010 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  integer :: cutCb_N
  real*8 :: cutCb_lnr0, cutCb_lnk0, cutCb_lnrange
  real*8 :: cutCb_lnwidth, cutCb_ilogcut
  real*8, allocatable :: cutCb_vvk(:)

contains

  !----------------------------------------------------------------------------
  !****s* cut_coulomb_operator/initialize_cut_coulomb
  !  NAME
  !    initialize_cut_coulomb
  !  SYNOPSIS

  subroutine initialize_cut_coulomb(N, lnr0, lnk0, lnrange, widthfac, rcut)

    !  PURPOSE
    !
    !    Initialize the cut Coulomb operator.
    !
    !  USES

    use arch_specific, only: arch_erf, arch_erfc
    use constants, only: bohr, pi
    use localorb_io, only: localorb_info, output_priority, OL_norm
    use logsbt, only: logsbt_ilog2r, logsbt_r2ilog, logsbt_multi_double_driver
    use mpi_tasks, only: check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    real*8, intent(IN) :: widthfac, rcut

    !  INPUTS
    !    o N, lnr0, lnk0, lnrange -- sbtgrid parameter (see logsbt.f90).
    !    o lnwidth -- width of cutting erfc() in logspace.
    !    o rcut -- real space cutting radius
    !  OUTPUTS
    !    none [module variables]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer, parameter :: one_zeros(1) = 0
    real*8, external :: radial_realspace_hse, radial_fourier_hse
    real*8, allocatable :: vvk_hse(:), vvk_rest(:)
    real*8 :: lnwidth
    integer :: i
    real*8 :: r, k, dln
    real*8 :: omega, damp
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'initialize_cut_coulomb'

    lnwidth = log(widthfac)
    dln = lnrange / N

    if(output_priority .le. OL_norm ) then
       write(info_str, "(2X,A)") 'Initializing cut Coulomb operator.'
       call localorb_info(info_str)
       write(info_str, "(2X,'| ',A,F10.4,A,F10.4,A,I5,A)") &
       & 'Cut [erfc(-/+1)/2] from', rcut/widthfac*bohr, ' A to', &
       & rcut*widthfac*bohr, ' A within', &
       & floor(2*lnwidth/dln), ' sbt_grid points.'
       call localorb_info(info_str)
       call localorb_info('')
    endif

!!$    write(info_str, "(2X,'| ',A,F10.4,A)") &
!!$    & 'Initializing cut Coulomb operator to range', rcut*bohr, ' A'
!!$    call localorb_info(info_str)
!!$    write(info_str, "(2X,'| ',A,F10.4,A,F10.4,A)") &
!!$    & 'and a 0.5*erfc(-/+1) smearing between ', &
!!$    & rcut/widthfac*bohr, ' A and', rcut*widthfac*bohr, ' A.'
!!$    call localorb_info(info_str)
!!$    write(info_str, "(2X,'| ',A,I6,A)") &
!!$    & 'Within this region, there are at least', floor(2*lnwidth/dln), &
!!$    & ' sbt_grid points.'
!!$    call localorb_info(info_str)

    cutCb_N = N
    cutCb_lnr0 = lnr0
    cutCb_lnk0 = lnk0
    cutCb_lnrange = lnrange
    cutCb_lnwidth = lnwidth
    cutCb_ilogcut = logsbt_r2ilog(rcut, lnr0, dln)
    omega = 5.d0 / rcut    ! Do not extent range
    allocate(vvk_hse(N), vvk_rest(N), stat=info)
    call check_allocation(info, 'vvk_hse, vvk_rest', func)

    ! --- Short range part directly in k-space

    do i = 1, N
       k = logsbt_ilog2r(dble(i), lnk0, dln)
       vvk_hse(i) = radial_fourier_hse(k, omega)
    end do

    ! --- Long range correction in r-space

    do i = 1, N
       r = logsbt_ilog2r(dble(i), lnr0, dln)
       damp = 0.5d0 * arch_erfc(dln*(dble(i) - cutCb_ilogcut) / cutCb_lnwidth)
       ! 1 <= damp <= 0
       vvk_rest(i) = damp * arch_erf(omega * r) / r
    end do

    call logsbt_multi_double_driver(N, lnr0, lnk0, lnrange, &
    &                               1, one_zeros, vvk_rest)

    ! --- Put together

    if (allocated(cutCb_vvk)) deallocate(cutCb_vvk)
    allocate(cutCb_vvk(N), stat=info)
    call check_allocation(info, 'cutCb_vvk', func)

    cutCb_vvk = (2*pi)**(-1.5d0) * vvk_hse + vvk_rest

    deallocate(vvk_hse, vvk_rest)

  end subroutine initialize_cut_coulomb
  !******
  !----------------------------------------------------------------------------
  !****s* cut_coulomb_operator/deallocate_cut_coulomb
  !  NAME
  !    deallocate_cut_coulomb
  !  SYNOPSIS

  subroutine deallocate_cut_coulomb()

    !  PURPOSE
    !    Tidy up module data
    !  USES

    implicit none

    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    character(*), parameter :: func = 'deallocate_cut_coulomb'

    cutCb_N = 0
    cutCb_lnr0 = 0.d0
    cutCb_lnk0 = 0.d0
    cutCb_lnrange = 0.d0
    cutCb_lnwidth = 0.d0
    cutCb_ilogcut = 0.d0
    if (allocated(cutCb_vvk)) deallocate(cutCb_vvk)

  end subroutine deallocate_cut_coulomb
  !******
  !----------------------------------------------------------------------------
  !****s* cut_coulomb_operator/cut_coulomb_logsbt_integrator
  !  NAME
  !    cut_coulomb_logsbt_integrator
  !  SYNOPSIS

  subroutine cut_coulomb_logsbt_integrator(N, n_fn, ffr, fn_to_L, &
  &                                        lnr0, lnk0, lnrange)

    !  PURPOSE
    !
    !    Calculate the cut Coulomb field using logSBT.
    !
    !  USES

    use constants, only: pi
    use logsbt, only: logsbt_multi_double_driver
    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    integer, intent(IN) :: n_fn
    real*8, intent(INOUT) :: ffr(N, n_fn)
    integer, intent(IN) :: fn_to_L(n_fn)
    real*8, intent(IN) :: lnr0, lnk0, lnrange

    !  INPUTS
    !    o N, lnr0, lnk0, lnrange -- sbtgrid parameters (see logsbt.f90 header)
    !    o n_fn -- Number of radial parts
    !    o fn_to_L -- Angular momenta
    !    o ffr -- Radial part f(r)
    !  OUTPUTS
    !    o ffr -- Radial part of potential v(r)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: i
    real*8 :: fac
    character(*), parameter :: func = 'cut_coulomb_logsbt_integrator'

    if (N /= cutCb_N) call aims_stop('N mismatch', func)
    if (abs(lnr0 - cutCb_lnr0) > 1d-10) call aims_stop('lnr0 mismatch', func)
    if (abs(lnk0 - cutCb_lnk0) > 1d-10) call aims_stop('lnk0 mismatch', func)
    if (abs(lnrange - cutCb_lnrange) > 1d-10) then
       call aims_stop('lnrange mismatch', func)
    end if

    ! --- f(r) -> f~(k)

    call logsbt_multi_double_driver(N, lnr0, lnk0, lnrange, &
    &                               n_fn, fn_to_L, ffr)

    ! --- f~(k) -> V~(k)

    fac = (2*pi)**1.5d0
    do i = 1, N
       ffr(i,:) = ffr(i,:) * fac * cutCb_vvk(i)
    end do

    ! --- V~(k) -> V(r)

    call logsbt_multi_double_driver(N, lnk0, lnr0, lnrange, &
    &                               n_fn, fn_to_L, ffr)

  end subroutine cut_coulomb_logsbt_integrator
  !******
  !----------------------------------------------------------------------------
  !****s* cut_coulomb_operator/cut_coulomb_analytic
  !  NAME
  !    cut_coulomb_analytic
  !  SYNOPSIS

  subroutine initialize_cut_coulomb_analytic(N, lnr0, lnk0, lnrange, widthfac, rcut)

    !  PURPOSE
    !
    !    Initialize the cut Coulomb operator according to 
    !    J. Spencer, A. Alavi, PRB,77 193110 (2008)
    !
    !  USES

    use arch_specific, only: arch_erf, arch_erfc
    use constants, only: bohr, pi
    use localorb_io, only: localorb_info, output_priority, OL_norm
    use logsbt, only: logsbt_ilog2r, logsbt_r2ilog
    use mpi_tasks, only: check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    real*8, intent(IN) :: widthfac, rcut

    !  INPUTS
    !    o N, lnr0, lnk0, lnrange -- sbtgrid parameter (see logsbt.f90).
    !    o lnwidth -- width of cutting erfc() in logspace.
    !    o rcut -- real space cutting radius
    !  OUTPUTS
    !    none [module variables]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer, parameter :: one_zeros(1) = 0
    !real*8, external :: radial_realspace_hse, radial_fourier_hse
    !real*8, allocatable :: vvk_hse(:), vvk_rest(:)
    real*8 :: lnwidth
    integer :: i
    real*8 :: r, k, dln
    real*8 :: omega, damp
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'initialize_cut_coulomb_analytic'

    lnwidth = log(widthfac)
    dln = lnrange / N

    if(output_priority .le. OL_norm ) then
       write(info_str, "(2X,A)") 'Initializing analytic cut Coulomb operator.'
       call localorb_info(info_str)
       write(info_str, "(2X,'| ',A,F10.4,A)") &
       & 'Cut Coulomb start from', rcut*bohr, ' A'
       call localorb_info(info_str)
       call localorb_info('')
    endif

    cutCb_N = N
    cutCb_lnr0 = lnr0
    cutCb_lnk0 = lnk0
    cutCb_lnrange = lnrange
    cutCb_lnwidth = lnwidth
    cutCb_ilogcut = logsbt_r2ilog(rcut, lnr0, dln)

    if (allocated(cutCb_vvk)) deallocate(cutCb_vvk)
    allocate(cutCb_vvk(N), stat=info)
    call check_allocation(info, 'cutCb_vvk', func)

    do i = 1, N
       k = logsbt_ilog2r(dble(i), lnk0, dln)
       cutCb_vvk(i) = (2*pi)**(-1.5d0) * 4*pi / k**2 * &
       (1.d0 - cos(abs( k * rcut)))
    end do

  end subroutine initialize_cut_coulomb_analytic
  !----------------------------------------------------------------------------
  !****s* cut_coulomb_operator/initialize_ers
  !  NAME
  !    initialize_cut_coulomb
  !  SYNOPSIS

  subroutine initialize_ers(N, lnr0, lnk0, lnrange, widthfac, rcut, omega)

    !  PURPOSE
    !
    !    Initialize the cut Coulomb operator.
    !
    !  USES

    use arch_specific, only: arch_erf, arch_erfc
    use constants, only: pi
    use localorb_io, only: localorb_info, output_priority, OL_norm
    use logsbt, only: logsbt_ilog2r, logsbt_r2ilog, logsbt_multi_double_driver
    use mpi_tasks, only: check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    real*8, intent(IN) :: widthfac, rcut

    !  INPUTS
    !    o N, lnr0, lnk0, lnrange -- sbtgrid parameter (see logsbt.f90).
    !    o widthfac -- width of cutting erfc().
    !    o rcut -- real space cutting radius
    !  OUTPUTS
    !    none [module variables]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer, parameter :: one_zeros(1) = 0
    real*8, external :: radial_realspace_hse, radial_fourier_hse
    real*8, allocatable :: vvk_hse(:), vvk_rest(:)
    real*8 :: lnwidth
    integer :: i
    real*8 :: r, k, dln
    real*8 :: omega, damp
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'initialize_cut_coulomb'

    lnwidth = log(widthfac)
    dln = lnrange / N

    if(output_priority .le. OL_norm ) then
       write(info_str, "(2X,A)") 'Initializing cut Coulomb operator.'
       call localorb_info(info_str)
       !write(info_str, "(2X,'| ',A,F10.4,A,F10.4,A,I5,A)") &
       !& 'Cut [erfc(-/+1)/2] from', rcut/widthfac*bohr, ' A to', &
       !& rcut*widthfac*bohr, ' A within', &
       !& floor(2*lnwidth/dln), ' sbt_grid points.'
       !call localorb_info(info_str)
       !call localorb_info('')
    endif

    cutCb_N = N
    cutCb_lnr0 = lnr0
    cutCb_lnk0 = lnk0
    cutCb_lnrange = lnrange
    cutCb_lnwidth = lnwidth
    cutCb_ilogcut = logsbt_r2ilog(rcut, lnr0, dln)
    !omega = 5.d0 / rcut    ! Do not extent range
    allocate(vvk_hse(N), vvk_rest(N), stat=info)
    call check_allocation(info, 'vvk_hse, vvk_rest', func)

    ! --- Short range part directly in k-space

    do i = 1, N
       k = logsbt_ilog2r(dble(i), lnk0, dln)
       vvk_hse(i) = radial_fourier_hse(k, omega)
    end do

    ! --- Long range correction in r-space

    do i = 1, N
       r = logsbt_ilog2r(dble(i), lnr0, dln)
       !damp = 0.5d0 * arch_erfc(dln*(dble(i) - cutCb_ilogcut) / cutCb_lnwidth)
       damp = 0.5d0 * arch_erfc((r - rcut) * widthfac)
       ! 1 <= damp <= 0
       vvk_rest(i) = damp * arch_erf(omega * r) / r
    end do

    call logsbt_multi_double_driver(N, lnr0, lnk0, lnrange, &
    &                               1, one_zeros, vvk_rest)

    ! --- Put together

    if (allocated(cutCb_vvk)) deallocate(cutCb_vvk)
    allocate(cutCb_vvk(N), stat=info)
    call check_allocation(info, 'cutCb_vvk', func)

    cutCb_vvk = (2*pi)**(-1.5d0) * vvk_hse + vvk_rest

    deallocate(vvk_hse, vvk_rest)

  end subroutine initialize_ers
  !******
end module cut_coulomb_operator
!******
