!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Test the radial extent of each Gaussian orbital and if any of them
!!  are larger than free_r_cut then set free_r_cut to that value. The
!!  radial extent is defined as the distance at which the radial
!!  function becomes less than wave_threshold.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
!!  COMMENTS
!!
!!  Here by radial function we mean u(r) and not u(r)/r.
!!
subroutine set_gaussian_free_cut(i_species, wave_threshold, free_r_cut)

  use constants,    only: bohr, pi
  use grids,        only: r_grid, n_grid, get_logarithmic_grid
  use localorb_io,  only: localorb_multi
  use mpi_tasks,    only: aims_stop
  use species_data, only: gaussian_alpha, gaussian_coeff, gaussian_n, &
       & gaussian_n_contr, n_gaussian, species_name
  use types,        only: dp

  implicit none

  interface
     real(dp) elemental function factorial(n)
       import dp
       integer, intent(in) :: n
     end function factorial
  end interface

  integer, intent(in) :: i_species
  real(dp), intent(in) :: wave_threshold
  real(dp), intent(in out) :: free_r_cut
  real(dp), allocatable :: coeffs(:), alphas(:), tmp(:)
  real(dp) :: total_norm
  real(dp) :: r_extents(n_gaussian(i_species))
  character(60) :: info_str
  integer :: i_gaussian, l_channel, ig, i_radial(1)

  ! The logarithmic grids have not been set yet. This subroutine will
  ! be called again later, which is fine.
  call get_logarithmic_grid(i_species)

  do i_gaussian = 1, n_gaussian(i_species)
     ! STEP 1 - Build the gaussian orbital
     l_channel = gaussian_n(i_species,i_gaussian)
     allocate(alphas(gaussian_n_contr(i_species, i_gaussian)))
     allocate(coeffs(gaussian_n_contr(i_species, i_gaussian)))
     alphas = gaussian_alpha(i_species, i_gaussian, :size(alphas))
     coeffs = gaussian_coeff(i_species, i_gaussian, :size(coeffs))
     ! Normalize the coefficients by the norms of the elemental
     ! Gaussian functions
     coeffs = coeffs/sqrt(norm(alphas))
     ! Compute the norm of the complete orbital
     total_norm = sum(coeffs**2*norm(alphas))
     do ig = 1, size(alphas)
        total_norm = total_norm + sum(2*coeffs(ig)*coeffs(ig+1:)* &
             & norm((alphas(ig)+alphas(ig+1:))/2))
     end do
     ! STEP 2 - Test wave_threshold
     ! First determine the position of the first maximum
     i_radial = maxloc(f_gauss(r_grid(:n_grid(i_species), i_species)))
     allocate(tmp(n_grid(i_species)))
     ! Find where the radial function becomes less than wave_threshold
     ! only for points that come after the first maximum
     tmp = f_gauss(r_grid(:n_grid(i_species), i_species)) - wave_threshold
     i_radial = i_radial-1 + minloc(tmp(i_radial(1):), tmp(i_radial(1):) > 0d0)
     ! The position just determined becomes the radial extent
     r_extents(i_gaussian) = r_grid(i_radial(1),i_species)
     write(info_str, '(a, i2, a, i1, a, f9.6, a)') 'Species: '// &
          & trim(species_name(i_species))//', Function: ', i_gaussian, &
          & ', L = ', l_channel, ', extent = ', r_extents(i_gaussian)*bohr, &
          & ' AA'
     call localorb_multi(info_str, format='(4x, a)')
     deallocate(alphas, coeffs, tmp)
  end do
  ! STEP 3 - Update free_r_cut if any radial extent is greater than
  !          its current value.
  if (maxval(r_extents) > free_r_cut) then
     free_r_cut = maxval(r_extents)
     write(info_str, '(f9.6)') free_r_cut*bohr
     call localorb_multi( &
          & 'Species: '//trim(species_name(i_species))//', new cutoff onset &
          &for free atom density:', '                             free_r_cut &
          &= '//trim(info_str)//' AA', format='(4x, a)')
  end if

contains
  ! Norm of an elemental Gaussian orbital.
  elemental real(dp) function norm(alpha) result(y)
    real(dp), intent(in) :: alpha
    y = factorial(2*l_channel+2)/(factorial(l_channel+1)* &
         & 2d0**(2*l_channel+3))*sqrt(pi/((2*alpha)**(2*l_channel+3)))
  end function norm

  ! Evaluates the value of a complete Gaussian orbital
  elemental real(dp) function f_gauss(x) result(y)
    real(dp), intent(in) :: x
    y = sum(coeffs*x**(l_channel+1)*exp(-alphas*x**2))
    y = y/sqrt(total_norm)
  end function f_gauss
end subroutine set_gaussian_free_cut
