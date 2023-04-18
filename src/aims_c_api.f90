!****h* FHI-aims/aims_c_api
! PURPOSE
!   Implements C API to select aims functionality.
! AUTHOR
!   Jan Hermann
! CREATION DATE
!   2016-06-08
! COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e.V. Please note
!   that any use of the "FHI-aims-Software" is subject to the terms and
!   conditions of the respective license agreement.
!******
module aims_c_api

use iso_c_binding

implicit none

private

public :: c_evaluate_chi_0

contains

subroutine c_evaluate_chi_0(r_grid, n_grid, r_prime, u, chi_0) bind(c)
    use physics, only: &
        KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers
    use density_response, only: evaluate_chi_0

    integer(c_int), value :: n_grid
    real(c_double), intent(in) :: r_grid(n_grid, 3)
    real(c_double), intent(in) :: r_prime(3)
    complex(c_double_complex), intent(in) :: u
    real(c_double), intent(out) :: chi_0(n_grid)

    call evaluate_chi_0( &
        KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers, &
        r_grid, r_prime, u, chi_0 &
    )
end subroutine

end module
