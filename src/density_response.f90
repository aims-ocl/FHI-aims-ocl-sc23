!****h* FHI-aims/density_response
! PURPOSE
!   Collection of functions to study the density response function
! AUTHOR
!   Jan Hermann
! CREATION DATE
!   2016-06-08
! COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e.V. Please note
!   that any use of the "FHI-aims-Software" is subject to the terms and
!   conditions of the respective license agreement.
!******
module density_response

implicit none

private

public :: evaluate_chi_0

integer :: dummy_int

external :: &
    prune_basis_p0, prune_radial_basis_p0, tab_local_geometry_p0, &
    tab_trigonom_p0, tab_wave_ylm_p0, evaluate_radial_functions_p0, &
    evaluate_waves_p0

contains

subroutine evaluate_chi_0(c, c_cmplx, omega, f, r_grid, r_prime, u, chi_0)
    ! Expressing the Adler–Wiser formula in a form suitable for computation. ("$*​$" is component-wise multiplication.)
    ! $$
    ! \begin{aligned}
    ! \chi_0(\mathbf r,\mathbf r',\mathrm iu)
    ! &=\sum_{ab}\frac{f_a-f_b}{\omega_a-\omega_b+\mathrm iu}
    ! \psi_a^*(\mathbf r)\psi_b^*(\mathbf r)\psi_a(\mathbf r')\psi_b(\mathbf r')+\text{c.c.}\\
    ! &=\sum_{ab}\frac{f_a-f_b}{\omega_a-\omega_b+\mathrm iu}
    ! \left[\sum_{ij}c^*_{ia}c^*_{jb}\varphi_i(\mathbf r)\varphi_j(\mathbf r)\right]\!\!
    ! \left[\sum_{kl}c_{ka}c_{lb}\varphi_k(\mathbf r')\varphi_l(\mathbf r')\right]+\text{c.c.}\\
    ! &=\sum_{ij}\left[\sum_{ab}\underbrace{\frac{f_a-f_b}{\omega_a-\omega_b+\mathrm iu}}_\mathbf G
    ! c^*_{ia}c^*_{jb}\sum_{kl}c_{ka}c_{lb}\underbrace{\varphi_k(\mathbf r')\varphi_l(\mathbf r')}_{\mathbf\Phi(\mathbf r')}+\text{c.c.}\right]
    ! \varphi_i(\mathbf r)\varphi_j(\mathbf r)\\
    ! &=\boldsymbol\varphi^\mathrm T(\mathbf r)\underbrace{\left[
    ! \mathbf c^*\big(\mathbf G*\mathbf c^\mathrm T\mathbf\Phi(\mathbf r')\mathbf c\big)\mathbf c^\dagger+\text{c.c.}
    ! \right]}_{\mathbf Q(\mathbf r')}\boldsymbol\varphi(\mathbf r)
    ! \end{aligned}
    ! $$
    ! Pseudocode.
    ! $$
    ! \begin{aligned}
    ! &\text{with given $\mathbf r'$, $u$:}\\
    ! &G_{ab}:=\frac{f_a-f_b}{\omega_a-\omega_b+\mathrm iu}\\
    ! &\Phi_{kl}(\mathbf r'):=\varphi_k(\mathbf r')\varphi_l(\mathbf r')\\
    ! &\mathbf Q(\mathbf r'):=\mathbf c^*\big(\mathbf G*\mathbf c^\mathrm T\mathbf\Phi(\mathbf r')\mathbf c\big)\mathbf c^\dagger+\text{c.c.}\\
    ! &\text{for each }\mathbf r:\\
    ! &\quad\chi_0(\mathbf r,\mathbf r',\mathrm iu):=\boldsymbol\varphi^\mathrm T(\mathbf r)\mathbf Q(\mathbf r')\boldsymbol\varphi(\mathbf r)
    ! \end{aligned}
    ! $$
    use dimensions, only: &
        n_basis, n_states, n_spin, n_k_points_task, n_k_points, &
        n_max_compute_ham, n_cbasis
    use pbc_lists, only: &
        cbasis_to_basis

    real(8), intent(in) :: c(n_basis, n_states, n_spin, n_k_points_task)
    complex(8), intent(in) :: c_cmplx(n_basis, n_states, n_spin, n_k_points_task)
    real(8), intent(in) :: omega(n_states, n_spin, n_k_points)
    real(8), intent(inout) :: f(n_states, n_spin, n_k_points)
    real(8), intent(in) :: r_grid(:, :)
    real(8), intent(in) :: r_prime(3)
    complex(8), intent(in) :: u
    real(8), intent(out) :: chi_0(size(r_grid, 1))

    integer :: i_basis(n_cbasis)
    real(8) :: varphi(n_max_compute_ham)
    integer :: n_compute, i_point
    complex(8), dimension(n_states, n_states) :: G
    real(8), dimension(n_basis, n_basis) :: Phi, Q
    integer :: a, b, i, j, k, l, ip, jp, kp, lp

    G(:, :) = 0.d0
    do a = 1, n_states
        do b = 1, n_states
            if (abs(f(a, 1, 1)-f(b, 1, 1)) < 1d-10) cycle
            G(a, b) = (f(a, 1, 1)-f(b, 1, 1)) &
                / (omega(a, 1, 1)-omega(b, 1, 1)+(0.d0, 1.d0)*u)
        end do
    end do
    call evaluate_varphi(r_prime, varphi, i_basis, n_compute)
    Phi(:, :) = 0.d0
    do kp = 1, n_compute
        k = cbasis_to_basis(i_basis(kp))
        do lp = 1, n_compute
            l = cbasis_to_basis(i_basis(lp))
            Phi(k, l) = varphi(kp)*varphi(lp)
        end do
    end do
    Q = 2*real(matmul( &
        c(:, :, 1, 1), matmul( &
        G*matmul( &
            transpose(c(:, :, 1, 1)), matmul( &
            Phi, &
            c(:, :, 1, 1))), &
        transpose(c(:, :, 1, 1))) &
    ))
    chi_0(:) = 0.d0
    do i_point = 1, size(r_grid, 1)
        call evaluate_varphi(r_grid(i_point, :), varphi, i_basis, n_compute)
        do ip = 1, n_compute
            i = cbasis_to_basis(i_basis(ip))
            do jp = 1, n_compute
                j = cbasis_to_basis(i_basis(jp))
                chi_0(i_point) = chi_0(i_point) &
                    + Q(i, j)*varphi(ip)*varphi(jp)
            end do
        end do
    end do
end subroutine

subroutine evaluate_varphi(r, varphi, i_basis, n_compute)
    use dimensions, only: &
        n_centers_integrals, &
        n_centers, n_basis_fns, l_wave_max, &
        n_max_compute_atoms, n_cbasis, n_max_compute_fns_ham, &
        n_max_compute_ham
    use pbc_lists, only: &
        centers_basis_integrals, coords_center, &
        inv_centers_basis_integrals
    use basis, only: basis_wave_ordered
    use species_data, only: l_shell_max

    real(8), intent(in) :: r(3)
    real(8), intent(out) :: varphi(n_max_compute_ham)
    integer, intent(out) :: i_basis(n_cbasis)
    integer, intent(out) :: n_compute

    integer :: &
        i_center, i_center_integrals, &
        n_compute_fns, n_compute_atoms, &
        i_index, i_l, i_m
    integer :: &
        i_basis_fns(n_basis_fns, n_centers_integrals), &
        i_basis_fns_inv(n_basis_fns, n_centers), &
        i_atom_fns(n_basis_fns, n_centers_integrals), &
        atom_index(n_centers_integrals), &
        atom_index_inv(n_centers), &
        spline_array_start(n_centers_integrals), &
        spline_array_end(n_centers_integrals), &
        index_lm(-l_wave_max:l_wave_max, 0:l_wave_max)
    real(8) :: &
        r_diff(3), &
        dist_tab(n_centers_integrals), &
        dist_tab_sq(n_centers_integrals), &
        i_r(n_max_compute_atoms), &
        dir_tab(3, n_centers_integrals), &
        trigonom_tab(4, n_max_compute_atoms), &
        radial_wave(n_max_compute_fns_ham), &
        ylm_tab((l_wave_max+1)**2, n_max_compute_atoms)

    i_index = 0
    do i_l = 0, l_wave_max
        do i_m = -i_l, i_l
            i_index = i_index+1
            index_lm(i_m, i_l) = i_index
        end do
    end do
    n_compute = 0
    n_compute_atoms = 0
    n_compute_fns = 0
    varphi(:) = 0.d0
    do i_center_integrals = 1, n_centers_integrals
        i_center = centers_basis_integrals(i_center_integrals)
        r_diff(:) = r(:)-coords_center(:, i_center)
        dir_tab(:, i_center_integrals) = r_diff(:)
        dist_tab_sq(i_center_integrals) = sum(r_diff**2)
    enddo
    call prune_basis_p0( &
        dist_tab_sq, dummy_int, n_compute, i_basis, &
        n_cbasis, n_centers_integrals, inv_centers_basis_integrals &
    )
    call prune_radial_basis_p0( &
        dist_tab_sq, dist_tab, dir_tab, &
        n_compute_atoms, atom_index, atom_index_inv, n_compute_fns, &
        i_basis_fns, i_basis_fns_inv, i_atom_fns, spline_array_start, &
        spline_array_end, n_centers_integrals, centers_basis_integrals &
    )
    call tab_local_geometry_p0( &
        dist_tab_sq, n_compute_atoms, atom_index, dir_tab, dist_tab, i_r &
    )
    call tab_trigonom_p0(n_compute_atoms, dir_tab, trigonom_tab)
    call tab_wave_ylm_p0( &
        n_compute_atoms, atom_index, trigonom_tab, l_shell_max, &
        l_wave_max, ylm_tab &
    )
    call evaluate_radial_functions_p0( &
        spline_array_start, spline_array_end, n_compute_atoms, &
        n_compute_fns, dist_tab, i_r, atom_index, i_basis_fns_inv, &
        basis_wave_ordered, radial_wave, .false., n_compute, &
        n_max_compute_fns_ham &
    )
    call evaluate_waves_p0( &
        l_wave_max, ylm_tab, dist_tab, index_lm, n_compute, i_basis, &
        radial_wave, varphi, n_compute_atoms, atom_index_inv, &
        n_compute_fns, i_basis_fns_inv, n_max_compute_fns_ham &
    )
end subroutine

end module
