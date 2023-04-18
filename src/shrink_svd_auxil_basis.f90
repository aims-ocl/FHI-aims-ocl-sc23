!****s* FHI-aims/shrink_svd_auxil_basis
!  NAME
!    shrink_svd_auxil_basis
!  SYNOPSIS

subroutine shrink_svd_auxil_basis()

  !  PURPOSE
  !
  !   Generates one-particle basis function onsite products and separates
  !   these products into different angular momentum channels.  For each
  !   L-channel, perform a singular value decomposition with respect to the
  !   Coulomb metric (well, actually a diagonalization) to generate the
  !   product basis (called "prodbas" or "basbas").
  !
  !   The idea is to minimize the sum over the squared norm (Coulomb
  !   self-repulsion) of all residues of pairs of basis functions:
  !
  !     I = \sum_ab ( \delta\rho_ab | \delta\rho_ab )
  !
  !   with
  !
  !      (a|b) = \int_R^3  [ a(vecr) b(vecr') ] / |vecr - vecr'|
  !
  !      \delta\rho_ab(vecr) = f_a * f_b - \sum_\mu C^ab_\mu P_\mu
  !                          =  \rho_ab   - \tilde\rho_ab.
  !
  !   As log as the coefficients C^ab_\mu are choosen optimally, this is
  !   equivalent to maximizing
  !
  !     J = \sum_ab ( \tilde\rho_ab | \tilde\rho_ab ).
  !
  !   This way, product basis functions are kept according to their benefit
  !   for describing onsite products.  In principle, this procedure can be
  !   straightforwardly generalized to also take into account some offsite
  !   products.
  !
  !  USES

  use debug_output
  use species_data


  use dimensions
  use runtime_choices
  use grids
  use geometry
  use basis
  use species_data
  use spline
  use prodbas
  use basbas_fn_coulomb, only: localize_basbas_fn_coulomb, &
      sort_species_basbas_fn
  use mpi_tasks
  use constants
  use localorb_io
  use numerical_utilities
  use triple_Y
  implicit none

  !  ARGUMENTS
  !    none
  !  INPUTS
  !    none [Things from basis.f90, etc.]
  !  OUTPUTS
  !    none [Things in prodbas.f90, etc.]
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
  !    Release version, FHI-aims (2010).
  !  SOURCE

  logical :: use_coulomb_metric, use_hse_localization
  integer :: l_prod_max, n_max_prod_fns, n_prods, i_prod_1, i_prod_2
  integer :: n_basbasfn_uptonow, i_basbasfn_start_L, i_basbasfn_start_species
  integer :: n_basbasfn_species, n_reject_uptonow, n_basbas_uptonow
  integer :: i_species, L, i_fn_1, i_fn_2, l1, l2, i_prodbas_fn, i_prod
  integer :: i_basbas_fn
  integer, allocatable :: n_prods_uptonow(:)
  real*8 :: onsite_fac, Y_norm
  real*8, allocatable :: tripY(:,:)
  real*8, allocatable :: metric(:,:), eigens(:,:), diags(:)
  real*8, allocatable :: prod_wave(:,:,:)
  real*8, allocatable :: basbasfn_wave(:,:)
  integer, allocatable :: fn2species(:)
  integer, allocatable :: fn2L(:)

  real*8 :: int_log_mesh, int_log_coulomb_metric

  integer :: info
  character*150 :: info_str
  character(*), parameter :: func = 'shrink_svd_auxil_basis'

  ! --- Initialize

  call localorb_info('')
  call localorb_info("--------------------------------------------")
  call localorb_info("Constructing auxiliary basis (full product) ...")
  call localorb_info('')

  use_coulomb_metric = .not. (RI_type == RI_SVS)
  l_prod_max = 2*l_wave_max
  n_max_prod_fns = n_max_basis_fns**2
  use_hse_localization = (use_hse .and. hse_omega_hf /= 0.d0 .and. .not. use_gw_and_hse &
      .and. .not. use_dftpt2_and_hse)

  ! --- Per species buffers

  allocate(tripY(-l_prod_max:l_prod_max, -l_prod_max:l_prod_max), stat=info)
  call check_allocation(info, 'tripY', func)
  allocate(prod_wave(n_max_grid, n_max_prod_fns, 0:l_prod_max), stat=info)
  call check_allocation(info, 'prod_wave', func)
  allocate(n_prods_uptonow(0:l_prod_max), stat=info)
  call check_allocation(info, 'n_prods_uptonow', func)

  ! --- Global buffers (input for analyze_basbas).

  allocate(basbasfn_wave(n_max_grid, n_max_prod_fns), stat=info)
  call check_allocation(info, 'basbasfn_wave')
  allocate(fn2species(n_max_prod_fns), stat=info)
  call check_allocation(info, 'fn2species')
  allocate(fn2L(n_max_prod_fns), stat=info)
  call check_allocation(info, 'fn2L')
  n_basbasfn_uptonow = 0
  n_reject_uptonow = 0
  
  ! --- Main species loop

  SPECIES_L: do i_species = 1, n_species
     
     ! --- Generate products

     ! F(rvec) = f1(rvec) * f2(rvec) = f1(r) * Y_l1m1 * f2(r) * Y_l2m2
     !         = f1(r) * f2(r) * \sum_L tripY(m1, m2) * Y_L(m1+m2)

     ! Loop over all pairs of onsite functions
     n_prods_uptonow = 0
     FN_LOOP: do i_fn_1 = 1, n_basis_fns
        if (basisfn_species(i_fn_1) /= i_species) cycle
        l1 = basisfn_l(i_fn_1)
        do i_fn_2 = 1, i_fn_1
           if (basisfn_species(i_fn_2) /= i_species) cycle
           l2 = basisfn_l(i_fn_2)

           if (i_fn_1 == i_fn_2) then
              onsite_fac = 1.d0
           else
              ! There are two, and we take only one.  Should make things
              ! invariant under unitary transformation of the one-particle
              ! basis.
              onsite_fac = sqrt(2.d0)
           end if


           ! Loop over all symmetry-allowed product angular momenta:
           !   l1=2, l2=3   ->   L = 1, 3, 5
           do L = abs(l1 - l2), min(l1 + l2, max_l_prodbas(i_species)), 2
              ! Get coefficients in Y_l1m1 * Y_l2m2 = \sum_L C_L Y_L(m1+m2).
              call triple_Y_cmplx(l1, l2, L, tripY(-l1:l1, -l2:l2))
              Y_norm = sum(tripY(-l1:l1, -l2:l2)**2)

              n_prods_uptonow(L) = n_prods_uptonow(L) + 1
              prod_wave(:, n_prods_uptonow(L), L) = &
              & basis_wave_spl( 1, :, i_fn_1) * &
              & basis_wave_spl( 1, :, i_fn_2) * &
              & Y_norm * onsite_fac / r_grid(:, i_species) / (2*L+1)
              ! Divide by r_grid because all quantities are u(r)=r*f(r).
              ! Divide by (2*L+1) because what costs is the number of
              ! product basis functions, not the number of radial parts.
              ! Say it that way: We throw out a product basis functions
              ! if it contributes less than prodbas_acc.
           end do
        end do
     end do FN_LOOP
     
     ! --- Diagonalization of metric

     i_basbasfn_start_species = n_basbasfn_uptonow + 1
     ADD_L_LOOP: do L = 0, l_prod_max
        n_prods = n_prods_uptonow(L)
        if (n_prods == 0) cycle
        allocate(metric(n_prods, n_prods), stat=info)
        call check_allocation(info, 'metric', func)
        allocate(eigens(n_prods, n_prods), diags(n_prods), stat=info)
        call check_allocation(info, 'eigens, diags', func)

        do i_prod_1 = 1, n_prods
           do i_prod_2 = 1, n_prods
              if (use_coulomb_metric) then
                 metric(i_prod_1, i_prod_2) = &
                 & int_log_coulomb_metric(L, &
                 &    prod_wave(:, i_prod_1, L), prod_wave(:, i_prod_2, L), &
                 &    n_grid(i_species), r_grid(:, i_species))
              else
                 metric(i_prod_1, i_prod_2) = &
                 & int_log_mesh(L, &
                 &    prod_wave(:, i_prod_1, L), prod_wave(:, i_prod_2, L), &
                 &    n_grid(i_species), r_grid(:, i_species))
              end if
           end do
        end do
        eigens = metric
        call diagonalize_rmatrix(n_prods, eigens, diags, .true.)
        ! metric * eigens = eigens * diags

        ! Keep all linear combinations with an eigenvalue larger than the
        ! threshold.  Do *not* normalize them, so that the global threshold
        ! cuts according to the weights from the onsite products.
        i_basbasfn_start_L = n_basbasfn_uptonow + 1
        do i_prodbas_fn = 1, n_prods
           ONSITE_THRES: if (diags(i_prodbas_fn) > prodbas_acc(i_species)) then
              n_basbasfn_uptonow = n_basbasfn_uptonow + 1
              fn2species(n_basbasfn_uptonow) = i_species
              fn2L(n_basbasfn_uptonow) = L
              basbasfn_wave(:, n_basbasfn_uptonow) = 0.d0
              do i_prod = 1, n_prods
                 basbasfn_wave(:, n_basbasfn_uptonow) = &
                 & basbasfn_wave(:, n_basbasfn_uptonow) + &
                 & prod_wave(:, i_prod, L) * eigens(i_prod, i_prodbas_fn)
              end do
           else
              n_reject_uptonow = n_reject_uptonow + 1
           end if ONSITE_THRES
        end do

        ! Unitary transform (take linear combinations) the radial parts
        ! to ensure that only one takes a finite multipole moment.
        if (.not. (use_hse .and. hse_omega_hf /= 0.d0).or.use_gw_and_hse&
            .or.use_dftpt2_and_hse) then
           call localize_basbas_fn_coulomb(L, &
           & i_basbasfn_start_L, n_basbasfn_uptonow, basbasfn_wave, &
           & n_grid(i_species), r_grid(:, i_species), &
           & Adams_Moulton_integrator, use_hse_localization)
        end if

        deallocate(metric, eigens, diags)

     end do ADD_L_LOOP

     ! Sort basbas_fn within one species by extent.
     if (.not. (use_hse .and. hse_omega_hf /= 0.d0) .or. use_gw_and_hse &
         .or. use_dftpt2_and_hse) then
        n_basbasfn_species = n_basbasfn_uptonow - i_basbasfn_start_species + 1
        call sort_species_basbas_fn(n_basbasfn_species, &
        & n_grid(i_species), r_grid(:, i_species), wave_threshold, &
        & basbasfn_wave(:, i_basbasfn_start_species:n_basbasfn_uptonow), &
        & fn2L(i_basbasfn_start_species:n_basbasfn_uptonow))
     end if

  end do SPECIES_L
  
  ! --- Set dimensions

  n_basbas_fns = n_basbasfn_uptonow
  n_basbas_uptonow = 0
  do i_basbas_fn = 1, n_basbas_fns
     i_species = fn2species(i_basbas_fn)
     L = fn2L(i_basbas_fn)
     n_basbas_uptonow = n_basbas_uptonow + &
     &                  (2*L+1) * atoms_in_structure(i_species)
  enddo
  n_basbas = n_basbas_uptonow
  n_basbas_supercell = n_basbas

  if(myid.eq.0) then
     write(use_unit, '(2X, A,1X, A, I8, 1X, A, I8, 1X, A)') &
     "| Shrink_svd_auxil_basis :", "there are ", n_basbas_fns, &
     " radial auxiliary wave functions"
     write(use_unit, '(27X,A,I8,A)') &
     " accepted and", n_reject_uptonow, " rejected."
     write(use_unit, '(2X,A,I8,A)') &
     "| Shrink_svd_auxil_basis : there are totally", n_basbas, &
     " partial auxiliary wave functions."
  endif

  ! --- Allocate and fill product basis arrays arrays

  call get_bas_dimensions(n_basbas_fns, &
  &                       fn2L(1:n_basbas_fns), fn2species(1:n_basbas_fns), &
  &                       max_basbas_L, max_n_basbas_fnLsp, n_basbas)

  call allocate_basbas()

  basbasfn_l = fn2L(1:n_basbas_fns)
  basbasfn_species = fn2species(1:n_basbas_fns)

  call generate_basbasfn(basbasfn_wave, basbasfn_l, basbasfn_species)

  call generate_full_bas(n_basbas_fns, basbasfn_l, basbasfn_species, &
  &                      max_basbas_L, max_n_basbas_fnLsp, n_basbas, &
  &                      basbas_atom, basbas_l, basbas_m, basbas_fn, &
  &                      Lsp2n_basbas_fnLsp, Lsp2basbas_fn, Lsp2basbas_sp, &
  &                      atom2basbas_off, sp2n_basbas_sp)
  max_n_basbas_sp = maxval(sp2n_basbas_sp)

  call prodbas_tasks_distribution()

  deallocate(tripY, prod_wave, n_prods_uptonow)
  deallocate(basbasfn_wave, fn2species, fn2L)

end subroutine shrink_svd_auxil_basis
!******
