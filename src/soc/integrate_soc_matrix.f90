!****f* FHI-aims/integrate_soc_matrix
!*  NAME
!*    integrate_soc_matrix
!*  SYNOPSIS
subroutine integrate_soc_matrix(rho_std, hartree_potential_std, &
                                partition_tab_std, soc_matrix)
!*  PURPOSE
!*    Integrates the real-space SOC matrix elements.  The SOC operator has the
!*    form
!*        V_SOC \prop sigma . p x Vp
!*    and we calculate the three components of the vector operator (p x Vp)
!*    here.  The contribution due to the Pauli matrix will come later when
!*    constructing the SOC Hamiltonian.
!*
!*    Because (p x Vp) is a cross product, there are six terms that we need to
!*    calculate.  We index them in terms of the Cartesian coordinate (x,y,z)
!*    which they contribute to, and the sign (+1, -1) in front of them.
!*  USES
  use dimensions, only: n_spin, n_full_points, n_hamiltonian_matrix_size, &
      n_centers_integrals, n_centers, n_max_compute_fns_ham, &
      n_max_compute_atoms, include_pw_lda_in_v_soc, l_wave_max,  n_basis_fns, &
      n_centers_basis_I, n_my_batches, n_max_batch_size, n_max_compute_ham, &
      n_periodic
  use constants, only: light_speed_sq
  use basis, only: basis_wave_ordered, basis_deriv_ordered
  use runtime_choices, only: prune_basis_once, use_local_index
  use mpi_tasks, only: aims_stop
  use grids, only: batch_of_points, batches
  use mpi_utilities, only : batch_task_list
  use pbc_lists, only: centers_basis_integrals, inv_centers_basis_integrals
  use synchronize_mpi_basic, only: sync_vector
  use localorb_io, only: localorb_info
  use cartesian_ylm, only: n_max_cartesian, evaluate_cartesians, &
      initialize_cartesian_ylm, evaluate_cartesian_hessian_and_gradient_terms_p2
  use load_balancing, only: use_batch_permutation, batch_perm, &
      permute_point_array
  use species_data, only: species_name,l_shell_max
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use xc, only: xcpot_pw91_lda
  implicit none
!*  ARGUMENTS
  real*8, target, intent(in)  :: rho_std(n_spin, n_full_points)
  real*8, target, intent(in)  :: hartree_potential_std(n_full_points)
  real*8, target, intent(in)  :: partition_tab_std(n_full_points)
  ! when this routine is called, soc_matrix has either the dimension
  ! (n_hamiltonian_matrix_size, n_spin) or (n_local_matrix_size, n_spin)
  ! so we declare it here as a 1D assumed size array
  real*8, intent(out) :: soc_matrix(*)
!*  INPUTS
!*    o rho_std - electron density for all points assigned to current MPI task
!*    o hartree_potential_std - Hartree potential for all points assigned to current MPI task
!*    o partition_tab_std - Integration weights for all points assigned to the current MPI tas
!*  OUTPUT
!*    o soc_matrix            - the real-space SOC matrix elements for interaction strength between
!*                              the basis functions
!*  AUTHORS
!*    William Huhn and Matthias Gramzow
!*  NOTES
!*    The bulk of this code was lifted from integrate_hamiltonian_matrix_p2.
!*
!*    At present, only the Hartree term enters into the potential.  This is done
!*    to considerably simplify the calculation, as semi-local and especially hybrid
!*    functionals would complicate this subroutine (and in the case of hybrids,
!*    significantly increase the runtime) without notably changing the results, as
!*    the SOC interaction strength dominates in the nuclear region where the relative
!*    contribution from the XC function is minimal.  There is a flag
!*    include_pw_lda_in_v_soc which will calculate the PW91 XC contribution and add it
!*    to the potential, yielding the true effective potential for an PW91 calculation.
!*    As is documented in the paper (see below), this changes values by 1%.
!*
!*    The implementation of second-variational SOC in FHI-aims is published in
!*      Huhn and Blum, Phys. Rev. Materials 1, 033803 (2017)
!*      https://dx.doi.org/10.1103/PhysRevMaterials.1.033803
!*
!*    This subroutine implements step 3 in Section III.3 from said paper.
!*  SEE ALSO
!*    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!*    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!*    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!*    Computer Physics Communications (2008), submitted.
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!*    the terms and conditions of the respective license agreement."
!*  SOURCE

! local_variables

  real*8, dimension(:), allocatable ::  potential_shell!, xc_shell
  real*8, dimension(:,:,:), allocatable :: cartesians
  real*8, dimension(:,:,:), allocatable :: sum_gradient
  real*8, dimension(:,:,:), allocatable :: sum_hessian

  real*8,dimension(4) :: temp

  integer :: l_ylm_max
  real*8, dimension(:,:), allocatable :: ylm_tab

  real*8 coord_current(3)
  real*8, dimension(:,:), allocatable :: dist_tab
  real*8, dimension(:,:), allocatable :: dist_tab_sq
  real*8 i_r(n_max_compute_atoms)
  real*8, dimension(:,:,:), allocatable :: dir_tab

  real*8,dimension(:)  ,allocatable:: radial_wave
  real*8,dimension(:)  ,allocatable:: radial_wave_deriv
  real*8,dimension(:), allocatable  ::  soc_shell

  real*8 :: partition(n_max_batch_size)

  real*8, dimension(:,:,:), allocatable :: gradient_basis_wave

  integer :: n_compute_c
  integer,dimension(:),allocatable :: i_basis

  integer :: n_points

  integer :: n_compute_fns
  integer, dimension(:), allocatable :: i_basis_fns
  integer, dimension(:,:), allocatable :: i_basis_fns_inv
  integer, dimension(:), allocatable :: i_atom_fns

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_centers_integrals)
  integer :: spline_array_end(n_centers_integrals)

! VB - renewed index infrastructure starts here

  real*8 one_over_dist_tab(n_max_compute_atoms)

  ! indices for basis functions that are nonzero at current point

  integer :: rad_index(n_max_compute_atoms)
  integer :: wave_index(n_max_compute_fns_ham)
  integer :: l_index(n_max_compute_fns_ham)
  integer :: l_count(n_max_compute_fns_ham)
  integer :: fn_atom(n_max_compute_fns_ham)

  ! indices for known zero basis functions at current point
  integer :: n_zero_compute
  integer :: zero_index_point(n_max_compute_ham)

  ! active atoms in current batch
  integer :: n_batch_centers
  integer :: batch_center(n_centers_integrals)

  !  counters
  integer i_basis_1, i_basis_2, i_spin
  integer i_index
  integer i_coord
  integer :: i, j, i_off, i_pauli, i_sign

  integer i_point
  integer :: i_full_points
  integer :: i_full_points_2

  integer :: i_fn_1, i_fn_2, i_cell, i_size, i_index_real
  character l_to_str ! Function found in convert_l_str.f90
  character l_char_1, l_char_2

  character*100 :: info_str

  integer :: i_my_batch

  integer :: info

  !  Quantities for calculation of the XC potential
  !  Currently NOT using, but leaving for future
  !  improvements
  !  WELCOME TO THE FUTURE
  !  It's a pretty terrible place
  !  Trump is president, and people eat kale
  real*8 :: pot_xc_pw_lda(n_spin)
  real*8 :: en_density_xc_pw_lda
  real*8 :: en_density_x_pw_lda(n_spin)
  real*8 :: en_density_c_pw_lda
  real*8 :: temp_rho(n_spin)

  ! Load balancing stuff
  integer :: ld_soc_matrix    ! Leading dimension of soc_matrix
  integer :: n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches used
  integer :: n_bp

  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)
  real*8, pointer :: rho(:,:)
  real*8, pointer :: hartree_potential(:)

  character(*), parameter :: func = 'integrate_soc_matrix'

  !  begin work

  write(info_str,'(2X,A)') "Integrating SOC matrix: batch-based integration p2 version."
  call localorb_info( info_str )

  l_ylm_max = l_wave_max

  call aims_allocate(dist_tab, n_centers_integrals, n_max_batch_size, "+dist_tab")
  call aims_allocate(dist_tab_sq, n_centers_integrals, n_max_batch_size, "+dist_tab_sq")
  call aims_allocate(dir_tab, 3, n_centers_integrals, n_max_batch_size, "+dir_tab")
  call aims_allocate(i_basis_fns, n_basis_fns*n_centers_integrals, "+i_basis_fns")
  call aims_allocate(i_basis_fns_inv, n_basis_fns, n_centers, "+i_basis_fns_inv")
  call aims_allocate(i_atom_fns, n_basis_fns*n_centers_integrals, "+i_atom_fns")
  call aims_allocate(gradient_basis_wave, n_max_compute_ham, 3, n_max_batch_size, "+gradient_basis_wave")
  gradient_basis_wave = 0

  call aims_allocate(ylm_tab, (l_ylm_max+1)**2, n_max_compute_atoms, "+ylm_tab" )
  call aims_allocate(radial_wave, n_max_compute_fns_ham, "+radial_wave")
  call aims_allocate(radial_wave_deriv, n_max_compute_fns_ham, "+radial_wave_deriv")
  call aims_allocate(i_basis, n_centers_basis_I, "+i_basis")
  call aims_allocate(soc_shell, n_max_compute_ham*n_max_compute_ham, "+soc_shell")

  ! Because cartesian uses non-standard indexing in the second argument,
  ! requires all indices
  call aims_allocate(cartesians, 1, n_max_cartesian, 0, l_wave_max, 1, &
       n_max_compute_atoms, "+cartesians")
  call aims_allocate(sum_gradient, 3, (l_wave_max+1) ** 2, n_max_compute_atoms,&
       "+sum_gradient")
  call aims_allocate(sum_hessian, 6, (l_wave_max+1) ** 2, n_max_compute_atoms, &
       "+sum_hessian")
  call aims_allocate(potential_shell, n_max_batch_size, "+potential_shell")

  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then
    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab

    allocate(rho(n_spin,batch_perm(n_bp)%n_full_points))
    call permute_point_array(n_bp,n_spin,rho_std,rho)

    allocate(hartree_potential(batch_perm(n_bp)%n_full_points))
    call permute_point_array(n_bp,1,hartree_potential_std,hartree_potential)

    ld_soc_matrix    =  batch_perm(n_bp)%n_local_matrix_size
  else
    n_my_batches_work =  n_my_batches
    batches_work      => batches
    partition_tab     => partition_tab_std
    rho               => rho_std
    hartree_potential => hartree_potential_std

    ld_soc_matrix    =  n_hamiltonian_matrix_size
  end if

  call initialize_cartesian_ylm(l_wave_max)

  i_basis_fns_inv                = 0
  i_full_points                  = 0
  i_full_points_2                = 0

  soc_matrix(1:3*ld_soc_matrix) = 0.0

  do i_my_batch = 1, n_my_batches_work, 1
    n_compute_c         = 0
    i_basis             = 0
    i_point             = 0

    gradient_basis_wave = 0.0d0
    partition           = 0.0d0
    potential_shell     = 0.0d0
    soc_shell           = 0.0d0

    do i_index = 1, batches_work(i_my_batch)%size, 1

      i_full_points_2 = i_full_points_2 + 1
      if (partition_tab(i_full_points_2).gt.0.d0) then
        i_point = i_point+1
        ! get current integration point coordinate
        coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)

        if(n_periodic > 0)then
          call map_to_center_cell(coord_current(1:3) )
        end if

        ! compute atom-centered coordinates of current integration point,
        ! as viewed from all atoms
        call tab_atom_centered_coords_p0 &
              ( coord_current,  &
              dist_tab_sq(1,i_point),  &
              dir_tab(1,1,i_point), &
              n_centers_integrals, centers_basis_integrals )

        ! determine which basis functions are relevant at current integration
        ! point, and tabulate their indices

        ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are
        ! actually needed
        if (.not.prune_basis_once) then
          call prune_basis_p2( &
                 dist_tab_sq(1,i_point), &
                 n_compute_c, i_basis,  &
                 n_centers_basis_I, n_centers_integrals, &
                 inv_centers_basis_integrals)
        end if
      end if !partition_tab .gt. 0
    end do ! i_index loop over batch

    if (prune_basis_once) then
      n_compute_c = batches_work(i_my_batch)%batch_n_compute
      i_basis(1:n_compute_c) = batches_work(i_my_batch)%batch_i_basis
    end if

    ! from list of n_compute active basis functions in batch, collect all atoms
    ! that are ever needed in batch.
    call collect_batch_centers_p2( n_compute_c, i_basis, n_centers_basis_I, &
           n_centers_integrals, inv_centers_basis_integrals, n_batch_centers, &
           batch_center)

    n_points = i_point
    if (n_compute_c.gt.0) then
      i_point = 0
      ! As the SOC operator has two gradients, we need to calculate the
      ! gradient of the basis element gradients
      ! This is mostly taken from integrate_hamiltonian, specifically the
      ! portion where the basis element gradients are calculated as input
      ! into the GGA functional
      do i_index = 1, batches_work(i_my_batch)%size, 1
        i_full_points = i_full_points + 1
        if (partition_tab(i_full_points).gt.0.d0) then
          i_point = i_point+1

          partition(i_point) = partition_tab(i_full_points)
          potential_shell(i_point) = hartree_potential(i_full_points)

          ! The following includes the PW-LDA XC potential evaluated at the
          ! current electron density in the effective potential entering into
          ! the SOC operator.
          ! This is intended for the high density near the nucleus
          ! At present, we do not make this a default because we are unclear
          ! how to properly include spin-polarized potentials, because there
          ! is no notion of "the up potential" or "the down potential" in the
          ! SOC operator; it is non-collinear and couples spin channels
          if (include_pw_lda_in_v_soc) then
            do i_spin = 1, n_spin
              ! Something silly that needs to be done due to usage of pointers
              temp_rho(i_spin) = rho(i_spin,i_full_points)
            end do
            call xcpot_pw91_lda( temp_rho, pot_xc_pw_lda, en_density_xc_pw_lda,&
                 en_density_x_pw_lda, en_density_c_pw_lda)

            if (n_spin.ne.1) then
              call aims_stop('include_pw_lda_in_v_soc only supports spin-non-polarized calculations.',func)
            end if
            potential_shell(i_point) = potential_shell(i_point) + &
                 pot_xc_pw_lda(1)
          end if

          n_compute_atoms = 0
          n_compute_fns = 0

          ! All these function calls generate inputs into basis function
          ! gradients
          call prune_radial_basis_p2 &
                 ( n_max_compute_atoms, n_max_compute_fns_ham, &
                   dist_tab_sq(1,i_point), dist_tab(1,i_point), dir_tab(1,1,i_point), &
                   n_compute_atoms, atom_index, atom_index_inv, &
                   n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                   i_atom_fns, spline_array_start, spline_array_end, &
                   n_centers_integrals, centers_basis_integrals, n_compute_c, i_basis, &
                   n_batch_centers, batch_center, &
                   one_over_dist_tab, rad_index, wave_index, l_index, l_count, &
                   fn_atom, n_zero_compute, zero_index_point &
                 )
          call tab_local_geometry_p2 &
                 ( n_compute_atoms, atom_index, &
                   dist_tab(1,i_point), i_r )

          ! I think this may not be needed, since this is intended for
          ! calculation of the wave functions (which we don't use here, only the
          ! gradients)
          call evaluate_radial_functions_p0  &
                 ( spline_array_start, spline_array_end,  &
                   n_compute_atoms, n_compute_fns,   &
                   dist_tab(1,i_point), i_r,  &
                   atom_index, i_basis_fns_inv,  &
                   basis_wave_ordered, radial_wave,  &
                   .false. , n_compute_c, n_max_compute_fns_ham )

          ! Here is where it starts to deviate from integrate_hamiltonian.
          ! Best I can tell, this is because we're using the subroutine
          ! evaluate_wave_gradient_cartesian_p2, whereas integrate_hamiltonian
          ! uses evaluate_wave_gradient_p2.
          ! This uses basis_deriv_ordered, and uses
          ! radial_wave_derive(1) instead of kinectic_wave(1), and has .true.
          ! instead of .false.
          call evaluate_radial_functions_p0  &
               ( spline_array_start, spline_array_end,  &
                 n_compute_atoms, n_compute_fns,   &
                 dist_tab(1,i_point), i_r,  &
                 atom_index, i_basis_fns_inv,  &
                 basis_deriv_ordered,   &
                 radial_wave_deriv(1), .true.,  &
                 n_compute_c, n_max_compute_fns_ham )

          call evaluate_cartesians &
               ( dir_tab(1,1,i_point), l_shell_max, l_wave_max,  &
                 atom_index, n_compute_atoms, cartesians(1,0,1) )

          call evaluate_cartesian_hessian_and_gradient_terms_p2 &
               ( l_shell_max, l_wave_max, cartesians(1,0,1), &
                 atom_index, n_compute_atoms,  &
                 ylm_tab,  &
                 sum_gradient,  &
                 sum_hessian(1,1,1))

          ! Here we calculate the basis element gradients
          call evaluate_wave_gradient_cartesian_p2 &
               ( n_compute_c, n_compute_atoms, n_compute_fns, &
                 one_over_dist_tab, dir_tab(1,1,i_point), l_wave_max,  &
                 ylm_tab, radial_wave, radial_wave_deriv, &
                 sum_gradient, gradient_basis_wave(1,1,i_point), &
                 rad_index, wave_index, l_index, l_count, fn_atom, &
                 n_zero_compute, zero_index_point,n_max_compute_ham &
               )

          i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0

        end if ! partition_tab .gt. 0
      end do ! i_index loop over batch

      ! Calculate matrix elements <dbra/dx_{i}|1/r|dket/dx_{j}>
      do i_pauli = 1, 3
        do i_sign = 1, -1, -2
          soc_shell = 0.0d0
          call evaluate_soc_shell(n_points, partition, n_compute_c, &
               potential_shell, gradient_basis_wave, soc_shell, i_pauli, i_sign)

          ! Calculates the anti-symmetric combinations of
          !       <dbra/dx_{i}|1/r|dket/dx_{j}>
          ! (which is what enters into SOC_Hamiltonian due to cross product in
          ! SOC term), see header for more information
          !
          ! I've chosen to move the update for local indexing + load balancing
          ! inside this subroutine, since it does belong in there, unlike other
          ! parts of the code where it is done outside (I suspect for legacy
          ! reasons)
          call update_soc_matrix(n_compute_c, i_basis, soc_shell, &
               i_pauli, i_sign, ld_soc_matrix, soc_matrix)
        end do
      end do
    else
      i_full_points = i_full_points + batches_work(i_my_batch)%size
    end if ! n_compute_c .gt. 0
  end do !loop over batches_work, i_my_batch

  if(.not. use_local_index) call sync_vector(soc_matrix,ld_soc_matrix*3)

  !  Apply the prefactor in Hartree atomic units
  soc_matrix(1:3*ld_soc_matrix) = &
       soc_matrix(1:3*ld_soc_matrix)/(4.0*light_speed_sq)

  ! Deallocate tracked arrays
  call aims_deallocate(dist_tab, "dist_tab")
  call aims_deallocate(dist_tab_sq, "dist_tab_sq")
  call aims_deallocate(dir_tab, "dir_tab")
  call aims_deallocate(i_basis_fns, "i_basis_fns")
  call aims_deallocate(i_basis_fns_inv, "i_basis_fns_inv")
  call aims_deallocate(i_atom_fns, "i_atom_fns")
  call aims_deallocate(gradient_basis_wave, "gradient_basis_wave")
  call aims_deallocate(ylm_tab, "ylm_tab")
  call aims_deallocate(radial_wave, "radial_wave")
  call aims_deallocate(radial_wave_deriv, "radial_wave_deriv")
  call aims_deallocate(i_basis, "i_basis")
  call aims_deallocate(soc_shell, "soc_shell")
  call aims_deallocate(cartesians, "cartesians")
  call aims_deallocate(sum_gradient, "sum_gradient")
  call aims_deallocate(sum_hessian, "sum_hessian")
  call aims_deallocate(potential_shell, "potential_shell")

  if(use_batch_permutation > 0) then
    ! In the load balancing path, these are dynamically allocated pointers, so
    ! we must deallocate them to avoid memory leaks
    deallocate(rho)
    deallocate(hartree_potential)
  end if
end subroutine integrate_soc_matrix
!******
