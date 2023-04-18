!****s* FHI-aims/get_real_space_density
!  NAME
!    get_real_space_density
!  SYNOPSIS

subroutine get_real_space_density &  
     ( KS_eigenvector, KS_eigenvector_complex, occ_numbers, partition_tab_std,  &
     hartree_partition_tab_std, rho_std, rho_gradient_std,  basis_l_max, &
     delta_rho_KS_std, delta_rho_gradient_std, rho_change,  density_matrix_sparse )

!  PURPOSE
!
!  The subroutine  obtains the new KS density using the density
!  matrix.
!
!  For the benefit of Pulay mixing later on, we obtain the change between
!  the input and output density, delta_rho_KS = rho_KS - rho_in,
!  rather than overwriting rho itself. So the results are in variable delta_rho_KS.
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use basis
  use cartesian_ylm
  use constants
  use analyze_arrays
  use scalapack_wrapper
  use wf_extrapolation
  use KH_core_states
  use plus_u
  use density_matrix_evaluation
  use load_balancing
  use gt
  implicit none

!  ARGUMENTS

  real*8,     dimension(n_basis, n_states, n_spin,n_k_points_task):: KS_eigenvector
  complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task):: KS_eigenvector_complex
  real*8,     dimension(n_states, n_spin, n_k_points)             :: occ_numbers

  real*8, target, dimension(n_full_points)            :: partition_tab_std
  real*8, target, dimension(n_full_points)            :: hartree_partition_tab_std
  real*8, target, dimension(n_spin, n_full_points)    :: rho_std
  real*8, target, dimension(3, n_spin, n_full_points) :: rho_gradient_std
  integer, dimension(n_species)                 :: basis_l_max
  real*8, target, intent(OUT) :: delta_rho_KS_std(n_full_points, n_spin)
  real*8, target, intent(OUT) :: delta_rho_gradient_std(3, n_full_points, n_spin)
  real*8,  dimension(n_spin)                    :: rho_change
  ! when this routine is called, density_matrix_sparse has either the dimension
  ! (n_hamiltonian_matrix_size) or (n_local_matrix_size)
  ! so we declare it here as a 1D assumed size array
  real*8,  dimension(*)                          :: density_matrix_sparse

! INPUTS
! o KS_eigenvector -- eigenvectors if real eigenvectors are in use
! o KS_eigenvector_complex -- eigenvectors is complex eigenvectors are in use
! o occ_numbers -- occupation numbers of states
! o partition_tab -- partition function values
! o hartree_partition_tab -- hartree potentials partition function values
! o rho -- electron density, rho is only input; what we store is the density residual (i.e. the change in the density)
! o rho_gradient -- gradient of the electron density
! o basis_l_max -- maximum basis l value of the basis functions
! o density_matrix_sparse -- this is the works space for density matrix, it is used if packed matrixes are in use.
!
! OUTPUT
!   o delta_rho_KS -- calculated charge density  -  old rho (input)
!   o delta_rho_gradient -- calculated gradient  -  old rho_gradient (input)
! o rho_change -- how much the electron density changed during the update.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE





  real*8, dimension(:),     allocatable:: dist_tab
  real*8, dimension(:,:),   allocatable:: dist_tab_sq
  real*8, dimension(:),     allocatable:: i_r
  real*8, dimension(:,:,:), allocatable:: dir_tab
  real*8, dimension(:,:),   allocatable:: trigonom_tab

  real*8,dimension(:),  allocatable:: radial_wave
  real*8,dimension(:),  allocatable:: radial_wave_deriv
  real*8,dimension(:),  allocatable:: radial_wave_2nd_deriv
  real*8,dimension(:,:),allocatable:: wave
  real*8,dimension(:),  allocatable:: kinetic_wave

  integer,dimension(:),allocatable :: i_basis
  
  integer :: n_compute_c

  real*8 coord_current(3)

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_centers_integrals)
  integer :: spline_array_end(n_centers_integrals)

! VB - renewed index infrastructure starts here

  real*8 :: one_over_dist_tab(n_max_compute_atoms)

      ! indices for basis functions that are nonzero at current point

      integer :: rad_index(n_max_compute_atoms)
      integer :: wave_index(n_max_compute_fns_dens)
      integer :: l_index(n_max_compute_fns_dens)
      integer :: l_count(n_max_compute_fns_dens)
      integer :: fn_atom(n_max_compute_fns_dens)

      ! indices for known zero basis functions at current point
      integer :: n_zero_compute
      integer :: zero_index_point(n_max_compute_dens)

      ! active atoms in current batch
      integer :: n_batch_centers
      integer :: batch_center(n_centers_integrals)

  !     other local variables

  integer :: l_ylm_max
  integer :: n_points

  ! (n_basis, n_states, n_spin) on purpose
  ! beware: inside the subroutines (evaluate_KS_density_v1 for instance)
  ! dimension is (n_compute, max_occ_number)
  ! that's a trick to get a continuous data flow and not
  ! bad programming

  real*8 :: temp_rho(n_max_batch_size,n_spin)
  real*8 :: temp_rho_gradient(3,n_max_batch_size,n_spin)
  real*8,dimension(:,:),allocatable :: temp_rho_small

  integer, dimension(:,:), allocatable :: index_lm
  real*8,  dimension(:,:), allocatable :: ylm_tab

  !     only allocated and referenced for gradient functionals
  real*8, dimension(:,:),  allocatable :: dylm_dtheta_tab
  real*8, dimension(:,:),  allocatable :: scaled_dylm_dphi_tab
  real*8, dimension(:,:,:),allocatable :: gradient_basis_wave

  real*8, dimension(:,:),  allocatable :: density_matrix
  real*8, dimension(:),  allocatable :: density_matrix_con
  real*8, dimension(:,:),  allocatable :: work
  real*8, dimension(:,:),  allocatable :: grad_work
  complex*16, dimension(:,:),allocatable :: work_complex
!  real*8, allocatable :: density_matrix_extrapolated(:,:,:)
!  real*8, allocatable :: density_matrix_sparse_extrapolated(:,:)

  !     counters
  integer :: i_l
  integer :: i_m
  integer :: i_coord
  integer :: i_state
  integer :: i_point

  integer :: i_full_points
  integer :: i_full_points_2
  integer :: i_full_points_3

  integer :: i_spin,  i_spin_2
  integer :: i_index, i_k_point
  integer :: i_basis_index, j_basis_index
  integer :: i_my_batch
  integer :: info
  logical :: is_restarted
  character*120 :: info_str

  ! Load balancing stuff

  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)
  real*8, pointer :: hartree_partition_tab(:)
  real*8, pointer :: rho(:,:)
  real*8, pointer :: rho_gradient(:,:,:)
  real*8, pointer :: delta_rho_KS(:,:)
  real*8, pointer :: delta_rho_gradient(:,:,:)

  integer i, j, i_off, n_bp
  integer, allocatable :: ins_idx(:)

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all


  call localorb_info( &
       "Evaluating new KS density using the density matrix", use_unit,'(2X,A)', OL_norm )

  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)

  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then

    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab

    allocate(hartree_partition_tab(batch_perm(n_bp)%n_full_points))
    call permute_point_array(n_bp,1,hartree_partition_tab_std,hartree_partition_tab)

    allocate(rho(n_spin,batch_perm(n_bp)%n_full_points))
    call permute_point_array(n_bp,n_spin,rho_std,rho)

    nullify(rho_gradient)

    allocate(delta_rho_KS(batch_perm(n_bp)%n_full_points,n_spin))

    nullify(delta_rho_gradient)

    allocate(ins_idx(batch_perm(n_bp)%n_basis_local))

  else

    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std
    hartree_partition_tab => hartree_partition_tab_std
    rho => rho_std
    rho_gradient => rho_gradient_std
    delta_rho_KS => delta_rho_KS_std
    delta_rho_gradient => delta_rho_gradient_std

  endif

  !-----------------------------------------------------------------------------

  ! initialize charge density convergence criterion
  rho_change = 0.d0

  delta_rho_KS = 0.d0

  if(packed_matrix_format == PM_none)then
     if(.not. allocated(density_matrix))then
        allocate(density_matrix(n_centers_basis_T, n_centers_basis_T),stat=info)
        call check_allocation(info, 'density_matrix                ')
     end if
  else
     ! Allocate dummy since otherways compiling with checking doesn't work
     allocate(density_matrix(1,1))
  end if

  ! --- start actual calculation

  time_work = 0
  time_all  = 0

  do i_spin = 1, n_spin

     ! --- density matrix construction/administration

     is_restarted = .false.
     if (use_scalapack .and. restart_read .and. .not.use_local_index) then
        call restart_scalapack_read(density_matrix_sparse, i_spin, is_restarted)
     end if

     if (.not. is_restarted) then
           call evaluate_densmat(KS_eigenvector, KS_eigenvector_complex, occ_numbers,&
           &                     density_matrix, density_matrix_sparse, i_spin, .false.)
     end if

     if (use_scalapack .and. restart_write .and. .not.use_local_index) then
        call restart_scalapack_write(density_matrix_sparse, i_spin)
     end if

     density_matrix (:,:) = green_fn(:,:)

      i_index=0
      density_matrix_sparse (1:n_hamiltonian_matrix_size) = 0.d0
      do i_basis_index = 1, n_basis, 1
        do j_basis_index =1, i_basis_index, 1
          i_index = i_index+1
          density_matrix_sparse (i_basis_index) = density_matrix(i_basis_index, j_basis_index)
        enddo
      enddo

     ! ---------------- end construct density matrix -------------------------

     if(.not. allocated(dist_tab))then
        allocate(dist_tab(n_centers_integrals),stat=info)
        call check_allocation(info, 'dist_tab                      ')
     end if
     if(.not. allocated(dist_tab_sq))  then
        allocate(dist_tab_sq(n_centers_integrals, n_max_batch_size),stat=info)
        call check_allocation(info, 'dist_tab_sq                   ')
     end if
     if(.not. allocated(i_r))then
        allocate(i_r(n_max_compute_atoms),stat=info)
        call check_allocation(info, 'i_r                           ')
     end if
     if(.not. allocated(dir_tab)) then
        allocate(dir_tab(3, n_centers_integrals, n_max_batch_size),stat=info)
        call check_allocation(info, 'dir_tab                       ')
     end if
     if(.not. allocated(trigonom_tab)) then
        allocate(trigonom_tab(4, n_max_compute_atoms),stat=info)
        call check_allocation(info, 'trigonom_tab                  ')
     end if
     if(.not. allocated(density_matrix_con))then
        allocate(density_matrix_con(n_max_compute_dens*n_max_compute_dens),stat=info)
        call check_allocation(info, 'density_matrix_con            ')
     end if
     if(.not. allocated(work))then
        allocate(work(n_max_compute_dens, n_max_batch_size),stat=info)
        call check_allocation(info, 'work                          ')
     end if
     if(.not. allocated(grad_work))then
        allocate(grad_work(n_max_compute_dens, n_max_batch_size),stat=info)
        call check_allocation(info, 'grad_work                     ')
     end if
     if(.not. allocated(i_basis))then
        allocate(i_basis(n_centers_basis_T),stat=info)
        call check_allocation(info, 'i_basis                       ')
     end if
     if(.not.allocated(radial_wave))then
        allocate(radial_wave(n_max_compute_fns_dens),stat=info)
        call check_allocation(info, 'radial_wave                   ')
     end if
     if(.not. allocated(radial_wave_deriv))then
        allocate(radial_wave_deriv(n_max_compute_fns_dens),stat=info)
        call check_allocation(info, 'radial_wave_deriv             ')
     end if
     if(.not. allocated(radial_wave_2nd_deriv))then
        allocate(radial_wave_2nd_deriv(n_max_compute_fns_dens),stat=info)
        call check_allocation(info, 'radial_wave_2nd_deriv         ')
     end if
     if(.not. allocated(wave))then
        allocate(wave(n_max_compute_dens, n_max_batch_size),stat=info)
        call check_allocation(info, 'wave                          ')
     end if
     if(.not. allocated(kinetic_wave))then
        allocate(kinetic_wave(n_max_compute_fns_dens),stat=info)
        call check_allocation(info, 'kinetic_wave                  ')
     end if

     l_ylm_max = l_wave_max

     if(.not. allocated( ylm_tab))then
        allocate( ylm_tab( (l_ylm_max+1)**2,n_max_compute_atoms),stat=info )
        call check_allocation(info, 'ylm_tab                       ')
     end if
     if(.not. allocated( index_lm))then
        allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),stat=info)
        call check_allocation(info, 'index_lm                      ')
      end if

     ! initialize index_lm
     i_index = 0
     do i_l = 0, l_ylm_max, 1
        do i_m = -i_l, i_l
           i_index = i_index + 1
           index_lm(i_m, i_l) = i_index
        enddo
     enddo

     call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
     time0 = mpi_wtime()

     ! -------- go over the grids and construc the density using density matrix -------------

     i_basis_fns_inv = 0

     i_full_points = 0
     i_full_points_2 = 0
     i_full_points_3 = 0

     do i_my_batch = 1, n_my_batches_work, 1

           n_compute_c = 0
           i_basis = 0

           i_point = 0

           ! loop over one batch
           do i_index = 1, batches_work(i_my_batch)%size, 1

              i_full_points = i_full_points + 1

              if (max(partition_tab(i_full_points),&
                   hartree_partition_tab(i_full_points)).gt.0.d0) then

                 i_point = i_point+1

                 !     get current integration point coordinate
                 coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)

                 if(n_periodic>0)then
                    call map_to_center_cell( coord_current)
                 end if

                 ! compute atom-centered coordinates of current integration point,
                 ! as viewed from all atoms
                 call tab_atom_centered_coords_p0 &
                      ( coord_current,  &
                      dist_tab_sq(1,i_point),  &
                      dir_tab(1,1,i_point), &
                      n_centers_integrals, centers_basis_integrals )



                 !     determine which basis functions are relevant at current integration point,
                 !     and tabulate their indices

                 ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
                 if (.not.prune_basis_once) then
                   call prune_basis_p2 &
                   ( dist_tab_sq(1,i_point), &
                     n_compute_c, i_basis,  &
                     n_centers_basis_T, n_centers_basis_integrals, inv_centers_basis_integrals  )
                 end if
              end if
           enddo ! end loop over the angular shell

           if (prune_basis_once) then
              n_compute_c = batches_work(i_my_batch)%batch_n_compute
              i_basis(1:n_compute_c) = batches_work(i_my_batch)%batch_i_basis
           end if

           if(packed_matrix_format /= PM_none )then
              if(use_batch_permutation > 0) then
                do i=1,n_compute_c
                  ins_idx(i) = batch_perm(n_bp)%i_basis_glb_to_loc(i_basis(i))
                enddo
                density_matrix_con(1:n_compute_c*n_compute_c) = 0
                do i=1,n_compute_c
                  i_off = (ins_idx(i)*(ins_idx(i)-1))/2
                  do j=1,i
                    density_matrix_con(j+(i-1)*n_compute_c) = density_matrix_sparse(ins_idx(j)+i_off)
                  enddo
                enddo
              else
                call  prune_density_matrix_sparse(density_matrix_sparse, density_matrix_con, &
                     n_compute_c, i_basis)
              endif
           else
             call  prune_density_matrix(density_matrix, density_matrix_con, &
                  n_compute_c, i_basis)
           end if


           ! from list of n_compute active basis functions in batch, collect all atoms that are ever needed in batch.
           call collect_batch_centers_p2 &
           ( n_compute_c, i_basis, n_centers_basis_T, n_centers_basis_integrals, inv_centers_basis_integrals, &
             n_batch_centers, batch_center &
           )

           n_points = i_point

           if (n_compute_c.gt.0) then

              ! Determine all radial functions, ylm functions and their derivatives that
              ! are best evaluated strictly locally at each individual grid point.
              i_point = 0
              do i_index = 1, batches_work(i_my_batch)%size, 1

                 i_full_points_2 = i_full_points_2 + 1

                 if (max(partition_tab(i_full_points_2),&
                      hartree_partition_tab(i_full_points_2)).gt.0.d0) then

                    i_point = i_point+1
                    n_compute_atoms = 0
                    n_compute_fns = 0


                    ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                    ! Are stored in a compact spline array that can be accessed by spline_vector_waves,
                    ! without any copying and without doing any unnecessary operations.
                    ! The price is that the interface is no longer explicit in terms of physical
                    ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.


                    call prune_radial_basis_p2 &
                    ( n_max_compute_atoms, n_max_compute_fns_dens, &
                       dist_tab_sq(1,i_point), dist_tab, dir_tab(1,1,i_point), &
                       n_compute_atoms, atom_index, atom_index_inv, &
                       n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                       i_atom_fns, spline_array_start, spline_array_end, &
                       n_centers_basis_integrals, centers_basis_integrals, n_compute_c, i_basis, &
                       n_batch_centers, batch_center, &
                       one_over_dist_tab, rad_index, wave_index, l_index, l_count, &
                       fn_atom, n_zero_compute, zero_index_point &
                    )

                    ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                    ! for all atoms which are actually relevant
                    call tab_local_geometry_p2 &
                      ( n_compute_atoms, atom_index, &
                        dist_tab,  &
                        i_r &
                      )

                    ! Determine all needed radial functions from efficient splines

                    ! Now evaluate radial functions u(r) from the previously stored compressed
                    ! spline arrays
                    call evaluate_radial_functions_p0 &
                         (   spline_array_start, spline_array_end, &
                         n_compute_atoms, n_compute_fns,  &
                         dist_tab, i_r(1), &
                         atom_index, i_basis_fns_inv, &
                         basis_wave_ordered, radial_wave(1), &
                         .false., n_compute_c, n_max_compute_fns_dens   &
                         )


                    ! We compute the needed
                    ! ylm pieces separately from the cartesian evaluation that is done
                    ! for the hessians ...

                    ! compute trigonometric functions of spherical coordinate angles
                    ! of current integration point, viewed from all atoms

                       call tab_trigonom_p0 &
                            ( n_compute_atoms, dir_tab(1,1,i_point),  &
                            trigonom_tab(1,1) &
                            )

                         ! tabulate distance and Ylm's w.r.t. other atoms
                         call tab_wave_ylm_p0 &
                            ( n_compute_atoms, atom_index,  &
                            trigonom_tab(1,1), basis_l_max,  &
                            l_ylm_max, &
                            ylm_tab(1,1) )

                    ! tabulate total wave function value for each basis function in all cases -
                    ! but only now are we sure that we have ylm_tab ...

                    ! tabulate total wave function value for each basis function
                    call evaluate_waves_p2  &
                    ( n_compute_c, n_compute_atoms, n_compute_fns, &
                      l_ylm_max, ylm_tab, one_over_dist_tab,   &
                      radial_wave, wave(1,i_point), &
                      rad_index, wave_index, l_index, l_count, fn_atom, &
                      n_zero_compute, zero_index_point &
                    )

                    ! reset i_basis_fns_inv
                    i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0

                 end if ! end if (partition_tab.gt.0)
              enddo ! end loop over a batch
              ! - all quantities that are evaluated pointwise are now known ...


              if (n_points.gt.0) then
                 ! Now perform all operations which are done across the entire integration
                 ! shell at once.

                 ! New density is always evaluated

                 call evaluate_KS_density_densmat &
                      (  n_points, wave(1,1), n_compute_c,   &
                      temp_rho(1,i_spin), n_max_compute_dens, &
                      n_centers_basis_T, density_matrix_con, work(1,1) )


              end if ! (n_points>0)


              ! calculate change in electron density
              i_point = 0
              do i_index = 1, batches_work(i_my_batch)%size, 1

                 i_full_points_3 = i_full_points_3 + 1

                 if (max(partition_tab(i_full_points_3),&
                      hartree_partition_tab(i_full_points_3)).gt.0) then

                    i_point = i_point + 1

                    ! and local change in charge density and gradient for Pulay mixing

                    if (spin_treatment.eq.0) then

                       delta_rho_KS(i_full_points_3, 1) =   &
                            delta_rho_KS(i_full_points_3, 1) &
                            + temp_rho(i_point,1) !- rho(1, i_full_points_3

                       rho(1, i_full_points_3) = delta_rho_KS(i_full_points_3, 1) 

                    !NO SPIN FOR NOW
!                    elseif (spin_treatment.eq.1) then
!
!
!                       ! Note that delta_rho_KS mixes the spin components of the density - but, in the dm
!                       ! formalism, we do first an _outermost_ loop over spin up, and only then over spin down.
!                       if(i_spin == 1)then
!
!
!                          delta_rho_KS(i_full_points_3, 1) =  delta_rho_KS(i_full_points_3, 1)   &
!                               +( temp_rho(i_point, 1) )  !&
!                          !     - ( rho(1,i_full_points_3)  )
!
!
!                          delta_rho_KS(i_full_points_3, 2) =   delta_rho_KS(i_full_points_3, 2)   &
!                               +( temp_rho(i_point, 1) )  !&
!!                               - ( rho(1,i_full_points_3))
!
!                       else
!
!                          delta_rho_KS(i_full_points_3, 1) =  delta_rho_KS(i_full_points_3, 1)   &
!                               +( temp_rho(i_point, 2) )  !&
!!                               - ( rho(2,i_full_points_3) )
!
!                          delta_rho_KS(i_full_points_3, 2) =   delta_rho_KS(i_full_points_3, 2)   &
!                               -( temp_rho(i_point, 2) )  !&
!!                               + ( rho(2, i_full_points_3) )
!                       end if
                    end if

                    ! prepare charge density root-mean-square distance

                    if(i_spin == n_spin)then

                       do i_spin_2 = 1,n_spin

                          rho_change(i_spin_2) = rho_change(i_spin_2) +  &
                               partition_tab(i_full_points_3) *  &
                               delta_rho_KS(i_full_points_3, i_spin_2)**2
                       end do
                    end if
                 endif
              end do

           else
               ! Even if n_compute .eq. 0 for the entire current batch of grid points, we still need to
               ! make sure that the density _change_ at this point is ( zero minus previous density ).
               ! This ensures that. even for a zero KS density, a potentially non-zero initialization density
               ! is subtracted to properly account for the density change ....

               do i_index = 1, batches_work(i_my_batch)%size, 1

                  i_full_points_2 = i_full_points_2 + 1
                  i_full_points_3 = i_full_points_3 + 1

                  if (max(partition_tab(i_full_points_3),&
                        hartree_partition_tab(i_full_points_3)).gt.0.d0) then

                     ! only evaluate once both spin components are known
                     if (i_spin.eq.n_spin) then

                        ! local change in charge density and gradient for Pulay mixing

                        if (spin_treatment.eq.0) then

                           delta_rho_KS(i_full_points_3, 1) =   &
                                 delta_rho_KS(i_full_points_3, 1) !&
!                                 - rho(1, i_full_points_3)
                       rho(1, i_full_points_3) = delta_rho_KS(i_full_points_3, 1) 

!NO SPIN FOR NOW
!                        elseif (spin_treatment.eq.1) then
!
!                           delta_rho_KS(i_full_points_3, 1) =   &
!                                 delta_rho_KS(i_full_points_3, 1) !-  &
!!                                 ( rho(1,i_full_points_3) + rho(2,i_full_points_3) )
!
!                           delta_rho_KS(i_full_points_3, 2) =   &
!                                 delta_rho_KS(i_full_points_3, 2) !-  &
!!                                 ( rho(1,i_full_points_3) - rho(2, i_full_points_3) )
!
                        end if

                        ! prepare charge density root-mean-square distance
                        do i_spin_2 = 1, n_spin, 1

                           rho_change(i_spin_2) = rho_change(i_spin_2) +  &
                                 partition_tab(i_full_points_3) *  &
                                 delta_rho_KS(i_full_points_3, i_spin_2)**2

                        enddo
                     end if ! (i_spin.eq.n_spin)

                  endif !(partition_tab)

               enddo ! (loop over all grid points, for empty batches)

           end if  ! end if (n_compute.gt.0)
        ! end if ! end distribution over threads
     end do ! end loop over grid batches

     ! Get work time and total time after barrier
     time_work = time_work + mpi_wtime()-time0
     call mpi_barrier(mpi_comm_world,info)
     time_all  = time_all  +mpi_wtime()-time0

     !--------------------------------- end go over the grids ------------------------------
  end do ! end i_spin
  rho (:,:) = rho (:,:) *2.d0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for getting density from density matrix: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)

  ! broadcast the result to all threads
  call sync_density(rho_change)

  do i_spin = 1, n_spin, 1
     rho_change(i_spin) = sqrt(rho_change(i_spin))
  enddo

  if(use_batch_permutation > 0) then

    do i_spin = 1, n_spin
      call permute_point_array_back(n_bp,1,delta_rho_KS(:,i_spin),delta_rho_KS_std(:,i_spin))
    enddo

    deallocate(hartree_partition_tab)
    deallocate(rho)

    deallocate(delta_rho_KS)

    deallocate(ins_idx)

  endif

  !---------- finally, deallocate stuff -----------------

  if (allocated( temp_rho_small       )) deallocate( temp_rho_small       )
  if (allocated( scaled_dylm_dphi_tab )) deallocate( scaled_dylm_dphi_tab )
  if (allocated( dylm_dtheta_tab      )) deallocate( dylm_dtheta_tab      )
  if (allocated( gradient_basis_wave  )) deallocate( gradient_basis_wave  )
  if (allocated( index_lm             )) deallocate( index_lm             )
  if (allocated( ylm_tab              )) deallocate( ylm_tab              )
  if (allocated( kinetic_wave         )) deallocate( kinetic_wave         )
  if (allocated( wave                 )) deallocate( wave                 )
  if (allocated( radial_wave_2nd_deriv)) deallocate( radial_wave_2nd_deriv)
  if (allocated( radial_wave_deriv    )) deallocate( radial_wave_deriv    )
  if (allocated( radial_wave          )) deallocate( radial_wave          )
  if (allocated( i_basis              )) deallocate( i_basis              )
  if (allocated( work                 )) deallocate( work                 )
  if (allocated( grad_work            )) deallocate( grad_work            )
  if (allocated( density_matrix_con   )) deallocate( density_matrix_con   )
  if (allocated( trigonom_tab         )) deallocate( trigonom_tab         )
  if (allocated( dir_tab              )) deallocate( dir_tab              )
  if (allocated( i_r                  )) deallocate( i_r                  )
  if (allocated( dist_tab_sq          )) deallocate( dist_tab_sq          )
  if (allocated( dist_tab             )) deallocate( dist_tab             )
  if (allocated( work_complex         )) deallocate( work_complex         )
  if (allocated( density_matrix       )) deallocate( density_matrix       )

end subroutine get_real_space_density 
!******
