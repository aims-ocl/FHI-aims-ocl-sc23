subroutine integrate_fxc_matrix( partition_tab, rho, rho_gradient, basis_l_max, fxc_matr )

!  PURPOSE
!
!  This subroutine is to calculate the matrix <u|f_xc|v> in the auxiliary basis.
!  It is the link to libxc and does handle the relevant calls for specified
!  f_xc kernels from the overall possible functionals within libxc.
!  
!  USES
      use dimensions
      use runtime_choices
      use grids
      use geometry
      use basis
      use prodbas
      use constants
      use mpi_tasks
      use mpi_utilities
      use synchronize_mpi
      use pseudodata
      use xc_f03_lib_m

      implicit none

!  ARGUMENTS
      integer, dimension(n_species), intent(in)                  :: basis_l_max(n_species)
      real*8, dimension(n_full_points), intent(in)               :: partition_tab
      real*8, dimension(n_spin,n_full_points), intent(in)        :: rho
      real*8, dimension(3,n_spin, n_full_points), intent(in)     :: rho_gradient
      real*8, dimension(n_basbas, n_basbas, n_spin), intent(out) :: fxc_matr

! INPUTS
!  o partition_tab -- real array, the values of the partition function (for 
!           integration) at each spatial grid point    
!  o rho   -- real array, the electron density at each grid point and for each spin 
!           channel 
!  o rho_gradient -- real array, the electron density gradient
!  o basis_l_max -- the maximal angualr momentum number of the basis function for
!           each species  
! OUTPUTS
!  o fxc_matr -- the matrix elements <u|f_xc|v> within the auxiliary basis
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
!  local variables

      integer :: l_ylm_max
      integer, dimension(:,:), allocatable :: index_lm
      real*8, dimension(:,:,:), allocatable :: ylm_tab

!     flag_xc gets temporarily assigned the value of flag_fxc for TDDFT calculations
!     flag_xc_temp is used to re-assign original value to flag_xc at the end of the subroutine
      integer :: flag_xc_temp

!   only referenced and allocated if gradient correction requested
      real*8, dimension(:,:,:), allocatable :: dylm_dtheta_tab
      real*8, dimension(:,:,:), allocatable :: scaled_dylm_dphi_tab

      real*8 :: coord_current(3)
      real*8 :: dist_tab(n_atoms,n_max_batch_size)
      real*8 :: i_r(n_atoms, n_max_batch_size)
      real*8 :: dir_tab(3,n_atoms,n_max_batch_size)
      real*8 :: trigonom_tab(4,n_atoms, n_max_batch_size)

      real*8 :: fxc_times_psi(n_basbas, n_max_batch_size,n_spin)
      real*8 :: radial_wave(n_basbas, n_max_batch_size)
      real*8 :: wave(n_basbas, n_max_batch_size)

      integer :: n_compute
      integer :: i_basis(n_basbas)

!     optimal accounting for matrix multiplications: only use points with nonzero components
      integer :: n_points

!     and condensed version of partition_tabs on angular grids
      real*8 :: partition(n_max_batch_size)

      real*8, dimension(:,:,:), allocatable :: gradient_basis_wave

      real*8 :: en_density_xc
      real*8, dimension(n_spin) :: en_density_x
      real*8 :: en_density_c
      real*8 :: local_xc_density_deriv(n_spin)
      real*8 :: local_xc_gradient_deriv(3,n_spin)

      type(xc_f03_func_t) :: xc_func_x, xc_func_c 
      type(xc_f03_func_info_t) :: xc_info_x, xc_info_c
      real*8 :: f_xc(1), f_x(1), f_c(1)!, grad_contracted(1), dv2dnds(1), dv2ds2(1)
      ! dv2dnds and dv2ds2 are dummy variables for unused output from libxc (GGA's only)
      ! look them up at the website if you are interested in what they mean

!  counters
      integer :: i_my_batch
      integer :: i_index
      integer :: i_l, i_m
      integer :: i_coord
      integer :: i_spin
      integer :: i_compute
      integer :: i_point
      integer :: i_full_points
      integer :: i_full_points_2
      integer :: k_cnt, basbas_l_max(n_species)

      real*8, dimension(n_spin, n_full_points)  :: rho_inc_partialcore
      real*8, dimension(3, n_spin, n_full_points) :: rho_gradient_inc_partialcore

!  begin work

      if(myid.eq.0) then
        write(use_unit,'(2X,A,2X,A)') &
          "| Integrating  < u | f_xc | v >  in auxiliary basis..."
!        write(use_unit,'(2x,a)') '|'
      endif

      flag_xc_temp = flag_xc

      if(use_libxc_tddft) then
         call xc_f03_func_init(xc_func_x, libxc_tddft_x, XC_UNPOLARIZED)
         xc_info_x = xc_f03_func_get_info(xc_func_x)
         call xc_f03_func_init(xc_func_c, libxc_tddft_c, XC_UNPOLARIZED)
         xc_info_c = xc_f03_func_get_info(xc_func_c)
      else
        flag_xc = flag_fxc
      endif

!     begin with general allocations

      if ( use_density_gradient ) then
        l_ylm_max = l_wave_max
        allocate (gradient_basis_wave(n_basbas,3,n_max_batch_size))
        allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_atoms, &
             n_max_batch_size) )
        allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_atoms, &
             n_max_batch_size) )
      end if

      if(use_embedding_pp.and.use_nonlinear_core) then
         do i_spin = 1,n_spin
            rho_inc_partialcore(i_spin,:) = rho(i_spin,:) + partial_core_rho(:)
            rho_gradient_inc_partialcore(:,i_spin,:) = rho_gradient(:,i_spin,:)
            if(use_density_gradient) then
               rho_gradient_inc_partialcore(:,i_spin,:) = &
                  rho_gradient_inc_partialcore(:,i_spin,:) + partial_core_rho_grad(:,:)
            endif
         enddo
      endif

      ! l_max for ylm in product basis set to 2 * l_max for wave basis 
      l_ylm_max = 2 * l_wave_max
      allocate( ylm_tab( (l_ylm_max+1)**2, n_atoms, &
           n_max_batch_size) )
      allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max) )

!     initialize
      basbas_l_max = 2 * basis_l_max

      ! initialize f_xc(n_basbas,n_basbas,n_spin) with zeros for security
      fxc_matr(:,:,:) = 0d0

!     initialize index_lm
      i_index = 0
      do i_l = 0, l_ylm_max, 1
        do i_m = -i_l, i_l
          i_index = i_index+1
          index_lm(i_m,i_l) = i_index
        enddo
      enddo

      i_full_points_2 = 0
      i_full_points = 0
!     perform partitioned integration, atom by atom, and point by point
!     This will be the outermost loop, to save evaluations of the potential.
!     and the Y_lm functions

      do  i_my_batch = 1, n_my_batches, 1

           n_compute = 0
           i_basis = 0

           i_point = 0

           do i_index = 1, batches(i_my_batch)%size, 1

             i_full_points_2 = i_full_points_2 + 1

             if(partition_tab(i_full_points_2) .gt. 0.d0) then

              i_point = i_point + 1

!     get current integration point coordinate
              coord_current(:) = batches(i_my_batch) &
                                 % points(i_index) % coords(:)

!     compute atom-centered coordinates of current integration point,
!     as viewed from all atoms
              call tab_atom_centered_coords &
                   ( coord_current(1), &
                     dist_tab(1,i_point), i_r(1,i_point), &
                     dir_tab(1,1,i_point) &
                   )
!     TODO
!     determine which basis functions are relevant at current integration point,
!     and tabulate their indices --------- LATER

!             call prune_basis_v1(dist_tab(1,i_point), n_compute, &
!                      i_basis)


!call prune_prodbas_v1( myid+1, dist_tab, n_compute, i_basis)

             endif
           enddo

           n_points = i_point


!write(use_unit,*) n_compute, n_basbas

!           if (n_compute .gt. 0) then

           n_compute = n_basbas

           do k_cnt=1, n_basbas
              i_basis(k_cnt) = k_cnt
           enddo

           i_point = 0
           fxc_times_psi =0d0

           do i_index = 1, batches(i_my_batch)%size, 1

              i_full_points = i_full_points + 1
              if (partition_tab(i_full_points).gt.0d0) then
!     execute only if partition_tab.gt.0 here, i.e. if the integration point
!     makes sense
                 i_point = i_point + 1
                 partition(i_point) = partition_tab(i_full_points)
!     compute trigonometric functions of spherical coordinate angles
!     of current integration point, viewed from all atoms
                 call tab_trigonom &
                         ( dir_tab(1,1,i_point), &
                         trigonom_tab(1,1,i_point) &
                         )

!     tabulate distance and Ylm's w.r.t. other atoms
                 call tab_wave_ylm &
                         ( trigonom_tab(1,1,i_point), basbas_l_max, &
                         l_ylm_max, &
                         ylm_tab(1,1,i_point) )

!           tabulate total wave function value for each basis function
                 call evaluate_prod_waves &
                         (i_r(1,i_point), l_ylm_max, &
                         ylm_tab(1,1,i_point), &
                         dist_tab(1,i_point), index_lm, &
                         n_compute, &
                         i_basis, &
                         wave(1,i_point))

!                     if(use_embedding_pp.and.use_nonlinear_core) then

!                        call evaluate_xc &
!                            ( rho_inc_partialcore(1:n_spin,i_full_points), &
!                            rho_gradient_inc_partialcore(1:3,1:n_spin, i_full_points), &
!                           en_density_xc, &
!                           en_density_x, en_density_c,  &
!                           local_xc_density_deriv(1), &
!                           local_xc_gradient_deriv(1,1) &
!                           )

!                     else


                  if(use_libxc_tddft) then
                     call xc_f03_lda_fxc(xc_func_x, 1, rho(n_spin,i_full_points), f_x(1))
                     call xc_f03_lda_fxc(xc_func_c, 1, rho(n_spin,i_full_points), f_c(1))
                     f_xc(1) = f_x(1) + f_c(1)
                  else
                     call evaluate_xc_shanghui &
                            ( rho(1:n_spin,i_full_points), &
                            rho_gradient(1:3,1:n_spin, i_full_points), &
                           en_density_xc, &
                           en_density_x, en_density_c,  &
                           local_xc_density_deriv(1), &
                           local_xc_gradient_deriv(1,1), &
                           f_xc(n_spin))
                  endif


!     evaluate f_xc * psi in this particular point
                  do i_spin = 1, n_spin
                     fxc_times_psi(1:n_compute,i_point,i_spin) = f_xc(1) * wave(1:n_compute,i_point)
                  enddo

!     end if (partition_tab.gt.0)
             end if

!     end angular integration loop
          end do
!     add non-relativistic contributions to the matrix elements

          do i_spin = 1, n_spin
            call evaluate_fxc_matr_shell &
               ( i_point, partition(1:i_point), &
                 n_compute, i_basis(1:n_compute), &
                 fxc_times_psi(1,1,i_spin), &
                 wave(1,1), fxc_matr(1,1,i_spin) )
         end do

!       else
!           i_full_points = i_full_points + &
!                batches(i_my_batch)%size

!     end if (n_compute.gt.0) then
 !      end if


!     end integration loop over batches
      end do

      if(use_libxc_tddft) then
         call xc_f03_func_end(xc_func_x)
         call xc_f03_func_end(xc_func_c)
      endif

      do i_spin = 1, n_spin
         call sync_matrix(fxc_matr(:,:,i_spin),n_basbas,n_basbas)
      enddo

      flag_xc = flag_xc_temp

      if (allocated(ylm_tab)) then
         deallocate(ylm_tab)
      end if
      if (allocated(index_lm)) then
         deallocate(index_lm)
      end if

      return
end subroutine integrate_fxc_matrix





subroutine evaluate_fxc_matr_shell &
           ( n_points, partition, &
             n_compute, i_basis, fxc_times_psi, &
             wave, fxc_matr &
           )

!  PURPOSE
!  Subroutine evaluate_fxc_matr_shell evaluates the final matrix resulting from the 
!  integrate_fxc_matrix routine just above.
!
!  USES

      use dimensions
      use mpi_tasks
      implicit none

!  ARGUMENTS

      integer :: n_points

      integer n_compute
      integer i_basis(n_compute)

      real*8 partition(n_points)
      real*8 wave(n_basbas,n_points)
      real*8 fxc_times_psi(n_basbas,n_points)

      real*8 fxc_matr( n_basbas,n_basbas )

!  INPUTS
!  o n_points -- number of relevant points in the current integration shell 
!  o n_compute -- number of non-zero basis functions in the current integration shell
!  o i_basis -- array, specifies the non-zero basis functions  
!  o partition -- the values of the partition function (for integration) at each 
!         spatial grid point
!  o wave -- value of the basis functions
!  o fxc_times_psi -- the local f_xc times basis functions
!
!  OUTPUTS
!  o fxc_matr -- the full <u|f_xc|v> matrix in the auxiliary basis
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

      real*8,allocatable,dimension(:,:) :: aux_fxc_matrix, wave_compute

!     counters
      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_index
      integer :: i_compute
      integer :: i_compute_1
      integer :: i_index_real
      integer :: i_point

!     begin work

      allocate(wave_compute(n_compute,n_points),stat=i_point)
      call check_allocation(i_point, 'wave_compute')
      allocate(aux_fxc_matrix(n_compute,n_compute),stat=i_point)
      call check_allocation(i_point, 'aux_fxc_matrix')
      aux_fxc_matrix = 0d0

      do i_point = 1, n_points, 1

         wave_compute(1:n_compute, i_point) = &
              wave(1:n_compute, i_point)*partition(i_point)

      enddo

!     compute wave*(f_xc*psi) and add this to aux. overlap matrix
      call dgemm('N', 'T', n_compute, n_compute, &
           n_points, 1.0d0, &
           wave_compute, n_compute, fxc_times_psi, &
           n_basbas, 0.0d0, aux_fxc_matrix, &
           n_compute )

!  now add the aux. overlap matrix to the actual overlap matrix
!  this requires translating between the actually computed matrix elements and the
!  full overlap matrix ...
      if (use_density_gradient) then

        do i_compute = 1, n_compute
            do i_compute_1 = 1, i_compute -1

             aux_fxc_matrix(i_compute, i_compute_1) = &
             0.5d0 * ( aux_fxc_matrix(i_compute, i_compute_1) + &
                       aux_fxc_matrix(i_compute_1, i_compute) )


             aux_fxc_matrix(i_compute_1, i_compute) = aux_fxc_matrix(i_compute, i_compute_1)

            enddo
         enddo

      endif

      i_index = 0
      i_index_real = 0

      do i_compute =1, n_compute, 1

         i_basis_1=i_basis(i_compute)
         i_index=  (i_basis(i_compute)-1)*i_basis(i_compute)/2

         do i_compute_1 = 1, n_compute, 1

              i_basis_2=i_basis(i_compute_1)

              i_index_real =  i_index + &
                   i_basis (i_compute_1)

              fxc_matr(i_basis_1,i_basis_2) = &
                fxc_matr(i_basis_1,i_basis_2) &
               + aux_fxc_matrix(i_compute, i_compute_1)

!     end of i_compute_1
           enddo

!     end of i_compute
       enddo
end subroutine evaluate_fxc_matr_shell
