!****s*  FHI-aims/integrate_xc_matrix
!  NAME
!    integrate_xc_matrix
!  SYNOPSIS
     subroutine integrate_xc_matrix &
     ( partition_tab, &
       rho, rho_gradient, &
       kinetic_density, &
       basis_l_max, &
       xc_matr )

!  PURPOSE
!
!  This subroutine is intended to calculate the exchange-correlaation energy
!  matrix elements between every two basis functions, need for G0W0 quasiparticle
!  calculations
!  
!  USES

      use dimensions
      use runtime_choices
      use grids
      use geometry
      use basis
      use constants
      use mpi_tasks
      use mpi_utilities
      use synchronize_mpi
      use pseudodata

      implicit none

!  ARGUMENTS

      real*8, dimension(n_full_points) :: partition_tab
      real*8, dimension(n_spin,n_full_points) :: rho
      real*8, dimension(3,n_spin, n_full_points) :: rho_gradient
      real*8, dimension(n_spin,n_full_points) :: kinetic_density
      integer basis_l_max (n_species)
      real*8 xc_matr( n_basis,n_basis,n_spin)

! INPUTS
!  o partition_tab -- real array, the values of the partition function (for 
!           integration) at each spatial grid point    
!  o rho   -- real array, the electron density at each grid point and for each spin 
!           channel 
!  o rho_gradient -- real array, the electron density gradient
!  o kinetic_density -- real array, the kinetic-energy density
!  o basis_l_max -- the maximal angualr momentum number of the basis function for
!           each species  
! OUTPUTS
!  o xc_matr -- the matrix elements within the regular basis functions of the exchange
!           correlation potential
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

!   only referenced and allocated if gradient correction requested
      real*8, dimension(:,:,:), allocatable :: dylm_dtheta_tab
      real*8, dimension(:,:,:), allocatable :: scaled_dylm_dphi_tab


      real*8 coord_current(3)
      real*8 dist_tab(n_atoms,n_max_batch_size)
      real*8 i_r(n_atoms, n_max_batch_size)
      real*8 dir_tab(3,n_atoms,n_max_batch_size)
      real*8 trigonom_tab(4,n_atoms, n_max_batch_size)

      real*8 xc_times_psi(n_basis, n_max_batch_size,n_spin)
      real*8 radial_wave(n_basis, n_max_batch_size)
      real*8 wave(n_basis, n_max_batch_size)

      integer :: n_compute
      integer :: i_basis(n_basis)

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
      real*8 :: local_xc_tau_deriv(n_spin)

      real time_stamp
      real time_end

!  counters

      integer i_my_batch
      integer i_grid
      integer i_index
      integer i_l, i_m
      integer i_coord
      integer i_spin
      integer i_compute

      integer i_point

      integer i_full_points
      integer i_full_points_2

      integer i_basis_1
      integer i_basis_2

      real*8, dimension(:,:), allocatable  :: rho_inc_partialcore
      real*8, dimension(:,:,:), allocatable  :: rho_gradient_inc_partialcore

!  begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A,2X,A)') &
          "Integrating the xc potential matrix", &
          "for basis functions ..."
      endif

!     begin with general allocations
        l_ylm_max = l_wave_max

      allocate( ylm_tab( (l_ylm_max+1)**2, n_atoms, &
           n_max_batch_size) )
      allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max) )

      if ( use_density_gradient ) then
        l_ylm_max = l_wave_max
        allocate (gradient_basis_wave(n_basis,3,n_max_batch_size))
        allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_atoms, &
             n_max_batch_size) )
        allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_atoms, &
             n_max_batch_size) )
      end if

      if(use_embedding_pp.and.use_nonlinear_core) then
         allocate(rho_inc_partialcore(n_spin,n_full_points))
         allocate(rho_gradient_inc_partialcore(3,n_spin,n_full_points))
         do i_spin = 1,n_spin
           rho_inc_partialcore(i_spin,:) = rho(i_spin,:) + partial_core_rho(:)
           rho_gradient_inc_partialcore(:,i_spin,:) = rho_gradient(:,i_spin,:)
           if(use_density_gradient) then
              rho_gradient_inc_partialcore(:,i_spin,:) = &
                 rho_gradient_inc_partialcore(:,i_spin,:) + partial_core_rho_grad(:,:)
           endif
        enddo
      endif

!     initialize

      i_index = 0
!      do i_basis_1 = 1, n_basis, 1
!          do i_basis_2 = 1, i_basis_1, 1
!            i_index = i_index+1

!            xc_matr(i_index) = 0.
!        enddo
!      enddo
       xc_matr(:,:,:)=0.d0

!     initialize index_lm

      i_index = 0
      do i_l = 0, l_wave_max, 1
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
!      do i_atom = 1, n_atoms, 1

!test
!        write(use_unit,*) "i_atom: ", i_atom
!test end

!         do i_radial = 1, n_radial(species(i_atom)), 1

!  mpi task distribution
!          if (myid.eq.radial_task_list(i_radial,i_atom)) then

!test
!          write(use_unit,*) "  i_radial: ", i_radial
!test end

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

!     determine which basis functions are relevant at current integration point,
!     and tabulate their indices

             call prune_basis_v1(dist_tab(1,i_point), n_compute, &
                      i_basis)

!     evaluate the partition function for three atoms
!              call evaluate_partition_tab_2atoms
!     +             (i_atom, dist_tab(1:n_atoms,i_angular),
!     +              i_r(1:n_atoms,i_angular),
!     +              w_radial( i_radial, species (i_atom)),
!     +              w_angular( i_angular, i_radial, species (i_atom)),
!     +              partition_tab_2atoms(1, i_angular),
!     +              partition_type )
!
!
            endif
           enddo

           n_points = i_point
!test
!     if (i_radial.eq.140) then
!     write(use_unit,*) "here 2"
!     end if
!test end

           if (n_compute .gt. 0) then

              i_point = 0
              xc_times_psi =0.d0

!              if(use_density_gradient) then
!                 gradient_basis_wave(:,:,:) = 0.d0
!              endif

              do i_index = 1, batches(i_my_batch)%size, 1

                 i_full_points = i_full_points + 1
                 if (partition_tab(i_full_points).gt.0.d0) &
                      then
!     execute only if partition_tab.gt.0 here, i.e. if the integration point
!     makes sense
                    i_point = i_point + 1
                    partition(i_point) = &
                         partition_tab(i_full_points)
!     compute trigonometric functions of spherical coordinate angles
!     of current integration point, viewed from all atoms
                    call tab_trigonom &
                         ( dir_tab(1,1,i_point), &
                         trigonom_tab(1,1,i_point) &
                         )

!     tabulate distance and Ylm's w.r.t. other atoms
                    call tab_wave_ylm &
                         ( trigonom_tab(1,1,i_point), basis_l_max, &
                         l_ylm_max, &
                         ylm_tab(1,1,i_point) )

!           tabulate total wave function value for each basis function
                    call evaluate_waves_v0 &
                         (i_r(1,i_point), l_ylm_max, &
                         ylm_tab(1,1,i_point), &
                         dist_tab(1,i_point), index_lm, n_compute, &
                         i_basis, radial_wave(1,i_point), &
                         wave(1,i_point))


                    if(use_embedding_pp.and.use_nonlinear_core) then

                        call evaluate_xc &
                            ( rho_inc_partialcore(1:n_spin,i_full_points), &
                            rho_gradient_inc_partialcore(1:3,1:n_spin, i_full_points), &
                            kinetic_density(1:n_spin, i_full_points), &
                            en_density_xc, en_density_x, en_density_c,  &
                            local_xc_density_deriv(1), &
                            local_xc_gradient_deriv(1,1), &
                            local_xc_tau_deriv(1), &
                            .false. &
                            )

                     else

                        call evaluate_xc &
                            ( rho(1:n_spin,i_full_points), &
                            rho_gradient(1:3,1:n_spin, i_full_points), &
                            kinetic_density(1:n_spin, i_full_points), &
                            en_density_xc, en_density_x, en_density_c,  &
                            local_xc_density_deriv(1), &
                            local_xc_gradient_deriv(1,1), &
                            local_xc_tau_deriv(1), &
                            .false. &
                            )

                     endif


!     evaluate V_xc * psi in this particular point
                    do i_spin = 1, n_spin
                       xc_times_psi (1:n_compute,i_point,i_spin)= &
                         local_xc_density_deriv (i_spin) * &
                         wave(1:n_compute,i_point)
                    enddo


                   if (use_density_gradient) then

!     tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
                      call tab_gradient_ylm &
                         ( trigonom_tab(1,1,i_point), basis_l_max, &
                           l_ylm_max, ylm_tab(1,1,i_point), &
                           dylm_dtheta_tab(1,1,i_point), &
                           scaled_dylm_dphi_tab(1,1,i_point) &
                         )

                      call evaluate_wave_gradient_v1 &
                         ( dist_tab(1,i_point), i_r(1,i_point), &
                           dir_tab(1,1,i_point), &
                           trigonom_tab(1,1,i_point), &
                           l_ylm_max, ylm_tab(1,1,i_point), &
                           dylm_dtheta_tab(1,1,i_point), &
                           scaled_dylm_dphi_tab(1,1,i_point), &
                           index_lm, n_compute, &
                           i_basis(1:n_compute), &
                           radial_wave(1, i_point), &
                           gradient_basis_wave (1,1,i_point) &
                         )

! add gradient part to xc * psi

                      do i_spin = 1, n_spin, 1
                        do i_coord = 1, 3, 1
                          do i_compute = 1, n_compute, 1

                           xc_times_psi(i_compute, i_point, i_spin) = &
                            xc_times_psi(i_compute, i_point, i_spin) + &
                            4.0d0* &
                            local_xc_gradient_deriv(i_coord, i_spin) * &
!     +                      rho_gradient(i_coord, i_spin,i_full_points)*
                            gradient_basis_wave(i_compute,i_coord, &
                                                i_point)

                         enddo
                       enddo
                     enddo

!   end if (use_density_gradient)
                   end if

!     end if (partition_tab.gt.0)
             end if

!     end angular integration loop
          enddo
!     add non-relativistic contributions to the Hamiltonian matrix elements

          do i_spin = 1, n_spin
            call evaluate_xc_matr_shell &
               ( n_points, partition(1:n_points), &
                 n_compute, i_basis(1:n_compute), &
                 xc_times_psi(1,1,i_spin), &
                 wave(1,1), xc_matr(1,1,i_spin) )
         enddo

       else
           i_full_points = i_full_points + &
                batches(i_my_batch)%size

!     end if (n_compute.gt.0) then
       end if


!test
!            if (i_radial.eq.140) then
!            write(use_unit,*) "here 3"
!            end if
!test end

!test
!            if ( (i_atom.eq.1) .and.
!     +           ((i_radial.eq.47).or.(i_radial.eq.48)) .and.
!     +           (i_angular.eq.1) ) then
!              write(use_unit,*) "      n_compute: ", n_compute
!              write(use_unit,*) "      i_basis  : ",
!     +          (i_basis(i_compute), i_compute = 1,n_compute,1)
!            end if

!       end of mpi distribution
         ! endif
!test end

!     end integration loop over batches
      enddo

      do i_spin = 1, n_spin
        call sync_matrix &
            (xc_matr(:,:,i_spin),n_basis,n_basis)
      enddo
!test
!      i_index = 0
!      write(use_unit,'(2X,A,2X,A,7X,A)')
!     +          "i_basis_1", "i_basis_2", "kinetic energy matrix"
!
!      if(myid.eq.0) then
!       do i_basis_1 =1, n_basis
!         do i_basis_2 = 1, n_basis
!
!         write(use_unit, '(4X,I5,5X,I5,9X,f16.10,f16.10)') &
!                     i_basis_1, i_basis_2, &
!                     xc_matr(i_basis_1,i_basis_2, 1), &
!                     xc_matr(i_basis_2,i_basis_1, 1)
!
!         enddo
!              write(use_unit, *)
!       enddo
!      endif
!test

      if (allocated(ylm_tab)) then
        deallocate(ylm_tab)
      end if
      if (allocated(index_lm)) then
        deallocate(index_lm)
      end if

      if(allocated( rho_inc_partialcore  )) then 
         deallocate( rho_inc_partialcore  )
      endif

      if(allocated( rho_gradient_inc_partialcore )) then 
         deallocate( rho_gradient_inc_partialcore )
      endif


      return
      end subroutine integrate_xc_matrix

!----------------------------------------------------------------------
!******
