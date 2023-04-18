!****s* FHI-aims/broyden_mix_p1
!  NAME
!   broyden_mix_p1
!  SYNOPSIS

subroutine broyden_mix_p1 &
     ( partition_tab, hartree_partition_tab, &
     saved_iter, delta_rho_KS, rho_diff, &
     previous_rho_diff, previous_rho_error, &
     delta_rho_gradient, rho_gradient_diff, &
     previous_rho_gradient_diff, &
     previous_rho_gradient_error, &
     delta_kinetic_density, kinetic_density_diff, &
     previous_kinetic_density_diff, &
     previous_kinetic_density_error, &
     broyden_matrix, mixing_factor )

!  PURPOSE
!  Subroutine broyden_mix
!
!  Use broyden's algorithm to mix present density with the densities 
!  of the previous iterations. 
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use basis
  use species_data
  use xc
  use lapack_wrapper
  use localorb_io
  use mpi_utilities
  use synchronize_mpi
  use constants
  use precondition 
  use numerical_utilities
  use lpb_solver_utilities, only : atomic_MERM, mpb_solver_started
  implicit none

!  ARGUMENTS

  real*8, dimension(n_full_points) :: partition_tab
  real*8, dimension(n_full_points) :: hartree_partition_tab
  integer :: saved_iter

  real*8, dimension(n_full_points, n_spin)              :: delta_rho_KS
  real*8, dimension(n_full_points, n_spin)              :: rho_diff
  real*8, dimension(n_full_points, n_max_broyden, n_spin) :: previous_rho_diff
  real*8, dimension(n_full_points, n_max_broyden, n_spin) :: previous_rho_error
  
  real*8, dimension(3, n_full_points, n_spin)              :: delta_rho_gradient
  real*8, dimension(3, n_full_points, n_spin)              :: rho_gradient_diff
  real*8, dimension(3, n_full_points, n_max_broyden, n_spin) :: previous_rho_gradient_diff
  real*8, dimension(3, n_full_points, n_max_broyden, n_spin) :: previous_rho_gradient_error

  real*8, dimension(n_full_points, n_spin)              :: delta_kinetic_density
  real*8, dimension(n_full_points, n_spin)              :: kinetic_density_diff
  real*8, dimension(n_full_points, n_max_broyden, n_spin) :: previous_kinetic_density_diff
  real*8, dimension(n_full_points, n_max_broyden, n_spin) :: previous_kinetic_density_error

  real*8, dimension(n_max_broyden,n_max_broyden) :: broyden_matrix
  real*8, dimension(saved_iter)        :: mixing_factor

!  INPUTS
!  o partition_tab -- values of the partition function
!  o hartree_partition_tab -- values of the partition function for Hartree potential
!  o saved_iter -- number of saved iterations
!  o delta_rho_KS -- electron density - free atoms electron density
!  o rho_diff -- change in electron density during sc iteration
!  o previous_rho_diff -- previous delta charge density
!  o previous_rho_error -- ?????????????

!  Only referenced and allocated if gradient functional requested:
!  o delta_rho_gradient -- gradient of delta charge
!  o rho_gradient_diff -- change in gradient of electron density
!  o previous_rho_gradient_diff -- changes from previous iterations in gradient of electron density
!  o previous_rho_gradient_error -- ?????????????
!
!  Only referenced and allocated if meta-gga functional requested:
!  o delta_kinetic_density -- kinetic density
!  o kinetic_density_diff -- change in kinetic density
!  o previous_kinetic_density_diff -- changes from previous iterations in kinetic density 
!  o previous_kinetic_density_error -- ?????????????
!
!  OUTPUT
!  o broyden_matrix -- ????????
!  o mixing_factor -- ???????
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

  real*8, dimension(saved_iter+1) :: broyden_load
  real*8, dimension(saved_iter+1,saved_iter) :: broyden_load_basis

! Temp variables
    real*8, dimension(saved_iter,saved_iter+1) :: tmp_matrix_1
    real*8, dimension(saved_iter,saved_iter) :: tmp_matrix_2
    real*8, dimension(saved_iter) :: tmp_vector
    real*8 :: sigma

! counters

  integer :: i_index
  integer :: i_store, i_store_2
  integer :: i_coord
  integer :: i_offset
  integer :: i_my_batch
  integer :: i_spin

!  begin work

  sigma = relative_fp_charge_mix

  call localorb_info( &
       "Broyden mixing of updated and previous charge densities.", &
       use_unit,'(2X,A)', OL_norm )

!     initialize
  mixing_factor = 0.d0
  broyden_load = 0.d0
  broyden_load_basis = 0.d0
  tmp_matrix_1 = 0.d0
  tmp_matrix_2 = 0.d0
  tmp_vector = 0.d0

  broyden_load(1) = 1.d0
  do i_store = 1, saved_iter
      broyden_load_basis(i_store,i_store) = 1.d0
      broyden_load_basis(i_store+1,i_store) = -1.d0
  enddo

!    Prepare updated broyden matrix
  do i_store = saved_iter+1, 2, -1
     do i_store_2 = saved_iter+1, i_store, -1
        broyden_matrix(i_store, i_store_2) =  &
             broyden_matrix(i_store-1, i_store_2-1)
        broyden_matrix(i_store_2, i_store) =  &
             broyden_matrix(i_store, i_store_2)
     enddo
  enddo
!     initialize broyden matrix elements of present iteration
  do i_store = 1, saved_iter+1, 1
     broyden_matrix(1,i_store) = 0.d0      
     broyden_matrix(i_store,1) = 0.d0      
  enddo

!     calculate, use, and store density error err(rho_i) for present iteration i
!     between KS input and output densities from previous iterations n-i, n-j ...
  i_offset = 0

  do i_my_batch = 1, n_my_batches, 1

        do i_index = 1, batches(i_my_batch)%size, 1
           
           i_offset = i_offset + 1
!           execute only if partition_tab.gt.0 here, i.e. if the integration point
!           makes sense
           if (partition_tab(i_offset).gt.0.d0) then

!             determine broyden matrix elements: < err(rho_m+1) - err(rho_m) | err(rho_j+1) - err(rho_j) >
!             determine broyden vector: < err(rho_j+1) - err(rho_j) | err(rho_m) >
              
              do i_spin = 1, n_spin, 1
                  broyden_matrix(1,1) =  &
                      broyden_matrix(1,1) + &
                      partition_tab(i_offset) * &
                      delta_rho_KS(i_offset, i_spin) * &
                      delta_rho_KS(i_offset, i_spin)
              end do

              do i_store = 1, saved_iter, 1
                 
                 do i_spin = 1, n_spin, 1

                    broyden_matrix(1,i_store+1) =  &
                         broyden_matrix(1,i_store+1) + &
                         partition_tab(i_offset) * &
                         delta_rho_KS(i_offset, i_spin) * &
                         previous_rho_error(i_offset, i_store, i_spin)
                 
                 end do

                 broyden_matrix(i_store+1,1) = broyden_matrix(1,i_store+1)
                 
              end do
                 
           end if

           !         end loop over a batch
        end do
     !       end loop over batches
  end do

  call sync_broyden_matrix(broyden_matrix)

!     Determine mixing factor from equation system defined by broyden_matrix
!     use the lapack expert solver for a symmetric indefinite problem
  if (saved_iter.gt.0) then

     ! It is possible that two successive iterations have a zero residual.
     ! In this case, we simply remove the appropriate column(s) / row(s) 
     ! from the broyden matrix. This should give a zero mixing factor for the 
     ! density from that iteration. Is that the right thing to do?
     do i_store = 1, saved_iter, 1
        if (broyden_matrix(i_store,i_store).le.1.d-20) then
           broyden_matrix(:,i_store) = 0.d0
           broyden_matrix(i_store,:) = 0.d0
           broyden_matrix(i_store,i_store) = 1.d0
        end if
     enddo

  end if

  if (saved_iter .gt. 0) then
      tmp_matrix_1 = transpose(broyden_load_basis)
      tmp_matrix_1 = matmul(tmp_matrix_1, &
          broyden_matrix(1:(saved_iter+1),1:(saved_iter+1)))
      tmp_matrix_2 = matmul(tmp_matrix_1, broyden_load_basis)
      tmp_vector = matmul(tmp_matrix_1, broyden_load)
      call solve_LEQ("broyden_mix_p1", saved_iter, &
          tmp_matrix_2, tmp_vector, mixing_factor)
  end if

  ! put together the next density in two pieces:
  ! (1) the sum \Delta R_i terms 
  i_offset   = 0
  do i_my_batch = 1, n_my_batches, 1
        do i_index = 1, batches(i_my_batch)%size, 1
           i_offset = i_offset + 1
           if (partition_tab(i_offset).gt.0.d0.or.(atomic_MERM .and. mpb_solver_started)) then
              do i_spin = 1, n_spin, 1

              ! TODO Here is the linear mixer independent of the history.
              ! BUT: delta_rho_KS for the next iteration depends on what
              ! is in the history even though this rho_diff is independent
              ! of it.
                rho_diff(i_offset, i_spin) = sigma * charge_mix_param(i_spin) * &
                delta_rho_KS(i_offset, i_spin)

                ! - + -
                if (saved_iter .gt. 0) then
                    rho_diff(i_offset,i_spin) = rho_diff(i_offset,i_spin) - &
                    mixing_factor(1) * charge_mix_param(i_spin) * &
                    ( previous_rho_diff(i_offset,1,i_spin) + &
                    sigma * (delta_rho_KS(i_offset,i_spin) - &
                    previous_rho_error(i_offset,1,i_spin)))
                end if

                do i_store = 2, saved_iter, 1
                    rho_diff(i_offset, i_spin) = rho_diff(i_offset, i_spin) - &
                    mixing_factor(i_store) * charge_mix_param(i_spin) * &
                    ( previous_rho_diff(i_offset,i_store,i_spin) + &
                    sigma * (previous_rho_error(i_offset,i_store-1,i_spin) - &
                    previous_rho_error(i_offset,i_store,i_spin)))
                end do

              end do
           end if
        end do          ! end loop over a batch
  end do                ! end loop over batches

  ! ... These sum \Delta R terms should potentially be preconditioned after they've all been put together 
  if (use_kerker_preconditioner.and.kerker_preconditioner_on .and. &
  &   .not. precondition_before_mixer) then
     do i_spin = 1, n_spin
        call precondition_kerker(rho_diff(:,i_spin), hartree_partition_tab)
     end do
  end if

!    if required, play exactly the same mixing game for the gradients

!     FIXME: this whole sequence needs a subroutine of its own to
!      avoid the silly code duplication between density and gradient; 
!      i.e. we could split the
!      update of the broyden matrix and the actual broyden mixing into
!     completely separate subroutines ...

  if (use_density_gradient) then
      ! We need to store gradients for both residuals and positions for this to
      ! work.
     ! add together all pieces that might require preconditioning
     i_offset = 0
     do i_my_batch = 1, n_my_batches, 1
           do i_index = 1, batches(i_my_batch)%size, 1
              i_offset = i_offset + 1
              if (partition_tab(i_offset).gt.0.d0.or.(atomic_MERM .and. mpb_solver_started)) then
                 do i_coord = 1,3,1
                    do i_spin = 1, n_spin, 1
                       rho_gradient_diff(i_coord,i_offset, i_spin) = &
                            sigma * charge_mix_param(i_spin) * delta_rho_gradient(i_coord, i_offset, i_spin)
                       if (saved_iter.gt.0) then
                          rho_gradient_diff(i_coord, i_offset, i_spin) = &
                          rho_gradient_diff(i_coord, i_offset, i_spin) - &
                          mixing_factor(1) * charge_mix_param(i_spin) * &
                          ( previous_rho_gradient_diff(i_coord,i_offset,1,i_spin) + &
                          sigma * (delta_rho_gradient(i_coord,i_offset,i_spin) - &
                          previous_rho_gradient_error(i_coord,i_offset,1,i_spin)))
                          !( previous_rho_gradient_error(i_coord,i_offset,1,i_spin) + &
                          !sigma * (delta_rho_gradient(i_coord,i_offset,i_spin) - &
                          !previous_rho_gradient_diff(i_coord,i_offset,1,i_spin)))
                       end if
                       do i_store = 2, saved_iter, 1
                          rho_gradient_diff(i_coord, i_offset, i_spin) = &
                          rho_gradient_diff(i_coord, i_offset, i_spin) - &
                          mixing_factor(i_store) * charge_mix_param(i_spin) * &
                          ( previous_rho_gradient_diff(i_coord,i_offset,i_store,i_spin) + &
                          sigma * (previous_rho_gradient_error(i_coord,i_offset,i_store-1,i_spin) - &
                          previous_rho_gradient_error(i_coord,i_offset,i_store,i_spin)))
                          !( previous_rho_gradient_error(i_coord,i_offset,i_store,i_spin) + &
                          !sigma * (previous_rho_gradient_diff(i_coord,i_offset,i_store-1,i_spin) - &
                          !previous_rho_gradient_diff(i_coord,i_offset,i_store,i_spin)))
                       enddo
                    end do
                 enddo     !     next coordinate [of x,y,z]
              end if
           end do          !     end loop over a batch
     end do                !     end loop over batches

     ! precondition gradient, if required
     if (use_kerker_preconditioner.and.kerker_preconditioner_on .and. &
     &   .not.precondition_before_mixer) then
        do i_coord = 1, 3
           do i_spin = 1, n_spin
              call precondition_kerker(rho_gradient_diff(i_coord,:,i_spin), hartree_partition_tab)
           end do
        end do
     end if

     ! add remaining terms
!!$     i_offset = 0
!!$     do i_my_batch = 1, n_my_batches, 1
!!$           do i_index = 1, batches(i_my_batch)%size, 1
!!$              i_offset = i_offset + 1
!!$              if (partition_tab(i_offset).gt.0.d0) then
!!$                 do i_coord = 1,3,1
!!$                    do i_spin = 1, n_spin, 1
!!$                       if (saved_iter.gt.0) then
!!$                          rho_gradient_diff(i_coord, i_offset, i_spin) = &
!!$                               rho_gradient_diff(i_coord, i_offset, i_spin) + &
!!$                               mixing_factor(1) * &
!!$                               previous_rho_gradient_diff(i_coord, i_offset,1, i_spin)
!!$                       end if
!!$                       do i_store = 2, saved_iter, 1
!!$                          rho_gradient_diff(i_coord, i_offset, i_spin) = &
!!$                               rho_gradient_diff(i_coord, i_offset, i_spin) + &
!!$                               mixing_factor(i_store) * &
!!$                               previous_rho_gradient_diff(i_coord, i_offset,i_store, i_spin)
!!$                       enddo
!!$                    end do
!!$                 enddo     !     next coordinate [of x,y,z]
!!$              end if
!!$           end do          !     end loop over a batch
!!$        ! end if             !     end work distribution over threads
!!$     end do                !     end loop over batches

    ! AJL, Feb2018
    ! This should just use the pulay infrastructure, as the two are so similar.
    if (use_meta_gga) then
     ! We need to kinetic_density for both residuals and positions for this to
     ! work.
     ! add together all pieces that might require preconditioning
     i_offset = 0
     do i_my_batch = 1, n_my_batches, 1
           do i_index = 1, batches(i_my_batch)%size, 1
              i_offset = i_offset + 1
              if (partition_tab(i_offset).gt.0.d0.or.(atomic_MERM .and. mpb_solver_started)) then
                    do i_spin = 1, n_spin, 1
                       kinetic_density_diff(i_offset, i_spin) = &
                            sigma * charge_mix_param(i_spin) * delta_kinetic_density(i_offset, i_spin)
                       if (saved_iter.gt.0) then
                          kinetic_density_diff(i_offset, i_spin) = &
                          kinetic_density_diff(i_offset, i_spin) - &
                          mixing_factor(1) * charge_mix_param(i_spin) * &
                          ( previous_kinetic_density_diff(i_offset,1,i_spin) + &
                          sigma * (delta_kinetic_density(i_offset,i_spin) - &
                          previous_kinetic_density_error(i_offset,1,i_spin)))
                       end if
                       do i_store = 2, saved_iter, 1
                          kinetic_density_diff(i_offset, i_spin) = &
                          kinetic_density_diff(i_offset, i_spin) - &
                          mixing_factor(i_store) * charge_mix_param(i_spin) * &
                          ( previous_kinetic_density_diff(i_offset,i_store,i_spin) + &
                          sigma * (previous_kinetic_density_error(i_offset,i_store-1,i_spin) - &
                          previous_kinetic_density_error(i_offset,i_store,i_spin)))
                       enddo
                    end do
              end if
           end do          !     end loop over a batch
     end do                !     end loop over batches

     ! precondition gradient, if required
     if (use_kerker_preconditioner.and.kerker_preconditioner_on .and. &
     &   .not.precondition_before_mixer) then
           do i_spin = 1, n_spin
              call precondition_kerker(kinetic_density_diff(:,i_spin), hartree_partition_tab)
           end do
     end if
    endif ! use_meta_gga
   endif ! use_density_gradient

end subroutine broyden_mix_p1
!******		
