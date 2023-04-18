!****s* FHI-aims/pulay_mix_p1
!  NAME
!   pulay_mix_p1
!  SYNOPSIS

subroutine pulay_mix_p1 &
     ( partition_tab, hartree_partition_tab, &
     pulay_saved_iter, delta_rho_KS, rho_diff, &
     previous_rho_diff, previous_rho_error, &
     delta_rho_gradient, rho_gradient_diff, &
     previous_rho_gradient_diff, &
     previous_rho_gradient_error, &
     delta_kinetic_density, kinetic_density_diff, &
     previous_kinetic_density_diff, previous_kinetic_density_error, &
     pulay_matrix, mixing_factor )

!  PURPOSE
!  Subroutine pulay_mix
!
!  Use Pulay's algorithm to mix present density with the densities 
!  of the previous iterations. 
!
!  Notice that we employ the variant described by Kresse and Furthmueller,
!  Comp. Mat. Sci. 6, 15-50 (1996), on p. 34 of that paper. The recipe is 
!  described much clearer than the original version of Pulay. We implemented
!  the original Pulay recipe briefly, then switched to this one (equivalent?).
!
!  Gradients are pulled through in exactly the same way as the actual
!  density. If no density gradient is required, none of the relevant
!  arrays are actually allocated, i.e. never ever reference one of them.
!
!  Thanks go to Christoph Freysoldt, FHI, for an illuminating discussion
!  of the subject.
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
  use lpb_solver_utilities, only: atomic_MERM, mpb_solver_started
  implicit none

!  ARGUMENTS

  real*8, dimension(n_full_points) :: partition_tab
  real*8, dimension(n_full_points) :: hartree_partition_tab
  integer :: pulay_saved_iter

  real*8, dimension(n_full_points, n_spin)              :: delta_rho_KS
  real*8, dimension(n_full_points, n_spin)              :: rho_diff
  real*8, dimension(n_full_points, n_max_pulay, n_spin) :: previous_rho_diff
  real*8, dimension(n_full_points, n_max_pulay, n_spin) :: previous_rho_error
  
  real*8, dimension(3, n_full_points, n_spin)              :: delta_rho_gradient
  real*8, dimension(3, n_full_points, n_spin)              :: rho_gradient_diff
  real*8, dimension(3, n_full_points, n_max_pulay, n_spin) :: previous_rho_gradient_diff
  real*8, dimension(3, n_full_points, n_max_pulay, n_spin) :: previous_rho_gradient_error

  real*8, dimension(n_full_points, n_spin)              :: delta_kinetic_density
  real*8, dimension(n_full_points, n_spin)              :: kinetic_density_diff
  real*8, dimension(n_full_points, n_max_pulay, n_spin) :: previous_kinetic_density_diff
  real*8, dimension(n_full_points, n_max_pulay, n_spin) :: previous_kinetic_density_error

  real*8, dimension(n_max_pulay,n_max_pulay) :: pulay_matrix
  real*8, dimension(pulay_saved_iter)        :: mixing_factor

!  INPUTS
!  o partition_tab -- values of the partition function
!  o hartree_partition_tab -- values of the partition function for Hartree potential
!  o pulay_saved_iter -- number of saved iterations
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
!  o pulay_matrix -- ????????
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

  real*8, dimension(pulay_saved_iter) :: pulay_vector 

! counters
  integer :: i_index
  integer :: i_pulay_store, i_pulay_store_2
  integer :: i_coord
  integer :: i_offset
  integer :: i_my_batch
  integer :: i_spin

!  begin work
      
  call localorb_info( &
       "Pulay mixing of updated and previous charge densities.", &
       use_unit,'(2X,A)', OL_norm )

!     initialize
  pulay_vector = 0.d0
  mixing_factor = 0.d0

!    Prepare updated Pulay matrix
  do i_pulay_store = pulay_saved_iter, 2, -1
     do i_pulay_store_2 = pulay_saved_iter, i_pulay_store, -1
        pulay_matrix(i_pulay_store, i_pulay_store_2) =  &
             pulay_matrix(i_pulay_store-1, i_pulay_store_2-1)
     enddo
  enddo
!     initialize Pulay matrix elements of present iteration
  do i_pulay_store = 1, pulay_saved_iter, 1
     pulay_matrix(1,i_pulay_store) = 0.d0      
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

!             determine Pulay matrix elements: < err(rho_m+1) - err(rho_m) | err(rho_j+1) - err(rho_j) >
!             determine Pulay vector: < err(rho_j+1) - err(rho_j) | err(rho_m) >

              if (pulay_saved_iter.gt.0) then

                 do i_spin = 1, n_spin, 1
                    
                    pulay_matrix(1,1) = pulay_matrix(1,1) + &
                         partition_tab(i_offset) * &
                         ( delta_rho_ks(i_offset, i_spin) - previous_rho_error(i_offset,1, i_spin) )**2.d0
                    
                    pulay_vector(1) = pulay_vector(1) + &
                         partition_tab(i_offset) * &
                         (previous_rho_error(i_offset, 1, i_spin) - delta_rho_KS(i_offset, i_spin) ) * &
                         delta_rho_KS(i_offset, i_spin)
                    
                 end do
                 
              end if
              
              do i_pulay_store = 2, pulay_saved_iter, 1
                 
                 do i_spin = 1, n_spin, 1

                    pulay_matrix(1,i_pulay_store) =  &
                         pulay_matrix(1,i_pulay_store) + &
                         partition_tab(i_offset) * &
                         ( delta_rho_KS(i_offset, i_spin) - previous_rho_error(i_offset,1, i_spin) ) * &
                         ( previous_rho_error(i_offset, i_pulay_store-1, i_spin) - &
                         previous_rho_error(i_offset, i_pulay_store, i_spin) )
                    
                    pulay_vector(i_pulay_store) = &
                         pulay_vector(i_pulay_store) + &
                         partition_tab(i_offset) * &
                         (previous_rho_error(i_offset, i_pulay_store, i_spin) - &
                         previous_rho_error(i_offset, i_pulay_store-1, i_spin) ) * &
                         delta_rho_KS(i_offset, i_spin)
                 
                 end do
                 
              end do
                 
           end if

           !         end loop over a batch
        end do
        !     end work distribution over tasks
     ! end if
     !       end loop over batches
  end do

  if (pulay_saved_iter.gt.0) then
     call sync_pulay_matrix_first_row_and_vector( pulay_matrix, &
          pulay_vector, pulay_saved_iter )
  end if

!     Determine mixing factor from equation system defined by pulay_matrix
!     use the lapack expert solver for a symmetric indefinite problem
  if (pulay_saved_iter.gt.0) then

     ! It is possible that two successive iterations have a zero residual.
     ! In this case, we simply remove the appropriate column(s) / row(s) 
     ! from the Pulay matrix. This should give a zero mixing factor for the 
     ! density from that iteration. Is that the right thing to do?
     do i_pulay_store = 1, pulay_saved_iter, 1
        if (pulay_matrix(i_pulay_store,i_pulay_store).le.1.d-20) then
           pulay_matrix(:,i_pulay_store) = 0.d0
           pulay_matrix(i_pulay_store,:) = 0.d0
           pulay_matrix(i_pulay_store,i_pulay_store) = 1.d0
           pulay_vector(i_pulay_store) = 0.d0
        end if
     enddo

     call solve_pulay &
          ( n_max_pulay, pulay_saved_iter, pulay_matrix, pulay_vector, &
          mixing_factor )
  end if
      
!     now mix densities

!test
!  write(use_unit,*) "in pulay_mix_p2 (after solve_pulay) ..."

!  write(use_unit,*) "pulay-matrix..."

!  do i_pulay_store = pulay_saved_iter, 1, -1
!     do i_pulay_store_2 = pulay_saved_iter, i_pulay_store, -1
!        write(use_unit,'(E10.4, 2X)',advance='no') pulay_matrix(i_pulay_store, i_pulay_store_2) 
!     enddo
!     write(use_unit,*) 
!  enddo

!  write(use_unit,*) "pulay-vector..."
!  do i_pulay_store = pulay_saved_iter, 1, -1
!     write(use_unit,'(E10.4, 2X)') pulay_vector(i_pulay_store)
!  end do

!  write(use_unit,*) "mixing_factor..."
!  do i_pulay_store = pulay_saved_iter, 1, -1
!     write(use_unit,'(A, I2, 1X, I2, F20.4, 2X)') "mixing_factor ", i_spin, i_pulay_store, mixing_factor(i_pulay_store)
!  end do

!  write(use_unit,*) 
!test end

  ! put together the next density in two pieces:
  ! (1) the sum \Delta R_i terms 
  i_offset   = 0
  do i_my_batch = 1, n_my_batches, 1
        do i_index = 1, batches(i_my_batch)%size, 1
           i_offset = i_offset + 1
           if (partition_tab(i_offset).gt.0.d0 .or. (atomic_MERM .and. mpb_solver_started)) then
              do i_spin = 1, n_spin, 1
                 rho_diff(i_offset, i_spin) = charge_mix_param(i_spin) * delta_rho_KS(i_offset, i_spin)                           
                 if (pulay_saved_iter.gt.0) then
                    rho_diff(i_offset, i_spin) = rho_diff(i_offset, i_spin) + &
                         mixing_factor(1) * charge_mix_param(i_spin) * &
                         ( delta_rho_KS(i_offset, i_spin) - previous_rho_error(i_offset,1, i_spin) ) 
                 end if                 
                 do i_pulay_store = 2, pulay_saved_iter, 1
                    rho_diff(i_offset, i_spin) = rho_diff(i_offset, i_spin) + &
                         mixing_factor(i_pulay_store) * charge_mix_param(i_spin) * &
                         ( previous_rho_error(i_offset,i_pulay_store-1, i_spin) - &
                         previous_rho_error(i_offset,i_pulay_store, i_spin) ) 
                 enddo                 
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
  else if (use_dielectric_preconditioner) then
     do i_spin = 1, n_spin
        call precondition_dielectric(rho_diff(:, i_spin), partition_tab, hartree_partition_tab)
     end do
  end if

  ! Then add the actual 'density' terms - using the exact same looping structure
  i_offset   = 0
  do i_my_batch = 1, n_my_batches, 1
        do i_index = 1, batches(i_my_batch)%size, 1
           i_offset = i_offset + 1
           if (partition_tab(i_offset).gt.0.d0 .or. (atomic_MERM .and. mpb_solver_started)) then
              do i_spin = 1, n_spin, 1
                 if (pulay_saved_iter.gt.0) then
                    rho_diff(i_offset, i_spin) = rho_diff(i_offset, i_spin) + &
                         mixing_factor(1) * previous_rho_diff(i_offset,1, i_spin)
                 end if                 
                 do i_pulay_store = 2, pulay_saved_iter, 1
                    rho_diff(i_offset, i_spin) = rho_diff(i_offset, i_spin) + &
                         mixing_factor(i_pulay_store) * previous_rho_diff(i_offset,i_pulay_store, i_spin)
                 enddo                 
              end do
           end if
        end do        ! end loop over a batch
  end do              ! end loop over batches
 


!    if required, play exactly the same mixing game for the gradients

!     FIXME: this whole sequence needs a subroutine of its own to
!      avoid the silly code duplication between density and gradient; 
!      i.e. we could split the
!      update of the pulay matrix and the actual pulay mixing into
!     completely separate subroutines ...

  if (use_density_gradient) then
     ! add together all pieces that might require preconditioning
     i_offset = 0
     do i_my_batch = 1, n_my_batches, 1
           do i_index = 1, batches(i_my_batch)%size, 1
              i_offset = i_offset + 1
              if (partition_tab(i_offset).gt.0.d0 .or. (atomic_MERM .and. mpb_solver_started)) then
                 do i_coord = 1,3,1
                    do i_spin = 1, n_spin, 1
                       rho_gradient_diff(i_coord,i_offset, i_spin) = &
                            charge_mix_param(i_spin) * delta_rho_gradient(i_coord, i_offset, i_spin)
                       if (pulay_saved_iter.gt.0) then
                          rho_gradient_diff(i_coord, i_offset, i_spin) = &
                               rho_gradient_diff(i_coord, i_offset, i_spin) + &
                               mixing_factor(1) * &
                               charge_mix_param(i_spin) * &
                               ( delta_rho_gradient(i_coord, i_offset, i_spin) - &
                               previous_rho_gradient_error(i_coord,i_offset,1, i_spin)) 
                       end if
                       do i_pulay_store = 2, pulay_saved_iter, 1
                          rho_gradient_diff(i_coord, i_offset, i_spin) = &
                               rho_gradient_diff(i_coord, i_offset, i_spin) + &
                               mixing_factor(i_pulay_store) * &
                               charge_mix_param(i_spin) * &
                               ( previous_rho_gradient_error(i_coord, i_offset,i_pulay_store-1, i_spin) - &
                               previous_rho_gradient_error(i_coord, i_offset,i_pulay_store, i_spin) ) 
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
     else if (use_dielectric_preconditioner) then
        do i_coord = 1, 3
           do i_spin = 1, n_spin
              call precondition_dielectric(rho_gradient_diff(i_coord, :, i_spin), partition_tab, hartree_partition_tab)
           end do
        end do
     end if

     ! add remaining terms
     i_offset = 0
     do i_my_batch = 1, n_my_batches, 1
           do i_index = 1, batches(i_my_batch)%size, 1
              i_offset = i_offset + 1
              if (partition_tab(i_offset).gt.0.d0 .or. (atomic_MERM .and. mpb_solver_started)) then
                 do i_coord = 1,3,1
                    do i_spin = 1, n_spin, 1
                       if (pulay_saved_iter.gt.0) then
                          rho_gradient_diff(i_coord, i_offset, i_spin) = &
                               rho_gradient_diff(i_coord, i_offset, i_spin) + &
                               mixing_factor(1) * &
                               previous_rho_gradient_diff(i_coord, i_offset,1, i_spin)
                       end if
                       do i_pulay_store = 2, pulay_saved_iter, 1
                          rho_gradient_diff(i_coord, i_offset, i_spin) = &
                               rho_gradient_diff(i_coord, i_offset, i_spin) + &
                               mixing_factor(i_pulay_store) * &
                               previous_rho_gradient_diff(i_coord, i_offset,i_pulay_store, i_spin)
                       enddo
                    end do
                 enddo     !     next coordinate [of x,y,z]
              end if
           end do          !     end loop over a batch
     end do                !     end loop over batches

     if (use_meta_gga) then
     ! add together all pieces that might require preconditioning
     i_offset = 0
     do i_my_batch = 1, n_my_batches, 1
           do i_index = 1, batches(i_my_batch)%size, 1
              i_offset = i_offset + 1
              if (partition_tab(i_offset).gt.0.d0 .or. (atomic_MERM .and. mpb_solver_started)) then
                    do i_spin = 1, n_spin, 1
                       kinetic_density_diff(i_offset, i_spin) = &
                            charge_mix_param(i_spin) * delta_kinetic_density(i_offset, i_spin)
                       if (pulay_saved_iter.gt.0) then
                          kinetic_density_diff(i_offset, i_spin) = &
                               kinetic_density_diff(i_offset, i_spin) + &
                               mixing_factor(1) * &
                               charge_mix_param(i_spin) * &
                               ( delta_kinetic_density(i_offset, i_spin) - &
                               previous_kinetic_density_error(i_offset,1, i_spin))
                       end if
                       do i_pulay_store = 2, pulay_saved_iter, 1
                          kinetic_density_diff(i_offset, i_spin) = &
                               kinetic_density_diff(i_offset, i_spin) + &
                               mixing_factor(i_pulay_store) * &
                               charge_mix_param(i_spin) * &
                               ( previous_kinetic_density_error(i_offset,i_pulay_store-1, i_spin) - &
                               previous_kinetic_density_error(i_offset,i_pulay_store, i_spin) )
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
     else if (use_dielectric_preconditioner) then
           do i_spin = 1, n_spin
              call precondition_dielectric(kinetic_density_diff(:, i_spin), partition_tab, hartree_partition_tab)
           end do
     end if

     ! add remaining terms
     i_offset = 0
     do i_my_batch = 1, n_my_batches, 1
           do i_index = 1, batches(i_my_batch)%size, 1
              i_offset = i_offset + 1
              if (partition_tab(i_offset).gt.0.d0 .or. (atomic_MERM .and. mpb_solver_started)) then
                    do i_spin = 1, n_spin, 1
                       if (pulay_saved_iter.gt.0) then
                          kinetic_density_diff(i_offset, i_spin) = &
                               kinetic_density_diff(i_offset, i_spin) + &
                               mixing_factor(1) * &
                               previous_kinetic_density_diff(i_offset,1, i_spin)
                       end if
                       do i_pulay_store = 2, pulay_saved_iter, 1
                          kinetic_density_diff(i_offset, i_spin) = &
                               kinetic_density_diff(i_offset, i_spin) + &
                               mixing_factor(i_pulay_store) * &
                               previous_kinetic_density_diff(i_offset,i_pulay_store, i_spin)
                       enddo
                    end do
              end if
           end do          !     end loop over a batch
        ! end if             !     end work distribution over threads
     end do                !     end loop over batches
     end if ! use_meta_gga
   endif ! use_density_gradients

end subroutine pulay_mix_p1
!******		
