      module mixing_constraint

      ! MODULES - NOTE - FURTHER MODULES ARE USED DOWN BELOW BY IND. SUBROUTINES
      use dimensions
      use constraint
      use runtime_choices

      implicit none

!     global variable declarations - exported to other program parts

      real*8, dimension(:,:), allocatable :: potential_diff

!     private variable declarations - not exported to other program parts

      integer, private :: pulay_saved_iter_constraint

      real*8, dimension(:,:,:), allocatable, private :: &
           previous_potential_diff
      real*8, dimension(:,:,:), allocatable, private :: &
           previous_potential_error
      real*8, dimension(:,:,:,:), allocatable, private :: &
           pulay_matrix_constraint

      contains
!---------------------------------------------------------------------
!  Subroutine allocate_pulay allocates the necessary storage arrays
!  for Pulay mixing.

      subroutine allocate_pulay_constraint &
           ( &
           )

      implicit none

!       Local variables


!       counters

      integer :: i_store


      if (.not.allocated(potential_diff)) then
         allocate( potential_diff(n_region, n_spin) )
         potential_diff = 0.d0
      end if

      if (.not.allocated(pulay_matrix_constraint)) then
         allocate(pulay_matrix_constraint(n_max_pulay_constraint, &
              n_max_pulay_constraint,n_region,n_spin) )
         pulay_matrix_constraint = 0.0d0
      end if

      if (.not.allocated(previous_potential_diff)) then
         allocate(previous_potential_diff &
              (n_region, n_max_pulay_constraint, n_spin))
         previous_potential_diff = 0.d0
      end if
      if (.not.allocated(previous_potential_error)) then
         allocate(previous_potential_error &
              (n_region, n_max_pulay_constraint, n_spin))
         previous_potential_error = 0.d0
      end if

      pulay_saved_iter_constraint = 0

      end subroutine allocate_pulay_constraint
!---------------------------------------------------------------------
!  Subroutine prepare_pulay is only neede when we're in the linear mixing
!  phase, but need to store densities for later iterations with Pulay mixing.

      subroutine prepare_pulay_mixing_constraint &
           ( delta_potential_in_mixing )

      use mpi_tasks

!  Input variables

      real*8 :: delta_potential_in_mixing(n_region,n_spin)

!  Local variables

!     Counters

      integer :: i_spin
      integer :: i_region

!  begin work


      if (pulay_saved_iter_constraint.lt.(n_max_pulay_constraint)) then
!     we did not yet store the requested maximum number of densities
!     for the Pulay scheme
         pulay_saved_iter_constraint = pulay_saved_iter_constraint+1
      end if

      do i_spin = 1, n_spin, 1

         do i_region = 1, n_active_regions, 1

            call pulay_store_constraint &
                 ( potential_diff(i_region,i_spin), &
                 pulay_saved_iter_constraint, &
                 previous_potential_diff(i_region,:,i_spin) )

            call pulay_store_constraint &
                 ( delta_potential_in_mixing(i_region,i_spin), &
                 pulay_saved_iter_constraint, &
                 previous_potential_error(i_region,:,i_spin) )
         enddo

      enddo

      end subroutine prepare_pulay_mixing_constraint
!---------------------------------------------------------------------
!  Subroutine pulay_store stores the present density rho
!  in a packed array of previous densities for later use by the
!  Pulay mixer

      subroutine pulay_store_constraint &
           ( potential, pulay_saved_iter_constraint, &
           previous_potential &
           )

      use dimensions

      implicit none

!  imported variables

!     input

      real*8 :: potential
      integer :: pulay_saved_iter_constraint

!     output

      real*8, dimension(n_max_pulay_constraint) :: previous_potential

!     local variables

!     counters

      integer :: i_pulay_store

!  begin work

!     first, shift stored densities to make room for next density
      do i_pulay_store = pulay_saved_iter_constraint, 2, -1

         previous_potential(i_pulay_store) = &
              previous_potential(i_pulay_store-1)

      end do

!     next, store present density in first column of previous_rho
      previous_potential(1) = potential

!  end work

      end subroutine pulay_store_constraint

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!  Subroutine execute_pulay_mixing is a wrapper around the pulay
!  mixing utilities which executes the mixing itself.

      subroutine execute_pulay_mixing_constraint &
           ( number_constraint_iter, constraint_potential, &
           delta_potential_in_mixing )


!  Global input variables

        ! partition_tab is only needed to make sure that we don't store
        ! any quantities on irrelevant integration points
      integer :: number_constraint_iter
      real*8, dimension (n_region, n_spin) :: constraint_potential
      real*8, dimension (n_region, n_spin) :: delta_potential_in_mixing

!  Local variables


!     Counters

      integer :: i_spin
      integer :: i_region

!  begin work

      if ( ( number_constraint_iter.eq.(ini_linear_mixing_constraint+1)) &
           .and. ( number_constraint_iter.gt.1) ) then
!     This is the iteration where we switched from linear to
!     Pulay mixing. Must initialize the Pulay matrix.

         call initialize_pulay_matrix_constraint &
              ( )

      end if

!               now, actual Pulay mixing ...

!     pulay_mix only updates rho_diff and gradient_rho_diff -
!     but does not know what they contain,
!     and in particular does not know anything about spin polarisation
      call pulay_mix_constraint &
           ( delta_potential_in_mixing )

        ! simply update on all threads
      call pulay_update_constraint &
           ( constraint_potential )

      if (pulay_saved_iter_constraint.lt.(n_max_pulay_constraint)) then
!     we did not yet store the requested maximum number of densities
!     for the Pulay scheme
         pulay_saved_iter_constraint = pulay_saved_iter_constraint+1
      end if

      do i_spin = 1, n_spin, 1

         do i_region = 1, n_active_regions, 1

            call pulay_store_constraint &
                 ( potential_diff(i_region,i_spin), &
                 pulay_saved_iter_constraint, &
                 previous_potential_diff(i_region,:,i_spin) )

            call pulay_store_constraint &
                 ( delta_potential_in_mixing(i_region,i_spin), &
                 pulay_saved_iter_constraint, &
                 previous_potential_error(i_region,:,i_spin) )
         enddo

      enddo

      end subroutine execute_pulay_mixing_constraint
!----------------------------------------------------------------------
!  Subroutine initialize_pulay_matrix
!
!  constructs the Pulay matrix from the previously stored densities
!  from linear mixing.
!
!----------------------------------------------------------------------
        subroutine initialize_pulay_matrix_constraint &
      ( )

      use dimensions
      use mpi_tasks

      implicit none

!  imported variables

!  input

!  output

!  counters

      integer :: i_spin, i_region
      integer :: i_pulay_store, i_pulay_store_2

!  begin work


      if (myid.eq.0) then
         write(use_unit,'(2X,A,A)') &
              "Switching from linear to Pulay mixing: ", &
              "Initializing Pulay matrix."
      end if

      pulay_matrix_constraint = 0.0d0

!     calculate, use, and store density error err(rho_i) for present iteration i
!     between KS input and output densities from previous iterations n-i, n-j ...

!     determine Pulay matrix elements: < err(rho_i)-err(rho_i+1) | err(rho_j)-err(rho_j+1) >

      do i_spin = 1, n_spin, 1
         do i_region = 1, n_active_regions, 1

            do i_pulay_store = 1, pulay_saved_iter_constraint-1, 1
               do i_pulay_store_2 = &
                    i_pulay_store, pulay_saved_iter_constraint-1, 1
                  pulay_matrix_constraint(i_pulay_store, &
                       i_pulay_store_2,i_region,i_spin) = &
                       pulay_matrix_constraint(i_pulay_store, &
                       i_pulay_store_2,i_region,i_spin) + &
                       ( previous_potential_error(i_region, &
                       i_pulay_store,i_spin) - &
                       previous_potential_error(i_region, &
                       i_pulay_store+1,i_spin) )* &
                       ( previous_potential_error(i_region, &
                       i_pulay_store_2,i_spin) - &
                       previous_potential_error(i_region, &
                       i_pulay_store_2+1,i_spin) )
               enddo
            enddo

         enddo
      enddo

      end subroutine initialize_pulay_matrix_constraint

!----------------------------------------------------------------------
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
!----------------------------------------------------------------------
      subroutine pulay_mix_constraint &
           ( &
           delta_potential_in_mixing &
           )

      use dimensions
      use lapack_wrapper
      use localorb_io

      use constants
      implicit none


!  imported variables

!  input

      real*8, dimension(n_region,n_spin) :: delta_potential_in_mixing

!  local variables

      real*8, dimension(pulay_saved_iter_constraint) :: pulay_vector
      real*8, dimension(pulay_saved_iter_constraint) :: mixing_factor

!  counters

      integer :: i_spin, i_region
      integer :: i_pulay_store, i_pulay_store_2


!  begin work

!      call localorb_info(
!     +"Pulay mixing of updated and previous constraint potentials.",
!     +     6,'(2X,A)' )

      do i_spin = 1, n_spin, 1
         do i_region = 1, n_active_regions, 1

!     initialize
            pulay_vector = 0.d0
            mixing_factor = 0.d0

            do i_pulay_store = pulay_saved_iter_constraint, 2, -1
               do i_pulay_store_2 = &
                    pulay_saved_iter_constraint, i_pulay_store, -1
                  pulay_matrix_constraint(i_pulay_store, &
                       i_pulay_store_2,i_region,i_spin) = &
                       pulay_matrix_constraint(i_pulay_store-1, &
                       i_pulay_store_2-1,i_region,i_spin)
               enddo
            enddo
            do i_pulay_store = 1, pulay_saved_iter_constraint, 1
               pulay_matrix_constraint(1,i_pulay_store,i_region,i_spin) &
                    = 0.d0
            enddo

!     initialize Pulay matrix elements of present iteration

!     calculate, use, and store density error err(rho_i) for present iteration i
!     between KS input and output densities from previous iterations n-i, n-j ...

!     determine Pulay vector: < err(rho_j+1) - err(rho_j) | err(rho_m) >
            if (pulay_saved_iter_constraint.gt.0) then

               pulay_matrix_constraint(1,1,i_region,i_spin) = &
                    pulay_matrix_constraint(1,1,i_region,i_spin) + &
                    ( delta_potential_in_mixing(i_region,i_spin) - &
                    previous_potential_error(i_region,1,i_spin) )**2.d0

               pulay_vector(1) = &
                    pulay_vector(1) + &
                    ( previous_potential_error(i_region,1,i_spin) - &
                    delta_potential_in_mixing(i_region,i_spin) ) * &
                    delta_potential_in_mixing(i_region,i_spin)
            end if

            do i_pulay_store = 2, pulay_saved_iter_constraint, 1

               pulay_matrix_constraint(1,i_pulay_store, &
                    i_region,i_spin) = &
                    pulay_matrix_constraint(1,i_pulay_store, &
                    i_region,i_spin) + &
                    ( delta_potential_in_mixing(i_region,i_spin) - &
                    previous_potential_error(i_region,1,i_spin) ) * &
                    ( previous_potential_error(i_region, &
                    i_pulay_store-1,i_spin) - &
                    previous_potential_error(i_region, &
                    i_pulay_store,i_spin) )

               pulay_vector(i_pulay_store) = &
                    pulay_vector(i_pulay_store) + &
                    (previous_potential_error(i_region, &
                    i_pulay_store,i_spin) - &
                    previous_potential_error(i_region, &
                    i_pulay_store-1,i_spin) ) * &
                    delta_potential_in_mixing(i_region,i_spin)

            enddo

!     Determine mixing factor from equation system defined by pulay_matrix
!     use the lapack expert solver for a symmetric indefinite problem
            if (pulay_saved_iter_constraint.gt.0) then

        ! It is possible that two successive iterations have a zero residual.
        ! In this case, we simply remove the appropriate column(s) / row(s)
        ! from the Pulay matrix. This should give a zero mixing factor for the
        ! density from that iteration. Is that the right thing to do?
               do i_pulay_store = 1, pulay_saved_iter_constraint, 1
                  if (pulay_matrix_constraint(i_pulay_store, &
                       i_pulay_store,i_region,i_spin).le.1.d-20) then
                     pulay_matrix_constraint(:,i_pulay_store, &
                          i_region,i_spin) = 0.d0
                     pulay_matrix_constraint(i_pulay_store,:, &
                          i_region,i_spin) = 0.d0
                     pulay_matrix_constraint(i_pulay_store, &
                          i_pulay_store,i_region,i_spin) = 1.d0
                     pulay_vector(i_pulay_store) = 0.d0
                  end if
               enddo


               call solve_pulay &
                    ( n_max_pulay_constraint, &
                    pulay_saved_iter_constraint, &
                    pulay_matrix_constraint(:,:,i_region,i_spin), &
                    pulay_vector, &
                    mixing_factor &
                    )

            end if

!     now mix densities

            potential_diff(i_region,i_spin) = constraint_mix(1) * &
                 delta_potential_in_mixing(i_region,i_spin)

            if (pulay_saved_iter_constraint.gt.0) then
               potential_diff(i_region,i_spin) = &
                    potential_diff(i_region,i_spin) + &
                    mixing_factor(1) * &
                    ( previous_potential_diff(i_region,1,i_spin) + &
                    constraint_mix(1) * &
                    ( delta_potential_in_mixing(i_region,i_spin) - &
                    previous_potential_error(i_region,1,i_spin) ) &
                    )
            end if

            do i_pulay_store = 2, pulay_saved_iter_constraint, 1
               potential_diff(i_region,i_spin) = &
                    potential_diff(i_region,i_spin) + &
                    mixing_factor(i_pulay_store) * &
                    ( previous_potential_diff(i_region, &
                    i_pulay_store,i_spin) + &
                    constraint_mix(1) * &
                    ( previous_potential_error(i_region, &
                    i_pulay_store-1,i_spin) - &
                    previous_potential_error(i_region, &
                    i_pulay_store,i_spin) ) &
                    )
            enddo

!     end loop over active regions
         enddo
!     end loop over spin
      enddo

      end subroutine pulay_mix_constraint
!---------------------------------------------------------------------
!     Subroutine pulay_update uses previous_rho_[gradient]_diff to
!     update rho [, rho_gradient]
!
!                 Notice that previous_rho(i_spin) stores
!     (rho_up + rho_dn) if i_spin = 1
!     (rho_up - rho_dn) if i_spin = 2
!
!     i.e. for the spin-polarized case we need
!     rho(up) = 1/2 * ( previous_rho(1)+previous_rho(2) )
!     rho(dn) = 1/2 * ( previous_rho(1)-previous_rho(2) )

      subroutine pulay_update_constraint &
           ( constraint_potential &
           )

      implicit none

!  imported variables

      real*8, dimension (n_region,n_spin) :: constraint_potential

!  local variables


!     counters

      integer :: i_region, i_spin

!  begin work

      do i_spin = 1, n_spin, 1
         do i_region = 1, n_active_regions, 1

            constraint_potential (i_region,i_spin) = &
                 constraint_potential (i_region,i_spin) + &
                 potential_diff(i_region,i_spin)

         enddo
      enddo

      end subroutine pulay_update_constraint
!---------------------------------------------------------------------
!  Subroutine cleanup_pulay deallocates all storage to do with pulay mixing
!
      subroutine cleanup_pulay_constraint &
           ( &
           )

      implicit none

!  begin work

      if (allocated(potential_diff)) then
         deallocate( potential_diff )
      end if

      if (allocated(previous_potential_diff)) then
         deallocate( previous_potential_diff )
      end if
      if (allocated(previous_potential_error)) then
         deallocate( previous_potential_error )
      end if

      if (allocated(pulay_matrix_constraint)) then
         deallocate( pulay_matrix_constraint )
      end if

      end subroutine cleanup_pulay_constraint
!---------------------------------------------------------------------
        end module mixing_constraint
