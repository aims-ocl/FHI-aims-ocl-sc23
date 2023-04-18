!****h* FHI-aims/mixing
!  NAME
!    mixing
!  SYNOPSIS
      module mixing
!  PURPOSE
!    Module mixing handles the details of Pulay mixing etc.
!
!    IMPORTANT - SOME MIXING SUBROUTINES ARE NOT PART OF THE MODULE BECAUSE THEY WORK
!    AS GENERIC "TOOLS" WITH DIFFERENT INPUT VARIABLES - HENCE THEY NEED TO BE
!    SELF-CONTAINED ANYWAY AND IT MAKES NO SENSE TO ADD THEM TO THE MODULE
!
!    Subroutines inside the module:
!    * allocate_pulay
!    * prepare_pulay_mixing
!    * execute_pulay_mixing
!    * cleanup_pulay
!
!    Utility subroutines OUTSIDE THE MODULE (Separate files!!):
!    * pulay_store
!    * pulay_grad_store
!    * pulay_mix
!
!    FIXME: This file needs tidying - there are lots of repeating for
!    loops that can be concatenated. Also, is broyden_upate just a copy
!    of pulay_update? This could also be tidied. AJL, Feb 2018
!
!  USES
      use dimensions
      use runtime_choices
      use mpi_tasks
      use localorb_io
      use lpb_solver_utilities, only: atomic_MERM, mpb_solver_started
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

      implicit none

!     global variable declarations - exported to other program parts

!     rho_diff: Difference between input and output (mixed) density of present iteration:
!               rho_diff = rho_mix - rho_in
!     delta_rho_KS : If linear mixing and Pulay mixing requested later, store
!           rho_KS[i] - rho_in[i] for present generation in this array
!
!     Only if density gradients are needed:
!     rho_gradient_diff :  for our form of Pulay mixing, need
!              difference between present and previous charge density gradients
!     delta_rho_gradient : store grad ( rho_KS[i] - rho_in[i] ) for present generation
!     delta_rho_gradient_prec : these must be saved as well ...

      real*8, dimension(:,:), allocatable :: rho_diff
      real*8, dimension(:,:), allocatable :: delta_rho_KS

      real*8, dimension(:,:,:), allocatable :: rho_gradient_diff
      real*8, dimension(:,:,:), allocatable :: delta_rho_gradient

      ! AJL, Feb2018
      real*8, dimension(:,:), allocatable :: kinetic_density_diff
      real*8, dimension(:,:), allocatable :: delta_kinetic_density

!     These are to be used for pulay mixing of density matrix
!     see routine "density_matrix_mixing.f"
      integer :: pulay_saved_iter_denmat
      real*8, dimension(:), allocatable :: mixing_factor

!     stored plain mixing factor in case the preconditioner is used
      real*8 :: stored_charge_mix_param

!     private variable declarations - not exported to other program parts

!     pulay_saved_iter: In case of Pulay mixing, this is the number of
!        already stored previous densities
!
!     previous_rho_diff : Stored difference between input and output charge densities
!              from previous iterations n-i: rho_diff[i] = rho_mix[i] - rho_in[i]
!     previous_rho_error : Stored charge density errors delta_rho_KS from previous iterations
!     previous_rho_gradient_diff : Stored differences between input and output charge density gradients
!              from previous iterations n-i; see previous_rho_diff above. Only referenced if gradient
!              functionals required.
!     previous_rho_gradient_error : Stored charge density gradient errors delta_rho_gradient from previous iterations
!     previous_kinetic_density_diff: Stored difference between input and
!              output kinetic density
!     previous_kinetic_density_error:  Stored charge density gradient
!              errors (delta_kinetic_density) from previous iterations
!
!     pulay_matrix : Scalar products of density residuals for Pulay algorithm:
!                    int { r[rho_i](r)*r[rho_j](r) } d3r
!
!     broyden_matrix : Scalar products of density residuals for Pulay algorithm:
!                    int { r[rho_i](r)*r[rho_j](r) } d3r

      integer, private :: pulay_saved_iter

      real*8, dimension(:,:,:), allocatable, private :: previous_rho_diff
      real*8, dimension(:,:,:), allocatable, private :: previous_rho_error
      real*8, dimension(:,:), allocatable, private :: pulay_matrix
      real*8, dimension(:,:), allocatable, private :: broyden_matrix
       real*8, dimension(:,:,:,:), allocatable, private :: &
                previous_rho_gradient_diff
      real*8, dimension(:,:,:,:), allocatable, private :: &
                previous_rho_gradient_error
      ! AJL, Feb2018
      real*8, dimension(:,:,:), allocatable, private :: &
              previous_kinetic_density_diff
      real*8, dimension(:,:,:), allocatable, private :: &
              previous_kinetic_density_error


      ! if this flag is set, the charge density from the first iteration
      ! (scf initialization = superposition of free atom densities) is not
      ! included in the charge density mixing procedure
      logical :: first_iter_mixing

      ! Multipliers used to threshold charge density ... only needed in case of linear mixing
      ! and subsequent Pulay mixing ...
      real*8, dimension(:), allocatable, private :: charge_multiplier

      ! For localorb_io based output
      character*100, private :: info_str
!******
      contains
!----------------------------------------------------------------------------
!****s* mixing/allocate_mixing
!  NAME
!    allocate_mixing
!  SYNOPSIS
        subroutine allocate_mixing
!  PURPOSE
!    allocate variables necessary for mixing
!  INPUT
!  none
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
        implicit none


        if (.not.allocated(linear_mix_param)) then
          allocate(linear_mix_param(n_spin))
        end if
        if (.not.allocated(charge_mix_param)) then
          allocate(charge_mix_param(n_spin))
        end if
        if (.not.allocated(prec_mix_param)) then
          allocate(prec_mix_param(n_spin))
        end if

        if (use_mixer_threshold) then
          if (.not.allocated(max_rho_change)) then
            allocate(max_rho_change(n_spin))
          end if
          if (.not.allocated(charge_multiplier)) then
            allocate(charge_multiplier(n_spin))
          end if
        end if

        end subroutine allocate_mixing
!******
!----------------------------------------------------------------------------
!****s* mixing/allocate_delta_rho
!  NAME
!    allocate_delta_rho
!  SYNOPSIS
        subroutine allocate_delta_rho ( )
! PURPOSE
!  Subroutine allocate_delta_rho allocates the necessary storage arrays
!  for all mixing
!  INPUT
!  none
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
        implicit none



        if (.not.allocated(delta_rho_KS)) then
                allocate( delta_rho_KS(n_full_points, n_spin) )
        end if

        if (use_density_gradient) then
!         if we need density gradients, must store them for Pulay mixing, too
          if (.not.allocated(delta_rho_gradient)) then
           allocate( delta_rho_gradient(3,n_full_points,n_spin) )
          end if
        else
!         only add dummy allocation
          if (.not.allocated(delta_rho_gradient)) then
           allocate( delta_rho_gradient(1,1,1) )
          end if
        end if

        if (use_meta_gga) then
!         if we need kinetic_density, must store them for Pulay mixing, too
          if (.not.allocated(delta_kinetic_density)) then
           allocate( delta_kinetic_density(n_full_points,n_spin) )
          end if
        else
!         only add dummy allocation
          if (.not.allocated(delta_kinetic_density)) then
           allocate( delta_kinetic_density(1,1) )
          end if
        end if

        end subroutine allocate_delta_rho
!******
!----------------------------------------------------------------------------
!****s* mixing/allocate_pulay
!  NAME
!    allocate_pulay
!  SYNOPSIS
        subroutine allocate_pulay ( )
          !  PURPOSE
          !    Subroutine allocate_pulay allocates the necessary storage arrays
          !    for Pulay mixing.
          !  INPUT
          !  none
          !  OUTPUT
          !  none
          !  AUTHOR
          !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
          !  HISTORY
          !    Release version, FHI-aims (2008).
          !  SOURCE
          implicit none

          !       Local variables

          integer(kind=8) :: storage_size

          !       counters

          integer :: i_store



! when we are allocating the Pulay mixing arrays, must also reset pulay_saved_iter
          pulay_saved_iter = 0
          pulay_saved_iter_denmat = 0

!       The first few storage arrays can not easily be swapped out to disk -
!       we will have to live with them

          if (.not.allocated(rho_diff)) then
             allocate( rho_diff(n_full_points, n_spin) )
             rho_diff = 0.d0
          end if

          call allocate_delta_rho()

          if (.not.allocated(pulay_matrix)) then
             allocate(pulay_matrix(n_max_pulay,n_max_pulay))
          end if
          if (use_density_gradient) then
!            if we need density gradients, must store them for Pulay mixing, too
             if (.not.allocated(rho_gradient_diff)) then
                allocate( rho_gradient_diff(3,n_full_points,n_spin) )
                rho_gradient_diff = 0.d0
             end if
          else
             ! dummy allocation only
             if (.not.allocated(rho_gradient_diff)) then
                allocate( rho_gradient_diff(1,1,1) )
             end if
          end if

          if (use_meta_gga) then
!            if we need kinetic density, must store them for Pulay mixing, too
             if (.not.allocated(kinetic_density_diff)) then
                allocate( kinetic_density_diff(n_full_points,n_spin) )
                kinetic_density_diff = 0.d0
             end if
          else
             ! dummy allocation only
             if (.not.allocated(kinetic_density_diff)) then
                allocate( kinetic_density_diff(1,1) )
             end if
          end if

!       The following big storage arrays could be swapped out to disk if
!       really needed


          if (.not.allocated(previous_rho_diff)) then
             allocate(previous_rho_diff &
             (n_full_points, n_max_pulay, n_spin))
             previous_rho_diff = 0.d0
          end if
          if (.not.allocated(previous_rho_error)) then
             allocate(previous_rho_error &
             (n_full_points, n_max_pulay, n_spin))
             previous_rho_error = 0.d0
          end if

          if (use_density_gradient) then
             if (.not.allocated(previous_rho_gradient_diff)) then
                allocate( previous_rho_gradient_diff &
                (3, n_full_points, n_max_pulay, n_spin) )
             end if
             if (.not.allocated(previous_rho_gradient_error)) then
                allocate( previous_rho_gradient_error &
                (3, n_full_points, n_max_pulay, n_spin) )
             end if
          else
             ! dummy allocations only
             if (.not.allocated(previous_rho_gradient_diff)) then
                allocate( previous_rho_gradient_diff &
                (1, 1, 1, 1) )
             end if
             if (.not.allocated(previous_rho_gradient_error)) then
                allocate( previous_rho_gradient_error &
                (1, 1, 1, 1) )
             end if
          end if

          if (use_meta_gga) then
             if (.not.allocated(previous_kinetic_density_diff)) then
                allocate( previous_kinetic_density_diff &
                (n_full_points, n_max_pulay, n_spin) )
             end if
             if (.not.allocated(previous_kinetic_density_error)) then
                allocate( previous_kinetic_density_error &
                (n_full_points, n_max_pulay, n_spin) )
             end if
          else
             ! dummy allocations only
             if (.not.allocated(previous_kinetic_density_diff)) then
                allocate( previous_kinetic_density_diff &
                (1, 1, 1) )
             end if
             if (.not.allocated(previous_kinetic_density_error)) then
                allocate( previous_kinetic_density_error &
                (1, 1, 1) )
             end if
          end if

          if(.not.allocated(mixing_factor)) then
             allocate(mixing_factor(n_max_pulay))
             mixing_factor = 0.d0
          endif

        end subroutine allocate_pulay
!******
!----------------------------------------------------------------------------
!****s* mixing/allocate_broyden
!  NAME
!    allocate_broyden
!  SYNOPSIS
        subroutine allocate_broyden ( )
          !  PURPOSE
          !    Subroutine allocate_broyden allocates the necessary storage arrays
          !    for Broyden mixing.
          !  INPUT
          !  none
          !  OUTPUT
          !  none
          !  AUTHOR
          !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
          !  HISTORY
          !    Release version, FHI-aims (2008).
          !  SOURCE
          implicit none

          !       Local variables

          integer(kind=8) :: storage_size

          !       counters

          integer :: i_store



! when we are allocating the Broyden mixing arrays, must also reset pulay_saved_iter
          pulay_saved_iter = 0
          pulay_saved_iter_denmat = 0

!       The first few storage arrays can not easily be swapped out to disk -
!       we will have to live with them

          if (.not.allocated(rho_diff)) then
             allocate( rho_diff(n_full_points, n_spin) )
             rho_diff = 0.d0
          end if

          call allocate_delta_rho()

          if (.not.allocated(broyden_matrix)) then
             allocate(broyden_matrix(n_max_broyden,n_max_broyden))
             broyden_matrix = 0.d0
          end if
          if (use_density_gradient) then
!            if we need density gradients, must store them for Broyden mixing, too
             if (.not.allocated(rho_gradient_diff)) then
                allocate( rho_gradient_diff(3,n_full_points,n_spin) )
                rho_gradient_diff = 0.d0
             end if
          else
             ! dummy allocation only
             if (.not.allocated(rho_gradient_diff)) then
                allocate( rho_gradient_diff(1,1,1) )
             end if
          end if

          if (use_meta_gga) then
!            if we need density gradients, must store them for Pulay mixing, too
             if (.not.allocated(kinetic_density_diff)) then
                allocate( kinetic_density_diff(n_full_points,n_spin) )
                kinetic_density_diff = 0.d0
             end if
          else
             ! dummy allocation only
             if (.not.allocated(kinetic_density_diff)) then
                allocate( kinetic_density_diff(1,1) )
             end if
          end if

!       The following big storage arrays could be swapped out to disk if
!       really needed


          if (.not.allocated(previous_rho_diff)) then
             allocate(previous_rho_diff &
             (n_full_points, n_max_broyden, n_spin))
             previous_rho_diff = 0.d0
          end if
          if (.not.allocated(previous_rho_error)) then
             allocate(previous_rho_error &
             (n_full_points, n_max_broyden, n_spin))
             previous_rho_error = 0.d0
          end if

          if (use_density_gradient) then
             if (.not.allocated(previous_rho_gradient_diff)) then
                allocate( previous_rho_gradient_diff &
                (3, n_full_points, n_max_broyden, n_spin) )
             end if
             if (.not.allocated(previous_rho_gradient_error)) then
                allocate( previous_rho_gradient_error &
                (3, n_full_points, n_max_broyden, n_spin) )
             end if
          else
             ! dummy allocations only
             if (.not.allocated(previous_rho_gradient_diff)) then
                allocate( previous_rho_gradient_diff &
                (1, 1, 1, 1) )
             end if
             if (.not.allocated(previous_rho_gradient_error)) then
                allocate( previous_rho_gradient_error &
                (1, 1, 1, 1) )
             end if
          end if

          if (use_meta_gga) then
             if (.not.allocated(previous_kinetic_density_diff)) then
                allocate( previous_kinetic_density_diff &
                (n_full_points, n_max_pulay, n_spin) )
             end if
             if (.not.allocated(previous_kinetic_density_error)) then
                allocate( previous_kinetic_density_error &
                (n_full_points, n_max_pulay, n_spin) )
             end if
          else
             ! dummy allocations only
             if (.not.allocated(previous_kinetic_density_diff)) then
                allocate( previous_kinetic_density_diff &
                (1, 1, 1) )
             end if
             if (.not.allocated(previous_kinetic_density_error)) then
                allocate( previous_kinetic_density_error &
                (1, 1, 1) )
             end if
          end if

          if(.not.allocated(mixing_factor)) then
             allocate(mixing_factor(n_max_broyden))
             mixing_factor = 0.d0
          endif

        end subroutine allocate_broyden
!*****
!----------------------------------------------------------------------------
!****s* mixing/prepare_pulay_mixing
!  NAME
!    prepare_pulay_mixing
!  SYNOPSIS
        subroutine prepare_pulay_mixing ( )
          !  PURPOSE
          !    Subroutine prepare_pulay is only needed when we're in the linear mixing
          !    phase, but need to store densities for later iterations with Pulay mixing.
          !  INPUT
          !  none
          !  OUTPUT
          !  none
          !  USES
          use mpi_tasks
          !  AUTHOR
          !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
          !  HISTORY
          !    Release version, FHI-aims (2008).
          !  SOURCE

          !  Global input variables

          !  Local variables

          integer :: mpierr

          !     Counters

          integer :: i_spin
          integer :: i_store

          !  begin work




          if (pulay_saved_iter.lt.(n_max_pulay)) then
             !        we did not yet store the requested maximum number of densities
             !        for the Pulay scheme
             pulay_saved_iter = pulay_saved_iter+1
          end if

          do i_spin = 1, n_spin, 1

             ! VB: rho_diff is needed for later Pulay mixing stage.
             !     In the old linear mixing, it was calculated with the density itself;
             !     in the new variant, this is no longer the case. hence, calculate it outrightly here.
             rho_diff(:,i_spin) = linear_mix_param(i_spin) * &
             delta_rho_KS(:,i_spin)

             call pulay_store &
             ( rho_diff(:,i_spin), pulay_saved_iter, &
             previous_rho_diff(:,:,i_spin) )

             ! VB: Yuck. If charge density change was limited in limit_rho_change, must restore it here
             !     before storing it. What a complicated beast.
             if ( use_mixer_threshold ) then
                if ( charge_multiplier(i_spin).gt.0.d0 ) then
                   ! delta_rho_KS was throttled - restore
                   delta_rho_KS(:,i_spin) = &
                   delta_rho_KS(:,i_spin) / charge_multiplier(i_spin)
                end if
             end if

             call pulay_store &
             ( delta_rho_KS(:,i_spin), &
             pulay_saved_iter, &
             previous_rho_error(:,:,i_spin) )

             if (use_density_gradient) then

                ! VB: rho_diff is needed for later Pulay mixing stage.
                !     In the old linear mixing, it was calculated with the density itself;
                !     in the new variant, this is no longer the case. hence, calculate it outrightly here.
                rho_gradient_diff(:,:,i_spin) = &
                linear_mix_param(i_spin) * delta_rho_gradient(:,:,i_spin)

                call pulay_grad_store &
                ( rho_gradient_diff(:,:,i_spin), &
                pulay_saved_iter, &
                previous_rho_gradient_diff(:,:,:,i_spin) &
                )

                ! VB: Yuck. If charge density change was limited in limit_rho_change, must restore it here
                !     before storing it. What a complicated beast.
                if ( use_mixer_threshold ) then
                   if ( charge_multiplier(i_spin).gt.0.d0 ) then
                      ! delta_rho_KS was throttled - restore
                      delta_rho_gradient(:,:,i_spin) = &
                      delta_rho_gradient(:,:,i_spin) / &
                      charge_multiplier(i_spin)
                   end if
                end if

                call pulay_grad_store &
                ( delta_rho_gradient(:,:,i_spin), &
                pulay_saved_iter, &
                previous_rho_gradient_error(:,:,:,i_spin) &
                )

               if (use_meta_gga) then

                 ! VB: kinetic_density_diff is needed for later Pulay mixing stage.
                 kinetic_density_diff(:,i_spin) = linear_mix_param(i_spin) * &
                 delta_kinetic_density(:,i_spin)

                 call pulay_store &
                 ( kinetic_density_diff(:,i_spin), pulay_saved_iter, &
                 previous_kinetic_density_diff(:,:,i_spin) )

                 ! VB: Yuck. If kinetic density change was limited in limit_rho_change, must restore it here
                 !     before storing it. What a complicated beast.
                 if ( use_mixer_threshold ) then
                    if ( charge_multiplier(i_spin).gt.0.d0 ) then
                       ! delta_kinetic_density was throttled - restore
                       delta_kinetic_density(:,i_spin) = &
                       delta_kinetic_density(:,i_spin) / charge_multiplier(i_spin)
                    end if
                 end if

                 call pulay_store &
                 ( delta_kinetic_density(:,i_spin), &
                 pulay_saved_iter, &
                 previous_kinetic_density_error(:,:,i_spin) )

              endif !use_meta_gga
            endif !use_density_gradient

          enddo


        end subroutine prepare_pulay_mixing
!******
!----------------------------------------------------------------------------
!****s* mixing/linear_mix_p1
!  NAME
!    linear_mix_p1
!  SYNOPSIS

      subroutine linear_mix_p1 ( partition_tab, hartree_partition_tab, &
       rho, rho_gradient, kinetic_density )

! PURPOSE
!   Subroutine linear_mix executes a very simple version of linear charge density mixing
! USES

      use grids
      use geometry
      use mpi_utilities
      use runtime_choices
      use precondition
      implicit none

!  ARGUMENTS

      real*8, dimension(n_full_points) :: partition_tab
      real*8, dimension(n_full_points) :: hartree_partition_tab
      real*8, dimension (n_spin, n_full_points) :: rho
      real*8, dimension (3, n_spin, n_full_points) :: rho_gradient
      real*8, dimension (n_spin, n_full_points) :: kinetic_density

!  INPUTS
!    o partition_tab -- integration weights
!    o hartree_partition_tab -- integration weights for hartree potential, required for preconditioner
!    o rho -- density
!    o rho_gradient -- density gradient
!    o kinetic_density -- 2\tau
!  OUTPUTS
!    o rho -- mixed density
!    o rho_gradient -- mixed density gradient
!    o kinetic_density -- mixed 2\tau
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



      ! temporary variable to store delta_rho_KS AND delta_rho_gradient prior to any potential preconditioning
      real*8, dimension (:,:), allocatable :: rho_temp

      ! Local variables
      real*8 :: local_rho(n_spin)
      real*8 :: local_rho_gradient(3,n_spin)
      real*8 :: local_kinetic_density(n_spin)

      integer :: mixer_type

      integer :: i_offset
      integer :: i_my_batch
      integer :: i_index

      integer :: i_coord, i_spin

      ! Begin work

      call localorb_info( &
           "Linear mixing of updated and previous charge densities.", &
           use_unit,'(2X,A)',OL_norm)

      if (.not.allocated(rho_temp)) then
         allocate(rho_temp(n_full_points, n_spin))
      end if

!     If required, limit the allowed charge density step length
      if ( use_mixer_threshold ) then
         mixer_type = 0
         call limit_rho_change &
         ( partition_tab, mixer_type, delta_rho_KS, delta_rho_gradient, &
           delta_kinetic_density, charge_multiplier )
      end if

      ! store delta_rho_KS in rho_temp:
      do i_spin = 1, n_spin
         rho_temp(:,i_spin) = delta_rho_KS(:,i_spin)
      end do

      ! precondition charge density residual if requested ...
      if (use_kerker_preconditioner.and.kerker_preconditioner_on .and. &
      &   .not.precondition_before_mixer) then
         do i_spin = 1, n_spin
            call precondition_kerker(rho_temp(:,i_spin), &
                 hartree_partition_tab)
         end do
      else if (use_dielectric_preconditioner) then
         do i_spin = 1, n_spin
            call precondition_dielectric(rho_temp(:, i_spin), partition_tab, hartree_partition_tab)
         end do
      end if

      ! do charge density mixing with rho_temp, while PRESERVING delta_rho_KS
      ! for an eventual build of the proper pulay matrix
      i_offset   = 0
      do i_my_batch = 1, n_my_batches, 1
            do i_index = 1, batches(i_my_batch)%size, 1
               i_offset = i_offset + 1
               if (partition_tab(i_offset).gt.0.d0 .or. (atomic_MERM .and. mpb_solver_started)) then
                  if (spin_treatment.eq.0) then
                     i_spin = 1
                     rho (i_spin, i_offset) =  rho (i_spin, i_offset) + &
                          linear_mix_param(i_spin) * &
                          rho_temp(i_offset,i_spin)
                  else if (spin_treatment.eq.1) then
                     local_rho(1) = 0.5d0 * &
                      ( linear_mix_param(1) * rho_temp(i_offset,1) + &
                        linear_mix_param(2) * rho_temp(i_offset,2) )
                     local_rho(2) = 0.5d0 * &
                      ( linear_mix_param(1) * rho_temp(i_offset,1) - &
                        linear_mix_param(2) * rho_temp(i_offset,2))
                     do i_spin = 1, n_spin, 1
                        rho (i_spin, i_offset) = &
                             rho (i_spin, i_offset) + local_rho(i_spin)
                     enddo
                  end if
               end if
!           end loop over a batch
            end do
!         end work distribution over threads
         ! end if
!       end loop over batches
      end do

!       if required, play exactly the same mixing game for the gradients:
      if (use_density_gradient) then
         ! do the coordinate loop external to avoid allocating an array that's 3x rho_temp - use rho_temp instead
         do i_coord = 1, 3
            ! first, calculate rho_temp - as an exact copy of rho_gradient, which can be preconditioned
            do i_spin = 1, n_spin
               rho_temp(:,i_spin) = &
                    delta_rho_gradient(i_coord,:,i_spin)
            end do

            ! precondition charge density gradient if necessary:
            if (use_kerker_preconditioner.and. &
                    kerker_preconditioner_on .and. &
                      .not.precondition_before_mixer) then
               do i_spin = 1, n_spin
                  call precondition_kerker( &
                       rho_temp(:,i_spin), &
                       hartree_partition_tab)
               end do
            else if (use_dielectric_preconditioner) then
               do i_spin = 1, n_spin
                  call precondition_dielectric(rho_temp(:, i_spin), partition_tab, hartree_partition_tab)
               end do
            end if

            ! now do the actual mixing using rho_temp, while preserving the original delta_rho_gradient
            i_offset = 0
            do i_my_batch = 1, n_my_batches, 1
                  do i_index = 1, batches(i_my_batch)%size, 1
                     i_offset = i_offset + 1
                     if (partition_tab(i_offset).gt.0.d0 .or. (atomic_MERM .and. mpb_solver_started)) then
                        if (spin_treatment.eq.0) then
                           i_spin = 1
                           rho_gradient (i_coord, i_spin, i_offset) = &
                              rho_gradient(i_coord, i_spin, i_offset) + &
                                linear_mix_param(i_spin) * &
                              rho_temp(i_offset,i_spin)
                        else if (spin_treatment.eq.1) then
                           local_rho_gradient(i_coord,1) = &
                                0.5d0 * &
                                ( linear_mix_param(1) * &
                                rho_temp(i_offset,1) + &
                                linear_mix_param(2) * &
                                rho_temp(i_offset,2))
                           local_rho_gradient(i_coord,2) = &
                                0.5d0 * &
                                ( linear_mix_param(1) * &
                                rho_temp(i_offset,1) - &
                                linear_mix_param(2) * &
                                rho_temp(i_offset,2))

                           do i_spin = 1, n_spin, 1
                              rho_gradient (i_coord, i_spin, i_offset) = &
                              rho_gradient (i_coord, i_spin, i_offset) + &
                                   local_rho_gradient(i_coord, i_spin)
                           enddo
                        end if
                     end if
                  end do   ! end loop over a batch
               ! end if      ! end work distribution over threads
            end do         ! end loop over batches
         end do            ! loop over coordinates

         if (use_meta_gga) then
           do i_spin = 1, n_spin
             rho_temp(:,i_spin) = delta_kinetic_density(:,i_spin)
          end do

          ! precondition charge density residual if requested ...
          if (use_kerker_preconditioner.and.kerker_preconditioner_on .and. &
          &   .not.precondition_before_mixer) then
             do i_spin = 1, n_spin
                call precondition_kerker(rho_temp(:,i_spin), &
                     hartree_partition_tab)
             end do
          else if (use_dielectric_preconditioner) then
             do i_spin = 1, n_spin
                call precondition_dielectric(rho_temp(:, i_spin), partition_tab, hartree_partition_tab)
             end do
          end if

          ! do kinetic_density mixing with rho_temp, while PRESERVING
          ! delta_kinetic_density (is this necessary?)
          i_offset   = 0
          do i_my_batch = 1, n_my_batches, 1
                do i_index = 1, batches(i_my_batch)%size, 1
                   i_offset = i_offset + 1
                   if (partition_tab(i_offset).gt.0.d0 .or. (atomic_MERM .and. mpb_solver_started)) then
                      if (spin_treatment.eq.0) then
                         i_spin = 1
                         kinetic_density (i_spin, i_offset) = &
                            kinetic_density (i_spin, i_offset) + &
                            linear_mix_param(i_spin) * &
                            rho_temp(i_offset,i_spin)
                    else if (spin_treatment.eq.1) then
                       ! Reusing local_rho here
                       local_kinetic_density(1) = 0.5d0 * &
                        ( linear_mix_param(1) * rho_temp(i_offset,1) + &
                          linear_mix_param(2) * rho_temp(i_offset,2) )
                       local_kinetic_density(2) = 0.5d0 * &
                        ( linear_mix_param(1) * rho_temp(i_offset,1) - &
                          linear_mix_param(2) * rho_temp(i_offset,2) )
                       do i_spin = 1, n_spin, 1
                          kinetic_density (i_spin, i_offset) = &
                          kinetic_density (i_spin, i_offset) + &
                          local_kinetic_density(i_spin)
                       enddo
                    end if
                 end if
!             end loop over a batch
              end do
!           end work distribution over threads
           ! end if
!         end loop over batches
        end do

        endif ! use_meta_gga
      endif ! use_density_gradient

      if (allocated (rho_temp)) deallocate(rho_temp)

      end subroutine linear_mix_p1
!******
!------------------------------------------------------------------------------------------
!****s* mixing/execute_pulay_mixing
!  NAME
!    execute_pulay_mixing
!  SYNOPSIS
       subroutine execute_pulay_mixing &
        ( number_of_loops, partition_tab, hartree_partition_tab, &
         rho, rho_gradient, kinetic_density)
!  PURPOSE
!    main driver routine for the Pulay mixer
!  USES
        use mpi_tasks
        use synchronize_mpi
!  ARGUMENTS
        integer :: number_of_loops
        real*8, dimension(n_full_points) :: partition_tab
        real*8, dimension(n_full_points) :: hartree_partition_tab
        real*8, dimension (n_spin, n_full_points) :: rho
        real*8, dimension (3, n_spin, n_full_points) :: rho_gradient
        real*8, dimension (n_spin, n_full_points) :: kinetic_density
!  INPUTS
!    o number_of_loops -- number of Pulay mixing steps
!    o partition tab -- integration weights
!    o hartree_partition_tab -- integration weights for hartree potential (for preconditioner!)
!    o rho -- density
!    o rho_gradient -- density gradient
!    o kinetic_density -- 2\tau
!  OUTPUTS
!    o rho -- mixed density
!    o rho_gradient -- mixed density gradient
!    o kinetic_density -- \2tau
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      integer :: mixer_type
      real*8 :: diff_mixing
      integer :: mpierr

!     Counters

      integer :: i_spin
      integer :: i_store
      integer :: i_pulay_store
!  begin work

      if ( ( number_of_loops.eq.(ini_linear_mixing+1)) .and. &
           ( number_of_loops.gt.1) ) then
!                 This is the iteration where we switched from linear to
!                 Pulay mixing. Must initialize the Pulay matrix.



! -- FIXME : Why was this ever needed?                  pulay_saved_iter = max(pulay_saved_iter, n_max_pulay)

            call initialize_pulay_matrix_p1 &
                 ( partition_tab, pulay_saved_iter, previous_rho_error, &
                 pulay_matrix)


            call sync_pulay_matrix( pulay_matrix )

      end if

!               now, actual Pulay mixing ...

                  ! use old infrastructure for mixing, spin loops outside

         ! pulay for one matrix containing both charge and spin residua
         ! basically the matrix is just the sum of the former charge and spin pulay matrix
         i_spin = 1
         call pulay_mix_p1 &
              ( partition_tab, &           ! real*8 (n_full_points) &
              hartree_partition_tab, &      ! real*8 (n_full_points) &
              pulay_saved_iter,       &    ! integer &
              delta_rho_KS,            &   ! real*8 (n_full_points, n_spin) &
              rho_diff,                 &  ! real*8 (n_full_points, n_spin) &
              previous_rho_diff,        &  ! real*8 (n_full_points, n_max_pulay, n_spin) &
              previous_rho_error,       &  ! real*8 (n_full_points, n_max_pulay, n_spin) &
              delta_rho_gradient,       &  ! real*8 (3, n_full_points, n_spin) &
              rho_gradient_diff,         & ! real*8 (3, n_full_points, n_spin) &
              previous_rho_gradient_diff,& ! real*8 (3, n_full_points, n_max_pulay, n_spin) &
              previous_rho_gradient_error,&! real*8 (3, n_full_points, n_max_pulay, n_spin) &
              delta_kinetic_density,            &   ! real*8 (n_full_points, n_spin) &
              kinetic_density_diff,                 &  ! real*8 (n_full_points, n_spin) &
              previous_kinetic_density_diff,        &  ! real*8 (n_full_points, n_max_pulay, n_spin) &
              previous_kinetic_density_error,       &  ! real*8 (n_full_points, n_max_pulay, n_spin) &
              pulay_matrix,           &    ! real*8 (n_max_pulay,n_max_pulay) &
              mixing_factor &
              )


!     If required, limit the allowed charge density step length
      if ( use_mixer_threshold ) then
        mixer_type = 1
        call limit_rho_change &
        ( partition_tab, mixer_type, rho_diff, rho_gradient_diff, &
          kinetic_density_diff, charge_multiplier )
      end if

!     Pulay_update uses previous_rho_[gradient]_diff to
!     update rho [, rho_gradient]
!
!     Notice that previous_rho(i_spin) stores
!     (rho_up + rho_dn) if i_spin = 1
!     (rho_up - rho_dn) if i_spin = 2
!
!     i.e. for the spin-polarized case we need
!     rho(up) = 1/2 * ( previous_rho(1)+previous_rho(2) )
!     rho(dn) = 1/2 * ( previous_rho(1)-previous_rho(2) )

            call pulay_update_p1 &
                 ( partition_tab, rho, rho_gradient, kinetic_density )


!     Finally, store this iteration's densities for future use ...
                  ! use old infrastructure for mixing, spin loops outside

         pulay_saved_iter_denmat = pulay_saved_iter
         if (pulay_saved_iter.lt.(n_max_pulay)) then
!                   we did not yet store the requested maximum number of densities
!                   for the Pulay scheme
            pulay_saved_iter = pulay_saved_iter+1
         end if

         do i_spin = 1, n_spin, 1

            call pulay_store &
                 ( delta_rho_KS(:,i_spin), &
                 pulay_saved_iter, &
                 previous_rho_error(:,:,i_spin) )

            call pulay_store &
                 ( rho_diff(:,i_spin), &
                 pulay_saved_iter, &
                 previous_rho_diff(:,:,i_spin) )


            if (use_density_gradient) then

               call pulay_grad_store &
                    ( delta_rho_gradient(:,:,i_spin), &
                    pulay_saved_iter, &
                    previous_rho_gradient_error(:,:,:,i_spin) )

               call pulay_grad_store &
                    ( rho_gradient_diff(:,:,i_spin), &
                    pulay_saved_iter, &
                    previous_rho_gradient_diff(:,:,:,i_spin) )

               if (use_meta_gga) then

                 call pulay_store &
                      ( delta_kinetic_density(:,i_spin), &
                      pulay_saved_iter, &
                      previous_kinetic_density_error(:,:,i_spin) )

                 call pulay_store &
                      ( kinetic_density_diff(:,i_spin), &
                      pulay_saved_iter, &
                      previous_kinetic_density_diff(:,:,i_spin) )

               end if ! use_meta_gga
            endif ! use_density_gradient

         enddo

      end subroutine execute_pulay_mixing
!******
!------------------------------------------------------------------------------------------
!****s* mixing/execute_broyden_mixing
!  NAME
!    execute_broyden_mixing
!  SYNOPSIS
       subroutine execute_broyden_mixing &
        ( number_of_loops, partition_tab, hartree_partition_tab, &
         rho, rho_gradient, kinetic_density)
!  PURPOSE
!    main driver routine for the broyden mixer
!  USES
        use mpi_tasks
        use synchronize_mpi
!  ARGUMENTS
        integer :: number_of_loops
        real*8, dimension(n_full_points) :: partition_tab
        real*8, dimension(n_full_points) :: hartree_partition_tab
        real*8, dimension (n_spin, n_full_points) :: rho
        real*8, dimension (3, n_spin, n_full_points) :: rho_gradient
        real*8, dimension (n_spin, n_full_points) :: kinetic_density
!  INPUTS
!    o number_of_loops -- number of broyden mixing steps
!    o partition tab -- integration weights
!    o hartree_partition_tab -- integration weights for hartree potential (for preconditioner!)
!    o rho -- density
!    o rho_gradient -- density gradient
!    o kinetic_density -- 2\tau
!  OUTPUTS
!    o rho -- mixed density
!    o rho_gradient -- mixed density gradient
!    o kinetic_density -- mixed 2\tau
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      integer :: mixer_type
      real*8 :: diff_mixing
      integer :: mpierr

!     Counters
      integer :: i_spin
      integer :: i_store
      integer :: i_pulay_store
!  begin work

         ! broyden for one matrix containing both charge and spin residua
         ! basically the matrix is just the sum of the former charge and spin pulay matrix
         call broyden_mix_p1 &
              ( partition_tab, &           ! real*8 (n_full_points) &
              hartree_partition_tab, &      ! real*8 (n_full_points) &
              pulay_saved_iter,       &    ! integer &
              delta_rho_KS,            &   ! real*8 (n_full_points, n_spin) &
              rho_diff,                 &  ! real*8 (n_full_points, n_spin) &
              previous_rho_diff,        &  ! real*8 (n_full_points, n_max_broyden, n_spin) &
              previous_rho_error,       &  ! real*8 (n_full_points, n_max_broyden, n_spin) &
              delta_rho_gradient,       &  ! real*8 (3, n_full_points, n_spin) &
              rho_gradient_diff,         & ! real*8 (3, n_full_points, n_spin) &
              previous_rho_gradient_diff,& ! real*8 (3, n_full_points, n_max_broyden, n_spin) &
              previous_rho_gradient_error,&! real*8 (3, n_full_points, n_max_broyden, n_spin) &
              delta_kinetic_density,            &   ! real*8 (n_full_points, n_spin) &
              kinetic_density_diff,                 &  ! real*8 (n_full_points, n_spin) &
              previous_kinetic_density_diff,        &  ! real*8 (n_full_points, n_max_broyden, n_spin) &
              previous_kinetic_density_error,       &  ! real*8 (n_full_points, n_max_broyden, n_spin) &
              broyden_matrix,           & ! real*8 (n_max_broyden,n_max_broyden) &
              mixing_factor &
              )

!     Broyden_update uses previous_rho_[gradient]_diff to
!     update rho [, rho_gradient]
!
!     Notice that previous_rho(i_spin) stores
!     (rho_up + rho_dn) if i_spin = 1
!     (rho_up - rho_dn) if i_spin = 2
!
!     i.e. for the spin-polarized case we need
!     rho(up) = 1/2 * ( previous_rho(1)+previous_rho(2) )
!     rho(dn) = 1/2 * ( previous_rho(1)-previous_rho(2) )

            call broyden_update_p1 &
                 ( partition_tab, rho, rho_gradient, kinetic_density )

!     Finally, store this iteration's densities for future use ...
                  ! use old infrastructure for mixing, spin loops outside

         pulay_saved_iter_denmat = pulay_saved_iter
         if (pulay_saved_iter.lt.(n_max_broyden-1)) then
!                   we did not yet store the requested maximum number of densities
!                   for the Pulay scheme
            pulay_saved_iter = pulay_saved_iter+1
         end if

         do i_spin = 1, n_spin, 1

! Broyden_store is redundant, as it is the same subroutine
! as pulay_store, so why not use the same routine for both?
!            call broyden_store &
             call pulay_store &
                 ( delta_rho_KS(:,i_spin), &
                 pulay_saved_iter, &
                 previous_rho_error(:,:,i_spin) )

!            call broyden_store &
             call pulay_store &
                 ( rho_diff(:,i_spin), &
                 pulay_saved_iter, &
                 previous_rho_diff(:,:,i_spin) )


            if (use_density_gradient) then

               call pulay_grad_store &
                    ( delta_rho_gradient(:,:,i_spin), &
                    pulay_saved_iter, &
                    previous_rho_gradient_error(:,:,:,i_spin) )

               call pulay_grad_store &
                    ( rho_gradient_diff(:,:,i_spin), &
                    pulay_saved_iter, &
                    previous_rho_gradient_diff(:,:,:,i_spin) )

               if (use_meta_gga) then

                call pulay_store &
                   ( delta_kinetic_density(:,i_spin), &
                   pulay_saved_iter, &
                   previous_kinetic_density_error(:,:,i_spin) )

                call pulay_store &
                   ( kinetic_density_diff(:,i_spin), &
                   pulay_saved_iter, &
                   previous_kinetic_density_diff(:,:,i_spin) )

               endif ! use_meta_gga

            endif ! use_density_gradient

         enddo

      end subroutine execute_broyden_mixing
!******
!------------------------------------------------------------------------------------------
!****s* mixing/pulay_update_p1
!  NAME
!    pulay_update_p1
!  SYNOPSIS
      subroutine pulay_update_p1 ( partition_tab, rho, rho_gradient, &
                                   kinetic_density )
! PURPOSE
!  Subroutine pulay_update uses previous_rho_[gradient]_diff to
!  update rho [, rho_gradient]
!
!  Notice that previous_rho(i_spin) stores
!     (rho_up + rho_dn) if i_spin = 1
!     (rho_up - rho_dn) if i_spin = 2
!
!  i.e. for the spin-polarized case we need
!      rho(up) = 1/2 * ( previous_rho(1)+previous_rho(2) )
!      rho(dn) = 1/2 * ( previous_rho(1)-previous_rho(2) )
! USES
      use grids
      use geometry
      use mpi_utilities
      implicit none
!  ARGUMENTS
      real*8, dimension(n_full_points) :: partition_tab
      real*8, dimension (n_spin, n_full_points) :: rho
      real*8, dimension (3, n_spin, n_full_points) :: rho_gradient
      real*8, dimension (n_spin, n_full_points) :: kinetic_density
!  INPUTS
!   o partition_tab -- integration weights
!   o rho -- density
!   o rho_gradient -- density gradient
!   o kinetic_density -- 2\tau
!  OUTPUTS
!   o rho -- mixed density
!   o rho_gradient -- mixed density gradient
!   o kinetic_density -- 2\tau
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
      real*8, dimension(n_spin) :: local_rho
      real*8, dimension(3, n_spin) :: local_rho_gradient

!     counters
      integer :: i_index
      integer :: i_coord
      integer :: i_spin
      integer :: i_offset
      integer :: i_my_batch

!  begin work

      ! FIXME? Is there any reason these loops can't be concatenated
      ! into one? AJL, Feb2018

      i_offset = 0

      do i_my_batch = 1, n_my_batches, 1

            do i_index = 1, batches(i_my_batch)%size, 1

               i_offset = i_offset + 1
!           execute only if partition_tab.gt.0 here, i.e. if the integration point
!           makes sense
               if (partition_tab(i_offset).gt.0.d0.or.(atomic_MERM .and. mpb_solver_started)) then

                  if (spin_treatment.eq.0) then
                     i_spin = 1

!               do nothing special, just use i_spin = 1
                     rho (i_spin, i_offset) =  rho (i_spin, i_offset) + &
                          rho_diff(i_offset,i_spin)

                  else

                ! spin-up density difference
                     local_rho(1) = &
                          0.5d0*( rho_diff(i_offset,1) &
                          + rho_diff(i_offset,2) )
                ! spin-down density difference
                     local_rho(2) = &
                          0.5d0*( rho_diff(i_offset,1) &
                          - rho_diff(i_offset,2) )

                     do i_spin = 1, n_spin, 1
                        rho (i_spin, i_offset) = &
                             rho (i_spin, i_offset) + local_rho(i_spin)
                     enddo

                  end if

               end if

!     end loop over a batch
            end do
!     end work distribution over threads
         ! end if
!     end loop over batches
      end do

!     if required, play exactly the same mixing game for the gradients

      if (use_density_gradient) then
         i_offset = 0

         do i_my_batch = 1, n_my_batches, 1

               do i_index = 1, batches(i_my_batch)%size, 1

                  i_offset = i_offset + 1
!             execute only if partition_tab.gt.0 here, i.e. if the integration point
!             makes sense
                  if (partition_tab(i_offset).gt.0.d0.or.(atomic_MERM .and. mpb_solver_started)) then

                     if (spin_treatment.eq.0) then
                        i_spin = 1

!                 do nothing special, just use i_spin = 1
                        do i_coord = 1,3,1

                           rho_gradient(i_coord, i_spin, i_offset) = &
                                rho_gradient(i_coord, i_spin, i_offset) &
                                + rho_gradient_diff &
                                (i_coord, i_offset,i_spin)

!                 next coordinate [of x,y,z]
                        enddo

                     else

                        do i_coord = 1,3,1

                    ! spin-up density difference
                           local_rho_gradient(i_coord,1) = &
                                0.5d0*( rho_gradient_diff &
                                (i_coord,i_offset,1) &
                                + rho_gradient_diff(i_coord,i_offset,2))
                    ! spin-down density difference
                           local_rho_gradient(i_coord,2) = &
                                0.5d0*( rho_gradient_diff &
                                (i_coord,i_offset,1) &
                                - rho_gradient_diff(i_coord,i_offset,2))

                           do i_spin = 1, n_spin, 1
                              rho_gradient(i_coord,i_spin,i_offset) = &
                                   rho_gradient &
                                   (i_coord, i_spin, i_offset)+ &
                                   local_rho_gradient(i_coord,i_spin)
                           enddo

                        enddo

                     end if

                  end if

!     end loop over a batch
               end do
!     end work distribution over threads
            ! end if
!     end loop over batches
         end do

         if (use_meta_gga) then
         i_offset = 0

         do i_my_batch = 1, n_my_batches, 1

            do i_index = 1, batches(i_my_batch)%size, 1

               i_offset = i_offset + 1
!           execute only if partition_tab.gt.0 here, i.e. if the integration point
!           makes sense
               if (partition_tab(i_offset).gt.0.d0.or.(atomic_MERM .and. mpb_solver_started)) then

                  if (spin_treatment.eq.0) then
                     i_spin = 1

!               do nothing special, just use i_spin = 1
                     kinetic_density (i_spin, i_offset) = &
                     kinetic_density (i_spin, i_offset) + &
                          kinetic_density_diff(i_offset,i_spin)

                  else

                ! spin-up density difference
                     local_rho(1) = &
                          0.5d0*( kinetic_density_diff(i_offset,1) &
                          + kinetic_density_diff(i_offset,2) )
                ! spin-down density difference
                     local_rho(2) = &
                          0.5d0*( kinetic_density_diff(i_offset,1) &
                          - kinetic_density_diff(i_offset,2) )

                     do i_spin = 1, n_spin, 1
                        kinetic_density (i_spin, i_offset) = &
                             kinetic_density (i_spin, i_offset) + local_rho(i_spin)
                     enddo

                  end if

               end if

!     end loop over a batch
            end do
!     end loop over batches
         end do

         endif ! use_meta_gga

      endif ! use_density_gradient

      end subroutine pulay_update_p1
!******
!------------------------------------------------------------------------------------------
!****s* mixing/broyden_update_p1
!  NAME
!    broyden_update_p1
!  SYNOPSIS
      subroutine broyden_update_p1 ( partition_tab, rho, rho_gradient, &
                                     kinetic_density)
! PURPOSE
!  Subroutine broyden_update uses previous_rho_[gradient]_diff to
!  update rho [, rho_gradient]
!
!  Notice that previous_rho(i_spin) stores
!     (rho_up + rho_dn) if i_spin = 1
!     (rho_up - rho_dn) if i_spin = 2
!
!  i.e. for the spin-polarized case we need
!      rho(up) = 1/2 * ( previous_rho(1)+previous_rho(2) )
!      rho(dn) = 1/2 * ( previous_rho(1)-previous_rho(2) )
!
! AJL: Is this identical to pulay_update_p1? If so, can be removed...
! USES
      use grids
      use geometry
      use mpi_utilities
      implicit none
!  ARGUMENTS
      real*8, dimension(n_full_points) :: partition_tab
      real*8, dimension (n_spin, n_full_points) :: rho
      real*8, dimension (3, n_spin, n_full_points) :: rho_gradient
      real*8, dimension (n_spin, n_full_points) :: kinetic_density
!  INPUTS
!   o partition_tab -- integration weights
!   o rho -- density
!   o rho_gradient -- density gradient
!   o kinetic_density -- 2\tau
!  OUTPUTS
!   o rho -- mixed density
!   o rho_gradient -- mixed density gradient
!   o kinetic_density -- mixed 2\tau
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
      real*8, dimension(n_spin) :: local_rho
      real*8, dimension(3, n_spin) :: local_rho_gradient

!     counters
      integer :: i_index
      integer :: i_coord
      integer :: i_spin
      integer :: i_offset
      integer :: i_my_batch

!  begin work
      i_offset = 0

      do i_my_batch = 1, n_my_batches, 1

            do i_index = 1, batches(i_my_batch)%size, 1

               i_offset = i_offset + 1
!           execute only if partition_tab.gt.0 here, i.e. if the integration point
!           makes sense
               if (partition_tab(i_offset).gt.0.d0.or.(atomic_MERM .and. mpb_solver_started)) then

                  if (spin_treatment.eq.0) then
                     i_spin = 1

!               do nothing special, just use i_spin = 1
                     rho (i_spin, i_offset) =  rho (i_spin, i_offset) + &
                          rho_diff(i_offset,i_spin)

                  else

                ! spin-up density difference
                     local_rho(1) = &
                          0.5d0*( rho_diff(i_offset,1) &
                          + rho_diff(i_offset,2) )
                ! spin-down density difference
                     local_rho(2) = &
                          0.5d0*( rho_diff(i_offset,1) &
                          - rho_diff(i_offset,2) )

                     do i_spin = 1, n_spin, 1
                        rho (i_spin, i_offset) = &
                             rho (i_spin, i_offset) + local_rho(i_spin)
                     enddo

                  end if

               end if

!     end loop over a batch
            end do
!     end work distribution over threads
         ! end if
!     end loop over batches
      end do

!     if required, play exactly the same mixing game for the gradients

      if (use_density_gradient) then
         i_offset = 0

         do i_my_batch = 1, n_my_batches, 1

               do i_index = 1, batches(i_my_batch)%size, 1

                  i_offset = i_offset + 1
!             execute only if partition_tab.gt.0 here, i.e. if the integration point
!             makes sense
                  if (partition_tab(i_offset).gt.0.d0.or.(atomic_MERM .and. mpb_solver_started)) then

                     if (spin_treatment.eq.0) then
                        i_spin = 1

!                 do nothing special, just use i_spin = 1
                        do i_coord = 1,3,1

                           rho_gradient(i_coord, i_spin, i_offset) = &
                                rho_gradient(i_coord, i_spin, i_offset) &
                                + rho_gradient_diff &
                                (i_coord, i_offset,i_spin)

!                 next coordinate [of x,y,z]
                        enddo

                     else

                        do i_coord = 1,3,1

                    ! spin-up density difference
                           local_rho_gradient(i_coord,1) = &
                                0.5d0*( rho_gradient_diff &
                                (i_coord,i_offset,1) &
                                + rho_gradient_diff(i_coord,i_offset,2))
                    ! spin-down density difference
                           local_rho_gradient(i_coord,2) = &
                                0.5d0*( rho_gradient_diff &
                                (i_coord,i_offset,1) &
                                - rho_gradient_diff(i_coord,i_offset,2))

                           do i_spin = 1, n_spin, 1
                              rho_gradient(i_coord,i_spin,i_offset) = &
                                   rho_gradient &
                                   (i_coord, i_spin, i_offset)+ &
                                   local_rho_gradient(i_coord,i_spin)
                           enddo

                        enddo

                     end if

                  end if

!     end loop over a batch
               end do
!     end work distribution over threads
            ! end if
!     end loop over batches
         end do

        if (use_meta_gga) then

        i_offset = 0

        do i_my_batch = 1, n_my_batches, 1

              do i_index = 1, batches(i_my_batch)%size, 1

                 i_offset = i_offset + 1
!             execute only if partition_tab.gt.0 here, i.e. if the integration point
!             makes sense
                 if (partition_tab(i_offset).gt.0.d0.or.(atomic_MERM .and. mpb_solver_started)) then

                    if (spin_treatment.eq.0) then
                       i_spin = 1

!                 do nothing special, just use i_spin = 1
                       kinetic_density (i_spin, i_offset) = &
                       kinetic_density (i_spin, i_offset) + &
                            kinetic_density_diff(i_offset,i_spin)

                    else

                  ! spin-up density difference
                       local_rho(1) = &
                            0.5d0*( kinetic_density_diff(i_offset,1) &
                            + kinetic_density_diff(i_offset,2) )
                  ! spin-down density difference
                       local_rho(2) = &
                            0.5d0*( kinetic_density_diff(i_offset,1) &
                            - kinetic_density_diff(i_offset,2) )

                       do i_spin = 1, n_spin, 1
                          kinetic_density (i_spin, i_offset) = &
                          kinetic_density (i_spin, i_offset) + local_rho(i_spin)
                       enddo

                    end if

                 end if

!       end loop over a batch
              end do
!       end loop over batches
        end do

        endif ! use_meta_gga
      endif ! use_density_gradient

      end subroutine broyden_update_p1
!******
!------------------------------------------------------------------------------------------
!****s* mixing/cleanup_pulay
!  NAME
!    cleanup_pulay
!  SYNOPSIS
        subroutine cleanup_pulay ( )
!  PURPOSE
!  Subroutine cleanup_pulay deallocates all storage to do with pulay mixing
        implicit none
!  INPUTS
!  none
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
        if (allocated(rho_diff)) then
          deallocate( rho_diff )
        end if
        if (allocated(delta_rho_KS)) then
          deallocate( delta_rho_KS )
        end if
        if (allocated(previous_rho_diff)) then
          deallocate( previous_rho_diff )
        end if
        if (allocated(previous_rho_error)) then
          deallocate( previous_rho_error )
        end if

        if (allocated(rho_gradient_diff)) then
          deallocate( rho_gradient_diff )
        end if
        if (allocated(delta_rho_gradient)) then
          deallocate( delta_rho_gradient )
        end if
        if (allocated(previous_rho_gradient_diff)) then
          deallocate( previous_rho_gradient_diff )
        end if
        if (allocated(previous_rho_gradient_error)) then
          deallocate( previous_rho_gradient_error )
        end if

        if (allocated(kinetic_density_diff)) then
          deallocate( kinetic_density_diff )
        end if
        if (allocated(delta_kinetic_density)) then
          deallocate( delta_kinetic_density )
        end if
        if (allocated(previous_kinetic_density_diff)) then
          deallocate( previous_kinetic_density_diff )
        end if
        if (allocated(previous_kinetic_density_error)) then
          deallocate( previous_kinetic_density_error )
        end if

        if (allocated(pulay_matrix)) then
          deallocate( pulay_matrix )
        end if

        end subroutine cleanup_pulay
!******
!------------------------------------------------------------------------------------------
!****s* mixing/cleanup_broyden
!  NAME
!    cleanup_broyden
!  SYNOPSIS
        subroutine cleanup_broyden ( )
!  PURPOSE
!  Subroutine cleanup_broyden deallocates all storage to do with broyden mixing
        implicit none
!  INPUTS
!  none
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
        if (allocated(rho_diff)) then
          deallocate( rho_diff )
        end if
        if (allocated(delta_rho_KS)) then
          deallocate( delta_rho_KS )
        end if
        if (allocated(previous_rho_diff)) then
          deallocate( previous_rho_diff )
        end if
        if (allocated(previous_rho_error)) then
          deallocate( previous_rho_error )
        end if

        if (allocated(rho_gradient_diff)) then
          deallocate( rho_gradient_diff )
        end if
        if (allocated(delta_rho_gradient)) then
          deallocate( delta_rho_gradient )
        end if
        if (allocated(previous_rho_gradient_diff)) then
          deallocate( previous_rho_gradient_diff )
        end if
        if (allocated(previous_rho_gradient_error)) then
          deallocate( previous_rho_gradient_error )
        end if

        if (allocated(kinetic_density_diff)) then
          deallocate( kinetic_density_diff )
        end if
        if (allocated(delta_kinetic_density)) then
          deallocate( delta_kinetic_density )
        end if
        if (allocated(previous_kinetic_density_diff)) then
          deallocate( previous_kinetic_density_diff )
        end if
        if (allocated(previous_kinetic_density_error)) then
          deallocate( previous_kinetic_density_error )
        end if

        if (allocated(broyden_matrix)) then
          deallocate( broyden_matrix )
        end if

        end subroutine cleanup_broyden
!******
!------------------------------------------------------------------------------------------
!****s* mixing/cleanup_mixing
!  NAME
!    cleanup_mixing
!  SYNOPSIS
        subroutine cleanup_mixing
!  PURPOSE
!  Subroutine cleanup_pulay deallocates all storage to do with the remainder of the mixing module
!
!  INPUT
!  none
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        if (allocated(linear_mix_param)) then
          deallocate(linear_mix_param)
        end if
        if (allocated(charge_mix_param)) then
          deallocate(charge_mix_param)
        end if
        if (allocated(prec_mix_param)) then
          deallocate(prec_mix_param)
        end if
        if (allocated(max_rho_change)) then
          deallocate(max_rho_change)
        end if
        if (allocated(charge_multiplier)) then
          deallocate(charge_multiplier)
        end if
        if(allocated(mixing_factor)) then
             deallocate(mixing_factor)
        endif
        end subroutine cleanup_mixing
!******
      end module mixing
